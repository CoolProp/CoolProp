# Residual-Helmholtz SIMD — "go crazy" continuation report

**Date:** 2026-05-30
**Branch:** `ihb/helmholtz-simd-crazy`
**Parent:** `ihb/helmholtz-simd-investigate` (which has the full prior story
in `docs/PERF-helmholtz-simd.md`)
**Time spent:** ~2-3 hours on this continuation
**Hardware:** Apple Silicon (ARM64, NEON 2-wide doubles, ~4 NEON FMA pipes)

## TL;DR

Two changes shipped on this branch:

1. **`allFastNEON` is now the default on Apple ARM64** (env opt-out via
   `COOLPROP_HELMHOLTZ_FAST=0`).  ~17% end-to-end speedup, ULP-class
   numerics, validated.
2. **Standalone polynomial-exp + vvexp benches saved to `dev/bench/`** —
   future-reference tools for the "try a custom poly exp" path.

One thing attempted that **didn't ship** (and shouldn't):

3. **`allFastNEONfused`** — a fully-fused single-pass NEON loop with
   inline polynomial exp.  Failed for two interlocking reasons (detailed
   below); the diff was discarded.

## What did ship: default-on dispatch

`src/Helmholtz.cpp:148` — the live `all()` now opens with:

```cpp
static const bool use_fast = []() {
    const char* env = std::getenv("COOLPROP_HELMHOLTZ_FAST");
    if (env != nullptr) return std::string(env) != "0";
#if defined(__APPLE__) && defined(__aarch64__)
    return true;   // Apple Silicon: opt-out via COOLPROP_HELMHOLTZ_FAST=0
#else
    return false;  // Other platforms: opt-in via env var
#endif
}();
if (use_fast) { allFastNEON(tau, delta, derivs); return; }
```

Static-local-init runs once at first call (C++11 guarantees thread-safe
single execution).  Subsequent calls read a single bool: ~1 ns/call,
unobservable in the noise.

This means **on Apple ARM64**, every CoolProp `all()` call now goes
through the 2-wide NEON SIMD B-chain by default, with the standalone
opt-out for debugging / regression bisection.

## What didn't ship: `allFastNEONfused` — and why

The hypothesis: replace `vvexp()` with an inline NEON polynomial exp,
which would (a) remove the function-call sync point and (b) enable a
single fused loop (kill the SoA memory traffic for u-derivatives).
Estimated ceiling: 1.7-2× scalar (vs current allFastNEON 1.33×).

**Result: 22% slower than scalar, broken on R125 dilute-gas.**

### Why slower

Standalone bench (`dev/bench/poly_exp_bench.cpp`) at the time looked
great — the SIMD poly was 1.5-2× faster than vvexp:

| N | scalar std::exp | vvexp | SIMD poly | poly vs vvexp |
|---|---:|---:|---:|---:|
| 18 | 2.47 ns | 2.22 ns | 1.02 ns | 2.18× |
| 54 | 2.37 ns | 1.61 ns | 0.99 ns | 1.63× |
| 128 | 2.35 ns | 1.50 ns | 1.00 ns | 1.50× |

But the standalone bench measures the poly call **in isolation**, where
the prefetcher and pipeline are happy.  **Inside the fused inner loop**,
the degree-11 Horner evaluation creates an 11-cycle serial FMA dep
chain.  The B-chain math right after also has dep chains.  The two dep
chains can't overlap — the polynomial exp result feeds directly into
the next step.  Compiler can't ILP across the cycle barrier.

Meanwhile, `vvexp()` internally uses ILP across batched inputs (it
amortises the polynomial chain across many parallel exps).  Its
per-exp throughput as a batched function is what the standalone bench
showed.  In the fused-loop context, each `vexp_poly_neon` call has the
full 11-FMA latency.  Result: fused is slower despite "fewer passes".

### Why broken on R125 dilute-gas (τ=3.0, δ=0.01)

12/96 equivalence failures, worst case `inf` on R125 `d3ar_dt3`.

Root cause was almost certainly numerical: at dilute gas, `1/δ = 100`,
so `(1/δ)^3 = 10⁶` magnifies any small error in `d3u_ddelta3` to a
6-order-of-magnitude error in `d3alphar_ddelta3`.  My poly exp has
~4 ULP accuracy, scalar `std::exp` has ~1 ULP.  The 3-4 ULP gap times
the per-term derivative chain times 1e6 amplification = the difference
ends up registering as inf in the rel-diff display.

Could the poly degree be bumped to 13 or 15 to fix accuracy?  Yes, at
the cost of more FMA cycles, making the perf problem worse.  No win
without a fundamentally different approach (e.g., 4-term unrolling
with parallel poly chains for ILP).

### The diff was discarded

The fused code lived in `src/Helmholtz.cpp` for ~30 minutes; build
errors, then equivalence failures, then perf regression confirmed the
approach was wrong on this platform.  `git checkout HEAD -- src/`
reverted to the working allFastNEON state.  The standalone bench
(`dev/bench/poly_exp_bench.cpp`) is kept as evidence of the per-call
poly performance for future reference — useful if someone wants to
attempt ILP-based fusion (process 4 terms via 2 parallel SIMD
accumulator chains).

## Cross-platform follow-up (not done)

The current `allFastNEON` is `__APPLE__` + `__aarch64__` gated.  Linux
x86 still uses scalar.  Real ship-worthy work for that:

- AVX2 (4-wide doubles, 256-bit): 2× more SIMD width than NEON, likely
  ~2.6× scalar speedup vs NEON's ~1.33×.
- SLEEF as cross-platform SIMD-math dependency (or hand-rolled per-arch
  polynomial exp with ILP — the right way to do what fused tried).
- Build matrix for CI to validate.

Multi-day effort.  Punted to a follow-up branch; my recommendation is
**before** doing that, profile a Linux x86 build of the current state
to confirm where time goes there — the bottleneck distribution may be
different.

## PR #3034 default-on (orthogonal change pushed during this session)

While here, also flipped the PR #3034 mechanical-bypass default from
opt-in to opt-out on `ihb/pxflash-direct-eos` branch (`b6df95681`):

```cpp
static const bool pxflash_direct_eos = []() {
    const char* env = std::getenv("PXFLASH_DIRECT_EOS");
    if (env == nullptr) return true;
    return std::string(env) != "0";
}();
```

Validated end-to-end on PXcdj_sweep — 10-17% faster across 3 runs, zero
new failures, |dT|/T = 0, |dρ|/ρ = 2.84e-6 (ULP-class).

## Current speedup ceiling (Apple ARM64)

With both changes default-on (`ihb/pxflash-direct-eos` + `ihb/helmholtz-simd-crazy`
stacked on the same build):

| Configuration | vs master |
|---|---:|
| master scalar | 1.000× |
| PR #3034 default-on bypass only | 1.21× |
| allFastNEON default-on only | 1.17× |
| **Both stacked (the new default)** | **1.41×** |

## Files in this branch

- `src/Helmholtz.cpp` — `allFastNEON` default-on dispatch on Apple
  ARM64.
- `dev/bench/poly_exp_bench.cpp` — standalone NEON SIMD polynomial
  exp benchmark, useful for future ILP-based fusion attempts.
- `dev/bench/vvexp_bench.cpp` — Accelerate vvexp vs scalar exp bench
  from earlier investigation.
- `docs/PERF-helmholtz-simd-crazy.md` — this document.

## Reproducing

```bash
# Configure + build (Release, default-on now)
cmake -B build_catch_rel -S . -DCOOLPROP_CATCH_MODULE=ON \
    -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build_catch_rel --target CatchTestRunner -j8

# Microbench: scalar vs allFastNEON (still env-toggled via the bench harness)
./build_catch_rel/CatchTestRunner "Helmholtz inner-loop: scalar vs allFastNEON"

# Equivalence — should all pass at 1e-12
./build_catch_rel/CatchTestRunner "allFastNEON equivalence to scalar all()"

# End-to-end (default-on now; opt out with COOLPROP_HELMHOLTZ_FAST=0)
PXCDJ_NT=30 PXCDJ_NP=30 PXCDJ_CSV=/tmp/on.csv \
    ./build_catch_rel/CatchTestRunner "[PXcdj_sweep]"
COOLPROP_HELMHOLTZ_FAST=0 PXCDJ_NT=30 PXCDJ_NP=30 PXCDJ_CSV=/tmp/off.csv \
    ./build_catch_rel/CatchTestRunner "[PXcdj_sweep]"
```
