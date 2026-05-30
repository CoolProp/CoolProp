# Residual-Helmholtz SIMD investigation — full report

**Date:** 2026-05-30 (started)
**Branches:** `ihb/helmholtz-simd-investigate`, `ihb/stacked-measure`
**Hardware:** Apple Silicon (ARM64, NEON 2-wide doubles)
**Time budget:** 8 hours; this report covers what was tried, what worked, and
what's left.

## TL;DR — definitive headline

| Configuration | Sum μs (30×30 PXcdj_sweep, 362 004 successful rows) | vs master |
|---|---:|---:|
| master scalar | 9 189 045 | 1.000× |
| PR #3034 mechanical bypass only | 7 588 056 | **1.211× (21.1% faster)** |
| `allFastNEON` only (this branch) | 7 871 420 | **1.167× (16.7% faster)** |
| **BOTH stacked** | **6 499 808** | **1.414× (41.4% faster)** |

- **NEON vs master:** zero drift on T AND ρ across all 362 004 successful
  rows.  Bit-exact at the end-to-end flash level.
- **BOTH vs master:** zero drift on T, |Δρ|/ρ = 2.84 × 10⁻⁶ on ρ (inherited
  from the PR #3034 bypass path; NEON contributes none).
- **G3:** identical row passes across all 4 configurations.
- Optimizations compose **multiplicatively** (1.211 × 1.167 ≈ 1.413 ≈ 1.414
  observed).

## Headline microbench (per-call `all()` wall time on Apple Silicon)

24-fluid cross-section, 200 k reps × 7 trials, Release build.  Top winners
and one representative small fluid:

| Fluid | N | scalar (ns/call) | NEON (ns/call) | speedup |
|---|---:|---:|---:|---:|
| Methanol | 44 | 766 | 525 | 1.46× |
| Oxygen | 32 | 455 | 309 | 1.47× |
| Nitrogen | 36 | 564 | 396 | 1.43× |
| R152a | 40 | 528 | 372 | 1.42× |
| Water | 54 | 814 | 597 | 1.36× |
| n-Propane | 18 | 261 | 200 | 1.31× |
| Toluene | 12 | 172 | 136 | 1.26× |
| **Aggregate (24 fluids)** | | 9 200 | 6 900 | **1.33×** |

Every fluid is faster.  Small fluids (N≤14) benefit less because per-call
overhead (vvexp setup, NEON prelude) is amortised less.

## What this means structurally

Per-term cost on Apple Silicon ARM64 with scalar `all()`:
- ~2 `exp()` calls (`exp(l·log_δ)` for `delta_li_in_u` branch, main `exp(t·log_τ + d·log_δ + u)`)
- ~50-70 multiply-adds in the B-chain (B_δ, B_δ2, B_δ3, B_δ4, B_τ, B_τ2, B_τ3, B_τ4)
- 15 derivative accumulations (`derivs.X += ndteu * Bxxx`)

Scalar `exp()` on Apple Silicon is fast (~2.4 ns) — so `exp` is only ~30%
of per-term cost.  The B-chain is the dominant 60% of cost, and that's
where NEON 2-wide SIMD with FMAs wins big.

## Architecture: `allFastNEON`

Three-phase implementation, ~280 lines in `src/Helmholtz.cpp`:

**Phase A — scalar u-construction with batched delta_li exp().** Walk
elements once: write u_partial / 8 derivatives to stack SoA arrays
(`alignas(16) double u_partial[128]` and friends).  Collect delta_li exp
arguments into a separate batched array; `vvexp()` them once.  Add
batched results into u_partial / derivatives.

**Phase B — `vvexp()` the main exp arguments.** One more `vvexp()` call
for `exp(t·log_τ + d·log_δ + u_partial)` across all terms.

**Phase C — NEON 2-wide SIMD B-chain + 15 accumulators.** Process two
terms simultaneously with `float64x2_t` and `vfmaq_f64` FMAs.  All 15
derivative accumulators are also SIMD lanes; horizontal-reduce with
`vaddvq_f64` at the end.  Scalar tail handles odd N.

Local SoA copy of `n`/`d`/`t` from `elements[i]` — the conditional check
falls back to a stack copy when the class SoA isn't populated (mixture
excess GenExp; see "Bug caught" below).  Common-case (pure fluid) reads
from the class SoA directly.

**Constant-rate hot loop** — inlined arithmetic with FMAs, no function
calls except the two vvexp() invocations and the explicit `std::pow()`
in the rare `tau_mi_in_u` branch.

## Wiring: env-toggled dispatch

The live `all()` opens with:

```cpp
static const bool use_fast = (std::getenv("COOLPROP_HELMHOLTZ_FAST") != nullptr);
if (use_fast) { allFastNEON(tau, delta, derivs); return; }
// ... scalar body ...
```

`static const` evaluates `getenv()` once at first call.  Subsequent calls
read a single bool; cost is < 1 ns/call, lost in the noise.

To enable: `COOLPROP_HELMHOLTZ_FAST=1` in env.  Default off — zero risk
to existing builds.

## Bug caught: mixture excess GenExp doesn't call `finish()`

The original NEON implementation read directly from the class SoA
arrays:

```cpp
const float64x2_t v_ni = vld1q_f64(&n[i]);   // assumed n.size() >= N
```

On the R32 & R125 mixture's departure (excess) GenExp term — which is
constructed and added to `elements[]` but never has `finish()` called —
`n.size() == 0` and the SIMD load segfaulted.

Fix: conditional `n.size() == N ? n.data() : local_copy` at the top of
`allFastNEON`.  In the rare mixture-excess case we pay ~3*N stores to
build a stack-resident SoA copy; common case (pure fluid) pays only the
branch.

Caught by running `[mixture]` tests under the NEON toggle.  Would have
been a runtime crash in any production user.

## What was tried + verdicts

### A. `allEigen` revival (from prior recon `ihb/helmholtz-alleigen-recon`)

**Verdict: shelved.**  The commented-out Eigen-vectorized variant in
`src/Helmholtz.cpp:38-137` revived cleanly with minor API drift fixes
(`POW2(eigen_array)` → `.square()`, free-function `exp` →
`.exp()` member).  But:

- Slower than scalar on Apple Silicon NEON (0.5×-0.7× on the bench)
  due to per-call `Eigen::ArrayXd::Zero(N)` allocations.
- Missing 4th-order derivatives.
- `tau_mi_in_u` branch dead-commented (catastrophic on R125 & Methanol).

Could be revived for AVX2 platforms with pre-allocated working arrays,
but the Apple-Silicon overhead points the wrong way.  See
`docs/RECON-alleigen.md` from the prior recon branch.

### B. `allFastVDSP` — vvexp() only, no SIMD B-chain

**Verdict: superseded by `allFastNEON`.**  First attempt at SIMD: just
batch the per-term exp() calls into `vvexp()`, leaving the B-chain
scalar.  Stable measured speedup: ~3.5% aggregate (1.03× — 1.10× per
fluid).

The bench made the obvious conclusion: exp() is ~30% of per-term cost on
Apple Silicon (scalar exp is already very fast at ~2.4 ns).  Batching
exp recovers ~30% × 30% ≈ 9% of work, but the 3-pass refactor's memory
traffic eats most of that.  Net: ~3.5%.  Underwhelming.

`allFastVDSP` is kept in the codebase as a stepping-stone for the NEON
variant and for measurement (the env toggle can route to either).

### C. `allFastNEON` — vvexp + NEON 2-wide SIMD B-chain ✅

**Verdict: shipped behind env toggle.**  The big win.  See "Architecture"
above.  ~33% aggregate per-call, ~17% end-to-end.

### D. Custom polynomial exp (Estrin / Schraudolph)

**Not tried, deferred.**  Scalar `std::exp` on Apple Silicon at 2.4 ns
is already very fast — a polynomial approximation could only save ~1-1.5
ns per call.  With ~108 exp() calls per Water-class `all()` invocation,
that's ~150 ns savings on top of NEON's ~200 ns wins.  Material on small
fluids but accuracy risk (ULP-2 to ULP-5) needs careful validation
against CoolProp's test suite.  Punt.

### E. SLEEF SIMD library (cross-platform)

**Not tried, future work.**  SLEEF would give equivalent speedups on
Linux x86 (AVX2 4-wide doubles, AVX-512 8-wide).  Adding the dependency
is one CMake edit + one external_project.  Probably 2-3 hour follow-up;
gains depend on platform width.

## What's next

In rough order of value:

1. **Decide whether to merge** `allFastNEON` into PR #3034 or ship as a
   separate PR.  My recommendation: separate PR.  The PR #3034 mechanical
   bypass is a coherent story ("eliminate redundant `all()` calls"); this
   work is a coherent separate story ("make each `all()` call faster").
   They stack multiplicatively, but reviewers benefit from clear
   boundaries.
2. **Stop being env-toggled.**  Once tested in CI on real hardware,
   make `allFastNEON` the default `all()` body on Apple ARM64 (drop
   the env check, just unconditionally branch on the `#ifdef`).  Other
   platforms remain on scalar.
3. **Cross-platform via SLEEF.**  Linux x86 with AVX2 should see
   ~2× the speedup we got on NEON.  Add SLEEF, port the SIMD path to
   `__m256d` intrinsics or use SLEEF's helpers.
4. **Optional: custom polynomial exp** for very small fluids where vvexp
   call overhead dominates.

## Files in this branch

- `src/Helmholtz.cpp` — `allFastVDSP` (~210 LOC) and `allFastNEON` (~280
  LOC) implementations, env-toggled dispatch at start of `all()`.
- `include/Helmholtz.h` — `allFastVDSP` and `allFastNEON` declarations.
- `src/Tests/CoolProp-Tests-HelmholtzInnerBench.cpp` — equivalence tests
  + microbench harness for both fast paths.
- `src/Tests/CoolProp-Tests-PXcdj.cpp` — borrowed from PR #3034 for the
  end-to-end sweep measurement.
- `CMakeLists.txt` — registers new test sources, links Accelerate
  framework on macOS.
- `docs/PERF-helmholtz-simd.md` — this document.

## Reproducing the measurements

```bash
# Configure + build (Release)
cmake -B build_catch_rel -S . -DCOOLPROP_CATCH_MODULE=ON \
    -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build_catch_rel --target CatchTestRunner -j8

# Microbench (24 fluids, scalar baseline + NEON comparison)
./build_catch_rel/CatchTestRunner "[helmholtz_inner_bench]"
./build_catch_rel/CatchTestRunner "[helmholtz_inner_bench_neon]"

# Equivalence (NEON vs scalar, every fluid × 4 regimes × 10 fields)
./build_catch_rel/CatchTestRunner "allFastNEON equivalence to scalar all()"

# End-to-end PXcdj sweep — toggle via COOLPROP_HELMHOLTZ_FAST=1
PXCDJ_NT=30 PXCDJ_NP=30 PXCDJ_CSV=/tmp/off.csv \
    ./build_catch_rel/CatchTestRunner "[PXcdj_sweep]"
COOLPROP_HELMHOLTZ_FAST=1 PXCDJ_NT=30 PXCDJ_NP=30 PXCDJ_CSV=/tmp/on.csv \
    ./build_catch_rel/CatchTestRunner "[PXcdj_sweep]"
```

The CSVs have `success`, `T_solved`, `rho_solved`, `microseconds`
columns.  See the analysis Python in this commit's commit message for the
G3/G5/drift extraction.
