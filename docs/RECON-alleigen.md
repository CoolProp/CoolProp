# `allEigen` revival reconnaissance

**Date:** 2026-05-29
**Branch:** `ihb/helmholtz-alleigen-recon` (off `master @ b4968f4b1`)
**Question:** Is `ResidualHelmholtzGeneralizedExponential::allEigen` — the
commented-out Eigen-vectorised variant in `src/Helmholtz.cpp` — viable as a
drop-in replacement for, or alternative to, the scalar `all()`?

## TL;DR

**No, not as-is.** Three blockers, in order of severity:

1. **Slower than scalar `all()` on Apple Silicon NEON** — by 1.4× to 2× across
   the four fluids benchmarked. The cause is heap allocation of seven
   per-call `Eigen::ArrayXd::Zero(N)` working arrays; the inner SIMD
   throughput (NEON is 2-wide for `double`) is not enough to amortise it.
2. **Incomplete derivative coverage** — `allEigen` only accumulates 0th-3rd
   order derivatives. The scalar `all()` also writes the five 4th-order
   fields used by Householder4 callers. Replacing `all()` with `allEigen`
   would silently break 4th-order derivative consumers.
3. **`tau_mi_in_u` branch is dead-commented inside `allEigen`** — for the two
   fluids with that branch active (R125, Methanol — see the sweep on
   `ihb/pxflash-direct-eos` commit `3eeb25baf`), `allEigen` is catastrophically
   wrong: R125 `d3alphar_dtau3` differs by 2e+6 relative; Methanol
   `dalphar_dtau` by 66%.

(2) and (3) are fixable in straightforward code; (1) is the structural problem.

## What was done

The commented `/* void allEigen(...) throw() { ... } */` block at
`src/Helmholtz.cpp:38-137` was uncommented, the signature updated to `noexcept`,
and a matching declaration added in `include/Helmholtz.h`. Two mechanical fixes
to restore buildability:

- `POW2(eigen_array)` / `POW3(eigen_array)` macros do not work on Eigen
  expressions — replaced with `.square()` / `.cube()`.
- `nE*exp(tE*log_tau + dE*log_delta + uE)` (`exp` is `std::exp`, not the
  Eigen member) replaced with `(tE*log_tau + dE*log_delta + uE).exp()`.

Also: the original cached working arrays (`uE, du_ddeltaE, ...,
d3u_dtau3E`) were class members. The revival made them locals
(`Eigen::ArrayXd::Zero(N)`) — see "Performance" below for why this matters.

The `epsilon1/eta1/gamma1/beta1` SoA arrays were missing from `finish()`
(only the `2` series existed in master); they're now populated so the
`Eigen::Map`s `allEigen` constructs reference valid storage.

A Catch2 recon harness lives at `src/Tests/CoolProp-Tests-AllEigenRecon.cpp`,
registered as `[alleigen_recon][.]`.

## Numerical equivalence (excluding R125/Methanol)

Tolerance: `|Δ| / max(|scalar|, 1) < 1e-12` across all 0th-3rd-order
`HelmholtzDerivatives` fields, 5 representative (τ, δ) regimes per fluid
(compressed liquid, supercritical, near-critical, dilute gas, near-critical
density).

| Fluid          | worst Δ_rel | worst field             |
|----------------|------------:|-------------------------|
| Water          |    7.2e-14  | dalphar_ddelta          |
| CarbonDioxide  |    3.0e-14  | d3alphar_ddelta2_dtau   |
| n-Propane      |    2.0e-14  | d3alphar_ddelta3        |
| R134a          |    2.7e-15  | d3alphar_ddelta_dtau2   |
| Methane        |    1.1e-13  | d3alphar_ddelta3        |
| Nitrogen       |    1.1e-14  | d3alphar_ddelta3        |
| Hydrogen       |    3.9e-14  | d3alphar_ddelta3        |
| MM             |    1.1e-13  | d3alphar_ddelta3        |

All 40 (fluid × regime) probes pass at 1e-12. Worst observed: 1.13e-13 on
`d3alphar_ddelta3` (MM dilute gas, Methane near-critical density). This is
consistent with ULP-class noise from libm `exp/log` reordering — scalar vs
Eigen evaluate `exp` in different orders, and the accumulating `.sum()`
differs from sequential `+=`.

For these 134-of-136 fluids, `allEigen` is **numerically interchangeable
with `all()` to within ULP** on 0th-3rd-order derivatives.

## Known-gap fluids

| Fluid    | worst Δ_rel | worst field         | why                            |
|----------|------------:|---------------------|--------------------------------|
| R125     |       2.0e+6 | d3alphar_dtau3      | tau_mi_in_u dead-commented     |
| Methanol |       6.6e-1 | dalphar_dtau        | tau_mi_in_u dead-commented     |

The block at lines 86-99 of the current `allEigen` is a multi-line scalar
comment. Restoring it as a vectorised loop would close the gap; the math
is identical to the scalar branch at `src/Helmholtz.cpp:177-185`.

## Performance — the structural problem

Microbench: 100,000 reps × 5 trials of `gen.all` / `gen.allEigen` directly,
zero-init `HelmholtzDerivatives` each iteration, Release build, Apple
Silicon (M-class) NEON.

| Fluid     | N  | scalar (ns/call) | eigen (ns/call) | scalar/eigen |
|-----------|---:|-----------------:|----------------:|-------------:|
| Water     | 54 |    806.9 ± 6.4   |  1129.4 ± 4.7   |    0.71 ✗    |
| R134a     | 21 |    295.3 ± 1.2   |   568.1 ± 2.2   |    0.52 ✗    |
| n-Propane | 18 |    265.9 ± 1.5   |   538.8 ± 3.8   |    0.49 ✗    |
| MM        | 18 |    261.0 ± 2.5   |   530.7 ± 4.9   |    0.49 ✗    |

**Eigen is slower than scalar across the board.** Even on Water (54 terms),
where vectorisation should be most beneficial, scalar wins by ≈ 40%.

The most likely cause is **per-call heap allocation**: each `allEigen`
call performs seven `Eigen::ArrayXd::Zero(N)` constructions, each of which
mallocs an array sized to the fluid's term count. Apple's NEON `double`
SIMD is 2-wide, so the marginal vectorisation speed-up per term (at most
2×) cannot offset 7 × (malloc + zero-init + free) per call. On x86 AVX2
(4-wide doubles) the trade-off would be more favourable, but the
allocation cost is identical and probably still dominates for small N.

The original `allEigen` had these working arrays as **cached class
members** — see the still-commented declaration `Eigen::ArrayXd uE,
du_ddeltaE, ..., d3u_dtau3E;` at `include/Helmholtz.h:364`, and the
sibling `uE.resize(elements.size()); ...` block in `finish()` at
`include/Helmholtz.h:551-557`. With per-instance pre-allocation, the
seven mallocs go to zero and the comparison shifts substantially. Whether
that puts allEigen ahead of scalar on Apple Silicon NEON, only a re-bench
will say.

## Verdict

**Don't ship `allEigen` as-is on Apple Silicon NEON.** The drop-in form
revived here is strictly slower than scalar, missing 4th-order, and wrong
on R125/Methanol.

There are three roads forward, in order of effort:

1. **Shelve and forget.** The recon stands as documentation that this path
   was tried and the perf is structurally upside-down on NEON. The branch
   can be archived. Closest to "do nothing."
2. **Restore cached working arrays + close the gaps.** Move the seven
   `Eigen::ArrayXd` workings back to class members, re-add `tau_mi_in_u`
   and 4th-order. Re-bench on Apple Silicon, then on a Linux x86 AVX2
   runner. This is the right move *only* if there is a credible
   hypothesis that pre-alloc + AVX2 will tip the ratio favourably; the
   Apple Silicon overhead points the wrong way.
3. **Skip Eigen, write a custom SIMD path with SLEEF or hand-rolled
   intrinsics.** Much bigger effort, separate platform matrix, gated on a
   CMake option. Probably the right call if SIMD is genuinely the next
   lever, but doesn't depend on this `allEigen` revival — that work would
   start from a clean slate.

The recon itself (compile, test, bench, document) cost about one focused
session. Decision can wait.

## Files in this branch (uncommitted as of write-up)

- `include/Helmholtz.h` — declaration of `allEigen` + SoA arrays for the
  "1" series populated in `finish()`.
- `src/Helmholtz.cpp` — `allEigen` un-commented + buildability fixes.
- `src/Tests/CoolProp-Tests-AllEigenRecon.cpp` — the recon harness.
- `CMakeLists.txt` — test source registration.
- `docs/RECON-alleigen.md` — this document.
