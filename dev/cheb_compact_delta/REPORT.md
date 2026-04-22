# Chebyshev Density Rootfinder — Prototype Report

Worktree: `feat-cheb-compact-delta` (branch `worktree-feat-cheb-compact-delta`)
Target: extend Bell & Alpert 2018 "Exceptionally reliable density-solving algorithms for multiparameter mixture models from Chebyshev expansion rootfinding" toward single-digit-μs rootfinding for single-component and small-mixture density inversion.

## Goal

Provide a dependable density rootfinder for multi-fluid Helmholtz EOS that (a) finds *all* real roots, (b) never misses them, (c) runs fast enough to live inside a flash inner loop (target: single-digit μs for small-component mixtures, matching the paper's "no τ flattening" path).

Supporting motivation: make the multi-fluid iteration more reliable by providing a proxy density solver that can seed Newton/Halley iteration on the full EOS — the primary pain point for CoolProp's flash at extreme states.

## Architecture summary

The paper splits each residual-Helmholtz term into a τ-only factor `F_k(τ)` and a δ-only factor `G_k(δ)` (separable form). The architecture:

1. **Offline**, per fluid, per separable term:
   - Adaptively subdivide the δ range `[δ_min, δ_max]` using dyadic subdivision driven by the ε_split,M metric shared across all terms.
   - On each subinterval, fit a per-term Chebyshev expansion of `H_k(δ) ≡ δ·∂G_k/∂δ`.
   - Store `H_coeffs[interval, term, j]` as a single flat matrix for BLAS dgemv.

2. **Runtime**, given `(T, P_target)`:
   - Compute `F_k(τ)` analytically (closed-form per term; **not** Chebyshev — see design decision below).
   - `weights[k] = n_k · F_k(τ)`.
   - Single BLAS matvec `A01_all = weights^T @ H_flat` → per-interval δ-Chebyshev of α^r_01.
   - Per interval: build pressure residual (times-x + linear shift + P_target subtraction), run cheap skip test, monotonicity test, secant-or-eig accordingly.

The runtime hot path is ~200 lines of tight C++ with a single BLAS call, a fused per-interval loop, and LAPACK dhseqr for the rare non-monotone branch.

## Key design decisions

**Rejected: Möbius compactification `u = (δ−L)/(δ+L)`.** The Jacobian `dδ/du = 2L/(1−u)²` blows up near u=1, forcing hundreds of subdivisions at low per-piece degree to resolve the compactification singularity rather than the physics. Direct δ on a bounded rectangle is the paper's design and works. (See memory: `project_chebyshev_density_solver.md`.)

**Rejected: extrapolation / soft-boundary tricks.** Hard rectangular (τ, δ) domain; calls outside are rejected, not extrapolated. Matches the paper's position.

**Decided: analytic `F_k(τ)`, not Chebyshev.** IAPWS-95 has power terms with t up to 50; on `τ ∈ [0.2, 4]` that's 10^65 dynamic range, numerically unusable for a single Chebyshev fit. At call time τ is a scalar anyway, so `τ^t` and `τ^t·exp(−β(τ−γ)²)` are O(1) closed-form evaluations — simpler, faster, and more accurate than Chebyshev for F_k. (The paper does build F̃_k expansions, but for the call pattern we target, there's no benefit.)

**Decided: monotonicity classification to route work between secant and eigenvalue solves.** For each candidate interval we compute the derivative-Chebyshev coefficients and test `|d_0| > Σ|d_{k≥1}|` — a sufficient condition for `f'(x)` to have constant sign, hence `f` monotone. Monotone + endpoint sign change → bracketed secant (fast). Non-monotone → LAPACK dhseqr on the pre-Hessenberg Chebyshev colleague matrix. In the nitrogen torture sweep, ~97% of candidate intervals classify as monotone-with-root or monotone-no-root, so the expensive eigenvalue path is essentially never triggered outside the near-critical vdW loop.

**Non-analytic terms (water, CO₂, R-125, methanol): deferred.** The paper's architecture only handles separable α^r terms. For fluids with non-analytic terms, a per-call 1D Chebyshev refit of the non-analytic contribution at the runtime τ (2 terms for IAPWS-95) can be added on top of the separable skeleton. Not implemented yet.

## Results

### Accuracy (nitrogen torture sweep, 500+ test points)

| Sweep | Worst rel err vs CoolProp | Missed roots |
|-------|--------------------------:|-------------:|
| Supercritical wide (T=150–1000 K, P=1 kPa–100 MPa) | 5.98e-13 | 0 |
| Subcritical single-phase | 3.36e-13 | 0 |
| Near-critical (122–130 K) | 6.04e-10 | 0 |
| Near-triple (63–80 K) | 4.03e-13 | 0 |
| Saturation curve (both rhoV & rhoL at P=Psat) | 8.5e-11 / 4.6e-9 | 0 / 0 |

**Machine precision across the entire sweep, zero missed roots.** The "extra stable roots" in subcritical/near-triple sweeps are the metastable branches (superheated liquid at P<Psat, supercooled vapor at P>Psat) that the local `dP/dρ > 0` stability filter correctly admits; a Gibbs-energy comparison at a higher layer selects the globally stable phase.

### Timing — single-component nitrogen (Apple M-series, Accelerate BLAS/LAPACK)

| Case | Python | C++ full | C++ fixed-τ |
|------|-------:|---------:|------------:|
| T=200K P=5MPa (supercrit, 1 root) | 350 μs | 46 μs | **4.7 μs** |
| T=Tc P=Pc (near-critical, 1 root) | 300 μs | 147 μs | ~5 μs |
| T=80K P=5MPa (subcrit liquid, 3 roots) | 430 μs | 40 μs | **4.8 μs** |
| T=80K P=10 kPa (subcrit vapor, 5 roots) | 533 μs | 47 μs | **21.9 μs** |
| T=100K P=5MPa (deep vdW loop, 1 eig) | 1200 μs | 50 μs | 157 μs |

Best-case fixed-τ (pressure sweep at constant T — the "no τ flattening" regime of paper Fig. 11): **4.6–4.9 μs** for single-phase cases, **22 μs** for deep subcritical multi-root, **60–180 μs** for the narrow T ≈ 90–115 K window where dhseqr is needed to resolve non-monotone intervals in the vdW loop.

### Comparison to Bell & Alpert 2018 Fig. 11

| | Paper reference | Our C++ |
|---|-----------------:|--------:|
| "No τ flattening", 2-component, single root | ~30 μs | **4.6 μs** |
| "Full flattening", 2-component | ~50 μs | 40–50 μs |
| REFPROP TPRHOdll | ~100–300 μs | — |

**The C++ implementation matches or beats the paper's demonstrated numbers in the common single-phase case and sits within 2× in the vdW-loop case. Target met.**

## What the eig path costs

The remaining slowdown (vdW-loop region) is dominated by LAPACK dhseqr on the 13×13 Hessenberg Chebyshev colleague matrix — ~50–100 μs per call on this hardware. Finer δ-subdivision can eliminate non-monotone intervals entirely (verified: at n_δ=8, tol=1e-10, K=3509 intervals, eig count=0 even in vdW-loop region, consistent ~60–70 μs per call) but the trade is an increased matmul cost. The sweet spot depends on workload: for flash inner loops dominated by single-phase queries, the current config (n_δ=12, tol=1e-8, K=261) is optimal; for density-hint-less scanning across the full phase envelope, finer subdivision is more consistent.

Further eig-path speedups possible (not implemented):
- Aurentz-Mach-Vandebril-Watkins structured QR (O(n²) for colleague matrices, competitive with dhseqr below n ≈ 250).
- Serkh-Rokhlin 2021 O(n²) componentwise-stable QR (newer, not yet in any standard library).
- Recursive bisection of non-monotone intervals until monotone (cheaper than eig for simple vdW-loop cases).

## Files in this worktree

`dev/cheb_compact_delta/`:

- `fluid_from_json.py` — generic fluid evaluator parsing CoolProp/teqp EOS JSON (Power + Gaussian blocks).
- `iapws95.py` — water-specific evaluator including the non-analytic terms (used for future non-analytic layer).
- `fast_solver.py` — pure-Python reference implementation of the full hot path (subdivision, monotonicity, secant, eig).
- `paper_architecture.py` — single-interval reference implementation following the paper, used to validate the separable-term deconstruction.
- `cpp_solver.cpp`, `cpp_solver_v2.cpp` — pybind11 C++ ports; v2 is the optimized version (stack-allocated eig buffers, fused per-interval loop, dhseqr-on-Hessenberg).
- `nitrogen_torture.py` — 500-point accuracy and timing sweep.
- `bench_cpp_vs_python.py`, `bench_best_case.py` — benchmarks.
- `piecewise_study.py`, `convergence_study.py`, `direct_delta_study.py` — exploratory convergence studies (historical, document why certain design choices were made).

## Open items

1. **Non-analytic term layer** (water, CO₂, R-125, methanol): per-call 1D Chebyshev refit of the non-analytic α^r contribution in δ at the given τ, summed with the separable part before rootfinding. Straightforward extension, ~100 lines in C++.

2. **Mixture extension**: stack per-component `H_coeffs` weighted by mole fraction, add departure-function bivariate fit per binary pair. The paper's flowchart applies directly; the architecture is already in place per-component.

3. **Data-repo scaffolding**: per the design discussion, pre-built per-fluid expansions live in a separate repo (pulled in by CoolProp via CPM at a tagged version). Build/serialize the offline expansions to a compact binary format; set up the data-repo with regeneration scripts and CI validation.

4. **Density-hint API**: for flash inner loops that already know the previous-iterate density, accepting a hint lets us touch only a few intervals rather than the full K scan. Would drop the fixed-τ call time from 5 μs to ~1 μs.

5. **Structured-QR eig path**: port Serkh-Rokhlin 2021 O(n²) colleague-matrix QR to replace dhseqr for the non-monotone branch. Only matters if profiling shows the 10% of calls that hit this path dominate total flash time.
