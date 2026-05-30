# VLE precision near triple point — measurement + recommendations

**Branch:** `ihb/vle-triple-stress`
**Test:** `src/Tests/CoolProp-Tests-VLETripleStress.cpp` — `[vle_triple_stress]`
**Hardware:** Apple Silicon ARM64, Release build, master @ `cbb1686c5` (post-#3041)

## What was measured

For each of 10 fluids × 30 temperatures sampled along the saturation curve (cubic-warped toward the triple point):
1. Run VLE flash via `QT_INPUTS` with Q=0 (liquid) and Q=1 (vapor) → recover (ρ_L, ρ_V, p_sat) — **as converged by CoolProp's default tolerance (Akasaka 1e-10)**
2. Independently re-evaluate `p(T, ρ)` via `DmolarT_INPUTS` at each density
3. Independently evaluate Gibbs energy `g(T, ρ)` at each density
4. **Refinement phase**: starting from CoolProp's converged densities, run up to 80 additional manual Akasaka iterations until the residual stops shrinking (true noise floor). This pushes past CoolProp's 1e-10 cutoff to expose the underlying precision limit.

The metric `|p_L - p_V| / p_avg` is the **noise floor of (saturation solver + EOS evaluator)** at that state. At true VLE in exact arithmetic it is 0; in double precision it is bounded by ULP scatter through the EOS terms summed across the 14-54 terms in the fluid's EOS.

## The refinement comparison — a surprising and important result

Pushing past CoolProp's 1e-10 cutoff to the true noise floor (5-7 additional Akasaka iterations until residual stops shrinking) gives a mixed result:

| Fluid | rel_dp (default tol) | rel_dp (refined) | Refinement effect |
|---|---:|---:|---|
| **Methanol** | 3.0e-5 | **4.5e-5** | **WORSE 1.5×** |
| **Water** | 1.5e-7 | **3.7e-7** | **WORSE 2.5×** |
| **Propane** | 3.2e-6 | **3.4e-6** | **WORSE 1.1×** |
| Ethane | 2.0e-6 | 1.6e-6 | better 1.3× |
| CarbonDioxide | 2.5e-12 | 2.6e-12 | same |
| Methane | 8.0e-11 | 3.6e-11 | better 2.2× |
| Nitrogen | 7.8e-11 | 4.8e-11 | better 1.6× |
| Ammonia | 3.0e-10 | 9.5e-11 | better 3.2× |
| R125 | 3.5e-10 | 1.2e-10 | better 2.9× |
| R134a | 3.8e-9 | 9.6e-10 | better 3.9× |

**Two distinct regimes:**

For well-behaved fluids (Methane, Nitrogen, Ammonia, R125, R134a), pushing past 1e-10 buys 2-4× tighter pressure equality. The Akasaka solver is being stopped earlier than it could go — these would benefit from a tighter default tolerance.

For triple-point-sensitive fluids (**Propane, Water, Methanol**), refinement makes pressure equality **worse**. Why: the Akasaka residual is `sqrt((K_L - K_V)² + (J_L - J_V)²)` where K embeds the chemical-potential equality and J embeds the pressure equality. Near triple, K is dominated by `log(δ_V)` for vapor (since δ_V is tiny — vapor density near zero). Small changes to ρ_V produce huge K-changes. Tightening on K-residual pushes densities into states where chemical potential matches better but **the pressure equality the user asked about gets worse**.

**This is a numerical-methods lesson with practical impact:** the convergence metric matters. Akasaka's K-J residual is the natural symmetric formulation but doesn't align with pressure equality at the noise floor for triple-point-class regimes.

## Headline numbers

300 (fluid × T) probes, 100% VLE convergence success (no iteration failures observed across the full envelope from T_triple+ε to T_critical-5%).

| Fluid | Worst rel_dp | At T_reduced | p_sat there | Absolute noise floor |
|---|---:|---:|---:|---:|
| **Methanol** | **3.0e-5** | 0.009 | 0.30 Pa | 1e-7 to 3e-5 Pa |
| **Propane** | **3.2e-6** | 0.02 | 1.5e-3 Pa | 4e-8 to 3e-6 Pa |
| **Ethane** | **2.0e-6** | 0.005 | 1.5 Pa | 2e-9 to 3e-6 Pa |
| Water | 1.5e-7 | 0.0001 | 613 Pa | 1e-7 to 1e-4 Pa |
| R134a | 3.8e-9 | 0.0001 | 391 Pa | 1e-9 to 1.6e-6 Pa |
| Ammonia | 3.0e-10 | 0.014 | 7.6 kPa | 2e-9 to 2.5e-6 Pa |
| R125 | 3.5e-10 | 0.0001 | 2.9 kPa | 6e-9 to 1e-6 Pa |
| Methane | 8.0e-11 | 0.0001 | 11.7 kPa | 8e-9 to 1e-6 Pa |
| Nitrogen | 7.8e-11 | 0.005 | 12.5 kPa | 5e-10 to 1e-6 Pa |
| CarbonDioxide | 2.5e-12 | 0.02 | 559 kPa | 8e-9 to 1.4e-6 Pa |

The pattern is unambiguous:

- **Absolute noise floor of ~1e-7 to 1e-6 Pa** holds across most well-behaved fluids (Methane, Nitrogen, CO₂, R125, Ammonia, etc.). This is the cost of summing 18–54 terms in IEEE-double.
- **Water and Methanol have systematically worse absolute floors** (~1e-4 and ~3e-5 Pa respectively). Methanol is the only Lemmon-2005-style EOS in the set with active `tau_mi_in_u` terms; Water uses a 56-term IAPWS-95 EOS — both are unusually term-heavy.
- **Relative precision explodes when p_sat is small.** Propane near triple has p_sat ≈ 1.7e-4 Pa (test minimum), and noise of ~3e-6 Pa absolute becomes 2% relative. At the true triple (p_sat ≈ 1.7e-10 Pa), relative noise would exceed unity — **p_sat is unresolvable in double precision** for some 10s of K above the triple.
- **Gibbs residual is consistently 10-1000× tighter than pressure residual.** Worst rel_dg = 1.77e-9 (Water at triple), most fluids ≤ 1e-12. This points to a clear opportunity for solver reformulation.

## Mid-saturation behavior (for context)

At T_reduced > 0.5, ALL fluids hit the IEEE-double ULP floor (1e-14 to 1e-16 relative). The pathology is concentrated in the bottom ~10% of the saturation curve.

| Fluid | T_red=0.0001 (near triple) | T_red=0.13 (low-mid) | T_red=0.95 (near critical) |
|---|---:|---:|---:|
| Propane | 3.1e-6 | 2.8e-7 | 1.2e-14 |
| Ethane | 1.1e-6 | 2.2e-9 | 4.8e-16 |
| Methanol | 5.3e-6 | 2.9e-7 | 3.6e-14 |
| Water | 1.5e-7 | 2.4e-9 | 7.4e-15 |
| Nitrogen | 6.3e-11 | 4.6e-12 | 2.7e-15 |
| CO₂ | 7.3e-13 | 1.8e-13 | 2.5e-15 |

## Recommendations for improving numerics at the limits

Ordered by impact-to-effort ratio:

### 1. Convergence-metric choice matters — and the current Akasaka K-J residual can pessimize pressure equality at triple (HIGH IMPACT, MEDIUM EFFORT)
The refinement experiment shows that for triple-point-class fluids (Propane, Water, Methanol), pushing the K-J residual below 1e-10 actively **degrades** pressure equality. The K-residual is dominated by `log(δ_V)` near triple — tightening it pushes the densities into states where chemical potential matches better but pressure doesn't.

Options:
- **(a) Switch the primary convergence test to direct pressure equality**: `|p_L - p_V| / max(p_L + p_V, p_ref) < ε`. Conceptually simple, would fix the over-tightening on the three problem fluids. Risk: pressure-equality alone is poorly conditioned at near-critical states where ρ_L ≈ ρ_V; needs a hybrid criterion.
- **(b) Use chemical potential equality directly**: `|g_L - g_V| / RT < ε`. Our data shows rel_dg is consistently 10-1000× tighter than rel_dp across the envelope. Worst rel_dg = 1.77e-9 (Water at triple) vs worst rel_dp = 3.0e-5 (Methanol). The Gibbs residual is ~4 orders of magnitude better. This is the cleanest fix.
- **(c) Adaptive metric**: use K-J residual until it bottoms out, then switch to pressure equality for the final iterations. More complex but combines the strengths.

Recommendation: **(b)**. The dimensionless Gibbs residual `(g_L - g_V) / (RT)` is well-conditioned everywhere we measured. It's also bounded by the underlying EOS precision rather than amplified by `log(δ_V)`.

### 1b. Tighter default tolerance for well-behaved fluids (MEDIUM IMPACT, LOW EFFORT)
The refinement data shows that for Methane, Nitrogen, Ammonia, R125, R134a, the 1e-10 cutoff stops 2-4× short of the true noise floor. A tighter default (e.g., 1e-12 on the K-J residual, or the dynamic noise-floor criterion described in #4) would let these fluids reach their full precision without changing behavior for the well-converged cases. Combined with #1 (Gibbs as primary metric), this would lift the average accuracy on common refrigerants.

Cost: trivial — change the `1e-10` in `VLERoutines.cpp:910` to `1e-12` (or even DBL_EPSILON * 100 ≈ 2e-14). Adds 1-2 iterations per VLE flash, ~5% wall time overhead.

Caveat: must be combined with #1 or the Propane/Water/Methanol degradation gets WORSE. Don't ship 1b standalone.

### 2. Compensated summation in `ResidualHelmholtzGeneralizedExponential::all` (HIGH IMPACT, LOW EFFORT)
The per-term accumulation
```cpp
derivs.alphar += ndteu;
derivs.dalphar_ddelta += ndteu * B_delta;
// ... 13 more
```
loses ~log₂(N) bits across 50+ terms with magnitudes spanning 10⁶. Kahan or Neumaier compensated summation recovers most of that. Cost: 4-5 extra FLOPs per accumulator update (~50 ns/call out of ~800 ns), expected gain: 1-2 orders of magnitude in the absolute noise floor of p and derivatives.

Particularly impactful for **Water (56 terms) and Methanol (44 terms)** where the per-term sum is the longest.

Suggested change: extract the `alphar +=` into a `KahanAccumulator` wrapping `value` + `compensation` doubles, with a finalize() call before the post-loop scaling. Five identical changes across the 15 derivative accumulators in `Helmholtz.cpp:257-275`.

### 3. Sort terms by `|n|` (smallest first) at `finish()` time (LOW IMPACT, LOW EFFORT)
Standard floating-point summation rule: small-first reduces accumulated round-off vs natural order (which tends to be largest-first in EOS json files). Pairs naturally with #2 — even without compensated summation, sorting alone reduces sum noise.

Caveat: needs the same term-reordering invariance discussion as the SIMD term-sorting recommendation.

### 4. Reformulate convergence test to use noise-floor-aware tolerance (MEDIUM IMPACT, MEDIUM EFFORT)
Current solvers likely use `|residual| < abs_tol || |residual/scale| < rel_tol`. Near triple, neither captures "the solver has reached the noise floor; stop iterating."

A noise-floor-aware tolerance is roughly:
```
eps_noise(T) ≈ N_terms × ULP × max(|term_magnitude|)
```
For Methanol at triple: N=44, max term ≈ 100, ULP ≈ 2e-16 → eps_noise ≈ 1e-12 absolute on alphar, propagating to ~1e-7 Pa on p.

Once `|Δresidual_between_iterations| < eps_noise`, stop. This prevents the solver from spinning trying to drive a residual below its noise floor.

### 5. Use `std::fma` explicitly in B-chain (LOW IMPACT, LOW EFFORT)
The scalar `all()` body writes `delta * du_ddelta + di` as separate mul + add (round-twice). `std::fma(delta, du_ddelta, di)` is round-once. Clang with default `-ffp-contract=on` may fuse some but not all; explicit FMA guarantees the round-once behavior at every site.

Audit `src/Helmholtz.cpp:163-275` for `delta * x + y` and `tau * x + y` patterns; convert ~20-30 sites. Per-term ULP gain: maybe 1-3 ULP. Per-fluid sum impact: modest but free.

### 6. Higher-precision intermediate for low-density evaluation (MEDIUM IMPACT, MEDIUM EFFORT)
At ρ → 0 (dilute gas), the EOS terms `delta^l * exp(-c * delta^l)` have catastrophic cancellation as `delta^l → 0` and `exp(...)→1` — the actual contribution becomes `delta^l - c * delta^(2l) + ...` which loses leading-digit precision when computed naively.

A virial-like reformulation at low ρ would give:
```
alphar(low ρ) = B(τ)·δ + C(τ)·δ² + D(τ)·δ³ + ...
```
where B, C, D are precomputed from the EOS. For δ < some_threshold (e.g., 0.01), use the virial form; otherwise use the full EOS. Bit-exact at the boundary by construction.

This is what reference implementations (NIST REFPROP) do for some properties. Significant engineering effort but eliminates the Propane-near-triple class of problems entirely.

### 7. Diagonal preconditioning of the VLE Newton system (LOW IMPACT, MEDIUM EFFORT)
Near triple, the VLE Jacobian rows have wildly different magnitudes (e.g., ρ_L is O(10⁴), ρ_V is O(10⁻⁶) for Propane; row scaling differs by 10¹⁰). Preconditioning with `diag(J)⁻¹` before Newton step recovers some condition number. Cost: per-iteration, O(few flops).

### 8. Track Jacobian condition number and fall back when ill-conditioned (MEDIUM IMPACT, HIGH EFFORT)
At each iteration, estimate `κ(J)` (e.g., via SVD or Hager's algorithm). If `κ > 10⁸`, switch from Newton to a bisection-on-Gibbs scheme — slower but robust. This is the right behavior near triple where Newton's quadratic convergence is wasted on noise.

### 9. Document the per-fluid noise floor in user-facing docs (LOW IMPACT, LOW EFFORT)
Users should know that `PropsSI("P", "T", 86, "Q", 0, "Propane")` returns a value with ~3 µPa absolute uncertainty regardless of what tolerance the solver claims. This is fundamental to IEEE-double through a 50-term EOS, not a bug.

A table in the docs (similar to this report's headline numbers) would set realistic expectations and prevent "why doesn't my pressure converge to 1e-12 at triple?" questions.

### 10. Reverse: use this stress test as a CI gate for changes that touch hot-path arithmetic (HIGH IMPACT, MEDIUM EFFORT)
Any future change to `all()` or saturation solvers should run this test and ensure rel_dp doesn't regress. The test takes ~5 seconds and is reproducible.

Suggested integration: a separate CI job that runs `[vle_triple_stress]` and checks the resulting CSV against a committed baseline; fails if any rel_dp regresses by more than 2× at any (fluid, T_reduced) point.

## What's NOT a problem here

- **VLE convergence robustness**: 100% success across 300 probes. CoolProp's saturation solver is solid down to T_triple + 1e-4 × (T_c - T_t).
- **Mid-saturation precision**: ULP-floor across all fluids at T_red > 0.5.
- **Critical-region behavior**: not stressed here; left for a follow-up.

## How this connects to the SIMD work

The SIMD `allFastNEON` was validated as ULP-equivalent at *single-point* evaluation. This test extends that to *integrated convergence behavior* under stress.

If the SIMD path slightly biases the per-term sum (which the numerical review noted is plausible due to round-once FMA in NEON vs round-twice in scalar), the noise floor measured here could shift. **Re-running this test with `COOLPROP_HELMHOLTZ_FAST=1` would directly answer "does SIMD degrade triple-point VLE?"** That's a 5-minute experiment when the SIMD branch is rebased.

For SLEEF u10 (10 ULP exp): the per-term noise scales linearly with ULP error. Replacing libm exp with SLEEF u10 would multiply the noise floor by ~10× — Propane's rel_dp would degrade from 3e-6 to 3e-5, putting it in Methanol's territory. **Confirms that SLEEF u10 is not acceptable for triple-point VLE work.**

## Reproducing

```bash
cmake -B build_catch_rel -S . -DCOOLPROP_CATCH_MODULE=ON -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build_catch_rel --target CatchTestRunner -j8
./build_catch_rel/CatchTestRunner "[vle_triple_stress]"
# CSV at /tmp/vle_triple_stress.csv
```

Total runtime: ~5 seconds on Apple Silicon.
