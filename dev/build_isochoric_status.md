# `build_isochoric` Phase Envelope Tracer — Status & Handoff Notes

**Branch**: `feat/build-isochoric-phase-envelope`  
**Date stashed**: 2026-04-13  
**Author**: Ian Bell / Claude (claude-sonnet-4-6)

---

## Goal

Implement `PhaseEnvelopeRoutines::build_isochoric(HEOS)` that traces the **complete** phase
envelope — dew curve, critical point, and bubble curve — for a mixture using isochoric
thermodynamic derivatives (DmolarT_INPUTS), rather than the TP-flash inner loop used by the
existing `build()`.

The motivating application is speed and numerical robustness for cases where a TP flash is
expensive or ill-conditioned near the critical point.

## Benchmark

```
cmake .. -DCOOLPROP_MY_MAIN=../dev/bench_phase_envelope.cpp -DCOOLPROP_STATIC_LIBRARY=ON
ninja Main && ./Main
```

Five mixtures: CH4/C2H6 50/50, CH4/nC4 80/20, CH4/C3H8/N2, CH4/C2H6/C3H8, CH4/CO2/N2.

## Algorithmic History (what was tried)

### 1. Pseudo-arclength continuation (first approach, abandoned)

State vector **s** = (x[0..N-2], T, ρ_L, ρ_V).  
Equilibrium residuals **F** (N+1): fugacity balance + pressure balance.  
Predictor: `s_pred = s + h · v` where `v = null(J)` via SVD.  
Corrector: augmented (N+2)×(N+2) Newton system:
```
  [ J(s)  ] [Δs] = [ -F(s)              ]
  [ v^T   ]        [ -v·(s - s_pred)    ]
```

**Why abandoned**: The cricondentherm (temperature maximum, distinct from the critical point)
is a turning point in T. After the cricondentherm the SVD null vector's T-component must flip
sign to follow the descending branch; however the algorithm found a spurious branch where T
continued to increase toward the pure-component critical temperature (~850 K), producing
thousands of tiny steps and an 18-minute hang. Detection of the sign flip requires already
being on the correct branch — a chicken-and-egg problem.

### 2. ρ_V-imposed Newton marching (current approach, partially working)

**Key insight**: ρ_V increases monotonically along the *entire* phase envelope — dew curve
outbound, through the critical point (ρ_V_crit = ρ_L_crit), and along the bubble curve
inbound. This is the natural parameter that sidesteps the cricondentherm problem: there is no
turning point in ρ_V.

For each step:
1. **Predict**: multiply ρ_V by an adaptive factor (1.001 … 1.08).
2. **Correct**: call `newton_raphson_saturation` with `RHOV_IMPOSED`, solving for (x[0..N-2],
   T, ρ_L) with ρ_V fixed.
3. **Accept/reject** based on sanity checks (ρ_L ≈ ρ_V, ΔT > 100 K, negative p).
4. **Terminate** when p < p_start (loop closed) or a pure component is approached.

The initial guess for each step uses linear extrapolation in ρ_V from the last two stored
points. Near the cricondentherm the extrapolated T overshoots the temperature maximum, so we
guard: if `T_extrap > env.T.back() + 5 K` we fall back to the last accepted state.

**Inner solver**: `SaturationSolvers::newton_raphson_saturation::call(HEOS, y, x, IO)` with
`RHOV_IMPOSED`. This is exactly what `build()` uses, so we are sure the inner solver is
correct — only the outer loop is new.

## Current Results (as of stash date)

```
mixture                  method           N_pts    time_ms  built  T_max_K  p_max_bar
CH4/C2H6 50/50           build              217      27.3      1   267.35     67.77
CH4/C2H6 50/50           build_isochoric    203     310.5      0   267.35     67.77  ✗
CH4/nC4 80/20            build              221       9.1      1   336.49    142.39
CH4/nC4 80/20            build_isochoric    234      38.5      1   336.49    142.39  ✓
CH4/C3H8/N2              build              ---       ---      0      ---       ---   (both fail)
CH4/C2H6/C3H8            build              226      13.2      1   270.75     82.73
CH4/C2H6/C3H8            build_isochoric    220     380.1      0   270.75     82.73  ✗
CH4/CO2/N2               build              ---       ---      0      ---       ---   (both fail)
```

**CH4/nC4**: fully working — correct T_max, p_max, and `built=1`.  
**CH4/C2H6 and CH4/C2H6/C3H8**: reach the cricondentherm correctly (T_max matches `build()`
exactly) but `newton_raphson_saturation` fails consistently for all guesses after the
cricondentherm. After `failure_max=8` rejections the loop exits with `built=0`.  
**CH4/C3H8/N2 and CH4/CO2/N2**: both methods fail; not a regression.

## Root Cause Analysis for Remaining Failures

### Why CH4/nC4 works but CH4/C2H6 does not

After the cricondentherm T must decrease. For CH4/nC4 the extrapolation guard fires early
enough that the fallback guess (previous accepted state) gives Newton a starting point inside
the convergence basin. For CH4/C2H6 and CH4/C2H6/C3H8 the extrapolated T exceeds the 5 K
guard threshold, the fallback guess is used, but Newton still diverges — likely because the
convergence basin of `newton_raphson_saturation` with `RHOV_IMPOSED` is narrow near the
cricondentherm for these compositions.

### Jacobian normalization note

`envelope_jac_and_res` (the custom isochoric Jacobian, still present in the file but no longer
used as the inner corrector) normalizes the pressure row by p_V:

```cpp
F(N) = (p_L - p_V) / p_V;      // dimensionless
```

`build()`'s `build_arrays` does NOT normalize:
```cpp
r(N) = p_L - p_V;               // in Pa
```

This conditioning difference causes `envelope_jac_and_res`-based Newton to fail past the
cricondentherm for CH4/nC4 even when NR.call succeeds. If the isochoric Jacobian is ever
reinstated as the inner solver, try removing the `/p_scale` normalization (match `build_arrays`
conditioning) or use row scaling that balances the fugacity and pressure residuals differently.

## What a Future Agent Should Try

### Option A: Larger failure_max + smarter cricondentherm detection
Increase `failure_max` from 8 to 20+. Near the cricondentherm, detect that T is near its
maximum (e.g. `env.T.back() < env.T[env.T.size()-2]` — T started decreasing) and explicitly
set `T_c = env.T.back() - ΔT` instead of extrapolating. This is the most targeted fix.

### Option B: Hybrid switching at cricondentherm
Detect the cricondentherm (dT/dρ_V changes sign between successive accepted steps). Switch
parametrisation to T-imposed (solve for ρ_L, ρ_V with T fixed and decreasing) for the
descending branch, then back to ρ_V-imposed once ρ_V starts decreasing again.

### Option C: Full arc-length with better branch tracking
Return to pseudo-arclength but add explicit cricondentherm detection: when the stored T
sequence forms a maximum (`env.T.back() < env.T[env.T.size()-2]`), force the T-component of
the tangent to be negative and re-orthogonalize via Gram-Schmidt. This is more complex but
avoids the parametrisation switch.

### Option D: Match build() exactly for the outer loop
Copy `build()`'s outer loop verbatim (including its `factor` adaptation and stopping logic)
but replace `build_arrays` with `envelope_jac_and_res` as the inner Jacobian. The current
implementation already uses NR.call as the inner corrector — this is identical to `build()`.
The only difference left would be the outer stepping variable (ρ_V vs. whatever `build()` uses
internally). Audit `build()` carefully: it also steps ρ_V with RHOV_IMPOSED. If it works for
CH4/C2H6, the difference must be in the initial guess or the stopping/retry logic.

## Files Modified

- `src/Backends/Helmholtz/PhaseEnvelopeRoutines.cpp` — Primary implementation.
  - Anonymous namespace helpers (lines ~849–931):
    - `envelope_jac_and_res`: computes F (N+1) and J (N+1)×(N+2) using isochoric derivatives.
      Present but **not used** in the current outer loop (uses NR.call instead).
    - `pack_arclen`, `unpack_arclen`, `tangent_from_J`: arc-length helpers, **not used**,
      retained for Option C above.
  - `build_isochoric` function (lines ~933–1169): ρ_V-imposed marching loop.
- `src/Backends/Helmholtz/PhaseEnvelopeRoutines.h` — Declaration of `build_isochoric`.
- `dev/bench_phase_envelope.cpp` — Benchmark driver; `fflush(stdout)` after each row.

## Files NOT Modified by This Work (but dirty in the worktree)

- `src/Backends/Helmholtz/ExcessHEFunction.h` — `phi.finish()` fix (prior session, separate).
- `src/Backends/Helmholtz/FlashRoutines.cpp`, `FlashRoutines.h` — GDEM/HELD flash work
  (committed on `feat/gdem-held-flash`, but worktree still has uncommitted tweaks).
