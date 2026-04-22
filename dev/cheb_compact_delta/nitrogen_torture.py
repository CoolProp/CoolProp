"""Comprehensive torture test for nitrogen density rootfinder.

Stresses the architecture across:
  - Near-critical (T ≈ Tc = 126.192 K) where dP/dρ → 0
  - Deep subcritical two-phase (T = 70-110K) with 3 physical roots
  - Triple-point vicinity (T ≈ 63.151 K)
  - Extreme pressures (> 100 MPa, pushing δ → δ_max)
  - Very low P (ideal-gas limit, δ → 0)
  - Supercritical (T up to 1000 K)

For each (T, P) case we compare to CoolProp's reference root(s) and tally errors.
Reports: max error across sweep, cases where we miss a root, cases where we
return extra spurious roots, worst-case timing.
"""

import time
import numpy as np
import CoolProp.CoolProp as CP

from fluid_from_json import FluidFromJSON
from fast_solver import build_fluid_fast, density_roots_vectorized


def stable_filter(fluid, T, roots):
    """Keep only roots where dP/drho > 0 (thermodynamically stable)."""
    if len(roots) == 0:
        return roots
    tau = fluid.Tc / T
    stable = []
    for rho in roots:
        delta = rho / fluid.rhoc
        # Numerical dP/dρ from direct EOS evaluation: P = ρRT(1 + α_01)
        # dP/drho = RT * (1 + α_01) + ρ R T * dα_01/drho
        # = RT * (1 + α_01 + δ α'_01 / rho_c * rho_c) = RT * (1 + 2 α_01 + α_02)
        # Use finite-difference in rho
        h = max(1e-3, abs(rho) * 1e-6)
        P_p = fluid.pressure_direct(T, rho + h)
        P_m = fluid.pressure_direct(T, max(rho - h, 1e-12))
        dPdrho = (P_p - P_m) / (2 * h)
        if dPdrho > 0:
            stable.append(rho)
    return np.array(stable)


def coolprop_reference(T, P, fluid):
    """Return CoolProp reference densities at (T, P). For a given (T, P) there
    is at most one stable root per phase. Try to return (vapor, liquid, single)
    depending on where the point sits relative to the saturation curve."""
    M = fluid.M
    T_c = fluid.Tc
    try:
        if T < T_c:
            Psat = CP.PropsSI('P', 'T', float(T), 'Q', 0.0, fluid.name)
            if P < Psat:
                return [CP.PropsSI('D', 'T', float(T), 'P', float(P), fluid.name) / M]
            elif P > Psat:
                return [CP.PropsSI('D', 'T', float(T), 'P', float(P), fluid.name) / M]
            else:  # exactly at Psat: both phases
                rhoL = CP.PropsSI('D', 'T', float(T), 'Q', 0.0, fluid.name) / M
                rhoV = CP.PropsSI('D', 'T', float(T), 'Q', 1.0, fluid.name) / M
                return [rhoV, rhoL]
        else:
            return [CP.PropsSI('D', 'T', float(T), 'P', float(P), fluid.name) / M]
    except Exception:
        return None


def run_sweep(fluid, state, T_vals, P_vals, label=""):
    """Run rootfinder at all (T, P) combinations and tally errors against CoolProp."""
    print(f"\n=== Sweep: {label} (T points={len(T_vals)}, P points={len(P_vals)}) ===")
    results = []
    worst_rel = 0.0
    worst_rel_case = None
    missed = 0
    extra = 0
    failed = 0
    timings = []

    for T in T_vals:
        for P in P_vals:
            try:
                t0 = time.perf_counter()
                roots, stats = density_roots_vectorized(state, T, P)
                dt = time.perf_counter() - t0
                timings.append(dt * 1e6)
            except Exception as e:
                failed += 1
                continue

            ref = coolprop_reference(T, P, fluid)
            if ref is None:
                continue  # CoolProp didn't give a reference

            # Filter our roots by stability
            stable = stable_filter(fluid, T, roots)
            stable = np.sort(stable)
            ref_arr = np.sort(np.array(ref))

            # Match each reference root to nearest in stable
            for rho_ref in ref_arr:
                if len(stable) == 0:
                    missed += 1
                    continue
                i = np.argmin(np.abs(stable - rho_ref))
                rel_err = abs(stable[i] - rho_ref) / max(abs(rho_ref), 1e-10)
                if rel_err > worst_rel:
                    worst_rel = rel_err
                    worst_rel_case = (T, P, rho_ref, stable[i])
                # Tolerance: allow 1% (seed quality)
                if rel_err > 1e-2:
                    missed += 1

            # Count stable roots that don't match any reference
            for rho_ours in stable:
                if len(ref_arr) == 0 or np.min(np.abs(ref_arr - rho_ours)) / max(abs(rho_ours), 1e-10) > 1e-2:
                    extra += 1

            results.append((T, P, stats, roots, stable, ref_arr))

    timings = np.array(timings)
    print(f"  cases: {len(results)} ok, {failed} failed")
    print(f"  worst rel err vs CP: {worst_rel:.2e}")
    if worst_rel_case:
        T, P, rho_ref, rho_ours = worst_rel_case
        print(f"     at T={T:.2f}, P={P:.3e}, ref={rho_ref:.4e}, ours={rho_ours:.4e}")
    print(f"  missed (>1% off): {missed},   extra stable roots: {extra}")
    print(f"  timing: median={np.median(timings):.1f}, "
          f"p95={np.percentile(timings, 95):.1f}, "
          f"max={np.max(timings):.1f} us")
    return results


def main():
    path = "/Users/ianbell/miniforge3/lib/python3.9/site-packages/teqp/fluiddata/dev/fluids/Nitrogen.json"
    fluid = FluidFromJSON(path)
    print(f"Nitrogen: Tc={fluid.Tc} K, rhoc={fluid.rhoc:.3f} mol/m^3, M={fluid.M}")
    Ttriple = 63.151

    t0 = time.perf_counter()
    state = build_fluid_fast(fluid, 1e-10, 5.0, n_delta=12, tol=1e-10)
    print(f"Built state: K={state['K']} intervals in {(time.perf_counter()-t0)*1e3:.1f} ms")

    # SWEEP 1: Supercritical wide pressure
    T_sup = np.linspace(150, 1000, 18)
    P_sup = np.logspace(3, 8, 16)  # 1 kPa to 100 MPa
    run_sweep(fluid, state, T_sup, P_sup, "supercritical wide")

    # SWEEP 2: Subcritical, away from Psat
    T_sub = np.linspace(70, 125, 12)
    P_sub = np.logspace(2, 7, 16)  # 100 Pa to 10 MPa
    run_sweep(fluid, state, T_sub, P_sub, "subcritical single-phase")

    # SWEEP 3: Near-critical (|T - Tc| < 5K) where the roots are numerically hard
    T_nc = np.linspace(122, 130, 9)
    P_nc = np.linspace(2.5e6, 4.5e6, 11)  # around Pc ≈ 3.39 MPa
    run_sweep(fluid, state, T_nc, P_nc, "near-critical")

    # SWEEP 4: Very low temperature near triple (tests δ_max boundary for dense liquid)
    T_lo = np.linspace(Ttriple + 0.1, 80, 10)
    P_lo = np.logspace(3, 8, 10)
    run_sweep(fluid, state, T_lo, P_lo, "near-triple")

    # SWEEP 5: Saturation curve sampling — explicit P=Psat test (3 roots expected)
    T_sat = np.linspace(70, 125, 20)
    print(f"\n=== Saturation-curve test (P = Psat exactly, expect both rhoL and rhoV) ===")
    n_matched = n_missed_V = n_missed_L = 0
    worst_V = 0.0; worst_L = 0.0
    for T in T_sat:
        try:
            Psat = CP.PropsSI('P', 'T', float(T), 'Q', 0.0, fluid.name)
            rhoL = CP.PropsSI('D', 'T', float(T), 'Q', 0.0, fluid.name) / fluid.M
            rhoV = CP.PropsSI('D', 'T', float(T), 'Q', 1.0, fluid.name) / fluid.M
        except Exception:
            continue
        roots, stats = density_roots_vectorized(state, T, Psat)
        stable = stable_filter(fluid, T, roots)
        if len(stable) < 2:
            # Did we find vapor at least?
            if len(stable) == 1:
                rel = abs(stable[0] - rhoV) / rhoV
                if rel < 1e-2: n_missed_L += 1
                else: n_missed_V += 1
            continue
        # match to rhoV (min) and rhoL (max)
        rhoV_ours = stable.min()
        rhoL_ours = stable.max()
        rv = abs(rhoV_ours - rhoV) / rhoV
        rl = abs(rhoL_ours - rhoL) / rhoL
        if rv < 1e-2 and rl < 1e-2: n_matched += 1
        else:
            if rv > 1e-2: n_missed_V += 1
            if rl > 1e-2: n_missed_L += 1
        worst_V = max(worst_V, rv); worst_L = max(worst_L, rl)
    print(f"  matched both phases: {n_matched}/{len(T_sat)}")
    print(f"  missed rhoV: {n_missed_V}, missed rhoL: {n_missed_L}")
    print(f"  worst rhoV err: {worst_V:.2e}, worst rhoL err: {worst_L:.2e}")


if __name__ == "__main__":
    main()
