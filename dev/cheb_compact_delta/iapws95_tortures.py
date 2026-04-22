"""
IAPWS-95 torture tests for compactified piecewise Chebyshev rootfinding.

For each (T, P_target) point:
  - construct f(u) = P_eos(delta(u)) - P_target
  - piecewise-adaptively approximate on u in [-1+eps, 1-eps]
  - find all roots; map back to delta and then rho [mol/m^3]
  - plot the piecewise approximation overlaid on the raw residual

teqp is used solely as the alpha^r evaluator.
"""

import os
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import teqp

from compact_cheb import (
    delta_of_u,
    u_of_delta,
    all_roots_piecewise,
    build_piecewise,
)

R_GAS = 8.314462618  # J/(mol K); teqp get_R segfaults in our env, use CODATA


class IAPWS95:
    def __init__(self):
        self.model = teqp.build_multifluid_model(["Water"], teqp.get_datapath())
        z = np.array([1.0])
        self.z = z
        self.Tc = self.model.get_Tr(z)  # 647.096
        self.rhoc = self.model.get_rhor(z)  # 17873.73 mol/m^3
        self.R = R_GAS

    def pressure(self, T, rho):
        """Raw EOS pressure, Pa. Scalar T, scalar or array rho."""
        if np.isscalar(rho):
            Ar01 = self.model.get_Ar01(float(T), float(rho), self.z)
            return rho * self.R * T * (1.0 + Ar01)
        out = np.empty_like(rho, dtype=float)
        for i, r in enumerate(rho):
            Ar01 = self.model.get_Ar01(float(T), float(r), self.z)
            out[i] = r * self.R * T * (1.0 + Ar01)
        return out

    def dP_drho(self, T, rho, h=None):
        """Finite-difference dP/drho in log space for smoothness reports."""
        # Used only for diagnostics, so FD is fine.
        if h is None:
            h = max(1e-3, abs(rho) * 1e-6)
        return (self.pressure(T, rho + h) - self.pressure(T, rho - h)) / (2 * h)


def residual_in_u(eos, T, P_target, L=1.0):
    """Return a scalar function f(u) = P_eos(delta(u)*rho_c, T) - P_target [Pa]."""

    def f(u):
        delta = delta_of_u(u, L=L)
        rho = delta * eos.rhoc
        return eos.pressure(T, rho) - P_target

    return f


def run_case(eos, T, P_target, eps=1e-4, tol=1e-10, n_init=16, n_max=128,
             min_width=1e-6, plot_path=None, case_label=""):
    f = residual_in_u(eos, T, P_target)
    a, b = -1.0 + eps, 1.0 - eps

    t0 = time.perf_counter()
    roots_u, pieces = all_roots_piecewise(
        f, a, b, tol=tol, n_init=n_init, n_max=n_max, min_width=min_width
    )
    dt = time.perf_counter() - t0

    # map back
    delta_roots = np.array([delta_of_u(u) for u in roots_u])
    rho_roots = delta_roots * eos.rhoc

    # sanity: compute residual at each root
    residuals = np.array([f(u) for u in roots_u])

    print(f"--- {case_label}: T={T:.2f} K (tau={eos.Tc/T:.3f}), P={P_target:.3e} Pa ---")
    print(f"  pieces={len(pieces)}  roots_u={len(roots_u)}  solve_time={dt*1e3:.2f} ms")
    for i, (u, d, rho, res) in enumerate(zip(roots_u, delta_roots, rho_roots, residuals)):
        # slope for phase classification
        dP = eos.dP_drho(T, rho) if rho > 0 else np.nan
        print(f"    root {i}: u={u:+.6f} delta={d:.6f} rho={rho:11.4f} mol/m^3   "
              f"P_res={res:+.3e} Pa  dP/drho={dP:+.3e}  "
              f"{'[stable]' if dP > 0 else '[unstable]'}")

    if plot_path is not None:
        plot_case(eos, T, P_target, pieces, roots_u, plot_path, case_label)

    return dict(T=T, P=P_target, roots_u=roots_u, roots_rho=rho_roots,
                pieces=pieces, time_ms=dt*1e3)


def plot_case(eos, T, P_target, pieces, roots_u, path, case_label):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 9))

    # Dense raw residual in u
    u_dense = np.linspace(-0.9999, 0.9999, 6000)
    rho_dense = delta_of_u(u_dense) * eos.rhoc
    P_dense = eos.pressure(T, rho_dense)
    res_dense = P_dense - P_target

    # clip for plotting (extreme values)
    ylim = 3 * abs(P_target) if P_target != 0 else 1e7
    ylim = max(ylim, 1e7)
    res_clipped = np.clip(res_dense, -ylim, ylim)

    ax1.plot(u_dense, res_clipped, 'k-', lw=0.5, label='EOS residual (clipped)')
    # overlay piecewise approximation
    for (aa, bb, coeffs) in pieces:
        uu = np.linspace(aa, bb, 64)
        xx = 2 * (uu - aa) / (bb - aa) - 1.0
        from numpy.polynomial import chebyshev as npcheb
        vv = npcheb.chebval(xx, coeffs)
        vv_clipped = np.clip(vv, -ylim, ylim)
        ax1.plot(uu, vv_clipped, lw=1.0, alpha=0.8)
        ax1.axvline(aa, color='gray', ls=':', lw=0.4, alpha=0.3)
        ax1.axvline(bb, color='gray', ls=':', lw=0.4, alpha=0.3)
    for ur in roots_u:
        ax1.axvline(ur, color='red', ls='--', lw=0.7)
    ax1.axhline(0, color='k', lw=0.3)
    ax1.set_xlabel('u = (delta-1)/(delta+1)')
    ax1.set_ylabel('P_eos - P_target  [Pa]  (clipped)')
    ax1.set_title(f'{case_label}: T={T:.2f} K, P_target={P_target:.3e} Pa, '
                  f'{len(pieces)} pieces, {len(roots_u)} roots')
    ax1.legend(fontsize=8)
    ax1.grid(alpha=0.3)

    # Delta-space view with log-symlog pressure
    ax2.set_xscale('log')
    rho_plot = np.logspace(-3, np.log10(80_000), 5000)
    P_plot = eos.pressure(T, rho_plot)
    ax2.plot(rho_plot, np.clip(P_plot - P_target, -ylim, ylim), 'k-', lw=0.5,
             label=f'P_eos - P_target (clipped to ±{ylim:.1e})')
    for ur in roots_u:
        rho_r = delta_of_u(ur) * eos.rhoc
        ax2.axvline(rho_r, color='red', ls='--', lw=0.7)
    ax2.axhline(0, color='k', lw=0.3)
    ax2.set_xlabel('rho [mol/m^3]')
    ax2.set_ylabel('P_eos - P_target [Pa]')
    ax2.grid(alpha=0.3, which='both')
    ax2.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(path, dpi=120)
    plt.close(fig)


def reference_saturation(eos, T):
    """Use CoolProp's IAPWS-95 saturation as reference (hard to beat)."""
    import CoolProp.CoolProp as CP
    try:
        if T >= eos.Tc:
            return None, None, None
        Psat = CP.PropsSI("P", "T", float(T), "Q", 0.0, "Water")  # Pa
        rhoL_mass = CP.PropsSI("D", "T", float(T), "Q", 0.0, "Water")
        rhoV_mass = CP.PropsSI("D", "T", float(T), "Q", 1.0, "Water")
        M = 18.015268e-3
        return Psat, rhoL_mass / M, rhoV_mass / M
    except Exception:
        return None, None, None


def main():
    eos = IAPWS95()
    outdir = os.path.dirname(os.path.abspath(__file__))
    print(f"IAPWS-95: Tc={eos.Tc} K, rhoc={eos.rhoc:.3f} mol/m^3")

    # Build a set of torture cases.
    # Strategy: pick T spanning triple to supercritical; for each, test
    # P_target at Psat (both vapor and liquid roots should appear), and
    # above/below Psat.
    cases = []
    for T in [280.0, 400.0, 500.0, 600.0, 640.0, 647.0, 700.0, 900.0]:
        Psat, rhoL, rhoV = reference_saturation(eos, T)
        if Psat is not None and Psat > 0:
            cases.append((T, Psat * 0.5, f"T={T}_belowsat"))
            cases.append((T, Psat, f"T={T}_atsat"))
            cases.append((T, Psat * 2.0, f"T={T}_abovesat"))
        else:
            # supercritical or saturation failed: sweep a few P
            cases.append((T, 1e6, f"T={T}_1MPa"))
            cases.append((T, 1e7, f"T={T}_10MPa"))
            cases.append((T, 1e8, f"T={T}_100MPa"))

    results = []
    for T, P, label in cases:
        Psat, rhoL, rhoV = reference_saturation(eos, T)
        if Psat is not None:
            print(f"[ref] T={T}: Psat={Psat:.4e} Pa, rhoL={rhoL:.2f}, rhoV={rhoV:.4f} mol/m^3")
        plot_path = os.path.join(outdir, f"torture_{label}.png")
        try:
            r = run_case(eos, T, P, plot_path=plot_path, case_label=label)
            results.append(r)
        except Exception as e:
            print(f"  FAILED: {e}")
            results.append(dict(T=T, P=P, error=str(e)))

    # Summary
    print("\n=== Summary ===")
    print(f"{'Case':<26} {'T [K]':>8} {'P [Pa]':>12} {'pieces':>7} "
          f"{'roots':>6} {'time [ms]':>10}")
    for c, r in zip(cases, results):
        _, _, label = c
        if 'error' in r:
            print(f"{label:<26} {r['T']:>8.2f} {r['P']:>12.3e}   ERROR: {r['error']}")
        else:
            print(f"{label:<26} {r['T']:>8.2f} {r['P']:>12.3e} "
                  f"{len(r['pieces']):>7d} {len(r['roots_u']):>6d} {r['time_ms']:>10.2f}")


if __name__ == "__main__":
    main()
