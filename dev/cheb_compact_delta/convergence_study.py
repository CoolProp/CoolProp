"""Per-term Chebyshev-in-u convergence study for IAPWS-95.

For each term family (power l=0, power l>0, Gaussian bell, non-analytic),
we build the pressure-contribution quantity at fixed tau:
    f_i(delta) = delta * d(alpha^r_i)/d(delta)     (i.e. alphar_01 contribution)
then study Chebyshev-in-u approximation convergence, where
    u = (delta - 1) / (delta + 1)
on the interval [-1+eps, 1-eps].

Outputs:
  1. Coefficient-decay plots per family for a few representative tau values.
  2. Max-error-vs-degree curve per family.
  3. A scan over all 56 IAPWS-95 terms giving degree needed to reach 1e-10
     relative max error on u in [-1+1e-4, 1-1e-4].
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from numpy.polynomial import chebyshev as npcheb

from iapws95 import IAPWS95Terms
from compact_cheb import (
    chebyshev_lobatto_nodes,
    values_to_coeffs,
    delta_of_u,
)


OUTDIR = os.path.dirname(os.path.abspath(__file__))


def fit_cheb_at_degree(f, a, b, n):
    """Fit f on [a, b] at Chebyshev-Lobatto nodes, return coefficients."""
    nodes = chebyshev_lobatto_nodes(n, a, b)
    vals = np.array([f(x) for x in nodes])
    return values_to_coeffs(vals), vals


def max_rel_error(f, a, b, coeffs, n_test=4000):
    """Max |f - cheb_eval| / max(|f|, abstol) on a dense grid inside (a,b)."""
    xs = np.linspace(a + 1e-10, b - 1e-10, n_test)
    # map to [-1, 1] for chebval
    xs_std = 2 * (xs - a) / (b - a) - 1.0
    fs = np.array([f(x) for x in xs])
    ps = npcheb.chebval(xs_std, coeffs)
    fmax = np.max(np.abs(fs))
    abstol = max(fmax * 1e-16, 1e-300)
    return np.max(np.abs(fs - ps)) / max(fmax, abstol), np.max(np.abs(fs - ps)), fmax


def study_one_term(label, f_term, a=-0.9999, b=0.9999,
                   degrees=(8, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512),
                   target_rel=1e-10):
    """Return (coeffs at highest degree, error table list)."""
    table = []
    final_coeffs = None
    for n in degrees:
        coeffs, vals = fit_cheb_at_degree(f_term, a, b, n)
        rel, absv, fmax = max_rel_error(f_term, a, b, coeffs, n_test=2000)
        table.append((n, rel, absv, fmax))
        final_coeffs = coeffs
        if rel < target_rel:
            break
    return final_coeffs, table


def family_convergence_plots(eos, tau_values, eps=1e-4):
    """For each representative term, produce a coefficient-magnitude plot and
    an error-vs-degree plot."""
    # Representative term selections (by index into each family's arrays):
    reps = [
        ("power_l0_idx5", "power term #5 (polynomial l=0)",
         lambda eos: _single_power_term_fn(eos, 5)),
        ("power_l1_idx10", "power term #10 (l=1, poly-exp)",
         lambda eos: _single_power_term_fn(eos, 10)),
        ("power_l2_idx28", "power term #28 (l=2, poly-exp)",
         lambda eos: _single_power_term_fn(eos, 28)),
        ("power_l6_idx50", "power term #50 (l=6, poly-exp)",
         lambda eos: _single_power_term_fn(eos, 50)),
        ("gauss_idx0", "Gaussian bell term #0",
         lambda eos: _single_gauss_term_fn(eos, 0)),
        ("gauss_idx2", "Gaussian bell term #2",
         lambda eos: _single_gauss_term_fn(eos, 2)),
        ("nonan_idx0", "non-analytic term #0",
         lambda eos: _single_nonan_term_fn(eos, 0)),
        ("nonan_idx1", "non-analytic term #1",
         lambda eos: _single_nonan_term_fn(eos, 1)),
    ]

    fig_c, axes_c = plt.subplots(4, 2, figsize=(14, 14))
    fig_e, axes_e = plt.subplots(4, 2, figsize=(14, 14))
    axes_c = axes_c.flatten()
    axes_e = axes_e.flatten()

    for idx, (key, desc, make_fn) in enumerate(reps):
        ax_c = axes_c[idx]
        ax_e = axes_e[idx]
        for tau in tau_values:
            f_fixed = make_fn(eos)  # returns fn(u) at this particular term, then we parameterize by tau
            def f_of_u_at_tau(u, _tau=tau, _f=f_fixed):
                return _f(u, _tau)
            # Single high-degree fit for coefficient plot
            coeffs, vals = fit_cheb_at_degree(f_of_u_at_tau, -1 + 1e-4, 1 - 1e-4, 256)
            ax_c.semilogy(np.abs(coeffs) + 1e-300, label=f'tau={tau:.3f}', lw=0.8)

            # Error-vs-degree sweep
            _, tbl = study_one_term(key, f_of_u_at_tau, a=-1+1e-4, b=1-1e-4,
                                    degrees=(8, 16, 32, 48, 64, 96, 128, 192, 256, 384))
            ns = [row[0] for row in tbl]
            rels = [row[1] for row in tbl]
            ax_e.semilogy(ns, rels, 'o-', label=f'tau={tau:.3f}', lw=0.8, ms=3)

        ax_c.set_title(f'{desc}: |coef|', fontsize=9)
        ax_c.set_xlabel('Chebyshev index')
        ax_c.set_ylabel('|c_k|')
        ax_c.grid(alpha=0.3, which='both')
        ax_c.legend(fontsize=7)
        ax_c.set_ylim(1e-20, 1e30)

        ax_e.set_title(f'{desc}: max rel error', fontsize=9)
        ax_e.set_xlabel('degree n')
        ax_e.set_ylabel('max rel err')
        ax_e.grid(alpha=0.3, which='both')
        ax_e.legend(fontsize=7)
        ax_e.axhline(1e-10, color='r', lw=0.5, ls='--')

    fig_c.suptitle('Per-term Chebyshev coefficient magnitudes in u (compactified delta)')
    fig_e.suptitle('Per-term Chebyshev approx error vs degree in u (compactified delta)')
    fig_c.tight_layout()
    fig_e.tight_layout()
    fig_c.savefig(os.path.join(OUTDIR, "term_coeffs.png"), dpi=120)
    fig_e.savefig(os.path.join(OUTDIR, "term_errors.png"), dpi=120)
    plt.close(fig_c)
    plt.close(fig_e)


def _single_power_term_fn(eos, idx):
    """Return f(u, tau) for the i-th power term's delta * d(alphar)/d(delta)."""
    n = eos.n_p[idx]; d = eos.d_p[idx]; t = eos.t_p[idx]; l = eos.l_p[idx]

    def f(u, tau):
        delta = delta_of_u(u)
        tp = tau ** t
        dp = delta ** d
        ep = np.exp(-delta ** l) if l > 0 else 1.0
        factor = d - l * (delta ** l)
        return n * dp * tp * ep * factor
    return f


def _single_gauss_term_fn(eos, idx):
    n = eos.n_g[idx]; d = eos.d_g[idx]; t = eos.t_g[idx]
    eta = eos.eta_g[idx]; eps = eos.eps_g[idx]
    beta = eos.beta_g[idx]; gamma = eos.gamma_g[idx]

    def f(u, tau):
        delta = delta_of_u(u)
        arg = -eta * (delta - eps) ** 2 - beta * (tau - gamma) ** 2
        inner = d - 2 * eta * delta * (delta - eps)
        return n * delta ** d * tau ** t * np.exp(arg) * inner
    return f


def _single_nonan_term_fn(eos, idx):
    """Return f(u, tau) for the i-th non-analytic term's δ · dα^r/dδ."""
    n = eos.n_na[idx]; a = eos.a_na[idx]; b = eos.b_na[idx]
    beta = eos.beta_na[idx]; A = eos.A_na[idx]; B = eos.B_na[idx]
    C = eos.C_na[idx]; D = eos.D_na[idx]

    def f(u, tau):
        delta = delta_of_u(u)
        dm1 = delta - 1.0
        dm1sq = dm1 * dm1
        tiny = 1e-300
        dm1sq_safe = dm1sq if dm1sq > 0 else tiny
        theta = (1.0 - tau) + A * np.power(dm1sq, 1.0 / (2.0 * beta)) if dm1sq > 0 else (1.0 - tau)
        Delta = theta ** 2 + B * np.power(dm1sq, a) if dm1sq > 0 else theta ** 2
        psi = np.exp(-C * dm1sq - D * (tau - 1.0) ** 2)
        # derivatives
        if dm1sq == 0:
            dtheta_dd = 0.0
            dDelta_dd = 0.0
        else:
            dtheta_dd = (A / beta) * np.power(dm1sq, 1.0 / (2.0 * beta) - 1.0) * dm1
            dDelta_dd = 2 * theta * dtheta_dd + \
                2 * a * B * np.power(dm1sq, a - 1.0) * dm1
        dpsi_dd = -2 * C * dm1 * psi

        term1 = np.power(Delta, b) * (psi + delta * dpsi_dd)
        term2 = b * np.power(Delta, b - 1.0) * dDelta_dd * delta * psi
        dalphar_dd = n * (term1 + term2)
        return delta * dalphar_dd
    return f


def full_term_scan(eos, tau_values=(1.5, 1.0, 0.8), eps=1e-4,
                   target_rel=1e-10, n_max=512):
    """Report, for every single term across all families, the degree
    needed to reach `target_rel` max relative error, or 'FAIL' if n_max
    isn't enough. Returns a list of dicts."""
    out = []
    degrees_try = [8, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512]
    for idx in range(len(eos.n_p)):
        term_name = f"power_l{int(eos.l_p[idx])}_i{idx}"
        fn = _single_power_term_fn(eos, idx)
        out.append(_scan_term(term_name, fn, tau_values, degrees_try, target_rel))
    for idx in range(len(eos.n_g)):
        term_name = f"gauss_i{idx}"
        fn = _single_gauss_term_fn(eos, idx)
        out.append(_scan_term(term_name, fn, tau_values, degrees_try, target_rel))
    for idx in range(len(eos.n_na)):
        term_name = f"nonan_i{idx}"
        fn = _single_nonan_term_fn(eos, idx)
        out.append(_scan_term(term_name, fn, tau_values, degrees_try, target_rel))
    return out


def _scan_term(term_name, fn, tau_values, degrees_try, target_rel, eps=1e-4):
    res = dict(term=term_name, per_tau={})
    for tau in tau_values:
        def f_of_u(u, _tau=tau, _f=fn):
            return _f(u, _tau)
        best = None
        for n in degrees_try:
            coeffs, _ = fit_cheb_at_degree(f_of_u, -1 + eps, 1 - eps, n)
            rel, absv, fmax = max_rel_error(f_of_u, -1 + eps, 1 - eps, coeffs, n_test=1000)
            if rel < target_rel:
                best = (n, rel, absv, fmax)
                break
        if best is None:
            # report worst-seen
            best = (degrees_try[-1], rel, absv, fmax)
        res["per_tau"][tau] = best
    return res


def main():
    eos = IAPWS95Terms()
    print(f"IAPWS-95 loaded: Tc={eos.Tc}, rhoc={eos.rhoc:.4f} mol/m^3")
    print(f"  power terms: {len(eos.n_p)}  (l=0: {sum(eos.l_p == 0)}, "
          f"l>0: {sum(eos.l_p > 0)})")
    print(f"  gaussian terms: {len(eos.n_g)}")
    print(f"  non-analytic terms: {len(eos.n_na)}")

    tau_vals = [1.5, 1.0, 0.8]  # subcritical (T=431K), critical, supercritical (T=809K)

    print("\n=== Generating per-term coefficient and error plots ===")
    family_convergence_plots(eos, tau_vals)
    print(f"  -> {OUTDIR}/term_coeffs.png")
    print(f"  -> {OUTDIR}/term_errors.png")

    print("\n=== Full per-term degree scan (target rel err = 1e-10) ===")
    results = full_term_scan(eos, tau_values=tau_vals, target_rel=1e-10)
    # print table
    print(f"\n{'term':<20} " + "  ".join([f"tau={tv:<5.3f} (n / rel)"
                                          for tv in tau_vals]))
    for r in results:
        row = f"{r['term']:<20} "
        for tv in tau_vals:
            n, rel, absv, fmax = r['per_tau'][tv]
            row += f"  {n:>4d} / {rel:.1e}   "
        print(row)


if __name__ == "__main__":
    main()
