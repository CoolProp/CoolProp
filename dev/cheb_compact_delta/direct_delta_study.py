"""Direct-variable piecewise Chebyshev study for IAPWS-95 per-term alpha^r_01
contributions.

For each separable term (power + Gaussian: 54 terms), the alpha^r_01 factors as
  S_i(tau) * g_i(delta)
so we fit only g_i(delta) on delta in [0, delta_max] with adaptive subdivision
at bounded degree n_max. For non-analytic terms (2 terms), the tau and delta
parts are coupled -- we study both a per-tau 1D fit and a 2D tensor fit.

Report: pieces per term at n_max in {6, 8, 10, 12}, with tol 1e-12 tail.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from iapws95 import IAPWS95Terms
from compact_cheb import chebyshev_lobatto_nodes, values_to_coeffs


OUTDIR = os.path.dirname(os.path.abspath(__file__))
DELTA_MAX_DEFAULT = 5.0
TAU_MIN_DEFAULT = 0.3
TAU_MAX_DEFAULT = 3.0


def adaptive_subdivide(f, a, b, n_max, tol_rel=1e-12, tol_abs=1e-300,
                       min_width=1e-10, max_pieces=5000):
    stack = [(a, b)]
    pieces = []
    while stack:
        if len(pieces) + len(stack) > max_pieces:
            raise RuntimeError(f"Exceeded {max_pieces} pieces")
        aa, bb = stack.pop()
        nodes = chebyshev_lobatto_nodes(n_max, aa, bb)
        try:
            vals = np.array([f(x) for x in nodes])
        except Exception:
            vals = np.full(n_max + 1, np.nan)
        if not np.all(np.isfinite(vals)):
            if bb - aa < min_width:
                pieces.append((aa, bb, np.full(n_max + 1, np.nan)))
                continue
            mid = 0.5 * (aa + bb)
            stack.append((mid, bb))
            stack.append((aa, mid))
            continue
        coeffs = values_to_coeffs(vals)
        abs_c = np.abs(coeffs)
        cmax = abs_c.max() if abs_c.max() > 0 else 1.0
        tail_len = max(2, (n_max + 1) // 4)
        tail = abs_c[-tail_len:]
        if np.max(tail) < tol_rel * cmax + tol_abs:
            pieces.append((aa, bb, coeffs))
            continue
        if bb - aa < min_width:
            pieces.append((aa, bb, coeffs))
            continue
        mid = 0.5 * (aa + bb)
        stack.append((mid, bb))
        stack.append((aa, mid))
    pieces.sort(key=lambda p: p[0])
    return pieces


# ---- 1D kernel for separable terms: g_i(delta) ----
def power_delta_kernel(eos, idx):
    """Return g_i(delta) for the i-th power term of alpha^r_01."""
    d = eos.d_p[idx]; l = eos.l_p[idx]

    def g(delta):
        if l == 0:
            return (delta ** d) * d  # factor = d - l*delta^l = d
        return (delta ** d) * np.exp(-delta ** l) * (d - l * (delta ** l))
    return g


def gauss_delta_kernel(eos, idx):
    d = eos.d_g[idx]
    eta = eos.eta_g[idx]
    eps = eos.eps_g[idx]

    def g(delta):
        inner = d - 2 * eta * delta * (delta - eps)
        return (delta ** d) * np.exp(-eta * (delta - eps) ** 2) * inner
    return g


def nonan_delta_kernel_at_tau(eos, idx, tau):
    """Return g(delta) = delta * d(alpha^r_na_i)/d(delta) AT a given tau
    (non-separable, so there is no tau-independent kernel)."""
    a = eos.a_na[idx]; b = eos.b_na[idx]
    beta = eos.beta_na[idx]; A = eos.A_na[idx]; B = eos.B_na[idx]
    C = eos.C_na[idx]; D = eos.D_na[idx]
    n = eos.n_na[idx]

    def fn(delta):
        dm1 = delta - 1.0
        dm1sq = dm1 * dm1
        if dm1sq == 0:
            theta = 1.0 - tau
            Delta = theta ** 2
            psi = np.exp(-D * (tau - 1.0) ** 2)
            dtheta_dd = 0.0
            dDelta_dd = 0.0
        else:
            theta = (1.0 - tau) + A * np.power(dm1sq, 1.0 / (2.0 * beta))
            Delta = theta ** 2 + B * np.power(dm1sq, a)
            psi = np.exp(-C * dm1sq - D * (tau - 1.0) ** 2)
            dtheta_dd = (A / beta) * np.power(dm1sq, 1.0 / (2.0 * beta) - 1.0) * dm1
            dDelta_dd = 2 * theta * dtheta_dd + \
                2 * a * B * np.power(dm1sq, a - 1.0) * dm1
        dpsi_dd = -2 * C * dm1 * psi
        term1 = np.power(Delta, b) * (psi + delta * dpsi_dd)
        term2 = b * np.power(Delta, b - 1.0) * dDelta_dd * delta * psi \
            if dm1sq != 0 else 0.0
        dalphar_dd = n * (term1 + term2)
        return delta * dalphar_dd

    return fn


def study_separable(eos, n_max_values=(6, 8, 10, 12), delta_max=DELTA_MAX_DEFAULT):
    print(f"\n### 1D piecewise fits in delta on [0, {delta_max}], tol 1e-12 ###")
    print(f"{'term':<18}  " + "  ".join(f"n={n:<3d}" for n in n_max_values))
    out = []
    # Power terms
    for i in range(len(eos.n_p)):
        g = power_delta_kernel(eos, i)
        label = f"power_l{int(eos.l_p[i])}_i{i}"
        row = f"{label:<18}"
        counts = {}
        for n_max in n_max_values:
            try:
                pieces = adaptive_subdivide(g, 0.0, delta_max, n_max=n_max)
                row += f"  {len(pieces):>4d}  "
                counts[n_max] = len(pieces)
            except Exception:
                row += f"  ERR   "
                counts[n_max] = -1
        print(row)
        out.append((label, counts))
    for i in range(len(eos.n_g)):
        g = gauss_delta_kernel(eos, i)
        label = f"gauss_i{i}"
        row = f"{label:<18}"
        counts = {}
        for n_max in n_max_values:
            try:
                pieces = adaptive_subdivide(g, 0.0, delta_max, n_max=n_max)
                row += f"  {len(pieces):>4d}  "
                counts[n_max] = len(pieces)
            except Exception:
                row += f"  ERR   "
                counts[n_max] = -1
        print(row)
        out.append((label, counts))
    return out


def study_nonan_per_tau(eos, tau_values=(0.3, 0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0),
                         n_max_values=(6, 8, 10, 12, 16),
                         delta_max=DELTA_MAX_DEFAULT):
    print(f"\n### Non-analytic terms: 1D slices at fixed tau (delta in [0, {delta_max}]) ###")
    print(f"{'term':<12} {'tau':>6}  " + "  ".join(f"n={n:<3d}" for n in n_max_values))
    for i in range(len(eos.n_na)):
        for tau in tau_values:
            fn = nonan_delta_kernel_at_tau(eos, i, tau)
            row = f"nonan_i{i}     {tau:>6.2f}  "
            for n_max in n_max_values:
                try:
                    pieces = adaptive_subdivide(fn, 0.0, delta_max, n_max=n_max)
                    row += f"  {len(pieces):>4d}  "
                except RuntimeError:
                    row += f"  >5k  "
                except Exception:
                    row += f"  ERR   "
            print(row)


def study_nonan_2d(eos, tau_min=TAU_MIN_DEFAULT, tau_max=TAU_MAX_DEFAULT,
                   delta_max=DELTA_MAX_DEFAULT, n_max=12, tol_rel=1e-8):
    """Quick-and-dirty 2D adaptive subdivision: split the (tau, delta) box
    into 4 equal tiles recursively until each tile has coefficient tail
    below tolerance in both dimensions."""
    print(f"\n### Non-analytic 2D piecewise fits on [tau,delta]=[{tau_min},{tau_max}]x[0,{delta_max}] ###")
    print(f"Target tail: {tol_rel}, per-tile degree: {n_max}x{n_max}")

    def fit_tile(fn, ta, tb, da, db, nmax):
        # sample on (nmax+1)^2 grid of tensor Chebyshev-Lobatto nodes
        tnodes = chebyshev_lobatto_nodes(nmax, ta, tb)
        dnodes = chebyshev_lobatto_nodes(nmax, da, db)
        V = np.empty((nmax + 1, nmax + 1))
        for i, tau in enumerate(tnodes):
            for j, delta in enumerate(dnodes):
                V[i, j] = fn(tau, delta)
        if not np.all(np.isfinite(V)):
            return None
        # 2D Chebyshev coeffs via row-then-column DCT
        Ct = np.array([values_to_coeffs(V[i, :]) for i in range(nmax + 1)])
        C = np.array([values_to_coeffs(Ct[:, j]) for j in range(nmax + 1)]).T
        return V, C

    def tail_ok(C, tol):
        cmax = np.max(np.abs(C))
        if cmax == 0:
            return True
        nmax = C.shape[0] - 1
        tail_len = max(2, (nmax + 1) // 4)
        # last 1/4 rows or last 1/4 cols
        tail_rows = np.max(np.abs(C[-tail_len:, :]))
        tail_cols = np.max(np.abs(C[:, -tail_len:]))
        return max(tail_rows, tail_cols) < tol * cmax

    def subdivide_2d(fn, ta, tb, da, db, nmax, tol, depth=0, max_depth=8,
                    min_w=1e-6):
        if depth >= max_depth or (tb - ta) < min_w or (db - da) < min_w:
            return 1
        res = fit_tile(fn, ta, tb, da, db, nmax)
        if res is None:
            # bisect into 4
            tm = 0.5 * (ta + tb); dm = 0.5 * (da + db)
            return (subdivide_2d(fn, ta, tm, da, dm, nmax, tol, depth+1) +
                    subdivide_2d(fn, tm, tb, da, dm, nmax, tol, depth+1) +
                    subdivide_2d(fn, ta, tm, dm, db, nmax, tol, depth+1) +
                    subdivide_2d(fn, tm, tb, dm, db, nmax, tol, depth+1))
        V, C = res
        if tail_ok(C, tol):
            return 1
        tm = 0.5 * (ta + tb); dm = 0.5 * (da + db)
        return (subdivide_2d(fn, ta, tm, da, dm, nmax, tol, depth+1) +
                subdivide_2d(fn, tm, tb, da, dm, nmax, tol, depth+1) +
                subdivide_2d(fn, ta, tm, dm, db, nmax, tol, depth+1) +
                subdivide_2d(fn, tm, tb, dm, db, nmax, tol, depth+1))

    for i in range(len(eos.n_na)):
        def fn(tau, delta, _i=i):
            return nonan_delta_kernel_at_tau(eos, _i, tau)(delta)
        tiles = subdivide_2d(fn, tau_min, tau_max, 0.0, delta_max,
                              nmax=n_max, tol=tol_rel)
        print(f"  nonan_i{i}: {tiles} tiles at n={n_max}x{n_max}, tol={tol_rel}")


def main():
    eos = IAPWS95Terms()
    out = study_separable(eos)
    study_nonan_per_tau(eos)
    study_nonan_2d(eos, n_max=12, tol_rel=1e-8)
    study_nonan_2d(eos, n_max=16, tol_rel=1e-8)
    # totals for separable terms
    print("\n### Totals for separable 54 terms ###")
    for nmax in (6, 8, 10, 12):
        total = sum(d[nmax] for _, d in out if d.get(nmax, -1) > 0)
        print(f"  n={nmax}: {total} total 1D pieces across all 54 separable terms")


if __name__ == "__main__":
    main()
