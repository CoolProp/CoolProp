"""
Piecewise Chebyshev rootfinder on a compactified delta domain for
pressure inversion of multi-parameter Helmholtz EOS.

Change of variable:  u = (delta - L) / (delta + L),  L = 1  (delta already reduced by rho_c)
  delta = 0 <-> u = -1     (vacuum)
  delta = 1 <-> u =  0     (critical density)
  delta -> inf <-> u -> +1

A single Chebyshev expansion on u in [-1+eps, 1-eps] covers the full
density range. Subdivision is adaptive (Chebfun-style tail-decay test).
Per-subinterval rootfinding uses the Chebyshev colleague matrix.
"""

import numpy as np
from numpy.polynomial import chebyshev as npcheb


# ---------------------------------------------------------------
# Compactification
# ---------------------------------------------------------------
def delta_of_u(u, L=1.0):
    return L * (1.0 + u) / (1.0 - u)


def u_of_delta(delta, L=1.0):
    return (delta - L) / (delta + L)


def ddelta_du(u, L=1.0):
    return 2.0 * L / (1.0 - u) ** 2


# ---------------------------------------------------------------
# Chebyshev basics (second-kind nodes, DCT-I coefficient extraction)
# ---------------------------------------------------------------
def chebyshev_lobatto_nodes(n, a=-1.0, b=1.0):
    """n+1 Chebyshev-Lobatto points on [a, b]."""
    k = np.arange(n + 1)
    x = np.cos(np.pi * k / n)
    return 0.5 * (a + b) + 0.5 * (b - a) * x


def values_to_coeffs(vals):
    """DCT-I: given function values at n+1 Chebyshev-Lobatto nodes on [-1,1]
    (ordered with x_0=1, x_n=-1), return Chebyshev series coefficients c_0..c_n
    such that f(x) ~ sum_k c_k T_k(x)."""
    n = len(vals) - 1
    v = np.concatenate([vals, vals[-2:0:-1]])
    F = np.real(np.fft.fft(v))
    c = F[: n + 1] / n
    c[0] *= 0.5
    c[-1] *= 0.5
    return c


def chebyshev_eval(coeffs, x):
    return npcheb.chebval(x, coeffs)


# ---------------------------------------------------------------
# Adaptive Chebyshev fit on one subinterval [a, b] in u
# ---------------------------------------------------------------
def fit_chebyshev_adaptive(f, a, b, tol=1e-10, n_init=16, n_max=256):
    """Try degrees n_init, 2*n_init, ..., n_max on [a, b] until the tail
    Chebyshev coefficients fall below tol * max(|c|). Return (coeffs, degree)
    or (None, -1) if no degree succeeded."""
    n = n_init
    cached_vals = {}

    def evaluate_at_nodes(nn):
        # reuse cached f-evaluations at the subset of Lobatto points that coincide
        nodes = chebyshev_lobatto_nodes(nn, a, b)
        vals = np.empty(nn + 1)
        for i, x in enumerate(nodes):
            key = (nn, i)
            if x in cached_vals:
                vals[i] = cached_vals[x]
            else:
                vals[i] = f(x)
                cached_vals[x] = vals[i]
        return vals

    while n <= n_max:
        vals = evaluate_at_nodes(n)
        if not np.all(np.isfinite(vals)):
            return None, -1
        coeffs = values_to_coeffs(vals)
        abs_c = np.abs(coeffs)
        cmax = abs_c.max() if abs_c.max() > 0 else 1.0
        # Chebfun-style "plateau" in the tail: last ~min(8, n//4) coeffs
        tail_len = max(4, n // 8)
        tail = abs_c[-tail_len:]
        if np.max(tail) < tol * cmax + tol * 1e-6:
            return coeffs, n
        n *= 2
    return None, -1


# ---------------------------------------------------------------
# Colleague matrix rootfinding for a Chebyshev series on [a, b]
# ---------------------------------------------------------------
def cheb_roots_on_interval(coeffs, a, b, drop_tail_tol=None):
    """Return all real roots of sum_k coeffs[k] T_k(x) on [a, b]."""
    if drop_tail_tol is not None:
        abs_c = np.abs(coeffs)
        cmax = abs_c.max() if len(abs_c) > 0 else 0.0
        if cmax > 0:
            keep = np.where(abs_c > drop_tail_tol * cmax)[0]
            if len(keep) > 0:
                coeffs = coeffs[: keep[-1] + 1]
    if len(coeffs) < 2:
        return np.array([])
    # numpy's chebroots uses the colleague matrix internally
    roots = npcheb.chebroots(coeffs)
    # keep real roots in [-1, 1] with small imag tolerance
    imag_tol = 1e-8
    real_roots = roots[np.abs(roots.imag) < imag_tol * (1 + np.abs(roots))].real
    # [-1, 1] with a small cushion for floating-point
    in_range = real_roots[(real_roots >= -1.0 - 1e-12) & (real_roots <= 1.0 + 1e-12)]
    in_range = np.clip(in_range, -1.0, 1.0)
    # map [-1, 1] back to [a, b]
    return 0.5 * (a + b) + 0.5 * (b - a) * in_range


# ---------------------------------------------------------------
# Piecewise adaptive: repeatedly bisect until each piece is happy
# ---------------------------------------------------------------
def build_piecewise(f, a, b, tol=1e-10, n_init=16, n_max=256,
                    min_width=1e-6, max_pieces=4096):
    """Return list of (a_i, b_i, coeffs_i, converged_i) pieces.
    converged_i is True if the Chebyshev tail-decay criterion was met; False
    if the piece bottomed out at min_width without happiness (in that case
    coeffs_i is still populated via the max-degree fit, but its roots should
    be treated as candidates only and validated against the raw function)."""
    stack = [(a, b)]
    pieces = []
    while stack:
        if len(pieces) + len(stack) > max_pieces:
            raise RuntimeError(f"Exceeded {max_pieces} pieces")
        aa, bb = stack.pop()
        coeffs, deg = fit_chebyshev_adaptive(
            lambda x: f(x), aa, bb, tol=tol, n_init=n_init, n_max=n_max
        )
        if coeffs is not None:
            pieces.append((aa, bb, coeffs, True))
            continue
        if bb - aa < min_width:
            nodes = chebyshev_lobatto_nodes(n_max, aa, bb)
            vals = np.array([f(x) for x in nodes])
            if np.all(np.isfinite(vals)):
                pieces.append((aa, bb, values_to_coeffs(vals), False))
            continue
        mid = 0.5 * (aa + bb)
        stack.append((mid, bb))
        stack.append((aa, mid))
    pieces.sort(key=lambda p: p[0])
    return pieces


# ---------------------------------------------------------------
# Top-level: all roots of f on [a, b]
# ---------------------------------------------------------------
def all_roots_piecewise(f, a, b, tol=1e-10, n_init=16, n_max=256,
                        min_width=1e-6, drop_tail_tol=1e-14,
                        validate_abstol=None, f_scale=1.0):
    """Return (roots, pieces). If validate_abstol is given, each candidate
    root is kept only if |f(root)| < validate_abstol (in whatever units
    f produces). Pieces that fail the happiness check are handled by
    scanning their Chebyshev nodes for sign changes and Brent-solving within
    each bracket -- this is robust to EOS pathology where the piecewise
    polynomial fit itself is unreliable."""
    from scipy.optimize import brentq

    pieces = build_piecewise(f, a, b, tol=tol, n_init=n_init,
                             n_max=n_max, min_width=min_width)
    roots = []
    for (aa, bb, coeffs, converged) in pieces:
        if converged:
            rr = cheb_roots_on_interval(coeffs, aa, bb, drop_tail_tol=drop_tail_tol)
            for r in rr:
                if aa - 1e-14 <= r <= bb + 1e-14:
                    if validate_abstol is None or abs(f(r)) < validate_abstol:
                        roots.append(r)
        else:
            # Unreliable piece: scan for sign changes on a fine grid of the raw
            # function and Brent-solve inside each bracket.
            n_scan = 257
            xs = np.linspace(aa, bb, n_scan)
            fs = np.array([f(x) for x in xs])
            for i in range(n_scan - 1):
                if np.isfinite(fs[i]) and np.isfinite(fs[i + 1]) and \
                        fs[i] * fs[i + 1] < 0:
                    try:
                        r = brentq(f, xs[i], xs[i + 1], xtol=1e-12)
                        if validate_abstol is None or abs(f(r)) < validate_abstol:
                            roots.append(r)
                    except Exception:
                        pass
    roots = np.array(sorted(roots))
    if len(roots) > 1:
        keep = [roots[0]]
        for r in roots[1:]:
            if r - keep[-1] > 1e-10:
                keep.append(r)
        roots = np.array(keep)
    return roots, pieces
