"""Piecewise Chebyshev density solver with small per-interval degree.

Architecture (following paper):
  - Offline: partition [delta_min, delta_max] into K subintervals; on each,
    per-term H̃_k(δ) Chebyshev coefficients of degree n_δ (small, e.g. 12).
    Subdivision: dyadic, driven by the ε_split,M metric (paper eq. 54) applied
    to the worst-case across all terms so the subdivision is shared.
  - Runtime: given (T, P_target), compute F_k(τ) analytically per term → scalar
    vector. Per-interval: weighted sum of per-term H_coeffs → local δ-Cheb.
    Build pressure residual r(δ) = ρ_r·R·T·δ·(1 + A_01) − P_target. Colleague
    matrix eig at small n.
  - Fast path: skip intervals whose ||coeffs||_1 bound shows the polynomial
    cannot cross zero within the interval (Chebyshev values are bounded by
    the sum of absolute coefficients).
"""

import numpy as np
from numpy.polynomial import chebyshev as npcheb

from paper_architecture import (
    chebyshev_lobatto_nodes, values_to_coeffs, eval_cheb_on,
    times_x_on_interval, add_coeffs,
)


def fit_cheb_on(f, n, a, b):
    nodes = chebyshev_lobatto_nodes(n, a, b)
    vals = np.array([f(x) for x in nodes])
    return values_to_coeffs(vals)


def epsilon_split_M(coeffs, M=3):
    """ε_split,M = sqrt(Σ c_{N-M+1..N}²) / sqrt(Σ c_0..M²). Paper eq. 54."""
    n = len(coeffs) - 1
    if n < 2 * M:
        return 0.0
    tail = coeffs[n - M + 1:]
    head = coeffs[: M + 1]
    num = np.sqrt(np.sum(tail * tail))
    den = np.sqrt(np.sum(head * head)) + 1e-300
    return num / den


def build_subdivision(terms, delta_min, delta_max, n_delta=12, tol=1e-10,
                       M=3, min_width=1e-8, max_intervals=4096):
    """Dyadic subdivision of [delta_min, delta_max] targeting per-interval
    ε_split,M < tol across all terms. Returns a list of (a, b, H_coeffs_matrix)
    where H_coeffs_matrix is (n_terms, n_delta+1)."""
    def fit_all_terms(a, b):
        H = np.zeros((len(terms), n_delta + 1))
        for k, term in enumerate(terms):
            H[k, :] = fit_cheb_on(term["H"], n_delta, a, b)
        return H

    def worst_eps_split(H):
        return max(epsilon_split_M(H[k, :], M) for k in range(H.shape[0]))

    pieces = []
    stack = [(delta_min, delta_max)]
    while stack:
        if len(pieces) + len(stack) > max_intervals:
            raise RuntimeError(f">{max_intervals} intervals")
        a, b = stack.pop()
        H = fit_all_terms(a, b)
        e = worst_eps_split(H)
        if e < tol or (b - a) < min_width:
            pieces.append((a, b, H))
            continue
        mid = 0.5 * (a + b)
        stack.append((mid, b))
        stack.append((a, mid))
    pieces.sort(key=lambda p: p[0])
    return pieces


def build_fluid_fast(fluid, delta_min, delta_max, n_delta=12, tol=1e-10):
    terms = fluid.separable_terms()
    n_values = np.array([t["n"] for t in terms])
    pieces = build_subdivision(terms, delta_min, delta_max,
                                n_delta=n_delta, tol=tol)
    # Vectorized layout: stack all intervals' H matrices together
    K = len(pieces)
    n_terms = len(terms)
    intervals_a = np.array([p[0] for p in pieces])
    intervals_b = np.array([p[1] for p in pieces])
    H_all = np.stack([p[2] for p in pieces], axis=0)  # (K, n_terms, n_delta+1)
    # Pre-flatten for BLAS matmul: H_flat has shape (n_terms, K*(n_delta+1)).
    # At runtime, weights @ H_flat is a single matrix-vector op → (K*(n_delta+1),) that
    # we reshape to (K, n_delta+1). Much faster than np.einsum in numpy.
    H_flat = np.ascontiguousarray(H_all.transpose(1, 0, 2).reshape(n_terms, -1))
    return dict(terms=terms, n_values=n_values, pieces=pieces,
                intervals_a=intervals_a, intervals_b=intervals_b,
                H_all=H_all, H_flat=H_flat, K=K, n_terms=n_terms,
                delta_min=delta_min, delta_max=delta_max, n_delta=n_delta,
                fluid=fluid)


def cheb_derivative_coeffs(c):
    """Given Chebyshev coefficients c of degree n (length n+1) on [-1, 1],
    return coefficients d of f'(x) of degree n-1 (length n). Vectorized over
    leading axis: c may be shape (..., n+1); output is (..., n)."""
    # Recurrence: d_{n-1} = 2n c_n; d_{n-2} = 2(n-1) c_{n-1};
    #             d_i = d_{i+2} + 2(i+1) c_{i+1} for i = n-3 .. 0
    #             then d_0 /= 2
    n = c.shape[-1] - 1
    if n < 1:
        return np.zeros_like(c[..., :0])
    d = np.zeros(c.shape[:-1] + (n,))
    d[..., n - 1] = 2 * n * c[..., n]
    if n >= 2:
        d[..., n - 2] = 2 * (n - 1) * c[..., n - 1]
    for i in range(n - 3, -1, -1):
        d[..., i] = d[..., i + 2] + 2 * (i + 1) * c[..., i + 1]
    d[..., 0] *= 0.5
    return d


def density_roots_vectorized(state, T, P_target):
    """Vectorized density rootfinder: builds all intervals' pressure-residual
    Chebyshev series at once, applies cheap skip test, then runs colleague-matrix
    rootfinding on the surviving intervals."""
    fluid = state["fluid"]
    tau = fluid.Tc / T
    # Per-term F_k(τ) scalars (vectorized analytic evaluation)
    F_vals = np.array([t["F"](tau) for t in state["terms"]])  # (n_terms,)
    weights = state["n_values"] * F_vals  # (n_terms,)

    # A01_all[k, j] = Σ_t weights[t] * H_all[k, t, j]  — shape (K, n_delta+1)
    # BLAS matvec on pre-flattened H is typically faster than einsum in numpy.
    A01_flat = weights @ state["H_flat"]  # shape (K*(n_delta+1),)
    A01_all = A01_flat.reshape(state["K"], state["n_delta"] + 1)

    # times_x on each interval: vectorized implementation of the multiply-by-x
    # identity (paper eq. 29) + interval rescaling.
    K = state["K"]
    n_delta = state["n_delta"]
    halfs = 0.5 * (state["intervals_b"] - state["intervals_a"])  # (K,)
    mids = 0.5 * (state["intervals_a"] + state["intervals_b"])    # (K,)
    # Multiply-by-x_std identity on each row A01_all[k, :]:
    #   c'_0 = 0.5 a_1
    #   c'_1 = a_0 + 0.5 a_2
    #   c'_i = 0.5 (a_{i-1} + a_{i+1}),   2 <= i <= n-1
    #   c'_n = 0.5 a_{n-1}
    #   c'_{n+1} = 0.5 a_n
    n = n_delta  # so A01_all has (n+1) cols
    c_prime = np.zeros((K, n + 2))
    if n >= 1:
        c_prime[:, 0] = 0.5 * A01_all[:, 1]
    if n >= 2:
        c_prime[:, 1] = A01_all[:, 0] + 0.5 * A01_all[:, 2]
    elif n == 1:
        c_prime[:, 1] = A01_all[:, 0]
    if n >= 3:
        c_prime[:, 2:n] = 0.5 * (A01_all[:, 1:n - 1] + A01_all[:, 3:n + 1])
    if n >= 1:
        c_prime[:, n] = 0.5 * A01_all[:, n - 1]
        c_prime[:, n + 1] = 0.5 * A01_all[:, n]

    # δ·A01 on each interval [a, b]:
    #   result = half[k] * c_prime + mid[k] * A01 (padded to n+2)
    dA01_all = halfs[:, None] * c_prime
    A01_padded = np.zeros((K, n + 2))
    A01_padded[:, : n + 1] = A01_all
    dA01_all[:, : n + 1] += mids[:, None] * A01_all
    # Now: ρ_r·R·T·(δ + δ·A01) − P_target
    factor = fluid.rhoc * fluid.R * T
    # δ on [a, b] is 2 coefficients: [mid, half]
    resid = factor * dA01_all  # (K, n+2)
    resid[:, 0] += factor * mids
    resid[:, 1] += factor * halfs
    resid[:, 0] -= P_target

    # Cheap "no root" test per interval: |c_0| > Σ |c_{k>=1}|
    # Proof: on [-1,1], |p(x_std)| ≤ Σ|c_k|, so if |c_0| > Σ|c_{k>=1}|, p has
    # constant sign on the interval and no zero.
    abs_res = np.abs(resid)
    no_root_mask = abs_res[:, 0] > np.sum(abs_res[:, 1:], axis=1)
    work_idxs = np.where(~no_root_mask)[0]

    rho_roots = []
    n_mono_root = 0
    n_eig = 0

    # Early exit if no intervals have potential roots
    if len(work_idxs) == 0:
        stats = dict(K=state["K"], no_root=state["K"],
                     mono_root=0, eig=0)
        return np.array([]), stats

    # Only compute derivative + monotonicity on non-skipped intervals (few).
    resid_work = resid[work_idxs, :]
    D_work = cheb_derivative_coeffs(resid_work)  # (|work|, n+1)
    abs_d = np.abs(D_work)
    mono_mask = abs_d[:, 0] > np.sum(abs_d[:, 1:], axis=1)

    # Endpoint values via alternating-sum / sum
    alternating = (-1.0) ** np.arange(resid.shape[1])
    f_at_a = np.sum(resid_work * alternating, axis=1)
    f_at_b = np.sum(resid_work, axis=1)
    sign_change = f_at_a * f_at_b < 0

    # Classify
    mono_with_root = mono_mask & sign_change
    needs_eig = ~mono_mask  # non-monotone → might have 0, 2, 3, ... roots

    # Monotone-with-root: inlined secant-with-bisection fallback.
    # Convergence criteria:
    #   (a) bracket width (b_s - a_s) < eps_bracket (1e-13 of [-1,1])
    #   (b) |fc| < eps_f * max(|fa|, |fb|, 1)  (relative f tolerance)
    #   (c) secant step lands on or outside the bracket (means |fc| is small
    #       in a scale-aware sense): accept c_new clamped to the bracket.
    eps_bracket = 1e-13
    eps_f = 1e-12
    for i in np.where(mono_with_root)[0]:
        k = work_idxs[i]
        c_k = resid[k, :]
        half = halfs[k]; mid = mids[k]
        a_s, b_s = -1.0, 1.0
        fa, fb = f_at_a[i], f_at_b[i]
        root_std = None
        for _ in range(12):
            c_new = b_s - fb * (b_s - a_s) / (fb - fa)
            # Secant lands on or outside the bracket → we've converged at
            # the endpoint where |f| is smaller; use that endpoint directly.
            if c_new <= a_s:
                root_std = a_s
                break
            if c_new >= b_s:
                root_std = b_s
                break
            fc = npcheb.chebval(c_new, c_k)
            if abs(fc) < eps_f * max(abs(fa), abs(fb), 1.0):
                root_std = c_new
                break
            if (b_s - a_s) < eps_bracket:
                root_std = c_new
                break
            if fa * fc < 0:
                b_s = c_new; fb = fc
            else:
                a_s = c_new; fa = fc
        if root_std is None:
            root_std = c_new
        delta_r = mid + half * root_std
        rho_roots.append(delta_r * fluid.rhoc)
        n_mono_root += 1

    # Non-monotone intervals → eigenvalue colleague matrix
    for i in np.where(needs_eig)[0]:
        k = work_idxs[i]
        coeffs = resid[k, :]
        roots_std = npcheb.chebroots(coeffs)
        real = roots_std[np.abs(roots_std.imag) < 1e-8].real
        in_range = real[(real >= -1.0 - 1e-10) & (real <= 1.0 + 1e-10)]
        if len(in_range) == 0:
            continue
        in_range = np.clip(in_range, -1.0, 1.0)
        delta_r = mids[k] + halfs[k] * in_range
        rho_roots.extend((delta_r * fluid.rhoc).tolist())
        n_eig += 1

    stats = dict(
        K=state["K"],
        no_root=int(no_root_mask.sum()),
        mono_root=n_mono_root,
        eig=n_eig,
    )
    return np.sort(np.asarray(rho_roots)), stats


def pressure_residual_coeffs_interval(state, tau, T, P_target, piece):
    """Return the Chebyshev coefficients on [a, b] of r(δ) = P(δ) - P_target."""
    a, b, H = piece
    # Evaluate F_k(tau) analytically per term
    F_vals = np.array([t["F"](tau) for t in state["terms"]])
    weights = state["n_values"] * F_vals
    A01_coeffs = weights @ H  # shape (n_delta+1,)
    # delta * A01 on [a, b]
    dA01 = times_x_on_interval(A01_coeffs, a, b)
    # delta on [a, b]
    half = 0.5 * (b - a); mid = 0.5 * (a + b)
    delta_cc = np.zeros(2)
    delta_cc[0] = mid; delta_cc[1] = half
    # pressure-prefactor
    factor = state["fluid"].rhoc * state["fluid"].R * T
    resid = add_coeffs(factor * delta_cc, factor * dA01)
    resid = add_coeffs(resid, np.array([-P_target]))
    return resid


def cheap_skip_test(coeffs):
    """Can this Chebyshev expansion possibly have a root in [-1, 1]?
    |p(x)| <= Σ |c_i| on [-1, 1]. If |c_0| > Σ_{i>=1} |c_i|, the polynomial
    has constant sign; no root in interval."""
    if len(coeffs) == 0:
        return True
    c0 = coeffs[0]
    tail = np.sum(np.abs(coeffs[1:]))
    return abs(c0) > tail


def density_roots_fast(state, T, P_target, report_spurious=False):
    """Find all real roots of P(delta) - P_target = 0 across all intervals.
    Returns rho roots sorted."""
    tau = state["fluid"].Tc / T
    rho_roots = []
    n_skip = 0; n_eig = 0
    for piece in state["pieces"]:
        a, b, _ = piece
        resid = pressure_residual_coeffs_interval(state, tau, T, P_target, piece)
        if cheap_skip_test(resid):
            n_skip += 1
            continue
        n_eig += 1
        roots_std = npcheb.chebroots(resid)
        real = roots_std[np.abs(roots_std.imag) < 1e-8].real
        in_range = real[(real >= -1.0 - 1e-10) & (real <= 1.0 + 1e-10)]
        if len(in_range) == 0:
            continue
        in_range = np.clip(in_range, -1.0, 1.0)
        half = 0.5 * (b - a); mid = 0.5 * (a + b)
        delta_piece = mid + half * in_range
        rho_piece = delta_piece * state["fluid"].rhoc
        rho_roots.extend(rho_piece.tolist())
    rho_roots = np.sort(np.asarray(rho_roots))
    if report_spurious:
        return rho_roots, (n_skip, n_eig)
    return rho_roots


if __name__ == "__main__":
    import time
    from fluid_from_json import FluidFromJSON
    import CoolProp.CoolProp as CP

    path = "/Users/ianbell/miniforge3/lib/python3.9/site-packages/teqp/fluiddata/dev/fluids/Nitrogen.json"
    fluid = FluidFromJSON(path)
    print(f"Fluid: {fluid.name} Tc={fluid.Tc} rhoc={fluid.rhoc:.4f}")

    for (n_delta, tol) in [(6, 1e-6), (8, 1e-8), (12, 1e-10), (16, 1e-12)]:
        t0 = time.perf_counter()
        state = build_fluid_fast(fluid, 1e-10, 5.0, n_delta=n_delta, tol=tol)
        dt_build = time.perf_counter() - t0
        K = len(state["pieces"])
        print(f"\nn_delta={n_delta}, tol={tol}: K={K} intervals, build={dt_build*1e3:.1f} ms")

        M = fluid.M
        cases = [
            (100.0, 1e5, "T=100K P=0.1MPa"),
            (126.0, 3.4e6, "T=126K (near Tc)"),
            (200.0, 5e7, "T=200K P=50MPa"),
            (80.0, 1e5, "T=80K (subcrit vapor)"),
            (80.0, 2e6, "T=80K (subcrit liquid)"),
        ]
        for T, P, label in cases:
            # warm-up
            density_roots_vectorized(state, T, P)
            t0 = time.perf_counter()
            N = 200
            for _ in range(N):
                roots, stats = density_roots_vectorized(state, T, P)
            dt = (time.perf_counter() - t0) / N
            try:
                rho_ref = CP.PropsSI('D', 'T', float(T), 'P', float(P), 'Nitrogen') / M
            except Exception:
                rho_ref = None
            print(f"  {label}: roots={roots}  cp_ref={rho_ref}  "
                  f"time={dt*1e6:.1f} us, K={stats['K']} "
                  f"no_root={stats['no_root']} mono_root={stats['mono_root']} "
                  f"eig={stats['eig']}")
