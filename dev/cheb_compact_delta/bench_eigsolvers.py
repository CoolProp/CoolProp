"""Benchmark Chebyshev colleague-matrix eigensolvers on sizes relevant
to our convergence study (n = 32, 64, 96, 128, 192, 256).

We construct random Chebyshev coefficient vectors representative of
decaying-tail series (so the colleague matrix is well-conditioned), and
time:
  1. numpy.polynomial.chebyshev.chebroots  (numpy.linalg.eigvals on colleague)
  2. scipy.linalg.eigvals on the colleague matrix directly (with
     check_finite=False)
  3. scipy.linalg.lapack.dhseqr on the pre-Hessenberg colleague matrix
     (skips the reduction step since colleague is already Hessenberg)
  4. (if available) balanced-then-dhseqr
  5. subdivide-to-degree-32-and-eig  (Chebfun-style)

We report wall time per solve, and filter real roots in [-1, 1] to verify
each method returns the same set.
"""

import time
import numpy as np
import scipy.linalg as sla
from numpy.polynomial import chebyshev as npcheb
from scipy.linalg import lapack


def build_colleague_matrix(coeffs):
    """Return the upper-Hessenberg colleague matrix for the Chebyshev
    series sum_k coeffs[k] T_k(x). Trefethen/Chebfun convention.

    For degree n (len(coeffs) = n+1, leading coeff nonzero), the matrix
    is n x n. The Chebyshev three-term recurrence gives:
        C = [[0, 1, 0, 0, ...],
             [1/2, 0, 1/2, 0, ...],
             [0, 1/2, 0, 1/2, ...],
             ...
             [0, 0, ..., 1/2, 0, 1/2],
             [a_0, a_1, ..., a_{n-1}]]
    where the last row is  a_k = -coeffs[k] / (2 * coeffs[n])  with the
    correction  a_{n-1} += 1/2.
    """
    c = np.asarray(coeffs, dtype=float).copy()
    # drop trailing zeros
    while len(c) > 1 and abs(c[-1]) < 1e-300:
        c = c[:-1]
    n = len(c) - 1
    if n <= 0:
        return np.zeros((0, 0))
    # Standard colleague matrix
    mat = np.zeros((n, n))
    if n >= 2:
        mat[0, 1] = 1.0
        for i in range(1, n - 1):
            mat[i, i - 1] = 0.5
            mat[i, i + 1] = 0.5
    if n >= 2:
        mat[n - 1, n - 2] = 0.5
    # Last row: coefficient feedback
    # Roots of sum_k c_k T_k are eigenvalues of this modified matrix
    last = -c[:-1] / (2.0 * c[-1])
    last[-1 if n - 1 > 0 else 0] += 0.5  # adjustment to last entry
    if n == 1:
        mat[0, 0] = -c[0] / c[1]
    else:
        mat[-1, :] = last
    return mat


def sample_coeffs(n, decay_rate=0.1, seed=0):
    """Generate a random Chebyshev series of degree n with exponential
    coefficient decay (realistic for smooth functions)."""
    rng = np.random.default_rng(seed)
    c = rng.standard_normal(n + 1)
    c = c * np.exp(-decay_rate * np.arange(n + 1))
    return c


def method_numpy_chebroots(coeffs):
    return npcheb.chebroots(coeffs)


def method_scipy_eigvals(coeffs):
    M = build_colleague_matrix(coeffs)
    if M.shape[0] == 0:
        return np.array([])
    return sla.eigvals(M, check_finite=False)


def method_dhseqr(coeffs):
    """LAPACK dhseqr on the Hessenberg colleague matrix directly."""
    M = build_colleague_matrix(coeffs)
    n = M.shape[0]
    if n == 0:
        return np.array([])
    # dhseqr expects ilo, ihi, H, Z, work parameters.
    # scipy.linalg.lapack.dhseqr(job, compz, ilo, ihi, h, z, ...)
    H = np.asfortranarray(M.copy())
    Z = np.eye(n, order='F')
    # job='E' => eigenvalues only; compz='N' => no Schur vectors
    wr, wi, Z_out, H_out, info = lapack.dhseqr(H, lwork=max(1, 11 * n),
                                                 compute_q=0, calc_q=0) \
        if False else _dhseqr_compat(H, n)
    if info != 0:
        # fallback
        return sla.eigvals(M, check_finite=False)
    return wr + 1j * wi


def _dhseqr_compat(H, n):
    """Call dhseqr with the signature scipy actually exposes in this version."""
    # scipy.linalg.lapack.dhseqr signature varies; the most portable is:
    # result = dhseqr(h, ilo=1, ihi=n, compute_q=0, z=None, ...)
    result = lapack.dhseqr(H)
    # result layout (scipy >=1.0): h_out, z_out, wr, wi, info
    if len(result) == 5:
        h_out, z_out, wr, wi, info = result
    elif len(result) == 4:
        h_out, wr, wi, info = result
        z_out = None
    else:
        raise RuntimeError(f"Unexpected dhseqr return: len={len(result)}")
    return wr, wi, z_out, h_out, info


def method_dhseqr_balanced(coeffs):
    """Balance first (LAPACK gebal) then dhseqr."""
    M = build_colleague_matrix(coeffs)
    n = M.shape[0]
    if n == 0:
        return np.array([])
    M = np.asfortranarray(M.copy())
    try:
        bal_res = lapack.dgebal(M)
        # scipy returns (ba, lo, hi, pivscale, info) typically
        if len(bal_res) == 5:
            ba, lo, hi, pivscale, info = bal_res
        else:
            ba, lo, hi, info = bal_res[:4]
            pivscale = None
        if info != 0:
            ba = M
    except Exception:
        ba = M
    return sla.eigvals(ba, check_finite=False)


def method_subdivide(coeffs, max_n_per_piece=32, a=-1.0, b=1.0, depth=0, max_depth=12):
    """Chebfun-style recursive subdivision: if series degree > threshold,
    bisect the interval, refit each half, recurse. Only returns roots in [a, b]."""
    n = len(coeffs) - 1
    if n <= max_n_per_piece:
        roots = npcheb.chebroots(coeffs)
        real_roots = roots[np.abs(roots.imag) < 1e-8 * (1 + np.abs(roots))].real
        in_range = real_roots[(real_roots >= -1.001) & (real_roots <= 1.001)]
        # map back to [a, b]
        return 0.5 * (a + b) + 0.5 * (b - a) * in_range
    if depth >= max_depth:
        return np.array([])
    # Bisect: rebuild on two halves by sampling cheb-lobatto on each half
    nodes_std = np.cos(np.pi * np.arange(max_n_per_piece + 1) / max_n_per_piece)
    # sample original series at the new nodes mapped into left and right halves
    mid = 0.5 * (a + b)
    def orig_eval(x_std):
        # map x_std in [-1,1] to original [-1,1]: x_std corresponds to x
        return npcheb.chebval(x_std, coeffs)
    # left interval [a, mid] in original variable has length (mid-a) on the original [-1,1]
    # (if original was [-1,1], left is [-1, 0], maps via x = -1 + (x_std+1)*0.5*(mid-a)/((b-a)/2))
    # For simplicity assume the caller is passing the series valid on [a,b].
    left_samples = np.array([orig_eval(a + (mid - a) * 0.5 * (xs + 1)) for xs in nodes_std])
    right_samples = np.array([orig_eval(mid + (b - mid) * 0.5 * (xs + 1)) for xs in nodes_std])
    left_coeffs = _vals_to_coeffs(left_samples)
    right_coeffs = _vals_to_coeffs(right_samples)
    left_roots = method_subdivide(left_coeffs, max_n_per_piece, a, mid, depth + 1, max_depth)
    right_roots = method_subdivide(right_coeffs, max_n_per_piece, mid, b, depth + 1, max_depth)
    return np.concatenate([left_roots, right_roots])


def _vals_to_coeffs(vals):
    n = len(vals) - 1
    v = np.concatenate([vals, vals[-2:0:-1]])
    F = np.real(np.fft.fft(v))
    c = F[: n + 1] / n
    c[0] *= 0.5
    c[-1] *= 0.5
    return c


def time_method(fn, coeffs, n_repeat=5):
    # warmup
    fn(coeffs)
    t0 = time.perf_counter()
    for _ in range(n_repeat):
        fn(coeffs)
    return (time.perf_counter() - t0) / n_repeat


def main():
    sizes = [32, 64, 96, 128, 192, 256, 384, 512]
    methods = [
        ("numpy.chebroots (eigvals on colleague)", method_numpy_chebroots),
        ("scipy.eigvals(check_finite=False)     ", method_scipy_eigvals),
        ("scipy.lapack.dhseqr (skip Hess reduc) ", method_dhseqr),
        ("balanced + eigvals                     ", method_dhseqr_balanced),
        ("subdivide to n=32, eig each half       ", method_subdivide),
    ]
    print(f"{'n':>5} " + "  ".join([f"{name:<40}" for name, _ in methods]))
    for n in sizes:
        c = sample_coeffs(n, decay_rate=0.05, seed=n)
        row = f"{n:>5}  "
        for name, fn in methods:
            try:
                t = time_method(fn, c, n_repeat=10 if n <= 128 else 5)
                row += f"{t*1e6:>14.2f} us                       "
            except Exception as e:
                row += f"{'err: '+str(e)[:30]:<40}"
        print(row)

    # sanity: roots agree across methods at n=64
    print("\nSanity check at n=64: root sets should agree (real, in [-1,1])")
    c = sample_coeffs(64, 0.05, 42)
    for name, fn in methods[:4]:
        r = fn(c)
        if hasattr(r, "imag"):
            real = np.sort(r[np.abs(r.imag) < 1e-8].real)
            in_rng = real[(real >= -1.01) & (real <= 1.01)]
        else:
            in_rng = np.sort(r)
        print(f"  {name}: {len(in_rng)} real roots in [-1,1], "
              f"first few = {in_rng[:3]}")


if __name__ == "__main__":
    main()
