"""Implementation of Bell & Alpert 2018 Chebyshev density rootfinding,
restricted to IAPWS-95 separable terms (skipping the 2 non-analytic terms).

Per-term Chebyshev expansions are built:
    - F̃_k(τ) of degree n_tau on [tau_min, tau_max]
    - H̃_k(δ) of degree n_delta on [delta_min, delta_max]
      where H_k(δ) = δ · dG_k/dδ is the pressure-contribution kernel
                    of the k-th term (eq. 47 of paper, integrated in)

At call time, given τ and x:
    Ã_01(δ) = Σ_k n_k F̃_k(τ) · H̃_k(δ)
This is one Chebyshev series in δ. Pressure residual r(δ) = ρRT δ (1 + Ã_01) − p_target
is another Chebyshev series; eigenvalue-rootfind via colleague matrix.
"""

import numpy as np
from numpy.polynomial import chebyshev as npcheb
from iapws95 import IAPWS95Terms


# ---------------------------------------------------------------
# Chebyshev utilities
# ---------------------------------------------------------------
def chebyshev_lobatto_nodes(n, a=-1.0, b=1.0):
    k = np.arange(n + 1)
    x = np.cos(np.pi * k / n)
    return 0.5 * (a + b) + 0.5 * (b - a) * x


def values_to_coeffs(vals):
    """DCT-I: n+1 values at Chebyshev-Lobatto nodes → n+1 Chebyshev coefficients."""
    n = len(vals) - 1
    v = np.concatenate([vals, vals[-2:0:-1]])
    F = np.real(np.fft.fft(v))
    c = F[: n + 1] / n
    c[0] *= 0.5
    c[-1] *= 0.5
    return c


def fit_cheb_on(f, n, a, b):
    """Fit f on [a, b] with degree n, return coefficients (length n+1)."""
    nodes = chebyshev_lobatto_nodes(n, a, b)
    vals = np.array([f(x) for x in nodes])
    return values_to_coeffs(vals)


def eval_cheb_on(coeffs, x, a, b):
    """Evaluate Σ c_k T_k((2x-a-b)/(b-a))."""
    x_std = (2.0 * np.asarray(x) - a - b) / (b - a)
    return npcheb.chebval(x_std, coeffs)


def times_x_on_interval(coeffs, a, b):
    """Given Chebyshev coefficients c of f(x) on [a, b], return coefficients of
    g(x) = x · f(x) on [a, b]. Uses the identity x = (b-a)/2 · x_std + (a+b)/2
    and paper eq. (29) for multiplication by x_std in the standard [-1, 1]."""
    half = 0.5 * (b - a)
    mid = 0.5 * (a + b)
    n = len(coeffs) - 1
    # First: multiply by x_std on [-1, 1]
    # Paper eq. 29:
    #   c'_0 = 0.5 * a_1
    #   c'_1 = a_0 + 0.5 * a_2
    #   c'_i = 0.5 * (a_{i-1} + a_{i+1}),   2 <= i <= n-1
    #   c'_n = 0.5 * a_{n-1}
    #   c'_{n+1} = 0.5 * a_n
    a_ = np.asarray(coeffs, dtype=float)
    c_prime = np.zeros(n + 2)
    c_prime[0] = 0.5 * a_[1] if n >= 1 else 0.0
    if n >= 2:
        c_prime[1] = a_[0] + 0.5 * a_[2]
    elif n == 1:
        c_prime[1] = a_[0]
    else:
        c_prime[1] = 0.0
    for i in range(2, n):
        c_prime[i] = 0.5 * (a_[i - 1] + a_[i + 1])
    if n >= 1:
        c_prime[n] = 0.5 * a_[n - 1]
        c_prime[n + 1] = 0.5 * a_[n]
    # Now x_real = half * x_std + mid, so: (half * x_std + mid) * f = half * (x_std * f) + mid * f
    result = np.zeros(max(n + 2, len(coeffs)))
    result[: n + 2] += half * c_prime
    result[: n + 1] += mid * a_
    return result


def add_coeffs(a, b):
    """Pad and add two Chebyshev coefficient vectors."""
    n = max(len(a), len(b))
    out = np.zeros(n)
    out[: len(a)] += a
    out[: len(b)] += b
    return out


# ---------------------------------------------------------------
# Per-term kernels: F_k(τ) and H_k(δ) = δ · dG_k/dδ
# ---------------------------------------------------------------
def iapws95_separable_terms(eos):
    """Return a list of dicts, one per separable term, with:
        - F: callable τ → F_k(τ)
        - H: callable δ → δ · dG_k/dδ
        - n: the coefficient n_k multiplying the product in α^r
    The list concatenates the 7 pure-polynomial, 44 poly-exp, and 3 Gaussian-bell terms.
    Non-analytic terms are skipped (to be added in a later layer).
    """
    terms = []
    # Power: F_k = τ^t_k, H_k = δ^d · exp(-δ^l) · (d - l·δ^l); for l=0, H = δ^d · d.
    for i in range(len(eos.n_p)):
        n = float(eos.n_p[i]); d = float(eos.d_p[i])
        t = float(eos.t_p[i]); l = float(eos.l_p[i])

        def Fk(tau, _t=t):
            return tau ** _t

        def Hk(delta, _d=d, _l=l):
            if _l == 0:
                return (delta ** _d) * _d
            return (delta ** _d) * np.exp(-delta ** _l) * (_d - _l * delta ** _l)

        terms.append(dict(F=Fk, H=Hk, n=n, family="power", index=i,
                          meta=dict(d=d, t=t, l=l)))

    # Gaussian bell: F_k = τ^t exp(-β(τ-γ)^2); H_k = δ^d exp(-η(δ-ε)^2) (d - 2η·δ(δ-ε))
    for i in range(len(eos.n_g)):
        n = float(eos.n_g[i]); d = float(eos.d_g[i])
        t = float(eos.t_g[i])
        eta = float(eos.eta_g[i]); eps = float(eos.eps_g[i])
        beta = float(eos.beta_g[i]); gamma = float(eos.gamma_g[i])

        def Fk(tau, _t=t, _b=beta, _g=gamma):
            return (tau ** _t) * np.exp(-_b * (tau - _g) ** 2)

        def Hk(delta, _d=d, _e=eta, _p=eps):
            return (delta ** _d) * np.exp(-_e * (delta - _p) ** 2) \
                * (_d - 2 * _e * delta * (delta - _p))

        terms.append(dict(F=Fk, H=Hk, n=n, family="gauss", index=i,
                          meta=dict(d=d, t=t, eta=eta, eps=eps,
                                    beta=beta, gamma=gamma)))

    return terms


# ---------------------------------------------------------------
# Per-term Chebyshev construction
# ---------------------------------------------------------------
class FluidCheb:
    """Per-fluid collection of F̃_k(τ) and H̃_k(δ) Chebyshev coefficient vectors."""

    def __init__(self, eos, tau_range, delta_range, n_delta):
        """F_k(τ) is evaluated analytically at call time — no τ-Chebyshev fit.
        High-t IAPWS-95 terms (t up to 50) have dynamic range on [0.2, 4] that
        exceeds 1e60, making Chebyshev fits of F_k numerically unusable without
        dyadic subdivision. Since τ is a scalar at call time, skipping the
        τ-Chebyshev layer and evaluating F_k analytically is both faster and
        more accurate."""
        self.eos = eos
        self.tau_min, self.tau_max = tau_range
        self.delta_min, self.delta_max = delta_range
        self.n_delta = n_delta
        self.terms = iapws95_separable_terms(eos)

        # Only H_k(δ) gets Chebyshev-fit. F_k(τ) remains a closed-form callable.
        self.H_coeffs = np.zeros((len(self.terms), n_delta + 1))
        self.n_values = np.zeros(len(self.terms))
        for k, term in enumerate(self.terms):
            self.H_coeffs[k, :] = fit_cheb_on(term["H"], n_delta,
                                                self.delta_min, self.delta_max)
            self.n_values[k] = term["n"]

    def eval_F_all(self, tau):
        """Evaluate F_k analytically at a single τ. Vector of length n_terms."""
        return np.array([term["F"](tau) for term in self.terms])

    def flatten_alphar01(self, tau):
        """Return Chebyshev coefficients of Ã_01(δ) at the given τ.
        Ã_01(δ) = Σ_k n_k F_k(τ) · H̃_k(δ)."""
        F_at_tau = self.eval_F_all(tau)
        weights = self.n_values * F_at_tau
        return weights @ self.H_coeffs

    def eval_alphar01(self, tau, delta):
        """Evaluate α^r_01 via the Chebyshev representation (for validation)."""
        coeffs = self.flatten_alphar01(tau)
        return eval_cheb_on(coeffs, delta, self.delta_min, self.delta_max)

    def pressure_residual_coeffs(self, T, P_target):
        """Return Chebyshev coefficients of r(δ) = P(T, δ) - P_target on the
        domain, where P = ρ_r · δ · R · T · (1 + α^r_01(τ, δ)). The result is
        a single polynomial in δ whose real roots in the domain are the
        density solutions to P = P_target.

        Math:
            r(δ) = ρ_r R T [δ + δ · α̃_01(δ)] - P_target
        The "δ · α̃_01" is computed via the Chebyshev-identity multiply-by-x
        on the fitted interval.
        """
        tau = self.eos.Tc / T
        A01_coeffs = self.flatten_alphar01(tau)
        # delta * A01 in Chebyshev coefficients on [delta_min, delta_max]
        dA01 = times_x_on_interval(A01_coeffs, self.delta_min, self.delta_max)
        # delta as a Chebyshev series on [a, b] is just two coefficients: [mid, half]
        half = 0.5 * (self.delta_max - self.delta_min)
        mid = 0.5 * (self.delta_min + self.delta_max)
        delta_coeffs = np.zeros(2)
        delta_coeffs[0] = mid
        delta_coeffs[1] = half
        factor = self.eos.rhoc * self.eos.R * T
        prefactor_delta = factor * delta_coeffs
        prefactor_dA01 = factor * dA01
        resid = add_coeffs(add_coeffs(prefactor_delta, prefactor_dA01),
                           np.array([-P_target]))  # subtract P_target (const)
        return resid

    def density_roots(self, T, P_target, eigvals_tol=1e-9):
        """Return list of real density roots rho [mol/m^3] at which P(T, rho) == P_target,
        found via eigenvalue colleague matrix on the pressure-residual Chebyshev series."""
        resid = self.pressure_residual_coeffs(T, P_target)
        roots_std = npcheb.chebroots(resid)
        # filter real roots in [-1, 1] (mapped back to [delta_min, delta_max])
        real = roots_std[np.abs(roots_std.imag) < eigvals_tol].real
        in_range = real[(real >= -1.0 - 1e-10) & (real <= 1.0 + 1e-10)]
        in_range = np.clip(in_range, -1.0, 1.0)
        half = 0.5 * (self.delta_max - self.delta_min)
        mid = 0.5 * (self.delta_min + self.delta_max)
        delta_roots = mid + half * in_range
        rho_roots = delta_roots * self.eos.rhoc
        return np.sort(rho_roots)


# ---------------------------------------------------------------
# Smoke test / validation against direct IAPWS-95 evaluation
# ---------------------------------------------------------------
def validation_report():
    eos = IAPWS95Terms()
    tau_range = (0.2, 4.0)
    delta_range = (1e-10, 5.0)
    n_delta = 120
    print(f"Building FluidCheb with n_delta={n_delta} (F_k evaluated analytically)")
    fc = FluidCheb(eos, tau_range, delta_range, n_delta)
    print(f"  built {len(fc.terms)} H̃_k(δ) expansions in {n_delta + 1} coeffs each")

    # Ground truth: sum of power + gaussian terms only (excluding non-analytic)
    def truth_alphar01_separable(tau, delta):
        power = eos.alphar_01_power_terms(tau, delta).sum()
        gauss = eos.alphar_01_gauss_terms(tau, delta).sum()
        return power + gauss

    print()
    print(f"{'T':>6} {'rho':>12} {'τ':>6} {'δ':>7} "
          f"{'separable truth':>18} {'FluidCheb':>18} "
          f"{'rel err':>10} {'nonan':>12}")
    Tc, rhoc = eos.Tc, eos.rhoc
    for T in [300.0, 500.0, 700.0, 1000.0]:
        for rho in [1.0, 1000.0, rhoc, 30000.0, 50000.0]:
            tau = Tc / T; delta = rho / rhoc
            truth = truth_alphar01_separable(tau, delta)
            cheb_val = fc.eval_alphar01(tau, delta)
            nonan = eos.alphar_01_nonan_terms(tau, delta).sum()
            rel = abs(cheb_val - truth) / max(abs(truth), 1e-30)
            print(f"{T:>6.1f} {rho:>12.2f} {tau:>6.3f} {delta:>7.4f} "
                  f"{truth:>+18.6e} {cheb_val:>+18.6e} "
                  f"{rel:>10.2e} {nonan:>+12.3e}")


if __name__ == "__main__":
    validation_report()
