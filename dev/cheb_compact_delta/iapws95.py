"""Minimal per-term IAPWS-95 residual Helmholtz evaluator, driven off the
CoolProp/teqp JSON. Exposes each of the four IAPWS-95 term families
separately so per-term Chebyshev-in-u convergence can be studied:

    family 'power'   : 51 terms, n δ^d τ^t exp(-δ^l) with l in {0,1,2,3,4,6}
    family 'gauss'   : 3 terms,  n δ^d τ^t exp(-η(δ-ε)^2 - β(τ-γ)^2)
    family 'nonan'   : 2 terms,  n Δ^b δ ψ, see IAPWS-95 release eq. 6.6

The quantity used for pressure inversion is
    alphar_01 := δ · ∂α^r/∂δ  (at fixed τ)
since  P = ρ R T (1 + alphar_01).

All functions return either α^r or alphar_01 per term, or summed by family.
"""

import json
import numpy as np

DEFAULT_WATER_JSON = "/Users/ianbell/miniforge3/lib/python3.9/site-packages/teqp/fluiddata/dev/fluids/Water.json"


class IAPWS95Terms:
    Tc = 647.096       # K
    rhoc = 322.0 / 18.015268e-3   # mol/m^3
    R = 8.314462618    # IAPWS-95 as shipped in teqp uses this; check below

    def __init__(self, json_path=DEFAULT_WATER_JSON):
        with open(json_path) as f:
            j = json.load(f)
        eos = j["EOS"][0]
        self.R = eos["gas_constant"]
        self.Tc = eos["STATES"]["reducing"]["T"]
        self.rhoc = eos["STATES"]["reducing"]["rhomolar"]
        alphar = eos["alphar"]
        # power block
        pwr = next(a for a in alphar if a["type"] == "ResidualHelmholtzPower")
        self.n_p = np.array(pwr["n"], dtype=float)
        self.d_p = np.array(pwr["d"], dtype=float)
        self.t_p = np.array(pwr["t"], dtype=float)
        self.l_p = np.array(pwr["l"], dtype=float)
        # gaussian
        g = next(a for a in alphar if a["type"] == "ResidualHelmholtzGaussian")
        self.n_g = np.array(g["n"], dtype=float)
        self.d_g = np.array(g["d"], dtype=float)
        self.t_g = np.array(g["t"], dtype=float)
        self.eta_g = np.array(g["eta"], dtype=float)
        self.eps_g = np.array(g["epsilon"], dtype=float)
        self.beta_g = np.array(g["beta"], dtype=float)
        self.gamma_g = np.array(g["gamma"], dtype=float)
        # non-analytic
        na = next(a for a in alphar if a["type"] == "ResidualHelmholtzNonAnalytic")
        self.n_na = np.array(na["n"], dtype=float)
        self.a_na = np.array(na["a"], dtype=float)
        self.b_na = np.array(na["b"], dtype=float)
        self.beta_na = np.array(na["beta"], dtype=float)
        self.A_na = np.array(na["A"], dtype=float)
        self.B_na = np.array(na["B"], dtype=float)
        self.C_na = np.array(na["C"], dtype=float)
        self.D_na = np.array(na["D"], dtype=float)

    # --------- Power: alpha^r term i = n_i δ^d_i τ^t_i exp(-δ^l_i) ----------
    def alphar_power_terms(self, tau, delta):
        """Return array of 51 per-term alpha^r values."""
        tp = tau ** self.t_p
        dp = delta ** self.d_p
        ep = np.exp(-delta ** self.l_p) if np.any(self.l_p > 0) else 1.0
        # handle l=0 (no exp)
        mask0 = self.l_p == 0
        res = self.n_p * dp * tp
        # multiply exp factor only where l > 0
        res_exp = res * np.exp(-delta ** np.where(self.l_p > 0, self.l_p, 0))
        # combine
        return np.where(mask0, res, res_exp)

    def alphar_01_power_terms(self, tau, delta):
        """Return array of 51 per-term δ · ∂α^r_i/∂δ values (pressure contribution)."""
        tp = tau ** self.t_p
        # d/d(delta) [ delta^d exp(-delta^l) ] = delta^(d-1) exp(-delta^l) * (d - l delta^l)
        # multiplied by delta : delta^d exp(-delta^l) * (d - l delta^l)
        dp = delta ** self.d_p
        ep = np.where(self.l_p > 0,
                      np.exp(-delta ** np.where(self.l_p > 0, self.l_p, 0)),
                      1.0)
        factor = self.d_p - self.l_p * (delta ** self.l_p)  # = d - l·δ^l
        return self.n_p * dp * tp * ep * factor

    # --------- Gaussian bell ------------------------------------------------
    def alphar_gauss_terms(self, tau, delta):
        arg = -self.eta_g * (delta - self.eps_g) ** 2 - self.beta_g * (tau - self.gamma_g) ** 2
        return self.n_g * delta ** self.d_g * tau ** self.t_g * np.exp(arg)

    def alphar_01_gauss_terms(self, tau, delta):
        # d/d(delta) [ delta^d * exp(-eta(delta-eps)^2 - beta(tau-gamma)^2) ]
        # = delta^(d-1) * (d - 2 eta delta (delta-eps)) * exp(...)
        # multiply by delta:
        arg = -self.eta_g * (delta - self.eps_g) ** 2 - self.beta_g * (tau - self.gamma_g) ** 2
        inner = self.d_g - 2 * self.eta_g * delta * (delta - self.eps_g)
        return self.n_g * delta ** self.d_g * tau ** self.t_g * np.exp(arg) * inner

    # --------- Non-analytic -------------------------------------------------
    def _Delta_theta_psi(self, tau, delta):
        # See IAPWS-95 release, eqs. for Δ, θ, ψ
        # theta = (1 - tau) + A_i * ((delta-1)^2)^(1/(2 beta_i))
        # Delta = theta^2 + B_i * ((delta-1)^2)^a_i
        # psi   = exp(-C_i (delta-1)^2 - D_i (tau-1)^2)
        dm1sq = (delta - 1.0) ** 2
        # ((delta-1)^2)^x needs care at delta=1: (0)^x = 0 for x>0
        # Use np.power for safety
        theta = (1.0 - tau) + self.A_na * np.power(dm1sq + 0.0, 1.0 / (2.0 * self.beta_na))
        Delta = theta ** 2 + self.B_na * np.power(dm1sq + 0.0, self.a_na)
        psi = np.exp(-self.C_na * dm1sq - self.D_na * (tau - 1.0) ** 2)
        return Delta, theta, psi

    def alphar_nonan_terms(self, tau, delta):
        Delta, theta, psi = self._Delta_theta_psi(tau, delta)
        # alpha^r_i = n_i Delta^b_i delta psi
        return self.n_na * np.power(Delta, self.b_na) * delta * psi

    def alphar_01_nonan_terms(self, tau, delta):
        """
        delta * d(alpha^r_i)/d(delta) for the non-analytic terms.
        Use IAPWS-95 release derivative formulas (eqs. 5.5, 5.6 of the 2018 update).
          alpha^r_i = n_i Delta^b_i delta psi
          d/d(delta) [Delta^b delta psi]
             = Delta^b (psi + delta dpsi/ddelta) + b Delta^{b-1} dDelta/ddelta * delta psi
        where
          dDelta/ddelta = 2 theta * dtheta/ddelta + 2 a_i B_i ((delta-1)^2)^(a_i - 1) (delta-1)
          dtheta/ddelta = A_i / beta_i * ((delta-1)^2)^(1/(2 beta_i) - 1) * (delta-1)
          dpsi/ddelta = -2 C_i (delta-1) psi
        """
        dm1 = delta - 1.0
        dm1sq = dm1 * dm1
        # Avoid 0^negative issues when delta exactly equals 1.
        # At delta==1 the contributions with (dm1sq)^(alpha-1) go to 0 for alpha>1,
        # and cusp-like for alpha<1. We regularize minimally.
        tiny = 1e-300
        dm1sq_safe = dm1sq if dm1sq > 0 else tiny
        Delta, theta, psi = self._Delta_theta_psi(tau, delta)

        dtheta_dd = (self.A_na / self.beta_na) * \
            np.power(dm1sq_safe, 1.0 / (2.0 * self.beta_na) - 1.0) * dm1
        dDelta_dd = 2 * theta * dtheta_dd + \
            2 * self.a_na * self.B_na * np.power(dm1sq_safe, self.a_na - 1.0) * dm1
        dpsi_dd = -2 * self.C_na * dm1 * psi

        # d/d(delta) [ Delta^b * delta * psi ]
        term1 = np.power(Delta, self.b_na) * (psi + delta * dpsi_dd)
        term2 = self.b_na * np.power(Delta, self.b_na - 1.0) * dDelta_dd * delta * psi
        dalphar_dd = self.n_na * (term1 + term2)
        return delta * dalphar_dd  # return delta * d(α^r)/d(δ)

    # --------- Summed forms -------------------------------------------------
    def alphar(self, tau, delta):
        return (self.alphar_power_terms(tau, delta).sum()
                + self.alphar_gauss_terms(tau, delta).sum()
                + self.alphar_nonan_terms(tau, delta).sum())

    def alphar_01(self, tau, delta):
        """delta * d(alpha^r)/d(delta)."""
        return (self.alphar_01_power_terms(tau, delta).sum()
                + self.alphar_01_gauss_terms(tau, delta).sum()
                + self.alphar_01_nonan_terms(tau, delta).sum())

    def pressure(self, T, rho):
        tau = self.Tc / T
        delta = rho / self.rhoc
        return rho * self.R * T * (1.0 + self.alphar_01(tau, delta))


def smoke_test_against_teqp():
    import teqp
    eos = IAPWS95Terms()
    model = teqp.build_multifluid_model(["Water"], teqp.get_datapath())
    z = np.array([1.0])
    Tc, rhoc = eos.Tc, eos.rhoc
    print(f"Local: Tc={Tc}, rhoc={rhoc:.4f}")
    print(f"teqp : Tc={model.get_Tr(z)}, rhoc={model.get_rhor(z):.4f}")
    print()
    print(f"{'T':>8} {'rho':>12} {'Ar01_teqp':>18} {'Ar01_local':>18} "
          f"{'rel_err':>10}")
    for T in [300.0, 500.0, 700.0, 1000.0]:
        for rho in [1.0, 1000.0, 17873.73, 30000.0, 50000.0]:
            a_t = model.get_Ar01(T, rho, z)
            a_l = eos.alphar_01(Tc / T, rho / rhoc)
            rel = abs(a_t - a_l) / max(abs(a_t), 1e-30)
            print(f"{T:>8.1f} {rho:>12.4f} {a_t:>+18.10e} {a_l:>+18.10e} "
                  f"{rel:>10.2e}")


if __name__ == "__main__":
    smoke_test_against_teqp()
