"""Generic fluid evaluator from a CoolProp/teqp EOS JSON. Supports
ResidualHelmholtzPower and ResidualHelmholtzGaussian blocks (the "simple"
separable case). Fluids with non-analytic or non-separable blocks will
error on load, so this is a clean test bed for the separable architecture.
"""

import json
import numpy as np


class FluidFromJSON:
    def __init__(self, json_path):
        with open(json_path) as f:
            j = json.load(f)
        eos = j["EOS"][0]
        self.R = eos["gas_constant"]
        self.Tc = eos["STATES"]["reducing"]["T"]
        self.rhoc = eos["STATES"]["reducing"]["rhomolar"]
        self.M = eos["molar_mass"]
        self.name = j.get("INFO", {}).get("NAME", "UNKNOWN")
        alphar = eos["alphar"]
        allowed = {"ResidualHelmholtzPower", "ResidualHelmholtzGaussian"}
        types_present = {a["type"] for a in alphar}
        unsupported = types_present - allowed
        if unsupported:
            raise NotImplementedError(
                f"Fluid {self.name} has unsupported alphar blocks: {unsupported}"
            )

        pwr = next((a for a in alphar if a["type"] == "ResidualHelmholtzPower"), None)
        if pwr is not None:
            self.n_p = np.array(pwr["n"], dtype=float)
            self.d_p = np.array(pwr["d"], dtype=float)
            self.t_p = np.array(pwr["t"], dtype=float)
            self.l_p = np.array(pwr["l"], dtype=float)
        else:
            self.n_p = self.d_p = self.t_p = self.l_p = np.array([])

        g = next((a for a in alphar if a["type"] == "ResidualHelmholtzGaussian"), None)
        if g is not None:
            self.n_g = np.array(g["n"], dtype=float)
            self.d_g = np.array(g["d"], dtype=float)
            self.t_g = np.array(g["t"], dtype=float)
            self.eta_g = np.array(g["eta"], dtype=float)
            self.eps_g = np.array(g["epsilon"], dtype=float)
            self.beta_g = np.array(g["beta"], dtype=float)
            self.gamma_g = np.array(g["gamma"], dtype=float)
        else:
            self.n_g = self.d_g = self.t_g = self.eta_g = self.eps_g = \
                self.beta_g = self.gamma_g = np.array([])

    def alphar_01_direct(self, tau, delta):
        """Direct summation of all terms' δ·∂α^r/∂δ. Ground truth for validation."""
        tot = 0.0
        if len(self.n_p) > 0:
            tp = tau ** self.t_p
            dp = delta ** self.d_p
            ep = np.where(self.l_p > 0,
                          np.exp(-delta ** np.where(self.l_p > 0, self.l_p, 0)),
                          1.0)
            factor = self.d_p - self.l_p * (delta ** self.l_p)
            tot += np.sum(self.n_p * dp * tp * ep * factor)
        if len(self.n_g) > 0:
            arg = -self.eta_g * (delta - self.eps_g) ** 2 \
                - self.beta_g * (tau - self.gamma_g) ** 2
            inner = self.d_g - 2 * self.eta_g * delta * (delta - self.eps_g)
            tot += np.sum(self.n_g * delta ** self.d_g * tau ** self.t_g
                          * np.exp(arg) * inner)
        return tot

    def pressure_direct(self, T, rho):
        tau = self.Tc / T
        delta = rho / self.rhoc
        return rho * self.R * T * (1.0 + self.alphar_01_direct(tau, delta))

    def separable_terms(self):
        """Return list of dicts like paper_architecture.iapws95_separable_terms."""
        terms = []
        for i in range(len(self.n_p)):
            n, d, t, l = float(self.n_p[i]), float(self.d_p[i]), \
                float(self.t_p[i]), float(self.l_p[i])

            def Fk(tau, _t=t):
                return tau ** _t

            def Hk(delta, _d=d, _l=l):
                if _l == 0:
                    return (delta ** _d) * _d
                return (delta ** _d) * np.exp(-delta ** _l) * (_d - _l * delta ** _l)
            terms.append(dict(F=Fk, H=Hk, n=n, family="power", index=i))
        for i in range(len(self.n_g)):
            n, d, t = float(self.n_g[i]), float(self.d_g[i]), float(self.t_g[i])
            eta, eps = float(self.eta_g[i]), float(self.eps_g[i])
            beta, gamma = float(self.beta_g[i]), float(self.gamma_g[i])

            def Fk(tau, _t=t, _b=beta, _g=gamma):
                return (tau ** _t) * np.exp(-_b * (tau - _g) ** 2)

            def Hk(delta, _d=d, _e=eta, _p=eps):
                return (delta ** _d) * np.exp(-_e * (delta - _p) ** 2) \
                    * (_d - 2 * _e * delta * (delta - _p))
            terms.append(dict(F=Fk, H=Hk, n=n, family="gauss", index=i))
        return terms


if __name__ == "__main__":
    import sys
    path = sys.argv[1] if len(sys.argv) > 1 else \
        "/Users/ianbell/miniforge3/lib/python3.9/site-packages/teqp/fluiddata/dev/fluids/Nitrogen.json"
    f = FluidFromJSON(path)
    print(f"Fluid {f.name}: Tc={f.Tc} K, rhoc={f.rhoc:.4f} mol/m^3, M={f.M} kg/mol")
    print(f"  power: {len(f.n_p)} terms (l set: {sorted(set(f.l_p.tolist()))})")
    print(f"  gaussian: {len(f.n_g)} terms")
    # smoke test vs teqp
    import teqp
    model = teqp.build_multifluid_model([f.name], teqp.get_datapath())
    z = np.array([1.0])
    T = 200.0
    for rho in [10.0, 1000.0, f.rhoc, 20000.0]:
        a_direct = f.alphar_01_direct(f.Tc / T, rho / f.rhoc)
        a_teqp = model.get_Ar01(T, rho, z)
        print(f"  T={T}, rho={rho}: direct={a_direct:+.10e}, teqp={a_teqp:+.10e}, "
              f"rel err {abs(a_direct-a_teqp)/max(abs(a_teqp), 1e-30):.2e}")
