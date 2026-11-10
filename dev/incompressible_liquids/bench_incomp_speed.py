"""Time INCOMP property calls on the classic polynomial fits.

Run with CoolProp installed:
    python dev/incompressible_liquids/bench_incomp_speed.py

Set plot = True in the main block to plot the results (needs matplotlib).
The timings include the per-call Python wrapper overhead. The bench_ prefix
(not test_) keeps pytest in the incomp-json-sanity job, which does not
install CoolProp, from collecting this file.
"""

import time

import CoolProp.CoolProp as CP
import numpy as np


def timed(label, func, values):
    start = time.perf_counter()
    out = [func(v) for v in values]
    print(f"{label} elapsed time: {time.perf_counter() - start} seconds")
    return out


if __name__ == "__main__":
    n = int(1e6)
    plot = False

    # CoolProp 8.0.0 has no such key and its caloric fits are classic already.
    if hasattr(CP, "INCOMP_PREFER_CHEBYSHEV_CALORIC"):
        CP.set_config_bool(CP.INCOMP_PREFER_CHEBYSHEV_CALORIC, False)

    for fld_string in ["MEG[0.35]", "Water"]:
        fld = fld_string.split("[")[0]
        fldconc = (
            [float(fld_string.split("[")[1].split("]")[0])]
            if "[" in fld_string
            else None
        )
        fluid = CP.AbstractState("INCOMP", fld)
        if fldconc:
            fluid.set_mass_fractions(fldconc)

        print(f"\nStarting CoolProp test: {fld_string}")
        print("-" * 40)

        T = np.linspace(300, 350, n)
        fluid.update(CP.PT_INPUTS, 1e5, T[0])
        h1, s1 = fluid.hmass(), fluid.smass()
        fluid.update(CP.PT_INPUTS, 1e5, T[-1])
        h2, s2 = fluid.hmass(), fluid.smass()
        h = np.linspace(h1, h2, n)
        s = np.linspace(s1, s2, n)

        def at_T(getter, p=1e5):
            def call(_T):
                fluid.update(CP.PT_INPUTS, p, _T)
                return getter()

            return call

        if fldconc:
            timed(
                "Tfreeze(T,p)", at_T(lambda: fluid.keyed_output(CP.iT_freeze), p=0.0), T
            )

        if fld == "Water":

            def psat(_T):
                fluid.update(CP.QT_INPUTS, 0.0, _T)
                return fluid.p()

            timed("psat(Q,T)", psat, T)

        mu = timed("mu(T,p)", at_T(fluid.viscosity), T)
        k = timed("k(T,p)", at_T(fluid.conductivity), T)
        d = timed("d(T,p)", at_T(fluid.rhomass), T)
        cp = timed("cp(T,p)", at_T(fluid.cpmass), T)

        def T_of_h(_h):
            fluid.update(CP.HmassP_INPUTS, _h, 1e5)
            return fluid.T()

        def T_of_s(_s):
            fluid.update(CP.PSmass_INPUTS, 1e5, _s)
            return fluid.T()

        Touth = timed("T(p,h)", T_of_h, h)
        Touts = timed("T(p,s)", T_of_s, s)

        if plot:
            import matplotlib.pyplot as plt

            fig, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
            panels = [
                (mu, "mu(T,p)", "Viscosity [Pa s]"),
                (k, "k(T,p)", "Conductivity [W/m/K]"),
                (d, "d(T,p)", "Density [kg/m3]"),
                (cp, "cp(T,p)", "Specific heat [J/kg/K]"),
                (Touth, "T(p,h)", "Temperature [K]"),
                (Touts, "T(p,s)", "Temperature [K]"),
            ]
            for ax, (values, title, ylabel) in zip(axes.flat, panels):
                ax.plot(T, values)
                ax.set_title(title)
                ax.set_xlabel("Temperature [K]")
                ax.set_ylabel(ylabel)
            fig.suptitle(f"Properties of {fld_string}")
            plt.show()
