"""Refit the rhosr-CS viscosity constants for R1233zd(E).

Reproduces the two constants in ``dev/fluids/R1233zd(E).json`` ->
``TRANSPORT.viscosity``: ``rhosr_critical`` and ``C``.  See CoolProp/CoolProp#3330
and CoolProp/CoolProp#1826.

``rhosr_critical`` is not fitted.  It is rho_c * R * (tau * dalphar_dtau - alphar)
evaluated at the EOS critical point -- the quantity that normalises the reduced
variable x = rho * s_r / rhosr_critical -- so it is a property of whichever
equation of state ships, and is simply recomputed here.

``C`` is the one fluid-specific fitted parameter.  Because

    eta = eta_dilute * (1 + C * (etastar_ref - 1))

is linear in C at fixed state, the objective is evaluated on a two-point basis
(C = 0 and C = 1) rather than by re-entering the backend for every trial value.

Data: Miyara, A., Alam, M. J. and Kariya, K., "Measurement of viscosity of
trans-1-chloro-3,3,3-trifluoropropene (R-1233zd(E)) by tandem capillary tubes
method", Int. J. Refrigeration 92 (2018) 86-93, doi:10.1016/j.ijrefrig.2018.05.021.
Table 3, liquid phase: 61 points, 313.67-433.32 K, 1.005-4.074 MPa, stated
combined standard uncertainty 3.0 percent.  The accepted manuscript carries the
notice "This manuscript is made available under the Elsevier user license".

Usage:  python dev/scripts/fit_R1233zdE_viscosity.py
"""

import copy
import json
import os
import statistics

import CoolProp.CoolProp as CP

# T [K], p [MPa], eta [micro Pa s] -- Miyara et al. (2018), Table 3, liquid phase
MIYARA_LIQUID = [
    ( 313.90,  1.005,  244.50),
    ( 353.70,  3.056,  168.80),
    ( 313.90,  1.005,  243.81),
    ( 353.68,  3.055,  169.19),
    ( 313.81,  1.010,  243.98),
    ( 353.63,  2.046,  166.86),
    ( 313.83,  1.938,  246.29),
    ( 353.66,  2.039,  166.83),
    ( 313.89,  1.907,  245.77),
    ( 353.63,  2.042,  167.15),
    ( 313.88,  2.032,  246.54),
    ( 374.31,  4.049,  142.14),
    ( 313.85,  2.032,  246.47),
    ( 374.31,  4.048,  142.11),
    ( 313.83,  2.024,  245.77),
    ( 374.31,  4.048,  142.09),
    ( 313.91,  3.033,  247.58),
    ( 374.24,  3.031,  138.63),
    ( 313.91,  3.027,  248.21),
    ( 374.22,  3.030,  138.63),
    ( 313.91,  3.027,  247.73),
    ( 374.23,  3.032,  138.47),
    ( 313.87,  4.025,  250.11),
    ( 374.24,  2.031,  135.40),
    ( 313.67,  4.028,  249.04),
    ( 374.24,  2.031,  135.45),
    ( 313.77,  4.023,  248.75),
    ( 374.23,  2.033,  135.41),
    ( 313.78,  4.022,  249.55),
    ( 393.52,  4.031,  116.79),
    ( 333.90,  4.037,  204.78),
    ( 393.52,  4.031,  116.49),
    ( 334.00,  4.041,  205.16),
    ( 393.52,  4.031,  116.57),
    ( 334.00,  4.044,  204.80),
    ( 393.50,  3.087,  112.17),
    ( 333.83,  3.050,  203.71),
    ( 393.54,  3.011,  112.10),
    ( 333.93,  3.011,  202.96),
    ( 393.57,  3.007,  112.17),
    ( 333.90,  3.015,  203.22),
    ( 393.57,  3.007,  111.76),
    ( 333.93,  2.050,  201.74),
    ( 413.49,  4.074,   90.54),
    ( 333.93,  2.048,  201.61),
    ( 413.49,  4.056,   90.47),
    ( 333.94,  2.050,  201.55),
    ( 413.55,  4.031,   91.80),
    ( 333.91,  1.033,  199.83),
    ( 412.90,  3.027,   87.03),
    ( 333.93,  1.037,  199.07),
    ( 413.31,  3.052,   86.48),
    ( 333.90,  1.037,  199.15),
    ( 413.45,  3.059,   86.11),
    ( 353.68,  4.043,  172.90),
    ( 433.32,  3.927,   62.85),
    ( 353.68,  4.037,  172.53),
    ( 433.32,  3.927,   62.49),
    ( 353.65,  4.027,  172.44),
    ( 433.32,  3.927,   62.67),
    ( 353.67,  3.054,  169.34),
]

FLUID_JSON = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "fluids", "R1233zd(E).json")


def probe(C, rhosr_critical):
    """Register a private clone of the fluid carrying the given constants."""
    with open(FLUID_JSON, encoding="utf-8") as fp:
        fluid = json.load(fp)
    fluid = copy.deepcopy(fluid[0] if isinstance(fluid, list) else fluid)
    probe.counter += 1
    name = "R1233zdE_fit_{0}".format(probe.counter)
    fluid["INFO"]["NAME"] = name
    fluid["INFO"]["ALIASES"] = []
    fluid["INFO"]["CAS"] = "999-00-{0}".format(probe.counter)
    fluid["INFO"]["REFPROP_NAME"] = "N/A"
    for key in ("INCHI_KEY", "INCHI_STRING", "SMILES", "CHEMSPIDER_ID", "2DPNG_URL"):
        fluid["INFO"].pop(key, None)
    fluid["TRANSPORT"]["viscosity"]["C"] = C
    fluid["TRANSPORT"]["viscosity"]["rhosr_critical"] = rhosr_critical
    CP.add_fluids_as_JSON("HEOS", json.dumps([fluid]))
    return CP.AbstractState("HEOS", name)


probe.counter = 0


def rhosr_critical_from_eos():
    """rho_c * R * (tau * dalphar_dtau - alphar) at the EOS critical point.

    Deliberately evaluated on a clone built from the repository's fluid file
    rather than on the built-in "R1233zd(E)", so that the constant printed here
    always belongs to the EOS in the tree.  Querying the installed CoolProp
    would silently report the wheel's EOS while C was being fitted against the
    repository's, the moment the two diverge.
    """
    state = probe(0.0, 0.0)  # transport constants are irrelevant to alphar
    Tc, rhoc = state.T_critical(), state.rhomolar_critical()
    state.update(CP.DmolarT_INPUTS, rhoc, Tc)
    return rhoc * state.gas_constant() * (state.tau() * state.dalphar_dTau() - state.alphar())


def main():
    """Recompute rhosr_critical, refit C, and print both with their fit statistics."""
    rhosr_critical = rhosr_critical_from_eos()
    print("rhosr_critical from the shipped EOS = {0:.17g}".format(rhosr_critical))

    # eta is linear in C, so two evaluations per state span the whole family.
    zero, one = probe(0.0, rhosr_critical), probe(1.0, rhosr_critical)
    basis = []
    for T, p_MPa, eta_uPas in MIYARA_LIQUID:
        zero.update(CP.PT_INPUTS, p_MPa * 1e6, T)
        dilute = zero.viscosity()
        one.update(CP.PT_INPUTS, p_MPa * 1e6, T)
        basis.append((dilute, one.viscosity() - dilute, eta_uPas * 1e-6))

    def aad(C):
        """Mean absolute relative deviation from the measurements, in percent."""
        return statistics.mean(abs((a + C * b) / e - 1) for a, b, e in basis) * 100

    # AAD(C) is a sum of terms |b_i/e_i| * |C - (e_i - a_i)/b_i|, i.e. piecewise
    # linear and convex with kinks only at those breakpoints.  The exact
    # minimiser is therefore one of them -- no grid, and no chance of silently
    # returning a scan boundary as if it were an optimum.
    candidates = [(e - a) / b for a, b, e in basis if b != 0.0]
    C_exact = min(candidates, key=aad)
    C = round(C_exact, 4)

    # The shipped constant is the rounded one, so report the rounded one; assert
    # the rounding costs nothing meaningful against a 3.0% measurement uncertainty.
    assert abs(aad(C) - aad(C_exact)) < 1e-3, "rounding C to 4 dp moved the AAD"

    deviations = [((a + C * b) / e - 1) * 100 for a, b, e in basis]
    print("C = {0:.4f}   (exact minimiser {1:.9f})".format(C, C_exact))
    print(
        "  N = {0}   AAD = {1:.2f}%   bias = {2:+.2f}%   range {3:+.1f}% .. {4:+.1f}%".format(
            len(basis), aad(C), statistics.mean(deviations), min(deviations), max(deviations)
        )
    )


if __name__ == "__main__":
    main()
