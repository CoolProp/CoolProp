"""Build the (fluid, T, p) grid used to verify the RES transport implementation.

The published sample sets carry ONE point per fluid, and several of those points sit almost
entirely in the dilute gas -- where the transport property barely depends on s_res at all, so a
large s_res error cannot show up.  That is what made the first round of HEOS-vs-REFPROP
exclusions unreliable.  This grid is built in REDUCED coordinates instead, so that every fluid is
sampled at states where the residual contribution actually dominates.

Phase ambiguity is designed out rather than handled: a (T, p) pair near the saturation line can
land HEOS and REFPROP on opposite sides of the dome, which inflates the deviation without saying
anything about the parameters.  So the grid uses supercritical isotherms, compressed liquid well
above p_sat, and superheated vapour well below it -- never a point close to saturation.

Binary mixtures are included as well.  They matter for a different reason than the pure fluids:
they are the only thing that exercises the Wilke rule for the dilute term and the mole-fraction
averaging of the residual coefficients, and the published tables carry just 17 viscosity and 8
conductivity binary points between them -- too few, and all at states the authors chose.

Writes dev/RES_comparison/grid_points.csv, consumed by the [RES_grid] Catch2 test.  The
`mole_fractions` column is empty for pure fluids and semicolon-separated for mixtures.
"""

import argparse
import csv
import os
import sys

import CoolProp.CoolProp as CP

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
OUT_DIR = os.path.join(SCRIPT_DIR, "RES_comparison")

# Supercritical block: no dome to straddle, and rho/rhoc up to 2.0 reaches the dense states the
# RES residual term is actually fitted for.  rho/rhoc = 0.05 is kept as a dilute control -- it is
# the regime the published samples over-represent, so it belongs in the report as a contrast.
SUPERCRIT_TR = (1.02, 1.10, 1.30, 1.60)
# NOT 2.00: that is exactly where Li's enhancement gate switches off (rho/rhoc >= 2), and the two
# backends have slightly different rhoc, so a point placed there lands inside the window on one and
# outside it on the other.  The resulting step is an artefact of the grid, not a property of the
# parameters, and it silently inflated the conductivity exclusion list by ~10 fluids.  1.80 sits
# clearly inside the window and 2.30 clearly outside it, for both backends.
SUPERCRIT_RHOR = (0.05, 0.50, 1.00, 1.50, 1.80, 2.30)

# Subcritical block, deliberately far from saturation in both directions.
LIQUID_TR = (0.60, 0.75, 0.90)  # evaluated at 3x p_sat, i.e. compressed liquid
VAPOUR_TR = (0.80, 0.95)  # evaluated at 0.3x p_sat, i.e. superheated vapour

# Binary mixtures, taken from the pairs the two papers tabulate so the components are certain to
# carry parameters.  Compositions are deliberately away from 0.5/0.5 where the paper used
# something else, since an even split can hide an error in the mole-fraction weighting.
MIXTURES = [
    ("METHANE", "ETHANE", 0.6, 0.4),
    ("BUTANE", "METHANE", 0.606, 0.394),
    ("ARGON", "NITROGEN", 0.5, 0.5),
    ("ARGON", "NEON", 0.3, 0.7),
    ("CO2", "ETHANE", 0.7, 0.3),
    ("CO2", "NITROGEN", 0.8, 0.2),
    ("DECANE", "METHANE", 0.4, 0.6),
    ("BENZENE", "C16", 0.5, 0.5),
]
# Mixtures are sampled only supercritically: below the mixture critical point a (T, p) pair can
# sit inside the two-phase envelope, and the reference implementation has no more claim to the
# right root there than we do.
MIX_TR = (1.05, 1.30)
MIX_RHOR = (0.30, 1.00, 1.80)


def refprop_name(fluid: str) -> str:
    return "REFPROP::" + fluid


def build_points(fluid: str) -> list[tuple[float, float, str]]:
    """Return [(T_K, p_Pa, region)] for one fluid, using REFPROP as the reference EOS."""
    rp = refprop_name(fluid)
    Tc = CP.PropsSI("Tcrit", "T", 0, "p", 0, rp)
    rhoc = CP.PropsSI("rhocrit", "T", 0, "p", 0, rp)
    Ttrip = CP.PropsSI("Ttriple", "T", 0, "p", 0, rp)
    pts = []

    for tr in SUPERCRIT_TR:
        T = tr * Tc
        for rr in SUPERCRIT_RHOR:
            rho = rr * rhoc
            try:
                p = CP.PropsSI("P", "T", T, "Dmass", rho, rp)
            except Exception:
                continue
            if p > 0:
                pts.append((T, p, f"sc_Tr{tr}_rhor{rr}"))

    for tr in LIQUID_TR:
        T = tr * Tc
        if T <= Ttrip:
            continue
        try:
            psat = CP.PropsSI("P", "T", T, "Q", 0, rp)
            p = 3.0 * psat
            # A pressure barely above p_sat is still ambiguous in practice; push well past it.
            p = max(p, psat + 0.2 * CP.PropsSI("Pcrit", "T", 0, "p", 0, rp))
        except Exception:
            continue
        pts.append((T, p, f"liq_Tr{tr}"))

    for tr in VAPOUR_TR:
        T = tr * Tc
        if T <= Ttrip:
            continue
        try:
            psat = CP.PropsSI("P", "T", T, "Q", 1, rp)
        except Exception:
            continue
        pts.append((T, 0.3 * psat, f"vap_Tr{tr}"))

    return pts


def flash_is_consistent(names: list[str], z: list[float], T: float, p: float) -> bool:
    """True when the PT flash lands on a state that actually corresponds to the pressure asked for.

    For a few mixtures it does not.  ARGON[0.3]&NEON[0.7] at T = 96.98 K is the clearest case: we
    ask REFPROP for the pressure at rho = 355.63 kg/m3, get 7.6565e6 Pa, flash back at that (p, T)
    and land on rho = 326.24 -- 8.3% away, and the EOS pressure there is 7.1646e6 Pa, not the
    7.6565e6 requested.  Nothing raises, because the state caches the pressure it was GIVEN rather
    than recomputing it, so p() agrees with the input by construction.

    Comparing the cached pressure against the EOS pressure at the returned density is what exposes
    it.  The identity p = rho*R*T*(1 + delta*alphar_delta) holds to ~1e-16 on a converged state.

    Only mixtures are affected: they are placed by reducing with CRITPdll's critical-point
    ESTIMATE, which for a pair like Ar/Ne is not a physical density (823 kg/m3), so the requested
    states can sit where the flash struggles.  All 3828 pure points pass this check.
    """
    try:
        AS = CP.AbstractState("REFPROP", "&".join(names))
        if len(names) > 1:
            AS.set_mole_fractions(z)
        AS.update(CP.PT_INPUTS, p, T)
        d = AS.delta()
        p_alphar = AS.rhomolar() * AS.gas_constant() * AS.T() * (1 + d * AS.dalphar_dDelta())
        return abs(p_alphar / AS.p() - 1) < 1e-9
    except Exception:
        return False


def build_mixture_points(a: str, b: str, xa: float, xb: float) -> list[tuple[float, float, str]]:
    """Reduced-coordinate points for one binary, using REFPROP's mixture critical estimate."""
    mix = f"REFPROP::{a}[{xa}]&{b}[{xb}]"
    Tc = CP.PropsSI("Tcrit", "T", 0, "p", 0, mix)
    rhoc = CP.PropsSI("rhocrit", "T", 0, "p", 0, mix)
    pts = []
    for tr in MIX_TR:
        T = tr * Tc
        for rr in MIX_RHOR:
            try:
                p = CP.PropsSI("P", "T", T, "Dmass", rr * rhoc, mix)
            except Exception:
                continue
            if p > 0 and flash_is_consistent([a, b], [xa, xb], T, p):
                pts.append((T, p, f"mix_Tr{tr}_rhor{rr}"))
    return pts


def main(fluids_arg: str | None) -> None:
    import json

    with open(os.path.join(SCRIPT_DIR, "res_transport_parameters.json"), encoding="utf-8") as fh:
        res = json.load(fh)

    if fluids_arg:
        fluids = [f.strip().upper() for f in fluids_arg.split(",") if f.strip()]
    else:
        # Only fluids that BOTH backends can build are useful for a transfer measurement, and only
        # fluids that still carry a HEOS key can be evaluated on HEOS at all.  Fluids currently in
        # the exclusion lists are included on purpose -- re-deriving those exclusions is the point
        # of this grid -- so the caller is expected to have run the converter with --keep-all-heos.
        fluids = sorted(set(res["viscosity"]) & set(res["conductivity"]))

    os.makedirs(OUT_DIR, exist_ok=True)
    out_path = os.path.join(OUT_DIR, "grid_points.csv")
    n_fluid = 0
    n_pts = 0
    skipped = []
    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["fluid", "T_K", "p_Pa", "region", "mole_fractions"])
        for fluid in fluids:
            try:
                pts = build_points(fluid)
            except Exception as exc:
                skipped.append(f"{fluid} ({str(exc).splitlines()[0][:50]})")
                continue
            if not pts:
                skipped.append(f"{fluid} (no usable points)")
                continue
            n_fluid += 1
            for T, p, region in pts:
                w.writerow([fluid, repr(T), repr(p), region, ""])
                n_pts += 1

        n_mix = 0
        if not fluids_arg:
            for a, bb, xa, xb in MIXTURES:
                try:
                    pts = build_mixture_points(a, bb, xa, xb)
                except Exception as exc:
                    skipped.append(f"{a}&{bb} ({str(exc).splitlines()[0][:50]})")
                    continue
                for T, p, region in pts:
                    w.writerow([f"{a}&{bb}", repr(T), repr(p), region, f"{xa};{xb}"])
                    n_pts += 1
                    n_mix += 1

    print(f"Wrote {out_path}")
    print(f"  {n_pts} points over {n_fluid} fluids ({n_mix} of them binary-mixture points)")
    if skipped:
        print(f"  {len(skipped)} fluids skipped (not constructible in REFPROP here):")
        for s in skipped:
            print(f"    {s}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fluids", default=None, help="comma-separated subset (default: every fluid in the JSON)")
    args = ap.parse_args()
    sys.exit(main(args.fluids))
