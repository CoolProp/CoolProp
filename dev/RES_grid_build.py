"""Build the (fluid, T, p) grid used by the Stage-6 RES verification.

The published sample sets carry ONE point per fluid, and several of those points sit almost
entirely in the dilute gas -- where the transport property barely depends on s_res at all, so a
large s_res error cannot show up.  That is what made the first round of HEOS-vs-REFPROP
exclusions unreliable.  This grid is built in REDUCED coordinates instead, so that every fluid is
sampled at states where the residual contribution actually dominates.

Phase ambiguity is designed out rather than handled: a (T, p) pair near the saturation line can
land HEOS and REFPROP on opposite sides of the dome, which inflates the deviation without saying
anything about the parameters.  So the grid uses supercritical isotherms, compressed liquid well
above p_sat, and superheated vapour well below it -- never a point close to saturation.

Writes dev/RES_comparison/grid_points.csv, consumed by the [RES_grid] Catch2 test.
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
        w.writerow(["fluid", "T_K", "p_Pa", "region"])
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
                w.writerow([fluid, repr(T), repr(p), region])
                n_pts += 1

    print(f"Wrote {out_path}")
    print(f"  {n_pts} points over {n_fluid} fluids")
    if skipped:
        print(f"  {len(skipped)} fluids skipped (not constructible in REFPROP here):")
        for s in skipped:
            print(f"    {s}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fluids", default=None, help="comma-separated subset (default: every fluid in the JSON)")
    args = ap.parse_args()
    sys.exit(main(args.fluids))
