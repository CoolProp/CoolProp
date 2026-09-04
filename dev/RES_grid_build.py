"""Build the (fluid, T, p) grid the RES validation harness evaluates.

The published sample sets carry ONE point per fluid, and several of those points sit almost
entirely in the dilute gas -- where the transport property barely depends on s_res at all, so a
large s_res error cannot show up.  That is what made the first round of HEOS exclusions
unreliable.  This grid is built in REDUCED coordinates instead, so that every fluid is sampled at
states where the residual contribution actually dominates.

Phase ambiguity is designed out rather than handled: a (T, p) pair near the saturation line can
land two equations of state on opposite sides of the dome, which inflates the deviation without
saying anything about the parameters.  So the grid uses supercritical isotherms, compressed
liquid well above p_sat, and superheated vapour well below it -- never a point close to
saturation.

Binary mixtures are included as well.  They matter for a different reason than the pure fluids:
they are the only thing that exercises the Wilke rule for the dilute term and the mole-fraction
averaging of the residual coefficients, and the published tables carry just 17 viscosity and 8
conductivity binary points between them -- too few, and all at states the authors chose.

Points are PLACED with HEOS, because HEOS is the backend under test; the reference values at
those points are produced separately by dev/RES_reference_run.py.

Placement is nonetheless LIMIT-AWARE, and that does need REFPROP -- to read its declared limits,
not to evaluate anything.  A point outside an equation of state's stated range is not a point:
REFPROP refuses it, and it refuses hardest at the dense end, which is exactly where the residual
term dominates and where the report is allowed to draw a conclusion.  So the ladder is clamped to
what both backends will answer, and where a rung is clipped the boundary state is substituted for
it, keeping the hottest and densest legal state of every fluid in the grid.

REFPROP answers up to 1.5x its declared Tmax and 2x its declared pmax and raises beyond.  Sampling
into that extrapolation band is deliberate: both sides run the SAME reference code at the SAME
(T, p), so what is compared there is still only the equation of state, which is the question the
grid exists to answer.  If REFPROP cannot be loaded, its limits simply do not constrain and the
grid falls back to HEOS's own.

Writes dev/RES_comparison/grid_points.csv.  The `mole_fractions` column is empty for pure fluids
and semicolon-separated for mixtures.
"""

import argparse
import collections
import csv
import json
import os

import CoolProp.CoolProp as CP

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
OUT_DIR = os.path.join(SCRIPT_DIR, "RES_comparison")

# Supercritical block: no dome to straddle, and rho/rhoc up to 2.3 reaches the dense states the
# RES residual term is actually fitted for.  rho/rhoc = 0.05 is kept as a dilute control -- it is
# the regime the published samples over-represent, so it belongs in the report as a contrast.
SUPERCRIT_TR = (1.02, 1.10, 1.30, 1.60)
# NOT 2.00: that is exactly where Li's enhancement gate switches off (rho/rhoc >= 2), and two
# equations of state have slightly different rhoc, so a point placed there lands inside the window
# on one and outside it on the other.  The resulting step is an artefact of the grid, not a
# property of the parameters, and it silently inflated an earlier conductivity exclusion list by
# ~10 fluids.  1.80 sits clearly inside the window and 2.30 clearly outside it, for both.
SUPERCRIT_RHOR = (0.05, 0.50, 1.00, 1.50, 1.80, 2.30)

# Subcritical block, deliberately far from saturation in both directions.
LIQUID_TR = (0.60, 0.75, 0.90)  # evaluated at 3x p_sat, i.e. compressed liquid
VAPOUR_TR = (0.80, 0.95)  # evaluated at 0.3x p_sat, i.e. superheated vapour

# Binary mixtures, taken from the pairs the two papers tabulate so the components are certain to
# carry parameters.  Compositions are deliberately away from 0.5/0.5 where the paper used
# something else, since an even split can hide an error in the mole-fraction weighting.
# Each pair must (a) be tabulated by one of the two papers, so both components certainly carry
# parameters, and (b) have a binary interaction pair in CoolProp.  ARGON&NEON and BENZENE&C16 are
# in Li's table but fail (b) -- no interaction pair, and C16 is not a CoolProp fluid -- so they
# are replaced by two other pairs from the same tables rather than left in to be skipped silently.
MIXTURES = [
    ("METHANE", "ETHANE", 0.6, 0.4),
    ("BUTANE", "METHANE", 0.606, 0.394),
    ("ARGON", "NITROGEN", 0.5, 0.5),
    ("ETHANOL", "WATER", 0.281, 0.719),
    ("CO2", "ETHANE", 0.7, 0.3),
    ("CO2", "NITROGEN", 0.8, 0.2),
    ("DECANE", "METHANE", 0.4, 0.6),
    ("ACETONE", "BENZENE", 0.092, 0.908),
]
# Mixtures are sampled only supercritically: below the mixture critical point a (T, p) pair can
# sit inside the two-phase envelope, and neither implementation has more claim to the right root
# there than the other.
MIX_TR = (1.05, 1.30)
MIX_RHOR = (0.30, 1.00, 1.80)


# REFPROP's documented extrapolation ceiling, confirmed by measurement (see the module docstring).
# The same factors are applied to HEOS: it does not enforce its limits on a PT flash -- it
# extrapolates silently, and occasionally returns NaN instead of raising -- so left alone it would
# never be the binding constraint even where its own range is the narrower of the two.
T_EXTRAP_FACTOR = 1.5
P_EXTRAP_FACTOR = 2.0
CAP_MARGIN = 0.98  # stay just inside the gate rather than exactly on it

# Reduced coordinates at which Li's enhancement switches off.  A substituted boundary point must
# not land on either, for the reason SUPERCRIT_RHOR gives for omitting 2.00: the two equations of
# state have slightly different Tc and rhoc, so a point placed on the gate is inside the window for
# one and outside it for the other.
GATE_RHOR, GATE_TR = 2.0, 1.4
GATE_KEEPOUT = 0.06


def heos(fluid):
    return "HEOS::" + fluid


def evaluation_caps(fluid):
    """Highest T and p, and lowest T, that EVERY backend the grid is evaluated on will accept.

    Returns (T_cap, p_cap, T_min).  A backend that cannot build the fluid at all does not
    constrain it: its points would fail for that reason whatever the grid did, and letting a
    failed lookup shrink the grid would silently punish every other backend for it.
    """
    T_cap = p_cap = float("inf")
    T_min = 0.0
    for backend in ("HEOS", "REFPROP"):
        name = "{}::{}".format(backend, fluid)
        try:
            Tmax = CP.PropsSI("Tmax", "T", 0, "p", 0, name)
            pmax = CP.PropsSI("pmax", "T", 0, "p", 0, name)
            Ttrip = CP.PropsSI("Ttriple", "T", 0, "p", 0, name)
        except Exception:
            continue
        T_cap = min(T_cap, CAP_MARGIN * T_EXTRAP_FACTOR * Tmax)
        p_cap = min(p_cap, CAP_MARGIN * P_EXTRAP_FACTOR * pmax)
        # REFPROP's triple point can sit ABOVE CoolProp's, which is the whole reason the liquid
        # block lost points at Tr = 0.6 for R1123 and R1224yd(Z).
        T_min = max(T_min, Ttrip)
    return T_cap, p_cap, T_min


def _isotherm(hp, T, rhoc, p_cap, label, capped):
    """One isotherm of the density ladder, clipped at p_cap, with the densest legal state added.

    Without the substitution a clipped isotherm ends early, and it ends early precisely at the
    dense end -- so the fluids whose EOS range is narrowest would contribute no dense evidence at
    all.
    """
    out = []
    reached = 0.0
    clipped = False
    for rr in SUPERCRIT_RHOR:
        try:
            p = CP.PropsSI("P", "T", T, "Dmass", rr * rhoc, hp)
        except Exception:
            continue
        if p <= 0:
            continue
        if p > p_cap:
            clipped = True
            continue
        out.append((T, p, "{}_rhor{}".format(label, rr)))
        reached = max(reached, rr)
    if not (clipped and out):
        return out
    try:
        rho_cap = CP.PropsSI("Dmass", "T", T, "P", p_cap, hp)
    except Exception:
        return out
    rr_cap = rho_cap / rhoc
    # Only worth a point if it is meaningfully denser than the last rung that fitted.
    if rr_cap <= 1.05 * reached:
        return out
    if abs(rr_cap - GATE_RHOR) < GATE_KEEPOUT:
        rr_cap = GATE_RHOR - GATE_KEEPOUT
        try:
            p_at = CP.PropsSI("P", "T", T, "Dmass", rr_cap * rhoc, hp)
        except Exception:
            return out
        if not (0 < p_at <= p_cap) or rr_cap <= 1.05 * reached:
            return out
        out.append((T, p_at, "{}_rhorCap".format(label)))
    else:
        out.append((T, p_cap, "{}_rhorCap".format(label)))
    capped["rhor"] += 1
    return out


def build_points(fluid, capped):
    """Return [(T_K, p_Pa, region)] for one pure fluid, placed in HEOS reduced coordinates."""
    hp = heos(fluid)
    Tc = CP.PropsSI("Tcrit", "T", 0, "p", 0, hp)
    rhoc = CP.PropsSI("rhocrit", "T", 0, "p", 0, hp)
    pcrit = CP.PropsSI("Pcrit", "T", 0, "p", 0, hp)
    T_cap, p_cap, T_min = evaluation_caps(fluid)
    pts = []

    reached_tr = 0.0
    clipped_tr = False
    for tr in SUPERCRIT_TR:
        T = tr * Tc
        if T > T_cap:
            clipped_tr = True
            continue
        pts.extend(_isotherm(hp, T, rhoc, p_cap, "sc_Tr{}".format(tr), capped))
        reached_tr = max(reached_tr, tr)

    # The hottest legal isotherm, for a fluid whose EOS range stops below the ladder.  Skipped when
    # it would sit on top of a rung that already fitted, or on the T/Tc = 1.4 enhancement gate.
    if clipped_tr:
        tr_cap = T_cap / Tc
        if abs(tr_cap - GATE_TR) < GATE_KEEPOUT:
            tr_cap = GATE_TR - GATE_KEEPOUT
        if tr_cap > 1.02 and tr_cap > 1.05 * reached_tr:
            pts.extend(_isotherm(hp, tr_cap * Tc, rhoc, p_cap, "sc_TrCap", capped))
            capped["tr"] += 1

    for tr in LIQUID_TR:
        T = tr * Tc
        if T <= T_min or T > T_cap:
            continue
        try:
            psat = CP.PropsSI("P", "T", T, "Q", 0, hp)
        except Exception:
            continue
        # A pressure barely above p_sat is still ambiguous in practice; push well past it.
        p = max(3.0 * psat, psat + 0.2 * pcrit)
        if p > p_cap:
            # Clamping to p_cap is only worth doing while it stays clear of the dome; otherwise
            # the point is exactly the ambiguous one this block exists to avoid.
            if p_cap < 1.5 * psat:
                continue
            p = p_cap
            capped["liq"] += 1
        pts.append((T, p, "liq_Tr{}".format(tr)))

    for tr in VAPOUR_TR:
        T = tr * Tc
        if T <= T_min or T > T_cap:
            continue
        try:
            psat = CP.PropsSI("P", "T", T, "Q", 1, hp)
        except Exception:
            continue
        pts.append((T, 0.3 * psat, "vap_Tr{}".format(tr)))

    return pts


def flash_is_consistent(names, z, T, p):
    """True when the PT flash lands on a state that corresponds to the pressure asked for.

    For a few mixtures it does not: we ask for the pressure at a given density, flash back at
    that (p, T), and land somewhere else entirely.  Nothing raises, because the state caches the
    pressure it was GIVEN rather than recomputing it, so p() agrees with the input by
    construction.

    Comparing the cached pressure against the EOS pressure at the returned density is what
    exposes it.  The identity p = rho*R*T*(1 + delta*alphar_delta) holds to ~1e-16 on a
    converged state.

    Only mixtures are affected: they are placed by reducing with a composition-averaged scale,
    which for a pair of very dissimilar components need not be near a physical critical density,
    so the requested states can sit where the flash struggles.
    """
    try:
        AS = CP.AbstractState("HEOS", "&".join(names))
        AS.set_mole_fractions(z)
        AS.update(CP.PT_INPUTS, p, T)
        d = AS.delta()
        p_alphar = AS.rhomolar() * AS.gas_constant() * AS.T() * (1 + d * AS.dalphar_dDelta())
        return abs(p_alphar / AS.p() - 1) < 1e-9
    except Exception:
        return False


def build_mixture_points(a, b, xa, xb):
    """Reduced-coordinate points for one binary, using HEOS's mixture REDUCING state.

    T_reducing / rhomolar_reducing, not T_critical: the mixture critical point costs a solve that
    can take ~half a second per call and does not always converge, and all this is used for is
    placing grid points -- any smooth, composition-aware scale does that job.
    """
    AS = CP.AbstractState("HEOS", "{}&{}".format(a, b))
    AS.set_mole_fractions([xa, xb])
    # HEOS has no density-temperature flash for mixtures ("DHSU_T_flash does not support
    # mixtures (yet)"), so the phase has to be imposed to reach the EOS directly.  That is not a
    # shortcut around an ambiguity: every mixture point below is supercritical by construction,
    # which is the whole reason mixtures are sampled only there.
    AS.specify_phase(CP.iphase_supercritical_gas)
    Tr = AS.T_reducing()
    rhor = AS.rhomolar_reducing() * AS.molar_mass()
    pts = []
    for tr in MIX_TR:
        T = tr * Tr
        for rr in MIX_RHOR:
            try:
                AS.update(CP.DmassT_INPUTS, rr * rhor, T)
                p = AS.p()
            except Exception:
                continue
            if p > 0 and flash_is_consistent([a, b], [xa, xb], T, p):
                pts.append((T, p, "mix_Tr{}_rhor{}".format(tr, rr)))
    AS.unspecify_phase()
    return pts


def main(fluids_arg):
    with open(os.path.join(SCRIPT_DIR, "res_transport_parameters.json"), encoding="utf-8") as fh:
        res = json.load(fh)

    if fluids_arg:
        fluids = [f.strip().upper() for f in fluids_arg.split(",") if f.strip()]
    else:
        # Fluids currently in the exclusion lists are included on purpose -- re-deriving those
        # exclusions is the point of this grid -- so the caller is expected to have run the
        # converter with --keep-all-heos.
        fluids = sorted(set(res["viscosity"]) & set(res["conductivity"]))

    os.makedirs(OUT_DIR, exist_ok=True)
    out_path = os.path.join(OUT_DIR, "grid_points.csv")
    n_fluid = 0
    n_pts = 0
    n_mix = 0
    capped = collections.Counter()
    skipped = []
    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["fluid", "T_K", "p_Pa", "region", "mole_fractions"])
        for fluid in fluids:
            try:
                pts = build_points(fluid, capped)
            except Exception as exc:
                skipped.append("{} ({})".format(fluid, str(exc).splitlines()[0][:60]))
                continue
            if not pts:
                skipped.append("{} (no usable points)".format(fluid))
                continue
            n_fluid += 1
            for T, p, region in pts:
                w.writerow([fluid, repr(T), repr(p), region, ""])
                n_pts += 1

        if not fluids_arg:
            for a, bb, xa, xb in MIXTURES:
                try:
                    pts = build_mixture_points(a, bb, xa, xb)
                except Exception as exc:
                    skipped.append("{}&{} ({})".format(a, bb, str(exc).splitlines()[0][:60]))
                    continue
                for T, p, region in pts:
                    w.writerow(["{}&{}".format(a, bb), repr(T), repr(p), region, "{};{}".format(xa, xb)])
                    n_pts += 1
                    n_mix += 1

    print("Wrote {}".format(out_path))
    print("  {} points over {} fluids ({} of them binary-mixture points)".format(n_pts, n_fluid, n_mix))
    print("  boundary states substituted for clipped rungs: {} density, {} isotherm, {} liquid"
          .format(capped["rhor"], capped["tr"], capped["liq"]))
    if skipped:
        print("  {} fluids skipped:".format(len(skipped)))
        for s in skipped:
            print("    {}".format(s))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fluids", default=None, help="comma-separated subset (default: every fluid in the JSON)")
    args = ap.parse_args()
    main(args.fluids)
