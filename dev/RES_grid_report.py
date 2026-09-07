"""Measure how well the REFPROP-fitted RES parameters transfer to CoolProp's HEOS backend.

Consumes the two CSVs written by dev/RES_reference_run.py, which evaluates the SAME code -- the
authors' own -- at the same (T, p) points on both equations of state.  One implementation, two
equations of state, so what is left is the equation of state and nothing else.  That is a
different and much sharper question than comparing against published values, which also carries
whatever an implementation gets wrong.

    python dev/RES_reference_run.py --eos REFPROP
    python dev/RES_reference_run.py --eos HEOS
    python dev/RES_grid_report.py

Output is the HEOS_TRANSFER_EXCLUDE set for dev/convert_RES_csv_to_json.py.

Two things make a naive comparison misleading, and both are handled here explicitly.

1. A transport deviation is only informative in proportion to how much of the property the
   residual term actually carries.  At a dilute-dominated state the dilute polynomial is identical
   on both sides, so the deviation is near zero however badly s_res disagrees.  Judging a fluid on
   such a point is what made the original single-sample-point exclusions unreliable.

   The dilute terms are not merely similar across the two equations of state, they are BIT-
   IDENTICAL -- the polynomial is a function of T alone, read from the same table, and measured
   over the whole grid the deviation between the two sides' eta0 and tc0 columns is exactly 0.
   So with s the residual share of the property, eta = eta0 + eta_res and s = eta_res/eta,

       d = |eta_H/eta_R - 1| = s * |eta_res,H/eta_res,R - 1| = s * e_res

   holds as an identity, not an approximation.  The criterion is therefore e_res = d / s: the
   error in the RESIDUAL TERM, which is the thing the coefficients actually determine, and which
   every point estimates on the same scale.

   Judging d directly instead would mean a different thing at every point: behind a share cutoff
   of 0.20 a 1% transport limit tests for a 5% residual error, at share 0.9 for 1.1%.  SHARE_FLOOR
   is only a numerical guard against dividing by a residual share indistinguishable from zero.  A
   fluid with no point above the floor is reported as UNJUDGED rather than silently passed.

2. Near the critical point the Olchowy-Sengers enhancement uses each side's OWN critical point
   and derivatives, so it moves for reasons that have nothing to do with the RES parameters.  Two
   criteria are reported side by side: NARROW judges only where the enhancement is inactive, so
   it isolates parameter transfer; BROAD judges everywhere, so it reflects what a user actually
   sees at a given (T, p).  NARROW drives the shipped list; BROAD is there so the cost of that
   choice stays visible.

The s_res clause is the second half of the criterion and is not redundant with the first: a fluid
can agree on the transport property at every sampled state while its residual entropy is badly
off, which means the sampled states happen to be insensitive rather than that the parameters
transfer.
"""

import argparse
import collections
import csv
import json
import math
import os
import statistics

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
COMP_DIR = os.path.join(SCRIPT_DIR, "RES_comparison")
REF_REFPROP = os.path.join(COMP_DIR, "grid_ref_refprop.csv")
REF_HEOS = os.path.join(COMP_DIR, "grid_ref_heos.csv")
JSON_PATH = os.path.join(SCRIPT_DIR, "res_transport_parameters.json")

# Numerical floor, NOT a physical criterion (see the module docstring).  e_res = d / share
# amplifies the deviation by 1/share, so a share indistinguishable from zero would amplify
# floating-point noise into a verdict.  At the floor the amplification is 50x, against a measured
# median d of 4e-7 -- five orders below DEV_LIMIT, so the floor costs nothing and guards the
# division.
SHARE_FLOOR = 0.02
THIN_POINTS = 5  # at or below this many judged points, say so -- do not let n hide in a table
DEV_LIMIT = 0.01  # transport deviation above which the parameters are judged not to transfer
SRES_LIMIT = 0.05  # s_res deviation above which the same conclusion is drawn directly


def enhancement_active(T, Tc, rho, rhoc):
    """Mirrors the gate Li 2024 applies, and CoolProp's RES conductivity with it."""
    return not (rho / rhoc >= 2.0 or T / Tc > 1.4)


def num(row, key):
    try:
        v = float(row.get(key, ""))
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def reldev(a, b):
    if a is None or b is None or b == 0:
        return None
    return abs(a / b - 1.0)


def share_of(total, dilute):
    """Fraction of `total` carried by the residual term, mapped into [0, 1).

    With eta = eta0 + eta_res this is exactly eta_res/eta, which is the factor by which an error
    in the residual term is attenuated in the property the user sees.  Dividing the observed
    deviation by it recovers the residual-term error itself.
    """
    if total is None or dilute is None or total == 0:
        return None
    x = abs(total / dilute - 1.0) if dilute else 1.0
    return x / (1.0 + x)


def read(path):
    if not os.path.exists(path):
        raise SystemExit("missing {} -- run dev/RES_reference_run.py for that equation of state first".format(path))
    with open(path, newline="", encoding="utf-8") as fh:
        return {(r["fluid"], r["T_K"], r["p_Pa"]): r for r in csv.DictReader(fh)}


def main(share_min, dev_limit, sres_limit, verbose, residual_criterion=True):
    rp_rows = read(REF_REFPROP)
    he_rows = read(REF_HEOS)
    residual_criterion = bool(residual_criterion)

    F = collections.defaultdict(lambda: collections.defaultdict(list))
    fail = collections.Counter()
    phase_mismatch = collections.Counter()
    gate_straddle = collections.Counter()

    # Seed F from every pure fluid in the grid, so a fluid whose points ALL fail still gets a
    # verdict row.  Without this it is simply absent from `result` -- not EXCLUDE, not UNJUDGED,
    # not NO-HEOS -- and silently keeps coefficients nothing ever measured.
    #
    # Spelled as an assignment rather than a bare `F[fluid]` subscript: the subscript does seed
    # the defaultdict, but it reads as a no-op and static analysers flag it as one.
    for key in rp_rows:
        fluid = key[0]
        if "&" not in fluid and fluid not in F:
            F[fluid] = collections.defaultdict(list)

    for key, rp in rp_rows.items():
        fluid = key[0]
        # Mixtures are skipped: the exclusion lists this drives are PER FLUID, so a binary has no
        # entry to withhold, and letting one in would colour a measurement about component
        # parameters with the mixture's own phase behaviour.  The mixing rules are checked
        # separately, against the C++ implementation rather than here.
        if "&" in fluid:
            continue
        he = he_rows.get(key)
        if rp["ok"] != "1":
            fail[(fluid, "REFPROP")] += 1
            continue
        if he is None or he["ok"] != "1":
            fail[(fluid, "HEOS")] += 1
            continue
        if num(rp, "phase") != num(he, "phase"):
            phase_mismatch[fluid] += 1

        T = float(key[1])
        rho, rhoc, Tc = num(rp, "rho"), num(rp, "rhoc"), num(rp, "Tc")
        rho_h, rhoc_h, Tc_h = num(he, "rho"), num(he, "rhoc"), num(he, "Tc")
        if None in (rho, rhoc, Tc, rho_h, rhoc_h, Tc_h) or 0 in (rhoc, Tc, rhoc_h, Tc_h):
            continue
        enh = enhancement_active(T, Tc, rho, rhoc)
        # Discard points where the two sides land on opposite sides of the enhancement gate.
        # Their rhoc differs slightly, so a point near rho/rhoc = 2 can have the enhancement on
        # for one and off for the other; the resulting step says nothing about the parameters.
        # The grid avoids the boundary, but the check stays -- a future grid could wander back.
        if enh != enhancement_active(T, Tc_h, rho_h, rhoc_h):
            gate_straddle[fluid] += 1
            continue

        sv = share_of(num(rp, "eta"), num(rp, "eta0"))
        sc = share_of(num(rp, "tc"), num(rp, "tc0"))
        d_eta = reldev(num(he, "eta"), num(rp, "eta"))
        d_tc = reldev(num(he, "tc"), num(rp, "tc"))
        d_s = reldev(num(he, "s_res"), num(rp, "s_res"))

        if d_eta is not None and sv is not None:
            F[fluid]["eta"].append((sv, d_eta, enh, d_s))
        if d_tc is not None and sc is not None:
            F[fluid]["tc"].append((sc, d_tc, enh, d_s))

    def judged(pairs, narrow):
        # Dividing by the residual share turns the observed transport deviation into the error in
        # the residual term itself.  The legacy branch does not divide; see the module docstring.
        sel = [((d / s if residual_criterion else d), ds)
               for s, d, enh, ds in pairs if s >= share_min and (not narrow or not enh)]
        if not sel:
            return None, None, 0
        ds_vals = [ds for _, ds in sel if ds is not None]
        # None here means "s_res could not be compared at any judged point", which is NOT the same
        # as "s_res agrees".  Returning it as None lets the caller report the fluid rather than
        # letting the criterion quietly collapse to its transport clause alone.
        return max(d for d, _ in sel), (max(ds_vals) if ds_vals else None), len(sel)

    with open(JSON_PATH, encoding="utf-8") as fh:
        shipped = json.load(fh)

    result = {}
    for fluid in sorted(F):
        for prop, key in (("eta", "viscosity"), ("tc", "conductivity")):
            pairs = F[fluid][prop]
            wn, dsn, nn = judged(pairs, True)
            wb, dsb, nb = judged(pairs, False)
            if fail[(fluid, "HEOS")] and not pairs:
                verdict = ("NO-HEOS", "NO-HEOS")
            else:
                # BOTH clauses, on both criteria.  NARROW and BROAD are meant to differ in exactly
                # one way -- whether points inside the enhancement window are judged -- so that the
                # gap between them measures the enhancement and nothing else.  Dropping the s_res
                # clause from BROAD alone would put a second, undocumented difference in that gap.
                def over(w, ds):
                    return (w is not None and w > dev_limit) or (ds is not None and ds > sres_limit)
                # Point COUNT does not enter the verdict.  It was tempting to demote a `keep` on
                # few points to UNJUDGED, but the fluids that would hit -- CYCLOPRO and RE143A,
                # both n=2 -- agree to 0.001%, four orders below the limit.  Calling that
                # "unjudged" overstates the doubt as badly as hiding it would understate it.  The
                # count is reported instead, in the table and in the thin-evidence list below, so
                # a reader can weigh it.
                verdict = (
                    "UNJUDGED" if nn == 0 else ("EXCLUDE" if over(wn, dsn) else "keep"),
                    "UNJUDGED" if nb == 0 else ("EXCLUDE" if over(wb, dsb) else "keep"),
                )
            result[(fluid, key)] = dict(narrow=verdict[0], broad=verdict[1], wn=wn, wb=wb, nn=nn, nb=nb, ds=dsn)

    # Name the column for what it holds.  Residual-term errors are LARGER than the transport
    # deviations under --legacy-share-criterion, and reading one as the other would badly misjudge
    # how close a fluid sits to the limit.
    col = "e_res" if residual_criterion else "dev"
    print("criterion: {}\n".format(
        "e_res = d / share, floor {:g} (residual-term error)".format(share_min) if residual_criterion
        else "d, share cutoff {:g} (legacy transport deviation)".format(share_min)))
    hdr = "{:<15} {:<4} {:>11} {:>4} {:>11} {:>4} {:>11}  narrow/broad".format(
        "fluid", "prop", "narrow " + col, "n", "broad " + col, "n", "max ds_res")
    print(hdr)
    print("-" * len(hdr))
    for (fluid, key), r in sorted(result.items()):
        if not verbose and r["narrow"] == "keep" and r["broad"] == "keep":
            continue
        def f(v):
            return "" if v is None else "{:10.3f}%".format(v * 100)
        print("{:<15} {:<4} {:>11} {:>4} {:>11} {:>4} {:>11}  {}/{}".format(
            fluid, key[:4], f(r["wn"]), r["nn"], f(r["wb"]), r["nb"], f(r["ds"]), r["narrow"], r["broad"]))

    print("\n" + "=" * 78)
    excluded_any = set()
    for key in ("viscosity", "conductivity"):
        cur = {n for n, e in shipped[key].items() if "HEOS" not in e}
        for crit in ("narrow", "broad"):
            ex = sorted(f for (f, k), r in result.items() if k == key and r[crit] == "EXCLUDE")
            print("{:<13} {:<7} {:>3} fluids".format(key.upper(), crit.upper(), len(ex)))
            print("  {{{}}}".format(", ".join('"' + f + '"' for f in ex)))
            if crit == "narrow":
                excluded_any |= set(ex)
        un = sorted(f for (f, k), r in result.items() if k == key and r["narrow"] == "UNJUDGED")
        if un:
            print("  unjudgeable (no residual-dominated point): {}".format(", ".join(un)))
        nh = sorted(f for (f, k), r in result.items() if k == key and r["narrow"] == "NO-HEOS")
        if nh:
            print("  not evaluable on HEOS at all: {}".format(", ".join(nh)))
        thin = sorted("{}(n={})".format(f, r["nn"]) for (f, k), r in result.items()
                      if k == key and r["narrow"] == "keep" and r["nn"] <= THIN_POINTS)
        if thin:
            print("  kept on few judged points: {}".format(", ".join(thin)))
        no_s = sorted(f for (f, k), r in result.items() if k == key and r["nn"] > 0 and r["ds"] is None)
        if no_s:
            print("  s_res clause could not be evaluated (transport clause only): {}".format(", ".join(no_s)))
        print("  currently shipped as excluded: {}".format(sorted(cur)))
        print()

    print("HEOS_TRANSFER_EXCLUDE = {")
    for f in sorted(excluded_any):
        print('    "{}",'.format(f))
    print("}")

    def med(xs):
        return statistics.median(xs) if xs else float("nan")

    all_eta = [d for f in F for _s, d, _e, _ds in F[f]["eta"]]
    all_tc = [d for f in F for _s, d, _e, _ds in F[f]["tc"]]
    print("\nmedian TRANSPORT deviation over all {} points: viscosity {:.5f} %, conductivity {:.5f} %".format(
        len(all_eta), med(all_eta) * 100, med(all_tc) * 100))
    print("fluids with a HEOS evaluation failure: {}".format(len({f for (f, b) in fail if b == "HEOS"})))
    print("fluids with a REFPROP evaluation failure: {}".format(len({f for (f, b) in fail if b == "REFPROP"})))
    print("points where the two sides disagreed on phase: {}".format(sum(phase_mismatch.values())))
    print("points discarded for straddling the enhancement gate: {}".format(sum(gate_straddle.values())))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--share-min", type=float, default=SHARE_FLOOR,
                    help="numerical floor on the residual share (legacy: the cutoff below which a "
                         "point was discarded)")
    ap.add_argument("--legacy-share-criterion", action="store_true",
                    help="judge the raw transport deviation against a hard share cutoff instead of "
                         "the residual-term error; with --share-min 0.2 this is the pre-2026-09 "
                         "criterion, for reproducing an older run")
    ap.add_argument("--dev-limit", type=float, default=DEV_LIMIT)
    ap.add_argument("--sres-limit", type=float, default=SRES_LIMIT)
    ap.add_argument("-v", "--verbose", action="store_true", help="list every fluid, not just the interesting ones")
    a = ap.parse_args()
    main(a.share_min, a.dev_limit, a.sres_limit, a.verbose, not a.legacy_share_criterion)
