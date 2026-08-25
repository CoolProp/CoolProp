"""Stage-6(b): measure how well the REFPROP-fitted RES parameters transfer to HEOS.

Consumes dev/RES_comparison/grid_cpp.csv, written by the [RES_grid] Catch2 harness, which
evaluates the SAME implementation on both backends at the same (T, p) points.  One implementation,
two equations of state, so what is left is the equation of state and nothing else -- unlike a
comparison against published values, which also carries whatever the implementation gets wrong.

Two things make a naive comparison misleading, and both are handled here explicitly.

1. A transport deviation is only informative where the property actually depends on the residual
   term.  At a dilute-dominated state the dilute polynomial is identical on both backends by
   construction, so the deviation is near zero however badly s_res disagrees.  Judging a fluid on
   such a point is what made the original single-sample-point exclusions unreliable.  Points are
   therefore weighted by residual share and fluids with no residual-dominated point are reported
   as unjudgeable rather than silently passed.

2. Near the critical point the Olchowy-Sengers enhancement uses each backend's OWN critical point
   and derivatives, so it moves for reasons that have nothing to do with the RES parameters.  Two
   criteria are reported side by side: NARROW judges only where the enhancement is inactive, so it
   isolates parameter transfer; BROAD judges everywhere, so it reflects what a user actually sees
   at a given (T, p).  They disagree, and which one should drive the shipped exclusion lists is a
   policy decision, not a measurement.

The [RES_grid] harness pins both backends to the same dilute-gas and enhancement-viscosity source,
without which this measures the Stage-4 source policy rather than the parameters.
"""

import argparse
import collections
import csv
import json
import math
import os
import statistics

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
GRID_CSV = os.path.join(SCRIPT_DIR, "RES_comparison", "grid_cpp.csv")
JSON_PATH = os.path.join(SCRIPT_DIR, "res_transport_parameters.json")

SHARE_MIN = 0.20  # a point counts only if this much of the property is non-dilute
DEV_LIMIT = 0.01  # deviation above which the parameters are judged not to transfer

# Mirrors the gate in RESTransport::conductivity, which follows Li 2024.
def enhancement_active(T, Tc, rho, rhoc):
    return not (rho / rhoc >= 2.0 or T / Tc > 1.4)


def num(r, k):
    try:
        f = float(r.get(k, ""))
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def reldev(a, b):
    if a is None or b is None or b == 0:
        return None
    return abs(a / b - 1.0)


def share_of(total, dilute):
    """Fraction of `total` that is not the dilute term, mapped into [0, 1)."""
    if total is None or dilute is None or total == 0:
        return None
    x = abs(total / dilute - 1.0) if dilute else 1.0
    return x / (1.0 + x)


def main(share_min, dev_limit, verbose):
    rows = collections.defaultdict(dict)
    with open(GRID_CSV, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            rows[(r["fluid"], r["T_K"], r["p_Pa"])][r["backend"]] = r

    F = collections.defaultdict(lambda: collections.defaultdict(list))
    fail = collections.Counter()
    phase_mismatch = collections.Counter()
    gate_straddle = collections.Counter()
    seen = collections.Counter()

    for (fluid, T_s, _p), byb in rows.items():
        # REFPROP_pinned, not REFPROP: the pinned column uses the same dilute term and
        # enhancement viscosity HEOS can, leaving the equation of state as the only difference.
        rp, he = byb.get("REFPROP_pinned"), byb.get("HEOS")
        seen[fluid] += 1
        if rp is None or rp["ok"] != "1":
            fail[(fluid, "REFPROP_pinned")] += 1
            continue
        if he is None or he["ok"] != "1":
            fail[(fluid, "HEOS")] += 1
            continue
        if num(rp, "phase") != num(he, "phase"):
            phase_mismatch[fluid] += 1

        T, rho, rhoc, Tc = float(T_s), num(rp, "rho"), num(rp, "rhoc"), num(rp, "Tc")
        rho_h, rhoc_h, Tc_h = num(he, "rho"), num(he, "rhoc"), num(he, "Tc")
        if None in (rho, rhoc, Tc, rho_h, rhoc_h, Tc_h) or 0 in (rhoc, Tc, rhoc_h, Tc_h):
            continue
        enh = enhancement_active(T, Tc, rho, rhoc)
        # Discard points where the two backends land on opposite sides of the enhancement gate.
        # Their rhoc differs slightly, so a point near rho/rhoc = 2 can have the enhancement on for
        # one and off for the other; the resulting step says nothing about the parameters.  The
        # grid avoids the boundary, but the check stays -- a future grid could wander back onto it.
        if enh != enhancement_active(T, Tc_h, rho_h, rhoc_h):
            gate_straddle[fluid] += 1
            continue

        # The dilute term the model actually subtracted, so the share is honest about it.
        eta0 = num(rp, "eta0_native")
        if eta0 is None:
            eta0 = num(rp, "eta0_poly")
        sv = share_of(num(rp, "eta"), eta0)
        sc = share_of(num(rp, "tc"), num(rp, "tc0_poly"))
        d_eta = reldev(num(he, "eta"), num(rp, "eta"))
        d_tc = reldev(num(he, "tc"), num(rp, "tc"))
        d_s = reldev(num(he, "s_res"), num(rp, "s_res"))

        if d_s is not None:
            F[fluid]["sres"].append(d_s)
        if d_eta is not None and sv is not None:
            F[fluid]["eta"].append((sv, d_eta, enh))
        if d_tc is not None and sc is not None:
            F[fluid]["tc"].append((sc, d_tc, enh))

    def worst(pairs, narrow):
        sel = [d for s, d, enh in pairs if s >= share_min and (not narrow or not enh)]
        return (max(sel) if sel else None), len(sel)

    with open(JSON_PATH, encoding="utf-8") as fh:
        shipped = json.load(fh)

    result = {}
    for fluid in sorted(F):
        for prop, key in (("eta", "viscosity"), ("tc", "conductivity")):
            pairs = F[fluid][prop]
            wn, nn = worst(pairs, True)
            wb, nb = worst(pairs, False)
            if fail[(fluid, "HEOS")] and not pairs:
                verdict = ("NO-HEOS", "NO-HEOS")
            else:
                verdict = (
                    "UNJUDGED" if nn == 0 else ("EXCLUDE" if wn > dev_limit else "keep"),
                    "UNJUDGED" if nb == 0 else ("EXCLUDE" if wb > dev_limit else "keep"),
                )
            result[(fluid, key)] = dict(narrow=verdict[0], broad=verdict[1], wn=wn, wb=wb, nn=nn, nb=nb,
                                        ds=max(F[fluid]["sres"]) if F[fluid]["sres"] else None)

    hdr = f"{'fluid':<15} {'prop':<4} {'narrow dev':>11} {'n':>4} {'broad dev':>11} {'n':>4} {'max ds_res':>11}  narrow/broad"
    print(hdr)
    print("-" * len(hdr))
    for (fluid, key), r in sorted(result.items()):
        if not verbose and r["narrow"] == "keep" and r["broad"] == "keep":
            continue
        f = lambda v: "" if v is None else f"{v * 100:10.3f}%"
        print(f"{fluid:<15} {key[:4]:<4} {f(r['wn']):>11} {r['nn']:>4} {f(r['wb']):>11} {r['nb']:>4} "
              f"{f(r['ds']):>11}  {r['narrow']}/{r['broad']}")

    print("\n" + "=" * 78)
    for key in ("viscosity", "conductivity"):
        cur = {n for n, e in shipped[key].items() if "HEOS" not in e and ("PR" in e or "SRK" in e)}
        for crit in ("narrow", "broad"):
            ex = sorted(f for (f_, k), r in result.items() if k == key and r[crit] == "EXCLUDE" for f in [f_])
            print(f"{key.upper():<13} {crit.upper():<7} {len(ex):>3} fluids")
            print(f"  {{{', '.join(chr(34) + f + chr(34) for f in ex)}}}")
        un = sorted(f for (f, k), r in result.items() if k == key and r["narrow"] == "UNJUDGED")
        if un:
            print(f"  unjudgeable (no residual-dominated point): {', '.join(un)}")
        print(f"  currently shipped as excluded: {sorted(cur)}")
        print()

    med = lambda xs: statistics.median(xs) if xs else float("nan")
    all_eta = [d for f in F for s, d, e in F[f]["eta"]]
    all_tc = [d for f in F for s, d, e in F[f]["tc"]]
    print(f"median deviation over all {len(all_eta)} points: viscosity {med(all_eta) * 100:.5f} %, "
          f"conductivity {med(all_tc) * 100:.5f} %")
    print(f"fluids with a HEOS evaluation failure: {len({f for (f, b) in fail if b == 'HEOS'})}")
    print(f"points where the backends disagreed on phase: {sum(phase_mismatch.values())}")
    print(f"points discarded for straddling the enhancement gate: {sum(gate_straddle.values())}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--share-min", type=float, default=SHARE_MIN)
    ap.add_argument("--dev-limit", type=float, default=DEV_LIMIT)
    ap.add_argument("-v", "--verbose", action="store_true")
    a = ap.parse_args()
    main(a.share_min, a.dev_limit, a.verbose)
