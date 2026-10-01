#!/usr/bin/env python3
"""Score CoolProp (one bench_ratio CSV per build) and REFPROP TPFLSH against the REFPROP-only per-state
reference from rp_truth_n.cpp, at EVERY state (not only where the codes disagree).

    python3 score_all.py <out_prefix> <build_label>:<dir_with_csvs> ... -- <truth_dir>

For mixture k the CSV is <dir>/<tag>_<k>.csv and the reference <truth_dir>/truth_<k>.txt.
Writes <out_prefix>_<build_label>.csv (plot verdicts: CoolProp markers, plus TPFLSH markers for the
first build only) and prints a per-mixture table.

Categories (density tolerances are relative, on the bulk molar density):
  ok              phase right and |rho/rho_ref - 1| <= 1e-4
  loose           phase right, 1e-4 < dev <= 1e-3                   -> cp_loose / rp_minor
  label only      two-phase reported, reference single phase, but density within 1e-4 -> cp_loose / rp_minor
  wrong           false split with wrong density, or dev > 1e-3     -> cp_rho / rp_wrong
  missed split    reference has a split (tm < -1e-8), code says single phase -> cp_wrong / rp_scope
  error / no ref  not marked
"""
import collections
import csv
import glob
import os
import sys

CATS = ["ok", "loose", "label only", "wrong", "missed split", "error", "no reference"]
CP_KEY = {"loose": "cp_loose", "label only": "cp_loose", "wrong": "cp_rho", "missed split": "cp_wrong"}
RP_KEY = {"loose": "rp_minor", "label only": "rp_minor", "wrong": "rp_wrong", "missed split": "rp_scope"}


def score(two, rho, err, ref):
    if err:
        return "error"
    if ref is None:
        return "no reference"
    tph, tr = ref
    if tph == "2ph" and not two:
        return "missed split"
    if tr is None:
        return "no reference"
    d = abs(rho / tr - 1)
    if tph == "1ph" and two:
        return "label only" if d <= 1e-4 else "wrong"
    return "ok" if d <= 1e-4 else "loose" if d <= 1e-3 else "wrong"


def main():
    args = sys.argv[1:]
    out = args[0]
    sep = args.index("--")
    builds = [a.split(":", 1) for a in args[1:sep]]
    tdir = args[sep + 1]
    for bi, (label, d) in enumerate(builds):
        rows = []
        print(f"\n=== {label}")
        print("%-28s %-9s" % ("mixture", "code") + "".join("%13s" % c for c in CATS))
        for fn in sorted(glob.glob(os.path.join(d, "*.csv"))):
            k = os.path.basename(fn).split("_")[-1].split(".")[0]
            tf = os.path.join(tdir, f"truth_{k}.txt")
            if not os.path.exists(tf):
                continue
            ref = {}
            for line in open(tf):
                f = line.split()
                if len(f) < 3 or f[1] == "NOREF":
                    continue
                r = float(f[2])
                ref[f[0]] = (f[1], r if r > 0 else None)
            R = list(csv.DictReader(open(fn)))
            name = R[0]["mixture"]
            cc, rc = collections.Counter(), collections.Counter()
            for r in R:
                t = ref.get(r["i"])
                c = score(0 < float(r["Q_cp"]) < 1, float(r["rho_cp"]), r["cp_fail"] == "1", t)
                cc[c] += 1
                if c in CP_KEY:
                    rows.append([name, r["i"], CP_KEY[c], c])
                q = float(r["q_rp"])
                c = score(0 <= q <= 1, float(r["rho_rp"]), int(r["ierr_rp"]) > 0, t)
                rc[c] += 1
                if bi == 0 and c in RP_KEY:
                    rows.append([name, r["i"], RP_KEY[c], c])
            print("%-28s %-9s" % (name[:28], "CoolProp") + "".join("%13d" % cc[c] for c in CATS))
            if bi == 0:
                print("%-28s %-9s" % ("", "TPFLSH") + "".join("%13d" % rc[c] for c in CATS))
        with open(f"{out}_{label}.csv", "w", newline="") as o:
            w = csv.writer(o)
            w.writerow(["mixture", "i", "verdict", "detail"])
            w.writerows(rows)


if __name__ == "__main__":
    main()
