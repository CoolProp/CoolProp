"""Stage-6(a): check the C++ RES implementation against the authors' own reference code.

Runs Martinek 2025 (viscosity) and Li 2024 (conductivity) over the same grid the C++ REFPROP+RES
harness used, and compares point by point.  Both sides evaluate on REFPROP, so the equation of
state is identical and any deviation is an implementation defect -- this is the half of Stage 6
that says whether the model is coded correctly, separate from whether its parameters transfer.

The published sample sets already give one point per fluid (the [RES] parity tests cover those).
This extends the comparison to states those points never reach: dense supercritical, compressed
liquid, and the near-critical region where the Olchowy-Sengers enhancement is active.

The reference code needs its own data files in the working directory, so point --ref-dir at a
copy of each SI code directory:

    python dev/RES_reference_check.py --vis-dir <...>/Martinek_2025_SI/.../code_SI \\
                                      --tc-dir  <...>/Li_2024_SI/code_SI

Both modules call pandas with delim_whitespace=, removed in pandas 2.2, so it is shimmed here
rather than editing the authors' source.
"""

import argparse
import collections
import csv
import importlib.util
import math
import os
import random
import statistics
import sys

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
GRID_CPP = os.path.join(SCRIPT_DIR, "RES_comparison", "grid_cpp.csv")


def shim_pandas():
    import pandas as pd

    orig = pd.read_csv

    def read_csv(*a, **k):
        if k.pop("delim_whitespace", None):
            k["sep"] = r"\s+"
        return orig(*a, **k)

    pd.read_csv = read_csv


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def read_cpp():
    out = {}
    with open(GRID_CPP, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            # The defaults column: the papers take the dilute term and the enhancement
            # viscosity from REFPROP's own transport model, so only this column is comparable.
            if r["backend"] != "REFPROP" or r["ok"] != "1":
                continue
            out[(r["fluid"], r["T_K"], r["p_Pa"])] = r
    return out


def main(vis_dir, tc_dir, per_fluid, seed, dev_limit):
    shim_pandas()
    cpp = read_cpp()
    by_fluid = collections.defaultdict(list)
    for k in cpp:
        by_fluid[k[0]].append(k)

    rng = random.Random(seed)
    selected = []
    for fluid, keys in sorted(by_fluid.items()):
        keys.sort()
        selected += keys if len(keys) <= per_fluid else rng.sample(keys, per_fluid)

    results = {"viscosity": [], "conductivity": []}
    failures = collections.Counter()
    cwd = os.getcwd()

    # Both reference functions return SI (Pa.s, W/m/K).  main.py scales to microPa.s only when
    # it writes its sample table, which is why the published column is in microPa.s and this is not.
    for label, path, prop, col in (
        ("Martinek", vis_dir, "viscosity", "eta"),
        ("Li", tc_dir, "conductivity", "tc"),
    ):
        if not path:
            print(f"-- {label}: skipped, no directory given")
            continue
        os.chdir(path)
        try:
            mod = load_module("ref_" + prop, os.path.join(path, "main.py" if prop == "viscosity" else "code_SI.py"))
            fn = mod.viscosity_RES if prop == "viscosity" else mod.TC_RES
            for key in selected:
                fluid, T_s, p_s = key
                try:
                    # Both modules return the value as a 1-element numpy array for pure fluids.
                    import numpy as np

                    val = float(np.ravel(fn(fluid, [1.0], float(p_s), float(T_s))[0])[0])
                except Exception:
                    failures[(label, "ref")] += 1
                    continue
                try:
                    ours = float(cpp[key][col])
                except (TypeError, ValueError):
                    failures[(label, "cpp")] += 1
                    continue
                if not (math.isfinite(val) and math.isfinite(ours)) or val == 0:
                    failures[(label, "nonfinite")] += 1
                    continue
                results[prop].append((fluid, cpp[key]["region"], abs(ours / val - 1.0)))
        finally:
            os.chdir(cwd)

    print(f"\ngrid points compared: {len(selected)} (up to {per_fluid} per fluid, seed {seed})")
    for prop in ("viscosity", "conductivity"):
        rows = results[prop]
        if not rows:
            continue
        devs = [d for _f, _r, d in rows]
        rows.sort(key=lambda t: -t[2])
        print(f"\n{prop.upper()}: n={len(rows)}  median={statistics.median(devs) * 100:.6f}%  "
              f"max={max(devs) * 100:.4f}%")
        over = [t for t in rows if t[2] > dev_limit]
        print(f"  points above {dev_limit * 100:g}%: {len(over)}")
        # Near the critical point the enhancement contains 1/(dp/drho), which is singular there,
        # so rho/rhoc = 1 at T/Tc = 1.02 amplifies any difference without indicating a defect.
        # Report the two populations separately rather than letting the singular one set the max.
        crit = [t for t in rows if "Tr1.02_rhor1.0" in t[1] or "Tr1.02_rhor1.5" in t[1]]
        rest = [t for t in rows if t not in crit]
        if crit and rest:
            print(f"  excluding the {len(crit)} near-critical points: "
                  f"median={statistics.median([d for _f, _r, d in rest]) * 100:.6f}%  "
                  f"max={max(d for _f, _r, d in rest) * 100:.4f}%")
        for fluid, region, d in rows[:12]:
            print(f"    {fluid:<14} {region:<22} {d * 100:9.4f}%")
    if failures:
        print("\nnot compared:", dict(failures))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vis-dir", default=None, help="Martinek 2025 code_SI directory (contains main.py)")
    ap.add_argument("--tc-dir", default=None, help="Li 2024 code_SI directory (contains code_SI.py)")
    ap.add_argument("--per-fluid", type=int, default=4, help="grid points sampled per fluid")
    ap.add_argument("--seed", type=int, default=20260824)
    ap.add_argument("--dev-limit", type=float, default=0.001)
    a = ap.parse_args()
    main(a.vis_dir, a.tc_dir, a.per_fluid, a.seed, a.dev_limit)
