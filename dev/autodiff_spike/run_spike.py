#!/usr/bin/env python3
"""
Drive the AD spike end to end:
  1. compile-time / object-size sweep per method, K = 1, 4, 16 model instantiations
  2. build + run bench (values and timings)
  3. 50-digit mpmath reference for every derivative, independent pure-Python PC-SAFT
  4. write results.json and print markdown tables

Usage: python3 run_spike.py [--outdir DIR] [--skip-compile-sweep]
Needs: a C++20 compiler (c++), mpmath; autodiff and mcx headers (defaults point at a
teqp checkout; override with AUTODIFF_DIR / MCX_DIR).
"""
import argparse
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).resolve().parent
METHODS = ["double", "ad_real", "ad_dual", "numdual", "mcx", "cstep", "fd"]
AUTODIFF_DIR = os.environ.get("AUTODIFF_DIR", str(Path.home() / "Code/teqp/externals/autodiff"))
MCX_DIR = os.environ.get("MCX_DIR", str(Path.home() / "Code/teqp/externals/mcx/multicomplex/include/multicomplex"))
CXX = os.environ.get("CXX", "c++")
INC = [f"-I{HERE}", f"-I{AUTODIFF_DIR}", f"-I{MCX_DIR}"]


def compile_obj(method, K, opt, out):
    cmd = [CXX, "-std=c++20", opt, *INC, f"-DKMODELS={K}", "-c", str(HERE / "methods" / f"m_{method}.cpp"), "-o", str(out)]
    t0 = time.perf_counter()
    r = subprocess.run(cmd, capture_output=True, text=True)
    dt = time.perf_counter() - t0
    if r.returncode != 0:
        sys.exit(f"compile failed: {' '.join(cmd)}\n{r.stderr}")
    return dt


def text_size(obj):
    # __TEXT,__text section size in bytes (machine code only)
    r = subprocess.run(["size", "-m", str(obj)], capture_output=True, text=True, check=True).stdout
    m = re.search(r"Section \(__TEXT, __text\): (\d+)", r)
    if m is None:
        sys.exit(f"could not parse `size -m` output for {obj}:\n{r}")
    return int(m.group(1))


def compile_sweep(outdir, repeats=3):
    res = []
    configs = [("-O2", 1), ("-O2", 4), ("-O2", 16), ("-O0", 1)]
    for method in METHODS:
        for opt, K in configs:
            obj = outdir / f"sweep_{method}_{opt[1:]}_{K}.o"
            times = [compile_obj(method, K, opt, obj) for _ in range(repeats)]
            res.append(dict(method=method, opt=opt, K=K, seconds=min(times), text_bytes=text_size(obj),
                            file_bytes=obj.stat().st_size))
            print(f"  {method:8s} {opt} K={K:2d}  {min(times):6.2f} s  text={res[-1]['text_bytes']:>9d} B", flush=True)
    return res


def run_bench(outdir):
    objs = []
    for method in METHODS:
        o = outdir / f"bench_{method}.o"
        compile_obj(method, 1, "-O2", o)
        objs.append(str(o))
    exe = outdir / "bench"
    subprocess.run([CXX, "-std=c++20", "-O2", *INC, str(HERE / "bench.cpp"), *objs, "-o", str(exe)], check=True)
    vals, ns = {}, {}
    for line in subprocess.run([str(exe)], capture_output=True, text=True, check=True).stdout.splitlines()[1:]:
        kind, method, state, q, idx, v = line.split(",")
        if kind == "val":
            vals[(method, state, q, int(idx))] = float(v)
        else:
            ns[(method, state, q)] = float(v)
    return vals, ns


# ------------------------------------------------------------------ mpmath reference
GS_A = [[0.9105631445, 0.6361281449, 2.6861347891, -26.547362491, 97.759208784, -159.59154087, 91.297774084],
        [-0.3084016918, 0.1860531159, -2.5030047259, 21.419793629, -65.255885330, 83.318680481, -33.746922930],
        [-0.0906148351, 0.4527842806, 0.5962700728, -1.7241829131, -4.1302112531, 13.776631870, -8.6728470368]]
GS_B = [[0.7240946941, 2.2382791861, -4.0025849485, -21.003576815, 26.855641363, 206.55133841, -355.60235612],
        [-0.5755498075, 0.6995095521, 3.8925673390, -17.215471648, 192.67226447, -161.82646165, -165.20769346],
        [0.0976883116, -0.2557574982, -9.1558561530, 20.642075974, -38.804430052, 93.626774077, -29.666905585]]


def alphar_mp(p, T, rho, x):
    # Constants enter as the exact binary64 values the C++ code sees (mpf(float) is exact).
    m, sig, eps = [[mp.mpf(v) for v in p[k]] for k in ("m", "sigma", "eps")]
    N = len(m)
    NA = mp.mpf(6.02214076e23) * mp.mpf(1e-30)
    pi = mp.mpf(3.14159265358979323846)
    d = [sig[i] * (1 - mp.mpf(0.12) * mp.exp(-3 * eps[i] / T)) for i in range(N)]
    rhoN = rho * NA
    z = [pi / 6 * rhoN * sum(x[i] * m[i] * d[i] ** n for i in range(N)) for n in range(4)]
    mbar = sum(x[i] * m[i] for i in range(N))
    om = 1 - z[3]
    ahs = (3 * z[1] * z[2] / om + z[2] ** 3 / (z[3] * om ** 2) + (z[2] ** 3 / z[3] ** 2 - z[0]) * mp.log(om)) / z[0]
    ahc = mbar * ahs
    for i in range(N):
        dij = d[i] / 2
        g = 1 / om + dij * 3 * z[2] / om ** 2 + dij ** 2 * 2 * z[2] ** 2 / om ** 3
        ahc -= x[i] * (m[i] - 1) * mp.log(g)
    eta = z[3]
    c1, c2 = (mbar - 1) / mbar, (mbar - 1) / mbar * (mbar - 2) / mbar
    I1 = sum((mp.mpf(GS_A[0][i]) + c1 * mp.mpf(GS_A[1][i]) + c2 * mp.mpf(GS_A[2][i])) * eta ** i for i in range(7))
    I2 = sum((mp.mpf(GS_B[0][i]) + c1 * mp.mpf(GS_B[1][i]) + c2 * mp.mpf(GS_B[2][i])) * eta ** i for i in range(7))
    s1 = s2 = 0
    for i in range(N):
        for j in range(N):
            sij = (sig[i] + sig[j]) / 2
            eT = mp.sqrt(eps[i] * eps[j]) / T
            pre = x[i] * x[j] * m[i] * m[j] * sij ** 3
            s1 += pre * eT
            s2 += pre * eT ** 2
    C1 = 1 / (1 + mbar * (8 * eta - 2 * eta ** 2) / om ** 4
              + (1 - mbar) * (20 * eta - 27 * eta ** 2 + 12 * eta ** 3 - 2 * eta ** 4) / (om * (2 - eta)) ** 2)
    return ahc - 2 * pi * rhoN * I1 * s1 - pi * rhoN * mbar * C1 * I2 * s2


STATES = {
    "propane_liq": (dict(m=[2.0020], sigma=[3.6184], eps=[208.11]), 300.0, 11000.0, [1.0]),
    "propane_gas": (dict(m=[2.0020], sigma=[3.6184], eps=[208.11]), 300.0, 100.0, [1.0]),
    "C1C2C3_dense": (dict(m=[1.0, 1.6069, 2.0020], sigma=[3.7039, 3.5206, 3.6184], eps=[150.03, 191.42, 208.11]),
                     250.0, 8000.0, [0.5, 0.3, 0.2]),
}


def reference():
    mp.mp.dps = 60
    ref = {}
    for name, (p, T, rho, x) in STATES.items():
        T, rho, x = mp.mpf(T), mp.mpf(rho), [mp.mpf(v) for v in x]
        ti = 1 / T
        fr = lambda r: alphar_mp(p, T, r, x)
        ft = lambda t: alphar_mp(p, 1 / t, rho, x)
        for n in range(4):
            ref[(name, "Ar0n", n)] = rho ** n * mp.diff(fr, rho, n)
        for n in range(3):
            ref[(name, "Arn0", n)] = ti ** n * mp.diff(ft, ti, n)
        ref[(name, "Ar11", 0)] = ti * rho * mp.diff(lambda t, r: alphar_mp(p, 1 / t, r, x), (ti, rho), (1, 1))
        for i in range(len(x)):
            def fx(xi, i=i):
                xx = list(x)
                xx[i] = xi
                return alphar_mp(p, T, rho, xx)
            ref[(name, "gradx", i)] = mp.diff(fx, x[i], 1)
    return ref


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=str(HERE / "build"))
    ap.add_argument("--skip-compile-sweep", action="store_true")
    a = ap.parse_args()
    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    cver = subprocess.run([CXX, "--version"], capture_output=True, text=True).stdout.splitlines()[0]
    print(f"compiler: {cver}")
    sweep = [] if a.skip_compile_sweep else (print("compile sweep:"), compile_sweep(outdir))[1]
    print("bench ...", flush=True)
    vals, ns = run_bench(outdir)
    print("mpmath reference ...", flush=True)
    ref = reference()

    # accuracy: max relative error over states and indices, per method x quantity
    acc = {}
    for (method, state, q, idx), v in vals.items():
        r = ref[(state, q, idx)]
        if v != v:  # NaN = not provided by this method
            continue
        err = float(abs((mp.mpf(v) - r) / r))
        key = (method, q)
        acc[key] = max(acc.get(key, 0.0), err)

    Q = ["Ar0n", "Arn0", "Ar11", "gradx"]
    print("\n## Max relative error vs 60-digit mpmath (all states, all orders in the group)\n")
    print("| method | " + " | ".join(Q) + " |\n|---|" + "---|" * len(Q))
    for m in METHODS[1:]:
        print(f"| {m} | " + " | ".join(f"{acc[(m, q)]:.1e}" if (m, q) in acc else "—" for q in Q) + " |")

    print("\n## Per-order relative error, Ar0n and Arn0, state C1C2C3_dense\n")
    print("| method | A01 | A02 | A03 | A10 | A20 |\n|---|---|---|---|---|---|")
    for m in METHODS[1:]:
        cells = []
        for q, i in [("Ar0n", 1), ("Ar0n", 2), ("Ar0n", 3), ("Arn0", 1), ("Arn0", 2)]:
            v = vals[(m, "C1C2C3_dense", q, i)]
            r = ref[("C1C2C3_dense", q, i)]
            cells.append("—" if v != v else f"{float(abs((mp.mpf(v) - r) / r)):.1e}")
        print(f"| {m} | " + " | ".join(cells) + " |")

    print("\n## Time per call [ns] (min of 7 x 20000); double = one alphar evaluation\n")
    for state in STATES:
        print(f"\n{state}\n\n| method | " + " | ".join(Q) + " |\n|---|" + "---|" * len(Q))
        for m in METHODS:
            print(f"| {m} | " + " | ".join(f"{ns[(m, state, q)]:.0f}" for q in Q) + " |")

    if sweep:
        print("\n## Compile time [s] / __text size [kB] per TU\n")
        cols = [("-O2", 1), ("-O2", 4), ("-O2", 16), ("-O0", 1)]
        print("| method | " + " | ".join(f"{o} K={k}" for o, k in cols) + " | s per extra model (O2) |\n|---|" + "---|" * (len(cols) + 1))
        by = {(r["method"], r["opt"], r["K"]): r for r in sweep}
        for m in METHODS:
            cells = [f"{by[(m, o, k)]['seconds']:.2f} / {by[(m, o, k)]['text_bytes'] / 1024:.0f}" for o, k in cols]
            slope = (by[(m, "-O2", 16)]["seconds"] - by[(m, "-O2", 1)]["seconds"]) / 15
            print(f"| {m} | " + " | ".join(cells) + f" | {slope:.2f} |")

    with open(outdir / "results.json", "w") as f:
        json.dump(dict(compiler=cver, sweep=sweep,
                       values=[dict(method=k[0], state=k[1], q=k[2], i=k[3], v=v) for k, v in vals.items()],
                       ns=[dict(method=k[0], state=k[1], q=k[2], ns=v) for k, v in ns.items()],
                       reference=[dict(state=k[0], q=k[1], i=k[2], v=mp.nstr(v, 30)) for k, v in ref.items()],
                       max_rel_err=[dict(method=k[0], q=k[1], err=v) for k, v in acc.items()]), f, indent=1)


if __name__ == "__main__":
    main()
