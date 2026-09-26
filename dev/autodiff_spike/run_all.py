#!/usr/bin/env python3
"""
Experiment 2: the whole derivative half-matrix A_ij, i + j <= N, in one call.
Compares a one-shot bivariate Taylor type (taylor2), polarization over N+1 univariate
autodiff::Real<N> passes (polar), and teqp's per-derivative scheme (teqp).

Usage: python3 run_all.py [--outdir DIR] [--skip-compile-sweep]
"""
import argparse
import json
import subprocess
import time
from pathlib import Path

import mpmath as mp

from run_spike import CXX, HERE, INC, STATES, alphar_mp, text_size

METHODS = ["taylor2", "polar", "teqp", "structured"]


def cc(src, out, *flags):
    t0 = time.perf_counter()
    subprocess.run([CXX, "-std=c++20", "-O2", *INC, *flags, "-c", str(src), "-o", str(out)], check=True)
    return time.perf_counter() - t0


def compile_sweep(outdir, repeats=3):
    res = []
    for m in METHODS:
        for n in (2, 4):
            for K in (1, 4, 16):
                o = outdir / f"sweep_{m}_N{n}_K{K}.o"
                t = min(cc(HERE / "methods_all" / f"a_{m}.cpp", o, f"-DKMODELS={K}", f"-DONLY_N={n}") for _ in range(repeats))
                res.append(dict(method=m, N=n, K=K, seconds=t, text_bytes=text_size(o)))
                print(f"  {m:8s} N={n} K={K:2d} {t:6.2f} s  text={res[-1]['text_bytes']:>9d} B", flush=True)
    return res


def reference():
    mp.mp.dps = 60
    ref = {}
    for name, (p, T, rho, x) in STATES.items():
        T, rho, x = mp.mpf(T), mp.mpf(rho), [mp.mpf(v) for v in x]
        ti = 1 / T
        f = lambda t, r: alphar_mp(p, 1 / t, r, x)
        for k in range(5):
            for j in range(k + 1):
                i = k - j
                ref[(name, i, j)] = ti ** i * rho ** j * mp.diff(f, (ti, rho), (i, j))
    return ref


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=str(HERE / "build_all"))
    ap.add_argument("--skip-compile-sweep", action="store_true")
    a = ap.parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    sweep = [] if a.skip_compile_sweep else (print("compile sweep:"), compile_sweep(out))[1]
    objs = []
    for m in METHODS:
        objs.append(out / f"a_{m}.o")
        cc(HERE / "methods_all" / f"a_{m}.cpp", objs[-1])
    for m in ("double", "numdual"):
        objs.append(out / f"m_{m}.o")
        cc(HERE / "methods" / f"m_{m}.cpp", objs[-1])
    exe = out / "bench_all"
    subprocess.run([CXX, "-std=c++20", "-O2", *INC, str(HERE / "bench_all.cpp"), *map(str, objs), "-o", str(exe)], check=True)
    vals, ns = {}, {}
    for line in subprocess.run([str(exe)], capture_output=True, text=True, check=True).stdout.splitlines()[1:]:
        kind, m, N, st, i, j, v = line.split(",")
        if kind == "val":
            vals[(m, int(N), st, int(i), int(j))] = float(v)
        else:
            ns[(m, int(N), st)] = float(v)
    print("mpmath reference ...", flush=True)
    ref = reference()

    print("\n## Max relative error by total order k = i + j (N = 4 run, all states)\n")
    print("| method | k=0 | k=1 | k=2 | k=3 | k=4 |\n|---|---|---|---|---|---|")
    for m in METHODS:
        worst = [0.0] * 5
        for (mm, N, st, i, j), v in vals.items():
            if mm == m and N == 4:
                r = ref[(st, i, j)]
                worst[i + j] = max(worst[i + j], float(abs((mp.mpf(v) - r) / r)))
        print(f"| {m} | " + " | ".join(f"{w:.1e}" for w in worst) + " |")

    print("\n## ns per call for the whole triangle i + j <= N\n")
    for st in STATES:
        print(f"\n{st}: value only {ns[('value_only', 0, st)]:.0f} ns; numdual A00..A03 {ns[('numdual_Ar0n', 3, st)]:.0f} ns;"
              f" numdual A11 {ns[('numdual_Ar11', 2, st)]:.0f} ns\n")
        print("| method | N=1 (3) | N=2 (6) | N=3 (10) | N=4 (15) |\n|---|---|---|---|---|")
        for m in METHODS:
            print(f"| {m} | " + " | ".join(f"{ns[(m, N, st)]:.0f}" for N in range(1, 5)) + " |")

    if sweep:
        print("\n## Compile time [s] / __text [kB], one TU, only order N instantiated\n")
        by = {(r["method"], r["N"], r["K"]): r for r in sweep}
        print("| method | N | K=1 | K=4 | K=16 | s / extra model | kB / extra model |\n|---|---|---|---|---|---|---|")
        for m in METHODS:
            for n in (2, 4):
                c = [by[(m, n, K)] for K in (1, 4, 16)]
                print(f"| {m} | {n} | " + " | ".join(f"{r['seconds']:.2f} / {r['text_bytes'] / 1024:.0f}" for r in c)
                      + f" | {(c[2]['seconds'] - c[0]['seconds']) / 15:.2f} | {(c[2]['text_bytes'] - c[0]['text_bytes']) / 1024 / 15:.0f} |")

    with open(out / "results_all.json", "w") as f:
        json.dump(dict(sweep=sweep, ns=[dict(method=k[0], N=k[1], state=k[2], ns=v) for k, v in ns.items()],
                       values=[dict(method=k[0], N=k[1], state=k[2], i=k[3], j=k[4], v=v) for k, v in vals.items()],
                       reference=[dict(state=k[0], i=k[1], j=k[2], v=mp.nstr(v, 30)) for k, v in ref.items()]), f, indent=1)


if __name__ == "__main__":
    main()
