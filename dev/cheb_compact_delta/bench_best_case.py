"""Stress-test benchmarks: paper's 'no τ flattening' best case.

Simulates a pressure sweep at fixed T (the common flash inner-loop pattern
where T is constant across many calls). set_tau() is called once, then
density_roots_fixed_tau() many times. This matches paper Fig. 11's
"No τ flattening" curve at ~30 μs.
"""

import time
import numpy as np

import cpp_solver_v2
from fluid_from_json import FluidFromJSON
from fast_solver import build_fluid_fast


def build_cpp_fluid_v2(fluid, state):
    terms = state["terms"]
    n_terms = state["n_terms"]
    family = np.zeros(n_terms, dtype=np.int32)
    t_vals = np.zeros(n_terms)
    beta_vals = np.zeros(n_terms)
    gamma_vals = np.zeros(n_terms)
    for k, term in enumerate(terms):
        if term["family"] == "power":
            family[k] = 0
            t_vals[k] = term["F"].__defaults__[0]
        else:
            family[k] = 1
            defaults = term["F"].__defaults__
            t_vals[k] = defaults[0]
            beta_vals[k] = defaults[1]
            gamma_vals[k] = defaults[2]
    return cpp_solver_v2.FluidChebV2(
        state["K"], n_terms, state["n_delta"],
        fluid.rhoc, fluid.R, fluid.Tc,
        state["intervals_a"].astype(np.float64),
        state["intervals_b"].astype(np.float64),
        family,
        t_vals, beta_vals, gamma_vals,
        state["n_values"].astype(np.float64),
        state["H_flat"].ravel().astype(np.float64),
    )


def bench_fixed_tau(cpp, T, P_vals, N_rep=20):
    """Simulate a pressure sweep at fixed T: set_tau once, solve at each P.
    set_tau is called ONCE outside the timed loop; only density_roots_fixed_tau
    is measured (true per-call cost at fixed τ)."""
    cpp.set_tau(T)
    # Warm-up
    for P in P_vals:
        cpp.density_roots_fixed_tau(T, P)
    t0 = time.perf_counter()
    for _ in range(N_rep):
        for P in P_vals:
            cpp.density_roots_fixed_tau(T, P)
    total = time.perf_counter() - t0
    return total / (N_rep * len(P_vals))


def bench_single_point_fixed_tau(cpp, T, P, N=10000):
    cpp.set_tau(T)
    cpp.density_roots_fixed_tau(T, P)
    t0 = time.perf_counter()
    for _ in range(N):
        cpp.density_roots_fixed_tau(T, P)
    return (time.perf_counter() - t0) / N


def bench_full(cpp, T, P, N=5000):
    cpp.density_roots(T, P)
    t0 = time.perf_counter()
    for _ in range(N):
        cpp.density_roots(T, P)
    return (time.perf_counter() - t0) / N


def main():
    path = "/Users/ianbell/miniforge3/lib/python3.9/site-packages/teqp/fluiddata/dev/fluids/Nitrogen.json"
    fluid = FluidFromJSON(path)
    print(f"Fluid: {fluid.name}  Tc={fluid.Tc}  rhoc={fluid.rhoc:.3f}")
    print(f"\nPaper (Fig. 11) reference timings for small-component mixtures:")
    print(f"  'No τ flattening' (fixed τ, pressure sweep) : ~30 μs")
    print(f"  'Full flattening' (τ and composition both)  : ~50–150 μs")
    print(f"  REFPROP TPRHOdll                            : ~100–300 μs")
    print()

    for n_delta, tol in [(8, 1e-8), (12, 1e-8), (12, 1e-10)]:
        state = build_fluid_fast(fluid, 1e-10, 5.0, n_delta=n_delta, tol=tol)
        K = state["K"]
        cpp = build_cpp_fluid_v2(fluid, state)
        print(f"=== n_delta={n_delta}, tol={tol:.0e}, K={K} intervals ===")


        # Specific representative points (single root / multi root)
        points = [
            (200.0, 5e6, "supercrit 1 root"),
            (200.0, 5e7, "supercrit dense"),
            (100.0, 5e4, "subcrit vapor"),
            (100.0, 5e6, "subcrit liquid"),
            (80.0, 1e4, "subcrit vapor deep"),
            (80.0, 5e6, "subcrit liquid deep"),
            (126.2, 3.4e6, "near-critical"),
        ]
        print(f"  {'Point':<25} {'full':>10} {'fixed-τ':>10}")
        for T, P, label in points:
            t_full = bench_full(cpp, T, P, N=5000)
            t_ft = bench_single_point_fixed_tau(cpp, T, P, N=20000)
            print(f"  {label:<25} T={T:5.1f}K P={P:7.1e}  "
                  f"full={t_full*1e6:6.1f} μs  fixed-τ={t_ft*1e6:6.1f} μs")
        print()


if __name__ == "__main__":
    main()
