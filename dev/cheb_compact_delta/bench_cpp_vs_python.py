"""Compare C++ density rootfinder against Python reference on nitrogen.

Benchmarks over (T, P) grids representative of typical flash inner-loop calls,
measures per-call wall time, verifies answers match the Python reference to
machine precision.

Reference for the single-digit μs target: Bell & Alpert 2018 Figure 11.
    "No τ flattening" (fixed τ, P sweep)  ~30 μs for 2–10 components
    "Full flattening"                      ~50–150 μs for 2–10 components
    REFPROP TPRHOdll                       ~100–300 μs
"""

import time
import numpy as np

import cpp_solver
import cpp_solver_v2
from fluid_from_json import FluidFromJSON
from fast_solver import build_fluid_fast, density_roots_vectorized


def build_cpp_fluid_v2(fluid, state):
    """Map Python state to optimized C++ FluidChebV2."""
    terms = state["terms"]
    n_terms = state["n_terms"]
    K = state["K"]
    n_delta = state["n_delta"]
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
        K, n_terms, n_delta,
        fluid.rhoc, fluid.R, fluid.Tc,
        state["intervals_a"].astype(np.float64),
        state["intervals_b"].astype(np.float64),
        family,
        t_vals, beta_vals, gamma_vals,
        state["n_values"].astype(np.float64),
        state["H_flat"].ravel().astype(np.float64),
    )


def build_cpp_fluid(fluid, state):
    """Map Python state to C++ FluidChebCpp."""
    terms = state["terms"]
    n_terms = state["n_terms"]
    K = state["K"]
    n_delta = state["n_delta"]

    # Term metadata for F_k(τ) eval
    family = np.zeros(n_terms, dtype=np.int32)
    t_vals = np.zeros(n_terms)
    beta_vals = np.zeros(n_terms)
    gamma_vals = np.zeros(n_terms)
    for k, term in enumerate(terms):
        if term["family"] == "power":
            family[k] = 0  # POWER
            # F_k(τ) = τ^t ; the index-into-power-coefficients is term["index"]
            t_vals[k] = term["F"].__defaults__[0]  # captured _t
            beta_vals[k] = 0.0
            gamma_vals[k] = 0.0
        else:  # gauss
            family[k] = 1  # GAUSS
            # F_k(τ) = τ^t exp(-β(τ-γ)^2); defaults order: _t, _b, _g
            defaults = term["F"].__defaults__
            t_vals[k] = defaults[0]
            beta_vals[k] = defaults[1]
            gamma_vals[k] = defaults[2]

    return cpp_solver.FluidChebCpp(
        K, n_terms, n_delta,
        fluid.rhoc, fluid.R, fluid.Tc,
        state["intervals_a"].astype(np.float64),
        state["intervals_b"].astype(np.float64),
        family,
        t_vals, beta_vals, gamma_vals,
        state["n_values"].astype(np.float64),
        state["H_flat"].ravel().astype(np.float64),
    )


def time_case(solver, T, P, n_rep=5000):
    # Warm-up
    solver.density_roots(T, P)
    t0 = time.perf_counter()
    for _ in range(n_rep):
        solver.density_roots(T, P)
    return (time.perf_counter() - t0) / n_rep


def time_python_case(state, T, P, n_rep=500):
    density_roots_vectorized(state, T, P)
    t0 = time.perf_counter()
    for _ in range(n_rep):
        density_roots_vectorized(state, T, P)
    return (time.perf_counter() - t0) / n_rep


def main():
    path = "/Users/ianbell/miniforge3/lib/python3.9/site-packages/teqp/fluiddata/dev/fluids/Nitrogen.json"
    fluid = FluidFromJSON(path)

    print(f"Fluid: {fluid.name}  Tc={fluid.Tc} K  rhoc={fluid.rhoc:.3f}  M={fluid.M}")
    for n_delta in (8, 12):
        for tol in (1e-8, 1e-10):
            state = build_fluid_fast(fluid, 1e-10, 5.0, n_delta=n_delta, tol=tol)
            K = state["K"]
            cpp = build_cpp_fluid(fluid, state)
            cpp_v2 = build_cpp_fluid_v2(fluid, state)
            print(f"\n{'='*80}\nn_delta={n_delta}, tol={tol:.0e}: K={K} intervals")

            cases = [
                (100.0, 1e5,  "T=100K P=0.1MPa     (supercrit vapor-ish)"),
                (126.2, 3.4e6, "T=Tc P=Pc            (critical point)"),
                (200.0, 5e7,  "T=200K P=50MPa       (supercrit dense)"),
                (80.0, 1e5,   "T=80K  P=0.1MPa      (subcrit vapor)"),
                (80.0, 2e6,   "T=80K  P=2MPa        (subcrit liquid)"),
                (75.0, 7.6e4, "T=75K  P=Psat        (both phases present)"),
            ]
            print(f"{'Case':<40} {'Py':>8} {'C++ v1':>8} {'C++ v2':>8} "
                  f"{'v1/v2':>6} {'match v2':>10}")
            for T, P, label in cases:
                py_roots, _ = density_roots_vectorized(state, T, P)
                cpp_roots = np.sort(np.array(cpp.density_roots(T, P)))
                cpp_v2_roots = np.sort(np.array(cpp_v2.density_roots(T, P)))
                if len(py_roots) != len(cpp_v2_roots):
                    match = f"DIFF {len(py_roots)}→{len(cpp_v2_roots)}"
                else:
                    rel = np.max(np.abs(py_roots - cpp_v2_roots) / np.maximum(np.abs(py_roots), 1e-10))
                    match = f"{rel:.1e}"
                t_py = time_python_case(state, T, P, 200)
                t_cpp = time_case(cpp, T, P, 3000)
                t_v2 = time_case(cpp_v2, T, P, 10000)
                print(f"  {label:<40} {t_py*1e6:6.0f}us {t_cpp*1e6:6.1f}us "
                      f"{t_v2*1e6:6.2f}us {(t_cpp/t_v2):>5.1f}x {match:>10}")


if __name__ == "__main__":
    main()
