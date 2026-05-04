# Mixture density rootfinder — baseline characterization

**Issue:** [CoolProp-aor]
**Code:** `dev/bench_density_baseline.cpp`, `dev/bench_density_baseline.py`
**Data:** `dev/bench_density_<system>.csv` (60×60 grid in T,p)
**Plots:** `dev/bench_density_<system>.png`

## Summary

| System | Cells | Flash err | Solver-isolated fail | Flash p50/p95/max (ms) | Solver p50/p95/max (μs) |
|---|---:|---:|---:|---|---|
| amarillo (10-comp AGA-8) | 3600 | 14 | **1399 (39%)** | 65.6 / 154.6 / 328.4 | 190.5 / 1376.5 / 5160.9 |
| CH4/H2S 70/30           | 3600 | 182 | **2214 (65%)** | 5.6 / 9.0 / 37.9 | 33.4 / 123.8 / 369.4 |
| CO2/N2 80/20            | 3600 | 15 | **1740 (49%)** | 13.3 / 43.5 / 413.6 | 74.4 / 193.6 / 938.9 |
| n-Decane/Methane 30/70  | 3600 | 48 | **895 (25%)**  | 3.9 / 8.2 / 18.3 | 15.7 / 54.7 / 265.6 |

"Solver-isolated fail" = `solver_rho_Tp(T, p, rho_guess=-1)` throws or returns non-finite at a cell where the full PT flash converged.

## Key findings

1. **The bare rootfinder is far more fragile than the flash-orchestrated path.** Across all four systems, 25–65 % of cells where the full flash succeeds cannot be re-solved by `solver_rho_Tp` called in isolation. Failures cluster along phase-envelope boundaries (visible in PNGs as red/white bands inside the bubble/dew curves).

2. **`solver_rho_Tp` cost scales steeply with component count.** Amarillo (10 components) has p95 = 1.4 ms and max = 5.2 ms per density rootfind. The binary mixtures sit at p95 ≈ 50–200 μs. Each Householder/Brent iteration evaluates `alphar` and its derivatives, all of which scale O(N²) with reducing-function combining rules.

3. **The full flash is dominated by orchestration + Newton outer loop, not the density rootfind itself.** Flash p95 / solver p95 ratios:
   - amarillo: 155 ms / 1.4 ms = ~110×
   - CH4/H2S: 9 ms / 0.12 ms = ~75×
   - CO2/N2: 44 ms / 0.19 ms = ~230×
   - c10c1: 8 ms / 0.05 ms = ~150×
   These ratios bound how much a faster rootfinder alone can speed up the user-visible cost — typically ≤ 2× without restructuring the flash.

4. **Near-critical and two-phase regions dominate flash tail latency.** CO2/N2 max flash time of 414 ms is in the immediate vicinity of the CO2 critical (304 K, 74 bar). Amarillo p95 of 155 ms tracks the high-pressure side of the dew curve.

## Implications for Gernert / Chebyshev assessment

- **Standalone robustness is the primary acceptance bar.** A replacement that matches `solver_rho_Tp_global` (bracketed Brent) on isolation-failure rate is more valuable than one that's faster but equally fragile.
- **Solver speed-up alone caps overall flash speedup.** To break the ~2× ceiling we'd need to also reduce the number of solver invocations per flash (separate scope).
- **N-scaling matters at large component counts.** Amarillo's 5 ms max should be a hard target; a Chebyshev pre-tabulation in (τ, δ) avoids the per-evaluation derivative cost but at memory cost — quantify in `CoolProp-3gt`.
- **Gas-side and near-critical regions are the priority hardpoints** for the robustness test set (`CoolProp-568`).

## Caveats

- "Solver-isolated fail" is a diagnostic, not a real-world failure mode — production callers reach `solver_rho_Tp` only via `update_TP`, which provides a phase index and an SRK-derived guess. The high failure count quantifies how much the flash orchestration is doing for the rootfinder, not user-visible breakage.
- Per-cell `flash_ms` and `solver_us` are best-of-N timings on a warm backend; cold-start (factory + first call) costs are not measured here.
- Tier breakdown (Wilson-SS / GDEM / HELD) is not captured because `get_last_pt_flash_tier()` is not implemented (filed as `CoolProp-29n`).

## Reproducing

```bash
mkdir -p build_bench
cd build_bench
cmake .. -G Ninja \
  -DCOOLPROP_MY_MAIN=../dev/bench_density_baseline.cpp \
  -DCOOLPROP_STATIC_LIBRARY=ON \
  -DCMAKE_BUILD_TYPE=Release
ninja Main

cd ..
for sys in amarillo ch4h2s co2n2 c10c1; do
  build_bench/Main "$sys" > "dev/bench_density_${sys}.csv" \
                          2> "dev/bench_density_${sys}.log"
done

python3 dev/bench_density_baseline.py
```

[CoolProp-aor]: ../.beads/  "bd show CoolProp-aor"
