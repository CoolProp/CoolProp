# Reference-data audit: incompressible fluids (2026-07-05)

Full physical-plausibility audit of every number feeding the incompressible
fits: 36 top-level `CPIncomp/data/*.txt` grids, ~180 SecCool
`xMass`/`xVolume`/`xPure` tables, 20 `xTables/xMass/*.csv`, the hardcoded
arrays in `PureFluids.py`/`SolutionFluids.py`/`ExampleObjects.py`, and the
committed fit residuals (NRMS) in all `json/*.json`.

**Headline: no wrong-unit slips exist in any production data source.** The
worst committed fit residual on a real fluid is 4.8% (TX22 saturation
pressure) — a rogue data point would blow that to tens of percent, and only
the synthetic `ExampleDigitalPure` (13.8%) does. Every value that trips a
generic liquid range check is either correct physics or a valid alternate
convention the loader already handles (see the don't-touch list).

## Real findings (fixed / to respect)

| # | Finding | Status |
|---|---|---|
| 1 | `data/SecCool/xPure/HFE-7100_{Rho,Cp,Mu,Cond}.txt` store T strictly **descending** (+64.27 → −80.33 °C) — the only non-ascending grids in the corpus. | The loader reversed fully-descending grids already; now **hardened to a real sort** (rows by T, columns by concentration) in both `SecCoolSolutionData.getFromFile` and `DigitalData.getFromFile`, and guarded by `test_data_sanity.py`. Raw files deliberately untouched (provenance). |
| 2 | Short grids that cannot support a temperature fit above order ~4: `IceNA` csv (4 T points), several `*_TFreeze` tables (4–5 points), `FRE_Tfreeze` (freeze curve, 1 T row). | Any fitter must cap the T-degree at `N_T − 1` and skip properties with < 3 usable points (the Chebyshev fitter does). |
| 3 | `data/SecCool/xMass/VDI, Methanol_*` grids are ~66% `-1`-sentinel — sparse but valid (data lives in a narrow T×X band). | Nothing to fix; NaN masking handles it. |
| 4 | Empty `*_Vol2Mass.txt` grids for mass-based fluids. | Expected — mass-based fluids carry no volume→mass table. |
| 5 | Nine **orphaned** `xTables/xMass/*.csv` files (`Freezium_{Cond,Cp,Mu}`, `Ice{EA,NA,PG}_{Cond,Mu}`) are **latin-1 encoded** (`·` in unit headers) and fail a default UTF-8 read. | Not read by the production pipeline (only the `Hfusion` csvs are loaded; Freezium reads the root `FRE_*.txt`). Anyone reviving them — e.g. to restore the ice-slurry conductivity/viscosity data (see the reproducibility follow-up) — must pass `encoding="latin-1"`. |

## Don't-touch list (looks wrong, is right)

| Value | Where | Why it's correct |
|---|---|---|
| density → 239 kg/m³, conductivity up to **87 W/m/K**, Prandtl 0.004–0.009 | `PureFluids.py` → `LiquidSodium` (LiqNa) | Textbook molten-sodium physics up to 2500 K. A "sanity fix" here would be the bug. |
| density 0.67–1.79 kg/m³, conductivity 0.018–0.041 W/m/K | `data/Air_*.txt` | Air is included as a **gas** reference fluid. |
| viscosity up to **59.3 Pa·s** | `xVolume/Zitrec LC_Mu.txt` at −50 °C / 70 vol-% | Concentrated glycol near its glass region; decays smoothly to ~0.5 mPa·s at 100 °C. Raw file already in Pa·s (`viscosityFactor=None` is correct). |
| freeze temperatures 214–263 (no °C offset) | `FRE_Tfreeze.txt`, `xTables/.../Freezium_TFreeze.csv` | These two tables are in **Kelvin** (csv header says so); the SecCool `_TFreeze.txt` files are in °C. Both conventions handled by their loaders. |
| freeze point −125.4 °C | `xVolume/Zitrec M_TFreeze.txt` at 100 vol-% | Monotonic freeze-depression curve endpoint (pure-glycol glass former), internally consistent. |
| cp dips to 437 J/kg/K; density rising with T | `ExampleObjects.py` → `DigitalExample` | Synthetic analytic test functions, not fluid data. |

## Unit conventions per source (for anyone touching loaders)

- Top-level `data/*.txt`: headerless; row 0 = concentration fraction
  (NaN corner), column 0 = T in **K**; property values in SI.
- `SecCool/x*/**.txt`: tab-separated `T\X` header; T in **°C**, concentration
  in **%**; property units vary per file and are normalized by the
  `densityFactor`/`heatFactor`/`conductivityFactor`/`viscosityFactor`
  arguments in `SecCoolFluids.py::factory()` — the factors were verified to
  match each file's raw units. `-1` is the out-of-range sentinel.
- `SecCool/xTables/**.csv`: 3 header rows including units; SI.
- Hardcoded arrays: SI, `np.nan` for missing points.
