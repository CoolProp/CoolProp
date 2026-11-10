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
| 5 | **23 of the 32** `xTables/*.csv` files were **latin-1 encoded** (`·` in unit headers) and failed a default UTF-8 read — the nine orphaned ones (`Freezium_{Cond,Cp,Mu}`, `Ice{EA,NA,PG}_{Cond,Mu}`) plus `AS{10,20,30,40,55}`, `HY{20,30,40,45,50}` and four `ASHRAE`/`CO2_*` tables. | The six `Ice{EA,NA,PG}_{Cond,Mu}` files are **now plain ASCII** (`W/(m K)`, `Pa s`) and back in the production pipeline: `SecCoolSolutionData.__init__` read them inside a bare `try/except`, so the `UnicodeDecodeError` silently dropped ice-slurry conductivity and viscosity from every refit (issue #3303). The refit now reproduces the committed coefficients to ~3e-7 relative. `test_data_sanity.py::test_loaded_data_files_are_ascii` pins the encoding for every file the loaders read, and `SecCoolIceData.getRequiredArray` refuses to construct the fluid when a grid is missing or unreadable (moving the reads out of the bare `try/except` alone was not enough: `getArray` returns `(None, None, None)` for a file it cannot find, rather than raising). The remaining 17 are still latin-1 and still orphaned; anyone reviving those must re-encode them the same way. |

## Don't-touch list (looks wrong, is right)

| Value | Where | Why it's correct |
|---|---|---|
| density → 239 kg/m³, conductivity up to **87 W/m/K**, Prandtl 0.004–0.009 | `PureFluids.py` → `LiquidSodium` (LiqNa) | Textbook molten-sodium physics up to 2500 K. A "sanity fix" here would be the bug. |
| density 0.67–1.79 kg/m³, conductivity 0.018–0.041 W/m/K | `data/Air_*.txt` | Air is included as a **gas** reference fluid. |
| viscosity up to **59.3 Pa·s** | `xVolume/Zitrec LC_Mu.txt` at −50 °C / 70 vol-% | Concentrated glycol near its glass region; decays smoothly to ~0.5 mPa·s at 100 °C. Raw file already in Pa·s (`viscosityFactor=None` is correct). |
| freeze temperatures 214–263 (no °C offset) | `FRE_Tfreeze.txt`, `xTables/.../Freezium_TFreeze.csv` | These two tables are in **Kelvin** (csv header says so); the SecCool `_TFreeze.txt` files are in °C. Both conventions handled by their loaders. |
| freeze point −125.4 °C | `xVolume/Zitrec M_TFreeze.txt` at 100 vol-% | Monotonic freeze-depression curve endpoint (pure-glycol glass former), internally consistent. |
| cp dips to 437 J/kg/K; density rising with T | `ExampleObjects.py` → `DigitalExample` | Synthetic analytic test functions, not fluid data. |
| no `T_freeze` at all | `Ice{EA,NA,PG}` | There is no `Ice*_TFreeze` table in the corpus, and none is missing. These are **ice slurries**: the composition axis is the *ice mass fraction*, so the fluid is already at solid-liquid equilibrium everywhere in its range. The equilibrium temperature is the `T` axis of the tables, not a separate curve over `x`. `T_freeze` stays `notdefined` and the backend raises, which is correct (issue #2567). |

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
