# SVDSBTL Heat-Exchanger Speedup Demo — Notebook Page Design

**Issue:** CoolProp-q00v (notebook extension)
**Date:** 2026-06-01
**Status:** Approved (brainstorming complete)
**Related:** `2026-05-31-svdsbtl-hx-speedup-demo-design.md` (the C++ `dev/demo_hx_speedup.cpp` demo). This document covers the **documentation page** counterpart and does not re-derive the moving-boundary math — see that spec and the recovered reference `dev/reference/HX.py`.

## Goal

A user-facing docs page that demonstrates the **value of the SVDSBTL backend** by reproducing the Section 3.2 "computational efficiency" result of

> I.H. Bell, S. Quoilin, E. Georges, J.E. Braun, E.A. Groll, W.T. Horton, V. Lemort,
> "A generalized moving-boundary algorithm to predict the heat transfer rate of
> counterflow heat exchangers for any phase configuration,"
> *Applied Thermal Engineering* 79 (2015) 192–201.

The page runs a faithful moving-boundary heat-exchanger model on the **water-heated n-Propane evaporator** of the paper's Table 3, once with the full Helmholtz EOS (`HEOS`) and once with `SVDSBTL&HEOS`, and shows that the two agree on the predicted heat-transfer rate Q while SVDSBTL is roughly an order of magnitude faster per run — the speedup being attributable to replacing the expensive `p,h → T,ρ` transformation (Water has the costliest EOS in the library) with a direct table lookup.

## Proof of concept (already validated)

A provider-based Python port of `dev/reference/HX.py` was run before writing this spec. At A = 4 m², N = 1500 runs on the dev machine (CoolProp 7.2.1dev):

| Backend | Q [W] | ε = Q/Qmax | ms/run | one-time table build |
|---|---|---|---|---|
| `HEOS` | 4577.2422 | 0.999070 | 5.38 | — |
| `SVDSBTL&HEOS` | 4577.2431 | 0.999070 | **0.26** | ~7 s (cached after) |

- Per-run speedup ≈ **20×** (stable across N = 300 and N = 1500 runs).
- Backend agreement `|ΔQ|/Q = 2.0e-7` — well within table tolerance; the effectiveness curves overlay.
- The Python port is bit-faithful to the oracle: `HEOS` Q matches `dev/reference/HX.py` (4577.242219 W) to `8e-11`.
- `SVDSBTL`'s `smass()` works after an `HmassP` lookup, so the entropy needed for cell boundaries is available.

The ~20× *exceeds* the paper's 10× because the paper compared against TTSE on 2011 hardware; SVDSBTL is a more aggressive lookup, and this ratio is measured *with* the shared pybind dispatch tax, which compresses it. The pure-C++ ratio from `demo_hx_speedup.cpp` will be larger.

## Success criteria

1. The page builds under nbsphinx (executed via `jupyter nbconvert --execute`, 3600 s timeout) without error.
2. An in-notebook assertion confirms `HEOS` and `SVDSBTL&HEOS` agree on Q to within table tolerance for the Table-3 case, and that `HEOS` Q matches the `dev/reference/HX.py` oracle — this is the regression guard.
3. The page renders an interactive plotly effectiveness-vs-area figure with the two backends overlaid (visually coincident), and a speedup table (ms/run, ratio, `|ΔQ|/Q`).
4. The reusable solver module is downloadable from the rendered HTML.

The 20× / 10× figures are **not** targets — the page reports the measured per-run ratio on the build machine alongside the paper's reference numbers and the quoted C++ ratio for context.

## Placement & build integration

- **Notebook:** `Web/coolprop/SVDSBTLHeatExchangerDemo.ipynb`, added to the `Web/coolprop/index.rst` toctree immediately after `SVDSBTLValidation.ipynb`.
- **Module:** `Web/coolprop/hx_moving_boundary.py` — the validated provider-based solver, sitting next to the notebook so it is importable when nbconvert executes the notebook (cwd = the notebook's directory) and is also the `:download:` target.
- **Execution:** nbsphinx is already configured (`conf.py` line ~106 runs `jupyter nbconvert --allow-errors --ExecutePreprocessor.timeout=3600 --to notebook --execute`). The one-time SVDSBTL table build (~7 s/fluid, cached under `~/.CoolProp/SVDTables/`) is comfortably inside the timeout.

## Components

Each unit has one purpose, a defined interface, and is independently testable.

### (a) `hx_moving_boundary.py` — reusable module (separate, importable)

**`PropertyProvider`** — the interface the solver talks to, so the *same* solver runs against either backend. Constructed per fluid. Mirrors the C++ spec's interface:

- `h_pT(p, T)` → hmass (via `PT_INPUTS`)
- `Ts_ph(p, h)` → (T, smass) from one `HmassP_INPUTS` update (the hot-loop transform; one update, two reads)
- `Tsat(p, Q)`, `hsat_TQ(T, Q)` → saturation (init only, not in the hot loop)

Implementations:

- **`HEOSProvider`** — one `HEOS` `AbstractState` for all calls.
- **`SVDSBTLProvider`** — a `SVDSBTL&HEOS` `AbstractState` for `h_pT` / `Ts_ph`; saturation routed through a **separate fast HEOS ancillary** `AbstractState` (SVDSBTL exposes no public PQ / T_sat path). **Both** providers source saturation through the identical HEOS ancillary route, so the only timed difference is the `p,h ↔ T,ρ` table lookup — the isolation approved in the C++ spec's "saturation sourcing" judgment.

**`HeatExchanger`** — a faithful Python port of `dev/reference/HX.py`'s `HeatExchanger`, with every property call routed through a `PropertyProvider` instead of hardcoded `PropsSI`:

- `external_pinching()` → Q_max,ext (eqns 4–7)
- `calculate_cell_boundaries(Q)` → sorted enthalpy vectors with phase-transition enthalpies inserted, plus per-boundary T and phase labels
- `internal_pinching(stream)` → reduced Q_max (eqns 9–16)
- `objective_function(Q)` = 1 − Σ wⱼ (eqn 28)
- `solve()` → `scipy.optimize.brentq(objective_function, 1e-5, Qmax-1e-10, rtol=1e-14, xtol=1e-10)`

Numerics are unchanged from the oracle (validated to `8e-11`). The module also provides a `make_evaporator(backend, A)` helper that constructs the Table-3 case.

### (b) `SVDSBTLHeatExchangerDemo.ipynb` — narrative notebook (self-contained narrative, imports the module)

Cells, in order:

1. **Intro (markdown)** — the 2015 paper, the §3.2 result, one-line "what SVDSBTL is," link to `SVDSBTL.rst`.
2. **Download link (raw reST cell)** — `:download:` role pointing at `hx_moving_boundary.py`, so the script is grabbable from the HTML.
3. **Imports + plotly setup (code)** — `from hx_moving_boundary import ...`; `pio.renderers.default = 'notebook_connected'` (matches `SVDSBTLValidation.ipynb`, renders interactively in HTML).
4. **Build (code)** — construct `HEOS` and `SVDSBTL&HEOS` providers for Water (hot) and n-Propane (cold); first SVDSBTL construction builds/caches tables, timed and reported *separately* (excluded from per-run timing, as in the paper).
5. **Correctness (code)** — solve once per backend at A = 4 m²; assert `|ΔQ|/Q` within table tolerance and assert `HEOS` Q matches the `dev/reference/HX.py` oracle (4577.242219 W); print Q, Qmax, ε.
6. **Effectiveness figure (code)** — ε = Q/Qmax vs area, log sweep 0.1–10 m², `HEOS` vs `SVDSBTL&HEOS` overlaid → visually coincident plotly traces.
7. **Timing (code)** — per-run timing at A = 4 m² averaged over `HX_REPEATS` (env-tunable; default 200, keeping doc-build timing under ~2 s given ~5 ms/run for HEOS) for both backends; report ms/run, ratio, `|ΔQ|/Q` in a small table.
8. **Closing (markdown)** — quote `demo_hx_speedup.cpp`'s pure-C++ ratio as the "no Python tax" library-level figure; caveat that absolute ms differ from the paper's 2011 hardware; the **Table-3 deviation note** (below).

## Fidelity & consistency

The notebook, `dev/reference/HX.py`, and `dev/demo_hx_speedup.cpp` all use the **same oracle case** so one set of numbers validates all three:

- Water (hot) ṁ = 0.1 kg/s, n-Propane (cold) ṁ = 0.01 kg/s — the **transposition** present in the original `PropaneEvaporatorPinching()` (Water given `mdot_c`, propane given `mdot_h`), opposite of the paper's printed Table 3 and the source of the evaporating cell structure.
- Two-phase α = 1000 W m⁻² K⁻¹, single-phase α = 100 W m⁻² K⁻¹ (as hard-coded in the original `objective_function`).
- Inlet states exactly as `PropaneEvaporatorPinching()`: p_Water = 101325 Pa, T_h,i = 330 K; n-Propane p_ref at T = 300 K / Q = 1, h_ref at T = 275 K.

The closing markdown explains these two deviations from the published Table 3 and why the page keeps them (they are what the recovered script — which produced the published evaporating-case figures — actually used; keeping them makes the page consistent with the C++ demo and the reference oracle).

## Plotting

Interactive **plotly** (`plotly.graph_objects` + `plotly.io`), `pio.renderers.default = 'notebook_connected'`, matching the nearest sibling `SVDSBTLValidation.ipynb`. The effectiveness figure overlays the two backends; hover shows (area, ε) per trace.

## Testing / validation

- **Primary guard:** the notebook's own correctness cell (criterion 2) fails the doc build if the backends disagree or the port drifts from the oracle.
- The C++ correctness gate (Catch2 `[SBTL]` test) is owned by the C++ demo spec, not duplicated here.
- The module is import-clean and side-effect-free at import (only class/function defs + `make_evaporator`), so it can be exercised directly from a Python REPL or a future pytest without running the notebook.

## Out of scope (YAGNI)

- Figs 8 (residual r(Q)) and 9 (Brent convergence) — the full-paper-reproduction option was declined; this page is the focused speedup story (effectiveness figure only).
- No new public SVDSBTL saturation API (saturation goes through the HEOS ancillary, as in the C++ spec).
- No mixtures; hot and cold are separate pure fluids.
- No corrected/as-published Table-3 variant alongside the oracle case (the "show both" option was declined).
- No matplotlib path; plotly only.
