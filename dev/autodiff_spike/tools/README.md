# Mixture-flash validation tools

Small drivers used to validate changes to the mixture PT / PQ / QT flash and the Chebyshev
density kernel (Linear: initiative *Mixture flash algorithms*; method write-up in the Linear doc
"How to validate a mixture-flash change"). They link statically against a Release CoolProp.

## Build

```bash
B=<release build dir>   # e.g. cmake -B $B -S . -DCMAKE_BUILD_TYPE=Release && cmake --build $B --target CoolProp
FL=(${=$(grep CXX_INCLUDES $(find $B/CMakeFiles/CoolProp.dir -name flags.make | head -1) | sed 's/^CXX_INCLUDES = //')})
clang++ -std=c++17 -O2 $FL -Isrc -Idev -Idev/autodiff_spike dev/autodiff_spike/tools/<tool>.cpp \
        -Wl,-force_load,$B/libCoolProp.a -o <tool>
```

(`${=...}` is zsh word splitting; in bash use an unquoted `$(...)`.) `refprop_tpflsh_check`
needs no CoolProp link; it dlopen()s REFPROP from `/Users/ianbell/REFPROP10/`.

Environment switches read by the library on this branch: `COOLPROP_CHEB_DENSITY=1` (kernel),
`SPIKE_SSSKIP=1`, `SPIKE_KDIRECT=1`, `SPIKE_ITER=1`, `SPIKE_DUMP=1` (spike instrumentation in
VLERoutines.cpp). Each is read once per process.

## Tools

| tool | what it does |
|---|---|
| `flash_dump.cpp` | 2000 (or `NS=`) states × 5 mixtures through `update(PT_INPUTS)`; one line per state (rho, Q) on stdout, median time per mixture on stderr. Run on two builds and diff per state: the equivalence check. |
| `tpd_counters.cpp` | Same states; median/mean time plus per-flash work counters (trial density solves, warm/global/kernel, SS steps, TPD iterations and exit reasons, split iterations, full αʳ derivative evaluations). Needs the `spike_counts` instrumentation on this branch. |
| `state_scan.cpp` | Same states, prints a `STATE` tag to stderr before each flash, so library-side tracers can be attributed to a state. |
| `pqqt_bench.cpp` | PQ and QT over 12 × 11 grids for 4 mixtures: time, failures, fallbacks (`warnstring`), mass-balance and equal-fugacity error of the published state; `PQ_DETAIL=1` lists offending states. |
| `split_gibbs_check.cpp` | In-model arbiter for disagreements with REFPROP: reads `idx T p` lines, flashes, and compares the published split's Gibbs energy with the kernel's **spinodal-branch selected** single-phase root. Negative Δg/RT = the split is right in-model. Arg `air` switches C1/H2S → humid air. Never use min-g over all roots as the reference: at low T the αʳ-well roots near δ ≈ 1 have spuriously low g. |
| `tpflsh_split_validity.cpp` | For states where TPFLSH says two-phase and CoolProp does not: evaluates TPFLSH's own split (x, y, ρL, ρV, q) in CoolProp's GERG-2008 — equal fugacities? both phases at the requested p? Δg vs CoolProp's single phase. Separates genuine CoolProp misses from TPFLSH returning non-equilibrium "splits" with ierr = 0. |
| `published_split_validity.cpp` | Runs a state sequence through ONE object and checks every published split through the public API only (equal fugacity, mass balance). `SET_Z_EACH=1` calls set_mole_fractions before every update. Use this, not the internal SatL/SatV objects, to judge what users get. |
| `refprop_tpflsh_check.cpp` | Reads `mixture idx T p rho_old Q_old rho_new Q_new` lines (the diff format) and adds REFPROP TPFLSH's answer, OK/x for each side. |
| `fd_lnphi_dnj.cpp` | Finite-difference check of n·∂lnφᵢ/∂nⱼ at constant T, p against the analytic forms (the #3357 Hessian building block). |
| `repro_*.cpp` | Single-state reproducers: the #3357 near-dew/near-bubble band, methanol/benzene at 308.15 K, the two Amarillo NaN-K-factor states. |
| `ssskip_audit.patch` | Audit build for #3427: whenever the SS-skip fires, also runs the old minimizer on a copy and logs its verdict (`SSSKIP_AUDIT=1`); `SSSKIP_OFF=1` reproduces pre-#3427 behaviour. Apply to master after #3427. |
| `verdict.cpp` | v1 per-state verdicts for CoolProp/TPFLSH **disagreements** in a `bench_ratio` CSV. Judges with CoolProp's own GERG code and the kernel's `select()` root, and never looks at states where the two codes agree. Superseded for the figures by the v2 pipeline below. |
| `rp_stability.cpp` | Brute-force tangent-plane stability test using **only REFPROP** (GERG mode: `TPRHO` + `FGCTY2` + `PRESS`). Five near-pure plus 20 random trial compositions, both density roots, SS to 1e-12. Args `'A.FLD\|B.FLD' 'z1,z2'`; reads `idx T p[Pa]`, prints `idx T p tm_min` (tm < 0: a split exists). Shares no code with CoolProp. |
| `rp_selfjudge.cpp` | Re-evaluates TPFLSH's own two-phase answers with REFPROP's own `FGCTY2`/`PRESS`: the max of \|Δln f\| and \|p_phase/p − 1\| per state. Reads `rpwrong_states2.tsv` (`mixture T p idx`). |
| `rp_truth_n.cpp` | REFPROP-only per-state reference for any N (`'A.FLD\|B.FLD' 'z1,z2' gerg\|default`). |
| `score_all.py` | Scores CoolProp builds and TPFLSH against the reference at every state (v3). |
| `build_verdicts_v2.py` | Builds the v2 verdicts (`v2/*.csv`) from the v1 files plus the two REFPROP-only judges. |

### Verdicts v3 (2026-10-01): score every state against a REFPROP-only reference

`rp_truth_n.cpp` gives every state its own answer from REFPROP routines alone: the lower-g single-phase root, a
brute-force stability test, and an independently converged split.  `score_all.py` then scores CoolProp (each
build) and TPFLSH against it at **every** state.  v1 and v2 judged only disagreements, or only states both codes
called single phase, and that hid CoolProp misses twice.  The figures are made with:
`score_all.py v3 on:<dir> master:<dir> -- <truth dir>`.  `tools/v3/` holds the verdicts.

### Verdicts v2 (2026-09-29): what changed and why

v1 flattered CoolProp. It only judged states where the two codes disagreed, and REFPROP never searches for liquid-liquid splits, so CoolProp LLE misses that REFPROP shared were invisible. It also used CoolProp code as the judge. v2 fixes both:

- **CoolProp wrong: 4 → 28.** `rp_stability` on all ~72k states that both codes call single phase finds 24 missed splits, all in the N2/C1/C2/nC4/nC5 LLE region (104–119 K, 1–28 MPa). It finds none in the other 8 GERG mixtures. The 4 v1 cases are confirmed.
- **TPFLSH wrong: 1024 → 767**, re-judged by REFPROP's own routines. The other 257 are `rp_minor`:
  - 219 where TPFLSH's density matches CoolProp's single phase and only the two-phase label is wrong;
  - loose splits off equilibrium by only 1e-5..1e-3;
  - 22 "not reproducible" re-calls. Those were made from 10-digit-rounded CSV inputs, the same artifact as `cp_history`.
- **TPFLSH misses a verified split: 1719.** All are confirmed by `rp_stability`.
- The R454B disagreements (2816) are not judged: the two sides use different models.

Validation of the judge: it flags all 328 splits CoolProp publishes in the N2-mixture LLE box. It needs a guard against 0·ln 0 when a trace component's trial fraction underflows; without it, humid air gave 14 false positives.

The benchmark and figure scripts are one level up: `bench_ratio.cpp` (per-state CoolProp vs
TPFLSH timing CSV) and `plot_ratio_pT.py` (the (T, p) ratio maps in `figs/`; `build:"..."` is required and names the CoolProp build in the title).  `figs/ratio_pT_vs_tpflsh_all_verdicts.png` is the UNMERGED kernel build with v2 verdicts; `figs/ratio_pT_vs_tpflsh_all_kernel_off.png` is the same states with the kernel off, which is close to what master ships.  Each panel title gives both the median per-state ratio and the total-time ratio. They differ a lot where either code has a heavy tail, e.g. Amarillo: median state 3.9× faster, total time 1.1× faster.

Machine discipline: timing runs only on a quiet machine, sequentially (or with each state timing
both codes back to back when run in parallel); min of 3 timings per state.
