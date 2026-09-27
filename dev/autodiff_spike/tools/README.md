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
| `split_gibbs_check.cpp` | In-model arbiter for disagreements with REFPROP: reads `idx T p` lines, flashes, and compares the published split's Gibbs energy with the best mechanically stable single-phase root (kernel all-roots). Negative Δg/RT = the split is right in-model. Arg `air` switches C1/H2S → humid air. |
| `refprop_tpflsh_check.cpp` | Reads `mixture idx T p rho_old Q_old rho_new Q_new` lines (the diff format) and adds REFPROP TPFLSH's answer, OK/x for each side. |
| `fd_lnphi_dnj.cpp` | Finite-difference check of n·∂lnφᵢ/∂nⱼ at constant T, p against the analytic forms (the #3357 Hessian building block). |
| `repro_*.cpp` | Single-state reproducers: the #3357 near-dew/near-bubble band, methanol/benzene at 308.15 K, the two Amarillo NaN-K-factor states. |
| `ssskip_audit.patch` | Audit build for #3427: whenever the SS-skip fires, also runs the old minimizer on a copy and logs its verdict (`SSSKIP_AUDIT=1`); `SSSKIP_OFF=1` reproduces pre-#3427 behaviour. Apply to master after #3427. |

The benchmark and figure scripts are one level up: `bench_ratio.cpp` (per-state CoolProp vs
TPFLSH timing CSV) and `plot_ratio_pT.py` (the (T, p) ratio maps in `figs/`).

Machine discipline: timing runs only on a quiet machine, sequentially (or with each state timing
both codes back to back when run in parallel); min of 3 timings per state.
