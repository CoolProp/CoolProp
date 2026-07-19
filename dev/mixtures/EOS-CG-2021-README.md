# EOS-CG-2021 binary mixture models — implementation notes

Adds/updates binary reducing parameters and departure functions from

> T. Neumann, S. Herrig, I. H. Bell, R. Beckmüller, E. W. Lemmon, M. Thol, R. Span,
> "EOS-CG-2021: A Mixture Model for the Calculation of Thermodynamic Properties of
> CCS Mixtures", Int. J. Thermophys. 44, 178 (2023). doi:10.1007/s10765-023-03263-6
> (BibTeX key `Neumann-IJT-2023`).

Scope: the **13 non-amine** components (CO2, H2O, N2, O2, Ar, CO, H2, CH4, H2S, SO2,
HCl, Cl2, NH3). The three amines (MEA, DEA, MDEA) are out of scope — CoolProp has no
pure-fluid EOS for them.

## Source data (provenance)

- `eoscg2021_table4.json` — Table 4 reducing parameters (β_T, γ_T, β_v, γ_v, F_ij),
  keyed by the paper's binary order `i-j`.
- `eoscg2021_table5.json` — Table 5 departure functions, each term tagged
  `pol`/`exp`/`spec`/`GBS`.

### Departure `type` mapping (Table 5 → CoolProp)

| paper term | meaning | CoolProp handling |
|---|---|---|
| `pol` | `n·δ^d·τ^t` | power term, `l=0` |
| `exp` | `·exp(−δ^l)` | power term, `l>0` |
| `spec` | `·exp(−η(δ−ε)²−β(δ−γ))` | GERG-2008 Gaussian (β linear in δ) |
| `GBS`  | `·exp(−η(δ−ε)²−β(τ−γ)²)` | modern Gaussian (β² in τ) |

`Npower` = count(`pol`)+count(`exp`); power terms precede Gaussian terms. CoolProp
`type` = `GERG-2008` (if Gaussian terms are `spec`), `Gaussian+Exponential` (if `GBS`),
or `Exponential` (no Gaussian terms).

## Redundancy — already in CoolProp, left untouched (cite original)

GERG-2008 (`Kunz-JCED-2012`) and EOS-CG-2016 (`Gernert-*`) pairs whose parameters
EOS-CG-2021 keeps unchanged: verified coefficient-identical and NOT modified, e.g.
N2-CO2, CH4-CO2, CH4-N2, CO2-H2O, N2/O2/CO-H2O (GeneralizedAirWater), and the 20
reducing-only GERG/EOS-CG pairs.

H2 pairs (CO2-H2, N2-H2, CO-H2, CH4-H2; Beckmüller `[7]`) are handled in the separate
Beckmüller-2021 PR (#3278) and are NOT touched here.

## Source-transcription caveats found

The paper's tables (transcribed by hand) contained two typos vs the validated
CoolProp/Beckmüller values — both in functions CoolProp already has correct:
- N2-H2 `gamma[7]`: 2.3 in the transcription vs 0.23 (validated to ppm in #3278).
- CH4-N2 `n[5]`: sign flip vs the GERG-2008 value already in CoolProp.
Because the seven genuinely-new departure functions cannot be cross-checked against
CoolProp, each is validated by tracing its VLE envelope with `pxy_teqp.py` and
comparing to the paper's figures / literature VLE data before the PR.

## Worklist (per-pair commits)

Reducing-only updates (Bell-2016 estimates → EOS-CG-2021 fits): CO2-SO2, SO2-O2,
SO2-CH4, SO2-HCl, Cl2-HCl. Reducing-only adds: SO2-H2O, SO2-N2, Cl2-SO2, CO-NH3,
N2-NH3, CH4-NH3, NH3-O2. Departure add/update: CO2-Ar (`GERG-2008`), CO2-CO,
CH4-H2O, H2S-H2O, NH3-H2O, Ar-NH3, H2-NH3 (`Gaussian+Exponential`/`Exponential`).
Low-priority predictive (linear combining rules): SO2-Ar, SO2-CO, SO2-H2, HCl-H2,
Cl2-H2, HCl-Ar, O2-HCl.

Validation tooling: `pxy_teqp.py` (teqp continuation VLE tracer reading these JSON
parameters).
