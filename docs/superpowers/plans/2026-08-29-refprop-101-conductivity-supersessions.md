# REFPROP 10.1 Thermal Conductivity Correlations That Supersede CoolProp's — Audit and Plan

**Date:** 2026-08-29
**Status:** Planning. No code changes proposed in this PR.
**Beads:** epic `CoolProp-rip4`
**Companion:** `2026-08-29-refprop-101-viscosity-supersessions.md` — same method,
same sources, same caveats. Read that document's "Source of truth and caveats"
section first; only the differences are restated here.

**Goal:** The viscosity audit, repeated for thermal conductivity.

## Method delta

- Reads REFPROP's **`#TCX`** blocks (dedicated conductivity correlations).
  70 of the 161 `.FLD` fluids in the 10.1 beta have one — nine more than have a
  dedicated `#ETA`.
- CoolProp side is `TRANSPORT.conductivity.BibTeX` in `dev/fluids/*.json`.
- Fluids whose primary transport model is `#TRN`/`@TRN` (ECS) are excluded as a
  scope decision, for the reason given in the companion document — there is no
  published correlation to adopt, *not* because the two ECS implementations
  agree. §A and §B likewise count only fluids whose **active** REFPROP model is
  a dedicated `#TCX` correlation. The eleven ECS-primary fluids that carry a `@TCX` secondary were
  checked individually: all are Chung-1988 or NIST14 estimates, so none is an
  upgrade over what CoolProp already does.

**Headline:** conductivity is in much better shape than viscosity. CoolProp
already carries the current Assael/Huber/Perkins reference correlations for
benzene, toluene, the xylenes, ethylbenzene, hexane, heptane, the pentanes,
methanol, ethanol, hydrogen, parahydrogen, CO₂, water, SF₆ and R-1234yf/ze(E).
Only 6 genuine supersessions, against 15 on the viscosity side.

## A. Supersessions — CoolProp has a model, REFPROP 10.1 has a newer published one

| Fluid | CoolProp today | REFPROP 10.1 (`#TCX`) | Open access |
|---|---|---|---|
| Ammonia | Tufeu et al. 1984 | `TC1` Monogenidou, Assael & Huber, *JPCRD* **47**:043101 (2018), [`10.1063/1.5053087`](https://doi.org/10.1063/1.5053087) | **Yes** — [NIST pub_id=925858](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=925858) |
| Heavy water | IAPWS 2007 | `TC0` Huber et al., *JPCRD* **51**:013102 (2022), [`10.1063/5.0084222`](https://doi.org/10.1063/5.0084222) | **Yes, via IAPWS** — paper paywalled, but [IAPWS R18-21](https://iapws.org/public/documents/SRPag/D2OThCond.pdf) carries the complete formulation, free |
| Nitrogen | Lemmon & Jacobsen 2004 | `TC7` Sotiriadou, Assael & Huber, *IJT* **46**:42 (2025), [`10.1007/s10765-025-03516-6`](https://doi.org/10.1007/s10765-025-03516-6) | **Yes** — Springer hybrid OA |
| Propylene | ECS, Huber et al. 2003 | `TC1` Assael, Koutian, Huber & Perkins, *JPCRD* **45**:033104 (2016), [`10.1063/1.4958984`](https://doi.org/10.1063/1.4958984) | **Yes** — [PMC5094801](https://pmc.ncbi.nlm.nih.gov/articles/PMC5094801/) |
| R-245fa | ECS, Huber et al. 2003 | `TC1` Perkins, Huber & Assael, *JCED* **61**:3286–3294 (2016), [`10.1021/acs.jced.6b00350`](https://doi.org/10.1021/acs.jced.6b00350) | **Yes** — [NIST pub_id=920723](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=920723) |
| R-143a | McLinden et al. 2000 | `TC1` Huber, **NISTIR 8209** (2018), [`10.6028/NIST.IR.8209`](https://doi.org/10.6028/NIST.IR.8209) — REFPROP labels this **"preliminary"** | **Yes** — [nvlpubs](https://nvlpubs.nist.gov/nistpubs/ir/2018/NIST.IR.8209.pdf) |

All six are freely obtainable. The R-245fa paper is the same one that supplies
the R-245fa *viscosity* supersession — one PDF covers both properties, so do
them together.

## B. Gaps — CoolProp has no conductivity model, REFPROP 10.1 has a published one

| Fluid | REFPROP 10.1 (`#TCX`) | Open access |
|---|---|---|
| Cyclohexane | `TC1` Koutian, Assael, Huber & Perkins, *JPCRD* **46**:013102 (2017), [`10.1063/1.4974325`](https://doi.org/10.1063/1.4974325) | **Yes** — [PMC5455799](https://pmc.ncbi.nlm.nih.gov/articles/PMC5455799/) |
| Ethylene | `TC1` Assael, Koutian, Huber & Perkins, *JPCRD* **45**:033104 (2016), [`10.1063/1.4958984`](https://doi.org/10.1063/1.4958984) | **Yes** — [PMC5094801](https://pmc.ncbi.nlm.nih.gov/articles/PMC5094801/) |
| Xenon | `TC7` Velliadou, Assael, Antoniadis & Huber, *IJT* **42**:51 (2021), [`10.1007/s10765-021-02803-2`](https://doi.org/10.1007/s10765-021-02803-2) | **Yes** — [NIST pub_id=931647](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=931647) |
| n-Undecane | `TC1` Assael, Papalas & Huber, *JPCRD* **46**:033103 (2017), [`10.1063/1.4996885`](https://doi.org/10.1063/1.4996885) | **Yes** — [PMC5721360](https://pmc.ncbi.nlm.nih.gov/articles/PMC5721360/) |
| R-161 | `TC1` Tsolakidou, Assael, Huber & Perkins, *JPCRD* **46**:023103 (2017), [`10.1063/1.4983027`](https://doi.org/10.1063/1.4983027) | **Yes** — [PMC5544035](https://pmc.ncbi.nlm.nih.gov/articles/PMC5544035/) |
| Novec 649 | `TC1` Perkins, Huber & Assael, *JCED* **63**:2783–2789 (2018), [`10.1021/acs.jced.8b00132`](https://doi.org/10.1021/acs.jced.8b00132) | **Yes** — [PMC6513336](https://pmc.ncbi.nlm.nih.gov/articles/PMC6513336/) |
| R-1233zd(E) | `TC1` Perkins, Huber & Assael, *JCED* **62**:2659–2665 (2017), [`10.1021/acs.jced.7b00106`](https://doi.org/10.1021/acs.jced.7b00106) | **Yes** — [PMC5721355](https://pmc.ncbi.nlm.nih.gov/articles/PMC5721355/) |
| R-1336mzz(Z) | `TC1` Perkins et al., *IJT* **41**:103 (2020), [`10.1007/s10765-020-02681-0`](https://doi.org/10.1007/s10765-020-02681-0) | **Yes** — [PMC7727279](https://pmc.ncbi.nlm.nih.gov/articles/PMC7727279/) |
| Tetrahydrofuran | `TC1` Sotiriadou et al., *IJT* **45**:123 (2024), [`10.1007/s10765-024-03415-2`](https://doi.org/10.1007/s10765-024-03415-2) | **Yes** — [NIST pub_id=958095](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=958095) |
| Methyl oleate | `TC1` Perkins & Huber, *Energy & Fuels* **25**:2383–2388 (2011), [`10.1021/ef200417x`](https://doi.org/10.1021/ef200417x) | **Yes** — [NIST pub_id=907397](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=907397) |
| Methyl linoleate | same Perkins & Huber (2011) paper | **Yes** — as above |
| Neon | `TC7` Assael, Sotiriadou, Thol & Huber, *IJT* **47**:126 (2026), [`10.1007/s10765-026-03782-y`](https://doi.org/10.1007/s10765-026-03782-y) | **Yes** — [NIST pub_id=961954](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=961954) (accepted manuscript, `.docx`) |

All 12 are freely obtainable as of 2026-08-29. Availability is a moving target
for the 2026 papers — neon's NIST deposit went up after this audit's first pass,
and the neon *viscosity* paper in the companion document still has no free copy.
Re-check before concluding anything is unobtainable.

Seven of these fluids appear in *both* audits — ethylene, xenon, neon,
n-undecane, R-161, THF and Novec 649 need viscosity *and* conductivity. For
n-undecane, R-161 and THF a single paper supplies both; **Novec 649 needs two**
(viscosity Wen et al. 2017, conductivity Perkins et al. 2018). R-245fa, in §A
rather than here, is the other single-paper pair. See the sequencing section —
the shared-paper fluids are not all equally cheap, because n-undecane and THF
carry `VS7` viscosity.

## C. Cautions — things that look like findings but are not

- **R-32 conductivity is not a supersession.** REFPROP's `#TCX` is
  "Perkins and Huber (2005) **(unpublished)**". CoolProp's ECS model is not
  obviously worse and there is nothing published to adopt. Note that this is the
  opposite of R-32 *viscosity*, which is a clean supersession.
- **Deuterium conductivity is not a published correlation.** REFPROP's block
  says "unpublished; based on scaling the Assael correlation" for normal
  hydrogen (*JPCRD* **40**:033101, 2011, [`10.1063/1.3606499`](https://doi.org/10.1063/1.3606499),
  free at [NIST](https://www.nist.gov/document/jpcrd402011corrpdf)). Adoptable,
  but it is a scaled H₂ correlation and must be documented as such — CoolProp
  already has the underlying Assael-2011 model for hydrogen.
- **Methyl linolenate, palmitate and stearate are estimates.** REFPROP cites
  NISTIR 8209 and states plainly that the correlation "is an estimation, based
  on results for methyl oleate, adjusted". Only methyl **oleate** and methyl
  **linoleate** have measured, published correlations. Treat the other three as
  a lower tier, or skip them.
- **R-143a is labelled "preliminary"** by REFPROP itself. It is still newer and
  published, but it is not a reference-quality correlation.
- **R-116** — REFPROP's `#TCX` is the R-134a model applied to R-116, not a
  dedicated correlation; CoolProp's ECS is equivalent. No action.
- Already matching, no action: argon, oxygen, benzene, toluene, o-/m-/p-xylene,
  ethylbenzene, hexane, heptane, octane, nonane, decane, dodecane, n-pentane,
  isopentane, cyclopentane, methane, ethane, propane, n-butane, isobutane,
  methanol, ethanol, hydrogen, parahydrogen, helium, CO₂, water, SF₆, R-123,
  R-125, R-134a, R-152a, R-1234yf, R-1234ze(E).

## D. Out of scope — REFPROP 10.1 fluids CoolProp does not have

1-hexene, n-hexadecane, methylcyclohexane, propylcyclohexane, ethylene glycol,
R-E347mcc, MIL-PRF-23699, POE5, POE7, POE9.

## Implementation notes

**Model-form distribution across §A + §B:** 14 × `TC1`, 3 × `TC7`, 1 × `TC0`.

- **`TC1`** is the dilute-gas ratio-of-polynomials (or η₀-and-polynomial) term
  plus a residual polynomial, optionally with a `simplified_Olchowy_Sengers`
  critical enhancement. `dilute_ratio_polynomials`, `dilute_eta0_and_poly`,
  `residual_polynomial` and `residual_polynomial_and_exponential` are all Tier-A
  forms in the transport expression DSL
  (`docs/superpowers/specs/2026-06-12-transport-expression-dsl-design.md`), so the
  bulk of this work is **data-only**. The critical enhancement is Tier B — it
  needs EOS-derived scalars and stays hardcoded, but CoolProp already has that
  routine and it is shared, not per-fluid.
- **`TC7`** (nitrogen, xenon, neon) is the RPN-encoded generic form, same
  situation as `VS7` in the viscosity document — see that document's
  implementation notes; a transpiler would serve both properties.
- **`TC0`** (heavy water) is hardcoded on both sides. CoolProp's existing
  `HeavyWater` conductivity routine would be rewritten against IAPWS R18-21.

**Do heavy water once, not twice.** IAPWS R17-20 (viscosity) and R18-21
(thermal conductivity) both supersede what CoolProp's single
`IAPWS-D2O-2007-Transport` entry covers, and both hardcoded routines live
together. One change, both properties.

## Suggested sequencing

1. **R-245fa and R-161** — one paper covers both properties *and* both models
   are `VS1`/`TC1`, so each fluid is a single data-only change closing two gaps.
   These two are the cheapest work in either audit; start here.
   Then, in decreasing order of convenience:
   - **n-undecane and THF** — one paper each covers both properties, but their
     viscosity is `VS7`, so only the `TC1` half is data-only today; the
     viscosity half waits on the `VS7` decision.
   - **Novec 649** — `VS1`/`TC1`, both data-only, but it needs *two* papers:
     viscosity from Wen et al., *JCED* **62**:3603 (2017) and conductivity from
     Perkins et al., *JCED* **63**:2783 (2018). Both are on PMC.
2. **Heavy water** — both IAPWS releases, one hardcoded pair.
3. **Ammonia, cyclohexane, ethylene, propylene, R-1233zd(E), R-1336mzz(Z),
   methyl oleate, methyl linoleate** — `TC1`, data-only, all open access.
   R-1233zd(E) is worth prioritising: PR #3335 restored and refit its
   **viscosity**, but `dev/fluids/R1233zd(E).json` still has no `conductivity`
   entry at all, so this is the one remaining transport gap for that fluid — and
   the Perkins/Huber/Assael correlation for it is on PMC.
4. **R-143a** — cheap, but label it preliminary in the JSON note.
5. **Nitrogen, xenon, neon** — blocked behind the `TC7`/`VS7` decision.
