# REFPROP 10.1 Viscosity Correlations That Supersede CoolProp's — Audit and Plan

**Date:** 2026-08-29
**Status:** Planning. No code changes proposed in this PR — this document is the
audit that the follow-on work is scoped from.
**Beads:** epic `CoolProp-rip4`

**Goal:** Enumerate every pure-fluid viscosity correlation shipped in the REFPROP
10.1 beta that is newer than, or fills a gap against, what `dev/fluids/*.json`
carries today — with a **verified access status for every reference**, and an
open-access route wherever one exists, so as much of the implementation work as
possible can start without a subscription.

## Source of truth and caveats

- **REFPROP fluid files:** `REFPROPv10.1 beta test invitation20260720045933/extracted/app/FLUIDS`,
  161 pure-fluid `.FLD` files plus 5 pseudo-pure `.PPF` files, dated
  2026-07-12/19. **This is a beta.** Every finding below
  should be re-confirmed against the 10.1 general release before its coefficients
  are transcribed into `dev/fluids/`.
- The 10.0 set at `~/REFPROP10/FLUIDS` (148 fluids, 2023-04) is materially
  different and produces a *much* shorter list. Do not audit against it.
- **CoolProp side:** `TRANSPORT.viscosity.BibTeX` in `dev/fluids/*.json`,
  resolved against `CoolPropBibTeXLibrary.bib`. Fluids matched by
  `INFO.REFPROP_NAME`, falling back to `INFO.NAME` / `INFO.ALIASES`.
- Only REFPROP's **`#ETA`** blocks count — the dedicated, NIST-recommended
  correlations. Fluids whose primary transport model is `#TRN`/`@TRN` (extended
  corresponding states) are excluded as a **deliberate scope decision**: ECS is
  unpublished fitting on both sides, so there is no published correlation to
  adopt. This is *not* a claim that the two implementations agree. They do not
  in general — `dev/fluids/R1132a.json` records that CoolProp's ECS implements
  neither REFPROP's BIG large-molecule dilute-gas correction (worked around
  there by folding a constant 0.934 into `sigma_eta`) nor the initial-density
  Rainwater–Friend correction, leaving residual deviations up to ~4 % at
  moderate vapour densities. Closing that gap is a separate piece of work from
  this audit and is not scoped here. 61 of the 161 `.FLD` fluids have an
  `#ETA` block. The
  `.PPF` files were checked separately: Air still carries Lemmon & Jacobsen
  (2004), unchanged from 10.0 and matching CoolProp.
- **Counting rule:** §A and §B count only fluids whose **active** REFPROP model
  is a dedicated `#ETA` correlation. A published correlation that REFPROP ships
  as an `@ETA` *alternative* to an ECS primary does not enter the counts — it is
  called out narratively instead (§C-bis is the one such case). Adopting a model
  REFPROP itself does not use as primary is a different, weaker argument, and
  mixing the two would make the totals mean two things at once.
- `#` marks the active model and `@` an alternative, and the two symbols can be
  swapped within a fluid file. The audit therefore reads the `#ETA` block, not
  the first `ETA` block in the file — for PROPANE and ETHANOL those are not the
  same block, and for PROPANE that changes the answer (see Cautions).
- Open-access status was checked per DOI against Crossref, OpenAlex, the PMC ID
  converter, and the NIST public-access repository. Every NIST-hosted PDF linked
  below was downloaded and its first page checked against the expected title;
  PMC IDs come from the NCBI ID converter; the IAPWS and Imperial Spiral records
  were fetched directly. All DOIs were verified against Crossref.

## A. Supersessions — CoolProp has a model, REFPROP 10.1 has a newer published one

| Fluid | CoolProp today | REFPROP 10.1 (`#ETA`) | Open access |
|---|---|---|---|
| Argon | Lemmon & Jacobsen 2004 | `VS7` Sotiriadou et al., *IJT* **46**:133 (2025), [`10.1007/s10765-025-03603-8`](https://doi.org/10.1007/s10765-025-03603-8) | **Yes** — Springer hybrid OA; [PMC12241276](https://pmc.ncbi.nlm.nih.gov/articles/PMC12241276/) |
| Nitrogen | Lemmon & Jacobsen 2004 | `VS7` Huber, Perkins & Lemmon, *IJT* **45**:146 (2024), [`10.1007/s10765-024-03440-1`](https://doi.org/10.1007/s10765-024-03440-1) | **Yes** — Springer hybrid OA; [PMC11466908](https://pmc.ncbi.nlm.nih.gov/articles/PMC11466908/) |
| Methane | Quiñones-Cisneros & Deiters 2006 (f-theory) | `VS7` Sotiriadou et al., *IJT* **47**:18 (2025), [`10.1007/s10765-025-03690-7`](https://doi.org/10.1007/s10765-025-03690-7) | **Yes** — Springer hybrid OA; [PMC12686087](https://pmc.ncbi.nlm.nih.gov/articles/PMC12686087/) |
| Ethanol | Kiselev et al. 2005 | `VS7` Sotiriadou et al., *IJT* **44**:40 (2023), [`10.1007/s10765-022-03149-z`](https://doi.org/10.1007/s10765-022-03149-z) | **Yes** — Springer hybrid OA; [NIST pub_id=935954](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=935954) |
| Heavy water | IAPWS 2007 (i.e. the 1994 formulation) | `VS0` Assael et al., *JPCRD* **50**:033102 (2021), [`10.1063/5.0048711`](https://doi.org/10.1063/5.0048711) | **Yes, via IAPWS** — paper is paywalled, but [IAPWS R17-20](https://iapws.org/technical-guidance/release/D2Ovisc) carries the complete formulation, free |
| Ammonia | Fenghour et al. 1995 | `VS1` Monogenidou, Assael & Huber, *JPCRD* **47**:023102 (2018), [`10.1063/1.5036724`](https://doi.org/10.1063/1.5036724) | **Yes** — [PMC6512859](https://pmc.ncbi.nlm.nih.gov/articles/PMC6512859/) |
| R-134a | Huber, Laesecke & Perkins 2003 | `VS7` Velliadou, Assael & Huber, *IJT* **43**:105 (2022), [`10.1007/s10765-022-03029-6`](https://doi.org/10.1007/s10765-022-03029-6) | **Yes** — [NIST pub_id=934107](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=934107) |
| R-32 | Bell ρ<sub>s</sub>r-CS (2016) + ECS | `VS7` Velliadou et al., *IJT* **43**:129 (2022), [`10.1007/s10765-022-03050-9`](https://doi.org/10.1007/s10765-022-03050-9) | **Yes** — [NIST pub_id=934726](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=934726) |
| Ethylbenzene | ECS, Huber NIST RP-912 (predictive) | `VS6` Meng, Cao, Wu & Vesovic, *JPCRD* **46**:013101 (2017), [`10.1063/1.4973501`](https://doi.org/10.1063/1.4973501) | **Yes** — accepted manuscript in [Imperial Spiral](https://spiral.imperial.ac.uk/handle/10044/1/43620) |
| R-1234yf | Bell ρ<sub>s</sub>r-CS + ECS | `VS7` Huber & Assael, *Int. J. Refrig.* **71**:39–45 (2016), [`10.1016/j.ijrefrig.2016.08.007`](https://doi.org/10.1016/j.ijrefrig.2016.08.007) | **Yes** — [PMC5103321](https://pmc.ncbi.nlm.nih.gov/articles/PMC5103321/) |
| R-1234ze(E) | Bell ρ<sub>s</sub>r-CS + ECS | same Huber & Assael (2016) paper | **Yes** — as above |
| R-245fa | Bell ρ<sub>s</sub>r-CS + ECS | `VS1` Perkins, Huber & Assael, *JCED* **61**:3286–3294 (2016), [`10.1021/acs.jced.6b00350`](https://doi.org/10.1021/acs.jced.6b00350) | **Yes** — [NIST pub_id=920723](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=920723) (accepted manuscript, `.docx`) |
| Cyclopentane | Chung et al. 1988 (generic) | `VS1` Tasidou, Huber & Assael, *JPCRD* **48**:043101 (2019), [`10.1063/1.5128321`](https://doi.org/10.1063/1.5128321) | **No** — publisher only; NIST landing page has no hosted copy |
| Ethane | Friend et al. 1991 | `VS7` Herrmann, Hellmann & Vogel, *JPCRD* **47**:023103 (2018), [`10.1063/1.5037239`](https://doi.org/10.1063/1.5037239) | **No** — publisher only (no NIST author) |
| n-Butane | Vogel et al. 1999 | `VS7` Herrmann & Vogel, *JPCRD* **47**:013104 (2018), [`10.1063/1.5020802`](https://doi.org/10.1063/1.5020802) | **No** — publisher only (no NIST author) |

12 of the 15 have a legitimate free route. Only cyclopentane, ethane and
n-butane need library access.

**Helium is deliberately excluded.** REFPROP 10.1 replaces Arp et al. (1998)
with "Huber, M.L., Friend, D.G., and Lemmon, E.W., *Reference Correlation for
the Viscosity of Helium*, Int. J. Thermophys., **for submission 2026**" — no
DOI, not published. Nothing to adopt yet; revisit when it appears.

## B. Gaps — CoolProp has no viscosity model, REFPROP 10.1 has a published one

| Fluid | REFPROP 10.1 (`#ETA`) | Open access |
|---|---|---|
| Krypton | `VS7` Polychroniadou, Antoniadis, Assael & **Bell**, *IJT* **43**:6 (2021), [`10.1007/s10765-021-02927-5`](https://doi.org/10.1007/s10765-021-02927-5) — entropy scaling | **Yes** — [NIST pub_id=933039](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=933039) |
| Xenon | `VS7` Velliadou et al., *IJT* **42**:74 (2021), [`10.1007/s10765-021-02818-9`](https://doi.org/10.1007/s10765-021-02818-9) | **Yes** — [PMC8356199](https://pmc.ncbi.nlm.nih.gov/articles/PMC8356199/) |
| Ethylene | `VS7` Sotiriadou et al., *IJT* **45**:87 (2024), [`10.1007/s10765-024-03378-4`](https://doi.org/10.1007/s10765-024-03378-4) | **Yes** — [NIST pub_id=957842](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=957842) |
| Deuterium | `VS7` Muzny, Huber & Kazakov, *JCED* **58**:969–979 (2013) — D₂ parameterization, [`10.1021/je301273j`](https://doi.org/10.1021/je301273j) | **Yes** — [NIST pub_id=912625](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=912625) |
| n-Undecane | `VS7` Assael, Papalas & Huber, *JPCRD* **46**:033103 (2017), [`10.1063/1.4996885`](https://doi.org/10.1063/1.4996885) | **Yes** — [PMC5721360](https://pmc.ncbi.nlm.nih.gov/articles/PMC5721360/); also [NIST SRD reprint](https://srd.nist.gov/jpcrdreprint/1.4996885.pdf) |
| R-161 | `VS1` Tsolakidou, Assael, Huber & Perkins, *JPCRD* **46**:023103 (2017), [`10.1063/1.4983027`](https://doi.org/10.1063/1.4983027) | **Yes** — [PMC5544035](https://pmc.ncbi.nlm.nih.gov/articles/PMC5544035/) |
| Novec 649 | `VS1` Wen, Meng, Huber & Wu, *JCED* **62**:3603–3609 (2017), [`10.1021/acs.jced.7b00572`](https://doi.org/10.1021/acs.jced.7b00572) | **Yes** — [PMC5755718](https://pmc.ncbi.nlm.nih.gov/articles/PMC5755718/) |
| Propylene glycol | `VS7` Velliadou et al., *IJT* **43**:42 (2022), [`10.1007/s10765-021-02970-2`](https://doi.org/10.1007/s10765-021-02970-2) | **Yes** — [NIST pub_id=933662](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=933662) |
| Tetrahydrofuran | `VS7` Sotiriadou et al., *IJT* **45**:123 (2024), [`10.1007/s10765-024-03415-2`](https://doi.org/10.1007/s10765-024-03415-2) | **Yes** — [NIST pub_id=958095](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=958095) |
| Neon | `VS7` Sotiriadou et al., *IJT* **47**:77 (2026), [`10.1007/s10765-026-03745-3`](https://doi.org/10.1007/s10765-026-03745-3) | **No copy found** — published 2026; NIST deposit likely not up yet, recheck later |

9 of 10 are freely obtainable today.

Note: `PropyleneGlycol.json` and `Tetrahydrofuran.json` both carry
`"REFPROP_NAME": "N/A"` even though REFPROP 10.1 has them as `PROPYLENEGLYCOL`
and `THF`. Worth fixing alongside — it also unblocks the REFPROP backend for
those fluids.

## C. Cautions — things that look like findings but are not

- **Propane is not a supersession.** REFPROP 10.0 had Vogel & Herrmann (2016) as
  the active model. In the 10.1 beta the markers are swapped: `#ETA` is Vogel et
  al. (1998) and the 2016 correlation sits at `@ETA`. CoolProp already uses
  Vogel et al. 1998, so it matches 10.1's active model. **Confirm this swap is
  intentional with NIST before acting on it either way** — it is exactly the kind
  of thing a beta gets wrong, and it is also exactly the kind of thing that gets
  deliberately reverted.
- **R-23 went backwards in REFPROP.** 10.0 had Shan et al. (2000) as `#ETA`;
  10.1 has only ECS. CoolProp's hardcoded Shan-2000 model is the better of the
  two. No action.
- **Hydrogen sulfide:** CoolProp's Quiñones-Cisneros et al. (2012) is newer than
  REFPROP's Schmidt et al. (2008). No action.
- **n-Nonane citation only.** REFPROP cites Huber et al., *FPE* (2005) where
  CoolProp cites the 2004 paper, but the coefficients are identical (checked
  term by term: `2.66987`, `1.32137`, `-0.0314367`, …). A bibliography fix at
  most.
- **Oxygen was not updated** while argon and nitrogen were — 10.1 still has
  Lemmon & Jacobsen (2004), same as CoolProp.
- Unchanged and already matching: CO₂, water, benzene, toluene, o-/m-/p-xylene,
  cyclohexane, hexane, heptane, octane, decane, dodecane, DME, methanol,
  hydrogen, parahydrogen, isobutane, R-123, R-125, R-116, SF₆.
- **R-113 is not a gap.** CoolProp has no viscosity for it, but REFPROP's `#ETA`
  block is explicitly "estimation based on the R-134a model of Huber et al.
  (2003), scaled" — an estimate, not a published correlation for R-113. Same for
  R-116, which CoolProp already covers with the equivalent ECS model.
- Everything CoolProp models with Chung-1988 or ECS that is *not* listed in §A
  is also ECS-or-unpublished in REFPROP 10.1 (R-11, R-12, R-13, R-14, R-22,
  R-124, R-141b, R-142b, R-143a, R-152a, R-218, R-227ea, R-236ea, R-236fa,
  RC318, isopentane, n-pentane, propylene). No published upgrade exists there.

## C-bis. R-1233zd(E) — a candidate the `#ETA` rule hides

REFPROP 10.1's *primary* transport model for R-1233zd(E) is ECS, so the fluid
does not appear in §B. But the file carries a published, measurement-based
viscosity correlation as an `@ETA` alternative: Meng, Wen & Wu,
*J. Chem. Thermodyn.* **123**:140–145 (2018),
[`10.1016/j.jct.2018.04.001`](https://doi.org/10.1016/j.jct.2018.04.001) — no
open-access copy found.

CoolProp's side of this moved while the audit was being written. PR #3335
restored the `TRANSPORT` block that #2768 had dropped and refit the ρ<sub>s</sub>r-CS
constants (`C = 0.8089`, `rhosr_critical` recomputed for the Akasaka & Lemmon
2022 EOS) against the 61 primary liquid points of Miyara, Alam & Kariya,
*Int. J. Refrig.* **92**:86–93 (2018). That refit reports AAD 2.62 % / bias
−0.07 % over 314–433 K, against 2.78 % for the REFPROP 10 ECS model on the same
states — so where data exist, CoolProp is now **at least as good as what REFPROP
actually ships as primary**, and this is no longer an obvious supersession.

It is still worth doing, for one specific reason. The refit note records that
over **243–313 K** — below Miyara's range, where only paywalled measurements
exist — the CoolProp model sits 6.7 % above REFPROP's ECS, and that saturated
vapour is ~10 % high regardless of `C` because the dilute-gas Lennard-Jones
parameters are Chung estimates. The Meng, Wen & Wu correlation is fitted to
measurements across 243–373 K and would replace a corresponding-states
extrapolation with data in exactly the range CoolProp is weakest. Treat this as
a **lower-priority, data-driven improvement**, not a supersession, and benchmark
against the Miyara points before switching — the current model is not broken.

## D. Out of scope — REFPROP 10.1 fluids CoolProp does not have

These carry `#ETA` correlations but there is no CoolProp fluid to attach them
to: 1-hexene, n-hexadecane, ethylene glycol, methyldiethanolamine, and the
lubricants MIL-PRF-23699, POE5, POE7, POE9.

## Implementation notes

**Model-form distribution across §A + §B:** 18 × `VS7`, 5 × `VS1`, 1 × `VS6`,
1 × `VS0`.

- **`VS1`** (ammonia, cyclopentane, R-161, Novec 649, R-245fa) is the
  collision-integral dilute gas + Rainwater–Friend second viscosity virial +
  modified Batschinski–Hildebrand residual. All three pieces are already Tier-A
  forms in the transport expression DSL
  (`docs/superpowers/specs/2026-06-12-transport-expression-dsl-design.md`, shipped
  in PR #3185), so these are **data-only additions** to `dev/fluids/*.json` — no
  C++ needed. Start here.
- **`VS7`** is REFPROP's generic form, encoded as an RPN token stream
  (`$DG`/`$VV`/`$RF`/`$CF` sections with `SUMLOGT`, `SUM:n`, `EXP` and a postfix
  operator stack). This is the same NIST RPN encoding the DSL design doc cites as
  its motivating prior art. The interesting option is a **`VS7` → CoolProp-DSL
  transpiler** rather than 18 hand-written model blocks; that decision deserves
  its own brainstorm before anyone starts transcribing.
- **`VS0`** (heavy water) is a hardcoded fluid-specific routine on both sides.
  CoolProp already has a `HeavyWater` hardcoded viscosity path; adopting IAPWS
  R17-20 means rewriting that routine, not adding data.
- **`VS6`** (ethylbenzene) is a one-off; check whether it reduces to a Tier-A
  form before adding a C++ routine.

**Validation.** Every fluid touched needs a `[refprop]`-tagged Catch2 comparison
against the REFPROP backend. REFPROP runs locally on this machine
(`COOLPROP_REFPROP_ROOT=/Users/ianbell/REFPROP10`) — but that env var points at
the **10.0** install, which carries *different* viscosity models for most of §A.

**These tests must assert the REFPROP version, not just its availability.** The
usual `REFPROP_supported()` / skip-if-unavailable guard fails open here: against
a 10.0 install the backend answers happily, from the older correlation, and the
test either passes for the wrong reason or fails and gets written off as noise.
Neither outcome validates anything. Gate on
`get_global_param_string("REFPROP_version")` (wired up at
`src/CoolProp.cpp:1202`) and **fail** — do not skip — when the loaded library is
not 10.1, since at that point REFPROP is present and the test simply cannot mean
what it claims. Keep the ordinary skip for REFPROP being absent entirely.

Note that this is currently blocked: the 10.1 beta ships Windows DLLs only, so a
macOS/Linux 10.1 shared library is a prerequisite for running any of it.

## Suggested sequencing

1. Confirm the propane marker swap and the helium submission status with NIST.
2. Land the five `VS1` fluids as pure JSON data + tests (ammonia, R-245fa,
   R-161, Novec 649, cyclopentane — cyclopentane last, it needs library access).
3. Fix `REFPROP_NAME` for `PropyleneGlycol` and `Tetrahydrofuran`.
4. Brainstorm the `VS7` question (transpiler vs. hand-written), then do the
   noble gases (krypton, xenon, neon) as the first batch — small, clean,
   entropy-scaling-based, and krypton is our own correlation.
5. Heavy water via IAPWS R17-20.
6. The rest of `VS7` in order of user impact: nitrogen, argon, methane, R-134a,
   R-32, ethane, n-butane, ethanol, ethylene.
