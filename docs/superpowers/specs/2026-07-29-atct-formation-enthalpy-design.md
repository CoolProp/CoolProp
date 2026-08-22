# Standard ideal-gas enthalpy of formation from ATcT — Design

**Date:** 2026-07-29
**Status:** Approved design, ready for implementation planning
**Epic:** `CoolProp-oett` (Standard-state thermochemical data)
**Issues:** `CoolProp-ijbw` (regenerator), `CoolProp-ihto` (JSON ingestion + C++ API)
**Upstream:** GitHub issue [#1360](https://github.com/CoolProp/CoolProp/issues/1360) ("Chemical reactions?", open since 2016-11-30)

---

## 1. Goal

Store the **standard molar enthalpy of formation of the ideal gas at 298.15 K**,
Δ<sub>f</sub>H°(298.15 K), as per-fluid metadata in `dev/fluids/*.json`, and expose it as
a trivial parameter through `AbstractState` / `PropsSI`.

Ship it with a **regeneration script that derives every value directly from the Active
Thermochemical Tables**, so no number in the repository rests on a hand-curated
spreadsheet.

This is Tier 1 of a larger arc. It deliberately does **not** compute heating values,
does **not** add a reference state, and does **not** store standard entropies (§8).

## 2. Why now

Issue #1360 has been open for ten years with exactly one blocker: nobody produced the
data. Maintainers requested it in 2016, 2019 and 2020, and the companion issue #58 was
bulk-closed for inactivity. On 2026-07-27 @Tobias-Reiter attached a curated CSV of 75
ATcT values plus a REFPROP cross-check program, which unblocks the data question.

We take the contributor's *dataset choice* (ATcT TN 1.220) and their REFPROP comparison
as the valuable contribution, but **regenerate the numbers from source** rather than
vendoring the CSV. Section 3 records why that distinction earns its keep.

## 3. Key facts (verified 2026-07-29)

### 3.1 ATcT is scrapeable, and is the more reliable source

- The whole TN 1.220 table (3,442 species) is a single 4.3 MB HTML page at
  `https://atct.anl.gov/Thermochemical%20Data/version%201.220/index.php`, with
  machine-readable rows carrying CAS, name, formula+phase, Δ<sub>f</sub>H(0 K),
  Δ<sub>f</sub>H(298.15 K), a 95 % uncertainty, and an ATcT species ID:

  ```
  <tr id="s1_6n4_1c0 i23 CAS74-82-8">
    Name=Methane  Formula="CH4  (g)"  DHf0=-66.543  DHf298=-74.513
    Uncert="± 0.043"  Units=kJ/mol  ATcTID=74-82-8*0
  ```

- **REFPROP's `HEATFRMdll` is a derived quantity, not an independent measurement.**
  Regressing REFPROP's returned `hFrm` against the gross heating value stored in each
  `.FLD` header (`:Heat:`) recovers, over 14 hydrocarbons with **max residual
  0.000 kJ/mol**:

  | implied constant | fitted from REFPROP | CODATA |
  |---|---|---|
  | Δ<sub>f</sub>H(CO₂)    | −393.5110 kJ/mol | −393.51 |
  | Δ<sub>f</sub>H(H₂O, l) | −285.8300 kJ/mol | −285.83 |

  H₂S separately implies Δ<sub>f</sub>H(SO₂) = −296.81 kJ/mol. `hFrm` is independent of
  both T and D (identical at 298.15 K, 500 K, 1000 K and at 150 K / 20 mol/L), and
  mixtures are a plain mole-fraction average.

- **REFPROP's stored heating values are almost entirely uncited.** Of 148 `.FLD` files,
  all carry the `:Heat:` slot; 82 hold `????`; 65 hold a bare number with no source; and
  exactly **one** (methane) cites ISO 6976:2016.

- **The two sources disagree beyond ATcT's stated uncertainty.** Over the 41 overlapping
  species: RMS 0.98 kJ/mol, median |diff| 0.51 kJ/mol, worst **cyclobutene 4.34 kJ/mol =
  7.2σ**. Only 8/41 agree within ATcT 1σ and 23/41 within 2σ. REFPROP's docs define
  `HEATFRM`'s standard state as "298.15 K for the ideal gas", identical to ATcT's, so
  this is a genuine disagreement and not a basis mismatch.

  **Conclusion:** ATcT is the reference. This comparison is *source-selection evidence
  only* — REFPROP is not consulted, called, or asserted against anywhere in this tier
  (see §7).

### 3.2 Three parsing traps, each a silent-failure mode

1. **Element rows differ in shape.** Ar, N₂, O₂, He, Ne, Kr, Xe, F₂, Cl₂ report
   uncertainty `exact` with an empty units cell. A `±`-anchored regex drops all nine
   without error.
2. **CAS is not a unique key.** Argon appears as `Ar (g)` *and* `Ar (aq, undissoc)`;
   water has **nine** rows including `(l)`, `(cr)`, `(cr,l)`, `(l, eq.press.)`,
   `(g, para)` and `(g, ortho)`.
3. **Phase qualifiers are load-bearing.** `cis-2-Butene` exists only as `(g, cis)`, so an
   exact-`(g)`-only filter silently drops it.

### 3.3 Coverage ceiling is real

R134a (811-97-2) and R1234yf (754-12-1) are genuinely **absent** from ATcT — this is not
a parser gap. No scraper improvement recovers the HFC/HFO refrigerants. Expected
coverage is ~75 of 137 fluids.

### 3.4 The hand-curated CSV contains at least one wrong value

It lists `OrthoHydrogen` and `ParaHydrogen` as `0.000 ± 0.000`. ATcT publishes
`H2 (g, para)` = **−0.058 kJ/mol** (ID `1333-74-0*2`) and a separate ortho row. This is
the concrete justification for regenerating rather than vendoring.

### 3.5 Build path

`dev/fluids/*.json` → `dev/package_json.py` → `dev/all_fluids.json` →
`dev/generate_headers.py` → `dev/all_fluids.cbor` → incbin header. Regenerating fluid
data therefore requires a rebuild to take effect. `dev/validate_fluid_schemas.py` covers
only pcsaft / cubics / mixtures, **not** the main fluid files, so no schema file needs
updating.

## 4. Regenerator — `dev/atct/fetch_atct_formation.py`

Single-purpose script. Touches the network only when explicitly run.

**Fetch.** ATcT index for a version supplied as `--version 1.220` (URL templated, never
hardcoded). Record the SHA-256 of the retrieved page.

**Parse.** Each `<tr id="...CAS...">` into
`(CAS, name, formula+phase, ΔfH₀, ΔfH₂₉₈, uncertainty, ATcT ID)`. Two rules get their own
unit tests because §3.2 shows both fail silently:

- *Elements*: `uncertainty == "exact"` with empty units parses as `0 ± 0`, not as a
  parse failure and not as a skip.
- *Phase*: prefer an exact `(g)` row. If none exists, accept a **single**
  `(g, <qualifier>)` row. More than one candidate after this rule is an error.

**Bind.** CoolProp fluid → ATcT species on `INFO.CAS`, gas rows only. CoolProp's
suffixed CAS for spin isomers (`1333-74-0p`, `1333-74-0o`) strips the suffix and requires
the corresponding `(g, para)` / `(g, ortho)` row, so para-H₂ receives −58 J/mol rather
than a fabricated zero. **0 or >1 candidates → hard error; nothing is written.**

**Fail-loud ledger.** `dev/atct/expected_coverage.json` is committed and lists every
CoolProp fluid with its expected state — `matched` plus its ATcT ID, or `absent` plus a
reason. The script diffs reality against the ledger and **exits non-zero on any
difference**. A future ATcT version that renames or drops a species becomes a failing
run, not a silent gap.

**Write.** `INFO.STANDARD_STATE` into `dev/fluids/*.json`, converting kJ/mol → J/mol.
Also emits `dev/atct/atct_report.csv` (committed audit trail: fluid, value, uncertainty,
ATcT ID, source page SHA-256). The 4.3 MB source page is **not** committed.

## 5. Fluid file schema

New `INFO.STANDARD_STATE` block, styled after the existing `STATES.critical` convention.
Shared state at the top; **provenance attached per quantity**, so a later standard-entropy
tier drops a `smolar` sibling in without a second migration through 137 files and without
inheriting the enthalpy's provenance:

```json
"STANDARD_STATE": {
  "T": 298.15, "T_units": "K",
  "p": 100000.0, "p_units": "Pa",
  "phase": "ideal_gas",
  "hmolar_formation": {
    "value": -74513.0, "units": "J/mol", "uncertainty": 43.0,
    "source": "ATcT", "version": "1.220", "id": "74-82-8*0"
  }
}
```

An absent block, or an absent `hmolar_formation` key, means the property is unavailable.

### 5.1 The block is shared; ownership is per quantity

*(Added post-review.)* `STANDARD_STATE` is deliberately shared — §8 plans a standard-entropy
tier from a different source, dropping a `smolar` sibling into the same object. Two
consequences the original spec did not draw out, both of which produced real defects:

- **The reference state (`T`, `T_units`, `p`, `p_units`, `phase`) is shared by every
  quantity in the block.** Writing ATcT's over another source's re-annotates their value
  with a reference state it was never published at — worse than deleting it, because the
  result still reads as valid. The regenerator therefore fills only *absent* scaffolding
  and refuses the run outright on a conflict. That refusal is **not** overridable by
  `--allow-removals`: there is no correct automatic answer when two sources disagree about
  the reference state.
- **Ownership is decided per quantity, by its `source`.** The regenerator only ever writes,
  overwrites or removes an `hmolar_formation` whose `source` is `ATcT` (or absent, which
  predates the field). Anything else belongs to another tier and is left alone.

### 5.2 Removal semantics

*(Added post-review; the original spec had none, and writing was matched-only.)* A fluid
that flips `matched → absent` between ATcT versions must lose its value, or it keeps
serving a number the ledger and report both deny, stamped with the old version, forever.

Removal is destructive and every path into `absent` is quiet — the parser skips rows with
no formula button or no `DHf298` span, the only floor is one surviving row, and the page
cache is written before any sanity check. A truncated download therefore reports almost
everything absent. So:

- Removals and foreign-data overwrites are **computed before anything is written** and
  **refused by default**, exiting non-zero with the affected fluids on stderr.
- `--allow-removals` authorizes them. The refusal covers the fluid files, the report **and**
  `expected_coverage.json` — the ledger is the tripwire the next run depends on, so a
  refused run must not leave it agreeing with a degraded page.
- Both `--write` and `--update-ledger` are gated; `--update-ledger` alone is not a
  read-only invocation.
- Removal takes only `hmolar_formation`, plus the reference-state scaffolding if no other
  quantity is still using it.

## 6. C++ ingestion and API

Mirrors the existing `iGWP100` / `iODP` path exactly — that pattern is the reason this
tier is low-risk.

- `struct FormationStruct { double hmolar = _HUGE, hmolar_uncertainty = _HUGE;
  std::string source, version, id; }` in `include/CoolProp/CoolPropFluid.h`, beside
  `EnvironmentalFactorsStruct`; held as a member of `CoolPropFluid`.
- `parse_standard_state()` in `src/Backends/Helmholtz/Fluids/FluidLibrary.h`, guarded by
  `contains("STANDARD_STATE")`, mirroring `parse_environmental`.
  *(Amended post-review — it does NOT mirror `parse_environmental`'s failure behaviour.)*
  A malformed block **declines the value, never the fluid**: all six keys of
  `hmolar_formation` must be present and correctly typed, and `units` must read `J/mol`;
  any miss leaves `hmolar` unset so `HFORMATION` throws its own "no value available"
  message. Throwing instead would abort `add_one` for the whole fluid — EOS, ancillaries
  and transport all unavailable over an optional annotation — and strand the fluid's name
  in `fluids_list` (bd CoolProp-dwuu). Nothing requires these keys: there is no HEOS fluid
  schema at all, so a third-party fluid registered through `add_fluids_as_JSON` with a
  partial or mistyped block must keep loading. This is fail-closed on the number that
  matters: a `kJ/mol` block is never read as J/mol. Our own data cannot rot unnoticed
  because `[formation]` pins the readable count at exactly 76.
- `iHmolar_formation` in the `parameters` enum (`include/CoolProp/DataStructures.h`),
  registered as **trivial**, short name `HFORMATION`, units `J/mol`.
- Dispatch in `AbstractState::trivial_keyed_output` → virtual `calc_Hmolar_formation()`;
  base throws `NotImplementedError`, `HelmholtzEOSMixtureBackend` implements it.

**Unavailable → `ValueError` naming the fluid and stating that no ATcT value exists.**
Never 0, never NaN. A silently-zero formation enthalpy is indistinguishable from a
legitimate element value and would corrupt any reaction sum built on it.

**Pure and pseudo-pure fluids only.** A mixture input throws, matching `calc_GWP100` and
`calc_ODP`.

The mole-fraction sum would be *exactly* linear for the ideal-gas standard state, and
REFPROP does return it (verified: 0.7 CH₄ / 0.3 C₂H₆ gives exactly
`0.7·(−74591) + 0.3·(−83822)`). It is nonetheless excluded, for three reasons:

- A mixture is not a compound. Its "formation" includes an arbitrary mixing step, so the
  term is being stretched beyond what it conventionally denotes.
- The caller already holds the composition. A weighted mean adds no capability CoolProp
  is uniquely placed to provide, while inviting the reading that the number accounts for
  real-mixture behaviour it does not.
- **It would not generalize to the entropy tier.** S° of a mixture carries an entropy of
  mixing, −R Σ xᵢ ln xᵢ, and is therefore *not* linear. Supporting mixtures here would
  make the two standard-state quantities behave differently in a way no user could
  predict from the API.

Mixture support remains a clean additive change if a concrete need appears.

## 7. Testing

| test | tag | what it catches |
|---|---|---|
| Known values: methane −74513, water −241808, CO₂ −393477 J/mol | `[formation]` | ingestion / unit errors |
| `PropsSI("HFORMATION", ..., "R134a")` throws | `[formation]` | absent-value handling regressing to 0/NaN |
| A mixture input throws | `[formation]` | pure-only contract (§6) |
| Every stored value within ±2000 kJ/mol | `[formation]` | a kJ↔J slip anywhere in the pipeline |
| Regenerator unit tests: element row, `(g, cis)` row, ambiguous-CAS row | pytest | the three §3.2 traps |
| Four-value spot-check (methane, water, CO2, para-hydrogen) | pytest + `[formation]` | regenerator output vs. known ATcT values |

Three deliberate choices:

- **Para-hydrogen is spot-checked precisely because the CSV is wrong there** (§3.4).
  ATcT publishes −0.058 kJ/mol; the hand-curated dataset lists 0.000. Asserting our
  value pins the disagreement rather than letting it pass unnoticed.
- **No REFPROP cross-check test.** §3.1 establishes that REFPROP's `HFRM` is a derived
  quantity disagreeing with ATcT by up to 7σ, so any threshold either encodes REFPROP's
  ISO-6976 provenance as truth or is too loose to catch anything the spot-check
  does not already catch. The REFPROP comparison did its job as source-selection
  evidence; it is not a regression gate. `HEATFRMdll` is not called anywhere in this
  tier.
- **The full 75-value CSV reproduction was descoped.** Reproducing the contributor's
  dataset entry by entry would require committing that CSV as a test fixture, which
  reintroduces the hand-curated artifact this design exists to eliminate. The
  four-value spot-check plus `expected_coverage.json` give the same protection: the
  spot-check pins the values most likely to expose a units or binding error
  (including para-hydrogen, where the CSV is wrong), and the coverage ledger fails
  the run if any fluid changes state.

## 8. Out of scope

Each of these gets its own spec:

- **Standard molar entropy S°(298.15 K).** ATcT publishes only enthalpies — the site
  offers a single data product and its species pages are enthalpy-provenance only — so
  S° requires a different source ([Burcat / ReSpecTh](https://respecth.elte.hu/burcat.php)
  is the leading candidate, as it already adopts ATcT enthalpies), a different version
  stamp, and an explicit decision on the nuclear-spin / isotopic-mixing convention. That
  convention cancels in balanced reactions but **not** for ortho/para hydrogen, which is
  exactly where §4 binds spin-resolved enthalpies. §5's per-quantity provenance reserves
  the slot.
- **Heating values (HHV/LHV).** Derivable from Δ<sub>f</sub>H plus the parsed `FORMULA`
  field, restricted to CHONS fluids. See §8.1 for how ISO 6976 relates.
- **A `"FORMATION"` reference state.** The trap: the standard state is the *hypothetical
  ideal gas*, so the offset must be set from `hmolar() - hmolar_residual()`, never from a
  real state at (298.15 K, 1 bar) — which for n-decane or water lands in the liquid.
- **Reaction equilibrium / K_eq solving.** Belongs downstream (Cantera), per
  @matthewzyates in the thread.

### 8.1 Recorded decision — how ISO 6976 is handled

Decided 2026-07-29, so it is not re-litigated when heating values are taken up.

ISO 6976 is **two separable things**, and they get opposite treatment:

1. **A calculation procedure** — mixture combination rules, the compression-factor
   correction via summation factors on a volumetric basis, molar/mass/volume bases,
   Wobbe index, relative density. This is method, it is not encumbered, and it is worth
   implementing.
2. **A normative table of per-component calorific values** — edition-locked constants
   whose purpose is reproducible custody-transfer billing. This is a legal reference,
   not a thermophysical property.

**ISO 6976 data must not live in `dev/fluids/*.json`.** Its values form a matrix over
combustion reference temperature (25/20/15/0 °C) × metering reference conditions, so a
fluid record — which holds one canonical value per quantity — cannot express it without
flattening. Editions differ numerically and contracts cite a specific one, and a fluid
file has no versioning axis. And it covers only the ~60 natural-gas components: a table
in a metering standard that happens to be indexed by species.

**Therefore:** implement the procedure as a self-contained module parameterized on a
component table; ship a CoolProp-derived table (ATcT Δ<sub>f</sub>H + CoolProp's own
cp₀, which generalizes to any temperature) as the default, named or documented so it
does **not** imply ISO conformance; expose a documented hook for users to supply the
normative table. On the REFPROP backend, forward `HG`/`HN` (and the `VOL`/`LQ`
variants) — those already are ISO 6976 and REFPROP owns that provenance. A single
heating-value call must never return an ATcT number on one backend and an ISO number on
another without the caller knowing which.

Two supporting facts. REFPROP already uses the anchor-plus-cp₀ factoring — one stored
25 °C value per fluid, generalized with its own ideal-gas heat capacity (methane `hg` =
890.58 kJ/mol at 298.15 K, 871.30 kJ/mol at 500 K) — so this is a proven architecture.
And CoolProp is MIT-licensed while ISO standards are copyrighted, which makes vendoring
the normative table a licensing question rather than a technical one. The design above
does not depend on how that question resolves, which is its main virtue.

## 9. Acceptance criteria

1. `dev/atct/fetch_atct_formation.py --version 1.220` runs from a clean checkout and
   reproduces the committed `dev/fluids/*.json` byte-for-byte in the
   `STANDARD_STATE` block.
2. The script exits non-zero if any fluid's coverage state differs from
   `dev/atct/expected_coverage.json`.
3. The four spot-checked values (methane −74513, water −241808, CO₂ −393477,
   para-hydrogen −58 J/mol) are reproduced exactly, and `expected_coverage.json`
   matches the run. Full CSV reproduction was descoped — see §7.
4. `PropsSI("HFORMATION", "", 0, "", 0, "Methane")` returns −74513 J/mol; the same call
   for `R134a` throws with a message naming the fluid.
5. `./dev/ci/preflight.sh` passes.
