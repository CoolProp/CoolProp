# GERG-2004/GERG-2008 developer tooling

Four standalone developer scripts supporting the strict GERG-2004/GERG-2008
backend (`src/Backends/GERG/`). None is invoked by CMakeLists.txt or CTest
— CoolProp's build has no dependency on `teqp` and never will; these are
one-off, by-hand tools you re-run after touching the GERG tables or after a
`teqp` upgrade.

## `verify_transcription.py` — table-family cross-check (Tasks 1-6)

Mechanically re-verifies all five *transcribed* coefficient-table families
(pure-fluid residual coefficients, ideal-gas coefficients, binary
reducing parameters, departure functions/F_ij, pure-fluid reducing/critical
parameters) transcribed into `src/Backends/GERG/GERGData.h` /
`GERGBackend.cpp` against a local `teqp` checkout's
`include/teqp/models/GERG/GERG.hpp`. It parses both sources with regexes (no
build/link against teqp — that would pull Eigen/Boost/autodiff into a
lookup-table sanity check) and diffs every extracted value.

```bash
python3 dev/gerg/verify_transcription.py                      # ~/Code/teqp by default
python3 dev/gerg/verify_transcription.py --teqp-root /path/to/teqp
```

plus a sixth, `GERGAcentric.h`, which is **derived rather than transcribed**
and is therefore verified by recomputation — see `compute_acentric.py` below.

```bash
python3 dev/gerg/verify_transcription.py --skip-acentric        # five families only
```

Exits 0 and prints "OK" for each family on success; exits 1 with a
diff-shaped mismatch report otherwise. Run this after editing any table in
`GERGData.h`/`GERGBackend.cpp` — a typo fix, a new fluid, a value correction.

Family 6 is the one check that needs `teqp`'s **python module** at run time
rather than just a source checkout, because there is no acentric-factor table
in `GERG.hpp` to diff against. A missing `teqp` module is therefore a
**failure**, not a silent skip; `--skip-acentric` is the explicit opt-out and
changes the summary line so "OK" never quietly means "five of six".

The fifth family, `pure_info` (Table A3.5's tabulated `Tc_K`, `rhoc_molm3`,
`M_kgmol`, i.e. `detail::pure_info_2004()` / `detail::pure_info_2008_overrides()`
in `GERGData.h`), was added after the other four and checked in the same
style: exact row counts asserted (18 GERG-2004, 5 GERG-2008 overrides) so a
dropped row fails loudly, and every row compared both as a parsed float
(`vec_close`, `1e-13` relative) and as a normalised raw literal token —
the raw-token compare catches a transcription slip a tolerance-based float
compare can paper over. It compares the RAW table literals (mol/dm^3, K,
kg/kmol) as written in source, not `get_pure_info()`'s converted output
(mol/m^3, kg/mol): both sides apply the identical `*1000`/`/1000`
conversion, so comparing pre-conversion literals against pre-conversion
literals is the actual transcription check, while comparing post-conversion
output would just re-verify that arithmetic against itself. Verified to
discriminate: corrupting one digit of methane's tabulated `Tc_K` in
`GERGData.h` makes the script exit 1 and name `methane` in both the
float-close and raw-literal mismatch lines; restoring it returns to exit 0
with no diff.

## `generate_reference_values.py` — EOS-level reference values (Task 7)

Generates `src/Backends/GERG/GERGReferenceValues.h`, the fixture Tasks 8-9
validate CoolProp's *assembled* GERG equation of state against. This is a
different, stronger claim than `verify_transcription.py`: that script only
proves CoolProp's *tables* match teqp's tables. This one proves CoolProp's
*evaluated EOS* (residual/ideal Helmholtz derivatives, pressure, c_v, speed
of sound) reproduces teqp's numbers at specific (T, rho, z) points.

### Why the header is committed, not built from teqp

CoolProp's build must not depend on teqp (a research-code Eigen/autodiff/
Boost-heavy C++ library) just to run its own test suite. So the reference
values are generated once, on a developer machine with a teqp venv, and the
resulting `GERGReferenceValues.h` is committed like any other source file.
Regenerating it (same teqp version, same script) is expected to reproduce a
byte-identical file unless teqp's own GERG implementation changes — if a
regeneration attempt produces a diff, that is itself worth investigating
before assuming the new values are "more correct."

### Environment (teqp 0.23.2)

The system `python3`'s teqp (0.15.3 as of this writing) predates the GERG
factory registration and fails with `Unknown kind:GERG2008resid`. Use an
isolated venv:

```bash
python3 -m venv /tmp/gergenv
/tmp/gergenv/bin/pip install "teqp>=0.18" numpy
/tmp/gergenv/bin/python -c "import teqp; print(teqp.__version__)"   # sanity check
```

### Regenerating the header

```bash
/tmp/gergenv/bin/python dev/gerg/generate_reference_values.py \
    > src/Backends/GERG/GERGReferenceValues.h
uvx clang-format@18.1.8 -i src/Backends/GERG/GERGReferenceValues.h
```

Run both commands with a CWD inside this repo, writing directly to
`src/Backends/GERG/GERGReferenceValues.h` (or another path under the repo
tree) — not to a scratch directory elsewhere on disk. clang-format finds
its style by walking UP from the target file looking for `.clang-format`;
outside the repo tree it silently falls back to a different default style
(e.g. Allman brace-wrapping becomes attached braces) with no warning, which
looks like nondeterminism but is really "formatted against the wrong
config." Verified: regenerating+reformatting in place under the repo
reproduces the committed header byte-for-byte; doing the same in `/tmp`
does not.

The clang-format pass is required, not cosmetic: the pre-PR gate runs
`clang-format --dry-run -Werror` against every changed `.h` file, and the
generator's own 4-space, one-point-per-line emission does not satisfy it.
Applying clang-format re-indents everything to 2 spaces and, for the AGA8
mixture rows (21-element `names`/`z` vectors per row), explodes each row
across roughly 20 lines — one element per line — to respect the column
limit. **This looks bad** (a single logical mixture-composition record is
no longer one line you can eyeball), but every value is complete and
unchanged, the file still compiles, and there is no `// clang-format off`
precedent anywhere in this repository to reach for instead. This was
flagged as an expected risk before generation, and it is exactly what
happened: the honest trade-off is "correct but visually noisy" over adding
a suppression convention that doesn't otherwise exist in this codebase.

The generator prints a coverage/cross-check summary to stderr (also
embedded verbatim, as `//`-prefixed lines, near the top of the generated
header) — check it after every regeneration:

```
pure_points_2004: 288 emitted (18 fluids x 16 grid pts, floor 14/16 finite alphar enforced per fluid), NaN'd-field counts: alphar=0, alphaig=0, p=0, cv=0, w=23
pure_points_2008: 336 emitted (21 fluids x 16 grid pts, floor 14/16 finite alphar enforced per fluid), NaN'd-field counts: alphar=0, alphaig=0, p=0, cv=0, w=27
mix_points_2004 (binary pairs only): 153 of 153 expected emitted (assert enforced), 19 with w = NaN
mix_points_2008 (binary pairs + AGA8): 397 emitted (210 of 210 expected pairs + 187 of 187 expected AGA8 gases, both counts assert-enforced), 40 pairs with w = NaN, 0 AGA8 with w = NaN
Mixture alphar/alphaig: 0 NaN on every emitted row of both vectors (assert-enforced). w is the ONLY mixture column that is ever NaN, so Task 9's per-field skip may fire on w and on nothing else.
AGA8 p_teqp vs validation_data.P_MPa: worst 6.770e-04 (gas 193), median 3.514e-12 ...
```

### Coverage is enforced, not just reported

`main()` asserts, and aborts with a descriptive message rather than
emitting a short header, if any of the following don't hold exactly:

- `len(mix_points_2004)` (binary pairs) `== 153` (`C(18,2)`)
- binary pairs in `mix_points_2008 == 210` (`C(21,2)`)
- AGA8 rows emitted `== len(VALIDATION_DATA) == 187`
- **every** emitted mixture row has a finite `alphar` *and* a finite
  `alphaig` (`nanalpha_*` / `aga8_nanalpha` must all be empty)

Why this matters more here than almost anywhere else in this backend: the
Task 5 review established that a single mistyped digit in one of the 225
binary reducing-parameter rows (153 GERG-2004 + 72 GERG-2008-override) is
guarded by **nothing** except these binary-pair reference points existing
and being numerically compared in Task 8/9. Before this assertion existed,
a future regeneration (a teqp version bump, a coefficient edit) that made
one pair evaluate to a caught exception instead of a finite-but-wrong
number would make that pair silently vanish from the header — removing the
only guard on that reducing-parameter row while every test stayed green.
Now that regeneration aborts loudly instead.

The pure-fluid grid gets a **per-fluid** floor instead of a global count:
each fluid must retain a finite `alphar` (teqp's `get_Ar00`) at least
`MIN_FINITE_ALPHAR_PER_FLUID = 14` of its 16 grid points. This is
per-fluid, not "at least N out of 288/336 total," specifically so one
fluid failing broadly cannot hide behind other fluids that have full
coverage — a global total can stay comfortably above threshold even if one
fluid degrades to near-zero usable points. In practice every fluid
currently has all 16/16 finite (see the next section for why alphar itself
essentially never goes non-finite); the floor is a tripwire for future
regressions, not a currently-binding constraint.

All three assertion classes were verified to actually fire: temporarily
forcing one binary pair, one AGA8 gas, and one fluid's `alphar` to drop out
(each on a throwaway copy of this script, never committed) independently
raised `AssertionError`, e.g.:

```
AssertionError: mix_points_2004 has 152 binary-pair rows, expected exactly 153 -- 1 pair(s)
dropped: [('methane', 'nitrogen')]. A missing pair removes the ONLY numeric guard on that
reducing-parameter row; investigate, do not ship.

AssertionError: AGA8 rows: 186 emitted from 187 validation_data entries, expected 187 of
both -- 1 gas(es) dropped: [2]

AssertionError: 2004/methane: only 0/16 grid points have a finite alphar -- this fluid is
failing broadly (not just isolated near-phase-boundary points); investigate before
shipping, do not silently ship a near-empty fluid.
```

### `validation_data.P_MPa` cross-check: why the worst-case row is 6.8e-4, not 1e-12

`_validation_data.py`/`_mixture_comps.py` are byte-for-byte transcriptions
of teqp's own `validation_data`/`mixture_comps` tables
(`~/Code/teqp/src/tests/catch_test_GERG.cxx:119-321,322-...`), verified by
diffing every row against the C++ source at generation time. All 187 AGA8
rows are generated with teqp directly (full `%.17g` double precision); the
table's own `P_MPa` column is used only as an independent cross-check, not
as the fixture.

When the mixture is built with the *correct* AGA8 component ordering (see
below), the median relative disagreement between teqp's computed pressure
and `validation_data.P_MPa` across all 187 rows is 3.5e-12 — machine-noise
level, confirming the ordering/composition wiring is right. But the *worst*
row (gas 193) disagrees by 6.8e-4 (41/187 rows exceed 1e-6). This is not a
composition/ordering bug:

- Deliberately rebuilding the model with the *wrong* ordering
  (`GERG2008::component_names` order instead of the AGA8 test-code order)
  changes the picture completely: median jumps to 7.1e-4, the worst row
  to 0.27 (27%), and only 7/187 rows stay below 1e-9 (down from 104/187).
  A real ordering bug produces that kind of uniform, everywhere-large
  disagreement — not a smooth distribution where most rows are at
  machine precision and a handful of unusual mixtures are worse.
- The worst-offending rows are consistently rich, atypical natural-gas
  blends far from GERG-2008's typical calibration composition — e.g. gas
  193 is 64.5% CO2, 21% methane, 7.7% N2, 5.1% H2S, 0.6% helium. GERG-2008
  is documented to have larger uncertainty for such acid-gas/CO2-dominant,
  multi-component blends than for typical pipeline-quality natural gas.
- teqp's own test (`catch_test_GERG.cxx:696-698`) computes this exact
  comparison and then *comments out* the `CHECK_THAT` against
  `validation_data[i].P_MPa` entirely, at a 1e-5 tolerance — i.e. teqp's own
  authors already know this comparison doesn't reliably pass and disabled
  it, rather than it being an oversight in this script.
- Independent confirmation from the Task 7 review: NIST's own bundled AGA8
  reference implementation, re-run directly (not through teqp) against the
  same gas-193 state, *also* disagrees with the published table value.
  That upgrades the conclusion from "teqp disagrees with the table" to
  "the reference C/FORTRAN code disagrees with its own published table" —
  i.e. this is a table-vs-reference-code discrepancy in the published AGA8
  validation data itself, not anything specific to teqp or to this
  generator.

Conclusion: `validation_data.P_MPa` is real published AGA8 measurement/
reference data that neither GERG-2008 (via teqp) nor NIST's own bundled
AGA8 reference code reproduces exactly for its most unusual compositions —
a property of the published reference table, not a transcription or
ordering bug anywhere in this pipeline. The reference values Tasks 8/9
validate against come from teqp directly, not from this column, so this
discrepancy does not propagate into the fixture.

### AGA8 component ordering (load-bearing, do not "simplify")

The columns of `mixture_comps` are **not** in `GERG2008::component_names`
order. They follow the AGA8 test-code ordering declared separately as
`components` (`catch_test_GERG.cxx:512`):

```
methane, nitrogen, carbondioxide, ethane, propane, isobutane, n-butane,
isopentane, n-pentane, n-hexane, n-heptane, n-octane, n-nonane, n-decane,
hydrogen, oxygen, carbonmonoxide, water, hydrogensulfide, helium, argon
```

Note isobutane precedes n-butane, isopentane precedes n-pentane, and
helium/argon are last — all three differ from `component_names` order.
Every `MixRefPoint` emitted by this generator carries its component
**names** alongside the mole fractions (not a bare `std::vector<double>`)
specifically so this ordering is self-describing in the fixture itself and
Task 9 cannot reintroduce this bug by re-deriving "which index is which
component" from a convention that lives only in a comment.

### Derivative-accessor naming (load-bearing)

teqp's Python accessors are `get_Ar<TAU order><DELTA order>` — e.g.
`get_Ar20` is the tau-tau second derivative, `get_Ar02` is delta-delta. An
earlier draft of this generator had these reversed, which still produced a
finite, plausible-looking number. Verified empirically: on the ideal-gas
model, `get_Ar20` gives `-cv_ideal/R` (3.303 for methane at 300 K, correct),
while `get_Ar02` gives exactly `-1.0` (a tell, not a real c_v). The p/c_v/w
expressions in `point()` mirror teqp's own test
(`catch_test_GERG.cxx:675-694`) exactly rather than being re-derived; any
future disagreement between this script and that test is a bug in this
script, not in teqp.

### Struct shapes consumed by Tasks 8/9

```cpp
namespace CoolProp { namespace GERG { namespace reference {

struct PureRefPoint {
    const char* name;                                       // GERG component name
    double T_K, rhomolar, alphar, alphaig, p_Pa, cvmolar, w; // ANY field may independently be NaN
};

struct MixRefPoint {
    std::vector<const char*> names;  // parallel to z -- component order is AGA8 order for
    std::vector<double> z;           // AGA8 rows, GERGData.h component_names() order for pair rows
    double T_K, rhomolar, alphar, alphaig, p_Pa, cvmolar, w;
    // w may independently be NaN (see below).  alphar/alphaig/p_Pa/cvmolar are never NaN --
    // for alphar/alphaig that is now assert-enforced at generation time, because those two
    // columns are the ones Task 9's gate depends on and a NaN in them would be silently
    // SKIPPED by the per-field comparison rather than failing.  (See the pure-fluid note
    // below for why the same guarantee is not claimed for PureRefPoint in general.)
};

extern const std::vector<PureRefPoint> pure_points_2004, pure_points_2008;
extern const std::vector<MixRefPoint> mix_points_2004, mix_points_2008;

}}}  // namespace CoolProp::GERG::reference
```

**Why `MixRefPoint` carries `alphar`/`alphaig` (added in Task 9).** The
struct originally had only `p_Pa`, `cvmolar` and `w`. Task 8's mutation
ledger established that `alphaig` is the *only* instrument in this suite
that can see an error in the ideal-gas integration constants below the
`h = s = 0` reference-state test's `1e-8` absolute floor: a `1e-9` additive
shift in `n0[1]` fails **624 `alphaig` assertions and nothing else** —
`alphar`, `p_Pa`, `cvmolar` and `w` all stay clean. Task 9's central hazard
(the mixture branch of `calc_alpha0_deriv_nocache` calling `set_Tred(Tr)`
while passing `tau_i = Tc_i/T`) is an ideal-gas error that *only* mixtures
expose, so without a mixture `alphaig` column the Task 9 gate could not see
the very bug it was written to catch. Re-running that same `n0[1] += 1e-9`
mutation with the column in place now fails **1174** assertions, every one
of them `alphaig`: 624 pure + 550 mixture. `alphar` was added alongside it
for symmetry with `PureRefPoint` and because it isolates the departure
function from the gas constant — a wrong `R` moves `p`, `cvmolar` and `w`
but leaves `alphar` exact.

**Per-FIELD nulling, not per-row dropping (load-bearing contract for Tasks
8-9).** Every field of both structs is computed and null'd to NaN
*independently*. No row is ever dropped because one field came out
non-finite. Concretely:

- `pure_points_2004`/`pure_points_2008` contain **exactly** 16 rows per
  fluid (4 T-factors x 4 rho-factors), always — nothing is ever omitted
  from the grid. In the currently-committed header, `alphar`, `alphaig`,
  `p_Pa`, and `cvmolar` are finite on **every single one** of the 624 pure
  rows; only `w` is ever NaN (23/288 rows in `pure_points_2004`, 27/336 in
  `pure_points_2008`), always at `T = 0.7 x Tc` combined with `rho = 0.5 x
  rhoc` or `1.5 x rhoc` — i.e. exactly where the 0.7-Tc isotherm is likely
  to cross the two-phase dome, putting that state on the
  mechanically-unstable branch where `w^2 < 0`, the same phenomenon
  described for mixtures below. An earlier version of this generator
  dropped the WHOLE row whenever ANY of the 5 fields was non-finite,
  which — since only `w` was ever actually the culprit — silently threw
  away ~50 perfectly good `p_Pa`/`cvmolar` points clustered in exactly the
  near-phase-boundary region where a wiring bug is most likely to show up.
  **Tasks 8/9 must check `std::isnan(...)` per field being compared, not
  skip the whole row when any single field is NaN.**
- `w == NaN` (`std::isnan(w)`) happens for 19/153 GERG-2004 and 40/210
  GERG-2008 binary-pair rows at the fixed (T=250K, rho=5000 mol/m^3,
  z=[0.4, 0.6]) state — that state sits on the mechanically-unstable
  (spinodal) branch of the single-phase EOS for some dissimilar-component
  pairs, where `w^2 < 0`. `p_Pa`/`cvmolar` remain finite and meaningful
  there; Task 9 must compare them unconditionally and skip only the `w`
  comparison when `std::isnan(w)` is true for that row.


## `fit_ancillaries.py` — saturation ancillaries (Task 11)

Traces the saturation curve of each of the **23 distinct pure GERG EOS** with
`teqp`, regresses CoolProp's standard ancillary forms against it, and emits
`src/Backends/GERG/GERGAncillaries.h`.

23, not 21: the GERG-2008 component list has 21 fluids, and GERG-2008 changed
the reducing parameters of **carbon monoxide** and **isopentane** relative to
GERG-2004 (`GERGData.h`, `pure_info_2008_overrides`). Changed reducing
parameters move the whole saturation curve, so those two need one fit per
model. The generated table is therefore a GERG-2004 base of 18 plus a
GERG-2008 override set of 5 — the same base-plus-overrides shape
`get_pure_info` uses.

### These are fitted against the GERG pure EOS, and MUST be refitted if a coefficient table changes

CoolProp already ships saturation ancillaries for every fluid GERG names, and
they would work perfectly well as VLE seeds — which is precisely the trap.
They belong to each fluid's **reference** equation of state (Setzmann-Wagner
for methane, Span-Wagner for carbon dioxide, IAPWS-95 for water), not to
GERG's shortened technical form. A backend named GERG-2004/GERG-2008 must not
carry data traceable to a different equation, with nothing anywhere to say so.

Consequently: **editing any pure-fluid `n`/`t`/`d`/`l` row, or any Table A3.5
reducing parameter, invalidates the corresponding ancillary.** Re-run this
script. Nothing in the build detects the staleness for you; what would catch it
is the `[GERG]` test *GERG ancillaries are close to the converged saturation
state*, which compares each fit against a freshly converged GERG saturation
state at five temperatures per fluid using that fit's own recorded
`max_rel_dev` as the tolerance.

### An ancillary is a seed, not an answer

`FlashRoutines::QT_flash` hands `rhoL`/`rhoV`/`pV` to
`SaturationSolvers::saturation_T_pure_Maxwell`, which iterates to the true GERG
saturation state; no ancillary value is ever returned to a caller. The
`[GERG]` test *GERG saturation is independent of the ancillary seed* pins this
by perturbing all three ancillaries by ±1% (about three times the worst fit
deviation in the shipped table) and requiring the converged answer not to move
beyond 1e-9 relative.

This is also why **no superancillary may ever be attached to a GERG fluid**.
`FlashRoutines::sat_superanc_path_applies` (`FlashRoutines.cpp:558`) routes
pure-fluid saturation straight to the Chebyshev expansion and *returns that as
the answer*, so a GERG fluid carrying CoolProp's blob would silently report
Setzmann-Wagner saturation densities labelled GERG-2008. `make_gerg_fluid`
never calls `set_superancillaries_str`; the test *GERG carries no
superancillary but does carry the saturation end states* pins it.

### Schema, and why every fit carries a `t = 0` term

The emitted rows use the same schema as `dev/fluids/*.json`'s `ANCILLARIES`
block, so `make_gerg_fluid` builds a plain CoolProp `SaturationAncillaryFunction`
from each one (through `SaturationAncillaryFunction::Values`, the same bundle
`cpjson::make_saturation_ancillary` fills from JSON) and the shipped evaluator
in `Ancillaries.cpp` is the only evaluator:

| slot   | `type`        | `using_tau_r` | form                                            |
|--------|---------------|---------------|-------------------------------------------------|
| `rhoL` | `rhoLnoexp`   | false         | `reducing_value * (1 + sum n_i theta^t_i)`       |
| `rhoV` | `rhoV`        | true          | `reducing_value * exp(T_r/T * sum n_i theta^t_i)`|
| `pV`   | `pV`          | true          | `reducing_value * exp(T_r/T * sum n_i theta^t_i)`|

with `theta = 1 - T/T_r`.

For most CoolProp fluids the ancillary reducing point **is** the EOS critical
point, so `rhoL(Tc) = rhoV(Tc) = rhoc` and `p(Tc) = pc` fall out of the form
for free (`theta = 0` kills every `t > 0` term). **That is not true for GERG.**
Table A3.5's reducing parameters are fitted quantities, and the true critical
point of the shortened form — located with `teqp`'s own critical solver — can
sit a long way from them:

| fluid       | `Tc(true) - T_reduce` | `rhoc(true)/rho_reduce - 1` |
|-------------|-----------------------|------------------------------|
| n-heptane   | +1.096 K              | -3.06%                       |
| n-butane    | +0.634 K              | -5.43%                       |
| isopentane (2004) | +0.381 K        | +3.35%                       |
| n-octane    | +0.250 K              | -3.10%                       |
| propane     | +0.114 K              | -3.25%                       |
| oxygen      | +0.114 K              | -3.73%                       |
| argon       | +0.104 K              | -2.81%                       |
| isobutane   | -0.067 K              | -3.16%                       |

So a `theta^0` term is included in every fit, leaving the value at `T = T_r`
a free fitted parameter instead of pinning it to the tabulated reducing value.

`T_r = max(Tc(true), T_reduce)`, and the maximum is load-bearing both ways:

- It can never be **below** the tabulated reducing temperature, because
  `SaturationAncillaryFunction::evaluate` returns NaN for `T > T_r` (guarding
  `pow(negative, fractional)`, GitHub #1611) while `QT_flash`'s upper
  saturation bound is the backend's `T_critical()`, i.e. the tabulated value.
- Where `Tc(true)` is the **higher** of the pair, using it puts `theta = 0`
  exactly at the end of the saturation curve, which is what makes the
  classical `theta^(1/2)` scaling representable. Pinning `theta = 0` 1.1 K
  short of the critical point instead leaves a square-root branch point inside
  the fit domain and costs an order of magnitude: propane fitted 7.3e-3 that
  way versus 1.4e-3 with `T_r = Tc(true)`.

`reducing_value` stays the **tabulated** reducing density, and for `pV` the
pressure evaluated from the EOS at the tabulated reducing state — the same
number `make_gerg_fluid` stores in `EOS.reduce.p` and the backend reports as
`p_critical()`.

### Method

1. **Trace.** From `theta = 1e-4` below the true critical temperature down to
   `max(0.30*Tc, T where p_sat < 1e-3 Pa)`, 250 points geometric in `theta`.
   Each step is seeded from the previous one (the technique
   `~/Code/fastchebpure`'s `fastcheb.cpp` uses), with a guess-free fallback
   that cannot be stranded by a bad seed: scan for the two spinodal densities,
   bisect on `p` between them with equal chemical potential as the residual,
   then polish with a damped 2-D Newton on `(rhoL, rhoV)`.
   The `1e-3 Pa` floor and the `0.30*Tc` floor are **numerical**, not
   statements about where the GERG EOS stops being valid.
2. **Regress.** Once the exponents are fixed each ancillary is linear in `n`,
   so the exponents come from greedy forward selection over a pre-spaced
   candidate ladder and `n` from a weighted linear least squares. Weights make
   the objective *relative* error in `rho`/`p`. A minimum-separation rule
   (`MIN_EXPONENT_GAP`/`MIN_EXPONENT_RATIO`) refuses two nearly-equal
   exponents: they span nearly the same function, so least squares can only
   distinguish them with huge opposite-signed coefficients — accurate on paper,
   but evaluated as the difference of two ~1e3 terms whose sum is ~1. An
   earlier unspaced ladder produced exactly that (coefficients to 2606 for
   methane's `rhoL`).
3. **Gate.** The worst relative deviation per fluid is printed and embedded in
   each row's `max_rel_dev`. `MAX_ACCEPTABLE_DEV = 0.005` **aborts** the script
   rather than shipping a bad fit; the `[GERG]` test independently rejects any
   row claiming worse than 1%.

### Running it

```bash
# reuse the existing venv from generate_reference_values.py (teqp 0.23.2)
/tmp/gergenv/bin/python dev/gerg/fit_ancillaries.py --report      # deviations per fluid
/tmp/gergenv/bin/python dev/gerg/fit_ancillaries.py \
    > src/Backends/GERG/GERGAncillaries.h
uvx clang-format@18.1.8 -i src/Backends/GERG/GERGAncillaries.h
```

`--report` also prints the reducing-versus-true-critical columns tabulated
above. Worst deviation across all 23 pure EOS as shipped: **3.822e-3**
(`rhoL`, n-hexane); worst `rhoV` 2.730e-3 (n-heptane), worst `pV` 1.180e-3
(n-decane).

### Saturation end states, and what `triple_liquid` means here

Each row also carries the saturation state at the low-temperature end of the
fitted range, traced rather than guessed. `make_gerg_fluid` writes it into
`EOS.sat_min_liquid`/`EOS.sat_min_vapor` — what `calc_Tmin_sat`/`calc_pmin_sat`
return, and hence what bounds `QT_flash` from below — and into
`CoolPropFluid::triple_liquid`/`triple_vapor`.

**`triple_*` is CoolProp's field name, not a claim about GERG.** GERG publishes
no triple point and `EOS.Ttriple` stays 0. `saturation_T_pure_Maxwell`
(`VLERoutines.cpp:965-980`) reads those two slots purely as "the
low-temperature end of the saturation curve", to sanity-band the ancillary seed
and to build a linear fallback seed. Left at `SimpleState`'s `_HUGE` default,
every seed would fail that band and every pure-fluid saturation call would take
the fallback path.

### Two consequences worth knowing about

- **Helium and hydrogen have no reachable saturation states.** `make_gerg_fluid`
  clamps `EOS.limits.Tmin` to `min(60 K, Tc)`, and for those two that *is* `Tc`
  (5.1953 K and 33.19 K), so the whole subcritical region is outside the
  backend's range. Ancillaries are still fitted and shipped for them — the
  table is complete, and `DONT_CHECK_PROPERTY_LIMITS` reaches them — but the
  `[GERG]` saturation sweep skips them, and asserts that it skips exactly 4
  fluid/model combinations so the skip cannot quietly widen.
- **The authoritative range is the mixture-model range.** GERG-2008's published
  envelope, 60-700 K and p <= 70 MPa (Kunz & Wagner 2012 §4.1), is stated for
  the mixture model as a whole. GERG publishes no per-component lower
  temperature limit, and this backend deliberately does not consult CoolProp's
  triple-point data. A pure-component `Tmin` below 60 K (helium's 5.1953 K) is
  the removal of a self-contradiction — `Tmin` above `Tc` — not a validity
  statement, and neither is an ancillary `Tmin` below 60 K (methane's is
  57.17 K).

## `compute_acentric.py` — acentric factors (bd CoolProp-deut)

Generates `src/Backends/GERG/GERGAcentric.h`, the 23-row table
`make_gerg_fluid` writes into `EquationOfState::acentric`.

### Why it exists

GERG publishes no acentric factor, so `make_gerg_fluid` used to leave
`EOS.acentric` at the `_HUGE` sentinel. That sentinel is read **directly**,
never through the overridable `calc_acentric_factor` accessor, by three
consumers that none of them range-check it:

- the Wilson K-factor seed in `SaturationSolvers` — every mixture QT/PQ flash
  and `build_phase_envelope()`;
- `FlashRoutines::T_DP_PengRobinson` — every gas-phase `DmolarP` flash;
- `HelmholtzEOSMixtureBackend::solver_rho_Tp_SRK`, which `solver_rho_Tp` seeds
  from.

Borrowing CoolProp's tabulated value for the same fluid would work numerically
and would break the property the whole backend exists to have. So it is
**derived from GERG's own equation** by Pitzer's defining relation:

```
omega = -1 - log10( p_sat(0.7 * T_c) / p_c )
```

23 rows for 21 components, because GERG-2008 moves the reducing parameters of
carbon monoxide and isopentane, which moves their saturation curve.

### The two judgement calls

1. **`p_sat(0.7 T_c)` is a CONVERGED solve, not the fitted `pV` ancillary.**
   An ancillary is a seed, never an answer; a shipped constant that depended on
   one would move whenever the ancillary was refit with a different term
   ladder. Reliability was measured, not assumed: the guess-free solve
   converges at `0.7*T_c` for **all 23** with a worst `|p'-p''|/p''` of
   5.5e-12, so there is no fallback path and none is implemented — a silent
   ancillary fallback is exactly the fail-open this backend refuses elsewhere.
   Measured size of the choice: the ancillary would move `omega` by at most
   9.1e-05 (nitrogen).
2. **`p_c` is the pressure at GERG's TABULATED reducing state**, i.e. the same
   number `make_gerg_fluid` stores in `EOS.reduce.p` and the backend reports
   from `p_critical()` — not the pressure at the true critical point of the
   shortened form. The reason is consistency with the consumer: Wilson's
   `ln K_i = ln(p_c,i/p) + 5.373 (1+omega_i)(1 - T_c,i/T)` reads `p_c,i` and
   `T_c,i` off the same fluid object that carries `omega_i`, so mixing a
   true-critical `p_c` into `omega` would put the inconsistency into the seed.
   `0.7*T_c` uses the tabulated `T_c` for the same reason. Measured size of
   the choice: up to 6.7e-03 in `omega` (n-heptane).

### Helium and hydrogen

Their enforced `Tmin` equals their reducing temperature (5.1953 K, 33.19 K), so
`0.7*T_c` — 3.6367 K and 23.2330 K — is below the range
`check_gerg_range_of_validity` lets a caller reach.

Computed there anyway, deliberately. `omega` is a definitional constant of the
equation, not a state a caller asks for: it is evaluated offline, once, and
only the scalar ships. The equations are perfectly well behaved there (the
`fit_ancillaries.py` trace reaches `0.30*T_c` for both). No run-time guard was
weakened, and the `[GERG]` self-consistency test reaches those two states
through CoolProp's own documented `DONT_CHECK_PROPERTY_LIMITS` opt-out rather
than around the guard.

### Sanity, against published values

Not a validation target — GERG's reducing point is a fitted quantity and these
come from the shortened technical form — but a useful net for gross errors.
Worst disagreement across all 23 rows: **3.1e-3** (n-octane). Methane
0.011385 (published 0.011), propane 0.152927 (0.152), n-butane 0.199237
(0.200), water 0.344969 (0.344), helium -0.385939 (-0.385), hydrogen
-0.218652 (-0.219).

### Running it

```bash
# reuse the existing venv from generate_reference_values.py (teqp 0.23.2)
/tmp/gergenv/bin/python dev/gerg/compute_acentric.py --report   # full table + comparisons
/tmp/gergenv/bin/python dev/gerg/compute_acentric.py \
    > src/Backends/GERG/GERGAcentric.h
uvx clang-format@18.1.8 -i src/Backends/GERG/GERGAcentric.h
```

Unlike `fit_ancillaries.py` (a least-squares solve, whose last two or three
digits move with the numpy version), this is a bisection plus a Newton polish
to a fixed residual and reproduces to ~1e-12 across environments;
`verify_transcription.py` compares with a 1e-9 absolute tolerance.

### Two verification layers, doing different jobs

- `verify_transcription.py` family 6 re-derives the table with **this
  script's** machinery and diffs it against the committed literals. That
  catches a stale table — someone edits a coefficient row and forgets to
  rerun — but reuses the same solver, so it is not an independent confirmation
  of the number.
- The Catch2 case *"GERG acentric factors satisfy their own definition"*
  re-derives all 39 (model, component) pairs through **CoolProp's own**
  saturation solver (`saturation_T_pure_Maxwell`, ancillary-seeded) and checks
  it reproduces the committed value. Nothing but the underlying equation is
  shared with teqp's path. Worst observed disagreement: 8.6e-12 in `omega`.

`GERG::get_acentric_factor` additionally rejects a non-finite table row on the
way out, through a **single exit point** — both the GERG-2004 base lookup and
the GERG-2008 override lookup pass through the same check. That structure is
load-bearing, not stylistic: an earlier version returned straight out of the
override branch, which left the guard inert for exactly the five rows in that
table. Injection-verified in both directions: an infinity in `n-nonane`
(override table) and in `water` (base table) both now throw
`The GERG acentric factor for [...] is not a finite number`; with the early
return restored, `n-nonane` built a fluid and answered `inf` while `water`
threw. The generator also refuses to emit a non-finite omega, but that does not
substitute for the guard — the generator cannot police a hand-edit made after
it ran, which is the case the guard exists for.

Mutation-verified: perturbing methane's `omega` from 0.011385 to 0.031385
fails both of those layers, and — the point of the exercise — moves the
*converged* mixture VLE answers by at most 4.2e-12 relative (bubble and dew
pressure for 90 % methane / 10 % ethane at 150 K, and the PQ round-trip
temperature), confirming that `omega` reaches those results only through the
initial guess.
