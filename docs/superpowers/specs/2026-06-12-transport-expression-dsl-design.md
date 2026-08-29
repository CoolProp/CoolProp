# Transport Property Expression DSL — Design

**Date:** 2026-06-12
**Status:** Implemented (PR #3185); this document is the reference for the
shipped behaviour, kept current through the review round that followed.
**Branch:** `worktree-exprtk`

## Problem

Adding a new transport-property (viscosity / thermal conductivity) correlation
form to CoolProp today requires editing the C++ core in four places: a new
`std::vector`-backed POD struct + enum value in `CoolPropFluid.h`, a new parser
branch in `FluidLibrary.h`, a new `case` in the dispatch switch in
`HelmholtzEOSMixtureBackend.cpp`, and a new static routine in
`TransportRoutines.cpp`. Every new functional form means a recompile of the
library core. This is the friction that makes implementing transport routines a
chore.

We want correlation forms expressible as **runtime-loaded formula strings** in
the fluid JSON, evaluated by an embedded expression interpreter, so that adding a
new *parametric* form is data, not code.

A prior NIST effort (Lemmon) used a reverse-polish-notation encoding for the same
goal; it worked but the RPN grammar is effectively unreadable. The explicit
design goal here is a small, **readable, infix, Python-flavored** grammar whose
formulas resemble the LaTeX/sigma notation in the source literature.

## Scope

### In scope (Tier A — pure parametric forms)

The DSL must be expressive enough to reproduce **every Tier-A correlation** — the
forms that are pure functions of `(T, ρ, τ, δ)` plus coefficient arrays:

- Viscosity: `dilute_powers_of_T`, `dilute_powers_of_Tr`,
  `dilute_collision_integral`, `dilute_collision_integral_powers_of_T`,
  `dilute_kinetic_theory`, `initial_density_dependence_Rainwater_Friend`,
  `initial_density_dependence_empirical`,
  `higher_order_modified_Batschinski_Hildebrand`.
- Conductivity: `dilute_ratio_polynomials`, `dilute_eta0_and_poly`,
  `residual_polynomial`, `residual_polynomial_and_exponential`.

Reproducing all of these to ULP is the **completeness proof** for the scope.

**Amendment (initial-density stage wired, PR follow-up to #3185).** Two of the
listed forms have no golden test and, as the stage is now wired, one of them
cannot have a like-for-like one:

- `initial_density_dependence_Rainwater_Friend` — the host arm scales the routine's
  return by `eta_dilute * rho`, whereas an `initial_density` expression block yields
  the stage's contribution **directly** (see §3a). So an expression block is not a
  drop-in replacement for an RF entry: to reproduce one it must compute
  `eta_0 * B_eta * rho` itself, recomputing the dilute term inside the block. That
  is what the shipped ethene and propylene-glycol correlations do, and it works —
  but it is a re-expression rather than a substitution, so "reproduces the RF
  routine to ULP" is not the right test and none is claimed.
- `dilute_kinetic_theory` — simply untested; no shipped fluid exercises it through
  the DSL yet.

The completeness claim therefore covers the eight golden-tested Tier-A forms plus
`initial_density_dependence_empirical`, not the full list above.

### Out of scope (v1 non-goals)

- **Tier B** — closed-form but requiring EOS-derived scalars:
  `simplified_Olchowy_Sengers` (critical enhancement), `friction_theory`,
  `Chung`. The derived-variable registry (§2b) is the designed extension path,
  but no Tier-B form is implemented in v1.
- **Tier C** — algorithmic / reference-fluid procedures: `viscosity_ECS`,
  `conductivity_ECS`. These stay hardcoded.
- **Fluid-specific `*_hardcoded` routines** (water, heavy water, helium, R23,
  xylenes, toluene, ethane, hydrogen, benzene, hexane, heptane, CO2 Laesecke,
  methanol, ammonia, methane …) — stay hardcoded, by explicit decision.
- No migration of existing fluids' JSON: the DSL is **opt-in per transport
  block**. Existing `type` values and their C++ routines are untouched.
- No autodiff, no mixture-specific transport handling, no bytecode VM.

## Decisions (locked during brainstorming)

| Axis | Decision |
|---|---|
| Boundary | Tier A only; Tier B/C and `*_hardcoded` stay in C++ |
| Integration | Additive — new `"type": "expression"` in existing transport sub-blocks |
| Implementation | Hand-rolled, zero third-party dependencies |
| Evaluator | Tree-walking AST, compiled once per fluid at load, evaluated many |
| Sum syntax | `sum(i: <body>)` explicit index; arrays subscripted `arr[i]` |
| Variables | Raw state + JSON-declared constants/arrays + `let` bindings; **all quantities base SI, always** |
| Thermodynamic inputs | ONE bucket, keyed by the existing `CoolProp::parameters` enum; the host fills every requested key with `AbstractState::keyed_output()`. Resolved per-block from a declared `state_variables` list via `is_valid_parameter()` — CoolProp's own names, no list kept here — minus a policy denylist (transport outputs; critical/reducing state) |

### Why tree-walking AST, not bytecode VM

A bytecode VM would shave ~15–20% off dispatch, but Tier-A evaluation is
dominated by `pow`/`exp`/`log`, so tree-walk overhead is negligible while the VM
doubles the code to own and test. Bytecode remains a localized later optimization
if an expression ever lands in a hot loop. YAGNI for v1.

## Architecture

Five small, independently testable units plus a thin host-integration layer.
The **evaluator is a pure function of its context** (scalar values + array spans)
and never touches `HelmholtzEOSMixtureBackend` — only the host layer does. This
preserves the unit boundary and lets golden tests inject values directly without
standing up an EOS.

```
formula string ──Lexer──▶ tokens ──Parser──▶ AST ──Binder──▶ compiled Program
                                                                    │
                                          (per eval) context ──Evaluator──▶ double
```

### 1. The language (grammar)

- **Program:** zero or more `let <name> = <expr>` bindings (newline- or
  `;`-separated), followed by exactly one final result expression. `let` is pure
  single-assignment name-binding — no reassignment, no control flow.
- **Operators:** `+ - * /`, `^` (power, **right-associative**), unary `-`.
  Standard precedence: `^` > unary `-` > `* /` > `+ -`.
- **Functions:** Tier A needs `exp, ln, log10, sqrt, abs, pow(x, y)`. Also
  registered now for cheap extensibility: `sinh, cosh, tanh, sin, cos, atan`.
- **Summation:** `sum(i: <body>)`. The index name (`i` here) is sum-local; arrays
  are referenced as `arr[i]`. The iteration length is the common length of all
  arrays subscripted by the index within the body, **validated equal at compile
  time** (mismatch is a compile error). Accumulation order is `i = 0 … n-1`,
  matching the existing C++ loops so results agree to ULP.
- **Literals:** decimal and scientific notation (`2.66958e-08`).

### 2. Symbol-table contract — name resolution (3 buckets)

At bind time each identifier resolves against, in order:

0. **`let` bindings** (highest precedence) — names introduced by `let <name> = …`
   earlier in the same formula. A `let` may therefore shadow an input/constant
   name within its formula (conventional local-shadows-global scoping); this is
   intentional and harmless. The remaining buckets follow.
1. **Declared state variables** (§2b) — the names in this block's
   `state_variables`, each bound to an existing `CoolProp::parameters` key by
   `is_valid_parameter()`: `T` (K), `Dmolar` (mol/m³), `Dmass` (kg/m³),
   `molar_mass` (kg/mol), `P` (Pa), `Smolar_residual` (J/mol/K) — CoolProp's own
   spellings, whatever `keyed_output()` supports. A name **not** declared here is
   never state, so an undeclared `p` is the author's to use. Checked *before*
   constants, and declaring a name that is also a `constants` or `arrays` key is a
   **compile error**: resolving it either way would give the author a number they
   did not write.
2. **Block-declared constants** — scalars from the JSON `constants` object
   (e.g. `T_reduce`, `epsilon_over_k`, `sigma_eta`, `C`).

Arrays live in their own namespace and resolve separately: a name is looked up in
the JSON `arrays` object only in subscripted form `name[index]` inside a `sum`.

A name found in none of the buckets is a descriptive compile error
(`unknown variable '<name>' at col <n>`).

**All inputs are exposed in base SI, always** (`T` in
K, `Dmolar` in mol/m³, `P` in Pa, result in Pa·s or W/m/K). The DSL imposes no
unit handling and there is no units metadata field — a formula yields base SI by
construction. Where the underlying physics needs a non-SI quantity (e.g. a
collision-integral correlation tabulated in nm and kg/kmol), the conversion
factor (`*1e9`, `*1000`) is written explicitly in the formula, exactly as the
current C++ does.

### 2b. Thermodynamic inputs — one bucket keyed by `CoolProp::parameters`

There is no separate "intrinsic" vs "derived" split, and no list of thermodynamic
names lives in this library at all. A block declares what it reads:

```json
"state_variables": ["T", "Dmolar"]
```

Each name is resolved by `CoolProp::is_valid_parameter()`, so **the DSL's
vocabulary is CoolProp's vocabulary** — every quantity `keyed_output()` can produce
is reachable, in CoolProp's own canonical spelling (`P`, `Dmolar`, `Dmass`,
`Smolar_residual`, `Bvirial`, ...). Adding a thermodynamic quantity to the language
costs nothing: no table row, no recompile.

Whether a key is free or costs an EOS call is not the DSL's business: the compiled
`Program` reports `requiredInputs()` as a `std::vector<CoolProp::parameters>` and
reads back a value array in that order. The host — `ExpressionCorrelation::eval`,
which holds the backend — fills it with one `HEOS.keyed_output(key)` per key. The
evaluator itself stays EOS-free, and there is no per-quantity `switch` on either
side to extend.

**Why declared rather than global.** The first design resolved a curated allowlist
(`T`, `rhomolar`, `rhomass`, `molar_mass`, `p`) in one namespace shared with the
author's own constants and arrays. That has two defects, and they compound.

*Every name the language knows is reserved everywhere.* Propylene glycol's dilute
term is η₀ = Σ nᵢ·Tr^pᵢ; its exponent array had to be renamed `np` because `p` was
pressure — in a block that never reads pressure. Across the eleven shipped blocks,
local identifiers outnumber state reads about 5:1 (median 11 against 2), and **not
one of them declares `P`**. The reservation was pure cost.

**Only the canonical spelling is accepted.** `is_valid_parameter()` also resolves
CoolProp's back-compat aliases and an upper-cased form of every name: `A` is the
speed of sound, `D` is `Dmass`, `M` is `molar_mass`, `TAU` is `Tau`. Harmless in
`PropsSI`; a trap here, because a DSL formula is full of single-letter coefficient
names. An author who forgets to declare a paper coefficient `A` in `constants` is
told to add it to `state_variables` — and doing so would silently bind the speed of
sound, with the compiler's own diagnostic having steered them there. So a
declaration must match `get_parameter_information(key, "short")` exactly. One name
per quantity, and still no list kept in this library.

*The DSL invented spellings CoolProp does not use.* `rhomolar`/`rhomass`/`p` are
not `get_parameter_index()` names; CoolProp says `Dmolar`/`Dmass`/`P`. That
divergence is what created the lowercase `p` that collided. The aliases are gone:
one name per quantity, and it is CoolProp's.

Declaration inverts both. A name the block did not declare is never state, however
well CoolProp knows it, so `p` is simply the author's. `requiredInputs()` is still
read off the AST walk, but the declaration now gates it: a name reaches the walk only
if the block asked for it. (The *order* remains first-reference, not declaration
order — the declaration says what may be read, the formula says in what order.)

Dropping the aliases costs nothing ergonomically, because **renaming is already in
the language**. A `let` binds any state variable to whatever reads best — including
the source paper's own symbol, which neither CoolProp's spelling nor a second
vocabulary maintained here would have given:

```
let rho = Dmolar
let tau = Tc/T
sum(i: n[i]*rho^d[i]*tau^t[i])
```

A `let` may not shadow a *declared* state variable at all: that is a compile error.
An earlier draft leaned on the "declared but never used" rule to catch it, but that
only fires when the declared name is read *nowhere* — `let a = T*2 / let T = 500 /
a + T` compiled, with `T` meaning the state on one line and 500 on the next. The
guard is explicit now rather than emergent.

**Compatibility.** This changed the fluid-JSON contract with no version gate and no
shim: a block written against the previous format does not compile. That is
deliberate — the retired spellings were the defect, so silently accepting them is the
one outcome worth avoiding — but it puts the whole burden on the failure being
legible. The three retired names are therefore kept in a table that produces *only
error text*:

```
unknown variable 'rhomolar' -- it was this DSL's own spelling for a thermodynamic
quantity before blocks declared their inputs; CoolProp calls it 'Dmolar', and it
must be listed in this block's "state_variables"
```

A version field was considered and rejected. It would only help if v1 were still
supported, which is exactly what we do not want; an *optional* version is
indistinguishable from a new block that omitted it; and a version on this block type
alone would be the only such field in CoolProp's fluid JSON. The failure above works
for every old block, including the ones that never carried a version. If schema
versioning is wanted later it belongs on the fluid file, not on one block type.

**What is still refused, and why opt-in is not enough.** Two classes are rejected at
compile time even when explicitly declared, because these are the cases where the
author's intent cannot be honoured:

| Refused | Reason |
|---|---|
| `V`, `viscosity`, `L`, `conductivity`, `Prandtl`, `surface_tension` | `keyed_output()` re-enters the correlation being defined — unbounded recursion, at eval time, in fluid-file data |
| `T_critical`, `P_critical`, `rhomolar_critical`, `rhomass_critical` | `calc_T_critical()`/`calc_rhomolar_critical()` return the *numerical* critical point under `ENABLE_SUPERANCILLARIES` (the default), so a correlation reducing on them changes answer with configuration. Xenon missed its reference values by 7e-5 exactly this way. Freeze the paper's own value as a constant instead. |
| the `*_reducing` family, **and `Tau`/`Delta`** | A different reason, worth stating separately: superancillaries do not touch `get_reducing_state()`. A correlation writing `T_reducing` means its *own* fitted reducing parameter, which is generally not the EOS's and moves when the EOS section is revised. `Tau` and `Delta` are that same state wearing another name — `keyed_output()` computes them as `_reducing.T/_T` and `_rhomolar/_reducing.rhomolar` — so refusing one while admitting the other would be a guard with a door next to it. |
| `Phase`, `Q`, `Qmass` | Not continuous state functions: `Phase` is an enum ordinal widened to `double` (and depends on any imposed phase), `Q` is `-1` outside the dome. Both would compile into plausible-looking arithmetic. |
| `T_freeze`, `fraction_min`, `fraction_max`, `GWP20/100/500`, `ODP`, `FH`, `HH`, `PH` | Fluid metadata, not functions of the current state; several are unimplemented for HEOS and throw from inside the evaluation, where fluid-file data has no business failing. |

Three further errors are compile-time, not runtime:

- declaring a name that is also a `constants` or `arrays` key (the collision, now
  raised only when the author actually asked for both);
- declaring a name the formula never reads (it would cost a `keyed_output()` per
  evaluation and misstate the block's dependencies);
- reading a valid parameter that was *not* declared — the message names
  `state_variables` rather than saying "unknown variable".

### 3. JSON schema (additive)

A transport sub-block (`dilute`, `initial_density`, `higher_order`, `residual`,
`critical`) may use `"type": "expression"`:

```json
"higher_order": {
  "type": "expression",
  "formula": "let delta = Dmolar/rhomolar_reduce\nlet tau = T_reduce/T\nsum(i: a[i]*delta^d1[i]*tau^t1[i]*exp(gamma[i]*delta^l[i]))",
  "state_variables": ["Dmolar", "T"],
  "constants": { "T_reduce": 132.6312, "rhomolar_reduce": 10447.7 },
  "arrays": { "a": [1.072e-05], "d1": [1], "t1": [0.2], "gamma": [0], "l": [0] }
}
```

- `formula` (string, required): the DSL source. Newlines separate `let`
  statements. Yields base SI by construction (Pa·s for viscosity, W/m/K for
  conductivity).
- `constants` (object, optional): name → scalar (base SI).
- `arrays` (object, optional): name → array of numbers.

There is no units field: every exposed quantity and the result are base SI
always (§2), so a units annotation would be redundant.

### 3a. What a block returns, per stage

**Every expression block yields its own stage's contribution, in base SI, ready to
be summed.** The host adds the stages; it applies no scaling of its own to an
expression result. That is the `EMPIRICAL` convention, and it is deliberately NOT
the `RAINWATER_FRIEND` one, whose hardcoded arm returns a second viscosity virial
`B_eta` that the host then multiplies by `eta_dilute * rho`.

The consequence is worth stating plainly: a formula cannot see any other stage's
output. `eta_dilute` is a within-correlation intermediate, not a thermodynamic
input, and exposing it would mean a stage reading another stage's result — the
same self-reference the transport-output denylist exists to prevent. A correlation whose
initial-density term is defined as `eta_1 = eta_0 * B_eta` therefore **recomputes
`eta_0` inside the initial-density block**. That duplicates coefficient data in the
fluid file, which is the price of not adding a code path; the two copies sit side
by side in one file, and a divergence shows up immediately against the
correlation's published verification values.

Stages must stay honest about what they report: `calc_viscosity_dilute()` is
consumed independently by conductivity models, and `viscosity_contributions()` is
public API. Collapsing a whole correlation into one block would be simpler to
write and wrong for both reasons.

### 4. Components

| Unit | File(s) | ~LOC | Responsibility |
|---|---|---|---|
| Lexer | `include/CoolProp/expression/Lexer.h` | 100 | source → token stream |
| AST | `include/CoolProp/expression/Ast.h` | 80 | node types: `Num, Var, Index, Unary, Binary, Call, Sum, Program` |
| Parser | `include/CoolProp/expression/Parser.h`, `src/expression/Parser.cpp` | 250 | Pratt/precedence-climbing expr parser + statement layer → AST |
| Binder | (in `Parser.cpp` or `Binder` unit) | 150 | 3-bucket name resolution; validate unknown-name, arity, equal sum-array lengths, index-used-only-as-subscript; cache variable slots; record the `CoolProp::parameters` inputs referenced; descriptive errors |
| Evaluator | `include/CoolProp/expression/ExpressionCorrelation.h` | 150 | recursive walk over eval context; `std::pow/exp/log`; sum accumulation `0…n-1` |
| Host correlation | `include/CoolProp/expression/ExpressionCorrelation.h`, `src/expression/ExpressionCorrelation.cpp` | 120 | owns compiled `Program` + bound constants/arrays; `eval(HEOS&)` fills context (one `keyed_output()` per required input) and evaluates |

Dispatch: new enum values (`VISCOSITY_*_EXPRESSION`, `CONDUCTIVITY_*_EXPRESSION`)
added to the existing transport switch in `HelmholtzEOSMixtureBackend.cpp`; each
calls `ExpressionCorrelation::eval`. A new `"type":"expression"` branch in
`FluidLibrary` parses the block and constructs the `ExpressionCorrelation` at
fluid load (compile-once).

No new build dependency; pure C++17; coexists with nlohmann/json under hidden
symbol visibility; WASM-clean by construction.

### 4b. Standalone / scripting path

`expression::compile_block(json_text, context)`
(`include/CoolProp/expression/ExpressionBlock.h`) is the **single** definition of
what a `"type": "expression"` block looks like. Both consumers go through it:

- `JSONFluidLibrary::parse_expression_block()` at fluid load, and
- `expression::ExpressionBlock`, a compiled block that can be evaluated against
  any `AbstractState` *without* being grafted into a fluid file first.

Both also share the host-side evaluation primitive
`expression::evaluate_at(const Program&, AbstractState&)` — fill the program's
required inputs with one `keyed_output()` per key, then evaluate. The
fluid-library correlation (`ExpressionCorrelation::eval(HEOS&)`) is a one-line
call to it, so a block evaluated standalone and the same block loaded into a
fluid cannot disagree.

`ExpressionBlock` takes JSON *text*, not an `nlohmann::json`, deliberately: the
header is reached from the scripting wrappers, and keeping it JSON-library-free
keeps nlohmann out of those translation units. The fluid library pays a `dump()`
round-trip per expression block at load, which is once-per-fluid and exact
(nlohmann round-trips doubles losslessly).

See "Authoring from Python" below for the Python surface built on it.

### 5. Error handling

- **Compile-time (at fluid load):** lex/parse/bind errors throw
  `CoolProp::ValueError` with the formula, the message, and the column of the
  offending token. A malformed formula fails the fluid load loudly; it never
  produces a silently-wrong correlation and never crashes.
- **Eval-time, inside the formula:** numeric domain results follow `std::pow/log`
  semantics exactly (e.g. `log` of a non-positive argument → NaN/-inf as in the
  current C++), so DSL output matches the hardcoded routines bit-for-bit in
  behavior, including at domain edges. No domain guards, no exceptions.
- **Eval-time, on the way in:** the *inputs* are guarded, which is a different
  thing. `evaluate_at()` throws `CoolProp::ValueError` naming the input if a
  `keyed_output()` comes back non-finite. This is the one exception the hot path
  can raise, and it costs two predictable compares per input (unmeasurable against
  the `pow`/`exp` calls that dominate). It exists because
  `AbstractState::p()` on a state that was never `update()`d returns `-_HUGE`
  rather than throwing, and the formula would propagate that into a
  plausible-looking `-inf` answer. A state reached through any completed
  `update()` already satisfies the guard — `HelmholtzEOSMixtureBackend::post_update()`
  enforces the same `ValidNumber` predicate on `_T`/`_p`/`_rhomolar` — so the
  guard fires only on a state that was never set.

  It is a guard, not a proof of freshness: `AbstractState::set_T()` mutates `_T`
  alone and leaves `_p`/`_rhomolar` at the previous state's (finite) values, and
  `molar_mass` is a trivial parameter that is finite even on an untouched state.
  Whoever owns the state object is responsible for it being the state they meant.

## Testing — also the completeness proof

New Catch2 tag `[expression]`.

1. **Unit tests** (no EOS): lexer tokenization incl. scientific notation; parser
   precedence and `^` right-associativity; `let` scoping; `sum` index/length
   semantics; and **failure modes** — sum array-length mismatch, unknown
   variable, unknown function, arity mismatch, malformed input → clean
   `ValueError`, never a crash. Evaluator correctness on hand-checked
   expressions, including derived-var injection via a directly-populated context
   (proves the evaluator is EOS-free). **Registry path:** a formula referencing
   `p` binds to the registry, and the host computes pressure from `HEOS` and
   injects it — assert the formula's value matches `HEOS.p()` at the state.
2. **Golden regression (the gate):** for **every Tier-A form**, author its DSL
   equivalent, then for several representative fluids compare
   `ExpressionCorrelation::eval` against the existing hardcoded C++ routine across
   a `(T, ρ)` grid spanning the fluid's transport validity range. Default gate is
   relative error `< 1e-14` (~tens of ULP). Because the DSL replicates the same
   library calls (`std::pow/exp/log`) in the same accumulation order, most forms
   match far tighter; the only expected divergence class is a hand-written
   `x*x`/`x*x*x` in some C++ routines vs the DSL's `x^2`/`x^3` (→ `std::pow`).
   **Where a form uses the identical operation sequence, assert bit-exact
   (0 ULP) per-form** and reserve the `1e-14` band only for the `pow`-vs-multiply
   forms. If all Tier-A forms reproduce within this gate, the DSL is provably
   complete for the scope.

Per project convention (`CLAUDE.md`), changes here run under `[SBTL]`-style
umbrella discipline only if they touch those paths; this feature is new code, so
the relevant local sweep is `[expression]` plus the existing transport tests to
prove no regression. `./dev/ci/preflight.sh` selects scope from changed paths.

## Authoring from Python

Correlation authoring and doc writing both want the same loop: type a formula,
evaluate it at a state, look at the number. `CoolProp.CoolProp.Expression` is that
loop, with no rebuild and no fluid-file edit in between:

```python
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import Expression

blk = """{
  "type": "expression",
  "formula": "sum(i: a[i]*T^t[i])",
  "arrays": {"a": [-4.87e-7, 3.29e-8], "t": [1.0, 1.5]}
}"""

e = Expression(blk)
e.required_inputs()          # ['T']

AS = CP.AbstractState("HEOS", "R123")
AS.update(CP.DmolarT_INPUTS, 1e4, 300.0)
e.evaluate(AS)
```

- The JSON text is exactly what goes into the fluid file, so a doc example and the
  fluid it documents cannot drift.
- `required_inputs()` reports the DSL names the formula actually references — the
  quickest way to see that `p` (say) pulled the EOS into the evaluation.
- `evaluate(AS)` reads whatever state `AS` is currently sitting at, via
  `keyed_output()` — the same path the production host takes. The caller sets the
  state through the ordinary `AbstractState` API, so **any** input pair, backend,
  or mixture composition works, and the block duplicates no state-setting logic of
  its own. Re-`update()` the state and the same compiled block follows it.
- A bad formula raises `ValueError` with the column of the offending token, at
  construction time.

## Risks / open trade-offs

- **`^` semantics:** chosen as `pow`. Forms whose C++ uses the identical
  operation sequence match bit-exact; the only divergence is `x^2`/`x^3` (→
  `std::pow`) vs a hand-written `x*x` in some routines. Golden gate is `1e-14`
  relative (~tens of ULP) with per-form bit-exact assertions where applicable.
  Accepted.
- **Per-block declared state variables:** the binder records the
  `CoolProp::parameters` keys a formula declares and reads, and the host fills them
  with `keyed_output()`. This replaced a curated five-entry allowlist, which made
  Tier B additive but only one row at a time — krypton's entropy scaling cost three
  rows and a recompile for quantities `keyed_output()` already knew, and its
  reserved names collided with coefficient arrays in blocks that never used them.
  Declaration makes Tier B free and the reservation local. `P` is still exercised
  end-to-end by a dedicated unit test, proving the EOS-backed path.
- **Evaluation requires an `AbstractState`,** even for a formula that reads only
  `T` and `Dmolar` (or nothing at all): the caller must stand up a backend for
  some fluid. An earlier draft took `(T, rhomolar, fluid_name)` and could evaluate
  with no fluid at all, which suited authoring a correlation for a fluid CoolProp
  does not ship yet. **Accepted deliberately** (PR #3185 review): this code lives
  *inside* CoolProp, so depending on CoolProp's own state type costs nothing and
  buys the whole existing API — any input pair, any backend, mixtures, and a
  single shared evaluation path with the fluid-library correlation. The calculus
  would be different for an API meant to be usable outside the library.
- **Performance unproven at production scale for `sum`:** Tier-A evals are
  `pow`/`exp`-bound and compiled once; expected tens-to-hundreds of ns. If
  profiling later shows an expression in a tight solver loop, bytecode is the
  localized fallback. Not pre-optimized.

## Future work (not v1)

- Register Tier-B derived variables (`dpdrho__constT`, `cpmolar`, `cvmolar`,
  correlation length, pressure parts) and express Olchowy–Sengers / friction
  theory / Chung as data.
- Optional bytecode compilation if profiling justifies it.
- Possible migration of existing Tier-A fluid JSON to `"type":"expression"` once
  the path is proven (separate, reversible effort).
