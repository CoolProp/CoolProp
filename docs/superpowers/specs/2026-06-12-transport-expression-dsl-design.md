# Transport Property Expression DSL — Design

**Date:** 2026-06-12
**Status:** Design approved, pending spec review
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
| Derived vars | 4-bucket name resolution incl. a derived-state registry, **seeded with one entry — `p` (pressure, Pa) — in v1** as the worked example |

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

### 2. Symbol-table contract — name resolution (4 buckets)

At bind time each identifier resolves against, in order:

0. **`let` bindings** (highest precedence) — names introduced by `let <name> = …`
   earlier in the same formula. A `let` may therefore shadow an intrinsic/constant
   name within its formula (conventional local-shadows-global scoping); this is
   intentional and harmless. The remaining buckets follow.
1. **Intrinsic state** (hardwired, always present, EOS-free): `T` (K),
   `rhomolar` (mol/m³), `rhomass` (kg/m³), `molar_mass` (kg/mol). These are the
   independent state variables and pure fluid metadata — available without
   evaluating the EOS. Pressure is **not** here: `p` requires the EOS at the
   current state, so it lives in the derived-state registry (§2b).
2. **Block-declared constants** — scalars from the JSON `constants` object
   (e.g. `T_reduce`, `epsilon_over_k`, `sigma_eta`, `C`).
3. **Block-declared arrays** — vectors from the JSON `arrays` object, usable only
   in subscripted form `name[index]` inside a `sum`.
4. **Derived-state registry** (§2b).

A name found in none of the four buckets is a descriptive compile error
(`unknown variable '<name>' at col <n>`).

**All intrinsic and derived quantities are exposed in base SI, always** (`T` in
K, `rhomolar` in mol/m³, `p` in Pa, result in Pa·s or W/m/K). The DSL imposes no
unit handling and there is no units metadata field — a formula yields base SI by
construction. Where the underlying physics needs a non-SI quantity (e.g. a
collision-integral correlation tabulated in nm and kg/kmol), the conversion
factor (`*1e9`, `*1000`) is written explicitly in the formula, exactly as the
current C++ does.

### 2b. Derived-state registry (one entry — `p` — in v1)

A host-side table mapping a canonical DSL name → a getter on
`HelmholtzEOSMixtureBackend`, for state-dependent quantities the EOS must compute
(e.g. a future `smolar_residual`, `cpmolar`, `dpdrho__constT`). Resolution and
dependency-tracking are built into the binder now (~30 lines). **In v1 the table
is seeded with exactly one entry — `p` (pressure, Pa)** — both because pressure
is a genuinely useful always-derivable quantity and because it exercises the
registry path end-to-end (a worked example) before any Tier-B input depends on
it. No Tier-A *completeness* form uses `p`; it is proven by a dedicated unit
test, not by the golden regression.

Mechanism: when the binder resolves an identifier to bucket 4, it records that
the compiled `Program` depends on that derived name. At eval time the host layer
(which holds `HEOS`) computes the recorded dependencies and writes them into the
eval context before the tree walk. The evaluator itself stays EOS-free.

Extension cost, made explicit:

- The v1 seed is one such line:
  ```cpp
  registry.add("p", [](HelmholtzEOSMixtureBackend& H){ return H.p(); });
  ```
  Adding any **new** derived quantity later is the same one-line pattern plus a
  CoolProp recompile, e.g.
  `registry.add("smolar_residual", [](auto& H){ return H.smolar_residual(); });`.
- After that, **every fluid and formula** can use that name as pure JSON data,
  no further code. The grammar/parser never changes.

This is the designed bridge to Tier B: each Tier-B input becomes a registry
one-liner rather than a new hardcoded routine + enum + switch case.

### 3. JSON schema (additive)

A transport sub-block (`dilute`, `initial_density`, `higher_order`, `residual`,
`critical`) may use `"type": "expression"`:

```json
"higher_order": {
  "type": "expression",
  "formula": "let delta = rhomolar/rhomolar_reduce\nlet tau = T_reduce/T\nsum(i: a[i]*delta^d1[i]*tau^t1[i]*exp(gamma[i]*delta^l[i]))",
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

### 4. Components

| Unit | File(s) | ~LOC | Responsibility |
|---|---|---|---|
| Lexer | `include/CoolProp/expression/Lexer.h` | 100 | source → token stream |
| AST | `include/CoolProp/expression/Ast.h` | 80 | node types: `Num, Var, Index, Unary, Binary, Call, Sum, Program` |
| Parser | `include/CoolProp/expression/Parser.h`, `src/expression/Parser.cpp` | 250 | Pratt/precedence-climbing expr parser + statement layer → AST |
| Binder | (in `Parser.cpp` or `Binder` unit) | 150 | 4-bucket name resolution; validate unknown-name, arity, equal sum-array lengths, index-used-only-as-subscript; cache variable slots; record derived-var deps; descriptive errors |
| Evaluator | `include/CoolProp/expression/ExpressionCorrelation.h` | 150 | recursive walk over eval context; `std::pow/exp/log`; sum accumulation `0…n-1` |
| Host correlation | `include/CoolProp/expression/ExpressionCorrelation.h`, `src/expression/ExpressionCorrelation.cpp` | 120 | owns compiled `Program` + bound constants/arrays; `eval(HEOS&)` fills context (intrinsics + derived deps) and evaluates |

Dispatch: new enum values (`VISCOSITY_*_EXPRESSION`, `CONDUCTIVITY_*_EXPRESSION`)
added to the existing transport switch in `HelmholtzEOSMixtureBackend.cpp`; each
calls `ExpressionCorrelation::eval`. A new `"type":"expression"` branch in
`FluidLibrary` parses the block and constructs the `ExpressionCorrelation` at
fluid load (compile-once).

No new build dependency; pure C++17; coexists with nlohmann/json under hidden
symbol visibility; WASM-clean by construction.

### 5. Error handling

- **Compile-time (at fluid load):** lex/parse/bind errors throw
  `CoolProp::ValueError` with the formula, the message, and the column of the
  offending token. A malformed formula fails the fluid load loudly; it never
  produces a silently-wrong correlation and never crashes.
- **Eval-time:** numeric domain results follow `std::pow/log` semantics exactly
  (e.g. `log` of a non-positive argument → NaN/-inf as in the current C++), so
  DSL output matches the hardcoded routines bit-for-bit in behavior, including at
  domain edges. No exceptions thrown on the hot path.

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

## Risks / open trade-offs

- **`^` semantics:** chosen as `pow`. Forms whose C++ uses the identical
  operation sequence match bit-exact; the only divergence is `x^2`/`x^3` (→
  `std::pow`) vs a hand-written `x*x` in some routines. Golden gate is `1e-14`
  relative (~tens of ULP) with per-form bit-exact assertions where applicable.
  Accepted.
- **Derived registry, one entry in v1:** ~30 lines of binder machinery plus a
  single registered getter (`p`). Accepted as cheap insurance that makes Tier B
  additive (a registration table) rather than a binder rewrite, and the `p` seed
  proves the path end-to-end.
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
