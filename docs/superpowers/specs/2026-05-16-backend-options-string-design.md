# Backend options via factory-string JSON

**Status:** draft — for review
**Owner:** Ian Bell
**Tracking:** CoolProp-vh0 (epic)
**Date:** 2026-05-16

## Problem

Several CoolProp backends need per-instance configuration that doesn't fit
the current model:

- **SVDSBTL** (the immediate driver): the critical-patch HEOS-fallback
  needs a per-fluid bounding box in (T, p). Today there's no per-fluid
  configuration, only Configuration globals; Ian flagged this in the
  earlier SVDSBTL review:

  > I don't want to hardcode. Something automatic is needed and/or
  > something the user can override/control at runtime.

  Auto-calibration covers "automatic". This spec covers "the override".

- **BICUBIC / TTSE**: grid resolution, axis ranges, and refinement
  policies are currently process-global via `Configuration`. Per-instance
  tweaks aren't expressible.

- **REFPROP**: `ALTERNATIVE_REFPROP_LIBRARY_PATH`, `_HMX_BNC_PATH`,
  `_PATH` are all Configuration globals. Switching libraries inside a
  single Python process is awkward / impossible.

The current `Configuration` mechanism is process-global and stringly-
typed; once set, every subsequent backend instantiation inherits it,
which makes mixing two SVDSBTL backends with different policies in the
same process effectively impossible.

## Goal

Add a single, **immutable**, **per-instance** option-passing mechanism
that:

1. Threads options through the factory entry-point in a way that
   composes cleanly with `PropsSI` one-liners.
2. Is **typed and nestable** — schema can grow without grammar evolution.
3. Is **strict** — unknown keys throw at construction time; no silent
   drift across CoolProp versions.
4. Is **canonical** — the same logical options produce the same bytes,
   so options participate in cache-key hashing.
5. Is **reproducible** — `AS->build_options_json()` returns the exact
   canonical form the instance was built with; copying it into a fresh
   `factory()` call replays the construction.
6. Is **backend-agnostic** — the mechanism lives in `AbstractState`
   plumbing, not in any single backend.

Explicit non-goals:

- **No mutators.** Options are construction-time only. Want different
  options? Construct a different `AbstractState`.
- **No Configuration globals as defaults.** The factory string is the
  single source of truth. Shared defaults live in a JSON file the user
  references via `?@path.json`.
- **No per-`PropsSI`-call options.** The factory call is the only place
  options enter.

## Design

### Grammar

The existing factory string accepts an optional `?<options>` suffix:

```text
SVDSBTL&HEOS::Water
SVDSBTL&HEOS::Water?{"critical_patch":"off"}
SVDSBTL&HEOS::Water?@/path/to/cfg.json
HEOS::CarbonDioxide?{"phase_envelope":{"max_steps":500}}
```

Parsing is intentionally dumb:

```cpp
auto pos = factory_string.find_first_of('?');
if (pos == std::string::npos) {
    options_json = "{}";
    cleaned    = factory_string;
} else {
    cleaned    = factory_string.substr(0, pos);
    auto tail  = factory_string.substr(pos + 1);
    if (tail.empty())            options_json = "{}";
    else if (tail[0] == '@')     options_json = read_file(tail.substr(1));
    else                         options_json = tail;            // inline JSON
}
```

`find_first_of('?')` is called **once**. Any further `?` characters are
inside the JSON value (or inside the file path for `@…`) and are handled
by the JSON parser's quoting rules, not by the factory-string parser.

### High-level interface support (PropsSI, FORTRAN, Excel)

The factory string is parsed in exactly one place
(`AbstractState::factory()`), and `PropsSI` / `Props1SI` /
`PhaseSI` etc. forward their `FluidName` argument through that same
parser unchanged. Consequence: the `?<options>` suffix works
identically for callers stuck on the high-level interface, with **no
new entry-points and no per-language wrapper changes**:

```python
# Python
PropsSI("D", "T", 300, "P", 1e5,
        'SVDSBTL&HEOS::Water?{"critical_patch":"off"}')
```

```fortran
! FORTRAN — escape quotes per host language rules
d = PropsSI('D'//c_null_char, 'T'//c_null_char, 300.0_dp, &
            'P'//c_null_char, 1.0e5_dp, &
            'SVDSBTL&HEOS::Water?@/opt/coolprop/h2o_no_patch.json'//c_null_char)
```

```excel
=PropsSI("D";"T";300;"P";1E5;"SVDSBTL&HEOS::Water?@C:\Configs\h2o.json")
```

```matlab
% MATLAB
d = py.CoolProp.CoolProp.PropsSI('D','T',300,'P',1e5, ...
        'SVDSBTL&HEOS::Water?@~/coolprop/opts.json');
```

For inline JSON the host language's string-quoting rules apply
(escape `"` as `\"` in C-family, doubled `""` in Excel, etc.). When
the quoting becomes painful — Excel cells especially are limited to
~32k chars and don't escape gracefully — the `?@path.json` form
sidesteps escaping entirely and is the recommended path for
high-level interfaces.

### Schema (initial, for SVDSBTL)

```jsonc
{
  "schema": 1,
  "pmin": 611.655,                         // absolute Pa; lower-pressure bound
                                           // for the PT/PH/PS surfaces.
                                           // Default = fluid p_triple; must be
                                           // >= p_triple (no sat boundary below).
  "grid": { "NT": 200, "NR": 800, "rank": 20 },
  "properties": {
    "transport": "auto"                    // auto | on | off
  },
  "critical_patch": {
    "mode": "auto",                        // auto | off | fixed
    "source": null,                        // null = same as truth source
    "tolerance": 1e-3,                     // auto-cal target residual
    "metric": "D",                         // D | A | H — auto-cal residual metric
    "bbox": null                           // [Tlo, Thi, plo, phi] when mode=="fixed"
  }
}
```

Nested objects let `critical_patch.tolerance` and `critical_patch.bbox`
coexist without flat-key collisions and leave room for `regions: {"R3":
{...}}` later. `schema` is the migration anchor — bump it, write a
migrator if a key moves.

Schema is published as a JSON Schema document (Draft 2020-12) in
`schemas/svdsbtl_options.schema.json`, rendered into the Sphinx docs.

### Validation

A shared helper `validate_against_schema(json_value, schema_document)`
runs at factory time. Default mode is **strict**:

- Unknown keys → throw.
- Type mismatch → throw.
- Required keys missing → throw.

No `strict=false` toggle in the initial release. If we want one later
it's a per-backend `"strict": false` in the JSON itself, not a global.

### Canonical serialization

Cache hashing and reproducibility both need a single canonical form per
logical-options-value. We sort object keys recursively, preserve
array order, and use RapidJSON's default number formatting. This
is *not* strict [RFC 8785 JCS][jcs] (no NFC string normalisation,
no ECMAScript number rounding) — but it's "deterministic within a
single CoolProp build", which is what cache hashing needs since the
canonical form never appears in any on-wire / cross-implementation
context.

The canonical bytes are the input to:

- **`AbstractState::build_options_json()`** — returns the canonical
  string for the instance.
- **Cache filename hashing** — for caching backends (SVDSBTL today),
  the cache filename gets a 16-hex-char [FNV-1a 64][fnv1a] prefix
  of the canonical bytes alongside fluid / source / input_pair.
  FNV-1a 64 is the same hash family CoolProp already uses for
  superancillary EOS-freshness verification (see the
  `Superancillary source_eos_hash` test in
  `src/Tests/CoolProp-Tests.cpp`) — no new dependencies, identical
  determinism guarantees across compilers. Drops the int-keyed
  input_pair fragility (CoolProp-b6v) at the same time.

[jcs]: https://datatracker.ietf.org/doc/html/rfc8785
[fnv1a]: http://www.isthe.com/chongo/tech/comp/fnv/

### API surface

```cpp
// AbstractState.h — new virtual, default returns empty.
virtual std::string build_options_json() const { return ""; }
```

```cpp
// AbstractStateGenerator interface — new virtual with default that
// rejects options.  Backwards-compatible: existing generators stay
// working untouched.
virtual AbstractState* get_AbstractState(
    const std::vector<std::string>& fluid_names,
    const rapidjson::Value& options) {
    if (options.IsNull() || (options.IsObject() && options.ObjectEmpty())) {
        return get_AbstractState(fluid_names);
    }
    throw NotImplementedError(
        "backend " + name() + " does not accept options");
}
```

```cpp
// Existing factory entry-point unchanged in signature; gains
// options-string parsing internally.
AbstractState* AbstractState::factory(const std::string& factory_string,
                                       const std::vector<std::string>& fluid_names);
```

`get_AbstractState(fluid_names)` (no-options) stays — backends that
ignore options need zero code changes.

### Per-backend opt-in

A backend opts in by:

1. Publishing a JSON Schema in `schemas/<backend>_options.schema.json`.
2. Overriding `get_AbstractState(fluid_names, options)` to validate
   against its schema and construct accordingly.
3. (Optional) Overriding `build_options_json()` to return the canonical
   form it was constructed with.

SVDSBTL ships as the first consumer in this PR. BICUBIC, TTSE, REFPROP
become follow-up patches; each one deprecates the corresponding
Configuration globals with a deprecation warning + a migration table in
the docs.

### `@path` indirection — file-path edge cases

When the suffix begins with `@`, the rest is a filesystem path read
verbatim:

```cpp
auto path = tail.substr(1);              // strip leading '@'
// NO further ? splitting on the path
```

Path may itself contain `?`, `&`, etc. Path resolution: relative paths
resolve against the process CWD; `~` expansion uses the same helper
already used for `ALTERNATIVE_REFPROP_PATH`.

If the file doesn't exist or isn't readable, throw a `FileIOError`
with the full path quoted, so misconfiguration surfaces immediately.

### Auto-calibration loop (referenced by `critical_patch.mode=auto`)

Runs at table-build time, once per (fluid, input_pair). Output is a
single `[Tlo, Thi, plo, phi]` tuple stored in the cache header.

Sketch:

1. Probe set: sample (T, p) on a hexagonal lattice in an oversized
   initial bbox (e.g. 0.85·T_c .. 1.20·T_c × 0.50·p_c .. 1.30·p_c).
2. For each candidate (Tlo, Thi, plo, phi), evaluate `metric` residual
   between SVDSBTL prediction and source backend over the probe set.
3. Binary-search-shrink along each axis independently from the
   oversized initial box. Stop when the residual hits `tolerance` or
   further shrinking blows past it.
4. Persist the final tuple in the cache header alongside the SVD
   blobs.

User who passes `critical_patch.mode="fixed"` + `bbox=[...]` skips the
loop and the supplied tuple goes straight into the cache. `mode="off"`
disables the patch entirely; the cache stores `null` for the bbox.

## Test contract

The parser is the riskiest piece. Tests live in
`src/Tests/CoolProp-Tests-FactoryOptions.cpp` (new file) and cover:

### Parsing edge cases

- No `?` → options = `{}`.
- Bare `?` → options = `{}`.
- `?{}` → options = `{}`.
- Inline JSON with nested objects round-trips.
- `@<path>` reads file; missing file throws `FileIOError`.
- `@<path>` where the path itself contains `?` and `&` — the path is
  read verbatim, no re-splitting.
- Garbage after `?` (non-JSON, non-`@`) → throws `ValueError` with a
  message that includes the offending tail.

### `?` characters inside JSON values

A mandatory subset (driven by Ian's earlier review note):

- `?{"hint":"what?"}` — `?` inside a string value parses correctly.
- `?{"q":"?"}` — `?` is the entire string value.
- `?{"meta":{"url":"http://x.com/path?id=1&q=2"}}` — `?` and `&`
  together inside a URL-style string.
- `?{"regex":"^[A-Z]+\\?$"}` — `?` after an escaped backslash.
- `?{"chain":"first?second?third"}` — multiple `?` characters inside
  one string.

Each test asserts the parsed value equals the expected Python
equivalent.

### Schema validation

- Unknown top-level key → throws with key name in message.
- Unknown nested key → throws with full dotted path
  (e.g. `critical_patch.unknown_field`).
- Type mismatch → throws with expected vs actual type.
- Required key missing → throws with key name.
- Default values fill in for unspecified keys.

### Canonical serialization round-trip

- Build instance with `?{"grid":{"NT":200,"rank":20,"NR":800}}` and
  with `?{"grid":{"rank":20,"NR":800,"NT":200}}`. Both must produce the
  same `build_options_json()` output (sorted keys) and the same
  FNV-1a 64 cache-filename prefix.
- Round-trip: `build_options_json()` output, fed back into the factory,
  produces an instance with identical `build_options_json()`.

### High-level-interface pass-through

- `PropsSI("D", "T", 300, "P", 1e5,
   'SVDSBTL&HEOS::Water?{"critical_patch":"off"}')` returns a finite
   density and matches the result of constructing the same backend
   via `AbstractState::factory()` with the same string.
- `PropsSI(..., 'SVDSBTL&HEOS::Water?@/tmp/cfg.json')` reads the file
   and applies the options. A missing file throws the same
   `FileIOError` PropsSI surfaces as a CoolProp error message — no
   silent fallback to defaults.
- A round-trip from `build_options_json()` through PropsSI's
   FluidName argument produces a backend with identical canonical
   options.

### SVDSBTL behavioural tests

- `critical_patch=off` → backend instance has the patch disabled;
  queries inside the would-be patch return NaN (no fallback).
- `critical_patch=fixed&bbox=[...]` → cache stores the supplied bbox
  verbatim; auto-cal loop is not invoked.
- `critical_patch=auto` (default) → auto-cal runs; resulting bbox is
  persistent in the cache header.

## Migration plan

| Generation | Backend | Action |
|---|---|---|
| Concurrent with this PR | SVDSBTL | First consumer; full schema, critical-patch + auto-cal landed |
| Follow-up #1 | BICUBIC | Deprecate `BICUBIC_*` Configuration globals; equivalent JSON keys land in the schema; emit deprecation warning when globals are set |
| Follow-up #2 | TTSE | Same as BICUBIC |
| Follow-up #3 | REFPROP | Deprecate `ALTERNATIVE_REFPROP_*` globals in favor of `?{"library_path":"…","hmx_path":"…"}` |
| Future | PCSAFT, SuperAncillary, HEOS | Per-backend judgment; many have nothing tunable today |

Each follow-up is one PR with its own migration notes and at least one
release-cycle of overlap before the Configuration globals are removed.

## Out-of-scope (deferred)

- Per-region SVDSBTL options (e.g. different rank in R3). Schema leaves
  room (`"regions": {...}`) but no plumbing.
- Programmatic options-builder API in language wrappers. The factory
  string is the universal entry-point; wrappers can wrap it but the
  spec doesn't require any per-wrapper helper.
- Read-only `Configuration` snapshot persisted into the cache. Cache
  already gets the options canonical form; that's sufficient for
  reproducibility.

## Decisions made during PR A

1. **JSON library: RapidJSON** (already vendored for the fluid catalog).
   nlohmann would offer a nicer API at a build-time cost; not worth it.
2. **Schema validator: RapidJSON's built-in** (`rapidjson::SchemaDocument`).
   Strict-mode = `additionalProperties:false` in each per-backend schema.
3. **Canonical-form hash: FNV-1a 64** (same hash family CoolProp already
   uses to stamp `source_eos_hash` for superancillary freshness; 16-hex
   prefix on cache filenames).
4. **Schema storage:** per-backend JSON Schemas live in
   `include/CoolProp/schemas/<backend>.schema.json` plus a generated
   `*.h` companion that exposes the schema as a string literal
   (`kSVDSBTLOptionsSchemaJson`), so the validator has no runtime file
   dependency.

## Open questions

1. **Deprecation cadence for Configuration globals** — One release with
   the JSON key + a deprecation warning when the global is set, then
   removal one release later? Leaving as TBD pending Ian's preference.

## Acceptance criteria

- [ ] Parser splits on first `?` only; handles inline JSON, `@path`,
      empty / missing tail.
- [ ] All `?`-in-JSON test cases pass.
- [ ] SVDSBTL `?{"critical_patch":"off"}` correctly disables the patch
      with no auto-cal invocation.
- [ ] SVDSBTL `?{"critical_patch":"fixed","bbox":[...]}` builds a cache
      with the supplied bbox in the header.
- [ ] SVDSBTL `?{"critical_patch":"auto"}` runs the auto-cal loop and
      persists the result.
- [ ] `AS->build_options_json()` round-trips: feed output back to
      `factory()`, get bit-identical canonical form.
- [ ] Cache filename hashing makes two logically-equal options produce
      the same cache filename; logically-different options produce
      different filenames.
- [ ] At least one BICUBIC / TTSE follow-up PR filed (not landed) to
      verify the mechanism generalizes.
- [ ] PropsSI / Props1SI / PhaseSI accept the `?<options>` suffix
      verbatim in their `FluidName` argument — both inline JSON and
      `@path` forms exercised in tests, Python wrapper through to the
      C++ entry-point.
