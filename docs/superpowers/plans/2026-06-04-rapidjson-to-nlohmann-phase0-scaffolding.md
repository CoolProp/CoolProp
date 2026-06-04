# RapidJSON → nlohmann/json — Phase 0 (Seam & Scaffolding) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stand up nlohmann/json + Valijson (both symbol-hidden), a CI symbol-leak gate, a Python build-time schema validator, and an nlohmann-backed `cpjson` helper header — additively, without migrating any loader, so the tree stays green.

**Architecture:** All three JSON libraries (RapidJSON during transition, nlohmann/json, Valijson) are fetched via CPM in `cmake/dependencies.cmake`. A single non-public wrapper header (`include/CoolProp/detail/json.h`) includes nlohmann under `#pragma GCC visibility push(hidden)` and re-implements the `cpjson` helpers on `nlohmann::json`. nlohmann's own versioned inline namespace (`nlohmann::json_abi_v3_x_x`) plus hidden visibility provide the no-leak guarantee; the spec's custom-namespace rename (§4 layer 2) is **dropped** because it breaks Valijson's nlohmann adapter, which references `nlohmann::json` literally.

**Tech Stack:** C++ (CoolProp), CMake + CPM.cmake, nlohmann/json, Valijson, Catch2, Python `jsonschema`, `nm`/`objdump`.

**Spec:** `docs/superpowers/specs/2026-06-04-rapidjson-to-nlohmann-migration-design.md`

**Scope note:** This is Phase 0 only. Phases 1…N (per-loader migrations) and Phase Final (delete RapidJSON) each get their own plan — see the roadmap at the end. Phase 0 adds infrastructure and leaves RapidJSON fully in place and in use.

**Prerequisites for every task:** configure once with
`cmake -B build_catch -S . -DCOOLPROP_CATCH_MODULE=ON -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release`
and build the runner with
`cmake --build build_catch --target CatchTestRunner -j8`.
New CPM packages and new test files require a **re-configure** (re-run the `cmake -B` command) before building.

---

## File Structure

| File | Responsibility | Action |
|------|----------------|--------|
| `cmake/dependencies.cmake` | CPM package declarations | Modify — add `nlohmann_json`, `valijson` |
| `CMakeLists.txt` | include dirs, test source list | Modify — add nlohmann/valijson include dirs; drop `externals/nlohmann-json`; register new test file |
| `include/CoolProp/detail/json.h` | nlohmann wrapper: hidden-visibility include + `cpjson` helpers on `nlohmann::json` + Valijson schema validation | Create |
| `src/Tests/CoolProp-Tests-JSONHelpers.cpp` | Catch2 tests for the new helpers + schema validation + no-escaping-exception contract | Create |
| `dev/ci/check-json-symbols.sh` | Post-build gate: fail if a shared product exports `nlohmann`/`rapidjson` symbols | Create |
| `dev/ci/preflight.sh` | Pre-push gate | Modify — invoke the symbol + python-schema checks |
| `dev/validate_fluid_schemas.py` | Build-time validation of `dev/*.json` against schemas via `jsonschema` | Create |
| `externals/nlohmann-json/` | Stale single-header placeholder | Delete |

---

## Task 1: Add nlohmann/json and Valijson via CPM

**Files:**
- Modify: `cmake/dependencies.cmake:33-38` (insert after the rapidjson block)
- Modify: `CMakeLists.txt:336` (add include dirs)

- [ ] **Step 1: Add the two CPM packages**

In `cmake/dependencies.cmake`, immediately after the existing `rapidjson` block (the `CPMAddPackage(NAME rapidjson …)` ending at line 38), insert:

```cmake
# nlohmann/json — replacement for rapidjson (GH: RapidJSON→nlohmann migration).
# Header-only; included via the hidden-visibility wrapper include/CoolProp/detail/json.h.
CPMAddPackage(
  NAME nlohmann_json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG        v3.12.0
  DOWNLOAD_ONLY  YES   # header-only; we only need the include dir
)

# Valijson — header-only JSON-Schema (draft-7) validator that validates an
# nlohmann::json instance directly via its bundled adapter. Used for runtime
# validation of user-supplied PC-SAFT / cubic fluids.
CPMAddPackage(
  NAME valijson
  GIT_REPOSITORY https://github.com/tristanpenman/valijson.git
  GIT_TAG        v1.0.6
  DOWNLOAD_ONLY  YES   # header-only; adapters/ + the core headers are include-only
)
```

- [ ] **Step 2: Add their include directories**

In `CMakeLists.txt`, immediately after line 336 (`list(APPEND APP_INCLUDE_DIRS "${rapidjson_SOURCE_DIR}/include")`), add:

```cmake
list(APPEND APP_INCLUDE_DIRS "${nlohmann_json_SOURCE_DIR}/include")
list(APPEND APP_INCLUDE_DIRS "${valijson_SOURCE_DIR}/include")
```

- [ ] **Step 3: Re-configure to verify the packages resolve**

Run:
```bash
cmake -B build_catch -S . -DCOOLPROP_CATCH_MODULE=ON -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
```
Expected: configure succeeds; output mentions fetching `nlohmann_json` and `valijson` (or reuses the CPM cache). No error.

- [ ] **Step 4: Confirm the headers are on disk**

Run:
```bash
ls build_catch/_deps/nlohmann_json-src/include/nlohmann/json.hpp \
   build_catch/_deps/valijson-src/include/valijson/adapters/nlohmann_json_adapter.hpp
```
Expected: both paths exist and are listed.

- [ ] **Step 5: Commit**

```bash
git add cmake/dependencies.cmake CMakeLists.txt
git commit -m "build(json): fetch nlohmann/json + valijson via CPM"
```

---

## Task 2: Create the hidden-visibility nlohmann wrapper header (smoke level)

This task creates `include/CoolProp/detail/json.h` with only the hidden-visibility include and a single trivial helper, plus a smoke test, to prove the include path and visibility pragma compile. Helpers are added in Task 3.

**Files:**
- Create: `include/CoolProp/detail/json.h`
- Create: `src/Tests/CoolProp-Tests-JSONHelpers.cpp`
- Modify: `CMakeLists.txt:2171` (register the test file)

- [ ] **Step 1: Write the failing test**

Create `src/Tests/CoolProp-Tests-JSONHelpers.cpp`:

```cpp
// Catch2 tests for the nlohmann-backed cpjson helpers and Valijson schema
// validation (RapidJSON→nlohmann migration, Phase 0).

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include <string>
#    include <vector>

#    include "CoolProp/detail/json.h"
#    include "CoolProp/Exceptions.h"

TEST_CASE("cpjson nlohmann wrapper parses a string", "[json]") {
    nlohmann::json j = cpjson::parse(R"({"a": 1.5})");
    REQUIRE(j.at("a").get<double>() == Approx(1.5));
}

#endif  // ENABLE_CATCH
```

- [ ] **Step 2: Register the test file in the Catch runner**

In `CMakeLists.txt`, immediately after the `CoolProp-Tests-TermCacheProfile.cpp` block (ending line 2171), add:

```cmake
  list(APPEND APP_SOURCES
       "${CMAKE_CURRENT_SOURCE_DIR}/src/Tests/CoolProp-Tests-JSONHelpers.cpp")
```

- [ ] **Step 3: Run the test to verify it fails to build**

Run:
```bash
cmake -B build_catch -S . -DCOOLPROP_CATCH_MODULE=ON -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build_catch --target CatchTestRunner -j8
```
Expected: compile error — `CoolProp/detail/json.h` not found (file does not exist yet).

- [ ] **Step 4: Create the wrapper header (minimal)**

Create `include/CoolProp/detail/json.h`:

```cpp
#ifndef COOLPROP_DETAIL_JSON_H
#define COOLPROP_DETAIL_JSON_H

// Internal, NON-INSTALLED wrapper around nlohmann/json.
//
// Symbol-leak strategy (see migration spec §4): nlohmann is included under
// hidden ELF/Mach-O visibility so none of its (weak/inline) symbols are
// exported from CoolProp's shared products. Combined with nlohmann's own
// versioned inline namespace (nlohmann::json_abi_v3_x_x), this prevents the
// cross-library ODR/symbol clashes that motivated the migration. Do NOT
// rename nlohmann's namespace here: Valijson's nlohmann adapter references
// nlohmann::json literally and a rename would break it.

#include "CoolProp/Exceptions.h"
#include "CoolProp/detail/tools.h"  // for CoolProp::format / CoolPropDbl

#include <string>
#include <string_view>
#include <vector>

#if defined(__GNUC__) || defined(__clang__)
#    pragma GCC visibility push(hidden)
#endif

#include <nlohmann/json.hpp>

#if defined(__GNUC__) || defined(__clang__)
#    pragma GCC visibility pop
#endif

namespace cpjson {

/// Parse a JSON-formatted string into an nlohmann::json document.
/// Throws CoolProp::ValueError (never a raw nlohmann exception) on failure.
inline nlohmann::json parse(std::string_view text) {
    try {
        return nlohmann::json::parse(text);
    } catch (const std::exception& e) {
        throw CoolProp::ValueError(std::string("Unable to parse JSON: ") + e.what());
    }
}

}  // namespace cpjson

#endif  // COOLPROP_DETAIL_JSON_H
```

- [ ] **Step 5: Build and run the smoke test**

Run:
```bash
cmake --build build_catch --target CatchTestRunner -j8
./build_catch/CatchTestRunner "[json]"
```
Expected: build succeeds; the `cpjson nlohmann wrapper parses a string` test PASSES.

- [ ] **Step 6: Commit**

```bash
git add include/CoolProp/detail/json.h src/Tests/CoolProp-Tests-JSONHelpers.cpp CMakeLists.txt
git commit -m "feat(json): add hidden-visibility nlohmann wrapper header + smoke test"
```

---

## Task 3: Re-implement the cpjson helpers on nlohmann::json

Port every helper the loaders use from `detail/rapidjson.h` to `detail/json.h`, keeping the **same `cpjson::` names and call signatures** so later loader migrations are a near-mechanical include swap. A given translation unit includes only one of the two headers.

**Files:**
- Modify: `include/CoolProp/detail/json.h`
- Modify: `src/Tests/CoolProp-Tests-JSONHelpers.cpp`

- [ ] **Step 1: Write the failing tests for the getters**

Append to `src/Tests/CoolProp-Tests-JSONHelpers.cpp`, before the final `#endif`:

```cpp
TEST_CASE("cpjson getters extract typed members", "[json]") {
    nlohmann::json j = cpjson::parse(R"({
        "i": 7, "x": 2.5, "b": true, "s": "hi",
        "darr": [1.0, 2.0, 3.0],
        "sarr": ["a", "b"],
        "d2d": [[1.0, 2.0], [3.0]]
    })");

    REQUIRE(cpjson::get_integer(j, "i") == 7);
    REQUIRE(cpjson::get_double(j, "x") == Approx(2.5));
    REQUIRE(cpjson::get_bool(j, "b") == true);
    REQUIRE(cpjson::get_string(j, "s") == "hi");

    std::vector<double> darr = cpjson::get_double_array(j, "darr");
    REQUIRE(darr.size() == 3);
    REQUIRE(darr[2] == Approx(3.0));

    std::vector<std::string> sarr = cpjson::get_string_array(j, "sarr");
    REQUIRE(sarr.size() == 2);
    REQUIRE(sarr[1] == "b");

    std::vector<std::vector<double>> d2d = cpjson::get_double_array2D(j.at("d2d"));
    REQUIRE(d2d.size() == 2);
    REQUIRE(d2d[0][1] == Approx(2.0));
}

TEST_CASE("cpjson getters throw ValueError on missing/mistyped members", "[json]") {
    nlohmann::json j = cpjson::parse(R"({"x": "not a number"})");
    REQUIRE_THROWS_AS(cpjson::get_double(j, "missing"), CoolProp::ValueError);
    REQUIRE_THROWS_AS(cpjson::get_double(j, "x"), CoolProp::ValueError);
}
```

- [ ] **Step 2: Run to verify failure**

Run:
```bash
cmake --build build_catch --target CatchTestRunner -j8
```
Expected: compile error — `get_integer`, `get_double`, etc. are not members of `cpjson`.

- [ ] **Step 3: Implement the helpers**

In `include/CoolProp/detail/json.h`, inside `namespace cpjson`, after `parse(...)`, add:

```cpp
/// Serialize an nlohmann::json value to a pretty-printed string.
inline std::string json2string(const nlohmann::json& v) {
    return v.dump(4);
}

inline int get_integer(const nlohmann::json& v, const std::string& m) {
    auto it = v.find(m);
    if (it == v.end()) throw CoolProp::ValueError(format("Does not have member [%s]", m.c_str()));
    if (!it->is_number_integer()) throw CoolProp::ValueError(format("Member [%s] is not an integer", m.c_str()));
    return it->get<int>();
}

inline double get_double(const nlohmann::json& v, const std::string& m) {
    auto it = v.find(m);
    if (it == v.end()) throw CoolProp::ValueError(format("Does not have member [%s]", m.c_str()));
    if (!it->is_number()) throw CoolProp::ValueError(format("Member [%s] is not a number", m.c_str()));
    return it->get<double>();
}

inline bool get_bool(const nlohmann::json& v, const std::string& m) {
    auto it = v.find(m);
    if (it == v.end()) throw CoolProp::ValueError(format("Does not have member [%s]", m.c_str()));
    if (!it->is_boolean()) throw CoolProp::ValueError(format("Member [%s] is not a boolean", m.c_str()));
    return it->get<bool>();
}

inline std::string get_string(const nlohmann::json& v, const std::string& m) {
    auto it = v.find(m);
    if (it == v.end()) throw CoolProp::ValueError(format("Does not have member [%s]", m.c_str()));
    if (!it->is_string()) throw CoolProp::ValueError(format("Member [%s] is not a string", m.c_str()));
    return it->get<std::string>();
}

inline std::vector<double> get_double_array(const nlohmann::json& v) {
    if (!v.is_array()) throw CoolProp::ValueError("input is not an array");
    std::vector<double> out;
    out.reserve(v.size());
    for (const auto& el : v) {
        if (!el.is_number()) throw CoolProp::ValueError("input is not a number");
        out.push_back(el.get<double>());
    }
    return out;
}

inline std::vector<double> get_double_array(const nlohmann::json& v, const std::string& m) {
    auto it = v.find(m);
    if (it == v.end()) throw CoolProp::ValueError(format("Does not have member [%s]", m.c_str()));
    return get_double_array(*it);
}

inline std::vector<CoolPropDbl> get_long_double_array(const nlohmann::json& v) {
    if (!v.is_array()) throw CoolProp::ValueError("input is not an array");
    std::vector<CoolPropDbl> out;
    out.reserve(v.size());
    for (const auto& el : v) {
        if (!el.is_number()) throw CoolProp::ValueError("input is not a number");
        out.push_back(static_cast<CoolPropDbl>(el.get<double>()));
    }
    return out;
}

inline std::vector<CoolPropDbl> get_long_double_array(const nlohmann::json& v, const std::string& name) {
    auto it = v.find(name);
    if (it == v.end()) throw CoolProp::ValueError(format("Does not have member [%s]", name.c_str()));
    return get_long_double_array(*it);
}

inline std::vector<std::vector<double>> get_double_array2D(const nlohmann::json& v) {
    if (!v.is_array()) throw CoolProp::ValueError("input is not an array");
    std::vector<std::vector<double>> out;
    for (const auto& row : v) {
        if (!row.is_array()) throw CoolProp::ValueError(format("input \"%s\" is not a 2D array", json2string(v).c_str()));
        out.push_back(get_double_array(row));
    }
    return out;
}

inline std::vector<std::vector<CoolPropDbl>> get_long_double_array2D(const nlohmann::json& v) {
    if (!v.is_array()) throw CoolProp::ValueError("input is not a 2D array");
    std::vector<std::vector<CoolPropDbl>> out;
    for (const auto& row : v) {
        if (!row.is_array()) throw CoolProp::ValueError("input is not a 2D array");
        out.push_back(get_long_double_array(row));
    }
    return out;
}

inline std::vector<std::string> get_string_array(const nlohmann::json& v) {
    if (!v.is_array()) throw CoolProp::ValueError("input is not an array");
    std::vector<std::string> out;
    out.reserve(v.size());
    for (const auto& el : v) {
        if (!el.is_string()) throw CoolProp::ValueError("input is not a string");
        out.push_back(el.get<std::string>());
    }
    return out;
}

inline std::vector<std::string> get_string_array(const nlohmann::json& v, const std::string& m) {
    auto it = v.find(m);
    if (it == v.end()) throw CoolProp::ValueError(format("Does not have member [%s]", m.c_str()));
    return get_string_array(*it);
}
```

> Note: the rapidjson-era `set_*` writers (`set_double_array`, `set_string`, etc.) and the `to_string`/`get_information` introspection helpers are intentionally NOT ported here — they are tied to rapidjson's allocator model and are replaced by idiomatic `nlohmann::json` assignment (`j["key"] = vec;`) at the call sites during each loader migration. Port any that a specific loader still needs when that loader's plan is written.

- [ ] **Step 4: Run to verify pass**

Run:
```bash
cmake --build build_catch --target CatchTestRunner -j8
./build_catch/CatchTestRunner "[json]"
```
Expected: all `[json]` tests PASS.

- [ ] **Step 5: Commit**

```bash
git add include/CoolProp/detail/json.h src/Tests/CoolProp-Tests-JSONHelpers.cpp
git commit -m "feat(json): port cpjson getter helpers onto nlohmann::json"
```

---

## Task 4: Add Valijson schema validation + exception-translation contract

**Files:**
- Modify: `include/CoolProp/detail/json.h`
- Modify: `src/Tests/CoolProp-Tests-JSONHelpers.cpp`

- [ ] **Step 1: Write the failing tests**

Append to `src/Tests/CoolProp-Tests-JSONHelpers.cpp`, before the final `#endif`:

```cpp
namespace {
constexpr const char* kJsonSchema = R"({
  "$schema": "http://json-schema.org/draft-07/schema#",
  "type": "object",
  "required": ["name", "T"],
  "properties": {
    "name": {"type": "string"},
    "T":    {"type": "number", "exclusiveMinimum": 0}
  }
})";
}

TEST_CASE("cpjson::validate_schema accepts a conforming document", "[json]") {
    std::string err;
    cpjson::schema_validation_code code =
        cpjson::validate_schema(kJsonSchema, R"({"name": "water", "T": 300.0})", err);
    REQUIRE(code == cpjson::SCHEMA_VALIDATION_OK);
    REQUIRE(err.empty());
}

TEST_CASE("cpjson::validate_schema rejects a non-conforming document with a message", "[json]") {
    std::string err;
    cpjson::schema_validation_code code =
        cpjson::validate_schema(kJsonSchema, R"({"name": "water", "T": -5.0})", err);
    REQUIRE(code == cpjson::SCHEMA_NOT_VALIDATED);
    REQUIRE_FALSE(err.empty());  // human-comprehensible error preserved
}

TEST_CASE("cpjson::validate_schema reports malformed input JSON, never leaks nlohmann exception", "[json]") {
    std::string err;
    cpjson::schema_validation_code code =
        cpjson::validate_schema(kJsonSchema, R"({not valid json)", err);
    REQUIRE(code == cpjson::INPUT_INVALID_JSON);
}
```

- [ ] **Step 2: Run to verify failure**

Run:
```bash
cmake --build build_catch --target CatchTestRunner -j8
```
Expected: compile error — `validate_schema` / `schema_validation_code` not members of `cpjson`.

- [ ] **Step 3: Implement Valijson-backed validation**

In `include/CoolProp/detail/json.h`, add the Valijson includes inside the existing hidden-visibility block (so Valijson's nlohmann adapter is also hidden). Replace the include block so it reads:

```cpp
#if defined(__GNUC__) || defined(__clang__)
#    pragma GCC visibility push(hidden)
#endif

#include <nlohmann/json.hpp>
#include <valijson/adapters/nlohmann_json_adapter.hpp>
#include <valijson/schema.hpp>
#include <valijson/schema_parser.hpp>
#include <valijson/validator.hpp>

#if defined(__GNUC__) || defined(__clang__)
#    pragma GCC visibility pop
#endif
```

Then, inside `namespace cpjson`, append:

```cpp
enum schema_validation_code
{
    SCHEMA_VALIDATION_OK = 0,
    SCHEMA_INVALID_JSON,
    INPUT_INVALID_JSON,
    SCHEMA_NOT_VALIDATED
};

/// Validate a JSON-formatted input string against a JSON-formatted draft-07
/// schema string. On a validation failure, `errstr` receives a
/// human-comprehensible description. Never propagates a raw nlohmann or
/// Valijson exception to the caller.
inline schema_validation_code validate_schema(std::string_view schemaJson, std::string_view inputJson, std::string& errstr) {
    nlohmann::json schemaDoc;
    try {
        schemaDoc = nlohmann::json::parse(schemaJson);
    } catch (const std::exception& e) {
        errstr = std::string("Invalid schema: ") + e.what();
        return SCHEMA_INVALID_JSON;
    }

    nlohmann::json inputDoc;
    try {
        inputDoc = nlohmann::json::parse(inputJson);
    } catch (const std::exception& e) {
        errstr = std::string("Invalid input json: ") + e.what();
        return INPUT_INVALID_JSON;
    }

    valijson::Schema schema;
    valijson::SchemaParser parser;
    valijson::adapters::NlohmannJsonAdapter schemaAdapter(schemaDoc);
    try {
        parser.populateSchema(schemaAdapter, schema);
    } catch (const std::exception& e) {
        errstr = std::string("Invalid schema: ") + e.what();
        return SCHEMA_INVALID_JSON;
    }

    valijson::Validator validator;
    valijson::ValidationResults results;
    valijson::adapters::NlohmannJsonAdapter inputAdapter(inputDoc);
    if (!validator.validate(schema, inputAdapter, &results)) {
        std::string msg;
        valijson::ValidationResults::Error error;
        while (results.popError(error)) {
            for (const std::string& ctx : error.context) msg += ctx;
            msg += ": " + error.description + "\n";
        }
        errstr = msg;
        return SCHEMA_NOT_VALIDATED;
    }
    return SCHEMA_VALIDATION_OK;
}
```

- [ ] **Step 4: Run to verify pass**

Run:
```bash
cmake --build build_catch --target CatchTestRunner -j8
./build_catch/CatchTestRunner "[json]"
```
Expected: all `[json]` tests PASS, including the three validation cases.

- [ ] **Step 5: Commit**

```bash
git add include/CoolProp/detail/json.h src/Tests/CoolProp-Tests-JSONHelpers.cpp
git commit -m "feat(json): add Valijson-backed cpjson::validate_schema + exception-translation tests"
```

---

## Task 5: Symbol-leak CI gate

Build CoolProp as a shared library and fail if any exported dynamic symbol matches `nlohmann` or `rapidjson`. Note: this Phase-0 wrapper header isn't yet referenced by library sources, so the check should currently PASS trivially — its value grows as loaders migrate. We still wire it now so the invariant is enforced from the first migration onward.

**Files:**
- Create: `dev/ci/check-json-symbols.sh`
- Modify: `dev/ci/preflight.sh`

- [ ] **Step 1: Create the check script**

Create `dev/ci/check-json-symbols.sh` (mode `0755`):

```bash
#!/usr/bin/env bash
# Fail if CoolProp's shared library exports any nlohmann or rapidjson symbol.
#
# Visibility attributes in a static .a only take effect once linked into a
# shared object, so this MUST inspect a shared product, never the archive.
set -euo pipefail

LIB="${1:-}"
if [[ -z "${LIB}" || ! -f "${LIB}" ]]; then
    echo "usage: $0 <path-to-shared-library (.so/.dylib)>" >&2
    exit 2
fi

case "$(uname -s)" in
    Darwin) SYMS="$(nm -gU "${LIB}" 2>/dev/null | c++filt || true)" ;;
    *)      SYMS="$(nm -D --defined-only "${LIB}" 2>/dev/null | c++filt || true)" ;;
esac

LEAKS="$(printf '%s\n' "${SYMS}" | grep -E 'nlohmann|rapidjson' || true)"
if [[ -n "${LEAKS}" ]]; then
    echo "FAIL: shared library exports JSON-library symbols:" >&2
    printf '%s\n' "${LEAKS}" | head -50 >&2
    exit 1
fi
echo "OK: no nlohmann/rapidjson symbols exported from ${LIB}"
```

- [ ] **Step 2: Make it executable and run it against a built shared library**

Run:
```bash
chmod +x dev/ci/check-json-symbols.sh
cmake -B build_shared -S . -DCOOLPROP_SHARED_LIBRARY=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build_shared -j8
SHARED_LIB="$(find build_shared -name 'libCoolProp.so' -o -name 'libCoolProp.dylib' | head -1)"
./dev/ci/check-json-symbols.sh "${SHARED_LIB}"
```
Expected: prints `OK: no nlohmann/rapidjson symbols exported …` (Phase 0: the wrapper isn't linked into any source yet, so this passes trivially).

- [ ] **Step 3: Wire the check into preflight**

Read `dev/ci/preflight.sh` to find where it builds and where it runs the final gates. After its existing build step, add a stanza that locates the shared library (building it if the preflight only built the Catch runner) and invokes `dev/ci/check-json-symbols.sh "${SHARED_LIB}"`, failing preflight on non-zero exit. Match the script's existing logging/step style (it uses numbered step echoes and a `--skip` mechanism — add a `--skip-json-symbols` flag consistent with the others).

- [ ] **Step 4: Run preflight to confirm the new gate executes**

Run:
```bash
./dev/ci/preflight.sh
```
Expected: preflight runs the JSON-symbol step and reports OK (or, if the shared build is heavy, the step is reached and passes). No failure attributable to the new gate.

- [ ] **Step 5: Commit**

```bash
git add dev/ci/check-json-symbols.sh dev/ci/preflight.sh
git commit -m "ci(json): gate exported nlohmann/rapidjson symbols on the shared library"
```

---

## Task 6: Python build-time schema validation

Validate the data files that ship embedded against their schemas, in Python, as the real correctness gate. This runs independently of the C++ build.

**Files:**
- Create: `dev/validate_fluid_schemas.py`
- Modify: `dev/ci/preflight.sh`

- [ ] **Step 1: Write the validator with a self-test**

Create `dev/validate_fluid_schemas.py`:

```python
#!/usr/bin/env python3
"""Validate CoolProp's source JSON data files against their JSON schemas.

This is the build-time correctness gate for embedded fluid data
(RapidJSON->nlohmann migration spec, section 5). Run from the repo root.
Exits non-zero on the first validation failure.
"""
import json
import sys
from pathlib import Path

try:
    import jsonschema
except ImportError:
    sys.exit("jsonschema is required: pip install jsonschema")

REPO = Path(__file__).resolve().parent.parent

# (data file, schema file) pairs. Each schema validates each top-level item
# in the corresponding data file's array/object as appropriate.
PAIRS = [
    (REPO / "dev/pcsaft/all_pcsaft_fluids.json", REPO / "dev/pcsaft/pcsaft_fluids_schema.json"),
    (REPO / "dev/cubics/all_cubic_fluids.json", REPO / "dev/cubics/cubic_fluids_schema.json"),
    (REPO / "dev/mixtures/mixture_departure_functions.json",
     REPO / "dev/mixtures/mixture_departure_functions_schema.json"),
]


def validate_pair(data_path: Path, schema_path: Path) -> int:
    if not data_path.exists() or not schema_path.exists():
        print(f"SKIP (missing): {data_path.name} / {schema_path.name}")
        return 0
    schema = json.loads(schema_path.read_text())
    data = json.loads(data_path.read_text())
    items = data if isinstance(data, list) else [data]
    failures = 0
    for i, item in enumerate(items):
        try:
            jsonschema.validate(instance=item, schema=schema)
        except jsonschema.ValidationError as e:
            failures += 1
            print(f"FAIL {data_path.name}[{i}]: {e.message}")
    if failures == 0:
        print(f"OK  {data_path.name} ({len(items)} items) vs {schema_path.name}")
    return failures


def main() -> int:
    total = sum(validate_pair(d, s) for d, s in PAIRS)
    if total:
        print(f"\n{total} schema validation failure(s)")
        return 1
    print("\nAll fluid data files validate against their schemas.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 2: Run it to verify it validates the real data**

Run:
```bash
chmod +x dev/validate_fluid_schemas.py
uvx --from jsonschema --with jsonschema python dev/validate_fluid_schemas.py || \
  python3 dev/validate_fluid_schemas.py
```
Expected: each existing (data, schema) pair prints `OK …`; the script exits 0. If a schema's array-vs-object shape differs from the assumption, adjust `validate_pair`'s `items` handling for that file and re-run until it reports OK against the current committed data.

- [ ] **Step 3: Wire into preflight (changed-path aware)**

In `dev/ci/preflight.sh`, in the section that selects checks by changed paths, add: when any `dev/**/*.json` or `dev/**/*_schema.json` file changed, run `dev/validate_fluid_schemas.py` (resolved via `uvx`/`python3` as the script already does for other Python tools) and fail preflight on non-zero exit. Add a `--skip-schema-validate` flag consistent with the existing skip flags.

- [ ] **Step 4: Run preflight to confirm**

Run:
```bash
./dev/ci/preflight.sh
```
Expected: the schema-validation step runs (or is correctly skipped when no `dev/*.json` changed) and passes.

- [ ] **Step 5: Commit**

```bash
git add dev/validate_fluid_schemas.py dev/ci/preflight.sh
git commit -m "ci(json): add Python build-time fluid-schema validation gate"
```

---

## Task 7: Delete the stale nlohmann placeholder

**Files:**
- Delete: `externals/nlohmann-json/` (the single-header placeholder)
- Modify: `CMakeLists.txt:331`

- [ ] **Step 1: Confirm nothing depends on the placeholder path**

Run:
```bash
grep -rn "externals/nlohmann-json" --include=*.txt --include=*.cmake --include=*.py . | grep -v build
grep -rn 'include <nlohmann/json' src include | grep -v _deps
```
Expected: the only `externals/nlohmann-json` reference is `CMakeLists.txt:331`; no `src/`/`include/` file relies on the old placeholder path (nlohmann now resolves via the CPM include dir from Task 1).

- [ ] **Step 2: Remove the include-dir line**

In `CMakeLists.txt`, delete line 331:
```cmake
list(APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/externals/nlohmann-json")
```

- [ ] **Step 3: Delete the placeholder directory**

Run:
```bash
git rm -r externals/nlohmann-json
```

- [ ] **Step 4: Re-configure, build, and run the JSON tests**

Run:
```bash
cmake -B build_catch -S . -DCOOLPROP_CATCH_MODULE=ON -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build_catch --target CatchTestRunner -j8
./build_catch/CatchTestRunner "[json]"
```
Expected: configure + build succeed (nlohmann resolves via CPM); all `[json]` tests PASS.

- [ ] **Step 5: Commit**

```bash
git add CMakeLists.txt
git commit -m "build(json): drop stale externals/nlohmann-json placeholder (now via CPM)"
```

---

## Task 8: Phase-0 verification sweep

- [ ] **Step 1: Full fast suite stays green**

Run:
```bash
cmake --build build_catch --target CatchTestRunner -j8
./build_catch/CatchTestRunner "[!slow]"
```
Expected: PASS. RapidJSON is untouched and all loaders still use it; nothing regressed.

- [ ] **Step 2: Symbol gate on the shared library**

Run:
```bash
cmake --build build_shared -j8
SHARED_LIB="$(find build_shared -name 'libCoolProp.so' -o -name 'libCoolProp.dylib' | head -1)"
./dev/ci/check-json-symbols.sh "${SHARED_LIB}"
```
Expected: `OK: no nlohmann/rapidjson symbols exported …`.

- [ ] **Step 3: Pre-PR adversarial review (per CLAUDE.md, MANDATORY before any PR)**

Invoke the `superpowers:code-reviewer` subagent against the diff vs `origin/master`. Address or justify each blocking finding.

- [ ] **Step 4: Preflight + push**

Run:
```bash
./dev/ci/preflight.sh
git push -u origin ihb/rapidjson-to-nlohmann-migration
```
Expected: preflight passes; push succeeds.

---

## Self-Review (completed during authoring)

- **Spec coverage (Phase 0 scope):** CPM deps → Task 1; localized hidden visibility + wrapper → Task 2; cpjson-on-nlohmann helpers → Task 3; Valijson runtime validation + exception-translation test → Task 4; symbol-leak gate on shared product → Task 5; Python build-time validation → Task 6; delete `externals/nlohmann-json` → Task 7. Spec §4 layer-2 (custom namespace) is intentionally dropped (breaks Valijson) and documented in the Architecture note. Per-loader migration and final RapidJSON deletion are out of Phase-0 scope (roadmap below).
- **Placeholder scan:** none — every code/step is concrete. The two preflight-edit steps (Tasks 5 & 6 Step 3) describe an edit to a file whose exact current contents must be read at execution time; they specify the required behavior, insertion point, and flag name precisely rather than guessing line numbers in a 15 KB script.
- **Type consistency:** `cpjson::parse`, `get_integer/get_double/get_bool/get_string`, `get_double_array`, `get_string_array`, `get_long_double_array`, `get_double_array2D`, `get_long_double_array2D`, `json2string`, `schema_validation_code` enum (`SCHEMA_VALIDATION_OK`, `SCHEMA_INVALID_JSON`, `INPUT_INVALID_JSON`, `SCHEMA_NOT_VALIDATED`), and `validate_schema(string_view, string_view, string&)` are used consistently across tasks and mirror the existing `detail/rapidjson.h` names so later loader swaps are mechanical.

---

## Roadmap: subsequent phases (each gets its own plan)

These are **not** part of this plan; listed so the whole migration is legible.

- **Phase 1 — Helmholtz loader** (`FluidLibrary.cpp/.h`, `fluids/Helmholtz.h`, `fluids/Ancillaries.h`): swap parse path to nlohmann, convert JSON-typed constructors/`to_json` to non-installed signatures, add `all_fluids.json` first-load benchmark.
- **Phase 2 — Incompressible loader** (`IncompressibleLibrary.cpp/.h`).
- **Phase 3 — PC-SAFT loader** (`PCSAFTLibrary.cpp/.h`, `fluids/PCSAFTFluid.h`): user-fluid validation now via `cpjson::validate_schema`.
- **Phase 4 — Cubics/UNIFAC loader** (`CubicsLibrary.cpp`, `UNIFACLibrary.h`).
- **Phase 5 — SVDSBTL options + Configuration** (`SVDSBTLBackend.cpp`, `Configuration.cpp/.h`, `SchemaValidation.cpp/.h`): public `SchemaValidation.h`/`Configuration.h` go string-in/string-out.
- **Phase Final — Remove RapidJSON:** drop the `rapidjson` CPM block + include dir, delete `include/CoolProp/detail/rapidjson.h` and `rapidjson_include.h`, add the `detail/` install `EXCLUDE` (now that no public header includes `detail/`), grep-clean any residual `rapidjson::`, confirm the symbol gate is green.
