# Exclude `detail/json.h` from installed headers — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stop shipping `include/CoolProp/detail/json.h` (which pulls `nlohmann/json.hpp` + valijson) in the installed header set, and add a fail-closed CI guard that it stays out.

**Architecture:** CMake-only change — add a precise `REGEX "detail/json\.h$" EXCLUDE` to the four header-install rules. No source edits (no installed header `#include`s `detail/json.h`; the two `cpjson` wrappers define the shared enum independently under a guard macro, with no include coupling). A new `dev/ci/check-installed-headers.sh` installs to a temp prefix and asserts `detail/json.h` is absent while `detail/rapidjson.h` is present, wired into preflight reusing the existing `build_shared`.

**Tech Stack:** CMake `install(DIRECTORY … FILES_MATCHING)`, Bash, the preflight gate (`dev/ci/preflight.sh`).

**Spec:** `docs/superpowers/specs/2026-06-07-exclude-detail-json-from-install-design.md`
**Issue:** `CoolProp-rxxx` (narrowed). **Branch:** `ihb/json-detail-install-exclude`.

---

## Task 1: Add the `detail/json.h` install exclude (4 rules)

**Files:**
- Modify: `CMakeLists.txt:631-632, 668-669, 857-858` (the three `install(DIRECTORY …/include/CoolProp …)` clauses)
- Modify: `wrappers/Python/CMakeLists.txt:348-355` (the `install(DIRECTORY "${ROOT_DIR}/include/" …)` clause)

Only the `install(DIRECTORY …/include/CoolProp …)` clause in each CMakeLists trio ships the `detail/` subtree; the companion `install(DIRECTORY …/include/ … PATTERN "CoolProp" EXCLUDE)` clause already skips it, so it needs no change.

- [ ] **Step 1: CMakeLists.txt — static_library clause (~631)**

Change:
```cmake
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/CoolProp
            DESTINATION static_library/include FILES_MATCHING PATTERN "*.h")
```
to:
```cmake
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/CoolProp
            DESTINATION static_library/include FILES_MATCHING PATTERN "*.h"
            REGEX "detail/json\\.h$" EXCLUDE)
```

- [ ] **Step 2: CMakeLists.txt — shared_library clause (~668)**

Change:
```cmake
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/CoolProp
            DESTINATION shared_library/include FILES_MATCHING PATTERN "*.h")
```
to:
```cmake
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/CoolProp
            DESTINATION shared_library/include FILES_MATCHING PATTERN "*.h"
            REGEX "detail/json\\.h$" EXCLUDE)
```

- [ ] **Step 3: CMakeLists.txt — cmake-install clause (~857)**

Change:
```cmake
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/CoolProp
          DESTINATION "${CMAKE_INSTALL_PREFIX}/usr/include" FILES_MATCHING PATTERN "*.h")
```
to:
```cmake
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/CoolProp
          DESTINATION "${CMAKE_INSTALL_PREFIX}/usr/include" FILES_MATCHING PATTERN "*.h"
          REGEX "detail/json\\.h$" EXCLUDE)
```

- [ ] **Step 4: wrappers/Python/CMakeLists.txt (~348)**

The clause currently reads:
```cmake
install(DIRECTORY "${ROOT_DIR}/include/"
    DESTINATION CoolProp/include
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.hpp"
    PATTERN "*_JSON.h" EXCLUDE
    PATTERN "*_JSON_z.h" EXCLUDE
)
```
Add the REGEX exclude before the closing paren:
```cmake
install(DIRECTORY "${ROOT_DIR}/include/"
    DESTINATION CoolProp/include
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.hpp"
    PATTERN "*_JSON.h" EXCLUDE
    PATTERN "*_JSON_z.h" EXCLUDE
    REGEX "detail/json\\.h$" EXCLUDE
)
```
(Read the actual file first to confirm exact surrounding lines; the `*.hpp` line may or may not be present — preserve whatever excludes already exist and just add the REGEX line.)

- [ ] **Step 5: Manually verify the exclude works (shared config)**

Run:
```bash
cmake -B build_shared -S . -DCOOLPROP_SHARED_LIBRARY=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build_shared -j8
rm -rf /tmp/cp-install-check && cmake --install build_shared --prefix /tmp/cp-install-check
echo "json.h present? ->"; find /tmp/cp-install-check -path '*/detail/json.h'
echo "rapidjson.h present? ->"; find /tmp/cp-install-check -path '*/detail/rapidjson.h'
```
Expected: the `json.h` find prints **nothing**; the `rapidjson.h` find prints **one path**.

- [ ] **Step 6: Commit**

```bash
git add CMakeLists.txt wrappers/Python/CMakeLists.txt
git restore --staged .beads/issues.jsonl 2>/dev/null || true
git commit -m "build(json): exclude detail/json.h from installed headers (CoolProp-rxxx)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: Add the fail-closed install-hygiene guard script

**Files:**
- Create: `dev/ci/check-installed-headers.sh`

- [ ] **Step 1: Write the script**

Create `dev/ci/check-installed-headers.sh` with mode `0755`:
```bash
#!/usr/bin/env bash
# Fail if forbidden internal headers are shipped in the installed tree.
#
# Enforced today: include/CoolProp/detail/json.h (which #includes
# nlohmann/json.hpp + valijson) must NOT be installed.  detail/rapidjson.h MUST
# still be installed (positive control — it is reachable via Configuration.h and
# is only removed at Phase Final; its presence proves the exclude is specific,
# not over-broad).
#
# Fail-closed: a failed `cmake --install`, an empty install tree, or a missing
# positive-control header is a FAILURE, never a vacuous pass (same discipline as
# check-json-symbols.sh — no `|| true` that masks a real failure).
set -euo pipefail

BUILD_DIR="${1:-build_shared}"
if [ ! -d "$BUILD_DIR" ]; then
    echo "FAIL: build dir '$BUILD_DIR' not found (configure+build a shared build first)" >&2
    exit 2
fi

PREFIX="$(mktemp -d)"
trap 'rm -rf "$PREFIX"' EXIT

if ! cmake --install "$BUILD_DIR" --prefix "$PREFIX" >/tmp/install-headers-check.log 2>&1; then
    echo "FAIL: cmake --install '$BUILD_DIR' failed (see /tmp/install-headers-check.log)" >&2
    exit 1
fi

# Guard against a vacuous pass: a real install ships many headers.
NHEADERS="$(find "$PREFIX" -name '*.h' | grep -c . || true)"
if [ "$NHEADERS" -eq 0 ]; then
    echo "FAIL: install tree under $PREFIX has no headers — wrong build dir or empty install? Cannot validate." >&2
    exit 1
fi

# Positive control: detail/rapidjson.h MUST ship (exclude must not be over-broad).
if ! find "$PREFIX" -path '*/detail/rapidjson.h' | grep -q .; then
    echo "FAIL: detail/rapidjson.h is absent from the install tree — exclude too broad or install layout changed." >&2
    exit 1
fi

# The assertion: detail/json.h must NOT ship.
LEAKED="$(find "$PREFIX" -path '*/detail/json.h' || true)"
if [ -n "$LEAKED" ]; then
    echo "FAIL: detail/json.h is shipped in the installed headers (it pulls nlohmann/json.hpp + valijson):" >&2
    printf '%s\n' "$LEAKED" >&2
    exit 1
fi

echo "OK: detail/json.h not installed; detail/rapidjson.h present (${NHEADERS} headers from ${BUILD_DIR})"

# TODO(superancillary de-leak): once superancillary.h no longer #includes
# nlohmann/json.hpp, broaden this to assert that NO shipped header pulls
# nlohmann/valijson (grep the installed tree), making it the install-side
# companion to the symbol-leak gate.
```

- [ ] **Step 2: Make it executable**

```bash
chmod +x dev/ci/check-installed-headers.sh
```

- [ ] **Step 3: Run it against the shared build from Task 1**

Run: `./dev/ci/check-installed-headers.sh build_shared`
Expected: `OK: detail/json.h not installed; detail/rapidjson.h present (N headers from build_shared)`.

- [ ] **Step 4: Prove the guard actually catches a regression (fail-closed check)**

Temporarily defeat the exclude to confirm the guard fails (do NOT commit this):
```bash
# Temporarily install WITHOUT the exclude by pointing at a stale config is hard;
# instead simulate a leak: copy json.h into a fresh prefix that already passed,
# then re-run the find logic by hand to confirm the detection branch works:
PREFIX=$(mktemp -d); cmake --install build_shared --prefix "$PREFIX" >/dev/null 2>&1
mkdir -p "$PREFIX/shared_library/include/CoolProp/detail" && cp include/CoolProp/detail/json.h "$PREFIX/shared_library/include/CoolProp/detail/json.h"
find "$PREFIX" -path '*/detail/json.h' && echo "DETECTION OK: find located the planted json.h" || echo "DETECTION BROKEN"
rm -rf "$PREFIX"
```
Expected: prints the planted path + `DETECTION OK`. (This validates the detection logic; the real script would `exit 1` on it.)

- [ ] **Step 5: Commit**

```bash
git add dev/ci/check-installed-headers.sh
git restore --staged .beads/issues.jsonl 2>/dev/null || true
git commit -m "ci(json): fail-closed guard that detail/json.h stays out of installed headers (CoolProp-rxxx)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: Wire the guard into preflight

**Files:**
- Modify: `dev/ci/preflight.sh` (add a step immediately after the `json-symbols` step, ~line 190, reusing `build_shared`; and add `install-headers` to the `--skip` usage examples near the top)

- [ ] **Step 1: Add the preflight step after the json-symbols block**

Immediately after the json-symbols `step`'s closing `fi` (the block that ends around line 190, just before `# ---------- check 3: Catch2 tests`), insert:
```bash
# ---------- check 2c: installed-header hygiene -----------------------
#
# Companion to the symbol-leak gate on the install side: assert that
# detail/json.h (which pulls nlohmann/json.hpp + valijson) is not shipped
# in the installed headers.  Reuses build_shared (built by the json-symbols
# step above).  Fail-closed lives in the script.
step "installed-header hygiene"
if skip_check install-headers; then
    skip "install-headers" "--skip=install-headers"
elif skip_check build; then
    skip "install-headers" "--skip=build (shared build needed)"
elif [ ! -d build_shared ]; then
    skip "install-headers" "build_shared not available (run without --skip=json-symbols)"
elif ./dev/ci/check-installed-headers.sh build_shared; then
    ok "install-headers (detail/json.h not shipped)"
else
    fail "install-headers (detail/json.h shipped, or install failed; see /tmp/install-headers-check.log)"
fi
```
(Match the surrounding `step`/`ok`/`skip`/`fail`/`skip_check` helper style exactly — read the json-symbols block first to confirm the helper names.)

- [ ] **Step 2: Add `install-headers` to the `--skip` usage comment**

Near the top of `dev/ci/preflight.sh` where `--skip=` examples are listed (around line 15-16, next to `json-symbols`), add `install-headers` to the documented skip targets so the flag is discoverable.

- [ ] **Step 3: Run preflight's new step in isolation to confirm it passes**

Run (reusing the already-built build_shared): `./dev/ci/check-installed-headers.sh build_shared && echo "STEP OK"`
Expected: the `OK:` line + `STEP OK`.

- [ ] **Step 4: Commit**

```bash
git add dev/ci/preflight.sh
git restore --staged .beads/issues.jsonl 2>/dev/null || true
git commit -m "ci(json): run install-header hygiene check in preflight (CoolProp-rxxx)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Final Verification (before PR)

- [ ] **Step 1: Full preflight**

Run: `./dev/ci/preflight.sh`
Expected: the new `installed-header hygiene` step passes; clang-format/build/json-symbols pass. (cppcheck/clang-tidy: this PR touches no `.cpp`/`.h` source, only CMake + a shell script, so the C++ linters have nothing in the diff to flag — they should skip or pass.)

- [ ] **Step 2: Sanity that the build/tests are unaffected**

This is an install-time + CI-script change only; no compiled code changed. If `build_catch/CatchTestRunner` exists, a quick `./build_catch/CatchTestRunner "[INCOMP]"` confirms nothing broke, but a full suite is not required for an install-rule change.

- [ ] **Step 3: bd bookkeeping**

```bash
# Narrow rxxx to the install-exclude:
bd update CoolProp-rxxx --title "Exclude detail/json.h from installed headers (decouple done; CMake exclude + hygiene guard)"
# File the superancillary de-leak as a new sibling under the epic, blocking Phase Final:
bd create --title "De-leak nlohmann::json from installed superancillary.h (explicit instantiation)" \
  --type=task --priority=2 --parent=CoolProp-xa8w \
  --description "superancillary.h is a header-only TEMPLATE (SuperAncillary<ArrayType>) that #includes nlohmann/json.hpp and exposes SuperAncillary(const nlohmann::json&) + SuperAncillary(const std::string&) (parses via nlohmann::json::parse) in its installed public API, reachable via CoolPropFluid.h -> superancillary.h. CoolPropFluid.h::get_superanc() constructs it INLINE, so the construction path must move out-of-line into libCoolProp via explicit template instantiation for the concrete ArrayType(s) (std::vector<double>, and verify Eigen::ArrayXd), keeping the perf-critical evaluation code inline. Then drop the nlohmann include + json ctor from the installed header. Blocks Phase Final."
# Link the new issue to block Phase Final (use the ID printed by bd create):
# bd dep <new-id> --blocks CoolProp-xa8w.3
```
(Capture the new issue's real ID from `bd create` output and run the `bd dep … --blocks CoolProp-xa8w.3` with it. Do NOT reuse a parent ID.)

- [ ] **Step 4: Adversarial review + finish**

Per CLAUDE.md, dispatch the review subagent against the diff vs `master` (focus: the script's fail-closed guarantees — is there any path where it prints OK / exits 0 while `detail/json.h` is actually shipped? e.g. wrong build dir, `cmake --install` partial failure swallowed, `grep -c` under `set -e`, the positive-control `find | grep -q` logic). Address findings, then push + `gh pr create`, then re-review any post-review commits.

---

## Self-Review (plan vs spec)

**Spec coverage:** §3 CMake excludes → Task 1 (4 rules); §4 fail-closed guard script → Task 2; preflight wiring → Task 3; §5 verification → Task 1 Step 5 (manual) + Task 2 Step 3-4 (guard) + Final Verification; §6 bd bookkeeping → Final Verification Step 3; §7 out-of-scope (superancillary, rapidjson.h, broader gate) → not implemented, the broader gate is a documented `TODO` in the script.

**Placeholder scan:** none — the one `TODO` in the script is an intentional, documented forward-pointer (not a plan gap). The bd `create` step requires capturing a generated ID at runtime (noted explicitly), which is inherent to `bd`.

**Consistency:** the matcher `REGEX "detail/json\\.h$" EXCLUDE` is identical across all 4 install rules and the script asserts exactly that file's absence + `detail/rapidjson.h`'s presence — the assertion matches the change.
