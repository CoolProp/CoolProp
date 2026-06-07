# De-leak nlohmann from installed superancillary.h — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove `#include "nlohmann/json.hpp"` and every `nlohmann::json` mention from the installed header `include/CoolProp/superancillary/superancillary.h`, so consumers no longer need nlohmann to compile against CoolProp's public headers.

**Architecture:** Keep all templated/perf-critical evaluation inline in the header; move only the JSON-consuming construction out-of-line into a new non-installed `src/superancillary.cpp`, explicitly instantiated for the one `ArrayType` ever used (`std::vector<double>`). Because `ChebyshevApproximation1D` is non-assignable (`const` member), the out-of-line string ctor delegates to a private `SuperAncillary(LoadedData&&)` ctor that member-initializes — a plain-typed `LoadedData` bundle carries the parsed pieces across the seam with no nlohmann in the header.

**Tech Stack:** C++17 templates + explicit instantiation, nlohmann/json (confined to the `.cpp`), Catch2 v3, CMake (`src/*.cpp` is auto-globbed — no CMake change).

**Spec:** `docs/superpowers/specs/2026-06-07-superancillary-deleak-nlohmann-design.md`
**Issue:** `CoolProp-xa8w.4`. **Branch:** `ihb/superancillary-deleak-nlohmann`.

---

## Task 1: Move JSON construction out-of-line; de-leak the header

This is one atomic change — the header change, the new `.cpp`, and the test retargeting must land together (the header *declares* a ctor the `.cpp` *defines*, and the tests use the removed json ctor). Commit once, at the end, after the build is green.

**Files:**
- Modify: `include/CoolProp/superancillary/superancillary.h` (drop nlohmann include L46; remove `loader` L883-897 + json ctor L938-954; add `LoadedData` + delegating ctor; declare string ctor L959)
- Create: `src/superancillary.cpp`
- Modify: `src/Tests/CoolProp-Tests.cpp:3399, 3456, 3773`

- [ ] **Step 1: Read the current header regions first**

Read `include/CoolProp/superancillary/superancillary.h` lines 44-48, 805-960 to confirm the exact current text of the include, the member declarations (the order: `m_rhoL, m_rhoV, m_p`, the `m_hL…m_invlnp` optionals, `m_Tmin, m_Tcrit_num, m_rhocrit_num, m_pmin, m_pmax, m_check_points, m_lazy_build_mutex`, then a build-stamp member), the public `CheckPoint` struct (~L813), the private `loader` (L883-897), the `SuperAncillary(const nlohmann::json&)` ctor (L938-954), and the `SuperAncillary(const std::string&)` inline ctor (L959).

- [ ] **Step 2: Remove the nlohmann include**

In `superancillary.h`, delete line 46: `#include "nlohmann/json.hpp"`. (Verify nothing else in the header needs it after Steps 3-5.)

- [ ] **Step 3: Delete the private `loader` member (L883-897)**

Remove the entire `auto loader(const nlohmann::json& j, const std::string& key) { … }` method (and its preceding doc comment block at ~L880-882). Its logic moves to `src/superancillary.cpp` (Step 6).

- [ ] **Step 4: Replace the json ctor with a `LoadedData` bundle + private delegating ctor**

Replace the public `SuperAncillary(const nlohmann::json& j) : … { … }` (L938-954, including its `/// Reading in a data structure…` doc comment at L936-937) with the following. Place the public `struct LoadedData` and the private delegating ctor; the `LoadedData` field types reference `ChebyshevExpansion<ArrayType>` and `CheckPoint`, both already defined above in this header:

```cpp
   public:
    /// Plain-typed bundle of already-parsed superancillary data. The non-installed
    /// src/superancillary.cpp parses the JSON into one of these and delegates to the
    /// private ctor below, keeping nlohmann::json out of this installed header.
    struct LoadedData
    {
        std::vector<ChebyshevExpansion<ArrayType>> rhoL, rhoV, p;
        double Tcrit_num = 0;
        double rhocrit_num = 0;
        std::vector<CheckPoint> check_points;
    };

   private:
    /// Member-initializing ctor from a parsed LoadedData bundle. Member-init (not
    /// assignment) is required because ChebyshevApproximation1D is non-assignable
    /// (it has a const member). Init-list order matches member declaration order.
    explicit SuperAncillary(LoadedData&& d)
      : m_rhoL(std::move(d.rhoL)),
        m_rhoV(std::move(d.rhoV)),
        m_p(std::move(d.p)),
        m_Tmin(m_p.xmin()),
        m_Tcrit_num(d.Tcrit_num),
        m_rhocrit_num(d.rhocrit_num),
        m_pmin(m_p.eval(m_p.xmin())),
        m_pmax(m_p.eval(m_p.xmax())),
        m_check_points(std::move(d.check_points)) {}

   public:
```

IMPORTANT: confirm the `private:`/`public:` transitions keep the rest of the class intact — the existing private member block (`m_rhoL` … the build-stamp) must remain private, and the existing public methods after the ctors must stay public. The cleanest edit is to insert the `LoadedData` struct into the existing `public:` section just above where the json ctor was, and the delegating ctor into the existing `private:` member section (or a private section adjacent to it). Match the surrounding access-specifier structure rather than blindly pasting the labels above; the goal is: `LoadedData` public, `SuperAncillary(LoadedData&&)` private.

- [ ] **Step 5: Turn the string ctor into a declaration (L959)**

Replace:
```cpp
    SuperAncillary(const std::string& s) : SuperAncillary(nlohmann::json::parse(s)) {};
```
with a declaration only (keep its doc comment):
```cpp
    explicit SuperAncillary(const std::string& s);
```

- [ ] **Step 6: Create `src/superancillary.cpp`**

```cpp
#include "CoolProp/superancillary/superancillary.h"
#include "nlohmann/json.hpp"

#include <string>
#include <vector>

namespace CoolProp {
namespace superancillary {
namespace {

template <typename ArrayType>
std::vector<ChebyshevExpansion<ArrayType>> load_expansions(const nlohmann::json& j, const std::string& key) {
    std::vector<ChebyshevExpansion<ArrayType>> buf;
    for (auto& block : j.at(key)) {
        buf.emplace_back(block.at("xmin"), block.at("xmax"), block.at("coef"));
    }
    return buf;
}

template <typename ArrayType>
typename SuperAncillary<ArrayType>::LoadedData parse_loaded(const std::string& s) {
    nlohmann::json j = nlohmann::json::parse(s);
    typename SuperAncillary<ArrayType>::LoadedData d;
    d.rhoL = load_expansions<ArrayType>(j, "jexpansions_rhoL");
    d.rhoV = load_expansions<ArrayType>(j, "jexpansions_rhoV");
    d.p = load_expansions<ArrayType>(j, "jexpansions_p");
    d.Tcrit_num = j.at("meta").at("Tcrittrue / K");
    d.rhocrit_num = j.at("meta").at("rhocrittrue / mol/m^3");
    if (j.contains("check_points")) {
        for (const auto& pt : j.at("check_points")) {
            d.check_points.push_back({pt.at("T / K").get<double>(), pt.at("p(mp) / Pa").get<double>(),
                                      pt.at("rho'(mp) / mol/m^3").get<double>(), pt.at("rho''(mp) / mol/m^3").get<double>(),
                                      pt.at("p(SA)/p(mp)").get<double>(), pt.at("rho'(SA)/rho'(mp)").get<double>(),
                                      pt.at("rho''(SA)/rho''(mp)").get<double>()});
        }
    }
    return d;
}

}  // namespace

template <typename ArrayType>
SuperAncillary<ArrayType>::SuperAncillary(const std::string& s) : SuperAncillary(parse_loaded<ArrayType>(s)) {}

// The one and only ArrayType used across the codebase. Construction of
// SuperAncillary<other types> from a string is intentionally not provided.
template SuperAncillary<std::vector<double>>::SuperAncillary(const std::string&);

}  // namespace superancillary
}  // namespace CoolProp
```

- [ ] **Step 7: Retarget the 3 tests to the string ctor**

In `src/Tests/CoolProp-Tests.cpp`:
- Line 3399: `superancillary::SuperAncillary<std::vector<double>> anc{json};` → `superancillary::SuperAncillary<std::vector<double>> anc{json.dump()};`
- Line 3456: same change (`anc{json}` → `anc{json.dump()}`).
- Line 3773: `superancillary::SuperAncillary<std::vector<double>> anc{jsuper};` → `superancillary::SuperAncillary<std::vector<double>> anc{jsuper.dump()};`

(These tests still parse/inspect json with nlohmann directly — that's fine, the test TU is non-installed. Only the *construction* changes to the string ctor.)

- [ ] **Step 8: Build**

Run: `cmake --build build_catch --target CatchTestRunner -j8`
Expected: clean build. The new `src/superancillary.cpp` is auto-globbed by `file(GLOB … src/*.cpp)`. If the linker complains about an undefined `SuperAncillary<std::vector<double>>::SuperAncillary(const std::string&)`, the explicit instantiation in Step 6 is missing/mismatched. If it complains about `SuperAncillary<OtherType>`, some other ArrayType is string-constructed somewhere — search for it and add a matching `template …` instantiation line (the spec expects only `std::vector<double>`).

- [ ] **Step 9: Header is nlohmann-free**

Run: `grep -nE "nlohmann|#include .*json" include/CoolProp/superancillary/superancillary.h`
Expected: no matches (doc-comment mentions of the word "JSON" in prose are acceptable; there must be no `#include` and no `nlohmann::` token).

- [ ] **Step 10: Behavior parity — superancillary suites**

Run: `./build_catch/CatchTestRunner "[superanc],[HSU_D_ptsweep],[SVDSBTL]" -s`
Expected: all pass — the retargeted `[superanc]` water tests, the H/S/U superancillary sweep, and the SVDSBTL tests (which use the fluid superancillary) load and evaluate to identical numbers. If `[HSU_D_ptsweep]` is hidden (`[.]`), run it explicitly: `./build_catch/CatchTestRunner "[HSU_D_ptsweep]"`.

- [ ] **Step 11: Symbol-leak gate stays clean**

Run: `cmake -B build_shared -S . -DCOOLPROP_SHARED_LIBRARY=ON -DCMAKE_BUILD_TYPE=Release >/tmp/sh.log 2>&1 && cmake --build build_shared -j8 >>/tmp/sh.log 2>&1; ./dev/ci/check-json-symbols.sh "$(find build_shared \( -name 'libCoolProp.so' -o -name 'libCoolProp.dylib' \) | head -1)"`
Expected: `OK: no symbols matching /nlohmann|valijson/`.

- [ ] **Step 12: Commit**

```bash
git add include/CoolProp/superancillary/superancillary.h src/superancillary.cpp src/Tests/CoolProp-Tests.cpp
git restore --staged .beads/issues.jsonl 2>/dev/null || true
git commit -m "refactor(json): de-leak nlohmann from installed superancillary.h (CoolProp-xa8w.4)

Move JSON construction out-of-line into src/superancillary.cpp (explicit
instantiation for std::vector<double>); the installed header drops the
nlohmann include and gains a plain LoadedData bundle + private delegating
ctor. Construction is now string-only; 3 tests use the string ctor.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```
Use `--no-verify` if the bd hook stages `.beads/issues.jsonl` (report it).

---

## Final Verification (before PR)

- [ ] **Step 1: Full fast suite (cross-backend regression)**

Run: `./build_catch/CatchTestRunner "~[slow]"`
Expected: 0 failures. Superancillary feeds Helmholtz saturation, SVDSBTL, region boundary curves, and HS routines, so the whole suite exercises it.

- [ ] **Step 2: Confirm no other installed header reintroduces nlohmann via superancillary**

Run: `grep -rn "nlohmann" include/CoolProp/superancillary/`
Expected: no matches anywhere under the superancillary include dir.

- [ ] **Step 3: Preflight**

Run: `./dev/ci/preflight.sh`
Expected: build + json-symbols + tests pass. clang-tidy diff-only: the new `src/superancillary.cpp` and the header edits — confirm no *new* gating findings on added lines (the cross-reference technique: only findings git-blamed to this branch's commit gate CI).

- [ ] **Step 4: Adversarial review (MANDATORY before `gh pr create`)**

Per CLAUDE.md, dispatch the review subagent vs `master`. Focus:
  - the delegating ctor's member-init reproduces the OLD json ctor field-by-field (rhoL/rhoV/p from the same keys; `m_Tmin = m_p.xmin()`; Tcrit/rhocrit from `meta`; `m_pmin`/`m_pmax` from `m_p.eval`; check_points loop identical) — a dropped/reordered field or a wrong `meta` key is a real bug;
  - init-list order matches member declaration order (no `-Wreorder`);
  - the explicit instantiation covers every `ArrayType` that is string-constructed anywhere (grep to be sure it's only `std::vector<double>`);
  - the lazy caloric builds (m_hL/m_sL/…) and `get_check_points()` still behave identically;
  - the header truly has zero nlohmann; `LoadedData` carries no json type.
Then push + `gh pr create`, then re-review any post-review commits.

- [ ] **Step 5: bd**

`bd update CoolProp-xa8w.4 --notes "PR #<n> open…"`. File a small follow-on issue (or note on `rxxx`/this) to **broaden `dev/ci/check-installed-headers.sh`** to "no installed header pulls nlohmann/valijson" once BOTH this and `CoolProp-rxxx` have merged (flips the TODO already in that script).

---

## Self-Review (plan vs spec)

**Spec coverage:** §3 header surgery → Task 1 Steps 2-5; §3 new `.cpp` + explicit instantiation → Step 6; §3 tests → Step 7; §6 verification (header-grep, parity, symbol gate, build) → Steps 9-11 + Final Verification; §7 guard-broadening follow-on → Final Verification Step 5.

**Placeholder scan:** none — all code is shown in full (the `LoadedData` struct, the delegating ctor with the exact init-list, the whole `.cpp`, the three exact test edits). No "TODO/handle X" steps.

**Type consistency:** `LoadedData` fields (`rhoL/rhoV/p` as `std::vector<ChebyshevExpansion<ArrayType>>`, `Tcrit_num`, `rhocrit_num`, `check_points`) match what `parse_loaded` populates and what the delegating ctor consumes (`std::move(d.rhoL)` → `m_rhoL`, `d.Tcrit_num` → `m_Tcrit_num`, etc.). The string ctor signature `SuperAncillary(const std::string&)` is identical in the header declaration (Step 5), the `.cpp` definition, and the explicit instantiation (Step 6). `load_expansions`/`parse_loaded` are defined before use in the `.cpp`.
