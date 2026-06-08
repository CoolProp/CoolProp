# De-leak `nlohmann::json` from installed `superancillary.h` — Design

**Date:** 2026-06-07
**Status:** Approved design, ready for implementation planning
**Epic:** `CoolProp-xa8w` (RapidJSON → nlohmann/json migration)
**Issue:** `CoolProp-xa8w.4` (blocks Phase Final `CoolProp-xa8w.3`)

---

## 1. Goal

Remove `#include "nlohmann/json.hpp"` (and every `nlohmann::json` mention) from the
installed header `include/CoolProp/superancillary/superancillary.h`, so a consumer
that includes `CoolProp.h` / `AbstractState.h` / `CoolPropFluid.h` (which transitively
reach `superancillary.h`) no longer needs `nlohmann/json.hpp` on its include path to
compile. This is a **compile-time / header-ABI** leak, not a symbol-export leak — the
symbol-leak gate is already green (no nlohmann symbol is exported), because the only
in-`libCoolProp` construction path (`CoolPropFluid.h::get_superanc()`) uses the
`std::string` ctor whose signature carries no nlohmann, and the `nlohmann::json` ctor
is instantiated only in the test binary. After this change, with `detail/json.h`
already install-excluded (`CoolProp-rxxx`), **no installed header pulls nlohmann/valijson.**

## 2. Key facts (verified)

- `superancillary.h` is **header-only** and its classes are **templates**:
  `ChebyshevExpansion<ArrayType>`, `ChebyshevApproximation1D<ArrayType>`,
  `SuperAncillary<ArrayType>`.
- **`SuperAncillary<ArrayType>` is only ever instantiated with `std::vector<double>`** —
  `CoolPropFluid.h` (`SuperAncillary_t`), `SVDSBTLBackend`, both `region/*BoundaryCurve.h`,
  and the tests all use `std::vector<double>`. (`MeltingCaloric` uses
  `ChebyshevApproximation1D<Eigen::ArrayXd>` directly, but **not** `SuperAncillary`.) So the
  explicit-instantiation set for `SuperAncillary` construction is exactly **one** type.
- **nlohmann reaches only into `SuperAncillary`'s construction**, not the evaluation code
  or the lower templates: `ChebyshevExpansion`'s ctor is `(double xmin, double xmax,
  const ArrayType& coeff)` — `ArrayType`, not json (the loader relies on nlohmann's
  implicit `json → ArrayType` conversion of `block.at("coef")`).
- The nlohmann surface in the header is exactly: the `#include` (line 46), the private
  `loader(const nlohmann::json&, key)` (883-897), the public
  `SuperAncillary(const nlohmann::json&)` ctor (938-954), and the inline
  `SuperAncillary(const std::string& s) : SuperAncillary(nlohmann::json::parse(s)) {}` (959).
- **`ChebyshevApproximation1D` is non-assignable** — it has a `const double thresh_imag`
  member (line 480). So `m_rhoL`/`m_rhoV`/`m_p` (of that type) must be **member-initialized
  (constructed)**, never assigned. This rules out an "assign-in-body" out-of-line ctor and
  dictates the delegating-ctor form below.
- Construction is reached **inline** from `CoolPropFluid.h::get_superanc()` (line 433-438,
  inline in an installed header) via `make_shared<SuperAncillary_t>(superancillaries_str)`.
- **Direct `nlohmann::json` construction exists only in 3 tests** (`CoolProp-Tests.cpp:3399,
  3456, 3773 — `SuperAncillary<…> anc{json}`). Production uses the string ctor.
  region/SVDSBTL/MeltingCaloric only *use* `SuperAncillary` (methods / shared_ptr), never
  construct it from json.

## 3. Design — explicit instantiation of the out-of-line string ctor

Keep every templated, perf-critical evaluation path inline in the header; move only the
JSON-consuming construction out-of-line into one `.cpp` compiled into `libCoolProp`, and
explicitly instantiate it for `std::vector<double>`.

**`include/CoolProp/superancillary/superancillary.h` (installed):**
- Remove `#include "nlohmann/json.hpp"` (line 46).
- Remove the private `loader(const nlohmann::json&, …)` member (883-897).
- Remove the public `SuperAncillary(const nlohmann::json&)` ctor (938-954).
- Add a **public nested POD `struct LoadedData`** (templated, nlohmann-free) holding the
  already-parsed pieces the ctor needs:
  ```cpp
  struct LoadedData {
      std::vector<ChebyshevExpansion<ArrayType>> rhoL, rhoV, p;  // loader outputs
      double Tcrit_num = 0, rhocrit_num = 0;
      std::vector<CheckPoint> check_points;
  };
  ```
- Add a **private delegating ctor** that member-initializes from it (inline, nlohmann-free —
  this preserves the exact member-init semantics of the old json ctor, including the derived
  `m_Tmin`/`m_pmin`/`m_pmax` from `m_p`):
  ```cpp
  explicit SuperAncillary(LoadedData&& d)
    : m_rhoL(std::move(d.rhoL)), m_rhoV(std::move(d.rhoV)), m_p(std::move(d.p)),
      m_Tmin(m_p.xmin()), m_Tcrit_num(d.Tcrit_num), m_rhocrit_num(d.rhocrit_num),
      m_pmin(m_p.eval(m_p.xmin())), m_pmax(m_p.eval(m_p.xmax())),
      m_check_points(std::move(d.check_points)) {}
  ```
  (Note: `m_rhoL` is `ChebyshevApproximation1D<ArrayType>` constructed from
  `std::vector<ChebyshevExpansion<ArrayType>>&&` — the same converting construction the old
  json ctor used via `loader(...)`.)
- Replace the inline string ctor (959) with a **declaration only**:
  `explicit SuperAncillary(const std::string& s);`.

**New `src/superancillary.cpp` (non-installed; picked up by `file(GLOB src/*.cpp)` at
CMakeLists.txt:297 — no CMake change):**
```cpp
#include "CoolProp/superancillary/superancillary.h"
#include "nlohmann/json.hpp"
namespace CoolProp { namespace superancillary {
namespace {
  template<typename ArrayType>
  std::vector<ChebyshevExpansion<ArrayType>> load_expansions(const nlohmann::json& j, const std::string& key) {
      std::vector<ChebyshevExpansion<ArrayType>> buf;
      for (auto& block : j.at(key)) {
          buf.emplace_back(block.at("xmin"), block.at("xmax"), block.at("coef"));
      }
      return buf;
  }
  template<typename ArrayType>
  typename SuperAncillary<ArrayType>::LoadedData parse_loaded(const std::string& s) {
      nlohmann::json j = nlohmann::json::parse(s);
      typename SuperAncillary<ArrayType>::LoadedData d;
      d.rhoL = load_expansions<ArrayType>(j, "jexpansions_rhoL");
      d.rhoV = load_expansions<ArrayType>(j, "jexpansions_rhoV");
      d.p    = load_expansions<ArrayType>(j, "jexpansions_p");
      d.Tcrit_num   = j.at("meta").at("Tcrittrue / K");
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
}
template<typename ArrayType>
SuperAncillary<ArrayType>::SuperAncillary(const std::string& s) : SuperAncillary(parse_loaded<ArrayType>(s)) {}
// The one and only instantiation actually used across the codebase.
template SuperAncillary<std::vector<double>>::SuperAncillary(const std::string&);
}}
```
The string ctor delegates to the private `LoadedData&&` ctor, so member-init (required by the
non-assignable `ChebyshevApproximation1D`) is preserved.

**Tests (`src/Tests/CoolProp-Tests.cpp:3399, 3456, 3773`):** change `SuperAncillary<…>
anc{json}` → `anc{json.dump()}` (string ctor). The tests already hold the json; the re-parse
cost is irrelevant in a test.

## 4. What stays the same

- All templated **evaluation** code (`eval_sat`, `get_T_from_p`, lazy caloric builds,
  `ChebyshevExpansion`/`ChebyshevApproximation1D`) stays inline in the header — **no
  de-inlining of the hot path.**
- `CoolPropFluid.h::get_superanc()` stays inline and unchanged — it calls the `std::string`
  ctor, now resolved against the explicitly-instantiated symbol in `libCoolProp`.
- The symbol-leak gate stays green (the out-of-line ctor's signature is `const std::string&`;
  its nlohmann use is internal to `superancillary.cpp` and inlined/static).

## 5. Trade-off

`SuperAncillary` construction is now **closed to `std::vector<double>`** (the only type ever
used). Constructing `SuperAncillary<SomeOtherType>` from a string would fail to link unless a
matching explicit instantiation is added. Acceptable — it is a CoolProp-internal performance
class, never instantiated by consumers with custom array types.

## 6. Testing & verification

- **Header is nlohmann-free:** `grep -n "nlohmann\|json" include/CoolProp/superancillary/superancillary.h`
  → no `#include`, no type usage (doc-comment mentions of "JSON format" are fine).
- **Behavior parity (the safety net):** `[HSU_D_ptsweep]` (superancillary H/S/U sweep),
  `[SVDSBTL]`, `[HS]`/`[HS_dbg]`, and the saturation tests must stay green — every fluid's
  superancillary must load and evaluate to identical numbers. The 3 retargeted tests must pass.
- **Symbol gate stays clean:** `./dev/ci/check-json-symbols.sh` on `build_shared` → still 0.
- **Build:** the new `src/superancillary.cpp` compiles and links; `make_shared<SuperAncillary_t>(str)`
  in consumer TUs links against its instantiated ctor (no nlohmann needed to compile a consumer).
- **Preflight + adversarial review** (focus: the delegating-ctor member-init reproduces the old
  json ctor's field-by-field semantics; the explicit instantiation covers every used `ArrayType`;
  no behavior change in the lazy caloric builds).

## 7. Out of scope / follow-on

- **Broaden `dev/ci/check-installed-headers.sh`** (the `rxxx` guard) from "`detail/json.h` not
  shipped" to "**no installed header pulls nlohmann/valijson**" — this becomes *true* only once
  BOTH `rxxx` (excludes `detail/json.h`) and this change (de-leaks `superancillary.h`) have
  merged. File/track as a small follow-on that lands after both, flipping the `TODO` already left
  in that script.
- `detail/rapidjson.h` still ships (rapidjson, not nlohmann) until Phase Final (`xa8w.3`).
- The broader consumer-compile gate (coordinated with the nanobind interface) remains deferred.
