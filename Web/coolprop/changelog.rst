Changelog for CoolProp
======================

8.0.0
-----

Highlights:

* Added the :doc:`SVDSBTL </coolprop/SVDSBTL>` SVD-compressed tabular lookup backend (factory string ``SVDSBTL&<source>``, where source ∈ {``HEOS``, ``REFPROP``, ``IF97``}).  Combines a region atlas with per-region SVD compression plus a critical-patch fallback to the source backend; produces sub-microsecond per-probe property evaluation in the batched ``fast_evaluate`` path at IAPWS-G13-15 ``T(p, h)`` conformance for water and a single-digit-percent accuracy ceiling for the multi-fluid HEOS-backed presets.  Disk footprint ~7-14 MB per (fluid, input pair, source backend), cached under ``~/.CoolProp/SVDTables/``.  Off by default in ``PropsSI`` — see ``ALLOW_SVDSBTL_IN_PROPSSI``.  See PRs `#2917`, `#2938-#2940`, `#2944-#2957`.
* Added mass-basis vapor quality (``Qmass``) support across HEOS and REFPROP backends, paralleling the existing molar quality (``Q``). Adds a new ``iQmass`` keyed parameter, a ``Qmass()`` accessor on ``AbstractState``, and 8 new ``Qmass``-bearing input pairs (``QmassT_INPUTS``, ``PQmass_INPUTS``, ``QmassSmolar_INPUTS``, ``QmassSmass_INPUTS``, ``HmolarQmass_INPUTS``, ``HmassQmass_INPUTS``, ``DmolarQmass_INPUTS``, ``DmassQmass_INPUTS``). Mixtures supported from day one — REFPROP uses its native ``kq=2`` flag in ``TQFLSHdll``/``PQFLSHdll`` for ``QmassT``/``PQmass``; the other 6 pairs and all HEOS pairs use a TOMS748 root-find on ``Qmolar`` (typically 5–8 iterations).
* Added Chebyshev superancillary functions for the SRK and Peng-Robinson cubic EOS, mirroring the multiparameter-EOS superancillaries from v7 — fast, robust saturation curves for the cubic backends.  The SRK/PR ``Omega`` coefficients are now stored at full precision rather than the historical rounded values.  See PRs `#2744`, `#2745`.
* Large expansion of the fluid and mixture libraries: new / updated pure-fluid EOS for Chlorine (Cl₂), R-1336mzz(Z), R-1130(E), R-1123, R-1224yd(Z), R-1132(E), THF, propylene glycol, vinyl chloride, the R-1233zd(E) international-standard EOS, perfluoro-n-butane / -pentane / -hexane, and 3rd-generation R-1243zf / R-1234yf EOS.  Many new mixture interaction pairs and predefined mixtures, including the Bell-JCED-2025, Bell-IJT-2020 and Bell-JPCRD-2022/2023 models, NIST IR 8570 pairs, JPCRD reference refrigerant mixtures, and the ASHRAE Standard 34 (2026) predefined blends.
* REFPROP backend: saturated-state properties and saturation / two-phase derivatives are now exposed (matching the HEOS surface); imposed-phase ``DmolarSmolar`` and the ``DmassQ`` / ``DmolarQ`` input pairs are wired through the corresponding REFPROP flash routines.  See PRs `#2990`, `#2865`, `#2866`.
* Thread-safety hardening for the native (HEOS / incompressible / HumidAir) code paths: fluid-library static initialization is now race-free, ``HumidAirProp`` uses per-thread Water/Air backends, and shared derivative counters are atomic.  The REFPROP backend remains **not** thread-safe — REFPROP itself is not reentrant, so calls into the ``REFPROP`` backend must still be serialized by the caller.  See PRs `#2800`, `#2831`, `#2855`.
* Tabular backends can now be evaluated directly (cache-bypassing batch ``fast_evaluate``), with new ``TABULAR_NX`` / ``TABULAR_NY`` configuration keys for grid resolution and a fix for a BICUBIC ``PT`` segfault near the saturation curve.  See PRs `#2920`, `#2894`, `#2891`.
* Numerous solver-robustness and graceful-error-handling fixes across the flash and density solvers (illegal quality inputs, zero-/one-length ``PropsSI`` arrays, incompressible molar requests, sub-``TminPsat`` saturation, ancillaries above the reducing temperature, and more).
* New native desktop GUI built with Tauri + React (`#2715`), plus a much-expanded Mathcad wrapper and interactive 3D molecule viewers on the fluid documentation pages.
* Build-system modernization: git submodules replaced by CPM.cmake (boost fetched as a trimmed subset from ``CoolProp/boost-headers``), Eigen bumped to 5.0.1, and a broad C++17 cleanup that also cuts compile times.
* Repository-wide code-quality program: enforced ``clang-format`` (pre-commit + CI), diff-only ``clang-tidy``, ``cppcheck``, CodeQL, include-what-you-use, and a single-script ``dev/ci/preflight.sh`` pre-push gate.
* **RapidJSON removed; JSON now on nlohmann/json + Valijson.**  All JSON parsing and JSON-Schema validation were migrated off RapidJSON to `nlohmann/json <https://github.com/nlohmann/json>`_ (with CBOR support) and `Valijson <https://github.com/tristanpenman/valijson>`_, and RapidJSON was deleted from the tree.  Crucially, **no nlohmann/valijson symbols are exported from any shipped binary and no JSON-library type appears in any installed public header**, so the JSON libraries are an internal implementation detail that cannot ODR-clash with a downstream consumer's own copy.  The symbol hiding is enforced at *link* time per shared product (ELF ``--version-script`` / Mach-O ``-unexported_symbols_list``), with CI gates over ``libCoolProp`` and the SWIG wrappers plus an installed-header hygiene check.  See PRs `#3094`, `#3100`, `#3104`, `#3106`, `#3108`, `#3111`, `#3112`, `#3116`, `#3117`.
* **Python wrapper rebuilt on nanobind** — a modern, lower-maintenance binding that ships as a single stable-ABI (``abi3``) wheel per platform on Python ≥ 3.12 (one wheel instead of one per Python version).  The public import tree is preserved, so it is a drop-in for almost all code; the legacy non-SI ``Props`` / ``HAProps`` are removed (use ``PropsSI`` / ``HAPropsSI``).  See the behavior-changes note below.

**Behavior changes (potentially breaking):**

* **WASM / JavaScript wrapper:** composition vectors and phase-envelope
  data now use **native JavaScript arrays** instead of the embind
  ``VectorDouble`` wrapper class. Setters take a JS array directly
  (``AS.set_mole_fractions([0.4, 0.6])``); getters return one
  (``var z = AS.get_mole_fractions(); console.log(z[0]);``);
  ``get_phase_envelope_data()`` returns a plain object whose fields
  (``T``, ``p``, ...) are JS arrays (``data.T[i]``). The
  ``VectorDouble``/``VectorString`` classes have been removed from
  the exposed surface.

* **SWIG-based language wrappers (C#, VB.NET, Java, R, PHP):** the
  generated shared library is now language-prefixed —
  ``CoolPropCsharp`` (C# and VB.NET), ``CoolPropJava``, ``CoolPropR``,
  ``CoolPropPHP`` — instead of ``CoolProp``. This removes the long-
  standing name collision with the conventional CoolProp shared
  library. The CoolProp API (classes, namespaces, functions) is
  unchanged; only the binary filename and the loader reference change.
  Update your ``DllImport`` (C#/VB.NET), ``System.loadLibrary`` (Java),
  or ``dyn.load`` (R) calls to the new name. See GitHub
  `#1674 <https://github.com/CoolProp/CoolProp/issues/1674>`_.

* **Python wrapper: the interface is now built on nanobind.** The hand-written
  Cython ``AbstractState`` interface has been replaced by a `nanobind
  <https://github.com/wjakob/nanobind>`_ core, and the parallel (and
  incomplete) pybind11 interface has been removed. The **public import tree is
  preserved** — ``CoolProp.CoolProp`` (now the nanobind core),
  ``CoolProp.AbstractState``, ``CoolProp.HumidAirProp``, ``CoolProp.State``,
  ``CoolProp.Plots``, ``CoolProp.BibtexParser`` and the rest resolve at exactly
  the same paths — so the large majority of user code needs no change.
  Downstream Cython that ``cimport``\ s CoolProp's ``State`` /
  ``AbstractState`` / ``constants_header`` surface (e.g. PDSim) also keeps
  working: a link-free, capsule-forwarding ``State`` shim preserves the
  cimportable contract. On Python ≥ 3.12 the wrapper ships as a single
  stable-ABI (``abi3``) wheel per platform instead of one wheel per Python
  version (Python 3.9–3.11 get per-version wheels, below nanobind's
  stable-ABI floor).

  ``AbstractState`` is now a genuine Python *type*: ``AbstractState("HEOS",
  "Water")`` still constructs a state, but ``isinstance(x, AbstractState)``
  now works too (previously the factory was a function and the check raised).

  The following **removals are breaking** (all long-deprecated or never part of
  the intended API):

  - the non-SI free function ``Props`` (deprecated behind a warning for years) —
    use ``PropsSI`` (SI units). The ``State`` class keeps its ``Props(key)``
    *method* (a keyed-output accessor); only the free function is gone.
  - the non-SI humid-air function ``HAProps`` — use ``HAPropsSI`` (SI units).
    ``CoolProp.HumidAirProp.HAProps`` now raises with a message pointing at
    ``HAPropsSI``; ``HAPropsSI``, ``HAProps_Aux`` and ``cair_sat`` are unchanged.
  - ``re_split`` (an accidental re-export of Python's ``re.split``) and
    ``rebuildState`` (an internal ``State``-unpickling helper) are no longer
    exported from ``CoolProp.CoolProp``.

* Default ``update`` for ``HmolarQ_INPUTS``, ``QSmolar_INPUTS``,
  ``DmolarQ_INPUTS`` (and their mass-input equivalents) on pure fluids
  now raises :cpapi:`CoolProp::MultipleSolutionsError` when the input
  value admits more than one temperature on the saturation curve. This
  affects water saturated-vapor enthalpy near the 540 K peak, water /
  D2O saturated-liquid density near the 4 °C / 11 °C maximum, and any
  similar non-monotonic input. Single-root inputs continue to work as
  before. To select a branch, use ``update_with_guesses`` with a
  ``guess.T`` set to the target temperature region. See GitHub
  `#2773 <https://github.com/CoolProp/CoolProp/issues/2773>`_.

* New input-pair dispatches in ``update_with_guesses``:
  ``HmolarQ_INPUTS``, ``HmassQ_INPUTS``, ``QSmolar_INPUTS``,
  ``QSmass_INPUTS``, ``DmolarQ_INPUTS``, ``DmassQ_INPUTS`` —
  use ``GuessesStructure.T`` to pick a branch. See
  ``LowLevelAPI.rst`` § *Disambiguating Multiple Saturation Roots*.

* New exception type :cpapi:`CoolProp::MultipleSolutionsError`
  (subclass of :cpapi:`CoolProp::ValueError`).

* **Multicomponent PT-flash results and** ``phase()`` **output change for
  mixtures (#3026).** The blind PT flash for mixtures (Cubic and HEOS
  backends) now runs a rigorous Michelsen TPD stability analysis followed by
  a Successive-Substitution / second-order Gibbs phase split, replacing the
  previous more limited blind-flash handling that could misclassify two-phase
  mixture states as single-phase. This changes computed compositions,
  densities and the reported phase for some mixture states. A new ``MIXTURE_STABILITY_ALGORITHM``
  config key selects the algorithm (``1`` = Michelsen, the new default;
  ``0`` = the legacy Gernert et al. 2014 path).

* **He–Ar mixture results change (#3118).** Helium–argon (and the new
  helium–neon and neon–argon) mixtures now use the Tkaczuk et al. (2020)
  reducing parameters and departure function instead of the previous
  ideal-mixing assumption.

* **Public C++ headers reorganized into an** ``include/CoolProp/`` **tier
  tree (GH #1280).** [BREAKING for C++ consumers, with back-compat shims.]
  The previously flat ``include/`` layout has been restructured into a
  namespaced tree rooted at ``include/CoolProp/``. The canonical include
  path for the primary API is now ``"CoolProp/CoolProp.h"`` (and
  ``"CoolProp/CoolPropLib.h"``, ``"CoolProp/AbstractState.h"``,
  ``"CoolProp/HumidAirProp.h"``, etc.). The tiers are: the top-level
  ``CoolProp/`` primary API; ``CoolProp/detail/`` (internal support);
  ``CoolProp/fluids/``, ``CoolProp/numerics/``, ``CoolProp/region/``,
  ``CoolProp/sbtl/``, ``CoolProp/svd/``, ``CoolProp/schemas/``,
  ``CoolProp/superancillary/``, ``CoolProp/plotting/`` and
  ``CoolProp/Backends/SVDSBTL/``.

  **Back-compat shims.** Every old flat path (``CoolProp.h``,
  ``AbstractState.h``, ``HumidAirProp.h``, ``MatrixMath.h``, …) still exists,
  but is now a one-line forwarding shim that ``#include``\ s the new location
  and emits a ``#pragma message`` on first include pointing at the canonical
  path. Existing C++ that includes the old paths keeps compiling unchanged; it
  just prints a deprecation notice. To silence the notices, compile with
  ``-DCOOLPROP_NO_DEPRECATED_HEADER_WARNINGS``. **The flat shims are slated for
  removal at v9** — migrate your includes before then.

  **Six headers were deleted outright** (no shim): ``miniz.h``,
  ``rapidjson_include.h``, ``crossplatform_shared_ptr.h``, ``Tests.h``,
  ``SpeedTest.h`` and ``TestObjects.h`` — internal third-party / compatibility
  wrappers with no public consumer (RapidJSON itself was removed) and test-only
  scaffolding that was never part of the shipped API.

  **What downstream C++ consumers must do.** Switch your includes from
  ``"CoolProp.h"`` (and friends) to ``"CoolProp/CoolProp.h"``. The install
  rules now ship the full ``include/CoolProp/`` tree (plus the flat shims), so
  a manual include of that directory gives usable include paths out of the box.
  Transitive dependencies exposed by the public headers: the ``fluids/``,
  ``numerics/`` and ``superancillary/`` tiers pull in **Eigen**, so Eigen must
  be on your include path if you include those tiers; and
  ``CoolProp/detail/strings.h`` includes **fmt** unless you compile with
  ``-DNO_FMTLIB``. (Boost is *not* required: the superancillary rootfinder is
  defined out-of-line in the compiled library.) (The generated database headers —
  ``*_JSON.h`` / ``*_CBOR.h`` — and the internal ``detail/json.h`` /
  ``detail/msgpack.h`` are intentionally *not* installed.) See PRs #3074,
  #3075 and #3078 (the #1280 phases), plus #3082
  (the related rapidjson / ``<filesystem>`` public-surface de-leak).

* **Build-system requirements (packagers).** Git submodules have been replaced
  by CPM.cmake, so a fresh configure now needs **network access** to fetch
  dependencies unless ``CPM_SOURCE_CACHE`` points at a populated cache. The
  minimum CMake is now **3.14** and a **C++17** compiler is required.

Issues closed:

* `#2254 <https://github.com/CoolProp/CoolProp/issues/2254>`_ : CSharp wrapper folder missing from download location since version 6.4.2
* `#1674 <https://github.com/CoolProp/CoolProp/issues/1674>`_ : Make the SWIG-generated DLL have a language prefix
* `#2326 <https://github.com/CoolProp/CoolProp/issues/2326>`_ : EntryPointNotFound after upgrade to newest version
* `#2189 <https://github.com/CoolProp/CoolProp/issues/2189>`_ : CoolProp C++ High level interface
* `#2193 <https://github.com/CoolProp/CoolProp/issues/2193>`_ : Handle PowerPC builds transparently
* `#2194 <https://github.com/CoolProp/CoolProp/issues/2194>`_ : Fix the builds for Windows on ARM and ARM64
* `#2195 <https://github.com/CoolProp/CoolProp/issues/2195>`_ : Wrong enthalpy outputs for illegal quality inputs
* `#2209 <https://github.com/CoolProp/CoolProp/issues/2209>`_ : Problem with zero Psat for incompressible mixture
* `#2215 <https://github.com/CoolProp/CoolProp/issues/2215>`_ : request  for r513a
* `#2217 <https://github.com/CoolProp/CoolProp/issues/2217>`_ : CoolProp SimpleCompressionCycle not close due to numerical tolerance(?)
* `#2226 <https://github.com/CoolProp/CoolProp/issues/2226>`_ : first_two_phase_deriv function is not implemented in CoolpropLib.h.
* `#2244 <https://github.com/CoolProp/CoolProp/issues/2244>`_ : PropsSI failed ungracefully
* `#2255 <https://github.com/CoolProp/CoolProp/issues/2255>`_ : HAPropsSI("Twb","T",200+273.15,"W",0.2,"P",1000E3) gives wrong result
* `#2308 <https://github.com/CoolProp/CoolProp/issues/2308>`_ : Issue calculating molar fractions of every component out the biphase zone for a custom fluid modelled with AbstractState
* `#2332 <https://github.com/CoolProp/CoolProp/issues/2332>`_ : CO2  at 297.45 K
* `#2339 <https://github.com/CoolProp/CoolProp/issues/2339>`_ : Comment creer diagramme des phases avec Excel ?
* `#2360 <https://github.com/CoolProp/CoolProp/issues/2360>`_ : Negative vapor and unsteady liquid quality caluculation for Helium
* `#2379 <https://github.com/CoolProp/CoolProp/issues/2379>`_ : Mixtures in C#
* `#2380 <https://github.com/CoolProp/CoolProp/issues/2380>`_ : Cannot use CoolProp in Excel anymore because add-in is not signed
* `#2385 <https://github.com/CoolProp/CoolProp/issues/2385>`_ : R1336mzz(Z) fluid
* `#2391 <https://github.com/CoolProp/CoolProp/issues/2391>`_ : CoolProp.dll has no versionnumber
* `#2417 <https://github.com/CoolProp/CoolProp/issues/2417>`_ : PropsSI array argument functionality gives unexpected results for arrays of length 1 or 0
* `#2427 <https://github.com/CoolProp/CoolProp/issues/2427>`_ : web site documentation for mixtures
* `#2433 <https://github.com/CoolProp/CoolProp/issues/2433>`_ : Mass fraction of saturated lithium bromide solution calculation question
* `#2434 <https://github.com/CoolProp/CoolProp/issues/2434>`_ : Discrepancy Between CoolProp DLL and Online Calculator Outputs
* `#2461 <https://github.com/CoolProp/CoolProp/issues/2461>`_ : Problem with displaying enthalpy diagram
* `#2533 <https://github.com/CoolProp/CoolProp/issues/2533>`_ : Official support for conda-forge
* `#2586 <https://github.com/CoolProp/CoolProp/issues/2586>`_ : Consistency issue when calculating enthalpy from entropy and density at dew line
* `#2602 <https://github.com/CoolProp/CoolProp/issues/2602>`_ : _crit needs to be synchronized with the superancillaries at load
* `#2654 <https://github.com/CoolProp/CoolProp/issues/2654>`_ : Github Release Tags missing?
* `#2655 <https://github.com/CoolProp/CoolProp/issues/2655>`_ : Wrappers need update.
* `#2658 <https://github.com/CoolProp/CoolProp/issues/2658>`_ : Data of Molar volume of saturated liquid water or ice (vbar_ws) [m^3/mol_H2O]
* `#2660 <https://github.com/CoolProp/CoolProp/issues/2660>`_ : Compatibility of CoolProp with MacBook Pro M4 chip
* `#2663 <https://github.com/CoolProp/CoolProp/issues/2663>`_ : Fluid mixtures not working for matlab wrapper
* `#2667 <https://github.com/CoolProp/CoolProp/issues/2667>`_ : MATLAB Wrapper
* `#2668 <https://github.com/CoolProp/CoolProp/issues/2668>`_ : Multithreaded Support with REFPROP Backend
* `#2670 <https://github.com/CoolProp/CoolProp/issues/2670>`_ : HAPropsSI with T_dp + R inputs returns 0 in v7.2.0 (worked in v6.3.0)
* `#2672 <https://github.com/CoolProp/CoolProp/issues/2672>`_ : REFPROP backend: mixture R1233zd(E) + Air cannot be loaded in CoolProp, while direct REFPROP calls work (with warning)
* `#2673 <https://github.com/CoolProp/CoolProp/issues/2673>`_ : DmassT and DmolarT inputs with Cubic backend not working in v7.2
* `#2680 <https://github.com/CoolProp/CoolProp/issues/2680>`_ : Address performance regression with superancillary functions
* `#2684 <https://github.com/CoolProp/CoolProp/issues/2684>`_ : DLL load failed while importing CoolProp
* `#2685 <https://github.com/CoolProp/CoolProp/issues/2685>`_ : CoolProp fails ungracefully when computing entropy through density+internal energy
* `#2691 <https://github.com/CoolProp/CoolProp/issues/2691>`_ : vcpkg port of coolprop out of date (6.4.3#3)
* `#2693 <https://github.com/CoolProp/CoolProp/issues/2693>`_ : PropsSImulti does not return _HUGE and does not trap Value Errors
* `#2694 <https://github.com/CoolProp/CoolProp/issues/2694>`_ : Download of REFPROP fails in CI build due to bad session key
* `#2696 <https://github.com/CoolProp/CoolProp/issues/2696>`_ : Issue template not being used since Nov. 2024
* `#2698 <https://github.com/CoolProp/CoolProp/issues/2698>`_ : Speed of sound shows numerical jump near vapor saturation line in two-phase HFRS model
* `#2701 <https://github.com/CoolProp/CoolProp/issues/2701>`_ : Equations Not Displaying on Ideal-Gas page of CoolProp.org
* `#2703 <https://github.com/CoolProp/CoolProp/issues/2703>`_ : HAPropsSI Patch (#2697) Created a Regression in the Docs/Workflows
* `#2711 <https://github.com/CoolProp/CoolProp/issues/2711>`_ : [ISSUE] Add table to docs with predefined mixtures
* `#2712 <https://github.com/CoolProp/CoolProp/issues/2712>`_ : [ISSUE] Add predefined mixtures from ASHRAE 2026 standard 34
* `#2714 <https://github.com/CoolProp/CoolProp/issues/2714>`_ : [ISSUE] Add Chlorine EOS
* `#2717 <https://github.com/CoolProp/CoolProp/issues/2717>`_ : [ISSUE] Get coverity to dump its output during build
* `#2718 <https://github.com/CoolProp/CoolProp/issues/2718>`_ : [ISSUE] Further speedup to caching
* `#2721 <https://github.com/CoolProp/CoolProp/issues/2721>`_ : [ISSUE] Implement reference refrigerant mixture models from JPCRD
* `#2724 <https://github.com/CoolProp/CoolProp/issues/2724>`_ : [ISSUE] Clang-format CI Misbehaving (silently) on Remote PRs
* `#2726 <https://github.com/CoolProp/CoolProp/issues/2726>`_ : [FEAT] Add even more mixture models
* `#2727 <https://github.com/CoolProp/CoolProp/issues/2727>`_ : Generate SUPERANCILLARY for R1234yf Lemmon-IJT-2022 EOS
* `#2738 <https://github.com/CoolProp/CoolProp/issues/2738>`_ : Critical density CO2
* `#2739 <https://github.com/CoolProp/CoolProp/issues/2739>`_ : [FEAT] Add superancillaries for cubic EOS
* `#2740 <https://github.com/CoolProp/CoolProp/issues/2740>`_ : [ISSUE] Mathematica won't build under new CMake Package Manager (CPM)
* `#2742 <https://github.com/CoolProp/CoolProp/issues/2742>`_ : Update cubic EOS coefficients
* `#2751 <https://github.com/CoolProp/CoolProp/issues/2751>`_ : Enforce clang-format
* `#2754 <https://github.com/CoolProp/CoolProp/issues/2754>`_ : Investigate the inclusion of teqp in CoolProp
* `#2755 <https://github.com/CoolProp/CoolProp/issues/2755>`_ : [ISSUE] Chlorine consistency plots are missing from docs
* `#2756 <https://github.com/CoolProp/CoolProp/issues/2756>`_ : [ISSUE] Wrong number of fluids in docs
* `#2762 <https://github.com/CoolProp/CoolProp/issues/2762>`_ : Tracking: bump fastchebpure pin after Akasaka/Lemmon EOS batch merges
* `#2763 <https://github.com/CoolProp/CoolProp/issues/2763>`_ : THF EOS (Fiedler et al. 2023): paper/FLD inconsistencies
* `#2764 <https://github.com/CoolProp/CoolProp/issues/2764>`_ : Propylene Glycol EOS (Eisenbach et al. 2021): critical-region instability
* `#2765 <https://github.com/CoolProp/CoolProp/issues/2765>`_ : R-1123 EOS (Akasaka et al., IJR 2020): Table 6 typo in Gaussian n_{14}
* `#2767 <https://github.com/CoolProp/CoolProp/issues/2767>`_ : Pure and mixture models in DOE report
* `#2769 <https://github.com/CoolProp/CoolProp/issues/2769>`_ : [ISSUE] CL2 does not precisely match REFPROP implementation
* `#2770 <https://github.com/CoolProp/CoolProp/issues/2770>`_ : HmolarQ calculation fails and rends AbstractState unstable
* `#2772 <https://github.com/CoolProp/CoolProp/issues/2772>`_ : Add at least one extended precision point for superancillaries to JSON
* `#2773 <https://github.com/CoolProp/CoolProp/issues/2773>`_ : [FEAT] Add means to specify which solution you want when multiple solutions are possible
* `#2774 <https://github.com/CoolProp/CoolProp/issues/2774>`_ : Move boost deps from archive to new repo
* `#2777 <https://github.com/CoolProp/CoolProp/issues/2777>`_ : fastchebpure: emit source_eos_hash in output/<fluid>_exps.json so CoolProp can verify SA freshness without git archaeology
* `#2779 <https://github.com/CoolProp/CoolProp/issues/2779>`_ : [ISSUE] eos_hash
* `#2787 <https://github.com/CoolProp/CoolProp/issues/2787>`_ : Fluid-library static init is not thread-safe (race during cold cache population)
* `#2815 <https://github.com/CoolProp/CoolProp/issues/2815>`_ : Dependabot should do actions too
* `#2825 <https://github.com/CoolProp/CoolProp/issues/2825>`_ : [ISSUE] CoolProp GUI: blank screen
* `#2828 <https://github.com/CoolProp/CoolProp/issues/2828>`_ : Windows Installer build: hashes.json decoding error from concurrent writes
* `#2830 <https://github.com/CoolProp/CoolProp/issues/2830>`_ : Doc builds are sometimes getting stuck in infinite loop
* `#2834 <https://github.com/CoolProp/CoolProp/issues/2834>`_ : Consider treating multi-valued saturation regions as an error in update() (follow-up to #2773)
* `#2869 <https://github.com/CoolProp/CoolProp/issues/2869>`_ : refactor: modernize-avoid-c-arrays — replace C arrays with std::array
* `#2870 <https://github.com/CoolProp/CoolProp/issues/2870>`_ : refactor: cppcoreguidelines-pro-type-vararg — replace remaining printf calls
* `#2871 <https://github.com/CoolProp/CoolProp/issues/2871>`_ : refactor: modernize-redundant-void-arg — drop (void) parameter lists
* `#2872 <https://github.com/CoolProp/CoolProp/issues/2872>`_ : refactor: modernize-make-shared — single-allocation shared_ptr construction
* `#2873 <https://github.com/CoolProp/CoolProp/issues/2873>`_ : refactor: bugprone-empty-catch — audit silent catch (...)  blocks
* `#2874 <https://github.com/CoolProp/CoolProp/issues/2874>`_ : refactor: cert-err58-cpp — static-init that can throw
* `#2885 <https://github.com/CoolProp/CoolProp/issues/2885>`_ : [ISSUE] Duplicated code at end of CoolPropLib.h in 7.2.0 shared_library SourceForge distribution
* `#2896 <https://github.com/CoolProp/CoolProp/issues/2896>`_ : R1132E should be R1132(E)
* `#2897 <https://github.com/CoolProp/CoolProp/issues/2897>`_ : Tabular methods should expose direct tabular evaluation
* `#2903 <https://github.com/CoolProp/CoolProp/issues/2903>`_ : CI contention on hashes.json
* `#2906 <https://github.com/CoolProp/CoolProp/issues/2906>`_ : HAPropsSI fails at T_wb = 0 °C across all RH (sub-freezing branch, #2690 follow-up)
* `#2913 <https://github.com/CoolProp/CoolProp/issues/2913>`_ : devdocs should only deploy on master
* `#2926 <https://github.com/CoolProp/CoolProp/issues/2926>`_ : clang-tidy sweep: actionable findings from files modified in the last month
* `#2973 <https://github.com/CoolProp/CoolProp/issues/2973>`_ : SVDSBTL: add DT-indexed surface to natively handle P(D, T) — supersedes BICUBIC inverter patching attempt (#1301)
* `#2988 <https://github.com/CoolProp/CoolProp/issues/2988>`_ : wheels not getting pushed to testpypi

Pull requests merged:

* `#2664 <https://github.com/CoolProp/CoolProp/pull/2664>`_ : Expose comprehensive AbstractState functionality and enums to WASM
* `#2669 <https://github.com/CoolProp/CoolProp/pull/2669>`_ : Document MEX function usage for MATLAB
* `#2677 <https://github.com/CoolProp/CoolProp/pull/2677>`_ : VB .Net README Update [skip ci]
* `#2679 <https://github.com/CoolProp/CoolProp/pull/2679>`_ : Fix HAPropsSI with T_dp+R or W+R inputs when T is unknown (issue #2670)
* `#2681 <https://github.com/CoolProp/CoolProp/pull/2681>`_ : Fix superancillary performance: store as shared_ptr instead of optional
* `#2682 <https://github.com/CoolProp/CoolProp/pull/2682>`_ : Some Refactoring of the Mathcad Wrapper for code succinctness and readability
* `#2683 <https://github.com/CoolProp/CoolProp/pull/2683>`_ : Update mixture tables for R1234ze(E) [skip ci]
* `#2687 <https://github.com/CoolProp/CoolProp/pull/2687>`_ : Install updated README.md to SourceForge with Mathcad binaries [skip ci]
* `#2688 <https://github.com/CoolProp/CoolProp/pull/2688>`_ : C++17 modernization: remove shims, stdlib idioms, attributes, filesystem
* `#2689 <https://github.com/CoolProp/CoolProp/pull/2689>`_ : Remove unused heavy headers to improve compile times
* `#2692 <https://github.com/CoolProp/CoolProp/pull/2692>`_ : make _binary arrays constexpr and use string_view instead of string
* `#2695 <https://github.com/CoolProp/CoolProp/pull/2695>`_ : Add PhaseSI and PropsSImulti to Mathcad Wrapper
* `#2697 <https://github.com/CoolProp/CoolProp/pull/2697>`_ : Fix HAPropsSI T_db from T_wb + low RelHum failing for narrow pressure band (issue #2690)
* `#2699 <https://github.com/CoolProp/CoolProp/pull/2699>`_ : fix(ci): skip REFPROP build for fork PRs where secrets are unavailable
* `#2700 <https://github.com/CoolProp/CoolProp/pull/2700>`_ : chore(deps): update miniz 3.0.2 → 3.1.1
* `#2702 <https://github.com/CoolProp/CoolProp/pull/2702>`_ : Add $$ to equations on Ideal Gas web page [skip ci]
* `#2704 <https://github.com/CoolProp/CoolProp/pull/2704>`_ : Revert "Fix HAPropsSI T_db from T_wb + low RelHum failing for narrow pressure band (issue #2690) (#2697)"
* `#2705 <https://github.com/CoolProp/CoolProp/pull/2705>`_ : test(humid_air): add comprehensive validation test suite
* `#2707 <https://github.com/CoolProp/CoolProp/pull/2707>`_ : feat(docs): add interactive 3D molecule viewers to fluid pages
* `#2708 <https://github.com/CoolProp/CoolProp/pull/2708>`_ : fix(docs): resolve all doxygen warnings
* `#2709 <https://github.com/CoolProp/CoolProp/pull/2709>`_ : Add Predefined Mixture Helper Functions to Mathcad Wrapper
* `#2710 <https://github.com/CoolProp/CoolProp/pull/2710>`_ : Skip REFPROP tests when not available during CI
* `#2713 <https://github.com/CoolProp/CoolProp/pull/2713>`_ : Add some more predefined mixtures from ASHRAE 34
* `#2715 <https://github.com/CoolProp/CoolProp/pull/2715>`_ : feat(GUI): add native Tauri/React desktop GUI for CoolProp
* `#2716 <https://github.com/CoolProp/CoolProp/pull/2716>`_ : perf(humid_air): cache virial/alpha0 coefficients to eliminate redundant EOS calls (26× speedup)
* `#2719 <https://github.com/CoolProp/CoolProp/pull/2719>`_ : docs: add predefined mixtures table to Mixtures page
* `#2722 <https://github.com/CoolProp/CoolProp/pull/2722>`_ : Mathcad wrapper compliance with CoolProp CLANG-format
* `#2723 <https://github.com/CoolProp/CoolProp/pull/2723>`_ : Mathcad wrapper additional error trapping for PropsSImulti
* `#2725 <https://github.com/CoolProp/CoolProp/pull/2725>`_ : Add refrigerant mixture models (Bell-JPCRD-2022/2023) and update R-1234yf EOS
* `#2728 <https://github.com/CoolProp/CoolProp/pull/2728>`_ : build: replace git submodules with CPM.cmake
* `#2729 <https://github.com/CoolProp/CoolProp/pull/2729>`_ : feat: add Chlorine (Cl2) fundamental equation of state with superancillaries
* `#2730 <https://github.com/CoolProp/CoolProp/pull/2730>`_ : Use requests with retry adapter for PubChem 3D/2D SDF downloads
* `#2731 <https://github.com/CoolProp/CoolProp/pull/2731>`_ : fix: repair clang-format CI for fork PRs and push events
* `#2732 <https://github.com/CoolProp/CoolProp/pull/2732>`_ : Add helper routine for HAPropsSI error handling in Mathcad wrapper
* `#2733 <https://github.com/CoolProp/CoolProp/pull/2733>`_ : fix(build): replace CPM.cmake dev snapshot with v0.40.7 release
* `#2734 <https://github.com/CoolProp/CoolProp/pull/2734>`_ : fix(cmake): change CACHE LIST to CACHE STRING for app sources and include dirs
* `#2735 <https://github.com/CoolProp/CoolProp/pull/2735>`_ : fix(build): default CPM_SOURCE_CACHE to enable incremental builds
* `#2741 <https://github.com/CoolProp/CoolProp/pull/2741>`_ : fix(cmake): Update FindMathematica dependency to latest version
* `#2744 <https://github.com/CoolProp/CoolProp/pull/2744>`_ : feat(cubics): add Chebyshev superancillaries for SRK and Peng-Robinson EOS
* `#2745 <https://github.com/CoolProp/CoolProp/pull/2745>`_ : fix(cubics): replace rounded SRK/PR Omega coefficients with exact values
* `#2747 <https://github.com/CoolProp/CoolProp/pull/2747>`_ : Mathematica wrapper web page .rst touchups [skip ci]
* `#2748 <https://github.com/CoolProp/CoolProp/pull/2748>`_ : Add mixture models from Bell-JCED-2025 and Bell-IJT-2020
* `#2749 <https://github.com/CoolProp/CoolProp/pull/2749>`_ : test(plot): regenerate isoline test data for exact SRK/PR coefficients
* `#2750 <https://github.com/CoolProp/CoolProp/pull/2750>`_ : fix(flash): return rhomolar_critical from PT_flash at the critical point (#2738)
* `#2752 <https://github.com/CoolProp/CoolProp/pull/2752>`_ : Ignore the .cache dir
* `#2758 <https://github.com/CoolProp/CoolProp/pull/2758>`_ : docs(index): auto-generate pure fluid count on landing page
* `#2759 <https://github.com/CoolProp/CoolProp/pull/2759>`_ : fix(docs): regenerate consistency plots when fluids are added
* `#2760 <https://github.com/CoolProp/CoolProp/pull/2760>`_ : test(plot): tolerate cross-platform libm jitter in isoline value checks
* `#2766 <https://github.com/CoolProp/CoolProp/pull/2766>`_ : fix(bib): correct DOIs for Bell-JPCRD-2022 and Bell-JPCRD-2023
* `#2768 <https://github.com/CoolProp/CoolProp/pull/2768>`_ : Batch update: 10 multiparameter EOS (R-1224yd(Z), R-1132(E), THF, PG, vinyl chloride, R-1123, n-C4/C5/C6 F, R-1233zd(E) intl std)
* `#2771 <https://github.com/CoolProp/CoolProp/pull/2771>`_ : fix(chlorine): align CL2 EOS with REFPROP reference implementation
* `#2775 <https://github.com/CoolProp/CoolProp/pull/2775>`_ : build: fetch boost subset via CPM from CoolProp/boost-headers
* `#2776 <https://github.com/CoolProp/CoolProp/pull/2776>`_ : test(superanc): extended-precision check points + EOS-freshness hash (#2772)
* `#2780 <https://github.com/CoolProp/CoolProp/pull/2780>`_ : Add R-1130(E); update R-1243zf to 3rd-gen EOS; harmonize source_eos_hash
* `#2781 <https://github.com/CoolProp/CoolProp/pull/2781>`_ : feat: NIST IR 8570 mixture pairs + R-1336mzz(Z) pure fluid
* `#2782 <https://github.com/CoolProp/CoolProp/pull/2782>`_ : fix(R1336mzzZ): inject SUPERANCILLARY block to fix sat-T-to-Tc test
* `#2784 <https://github.com/CoolProp/CoolProp/pull/2784>`_ : fix(R1336mzzZ): refresh hs_anchor h/s and SUPERANCILLARY after a1/a2 change
* `#2786 <https://github.com/CoolProp/CoolProp/pull/2786>`_ : ci(clang-format): drop push triggers, run PR-only (CoolProp-2uw.1)
* `#2788 <https://github.com/CoolProp/CoolProp/pull/2788>`_ : ci: add .pre-commit-config.yaml with clang-format hook (CoolProp-2uw.3)
* `#2789 <https://github.com/CoolProp/CoolProp/pull/2789>`_ : build(cmake): always export compile_commands.json (CoolProp-2uw.7)
* `#2790 <https://github.com/CoolProp/CoolProp/pull/2790>`_ : ci: skip CI on beads-only commits via shared commit-msg hook
* `#2791 <https://github.com/CoolProp/CoolProp/pull/2791>`_ : build(cmake): add format and format-check targets (CoolProp-2uw.2)
* `#2792 <https://github.com/CoolProp/CoolProp/pull/2792>`_ : ci: re-enable cppcheck as warning-only artifact (CoolProp-2uw.8)
* `#2793 <https://github.com/CoolProp/CoolProp/pull/2793>`_ : ci(codeql): run on PRs + master, bump to v3 actions (CoolProp-2uw.9)
* `#2794 <https://github.com/CoolProp/CoolProp/pull/2794>`_ : ci(coverity): schedule-only + dump JSON defects (CoolProp-2uw.10, closes #2717)
* `#2796 <https://github.com/CoolProp/CoolProp/pull/2796>`_ : docs: contributor code-quality workflow (CoolProp-2uw.4)
* `#2797 <https://github.com/CoolProp/CoolProp/pull/2797>`_ : ci: add IWYU report job as warning-only artifact (CoolProp-2uw.11)
* `#2798 <https://github.com/CoolProp/CoolProp/pull/2798>`_ : ci: add diff-only clang-tidy workflow on PRs (CoolProp-2uw.5)
* `#2799 <https://github.com/CoolProp/CoolProp/pull/2799>`_ : ci: add manual-stage clang-tidy pre-commit hook (CoolProp-2uw.6)
* `#2800 <https://github.com/CoolProp/CoolProp/pull/2800>`_ : fix: thread-safe lazy initialization of fluid libraries
* `#2801 <https://github.com/CoolProp/CoolProp/pull/2801>`_ : ci(clang-tidy): document informational-by-design
* `#2802 <https://github.com/CoolProp/CoolProp/pull/2802>`_ : ci(clang-tidy): switch to high-signal whitelist (CoolProp-2uw)
* `#2803 <https://github.com/CoolProp/CoolProp/pull/2803>`_ : style: whole-repo clang-format pass over src/ + include/ (CoolProp-2uw.12)
* `#2804 <https://github.com/CoolProp/CoolProp/pull/2804>`_ : chore: add .git-blame-ignore-revs for the whole-repo reformat (CoolProp-2uw.12)
* `#2805 <https://github.com/CoolProp/CoolProp/pull/2805>`_ : refactor: apply safe clang-tidy --fix passes (CoolProp-2uw.13, part 1)
* `#2806 <https://github.com/CoolProp/CoolProp/pull/2806>`_ : fix(cmake): IWYU wiring fires for both Main creation paths
* `#2807 <https://github.com/CoolProp/CoolProp/pull/2807>`_ : chore: blame-revs append + clang-tidy MathHeader (CoolProp-2uw.13b)
* `#2808 <https://github.com/CoolProp/CoolProp/pull/2808>`_ : refactor: cppcoreguidelines-init-variables --fix pass (CoolProp-2uw.13 part 2)
* `#2809 <https://github.com/CoolProp/CoolProp/pull/2809>`_ : chore: blame-ignore #2808 init-variables fix pass (CoolProp-2uw.13c)
* `#2810 <https://github.com/CoolProp/CoolProp/pull/2810>`_ : refactor: IWYU pass — drop 28 unused std-library includes (build speed)
* `#2814 <https://github.com/CoolProp/CoolProp/pull/2814>`_ : Expand Mathcad Wrapper Documentation [skip ci]
* `#2816 <https://github.com/CoolProp/CoolProp/pull/2816>`_ : build(deps): bump vite/vitest, override protocol-buffers-schema (combines #2811-2813)
* `#2817 <https://github.com/CoolProp/CoolProp/pull/2817>`_ : ci: enable Dependabot version updates for GitHub Actions
* `#2824 <https://github.com/CoolProp/CoolProp/pull/2824>`_ : ci(GUI): consolidate per-OS releases into one draft (idempotent release job)
* `#2826 <https://github.com/CoolProp/CoolProp/pull/2826>`_ : fix(GUI): unwrap react-plotly.js CJS default — fixes blank window on all platforms
* `#2827 <https://github.com/CoolProp/CoolProp/pull/2827>`_ : build(deps): bump Eigen 3.4.0 → 5.0.1
* `#2831 <https://github.com/CoolProp/CoolProp/pull/2831>`_ : fix(HumidAirProp): per-thread Water/Air backends via thread_local
* `#2832 <https://github.com/CoolProp/CoolProp/pull/2832>`_ : docs(pubchem): commit SDF cache + manifest-based freshness check
* `#2833 <https://github.com/CoolProp/CoolProp/pull/2833>`_ : fix(build): atomic dev/hashes.json write — fixes Windows installer race (#2828)
* `#2835 <https://github.com/CoolProp/CoolProp/pull/2835>`_ : feat(saturation): branch-selecting HQ/SQ/DQ flashes via guess.T (#2773)
* `#2837 <https://github.com/CoolProp/CoolProp/pull/2837>`_ : ci: deprecate Coverity Scan in favor of CodeQL + cppcheck
* `#2840 <https://github.com/CoolProp/CoolProp/pull/2840>`_ : Add Qmass (mass-basis vapor quality) support across CoolProp backends
* `#2841 <https://github.com/CoolProp/CoolProp/pull/2841>`_ : fix(cubic): add missing override on SRK/PR backend methods
* `#2842 <https://github.com/CoolProp/CoolProp/pull/2842>`_ : ci(clang-format): pin CI version via uvx to match pre-commit hook
* `#2843 <https://github.com/CoolProp/CoolProp/pull/2843>`_ : build(deps): cargo update — drop vulnerable rand 0.7.3/0.8.5 (Dependabot #8)
* `#2845 <https://github.com/CoolProp/CoolProp/pull/2845>`_ : ci(codeql): exclude vendored sources via paths-ignore config
* `#2846 <https://github.com/CoolProp/CoolProp/pull/2846>`_ : fix(pcsaft): widen num_sites multiply to size_t (CodeQL int-overflow)
* `#2847 <https://github.com/CoolProp/CoolProp/pull/2847>`_ : chore: retire generate_meta_info.py to .py.txt (orphan; closes 2 CodeQL jinja2 alerts)
* `#2848 <https://github.com/CoolProp/CoolProp/pull/2848>`_ : fix(plots): silence CodeQL uninitialized-locals in Plots/Common.py
* `#2849 <https://github.com/CoolProp/CoolProp/pull/2849>`_ : fix(cppcheck): four real bugs / suppressions across src/
* `#2850 <https://github.com/CoolProp/CoolProp/pull/2850>`_ : fix(HEOS): change_EOS throws on unknown EOS name (#1703)
* `#2851 <https://github.com/CoolProp/CoolProp/pull/2851>`_ : fix(HEOS): correct first_saturation_deriv for mixtures (#2091)
* `#2852 <https://github.com/CoolProp/CoolProp/pull/2852>`_ : fix(ancillaries): return NaN above the reducing T instead of UB / SIGFPE (#1611)
* `#2853 <https://github.com/CoolProp/CoolProp/pull/2853>`_ : fix(PropsSI): handle length-0 and length-1 array inputs (#2417)
* `#2854 <https://github.com/CoolProp/CoolProp/pull/2854>`_ : fix(INCOMP): throw a clean error for molar properties (#1908)
* `#2855 <https://github.com/CoolProp/CoolProp/pull/2855>`_ : fix(thread-safety): make deriv_counter atomic (#2844 race 1)
* `#2856 <https://github.com/CoolProp/CoolProp/pull/2856>`_ : test(water): regression tests for fixed flash bugs (#2079, #1730)
* `#2857 <https://github.com/CoolProp/CoolProp/pull/2857>`_ : fix(HEOS): re-evaluate alphar at converged T in HSU_D_flash (#1907)
* `#2858 <https://github.com/CoolProp/CoolProp/pull/2858>`_ : fix(PCSAFT): align gibbsmolar_residual with g_res = h_res - T*s_res (#1943)
* `#2859 <https://github.com/CoolProp/CoolProp/pull/2859>`_ : fix(HumidAir): Twb solver no longer brackets above Tsat(p) (#2255)
* `#2860 <https://github.com/CoolProp/CoolProp/pull/2860>`_ : fix(HEOS): validate Q before mutating state in update() (#2195)
* `#2861 <https://github.com/CoolProp/CoolProp/pull/2861>`_ : fix(INCOMP): psat throws below TminPsat instead of returning 0 silently (#2209)
* `#2862 <https://github.com/CoolProp/CoolProp/pull/2862>`_ : test(water): pin sub-zero water phase behaviour (#1098)
* `#2863 <https://github.com/CoolProp/CoolProp/pull/2863>`_ : build(MSVC): embed VERSIONINFO in CoolProp.dll (#2391)
* `#2864 <https://github.com/CoolProp/CoolProp/pull/2864>`_ : test: regression pins for fixed bugs (#2244, #2461)
* `#2865 <https://github.com/CoolProp/CoolProp/pull/2865>`_ : fix(REFPROP): DmolarSmolar honours imposed phase via DSFL1 (#2042)
* `#2866 <https://github.com/CoolProp/CoolProp/pull/2866>`_ : fix(REFPROP): wire DmassQ / DmolarQ inputs through DQFL2 (#1845)
* `#2875 <https://github.com/CoolProp/CoolProp/pull/2875>`_ : refactor(modernize-make-shared): single-allocation shared_ptr (#2872)
* `#2876 <https://github.com/CoolProp/CoolProp/pull/2876>`_ : refactor(vararg): replace direct printf/fprintf with iostream + format (#2870)
* `#2877 <https://github.com/CoolProp/CoolProp/pull/2877>`_ : refactor(modernize-avoid-c-arrays): convert local char arrays to std::array (#2869)
* `#2878 <https://github.com/CoolProp/CoolProp/pull/2878>`_ : fix(ci): use pull_request.head.sha for clang-format diff
* `#2886 <https://github.com/CoolProp/CoolProp/pull/2886>`_ : refactor(modernize-redundant-void-arg): drop (void) parameter lists (#2871)
* `#2887 <https://github.com/CoolProp/CoolProp/pull/2887>`_ : refactor(cert-err58-cpp): NOLINT static-init that may throw (#2874)
* `#2888 <https://github.com/CoolProp/CoolProp/pull/2888>`_ : fix(HEOS): mole_fractions_liquid/vapor reject single-phase states (#2308)
* `#2889 <https://github.com/CoolProp/CoolProp/pull/2889>`_ : feat(python): expose AbstractState::Qmass on Python interface
* `#2890 <https://github.com/CoolProp/CoolProp/pull/2890>`_ : ci(python): add PIP_RETRIES/PIP_TIMEOUT for cibw pip-install flakes
* `#2891 <https://github.com/CoolProp/CoolProp/pull/2891>`_ : fix(tabular): validate cell after sat-curve bump (BICUBIC PT segfault, #1950)
* `#2894 <https://github.com/CoolProp/CoolProp/pull/2894>`_ : feat(tabular): TABULAR_NX/NY config keys for grid resolution
* `#2895 <https://github.com/CoolProp/CoolProp/pull/2895>`_ : Fix equations and typos in MathcadWrappers documentation [skip ci]
* `#2904 <https://github.com/CoolProp/CoolProp/pull/2904>`_ : fix(fluids): rename R1132E to R1132(E); add R-1132E, R-1132(E) aliases (#2896)
* `#2905 <https://github.com/CoolProp/CoolProp/pull/2905>`_ : fix(build): retry hashes.json swap on Windows PermissionError (#2903)
* `#2907 <https://github.com/CoolProp/CoolProp/pull/2907>`_ : fix(HumidAir): T_db from (T_wb, RelHum, P) at low RH (#2690 modes B+C)
* `#2909 <https://github.com/CoolProp/CoolProp/pull/2909>`_ : fix(cmake): fetch Eigen via tarball to avoid Windows clone flakes
* `#2910 <https://github.com/CoolProp/CoolProp/pull/2910>`_ : ci(GUI): poll screenshot render check instead of fixed sleep
* `#2911 <https://github.com/CoolProp/CoolProp/pull/2911>`_ : ci: consolidate dev_* checks + gate wheel matrix on PRs
* `#2912 <https://github.com/CoolProp/CoolProp/pull/2912>`_ : fix(docs): pin 3Dmol.js download to jsDelivr CDN
* `#2914 <https://github.com/CoolProp/CoolProp/pull/2914>`_ : feat(svd,region): generic 2D rank-r SVD + curve-bounded region atlas
* `#2915 <https://github.com/CoolProp/CoolProp/pull/2915>`_ : feat(svd-2a): C++ multi-fluid SVD ρ(h,p) end-to-end validation harness
* `#2916 <https://github.com/CoolProp/CoolProp/pull/2916>`_ : fix(HEOS): RAII clear of imposed iphase_gas in HSU_D_flash::solver_resid
* `#2917 <https://github.com/CoolProp/CoolProp/pull/2917>`_ : feat(sbtl): SVDSurface adapter library + msgpack persistence (Phase 2b)
* `#2918 <https://github.com/CoolProp/CoolProp/pull/2918>`_ : ci(devdocs): drop pull_request trigger so devdocs deploys only on master
* `#2919 <https://github.com/CoolProp/CoolProp/pull/2919>`_ : feat(svdsbtl): SVDSBTL_BACKEND — full Phase 2c (backend + benchmark)
* `#2920 <https://github.com/CoolProp/CoolProp/pull/2920>`_ : feat(AbstractState): cache-bypassing batch fast_evaluate (Tabular + IF97)
* `#2921 <https://github.com/CoolProp/CoolProp/pull/2921>`_ : fix(region,clenshaw): drop redundant >= 1 + break in unsigned loop (CodeQL #3512)
* `#2922 <https://github.com/CoolProp/CoolProp/pull/2922>`_ : feat(svdsbtl): require explicit source-of-truth backend (HEOS/REFPROP/IF97)
* `#2923 <https://github.com/CoolProp/CoolProp/pull/2923>`_ : fix(if97): clear() must flush full cache; SVDSBTL viscosity + auto-conformance docs
* `#2924 <https://github.com/CoolProp/CoolProp/pull/2924>`_ : fix(cubic): pure-fluid DmolarT/DmassT via superancillary phase bracketing
* `#2925 <https://github.com/CoolProp/CoolProp/pull/2925>`_ : perf(HEOS): trim DHSU_T_flash imposed-phase overhead (#2718)
* `#2927 <https://github.com/CoolProp/CoolProp/pull/2927>`_ : Generic backend options via factory-string JSON (PR A — infrastructure)
* `#2928 <https://github.com/CoolProp/CoolProp/pull/2928>`_ : SVDSBTL opts in to factory-string options (PR B — schema + cache key)
* `#2929 <https://github.com/CoolProp/CoolProp/pull/2929>`_ : SVDSBTL critical-patch HEOS-fallback routing (PR C)
* `#2934 <https://github.com/CoolProp/CoolProp/pull/2934>`_ : feat(svdsbtl): Chebyshev η + Newton-refined IF97 + ε-band phase-pin
* `#2936 <https://github.com/CoolProp/CoolProp/pull/2936>`_ : test(svdsbtl): fix cache_present_for to match actual cache filenames
* `#2937 <https://github.com/CoolProp/CoolProp/pull/2937>`_ : ci(test_catch2): cache ~/.CoolProp/SVDTables across runs
* `#2938 <https://github.com/CoolProp/CoolProp/pull/2938>`_ : feat(svdsbtl): split SUPER at IF97 R2/R3 boundary + fix R3 sampling
* `#2940 <https://github.com/CoolProp/CoolProp/pull/2940>`_ : feat(svdsbtl): R1/R3 atlas split + TOMS748 h-inversion (combines reverted #2939 + #2940)
* `#2941 <https://github.com/CoolProp/CoolProp/pull/2941>`_ : docs(IF97): Python timing-profile script + figure for IF97 / HEOS / SVDSBTL
* `#2942 <https://github.com/CoolProp/CoolProp/pull/2942>`_ : feat(fast_evaluate): batch path for SVDSBTL + dome Q-blend for IF97/SVDSBTL
* `#2943 <https://github.com/CoolProp/CoolProp/pull/2943>`_ : fix(svdsbtl,docs): block SVDSBTL from PropsSI + speed up IF97Conformance.py ~3000x
* `#2944 <https://github.com/CoolProp/CoolProp/pull/2944>`_ : perf(svdsbtl): batched fast_evaluate — share locate+basis across surface outputs
* `#2945 <https://github.com/CoolProp/CoolProp/pull/2945>`_ : test(svdsbtl): close REFPROP-source coverage gaps (zx1)
* `#2946 <https://github.com/CoolProp/CoolProp/pull/2946>`_ : docs(svdsbtl,if97): new SVDSBTL page + restructure IF97 narrative (s60)
* `#2947 <https://github.com/CoolProp/CoolProp/pull/2947>`_ : feat(svdsbtl): auto-calibrate critical-patch bbox per fluid (5ni+dxd+229)
* `#2949 <https://github.com/CoolProp/CoolProp/pull/2949>`_ : fix(svdsbtl): skip patch polish for two-phase points
* `#2950 <https://github.com/CoolProp/CoolProp/pull/2950>`_ : test(svdsbtl): relax Q bit-exact compare to Approx for REFPROP fast_evaluate
* `#2951 <https://github.com/CoolProp/CoolProp/pull/2951>`_ : feat(svdsbtl): SuperAncillary-backed sat boundaries (CoolProp-8vg)
* `#2952 <https://github.com/CoolProp/CoolProp/pull/2952>`_ : feat(svdsbtl): fast eval surrogate for atlas curve_contains (~28% sub-critical speedup)
* `#2953 <https://github.com/CoolProp/CoolProp/pull/2953>`_ : feat(svdsbtl): parallel per-cell HEOS sampling during table build (CoolProp-43h)
* `#2954 <https://github.com/CoolProp/CoolProp/pull/2954>`_ : feat(ci): dev/ci/preflight.sh single-script pre-push gate (CoolProp-6r6)
* `#2955 <https://github.com/CoolProp/CoolProp/pull/2955>`_ : feat(ci): preflight clang-tidy works on macOS (per #2926)
* `#2956 <https://github.com/CoolProp/CoolProp/pull/2956>`_ : feat(ci): preflight gates on signal-filtered clang-tidy + code-reviewer banner (CoolProp-y51)
* `#2957 <https://github.com/CoolProp/CoolProp/pull/2957>`_ : feat(svdsbtl): extend SVDSBTL&IF97 coverage to IAPWS R5 (CoolProp-pd6)
* `#2958 <https://github.com/CoolProp/CoolProp/pull/2958>`_ : docs(if97-conformance): dense-knot R1/R3 isotherm + PDF emission + download link
* `#2959 <https://github.com/CoolProp/CoolProp/pull/2959>`_ : docs(svdsbtl): correct 19 doc bugs uncovered by audit (CoolProp-fhp)
* `#2960 <https://github.com/CoolProp/CoolProp/pull/2960>`_ : docs(if97,svdsbtl): correct 12 claims uncovered by IF97 audit
* `#2961 <https://github.com/CoolProp/CoolProp/pull/2961>`_ : feat(svdsbtl): implement ALTERNATIVE_SVDTABLES_DIRECTORY config key (CoolProp-fhp)
* `#2962 <https://github.com/CoolProp/CoolProp/pull/2962>`_ : fix(if97): correct n5 typo in R2/R3 boundary classifier (CoolProp-6zl)
* `#2963 <https://github.com/CoolProp/CoolProp/pull/2963>`_ : fix(svdsbtl): atomic write-temp + rename for cache files (CoolProp-4no.2)
* `#2964 <https://github.com/CoolProp/CoolProp/pull/2964>`_ : docs(svdsbtl): interactive SVDSBTL-vs-HEOS validation page (CoolProp-4no.4)
* `#2965 <https://github.com/CoolProp/CoolProp/pull/2965>`_ : feat(svdsbtl): SaturationSurrogate cubic-spline cache for REFPROP source (CoolProp-077)
* `#2966 <https://github.com/CoolProp/CoolProp/pull/2966>`_ : fix(svdsbtl): gate polish_patch_state_ to IF97 source only
* `#2967 <https://github.com/CoolProp/CoolProp/pull/2967>`_ : feat(svdsbtl): NC near-critical sub-regions with POWER(β=1/3) axis (CoolProp-4u9)
* `#2968 <https://github.com/CoolProp/CoolProp/pull/2968>`_ : feat(dev): SVDSBTL per-region sizing harness (CoolProp-6oe)
* `#2970 <https://github.com/CoolProp/CoolProp/pull/2970>`_ : fix(wasm): pin emsdk to CI sha + add mixture set_mole_fractions test
* `#2971 <https://github.com/CoolProp/CoolProp/pull/2971>`_ : feat(wasm): complete enum binding coverage + binding-surface gaps + tests
* `#2974 <https://github.com/CoolProp/CoolProp/pull/2974>`_ : fix(wasm): bypass register_optional<T> template to preserve _embind_register_optional import
* `#2976 <https://github.com/CoolProp/CoolProp/pull/2976>`_ : refactor(clang-tidy): modernize-use-override sweep (CoolProp-7ue)
* `#2977 <https://github.com/CoolProp/CoolProp/pull/2977>`_ : feat(wasm)!: native JS arrays via emscripten::val (replaces VectorDouble)
* `#2978 <https://github.com/CoolProp/CoolProp/pull/2978>`_ : refactor: clang-tidy bug-finder batch (CoolProp-m1q, #2926)
* `#2979 <https://github.com/CoolProp/CoolProp/pull/2979>`_ : fix(tests): per-PID tmpdir suffix on Catch2 fixtures (CoolProp-8ft, CoolProp-3t8)
* `#2980 <https://github.com/CoolProp/CoolProp/pull/2980>`_ : refactor: clang-tidy perf sweep + virtual-in-ctor + float-loops (CoolProp-kig, #2926)
* `#2981 <https://github.com/CoolProp/CoolProp/pull/2981>`_ : refactor(headers): break two header layering violations (CoolProp-6c1)
* `#2983 <https://github.com/CoolProp/CoolProp/pull/2983>`_ : refactor(clang-tidy): annotate documented silent catches in CoolPropLib batch APIs (#2873)
* `#2985 <https://github.com/CoolProp/CoolProp/pull/2985>`_ : refactor(clang-tidy): behaviour-preserving modernize sweep of CoolPropLib.cpp (CoolProp-jtl)
* `#2986 <https://github.com/CoolProp/CoolProp/pull/2986>`_ : fix(HumidAir): reject wet-bulb inputs in the water triple-point gap (#2906)
* `#2989 <https://github.com/CoolProp/CoolProp/pull/2989>`_ : ci: cancel superseded runs to cut CI backlog (concurrency groups)
* `#2990 <https://github.com/CoolProp/CoolProp/pull/2990>`_ : feat(REFPROP): saturated-state props + saturation/two-phase derivatives (#2013, supersedes #2016)
* `#2992 <https://github.com/CoolProp/CoolProp/pull/2992>`_ : refactor(svdsbtl): use IF97::Region23_p for B23 boundary, drop local copy
* `#2993 <https://github.com/CoolProp/CoolProp/pull/2993>`_ : fix(ci): publish unique TestPyPI wheel versions (#2988)
* `#2994 <https://github.com/CoolProp/CoolProp/pull/2994>`_ : docs: fix RST escaping and formatting issues in Web docs
* `#2995 <https://github.com/CoolProp/CoolProp/pull/2995>`_ : Superancillary happy-path for the D+{H,S,U} flashes (HEOS)
* `#2996 <https://github.com/CoolProp/CoolProp/pull/2996>`_ : ci: cancel a merged PR's in-flight CI runs
* `#2997 <https://github.com/CoolProp/CoolProp/pull/2997>`_ : style: clang-format whole codebase (pinned 18.1.8) [skip ci]
* `#2998 <https://github.com/CoolProp/CoolProp/pull/2998>`_ : ci(asan): drop detect_stack_use_after_return + add 90-min timeout (CoolProp-xuu)
* `#2999 <https://github.com/CoolProp/CoolProp/pull/2999>`_ : Superancillary happy-path for the H,S flash
* `#3000 <https://github.com/CoolProp/CoolProp/pull/3000>`_ : Correct D+{H,S,U} round-trip for sub-triple-point compressed liquid
* `#3001 <https://github.com/CoolProp/CoolProp/pull/3001>`_ : feat(plots/docs): consistency-plot failure reporting, per-panel timing, REFPROP support
* `#3003 <https://github.com/CoolProp/CoolProp/pull/3003>`_ : docs(virial): exact-virials/Axy derivations, plans & beads (#2991) [skip ci]
* `#3006 <https://github.com/CoolProp/CoolProp/pull/3006>`_ : fix(REFPROP): honor melting-line bound sentinels in calc_melting_line
* `#3009 <https://github.com/CoolProp/CoolProp/pull/3009>`_ : fix(HumidAir): drop a spurious 1e6 factor in the vbar_ws ice branch (#2657)
* `#3011 <https://github.com/CoolProp/CoolProp/pull/3011>`_ : fix(REFPROP): preserve the SATTP saturation pressure in the QT flash (#2671)
* `#3016 <https://github.com/CoolProp/CoolProp/pull/3016>`_ : fix: mask floating-point exceptions at the C-interface entry points
* `#3022 <https://github.com/CoolProp/CoolProp/pull/3022>`_ : fix: extend HEOS two-phase awareness to chemical_potential and fugacity_coefficient (follow-up to #2345)
* `#3026 <https://github.com/CoolProp/CoolProp/pull/3026>`_ : Multicomponent TP-flash: Michelsen stability analysis + phase-split for the Cubic and HEOS backends (changes mixture results / phase(); adds MIXTURE_STABILITY_ALGORITHM)
* `#3029 <https://github.com/CoolProp/CoolProp/pull/3029>`_ : fix: MPG2 incompressible viscosity was 10x too small (#1374)
* `#3031 <https://github.com/CoolProp/CoolProp/pull/3031>`_ : fix: HEOS melting-line PT guards honor DONT_CHECK_PROPERTY_LIMITS (#1936)
* `#3033 <https://github.com/CoolProp/CoolProp/pull/3033>`_ : Melting-line caloric seeding for the cold-compressed-liquid H,S flash
* `#3034 <https://github.com/CoolProp/CoolProp/pull/3034>`_ : perf: direct-EOS cache-bypass on warm probes in HSU_P_flash_singlephase_Brent
* `#3036 <https://github.com/CoolProp/CoolProp/pull/3036>`_ : Add the D2O (heavy water) melting-line ancillary (Herrig et al. 2018)
* `#3041 <https://github.com/CoolProp/CoolProp/pull/3041>`_ : fix: GERG-2008 / ExponentialDepartureFunction now calls phi.finish()
* `#3042 <https://github.com/CoolProp/CoolProp/pull/3042>`_ : feat(python): ship PEP 561 type stubs (CoolProp.pyi) for the Cython wrapper
* `#3081 <https://github.com/CoolProp/CoolProp/pull/3081>`_ : fix(REFPROP): correct the SATP phase flag / pressure / density in the PQ saturation fallback (GH #1502)
* `#3110 <https://github.com/CoolProp/CoolProp/pull/3110>`_ : feat(python): ship PEP 561 type stubs for the v8 wheel via nanobind stubgen
* `#3118 <https://github.com/CoolProp/CoolProp/pull/3118>`_ : feat(mixtures): add Tkaczuk et al. (2020) He-Ne, He-Ar, Ne-Ar cryogenic models (changes He-Ar results)
* `#3127 <https://github.com/CoolProp/CoolProp/pull/3127>`_ : ci: match the SONAME-versioned libCoolProp.so.<ver> in the symbol-leak gate
* `#3130 <https://github.com/CoolProp/CoolProp/pull/3130>`_ : fix(ci): name TestPyPI dev builds X.Y.Z.dev<ts> not .post<ts>
* `#3141 <https://github.com/CoolProp/CoolProp/pull/3141>`_ : perf(transport): cache the ECS reference fluid instead of rebuilding it per call
* `#3146 <https://github.com/CoolProp/CoolProp/pull/3146>`_ : fix(core): make debug_level + error/warning string globals thread-safe (per-thread error strings)

7.2.0
-----

Highlights:

* Added support for python 3.14
* Modernized the build system to use scikit-build-core, and supports uv now too 
* Added methods for ideal gas properties
* More tweaks to iterative routines

Issues closed:

* `#2589 <https://github.com/CoolProp/CoolProp/issues/2589>`_ : Properties of Air as an Ideal Gas
* `#2622 <https://github.com/CoolProp/CoolProp/issues/2622>`_ : R123 1phase PY issue

Pull requests merged:

* `#2579 <https://github.com/CoolProp/CoolProp/pull/2579>`_ : Improvements to mixture pair management
* `#2624 <https://github.com/CoolProp/CoolProp/pull/2624>`_ : Add Windows Instructions and example for Fluent Wrapper
* `#2625 <https://github.com/CoolProp/CoolProp/pull/2625>`_ : Faster alpha0
* `#2626 <https://github.com/CoolProp/CoolProp/pull/2626>`_ : Add Ideal gas methods
* `#2627 <https://github.com/CoolProp/CoolProp/pull/2627>`_ : Add fallback method for density solver for bad EOS
* `#2628 <https://github.com/CoolProp/CoolProp/pull/2628>`_ : CI add InnoSetup
* `#2630 <https://github.com/CoolProp/CoolProp/pull/2630>`_ : Python 3.14 support?
* `#2631 <https://github.com/CoolProp/CoolProp/pull/2631>`_ : RST->MD
* `#2632 <https://github.com/CoolProp/CoolProp/pull/2632>`_ : Modernize Python build system to use scikit-build-core
* `#2633 <https://github.com/CoolProp/CoolProp/pull/2633>`_ : Fix header gen on windows



7.1.0
-----

Highlights:

* Fixed performance regression for P+{H,S,U} and T+{H,S,U,D} caused by copies being made of the superancillary data structures
* Flash routine tuning to improve robustness in the critical region

Issues closed:

* `#2491 <https://github.com/CoolProp/CoolProp/issues/2491>`_ : Xenon properties near the critical point
* `#2582 <https://github.com/CoolProp/CoolProp/issues/2582>`_ : ExternalMedia CoolProp C02 : error in TestBasePropertiesTranscritical.mo (same problem with Helium)
* `#2587 <https://github.com/CoolProp/CoolProp/issues/2587>`_ : Handling of comma in 1,2dichloroethane in CP.get_aliases in Python wrapper
* `#2588 <https://github.com/CoolProp/CoolProp/issues/2588>`_ : gitrevision missing in docs
* `#2591 <https://github.com/CoolProp/CoolProp/issues/2591>`_ : Performance regression for P, {D,H,S,U} in two-phase
* `#2592 <https://github.com/CoolProp/CoolProp/issues/2592>`_ : Performance regression for T, {H,S,U,D} when two-phase
* `#2594 <https://github.com/CoolProp/CoolProp/issues/2594>`_ : Low-level interface CO2 transcritical phase problem
* `#2598 <https://github.com/CoolProp/CoolProp/issues/2598>`_ : PhaseSI broken in v7
* `#2606 <https://github.com/CoolProp/CoolProp/issues/2606>`_ : Cannot use 'import.meta' outside a module (at coolprop.js:1:474)
* `#2612 <https://github.com/CoolProp/CoolProp/issues/2612>`_ : Debugger crashing when breakpoints are set on Python 3.13

Pull requests merged:

* `#2601 <https://github.com/CoolProp/CoolProp/pull/2601>`_ : Always recalculate the critical point to be consistent with the use (…
* `#2604 <https://github.com/CoolProp/CoolProp/pull/2604>`_ : Fix/docs gitrevision
* `#2613 <https://github.com/CoolProp/CoolProp/pull/2613>`_ : Disable profiling by default with flag to enable
* `#2619 <https://github.com/CoolProp/CoolProp/pull/2619>`_ : Fix consistency plots


7.0.0
-----

Highlights:

* [BREAKING] Added superancillary functions for all pure fluids with multiparameter EOS; see docs: https://www.coolprop.org/coolprop/SuperAncillary.html . Huge speedup for some properties!
* [BREAKING] Changed EOS for methanol to match that in REFPROP 
* Finally(!) got the Catch2 tests to all pass
* Update docs to use the pydata sphinx theme 
* Improvements and changes to javascript wrapping
* Moved docs to be hosted on GitHub Pages (fixed https issues)
* Update rapidjson to final commit to ensure all libraries use the same symbols

Issues closed:

* `#2350 <https://github.com/CoolProp/CoolProp/issues/2350>`_ : coolprop.org doesn't support HTTPS
* `#2395 <https://github.com/CoolProp/CoolProp/issues/2395>`_ : Helium saturation routine bug
* `#2409 <https://github.com/CoolProp/CoolProp/issues/2409>`_ : Can't open `https://coolprop.org` redirected to `https://.sourceforge.net/` end in "Site not found"
* `#2451 <https://github.com/CoolProp/CoolProp/issues/2451>`_ : Glitches with cyclopentane in CoolProp 6.6.0
* `#2484 <https://github.com/CoolProp/CoolProp/issues/2484>`_ : Catch2 workflow Fails on PRs and Pushes
* `#2512 <https://github.com/CoolProp/CoolProp/issues/2512>`_ : CoolProp fails to import when cadquery is imported first
* `#2535 <https://github.com/CoolProp/CoolProp/issues/2535>`_ : PC-SAFT test failures for electrolytes
* `#2536 <https://github.com/CoolProp/CoolProp/issues/2536>`_ : Script to generate test points for plotting
* `#2538 <https://github.com/CoolProp/CoolProp/issues/2538>`_ : Switch Methanol to the same EOS as in REFPROP
* `#2540 <https://github.com/CoolProp/CoolProp/issues/2540>`_ : use std::array for cache in Helmholtz terms
* `#2541 <https://github.com/CoolProp/CoolProp/issues/2541>`_ : A way to get required buffer size for `get_fluid_param_string`
* `#2546 <https://github.com/CoolProp/CoolProp/issues/2546>`_ : Restore asan tests action
* `#2558 <https://github.com/CoolProp/CoolProp/issues/2558>`_ : querying `Tau` and `Delta` from multicomponent models fails
* `#2559 <https://github.com/CoolProp/CoolProp/issues/2559>`_ : Logo missing for docs
* `#2563 <https://github.com/CoolProp/CoolProp/issues/2563>`_ : Switch universal gas constant to the final number
* `#2568 <https://github.com/CoolProp/CoolProp/issues/2568>`_ : Performance Regression in C# Bindings: HAPropsSI and Refrigerant PropsSI in 6.8.0/6.8.1dev
* `#2570 <https://github.com/CoolProp/CoolProp/issues/2570>`_ : Fire devdocs on every commit to main/master
* `#2571 <https://github.com/CoolProp/CoolProp/issues/2571>`_ : Profile AbstractState construction and add ENVVAR to disable superancillary entirely

Pull requests merged:

* `#2511 <https://github.com/CoolProp/CoolProp/pull/2511>`_ : Superancillaries for pure fluids
* `#2544 <https://github.com/CoolProp/CoolProp/pull/2544>`_ : Implement the array-based caching in AbstractState.
* `#2548 <https://github.com/CoolProp/CoolProp/pull/2548>`_ : Add an action to run asan
* `#2551 <https://github.com/CoolProp/CoolProp/pull/2551>`_ : Modernize (and test!) embind wrapping for Javascript
* `#2562 <https://github.com/CoolProp/CoolProp/pull/2562>`_ : Switch to pydata theme
* `#2580 <https://github.com/CoolProp/CoolProp/pull/2580>`_ : Export additional JavaScript functions
* `#2585 <https://github.com/CoolProp/CoolProp/pull/2585>`_ : Update rapidjson to latest trunk to align with python-rapidjson

6.8.0
-----

Highlights:

* Implement the `COOLPROP_REFPROP_ROOT` environment variable to make REFPROP integration easier
* Allow the configuration variables to be set via environment variables

Issues closed:

* `#1829 <https://github.com/CoolProp/CoolProp/issues/1829>`_ : CoolProp can't load REFPROP v10 mixtures on Debian 9
* `#1870 <https://github.com/CoolProp/CoolProp/issues/1870>`_ : n-pentane inconsistency properties
* `#2400 <https://github.com/CoolProp/CoolProp/issues/2400>`_ : CoolProp not able to compute subcooled liquid MM  properties
* `#2445 <https://github.com/CoolProp/CoolProp/issues/2445>`_ : Incorrect Enthalpy Value Using CoolProp-REFPROP Wrapper
* `#2447 <https://github.com/CoolProp/CoolProp/issues/2447>`_ : Incorrect Vapor Pressure Output for INCOMP::S800 (Syltherm 800)
* `#2472 <https://github.com/CoolProp/CoolProp/issues/2472>`_ : Workflow Python cibuildwheel failing on upload to TestPyPi
* `#2477 <https://github.com/CoolProp/CoolProp/issues/2477>`_ : How do I programme an isentropic compression to 140°C condensation temperature with Coolprop? (got a ValueError)
* `#2485 <https://github.com/CoolProp/CoolProp/issues/2485>`_ : 11 fluids in Brines and Solutions have exclusive range check on fraction_max where an inclusive check would be an improvement.
* `#2488 <https://github.com/CoolProp/CoolProp/issues/2488>`_ : Incompressible fitting code is no longer working
* `#2497 <https://github.com/CoolProp/CoolProp/issues/2497>`_ : Unicode error in State class with coolprop 6.7.0
* `#2501 <https://github.com/CoolProp/CoolProp/issues/2501>`_ : Add options for REFPROP PATH to be set by environment variables
* `#2517 <https://github.com/CoolProp/CoolProp/issues/2517>`_ : Document COOLPROP_REFPROP_ROOT
* `#2518 <https://github.com/CoolProp/CoolProp/issues/2518>`_ : Fix gitrevision when using pipx run build
* `#2519 <https://github.com/CoolProp/CoolProp/issues/2519>`_ : Allow all config variables to be overwritten by environment variables at first load
* `#2520 <https://github.com/CoolProp/CoolProp/issues/2520>`_ : Fix Octave building
* `#2531 <https://github.com/CoolProp/CoolProp/issues/2531>`_ : Test that __gitrevision__ is being properly populated

Pull requests merged:

* `#2478 <https://github.com/CoolProp/CoolProp/pull/2478>`_ : Bump fmtlib to 11.1.3
* `#2482 <https://github.com/CoolProp/CoolProp/pull/2482>`_ : Fix a bunch of bugs with D,U inputs in REFPROP backend
* `#2496 <https://github.com/CoolProp/CoolProp/pull/2496>`_ : Add set_cubic_alpha_C function to python wrapper
* `#2532 <https://github.com/CoolProp/CoolProp/pull/2532>`_ : Overwrite config with env


6.7.0
-----

Highlights:

* Added wheels for Python 3.13
* Added EOS for R-1336mzz(E)

Issues closed:

* `#2307 <https://github.com/CoolProp/CoolProp/issues/2307>`_ : Feature Request: Implementation of the novel Refrigerant R1336mzz(e) into HEOS Backend
* `#2325 <https://github.com/CoolProp/CoolProp/issues/2325>`_ : HEOS::R1336MZZE does not work properly
* `#2430 <https://github.com/CoolProp/CoolProp/issues/2430>`_ : Support for python 3.13
* `#2462 <https://github.com/CoolProp/CoolProp/issues/2462>`_ : Python3.13 getting same error as #1876 Python 3.8 : Error in import #1876

Pull requests merged:

* `#2309 <https://github.com/CoolProp/CoolProp/pull/2309>`_ : Addition of EoS JSON of R1336mzz(E) from Akasaka-IJT-2023
* `#2329 <https://github.com/CoolProp/CoolProp/pull/2329>`_ : Fix the triple points densities for R-1336mzz(E)
* `#2341 <https://github.com/CoolProp/CoolProp/pull/2341>`_ : Update build process for LibreOffice wrapper
* `#2342 <https://github.com/CoolProp/CoolProp/pull/2342>`_ : Include stdbool for C interface
* `#2347 <https://github.com/CoolProp/CoolProp/pull/2347>`_ : Update CMakeLists.txt
* `#2351 <https://github.com/CoolProp/CoolProp/pull/2351>`_ : Fix bug introduced by 443a2fd
* `#2362 <https://github.com/CoolProp/CoolProp/pull/2362>`_ : update Excel-wrapper for Mac docs
* `#2396 <https://github.com/CoolProp/CoolProp/pull/2396>`_ : Add link to `coolprop-mat` repo
* `#2402 <https://github.com/CoolProp/CoolProp/pull/2402>`_ : Use temperature dependent hard sphere diameter for ion term in ePC-SAFT
* `#2404 <https://github.com/CoolProp/CoolProp/pull/2404>`_ : Improve Robustness of IF97 Reverse (P,H) and (P,S) Evaluations Along Saturation Curve
* `#2415 <https://github.com/CoolProp/CoolProp/pull/2415>`_ : Update index.rst
* `#2416 <https://github.com/CoolProp/CoolProp/pull/2416>`_ : Bypass Mathcad builder on pushes and PRs
* `#2418 <https://github.com/CoolProp/CoolProp/pull/2418>`_ : Update the Release Workflow for upload/download-artifact@v4 [skip ci]
* `#2419 <https://github.com/CoolProp/CoolProp/pull/2419>`_ : Fix libreoffice_builder
* `#2436 <https://github.com/CoolProp/CoolProp/pull/2436>`_ : Python 3.13 and replace distutils with setuptools PEP 632 – Deprecate distutils
* `#2439 <https://github.com/CoolProp/CoolProp/pull/2439>`_ : Fix np.NaN for numpy >=2
* `#2446 <https://github.com/CoolProp/CoolProp/pull/2446>`_ : Fix Plots of log-p-h diagrams
* `#2450 <https://github.com/CoolProp/CoolProp/pull/2450>`_ : Update msgpack-c and selectively add boost
* `#2463 <https://github.com/CoolProp/CoolProp/pull/2463>`_ : Support plotting in C++
* `#2471 <https://github.com/CoolProp/CoolProp/pull/2471>`_ : Get Mathcad Workflow Running Again

6.6.0
-----

Highlights:

* Added wheels for Python 3.12
* Added new functions to the library interface
* Include new binaries in the release workflow (Mathcad, Javascript)
* Fixed the base temperature bug in the compressible backend

Issues closed:

* `#1944 <https://github.com/CoolProp/CoolProp/issues/1944>`_ : coolprop in MathCAD - current version problem
* `#2198 <https://github.com/CoolProp/CoolProp/issues/2198>`_ : Binary folders have a "v" in their name
* `#2260 <https://github.com/CoolProp/CoolProp/issues/2260>`_ : manual install files not showing up on soundforge
* `#2278 <https://github.com/CoolProp/CoolProp/issues/2278>`_ : Javascript Wrapper unable to compiled for version coolprop 6.5.0
* `#2284 <https://github.com/CoolProp/CoolProp/issues/2284>`_ : Getting fugacity from the shared library
* `#2295 <https://github.com/CoolProp/CoolProp/issues/2295>`_ : Calculating enthalpy and entropy at exactly the middle value between min and max temperature does not work
* `#2310 <https://github.com/CoolProp/CoolProp/issues/2310>`_ : Wheels for Python 3.12

Pull requests merged:

* `#2213 <https://github.com/CoolProp/CoolProp/pull/2213>`_ : Use lazy initialization and avoid static objects
* `#2275 <https://github.com/CoolProp/CoolProp/pull/2275>`_ : CoolProp::apply_simple_mixing_rule missing from SWIG Wrapper
* `#2279 <https://github.com/CoolProp/CoolProp/pull/2279>`_ : Add javascript to the release script
* `#2286 <https://github.com/CoolProp/CoolProp/pull/2286>`_ : Add fugacity functions needed for compatability with `CoolProp.jl` pkg
* `#2291 <https://github.com/CoolProp/CoolProp/pull/2291>`_ : Added second_partial_deriv and  first_two_phase_deriv
* `#2294 <https://github.com/CoolProp/CoolProp/pull/2294>`_ : Update Mathcad docs for pre-compiled and discontinuation of Legacy Mathcad [skip-ci]
* `#2296 <https://github.com/CoolProp/CoolProp/pull/2296>`_ : Mathcad Wrapper README Formatting
* `#2317 <https://github.com/CoolProp/CoolProp/pull/2317>`_ : New C Interface Functions
* `#2320 <https://github.com/CoolProp/CoolProp/pull/2320>`_ : Food properties as incompressible liquids + ice
* `#2322 <https://github.com/CoolProp/CoolProp/pull/2322>`_ : Actions for Mathcad
* `#2323 <https://github.com/CoolProp/CoolProp/pull/2323>`_ : Incompresible versions of CoolProp fluids
* `#2324 <https://github.com/CoolProp/CoolProp/pull/2324>`_ : Fix base temperature and composition problems for incompressible fluids

6.5.0
-----

Highlights:

* Mostly small bugfixes and dependency updates
* Added ability to add predefined mixtures at runtime
* Updated transport models for CO2
* Fixed a bug in the hexane models

Issues closed:

* `#2051 <https://github.com/CoolProp/CoolProp/issues/2051>`_ : Cyclopentane EOS needs to be updated
* `#2142 <https://github.com/CoolProp/CoolProp/issues/2142>`_ : R1233zd does not work in the saturation region close to the bubble point
* `#2200 <https://github.com/CoolProp/CoolProp/issues/2200>`_ : CoolProp pure Hexane bug
* `#2201 <https://github.com/CoolProp/CoolProp/issues/2201>`_ : N-heptane has repeated IdealGasHelmholtzCP0AlyLee
* `#2205 <https://github.com/CoolProp/CoolProp/issues/2205>`_ : Python silently crashes when calling trivial_keyed_output on binary mixtures without specified mole fractions
* `#2251 <https://github.com/CoolProp/CoolProp/issues/2251>`_ : Unable to compile with fmt 10.0.0
* `#2265 <https://github.com/CoolProp/CoolProp/issues/2265>`_ : Sharp non-differentiable changes in thermal conductivity of CO2 and other gases
* `#2277 <https://github.com/CoolProp/CoolProp/issues/2277>`_ : Update State class

Pull requests merged:

* `#2203 <https://github.com/CoolProp/CoolProp/pull/2203>`_ : Provide better feedback for bad DQ inputs
* `#2207 <https://github.com/CoolProp/CoolProp/pull/2207>`_ : Verify that mole fractions are set before using them
* `#2214 <https://github.com/CoolProp/CoolProp/pull/2214>`_ : Change links from Google group to GitHub discussions
* `#2223 <https://github.com/CoolProp/CoolProp/pull/2223>`_ : Topic 2142
* `#2225 <https://github.com/CoolProp/CoolProp/pull/2225>`_ : update cyclopentane.json
* `#2230 <https://github.com/CoolProp/CoolProp/pull/2230>`_ : Topic-2200: Correct typo in n-Hexane rhoV auxilliary
* `#2238 <https://github.com/CoolProp/CoolProp/pull/2238>`_ : Incomp liqna
* `#2241 <https://github.com/CoolProp/CoolProp/pull/2241>`_ : Update index.rst
* `#2252 <https://github.com/CoolProp/CoolProp/pull/2252>`_ : Update fmt submodule to 10.0.0
* `#2261 <https://github.com/CoolProp/CoolProp/pull/2261>`_ : Create CITATION.bib
* `#2267 <https://github.com/CoolProp/CoolProp/pull/2267>`_ : implemented TCX Huber-JPCRD-2016 for CO2
* `#2268 <https://github.com/CoolProp/CoolProp/pull/2268>`_ : implemented VISC LAESECKE-JPCRD-2017-CO2


6.4.3
-----

Highlights:

* The first automated release that updates the homepage and all binaries

Issues closed:

* `#2196 <https://github.com/CoolProp/CoolProp/issues/2196>`_ : Automatically publish release binaries
* `#2197 <https://github.com/CoolProp/CoolProp/issues/2197>`_ : Add sdist for Python


6.4.2
-----

Highlights:

* The first release after 2 years
* Fixed the values in the vicinity of the critical point of ammonia
* Added Python wheels for Python 3.6 through 3.11 on many different architectures
* Added a reverse T(p,h) function to IF97
* Exposed more functions in the CoolPropLib interface
* Fixed a faulty density calculation for ice
* Added PC-SAFT as indepedent backend

Deprecated:

* Dropped support for Python 2.x

Issues Closed:

* `#1867 <https://github.com/CoolProp/CoolProp/issues/1867>`_ : TypeError after importing CoolProp / pip installation on Raspberry Pi
* `#1884 <https://github.com/CoolProp/CoolProp/issues/1884>`_ : Typo in enthalpy's unit of measure
* `#1962 <https://github.com/CoolProp/CoolProp/issues/1962>`_ : Ammonia (and maybe other?) calculations fail at the critical point
* `#1963 <https://github.com/CoolProp/CoolProp/issues/1963>`_ : Some examples don't work in docs
* `#1974 <https://github.com/CoolProp/CoolProp/issues/1974>`_ : Fix reducing density for Nitrogen
* `#1980 <https://github.com/CoolProp/CoolProp/issues/1980>`_ : Wrong alias in "R1243zf.json"
* `#1981 <https://github.com/CoolProp/CoolProp/issues/1981>`_ : Python CoolProp package doesn't work on Python 3.9.0 (32 bit and 64 bit)
* `#1992 <https://github.com/CoolProp/CoolProp/issues/1992>`_ : Installation errors with Python 3.9
* `#1999 <https://github.com/CoolProp/CoolProp/issues/1999>`_ :  PropsSI failed ungracefully with Water::IF97
* `#2003 <https://github.com/CoolProp/CoolProp/issues/2003>`_ : build error on MacOS 11.2 Big Sur
* `#2010 <https://github.com/CoolProp/CoolProp/issues/2010>`_ : cannot build the object library (COOLPROP_OBJECT_LIBRARY)
* `#2017 <https://github.com/CoolProp/CoolProp/issues/2017>`_ : I'm not able to install the coolprop with pip in python ...
* `#2020 <https://github.com/CoolProp/CoolProp/issues/2020>`_ : PC-SAFT integration
* `#2025 <https://github.com/CoolProp/CoolProp/issues/2025>`_ : Error in HAPropsSI when using enthalpy as an input (Excel VBA)
* `#2033 <https://github.com/CoolProp/CoolProp/issues/2033>`_ : Compatibility with Silicon chip in MacOS Big Sur 11.5.1
* `#2043 <https://github.com/CoolProp/CoolProp/issues/2043>`_ : Cannot create propertyplot for ammonia
* `#2049 <https://github.com/CoolProp/CoolProp/issues/2049>`_ : PropsSI("PHASE") calculate with ammonia, get error "options.p is not valid in saturation_T_pure_1D_P"
* `#2052 <https://github.com/CoolProp/CoolProp/issues/2052>`_ : How to install Coolprop in MacOS which has M1 chip?
* `#2053 <https://github.com/CoolProp/CoolProp/issues/2053>`_ : Small rounding issues for water
* `#2054 <https://github.com/CoolProp/CoolProp/issues/2054>`_ : Rounding for reducing density for R236ea
* `#2055 <https://github.com/CoolProp/CoolProp/issues/2055>`_ : Rounding for reducing density for nitrogen
* `#2067 <https://github.com/CoolProp/CoolProp/issues/2067>`_ : Adding a new fluid and compiled it. Not working when function is used.
* `#2073 <https://github.com/CoolProp/CoolProp/issues/2073>`_ : PHI0 density derivatives with REFPROP backend are wrong
* `#2078 <https://github.com/CoolProp/CoolProp/issues/2078>`_ : Python 3.8: Error in import
* `#2081 <https://github.com/CoolProp/CoolProp/issues/2081>`_ : Add support to release linux aarch64 wheels
* `#2095 <https://github.com/CoolProp/CoolProp/issues/2095>`_ : Issue when compiling shared library in docker on M1 - unrecognized command-line option ‘-m64’
* `#2100 <https://github.com/CoolProp/CoolProp/issues/2100>`_ : Cubic Mixtures: ideal gas contribution doesn't work properly (Rcomponent is wrong))
* `#2113 <https://github.com/CoolProp/CoolProp/issues/2113>`_ : Installation failed when using command: pip install coolprop
* `#2114 <https://github.com/CoolProp/CoolProp/issues/2114>`_ : Trouble installing MATLAB wrapper via Python
* `#2119 <https://github.com/CoolProp/CoolProp/issues/2119>`_ : Python bindings: Call for help from the community
* `#2126 <https://github.com/CoolProp/CoolProp/issues/2126>`_ : CoolProp 6.4.2dev0, MATLAB wrapper with Python 3.9
* `#2149 <https://github.com/CoolProp/CoolProp/issues/2149>`_ : Bug in the departure function parameters for GeneralizedHFC in CoolProp
* `#2178 <https://github.com/CoolProp/CoolProp/issues/2178>`_ : Please update github release
* `#2184 <https://github.com/CoolProp/CoolProp/issues/2184>`_ : CoolProp Online throwing internal error
* `#2186 <https://github.com/CoolProp/CoolProp/issues/2186>`_ : Ammonia critical point issue behaviour
* `#2187 <https://github.com/CoolProp/CoolProp/issues/2187>`_ : The online version of CoolProp cannot work
* `#2190 <https://github.com/CoolProp/CoolProp/issues/2190>`_ : Humid air property function HAPropsSI is not reversible
* `#2192 <https://github.com/CoolProp/CoolProp/issues/2192>`_ : Update the changelog for v6.4.2

Pull requests merged:

* `#1977 <https://github.com/CoolProp/CoolProp/pull/1977>`_ : Add Rust Wrapper
* `#1990 <https://github.com/CoolProp/CoolProp/pull/1990>`_ : Fix cxx17
* `#1993 <https://github.com/CoolProp/CoolProp/pull/1993>`_ : LibreOffice: Use pip for installing CoolProp python package
* `#2005 <https://github.com/CoolProp/CoolProp/pull/2005>`_ : Fix cxx17
* `#2008 <https://github.com/CoolProp/CoolProp/pull/2008>`_ : Fix build on macOS
* `#2011 <https://github.com/CoolProp/CoolProp/pull/2011>`_ : A minor correction in case of COOLPROP_OBJECT_LIBRARY=ON
* `#2050 <https://github.com/CoolProp/CoolProp/pull/2050>`_ : Update index.rst for the C# Wrapper
* `#2056 <https://github.com/CoolProp/CoolProp/pull/2056>`_ : Fix typo in iQ description
* `#2058 <https://github.com/CoolProp/CoolProp/pull/2058>`_ : IF97 Backend Q and Phase Patch
* `#2062 <https://github.com/CoolProp/CoolProp/pull/2062>`_ : Updated info for the C# Wrapper
* `#2076 <https://github.com/CoolProp/CoolProp/pull/2076>`_ : Included CoolPropJavascriptDemo
* `#2084 <https://github.com/CoolProp/CoolProp/pull/2084>`_ : Add functions to CoolPropLib
* `#2097 <https://github.com/CoolProp/CoolProp/pull/2097>`_ : Add github action to build python wheels (including python 3.9 and 3.10)
* `#2098 <https://github.com/CoolProp/CoolProp/pull/2098>`_ : Github Actions: add shared library and doxygen workflows.
* `#2101 <https://github.com/CoolProp/CoolProp/pull/2101>`_ : Fix Rcomponent in calc_alpha0_deriv_nocache
* `#2103 <https://github.com/CoolProp/CoolProp/pull/2103>`_ : Lint: use automated tooling to reformat C++ and CMakeLists files
* `#2105 <https://github.com/CoolProp/CoolProp/pull/2105>`_ : Bump Catch  1 to Catch v3.0.0-preview4
* `#2106 <https://github.com/CoolProp/CoolProp/pull/2106>`_ : Cppcheck workflow
* `#2107 <https://github.com/CoolProp/CoolProp/pull/2107>`_ : Add bound-check to setter and getter functions
* `#2108 <https://github.com/CoolProp/CoolProp/pull/2108>`_ : Format macros + strip trailing whitespaces
* `#2109 <https://github.com/CoolProp/CoolProp/pull/2109>`_ : Configure upload to pypi/testpypi
* `#2110 <https://github.com/CoolProp/CoolProp/pull/2110>`_ : Fix mac cibuildwheel
* `#2116 <https://github.com/CoolProp/CoolProp/pull/2116>`_ : Fix mac sed
* `#2118 <https://github.com/CoolProp/CoolProp/pull/2118>`_ : Python bindings upload to (test)pypi fixes
* `#2120 <https://github.com/CoolProp/CoolProp/pull/2120>`_ : Missing a py37 build for Windows x64 + fix py38 win32 and py39 win32
* `#2122 <https://github.com/CoolProp/CoolProp/pull/2122>`_ : Simplify CoolProp python bindings cibuildwheel
* `#2132 <https://github.com/CoolProp/CoolProp/pull/2132>`_ : Bump IF97 to included reverse T(P,H) patch [skip ci]
* `#2133 <https://github.com/CoolProp/CoolProp/pull/2133>`_ : New functions for CoolPropLib
* `#2134 <https://github.com/CoolProp/CoolProp/pull/2134>`_ : Add fluid_param_string and get_JSONstring to cubic backend
* `#2135 <https://github.com/CoolProp/CoolProp/pull/2135>`_ : AbstractState functions for CoolPropLib
* `#2143 <https://github.com/CoolProp/CoolProp/pull/2143>`_ : Corrected rho_ice route by replacing g_ice with dg_dp_Ice in Ice.cpp
* `#2146 <https://github.com/CoolProp/CoolProp/pull/2146>`_ : Bump FindMathematica to most recent version
* `#2161 <https://github.com/CoolProp/CoolProp/pull/2161>`_ : improve PC-SAFT flash
* `#2164 <https://github.com/CoolProp/CoolProp/pull/2164>`_ : Updated info about SharpProp (3-party wrapper for C#)
* `#2165 <https://github.com/CoolProp/CoolProp/pull/2165>`_ : Added info about PyFluids (3-party wrapper for Python)
* `#2173 <https://github.com/CoolProp/CoolProp/pull/2173>`_ : Prevent crashes near critical density due to saturation calc
* `#2176 <https://github.com/CoolProp/CoolProp/pull/2176>`_ : add PCSAFT page in docs
* `#2191 <https://github.com/CoolProp/CoolProp/pull/2191>`_ : Build the docs for v6.4.2


6.4.1
-----

Highlights:

* Fixed the nightly builds
* Added documentation for the NuGet packages by MadsKirkFoged

Deprecated:

* Removed the n-dimensional input for Python due to too many bugs

Issues Closed:

* `#1960 <https://github.com/CoolProp/CoolProp/issues/1960>`_ : Docs do not build
* `#1952 <https://github.com/CoolProp/CoolProp/issues/1952>`_ : CoolProp Module not Found
* `#1942 <https://github.com/CoolProp/CoolProp/issues/1942>`_ : Help with DLL hell
* `#1940 <https://github.com/CoolProp/CoolProp/issues/1940>`_ : CoolProp import doesn't work in v6.4.0

Pull Requests merged:

* `#1957 <https://github.com/CoolProp/CoolProp/pull/1957>`_ : Update BaseObjects.py
* `#1949 <https://github.com/CoolProp/CoolProp/pull/1949>`_ : fixed typo with units for sigma
* `#1947 <https://github.com/CoolProp/CoolProp/pull/1947>`_ : Update index.rst
* `#1932 <https://github.com/CoolProp/CoolProp/pull/1932>`_ : Allow nD-array input to PropsSI


6.4.0
-----

Highlights:

* Added a working version of PC-SAFT EOS (big thanks to Zach Baird)
* Updated EOS for neon, helium
* Added Python 3.8 wheels

Deprecated:

* Python 2.7 interface. This will be the last release supporting Python 2.7
* 32-bit builds on OSX

Issues Closed:

* `#1922 <https://github.com/CoolProp/CoolProp/issues/1922>`_ : Bug: The density of air is off by a magnitude of 10 given the P and T inputs specified below
* `#1915 <https://github.com/CoolProp/CoolProp/issues/1915>`_ : Error in Low Level Interface example code on coolprop.org
* `#1857 <https://github.com/CoolProp/CoolProp/issues/1857>`_ : calc_name is not implemented for REFPROP backend
* `#1856 <https://github.com/CoolProp/CoolProp/issues/1856>`_ : List of incompressible fluids/mixtures is missing
* `#1855 <https://github.com/CoolProp/CoolProp/issues/1855>`_ : Python AbstractState has no attribute 'compressibility_factor'
* `#1848 <https://github.com/CoolProp/CoolProp/issues/1848>`_ : PR density calc fail

Pull Requests merged:

* `#1921 <https://github.com/CoolProp/CoolProp/pull/1921>`_ : fix typos in pcsaft json
* `#1916 <https://github.com/CoolProp/CoolProp/pull/1916>`_ : Fix second_partial_deriv example
* `#1913 <https://github.com/CoolProp/CoolProp/pull/1913>`_ : Fix pcsaft flash
* `#1906 <https://github.com/CoolProp/CoolProp/pull/1906>`_ : Update index.rst
* `#1903 <https://github.com/CoolProp/CoolProp/pull/1903>`_ : Update index.rst adding PropiedadesDeFluidos Tool
* `#1901 <https://github.com/CoolProp/CoolProp/pull/1901>`_ : for CDNJS auto-update request
* `#1896 <https://github.com/CoolProp/CoolProp/pull/1896>`_ : Update init.py byte string in split for Py3.8
* `#1891 <https://github.com/CoolProp/CoolProp/pull/1891>`_ : Fix uppercase-only fluid naming
* `#1885 <https://github.com/CoolProp/CoolProp/pull/1885>`_ : fixed typo in getos() and else statement
* `#1881 <https://github.com/CoolProp/CoolProp/pull/1881>`_ : Fixed crash in PTflash_twophase::build_arrays
* `#1878 <https://github.com/CoolProp/CoolProp/pull/1878>`_ : Py38
* `#1877 <https://github.com/CoolProp/CoolProp/pull/1877>`_ : Adding PC-SAFT EOS
* `#1875 <https://github.com/CoolProp/CoolProp/pull/1875>`_ : Document apply_simple_mixing_rule initialisation behaviour
* `#1842 <https://github.com/CoolProp/CoolProp/pull/1842>`_ : fixes in doxygen docu for fundamental derivative
* `#1838 <https://github.com/CoolProp/CoolProp/pull/1838>`_ : Allow list delimiter to be changed
* `#1837 <https://github.com/CoolProp/CoolProp/pull/1837>`_ : Never allow two-phase SOS. 


6.3.0
-----

Highlights:

* The molar mass gets now reset properly - affected mixture calculations with changing concentrations.
* The humid air calculations check the inputs and outputs according to the limits from the publication.
* The isentropic expansion coefficient can now be accessed directly.
* ... and a lot of little bugfixes (see issues)

Issues Closed:

* `#1820 <https://github.com/CoolProp/CoolProp/issues/1820>`_ : Humid air example fails due to new limits
* `#1815 <https://github.com/CoolProp/CoolProp/issues/1815>`_ : molar_mass not getting cleared
* `#1811 <https://github.com/CoolProp/CoolProp/issues/1811>`_ : Humid air properties above 388 K
* `#1786 <https://github.com/CoolProp/CoolProp/issues/1786>`_ : Incompressible docs not building properly
* `#1784 <https://github.com/CoolProp/CoolProp/issues/1784>`_ : Sphinx builder still broken
* `#1782 <https://github.com/CoolProp/CoolProp/issues/1782>`_ : OSX 10.14 builds
* `#1778 <https://github.com/CoolProp/CoolProp/issues/1778>`_ : There are no nightly builds after 2018/11/04
* `#1777 <https://github.com/CoolProp/CoolProp/issues/1777>`_ : Building from the PyPI sdist on Python 3.7 results in compilation errors on macOS
* `#1775 <https://github.com/CoolProp/CoolProp/issues/1775>`_ : Tmin function cannot be evaluated at 0.0 concentration for incomp fluids
* `#1763 <https://github.com/CoolProp/CoolProp/issues/1763>`_ : Mathcad 15 binary builds dropped as of version 6.2
* `#1762 <https://github.com/CoolProp/CoolProp/issues/1762>`_ : IF97 Documentation Page Error
* `#1760 <https://github.com/CoolProp/CoolProp/issues/1760>`_ : Android Wrapper error at 6.2.0 and 6.2.2dev
* `#1759 <https://github.com/CoolProp/CoolProp/issues/1759>`_ : Memory leak in Mathematica interface
* `#1758 <https://github.com/CoolProp/CoolProp/issues/1758>`_ : Build AbstractState object from stored tabular data
* `#1756 <https://github.com/CoolProp/CoolProp/issues/1756>`_ : Issue with incompressible fluid in v6.2.1
* `#1753 <https://github.com/CoolProp/CoolProp/issues/1753>`_ : numpy.core.multiarray failed to import
* `#1752 <https://github.com/CoolProp/CoolProp/issues/1752>`_ : Add fluids to CoolProp if you are using matlab
* `#1748 <https://github.com/CoolProp/CoolProp/issues/1748>`_ : Apostrophe should be escaped in '...' strings or be used in "..." string
* `#1745 <https://github.com/CoolProp/CoolProp/issues/1745>`_ : Surface Tension calculation failing ungracefully
* `#1744 <https://github.com/CoolProp/CoolProp/issues/1744>`_ : Error value Excel on Mac
* `#1742 <https://github.com/CoolProp/CoolProp/issues/1742>`_ : r404a JT valve error
* `#1741 <https://github.com/CoolProp/CoolProp/issues/1741>`_ : Wrapper for OCTAVE-4.x.x 32 bit for Windows
* `#1734 <https://github.com/CoolProp/CoolProp/issues/1734>`_ : Can I access the older 'Props' functions in CoolProp 6
* `#1732 <https://github.com/CoolProp/CoolProp/issues/1732>`_ : Error 53 Excel Wrapper MacBook Issue
* `#1731 <https://github.com/CoolProp/CoolProp/issues/1731>`_ : Will CoolProp contain R513a refrigerant properties in the near future??
* `#1724 <https://github.com/CoolProp/CoolProp/issues/1724>`_ : REFPROP v10.0 enthalpy/pressure look-up bug
* `#1717 <https://github.com/CoolProp/CoolProp/issues/1717>`_ : Coolprop cannot work on local JS
* `#1708 <https://github.com/CoolProp/CoolProp/issues/1708>`_ : calc_phase missing for REFPROP backend.
* `#1707 <https://github.com/CoolProp/CoolProp/issues/1707>`_ : Inconsistency in incompressible backend when getting fluid name()
* `#1700 <https://github.com/CoolProp/CoolProp/issues/1700>`_ : HAPropsSi broken in 6.1.1dev for R(p,T,h)
* `#1618 <https://github.com/CoolProp/CoolProp/issues/1618>`_ : Excel and Python wrapper, HAPropsSI problem with P_w input
* `#1601 <https://github.com/CoolProp/CoolProp/issues/1601>`_ : ValueError: HAProps failed ungracefully with inputs: "W","T",2.7852603934336025e+02,"P",1.0132500000000000e+05,"H",1.4114309647910737e+04
* `#1465 <https://github.com/CoolProp/CoolProp/issues/1465>`_ : Humid air calculations when dry bulb is above freezing but wet bulb is below freezing

Pull Requests merged:

* `#1823 <https://github.com/CoolProp/CoolProp/pull/1823>`_ : Robustify humid air limit checks
* `#1814 <https://github.com/CoolProp/CoolProp/pull/1814>`_ : Feature/docs py3
* `#1813 <https://github.com/CoolProp/CoolProp/pull/1813>`_ : Added limits to the humid air properties, closes #1811
* `#1810 <https://github.com/CoolProp/CoolProp/pull/1810>`_ : Use std::shared_ptr if available in VS
* `#1797 <https://github.com/CoolProp/CoolProp/pull/1797>`_ : Set _phase in REFPROP backend and return it in calc_phase()
* `#1791 <https://github.com/CoolProp/CoolProp/pull/1791>`_ : Add isentropic expansion coefficient
* `#1787 <https://github.com/CoolProp/CoolProp/pull/1787>`_ : Add conda install instructions from conda-forge
* `#1783 <https://github.com/CoolProp/CoolProp/pull/1783>`_ : Feature/recent compilers
* `#1773 <https://github.com/CoolProp/CoolProp/pull/1773>`_ : Make wrapper for LibreOffice an extension
* `#1770 <https://github.com/CoolProp/CoolProp/pull/1770>`_ : Disambiguate IF97 Fluid Names - Water only
* `#1767 <https://github.com/CoolProp/CoolProp/pull/1767>`_ : Add IF97 links to CoolProp homepage and backends page
* `#1765 <https://github.com/CoolProp/CoolProp/pull/1765>`_ : Patch PropsSI for INCOMP fluid with zero mass fraction
* `#1761 <https://github.com/CoolProp/CoolProp/pull/1761>`_ : Free input strings in Mathematica interface
* `#1755 <https://github.com/CoolProp/CoolProp/pull/1755>`_ : Throw sensible error message if single-phase surface tension requested
* `#1751 <https://github.com/CoolProp/CoolProp/pull/1751>`_ : update syntax to julia 1.0

6.2.1
-----

Only a minor fix to the Javascript wrapper


6.2.0
-----

New Features:

* Added a new EoS for heavy water
* Added pre-built wheels for Python 3.6 and Python 3.7
* The MATLAB wrappers have been abandoned in favour of Python-based calls
* Add phase specification to high-level interface
* LabVIEW VIs can now call PropsSI and Props1SI
* Added a wrapper for Android
* ... and a lot of little bugfixes (see issues)

Issues Closed:

* `#1699 <https://github.com/CoolProp/CoolProp/issues/1699>`_ : Pip Install problem with Ubuntu 18.04 and Python 3.6.3
* `#1682 <https://github.com/CoolProp/CoolProp/issues/1682>`_ : Coolprop.Coolprop module does not exist
* `#1672 <https://github.com/CoolProp/CoolProp/issues/1672>`_ : In ODEintegrator, limits are wrong for backwards
* `#1662 <https://github.com/CoolProp/CoolProp/issues/1662>`_ : Wrong results when using INCOMP (incompressible) fluids using high-level interface after importing matplotlib.pyplot
* `#1661 <https://github.com/CoolProp/CoolProp/issues/1661>`_ : install fail with python 3.6 in anaconda for win10
* `#1659 <https://github.com/CoolProp/CoolProp/issues/1659>`_ : More reducing state resetting needed when fractions set
* `#1654 <https://github.com/CoolProp/CoolProp/issues/1654>`_ : Version 6.1 with Python 3.6.1 Anaconda (64-bit)
* `#1652 <https://github.com/CoolProp/CoolProp/issues/1652>`_ : Problem with saturated vapor internal energy calculations with quality/density inputs
* `#1649 <https://github.com/CoolProp/CoolProp/issues/1649>`_ : Cannot cimport CoolProp into cython in python 3.6
* `#1647 <https://github.com/CoolProp/CoolProp/issues/1647>`_ : Parsing of Mixtures depends on LOCALE
* `#1630 <https://github.com/CoolProp/CoolProp/issues/1630>`_ : Predefined mixture cannot have uppercase .MIX
* `#1629 <https://github.com/CoolProp/CoolProp/issues/1629>`_ : Phase envelopes fail for predefined mixtures with REFPROP backend
* `#1607 <https://github.com/CoolProp/CoolProp/issues/1607>`_ : Tabular Backend Fails with HmassP_INPUTS when iphase_twophase Imposed
* `#1604 <https://github.com/CoolProp/CoolProp/issues/1604>`_ : v6.2?
* `#1603 <https://github.com/CoolProp/CoolProp/issues/1603>`_ : Parse out Zero Mass Fraction Components in High-Level Interface
* `#1602 <https://github.com/CoolProp/CoolProp/issues/1602>`_ : hmass() gives strange result after calling update() with PQ_INPUTS in specific case
* `#1582 <https://github.com/CoolProp/CoolProp/issues/1582>`_ : Buildbot update
* `#1563 <https://github.com/CoolProp/CoolProp/issues/1563>`_ : Unify AbstractState's behavior when using HEOS or the tabular interpolations schemes
* `#1551 <https://github.com/CoolProp/CoolProp/issues/1551>`_ : Import of matplotlib.pyplot results in error for mixtures
* `#1530 <https://github.com/CoolProp/CoolProp/issues/1530>`_ : Catch tests failing
* `#1455 <https://github.com/CoolProp/CoolProp/issues/1455>`_ : apply_simple_mixing_rule broken
* `#1439 <https://github.com/CoolProp/CoolProp/issues/1439>`_ : Wrong dam_dtau for Twu
* `#1426 <https://github.com/CoolProp/CoolProp/issues/1426>`_ : UNIFAQ compile errors
* `#1422 <https://github.com/CoolProp/CoolProp/issues/1422>`_ : Ttriple wrong for REFPROP for water
* `#1406 <https://github.com/CoolProp/CoolProp/issues/1406>`_ : StateContainer print fails
* `#1396 <https://github.com/CoolProp/CoolProp/issues/1396>`_ : Formulas are wrong for dichloroethane and ethylene oxide
* `#1393 <https://github.com/CoolProp/CoolProp/issues/1393>`_ : Crash when set_mole_fractions() not called
* `#1381 <https://github.com/CoolProp/CoolProp/issues/1381>`_ : Calling acentric factor with cubic equation
* `#1372 <https://github.com/CoolProp/CoolProp/issues/1372>`_ : inconsistent result with mixture of "Ethylbenzene[0.5]&P-XYLENE[0.5]"<>"ethylbenzene[0.5]&P-XYLENE[0.5]"
* `#1371 <https://github.com/CoolProp/CoolProp/issues/1371>`_ : Get JSON string for fluid at runtime
* `#1369 <https://github.com/CoolProp/CoolProp/issues/1369>`_ : Return 'n/a'  REFPROP version when not loaded or supported.
* `#1368 <https://github.com/CoolProp/CoolProp/issues/1368>`_ : segmentation fault when calling get_global_param_string("REFPROP_version") from Python
* `#1366 <https://github.com/CoolProp/CoolProp/issues/1366>`_ : Allow fluids to be provided (and overwritten) at runtime; closes #1345
* `#1365 <https://github.com/CoolProp/CoolProp/issues/1365>`_ : SMath Wrapper refactoring
* `#1362 <https://github.com/CoolProp/CoolProp/issues/1362>`_ : LabVIEW VIs to call PropsSI and Props1SI
* `#1361 <https://github.com/CoolProp/CoolProp/issues/1361>`_ : Re-enable alpha0 mixture derivative tests for cubics
* `#1359 <https://github.com/CoolProp/CoolProp/issues/1359>`_ : Allow for cubic transformations in HEOS multi-fluid model
* `#1355 <https://github.com/CoolProp/CoolProp/issues/1355>`_ : SMath Wrapper refactoring
* `#1354 <https://github.com/CoolProp/CoolProp/issues/1354>`_ : splined properties, _fluid_type and _phase not cleared in AbstractState.h
* `#1352 <https://github.com/CoolProp/CoolProp/issues/1352>`_ : Faulty state update with DmassT_inputs in HEOS backend with specified phase
* `#1350 <https://github.com/CoolProp/CoolProp/issues/1350>`_ : Simulation error when using ExternalMedia in Dymola
* `#1348 <https://github.com/CoolProp/CoolProp/issues/1348>`_ : Allow alpha0 to be provided for cubic EOS
* `#1347 <https://github.com/CoolProp/CoolProp/issues/1347>`_ : Add ability to ignore setup errors for REFPROP mixtures
* `#1343 <https://github.com/CoolProp/CoolProp/issues/1343>`_ : Call git in the dev folder
* `#1339 <https://github.com/CoolProp/CoolProp/issues/1339>`_ : Set a standard departure function through the AbstractState
* `#1333 <https://github.com/CoolProp/CoolProp/issues/1333>`_ : Make it possible to use x[i]=0 in some alpha0 derivatives
* `#1329 <https://github.com/CoolProp/CoolProp/issues/1329>`_ : DO NOT allow for over-writing of departure functions when loading defaults
* `#1328 <https://github.com/CoolProp/CoolProp/issues/1328>`_ : Dmass wrong for saturated states for REFPROP
* `#1325 <https://github.com/CoolProp/CoolProp/issues/1325>`_ : Also export HAProps_Aux to pybind11 interface
* `#1324 <https://github.com/CoolProp/CoolProp/issues/1324>`_ : Figure out problem with linux wheels
* `#1323 <https://github.com/CoolProp/CoolProp/issues/1323>`_ : Added PQ and QT Input Pairs to provide Saturation Values to IF97 Backend
* `#1322 <https://github.com/CoolProp/CoolProp/issues/1322>`_ : Bigger buffer size for Julia wrapper
* `#1321 <https://github.com/CoolProp/CoolProp/issues/1321>`_ : Finally fix phase envelopes again
* `#1320 <https://github.com/CoolProp/CoolProp/issues/1320>`_ : Figure out why catch runs take forever on "*"nix
* `#1319 <https://github.com/CoolProp/CoolProp/issues/1319>`_ : Fix python windows builds
* `#1318 <https://github.com/CoolProp/CoolProp/issues/1318>`_ : Move entire Emscripten interface into its own file that is included separately
* `#1317 <https://github.com/CoolProp/CoolProp/issues/1317>`_ : Loading HMX.BNC through the DLL yields weird behavior
* `#1316 <https://github.com/CoolProp/CoolProp/issues/1316>`_ : Added configuration options for MSVCRT linking, changed the output di…
* `#1314 <https://github.com/CoolProp/CoolProp/issues/1314>`_ : Android Wrapper
* `#1312 <https://github.com/CoolProp/CoolProp/issues/1312>`_ : First step toward #1310
* `#1309 <https://github.com/CoolProp/CoolProp/issues/1309>`_ : version 6.1.0 not available from pypi
* `#1308 <https://github.com/CoolProp/CoolProp/issues/1308>`_ : Add Trivial Parameter calls to IF97 Backend
* `#1307 <https://github.com/CoolProp/CoolProp/issues/1307>`_ : get_config_string returns nothing in python
* `#1306 <https://github.com/CoolProp/CoolProp/issues/1306>`_ : Typo in CO2+Argon coefficients
* `#1305 <https://github.com/CoolProp/CoolProp/issues/1305>`_ : Fix some warnings in MSVC 2015
* `#1304 <https://github.com/CoolProp/CoolProp/issues/1304>`_ : Parse refprop HMX.BNC file and load coefficients
* `#1303 <https://github.com/CoolProp/CoolProp/issues/1303>`_ : call refprop from coolprop in scilab on linux
* `#1302 <https://github.com/CoolProp/CoolProp/issues/1302>`_ : Export cubic's alpha functions
* `#1300 <https://github.com/CoolProp/CoolProp/issues/1300>`_ : Add criticality_contour_values to pybind11 interface
* `#1299 <https://github.com/CoolProp/CoolProp/issues/1299>`_ : Add config keys to pybind11 interface
* `#1298 <https://github.com/CoolProp/CoolProp/issues/1298>`_ : HAPropsSI H, p, W lookups not working past 5.0.0
* `#1296 <https://github.com/CoolProp/CoolProp/issues/1296>`_ : Add phases enum to pybind11 interface
* `#1295 <https://github.com/CoolProp/CoolProp/issues/1295>`_ : Specify the minimum delta for spinodal tracer as config variable
* `#1294 <https://github.com/CoolProp/CoolProp/issues/1294>`_ : Add parser for HMX.BNC from REFPROP
* `#1292 <https://github.com/CoolProp/CoolProp/issues/1292>`_ : Source zip file on SourceForge is not correct again...
* `#1289 <https://github.com/CoolProp/CoolProp/issues/1289>`_ : Make triple point accessible in HEOS::get_fluid_constant
* `#1285 <https://github.com/CoolProp/CoolProp/issues/1285>`_ : Allow fluids to be overwritten
* `#1279 <https://github.com/CoolProp/CoolProp/issues/1279>`_ : Add phase specification to high-level interface
* `#1253 <https://github.com/CoolProp/CoolProp/issues/1253>`_ : Implement derivatives of alpha0 w.r.t. tau, delta
* `#1249 <https://github.com/CoolProp/CoolProp/issues/1249>`_ : IF97 Error in CoolProp Wrapping for SMath
* `#969 <https://github.com/CoolProp/CoolProp/issues/969>`_ : Support mixtures with component mole fractions of zero

Pull Requests merged:

* `#1675 <https://github.com/CoolProp/CoolProp/pull/1675>`_ : Let DARWIN build with libc++
* `#1666 <https://github.com/CoolProp/CoolProp/pull/1666>`_ : Make string->float conversion aware of the locale
* `#1665 <https://github.com/CoolProp/CoolProp/pull/1665>`_ : Patches PropsSI imposed phase for backends other than HEOS
* `#1660 <https://github.com/CoolProp/CoolProp/pull/1660>`_ : Update PropsSI() to Parse Imposed Phase Strings on Input Keys
* `#1656 <https://github.com/CoolProp/CoolProp/pull/1656>`_ : Mistake in function 'inline_label' in CoolProp/Plots/Common.py
* `#1645 <https://github.com/CoolProp/CoolProp/pull/1645>`_ : Provide return string from PhaseSI() if phase can't be determined.
* `#1609 <https://github.com/CoolProp/CoolProp/pull/1609>`_ : editorconfig
* `#1606 <https://github.com/CoolProp/CoolProp/pull/1606>`_ : Patch PT_flash() to update _phase with imposed phase, in case it changed
* `#1464 <https://github.com/CoolProp/CoolProp/pull/1464>`_ : Fix a few REFPROP functions; closes #1422
* `#1460 <https://github.com/CoolProp/CoolProp/pull/1460>`_ : Greatly improve the stability of REFPROP mixture calls at saturation …
* `#1459 <https://github.com/CoolProp/CoolProp/pull/1459>`_ : Call SATTP properly when QT inputs are given for REFPROP
* `#1458 <https://github.com/CoolProp/CoolProp/pull/1458>`_ : Actually set the Twu parameters if provided
* `#1457 <https://github.com/CoolProp/CoolProp/pull/1457>`_ : Add ierr checks to calls to SETKTV
* `#1450 <https://github.com/CoolProp/CoolProp/pull/1450>`_ : Fix typo in CoolPropLib.h
* `#1449 <https://github.com/CoolProp/CoolProp/pull/1449>`_ : Move F2K into emscripten_interface.cxx
* `#1448 <https://github.com/CoolProp/CoolProp/pull/1448>`_ : Update the ODE integrator to allow it to integrate backwards
* `#1376 <https://github.com/CoolProp/CoolProp/pull/1376>`_ : Update HumidAirProp.cpp

6.1.0
-----

New features:

* Windows installer for Microsoft Excel
* Added VTPR backend
* Twu and Mathias-Copeman attractive parameters can be set for PR and SRK
* Major improvements to Excel wrapper
* Added EOS for MDM of M. Thol
* Implemented first version of PT flash calculations for two-phase states
* Implemented PT flash for mixtures (not finished)
* Added a pybind11 module for CoolProp
* ... and a lot of little bugfixes (see issues)

Contributors to this release:
ibell, JonWel, jowr, babaksamareh, mikekaganski

* `#1290 <https://github.com/CoolProp/CoolProp/issues/1290>`_ : Catch runs should be Release builds
* `#1288 <https://github.com/CoolProp/CoolProp/issues/1288>`_ : Actually check if T < Tmelt for p > pmin
* `#1287 <https://github.com/CoolProp/CoolProp/issues/1287>`_ : Actually commit new pybind11 submodule
* `#1286 <https://github.com/CoolProp/CoolProp/issues/1286>`_ : in phase envelope construction, potential crash
* `#1284 <https://github.com/CoolProp/CoolProp/issues/1284>`_ : Make low-level interface accessible through high-level interface in FORTRAN
* `#1283 <https://github.com/CoolProp/CoolProp/issues/1283>`_ : Add pure fluid check to VTPR
* `#1282 <https://github.com/CoolProp/CoolProp/issues/1282>`_ : Correct typo, see #1270
* `#1281 <https://github.com/CoolProp/CoolProp/issues/1281>`_ : Add ability to add HEOS fluids as JSON at runtime
* `#1272 <https://github.com/CoolProp/CoolProp/issues/1272>`_ : Solves a bug in VTPR
* `#1271 <https://github.com/CoolProp/CoolProp/issues/1271>`_ : Remove possible division by 0, closes #1270
* `#1269 <https://github.com/CoolProp/CoolProp/issues/1269>`_ : SatL and SatV of type VTPR too
* `#1268 <https://github.com/CoolProp/CoolProp/issues/1268>`_ : Implement fluid_names for cubic backend
* `#1267 <https://github.com/CoolProp/CoolProp/issues/1267>`_ : PengRobinson doesn't pass alpha to SatL and SatV
* `#1266 <https://github.com/CoolProp/CoolProp/issues/1266>`_ : Small fixes for VTPR
* `#1264 <https://github.com/CoolProp/CoolProp/issues/1264>`_ : Update initialization for VTPR
* `#1262 <https://github.com/CoolProp/CoolProp/issues/1262>`_ : Set alpha function in JSON
* `#1261 <https://github.com/CoolProp/CoolProp/issues/1261>`_ : Update CMakeLists.txt
* `#1259 <https://github.com/CoolProp/CoolProp/issues/1259>`_ : Methanol-water mixture: strange results
* `#1258 <https://github.com/CoolProp/CoolProp/issues/1258>`_ : Solves a bug with cubic and mixtures
* `#1257 <https://github.com/CoolProp/CoolProp/issues/1257>`_ : Update iPhone compilation docs
* `#1255 <https://github.com/CoolProp/CoolProp/issues/1255>`_ : Allow ability to set Twu parameters for cubic EOS (from JSON)
* `#1252 <https://github.com/CoolProp/CoolProp/issues/1252>`_ : Implement set_double_array2D
* `#1250 <https://github.com/CoolProp/CoolProp/issues/1250>`_ : Implement coefficient derivatives of dYr_dxi in reducing function
* `#1248 <https://github.com/CoolProp/CoolProp/issues/1248>`_ : Problem with OSX compilation
* `#1240 <https://github.com/CoolProp/CoolProp/issues/1240>`_ : Make psi_plus public
* `#1239 <https://github.com/CoolProp/CoolProp/issues/1239>`_ : Shortcut VTPR when pure fluids, solves #1232
* `#1237 <https://github.com/CoolProp/CoolProp/issues/1237>`_ : Create an installer for selected Windows wrappers
* `#1235 <https://github.com/CoolProp/CoolProp/issues/1235>`_ : Excel 2016 Add-In Updates
* `#1234 <https://github.com/CoolProp/CoolProp/issues/1234>`_ : Add the ability to set limits in Consistency plots
* `#1232 <https://github.com/CoolProp/CoolProp/issues/1232>`_ : VTPR components with one group
* `#1230 <https://github.com/CoolProp/CoolProp/issues/1230>`_ : Allow ability to call REFPROP on OSX
* `#1229 <https://github.com/CoolProp/CoolProp/issues/1229>`_ : ConsistencyPlots updates
* `#1227 <https://github.com/CoolProp/CoolProp/issues/1227>`_ : Make all functions in DepartureFunction overrridable
* `#1226 <https://github.com/CoolProp/CoolProp/issues/1226>`_ : More critical point questions
* `#1222 <https://github.com/CoolProp/CoolProp/issues/1222>`_ : Critical point calc failure
* `#1221 <https://github.com/CoolProp/CoolProp/issues/1221>`_ : Take more steps in stability evaluator (at least 100)
* `#1220 <https://github.com/CoolProp/CoolProp/issues/1220>`_ : Add adaptive integrator code
* `#1219 <https://github.com/CoolProp/CoolProp/issues/1219>`_ : Double post_update in update_TP_guessrho
* `#1217 <https://github.com/CoolProp/CoolProp/issues/1217>`_ : Peng-Robinson issue with Hydrogen
* `#1215 <https://github.com/CoolProp/CoolProp/issues/1215>`_ : Vapour QT_INPUT with VTPR
* `#1214 <https://github.com/CoolProp/CoolProp/issues/1214>`_ : Refactor exceptions in CoolPropLib.cpp close #1200
* `#1213 <https://github.com/CoolProp/CoolProp/issues/1213>`_ : Add tests for Poling example with UNIFAC code
* `#1212 <https://github.com/CoolProp/CoolProp/issues/1212>`_ : Add derivatives of a*rho with respect to tau,delta,x
* `#1211 <https://github.com/CoolProp/CoolProp/issues/1211>`_ : Use aii_term and b0_ii from cubic
* `#1209 <https://github.com/CoolProp/CoolProp/issues/1209>`_ : Correct tau derivatives in VTPR
* `#1208 <https://github.com/CoolProp/CoolProp/issues/1208>`_ : Correct derivatives of am and test for VTPR
* `#1206 <https://github.com/CoolProp/CoolProp/issues/1206>`_ : Segmentation fault when calling get_mass_fractions() with SRK and PR
* `#1204 <https://github.com/CoolProp/CoolProp/issues/1204>`_ : Make all functions in reducing function const
* `#1203 <https://github.com/CoolProp/CoolProp/issues/1203>`_ : Allow VTPR to pass only names by setting default R_u value
* `#1202 <https://github.com/CoolProp/CoolProp/issues/1202>`_ : Better error message when UNIFAC component cannot be found
* `#1201 <https://github.com/CoolProp/CoolProp/issues/1201>`_ : Update MixtureDerivatives.cpp
* `#1199 <https://github.com/CoolProp/CoolProp/issues/1199>`_ : dalpha0_dxi is wrong
* `#1198 <https://github.com/CoolProp/CoolProp/issues/1198>`_ : Cubic CP
* `#1197 <https://github.com/CoolProp/CoolProp/issues/1197>`_ : Cubic QT_INPUTS
* `#1196 <https://github.com/CoolProp/CoolProp/issues/1196>`_ : Update CoolPropLib.def
* `#1195 <https://github.com/CoolProp/CoolProp/issues/1195>`_ : Merge VTPR
* `#1193 <https://github.com/CoolProp/CoolProp/issues/1193>`_ : REFPROP backend is missing acentric factor accessor
* `#1192 <https://github.com/CoolProp/CoolProp/issues/1192>`_ : Missing formulas for some HFO
* `#1191 <https://github.com/CoolProp/CoolProp/issues/1191>`_ : Linked states need to be updated in copy_k
* `#1190 <https://github.com/CoolProp/CoolProp/issues/1190>`_ : Problems running the VB.NET and C# wrappers
* `#1189 <https://github.com/CoolProp/CoolProp/issues/1189>`_ : Cubic backend broken for PQ calls
* `#1188 <https://github.com/CoolProp/CoolProp/issues/1188>`_ : Critical state not copying for cubics
* `#1187 <https://github.com/CoolProp/CoolProp/issues/1187>`_ : All critical points destroy density solver
* `#1185 <https://github.com/CoolProp/CoolProp/issues/1185>`_ : Add 4th order solver (Halley+)
* `#1184 <https://github.com/CoolProp/CoolProp/issues/1184>`_ : Add 4th order alphar derivatives to python
* `#1183 <https://github.com/CoolProp/CoolProp/issues/1183>`_ : QT/PQ inputs needs to polish with Newton-Raphson
* `#1182 <https://github.com/CoolProp/CoolProp/issues/1182>`_ : Add function to generate rapidjson instance from JSON string
* `#1181 <https://github.com/CoolProp/CoolProp/issues/1181>`_ : Add warning about T > Tmax for HS inputs
* `#1180 <https://github.com/CoolProp/CoolProp/issues/1180>`_ : CoolProp add-in for Excel not working on re-opened files
* `#1179 <https://github.com/CoolProp/CoolProp/issues/1179>`_ : Add derivatives of vr and Tr with respect to beta and gamma
* `#1178 <https://github.com/CoolProp/CoolProp/issues/1178>`_ : Android Wrapper Undefined Reference error with latest ndk
* `#1176 <https://github.com/CoolProp/CoolProp/issues/1176>`_ : [VTPR] mole fractions must be set before calling set_temperature
* `#1175 <https://github.com/CoolProp/CoolProp/issues/1175>`_ : Impose phase for REFPROP in low-level interface
* `#1174 <https://github.com/CoolProp/CoolProp/issues/1174>`_ : Update PHP module docs
* `#1172 <https://github.com/CoolProp/CoolProp/issues/1172>`_ : Please Help With Java Wrapper
* `#1170 <https://github.com/CoolProp/CoolProp/issues/1170>`_ : Incorrect InChI keys
* `#1169 <https://github.com/CoolProp/CoolProp/issues/1169>`_ : Issue with PropsSI on Methane-Ethane mixtures
* `#1168 <https://github.com/CoolProp/CoolProp/issues/1168>`_ : Volume translation for cubic
* `#1166 <https://github.com/CoolProp/CoolProp/issues/1166>`_ : Thermodynamic Properties of R1233zd(E)
* `#1165 <https://github.com/CoolProp/CoolProp/issues/1165>`_ : Not erroring if T < Tmin and p > ptriple
* `#1164 <https://github.com/CoolProp/CoolProp/issues/1164>`_ : REFPROP doesn't store mole fractions in phase envelope
* `#1161 <https://github.com/CoolProp/CoolProp/issues/1161>`_ : [VTPR] gE/RT needs to be multiplied by RT
* `#1158 <https://github.com/CoolProp/CoolProp/issues/1158>`_ : Retrieve phase envelope through high-level DLL
* `#1150 <https://github.com/CoolProp/CoolProp/issues/1150>`_ : IF97 backend: Wrong results for cvmass
* `#1148 <https://github.com/CoolProp/CoolProp/issues/1148>`_ : Add new EOS for MDM of Thol
* `#1146 <https://github.com/CoolProp/CoolProp/issues/1146>`_ : MEXW32 is actually 64-bit and crashes MATLAB
* `#1145 <https://github.com/CoolProp/CoolProp/issues/1145>`_ : Re-implement fundamental derivative of gas dynamics
* `#1144 <https://github.com/CoolProp/CoolProp/issues/1144>`_ : Repair use of spinodals and cubic backend
* `#1143 <https://github.com/CoolProp/CoolProp/issues/1143>`_ : PT inputs for cubics without phase specification
* `#1142 <https://github.com/CoolProp/CoolProp/issues/1142>`_ : PQ inputs very slow for cubic backends
* `#1141 <https://github.com/CoolProp/CoolProp/issues/1141>`_ : dichloroethane has the wrong CAS #
* `#1137 <https://github.com/CoolProp/CoolProp/issues/1137>`_ : Nonsensical results for mistaken inputs with INCOMP fluids
* `#1122 <https://github.com/CoolProp/CoolProp/issues/1122>`_ : Calculate density with PropsSi in Javascript
* `#1120 <https://github.com/CoolProp/CoolProp/issues/1120>`_ : Allow state generation from backend_name() return values
* `#1118 <https://github.com/CoolProp/CoolProp/issues/1118>`_ : Fix plots for cases with multiple critical points
* `#1114 <https://github.com/CoolProp/CoolProp/issues/1114>`_ : Export set_binary_interaction_double + Julia wrapper improvement
* `#1111 <https://github.com/CoolProp/CoolProp/issues/1111>`_ : Improvements to SMath wrapper error handling and some small tweaks
* `#1109 <https://github.com/CoolProp/CoolProp/issues/1109>`_ : SMath wrapper: update AssemblyInfo.cs.template
* `#1108 <https://github.com/CoolProp/CoolProp/issues/1108>`_ : SMath copyright year outdated
* `#1107 <https://github.com/CoolProp/CoolProp/issues/1107>`_ : Allow conditional build of SMath in source tree (fixes #1110)

Pull Requests merged:

* `#1283 <https://github.com/CoolProp/CoolProp/pull/1283>`_ : Add pure fluid check to VTPR
* `#1282 <https://github.com/CoolProp/CoolProp/pull/1282>`_ : Correct typo, see #1270
* `#1272 <https://github.com/CoolProp/CoolProp/pull/1272>`_ : Solves a bug in VTPR
* `#1271 <https://github.com/CoolProp/CoolProp/pull/1271>`_ : Remove possible division by 0, closes #1270
* `#1269 <https://github.com/CoolProp/CoolProp/pull/1269>`_ : SatL and SatV of type VTPR too
* `#1266 <https://github.com/CoolProp/CoolProp/pull/1266>`_ : Small fixes for VTPR
* `#1262 <https://github.com/CoolProp/CoolProp/pull/1262>`_ : Set alpha function in JSON
* `#1261 <https://github.com/CoolProp/CoolProp/pull/1261>`_ : Update CMakeLists.txt
* `#1258 <https://github.com/CoolProp/CoolProp/pull/1258>`_ : Solves a bug with cubic and mixtures
* `#1257 <https://github.com/CoolProp/CoolProp/pull/1257>`_ : Update iPhone compilation docs
* `#1239 <https://github.com/CoolProp/CoolProp/pull/1239>`_ : Shortcut VTPR when pure fluids, solves #1232
* `#1234 <https://github.com/CoolProp/CoolProp/pull/1234>`_ : Add the ability to set limits in Consistency plots
* `#1214 <https://github.com/CoolProp/CoolProp/pull/1214>`_ : Refactor exceptions in CoolPropLib.cpp close #1200
* `#1211 <https://github.com/CoolProp/CoolProp/pull/1211>`_ : Use aii_term and b0_ii from cubic
* `#1209 <https://github.com/CoolProp/CoolProp/pull/1209>`_ : Correct tau derivatives in VTPR
* `#1208 <https://github.com/CoolProp/CoolProp/pull/1208>`_ : Correct derivatives of am and test for VTPR
* `#1196 <https://github.com/CoolProp/CoolProp/pull/1196>`_ : Update CoolPropLib.def
* `#1195 <https://github.com/CoolProp/CoolProp/pull/1195>`_ : Merge VTPR
* `#1114 <https://github.com/CoolProp/CoolProp/pull/1114>`_ : Export set_binary_interaction_double + Julia wrapper improvement
* `#1111 <https://github.com/CoolProp/CoolProp/pull/1111>`_ : Improvements to SMath wrapper error handling and some small tweaks
* `#1109 <https://github.com/CoolProp/CoolProp/pull/1109>`_ : SMath wrapper: update AssemblyInfo.cs.template
* `#1107 <https://github.com/CoolProp/CoolProp/pull/1107>`_ : Allow conditional build of SMath in source tree (fixes #1110)
* `#1103 <https://github.com/CoolProp/CoolProp/pull/1103>`_ : One small tweak to Props1
* `#1101 <https://github.com/CoolProp/CoolProp/pull/1101>`_ : Add error handling to some functions, see #1096
* `#1100 <https://github.com/CoolProp/CoolProp/pull/1100>`_ : Allow cmake properly build SMath wrapper
* `#1097 <https://github.com/CoolProp/CoolProp/pull/1097>`_ : Set error string in get_parameter_information_string() and fix SMath wrapper : fixes #1096
* `#1093 <https://github.com/CoolProp/CoolProp/pull/1093>`_ : Revert part of 763d4ce to solve #1077

6.0.0
-----

New features:

* MathCAD wrapper working again (thanks to Jeff Henning)
* Added binary interaction parameters for more than 400 mixtures 
* Added a cubic backend supporting PR and SRK for some calculations
* Added new non-iterative viscosity model for a few refrigerants (especially R32 and R245fa)
* Implemented EOS for HCl, D4, ethylene oxide, and dichloroethane from M. Thol
* ... and a lot of little bugfixes (see issues)

Contributors to this release:
ibell, jowr, henningjp, bilderbuchi, dinojr, mapipolo, Mol3culo, stefann82, arashsk, pypamart, milesabarr, wahlenkus, saha84, EmiCas, Heathckliff, Tom0310, dizzux, davideziviani, paarfi

Issues Closed:

* `#1056 <https://github.com/CoolProp/CoolProp/issues/1056>`_ : Added "set_reference_state" wrapper for Mathcad and Updated Example Worksheets
* `#1053 <https://github.com/CoolProp/CoolProp/issues/1053>`_ : Align Tmax with REFPROP values
* `#1049 <https://github.com/CoolProp/CoolProp/issues/1049>`_ : apply_simple_mixing_rule should be implemented for HEOS instances
* `#1048 <https://github.com/CoolProp/CoolProp/issues/1048>`_ : Calling set_binary_interaction_double on AbstractState instance has no effect
* `#1047 <https://github.com/CoolProp/CoolProp/issues/1047>`_ : Mathcad Wrapper Updates for CoolProp 5.x and 6
* `#1044 <https://github.com/CoolProp/CoolProp/issues/1044>`_ : Manylinux build integration
* `#1041 <https://github.com/CoolProp/CoolProp/issues/1041>`_ : Fixed Minor MSVC Compiler Warnings
* `#1034 <https://github.com/CoolProp/CoolProp/issues/1034>`_ : Strange behaviour of densities at critical point
* `#1033 <https://github.com/CoolProp/CoolProp/issues/1033>`_ : Python builder issues
* `#1032 <https://github.com/CoolProp/CoolProp/issues/1032>`_ : Rewrite mixture derivatives tests to use new format
* `#1031 <https://github.com/CoolProp/CoolProp/issues/1031>`_ : Fixes STRING conflict between Mathcad library and cppformat
* `#1030 <https://github.com/CoolProp/CoolProp/issues/1030>`_ : Add pass-throughs for testing derivatives
* `#1029 <https://github.com/CoolProp/CoolProp/issues/1029>`_ : Sphinx builder
* `#1028 <https://github.com/CoolProp/CoolProp/issues/1028>`_ : ALTERNATIVE_REFPROP_PATH ignored for predefined mixtures
* `#1026 <https://github.com/CoolProp/CoolProp/issues/1026>`_ : Add REFPROP version to REFPROP comparison script
* `#1025 <https://github.com/CoolProp/CoolProp/issues/1025>`_ : Phase envelopes construction failing for example in docs 
* `#1024 <https://github.com/CoolProp/CoolProp/issues/1024>`_ : VLE calcs failing for SRK & PR backends
* `#1023 <https://github.com/CoolProp/CoolProp/issues/1023>`_ : AbstractState.update fails for mixtures containing specific refrigerants using REFPROP backend
* `#1020 <https://github.com/CoolProp/CoolProp/issues/1020>`_ : Add target_link_libraries to CMakeLists.txt
* `#1014 <https://github.com/CoolProp/CoolProp/issues/1014>`_ : Figure out how to make coolprop static library a clean cmake dependency
* `#1012 <https://github.com/CoolProp/CoolProp/issues/1012>`_ : Residual Helmholtz energy not work
* `#1011 <https://github.com/CoolProp/CoolProp/issues/1011>`_ : Update references
* `#1010 <https://github.com/CoolProp/CoolProp/issues/1010>`_ : Derivative of residual Helmholtz energy with delta
* `#1009 <https://github.com/CoolProp/CoolProp/issues/1009>`_ : Can't compute densities at the triple point
* `#1007 <https://github.com/CoolProp/CoolProp/issues/1007>`_ : 'error: key [Ar] was not found in string_to_index'
* `#1006 <https://github.com/CoolProp/CoolProp/issues/1006>`_ : Use c++14 when building on MINGW
* `#1005 <https://github.com/CoolProp/CoolProp/issues/1005>`_ : Derivative of the saturation enthalpy cair_sat = d(hsat)/dT
* `#1003 <https://github.com/CoolProp/CoolProp/issues/1003>`_ : Fix bug in Chung estimation model
* `#1002 <https://github.com/CoolProp/CoolProp/issues/1002>`_ : Add python 3.5 wheel
* `#1001 <https://github.com/CoolProp/CoolProp/issues/1001>`_ : DmolarP broken for Air
* `#1000 <https://github.com/CoolProp/CoolProp/issues/1000>`_ : Fix setting of BIP function
* `#999 <https://github.com/CoolProp/CoolProp/issues/999>`_ : Abbreviate all journal names
* `#998 <https://github.com/CoolProp/CoolProp/issues/998>`_ : Refine phase envelope better on liquid side
* `#997 <https://github.com/CoolProp/CoolProp/issues/997>`_ : Abbreviate IECR in CoolProp reference
* `#996 <https://github.com/CoolProp/CoolProp/issues/996>`_ : Update references for R245fa and R1234ze(E)
* `#995 <https://github.com/CoolProp/CoolProp/issues/995>`_ : Check double_equal in CPnumerics.h
* `#994 <https://github.com/CoolProp/CoolProp/issues/994>`_ : Find a way to simplify includes
* `#993 <https://github.com/CoolProp/CoolProp/issues/993>`_ : Test/Add example for DLL calling from C
* `#992 <https://github.com/CoolProp/CoolProp/issues/992>`_ : Fix reference for R1234ze(E) again
* `#987 <https://github.com/CoolProp/CoolProp/issues/987>`_ : Multiple EOS paper refs run together
* `#986 <https://github.com/CoolProp/CoolProp/issues/986>`_ : Air lookup in Excel v5.1.2
* `#982 <https://github.com/CoolProp/CoolProp/issues/982>`_ : Reorganize CoolPropTools.h into smaller modules
* `#981 <https://github.com/CoolProp/CoolProp/issues/981>`_ : Saturation states
* `#976 <https://github.com/CoolProp/CoolProp/issues/976>`_ : Add high-level functions to Julia wrapper
* `#975 <https://github.com/CoolProp/CoolProp/issues/975>`_ : Correct get_parameter_information_string, fixes #974
* `#973 <https://github.com/CoolProp/CoolProp/issues/973>`_ : Remove warnings when using Julia 0.4 realease
* `#971 <https://github.com/CoolProp/CoolProp/issues/971>`_ : Fix bug in PhaseEnvelopeRoutines::evaluate
* `#970 <https://github.com/CoolProp/CoolProp/issues/970>`_ : Props1SI function missing in Mathematica wrapper on OSX
* `#968 <https://github.com/CoolProp/CoolProp/issues/968>`_ : Update index.rst
* `#967 <https://github.com/CoolProp/CoolProp/issues/967>`_ : SO2 ancillaries broken
* `#964 <https://github.com/CoolProp/CoolProp/issues/964>`_ : Update index.rst
* `#963 <https://github.com/CoolProp/CoolProp/issues/963>`_ : Update index.rst
* `#962 <https://github.com/CoolProp/CoolProp/issues/962>`_ : Update sample.sce
* `#960 <https://github.com/CoolProp/CoolProp/issues/960>`_ : Update index.rst
* `#953 <https://github.com/CoolProp/CoolProp/issues/953>`_ : Remap CoolPropDbl to double
* `#952 <https://github.com/CoolProp/CoolProp/issues/952>`_ : Switch string formatting to use the cppformat library; see #907
* `#951 <https://github.com/CoolProp/CoolProp/issues/951>`_ : Allow gibbs as input to first_partial_deriv()
* `#950 <https://github.com/CoolProp/CoolProp/issues/950>`_ : Wrong units for residual entropy
* `#949 <https://github.com/CoolProp/CoolProp/issues/949>`_ : Fix {} in bibtex to protect title capitalization
* `#948 <https://github.com/CoolProp/CoolProp/issues/948>`_ : Update reference for  EOS-CG
* `#947 <https://github.com/CoolProp/CoolProp/issues/947>`_ : Add Fij to REFPROPMixtureBackend::get_binary_interaction_double
* `#945 <https://github.com/CoolProp/CoolProp/issues/945>`_ : Add EOS for R245ca
* `#944 <https://github.com/CoolProp/CoolProp/issues/944>`_ : Update reference for R1233ze(E)
* `#941 <https://github.com/CoolProp/CoolProp/issues/941>`_ : CoolProp returns same value for p_critical and p_triple for R503
* `#937 <https://github.com/CoolProp/CoolProp/issues/937>`_ : Allow ability to get refprop version
* `#934 <https://github.com/CoolProp/CoolProp/issues/934>`_ : Memory access violation on mixture update at very low pressures using tabular backend
* `#933 <https://github.com/CoolProp/CoolProp/issues/933>`_ : ValueError: Bad phase to solver_rho_Tp_SRK (CoolProp 5.1.2)
* `#932 <https://github.com/CoolProp/CoolProp/issues/932>`_ : Fix EOS reference for oxygen
* `#931 <https://github.com/CoolProp/CoolProp/issues/931>`_ : Remap CoolPropDbl to double permanently
* `#930 <https://github.com/CoolProp/CoolProp/issues/930>`_ : Phase envelopes should be begin at much lower pressure
* `#929 <https://github.com/CoolProp/CoolProp/issues/929>`_ : PT should start with Halley's method everywhere
* `#928 <https://github.com/CoolProp/CoolProp/issues/928>`_ : Add EOS for HCl, D4, ethylene oxide, and dichloroethane
* `#927 <https://github.com/CoolProp/CoolProp/issues/927>`_ : Add ability to use Henry's Law to get guesses for liquid phase composition
* `#926 <https://github.com/CoolProp/CoolProp/issues/926>`_ : hydrogen formula is wrong
* `#925 <https://github.com/CoolProp/CoolProp/issues/925>`_ : Fix HS inputs 
* `#921 <https://github.com/CoolProp/CoolProp/issues/921>`_ : Tabular calcs with mixtures often return Dew T< Bubble T using PQ input pair
* `#920 <https://github.com/CoolProp/CoolProp/issues/920>`_ : Can't find temperature at pressure and entropy
* `#917 <https://github.com/CoolProp/CoolProp/issues/917>`_ : Fix errors in docs
* `#907 <https://github.com/CoolProp/CoolProp/issues/907>`_ : Replace string formatting with C++ format library
* `#905 <https://github.com/CoolProp/CoolProp/issues/905>`_ : Using conda recipes
* `#885 <https://github.com/CoolProp/CoolProp/issues/885>`_ : Duplicate critical points found 
* `#854 <https://github.com/CoolProp/CoolProp/issues/854>`_ : Coolprop R448A, R449A or R450A
* `#816 <https://github.com/CoolProp/CoolProp/issues/816>`_ : Issue with viscosity of R245FA
* `#808 <https://github.com/CoolProp/CoolProp/issues/808>`_ : Implement tangent plane distance
* `#665 <https://github.com/CoolProp/CoolProp/issues/665>`_ : Viscosity convergence issue
* `#279 <https://github.com/CoolProp/CoolProp/issues/279>`_ : Rebuild MathCAD wrapper with v5 support
* `#186 <https://github.com/CoolProp/CoolProp/issues/186>`_ : Convert cubics to HE

Pull Requests merged:

* `#1062 <https://github.com/CoolProp/CoolProp/pull/1062>`_ : Export first_partial_deriv, see #946 #1062
* `#1056 <https://github.com/CoolProp/CoolProp/pull/1056>`_ : Added "set_reference_state" wrapper for Mathcad and Updated Example Worksheets
* `#1053 <https://github.com/CoolProp/CoolProp/pull/1053>`_ : Align Tmax with REFPROP values
* `#1047 <https://github.com/CoolProp/CoolProp/pull/1047>`_ : Mathcad Wrapper Updates for CoolProp 5.x and 6
* `#1041 <https://github.com/CoolProp/CoolProp/pull/1041>`_ : Fixed Minor MSVC Compiler Warnings
* `#1031 <https://github.com/CoolProp/CoolProp/pull/1031>`_ : Fixes STRING conflict between Mathcad library and cppformat
* `#1020 <https://github.com/CoolProp/CoolProp/pull/1020>`_ : Add target_link_libraries to CMakeLists.txt
* `#982 <https://github.com/CoolProp/CoolProp/pull/982>`_ : Reorganize CoolPropTools.h into smaller modules
* `#981 <https://github.com/CoolProp/CoolProp/pull/981>`_ : Saturation states
* `#976 <https://github.com/CoolProp/CoolProp/pull/976>`_ : Add high-level functions to Julia wrapper
* `#975 <https://github.com/CoolProp/CoolProp/pull/975>`_ : Correct get_parameter_information_string, fixes #974
* `#973 <https://github.com/CoolProp/CoolProp/pull/973>`_ : Remove warnings when using Julia 0.4 realease
* `#968 <https://github.com/CoolProp/CoolProp/pull/968>`_ : Update index.rst
* `#964 <https://github.com/CoolProp/CoolProp/pull/964>`_ : Update index.rst
* `#963 <https://github.com/CoolProp/CoolProp/pull/963>`_ : Update index.rst
* `#962 <https://github.com/CoolProp/CoolProp/pull/962>`_ : Update sample.sce
* `#960 <https://github.com/CoolProp/CoolProp/pull/960>`_ : Update index.rst
* `#953 <https://github.com/CoolProp/CoolProp/pull/953>`_ : Remap CoolPropDbl to double
* `#952 <https://github.com/CoolProp/CoolProp/pull/952>`_ : Switch string formatting to use the cppformat library; see #907

5.1.2
-----

New features:

* Android wrapper available
* Javascript interface extended to export AbstractState and some functions
* Fixed a wide range of issues with tables
* ... and a lot of little bugfixes (see issues)

Issues Closed:

* `#914 <https://github.com/CoolProp/CoolProp/issues/914>`_ : Tabular ammonia calc yields very different results using TTSE vs. bicubic, including non-physical and NaN quantities
* `#909 <https://github.com/CoolProp/CoolProp/issues/909>`_ : Fortran wrapper on Win...still unable to run it!
* `#906 <https://github.com/CoolProp/CoolProp/issues/906>`_ : Add DOI for Novec649
* `#904 <https://github.com/CoolProp/CoolProp/issues/904>`_ : Deuterium reference has wrong year
* `#903 <https://github.com/CoolProp/CoolProp/issues/903>`_ : Some BibTeX keys need updating
* `#902 <https://github.com/CoolProp/CoolProp/issues/902>`_ : Trap errors in get_BibTeXKey and throw
* `#901 <https://github.com/CoolProp/CoolProp/issues/901>`_ : Viscosity of some incompressibles off by a factor of 100 and 1000
* `#899 <https://github.com/CoolProp/CoolProp/issues/899>`_ : Cp, Cv, speed_sound cannot be calculated with QT inputs (Q=0 or 1) and tabular backends
* `#897 <https://github.com/CoolProp/CoolProp/issues/897>`_ : Update DEF for new AbstractState functions
* `#896 <https://github.com/CoolProp/CoolProp/issues/896>`_ : Tabular refactor
* `#894 <https://github.com/CoolProp/CoolProp/issues/894>`_ : License on homepage
* `#889 <https://github.com/CoolProp/CoolProp/issues/889>`_ :  MSVCP100.dll and MSVCR100.dll dependency issue...
* `#888 <https://github.com/CoolProp/CoolProp/issues/888>`_ : Multi-output library function
* `#886 <https://github.com/CoolProp/CoolProp/issues/886>`_ : ALTERNATE_REFPROP_PATH ignored in low-level interface
* `#882 <https://github.com/CoolProp/CoolProp/issues/882>`_ : Tabular backends and phase specification
* `#880 <https://github.com/CoolProp/CoolProp/issues/880>`_ : low-level interface MATLAB using shared library
* `#871 <https://github.com/CoolProp/CoolProp/issues/871>`_ : Issues with Cp, Cv, u, and viscosity with QT_INPUTS where Q=0 or 1 (xxx&REFPROP backend)
* `#869 <https://github.com/CoolProp/CoolProp/issues/869>`_ : Fix javascript builder on buildbot
* `#868 <https://github.com/CoolProp/CoolProp/issues/868>`_ : Fix fortran builds on buildbot
* `#865 <https://github.com/CoolProp/CoolProp/issues/865>`_ : Hide tabular generation outputs when debug_level=0
* `#859 <https://github.com/CoolProp/CoolProp/issues/859>`_ : Windows wrapper for Octave not working for v 4.0
* `#853 <https://github.com/CoolProp/CoolProp/issues/853>`_ : Problem with linking shared libraries using Code::Blocks and CoolProp
* `#849 <https://github.com/CoolProp/CoolProp/issues/849>`_ : Tidy up references in online docs
* `#848 <https://github.com/CoolProp/CoolProp/issues/848>`_ : PropsSImulti in Python
* `#845 <https://github.com/CoolProp/CoolProp/issues/845>`_ : Tabular calculations fail with message "Unable to bisect segmented vector slice..."
* `#844 <https://github.com/CoolProp/CoolProp/issues/844>`_ : failure in calculation enthalpy for water
* `#843 <https://github.com/CoolProp/CoolProp/issues/843>`_ : Calling AbstractState.update() using Dmass_P input pair causes stack overflow in tabular backends
* `#842 <https://github.com/CoolProp/CoolProp/issues/842>`_ : Wrong enthalpy calculation for SES36
* `#841 <https://github.com/CoolProp/CoolProp/issues/841>`_ : R1233zd(E) reference
* `#840 <https://github.com/CoolProp/CoolProp/issues/840>`_ : Failure to calculate any state using input pair QT_INPUTS with backend TTSE&REFPROP
* `#838 <https://github.com/CoolProp/CoolProp/issues/838>`_ : Request: implement a configuration variable to specify directory for tabular interpolation data
* `#837 <https://github.com/CoolProp/CoolProp/issues/837>`_ : Exceptions thrown when getting/setting MAXIMUM_TABLE_DIRECTORY_SIZE_IN_GB configuration setting
* `#835 <https://github.com/CoolProp/CoolProp/issues/835>`_ : Request: CoolProp.AbstractState.first_saturation_deriv wrapped in CoolPropLib.h
* `#831 <https://github.com/CoolProp/CoolProp/issues/831>`_ : Predefined mixtures fail for BICUBIC&REFPROP backend
* `#826 <https://github.com/CoolProp/CoolProp/issues/826>`_ : Unit conversion problem somewhere in Bicubic backend for enthalpy
* `#825 <https://github.com/CoolProp/CoolProp/issues/825>`_ : PQ_with_guesses assumes bubble point
* `#824 <https://github.com/CoolProp/CoolProp/issues/824>`_ : C-Sharp Wrapper AbstractState mole_fractions_liquid
* `#823 <https://github.com/CoolProp/CoolProp/issues/823>`_ : Documentation for use of static libraries is unclear
* `#822 <https://github.com/CoolProp/CoolProp/issues/822>`_ : Request: PropsSI Inputs of D and Q
* `#821 <https://github.com/CoolProp/CoolProp/issues/821>`_ : Fix pip command for nightly
* `#820 <https://github.com/CoolProp/CoolProp/issues/820>`_ : Add cmake option to generate Android .so library
* `#819 <https://github.com/CoolProp/CoolProp/issues/819>`_ : Expose phase envelope calculations in javascript
* `#814 <https://github.com/CoolProp/CoolProp/issues/814>`_ : saturated_liquid/vapor_keyed_output for tabular backend
* `#812 <https://github.com/CoolProp/CoolProp/issues/812>`_ : Add ability to retrieve mass fractions
* `#810 <https://github.com/CoolProp/CoolProp/issues/810>`_ : Python builds crash on Windows
* `#809 <https://github.com/CoolProp/CoolProp/issues/809>`_ : Implement fluid_param_string in python
* `#807 <https://github.com/CoolProp/CoolProp/issues/807>`_ : Return all critical points
* `#805 <https://github.com/CoolProp/CoolProp/issues/805>`_ : Coolprop function like Refprop Excel Fluidstring Function for mixtures
* `#804 <https://github.com/CoolProp/CoolProp/issues/804>`_ : Allow disabling parameter estimation in REFPROP
* `#802 <https://github.com/CoolProp/CoolProp/issues/802>`_ : Error with two-phase DT inputs for R134a
* `#800 <https://github.com/CoolProp/CoolProp/issues/800>`_ : Add access to contributions to viscosity and conductivity
* `#799 <https://github.com/CoolProp/CoolProp/issues/799>`_ : Add access to conformal state solver in AbstractState
* `#798 <https://github.com/CoolProp/CoolProp/issues/798>`_ : Add linear and Lorentz-Berthelot mixing rules
* `#796 <https://github.com/CoolProp/CoolProp/issues/796>`_ : Add SATTP guess implementation
* `#795 <https://github.com/CoolProp/CoolProp/issues/795>`_ : Provide swigged MATLAB wrapper code
* `#793 <https://github.com/CoolProp/CoolProp/issues/793>`_ : Set interaction parameters in REFPROP through CoolProp
* `#792 <https://github.com/CoolProp/CoolProp/issues/792>`_ : Allow possibility to set interaction parameters even if the mixture isn't already included
* `#789 <https://github.com/CoolProp/CoolProp/issues/789>`_ : Make sure all phases are calculated correctly for BICUBIC&HEOS backend
* `#788 <https://github.com/CoolProp/CoolProp/issues/788>`_ : Make sure all phases are calculated correctly for HEOS backend
* `#786 <https://github.com/CoolProp/CoolProp/issues/786>`_ : Implement conductivity for pentanes
* `#785 <https://github.com/CoolProp/CoolProp/issues/785>`_ : Implement viscosity for Toluene
* `#784 <https://github.com/CoolProp/CoolProp/issues/784>`_ : Add docs for get/set config functions
* `#783 <https://github.com/CoolProp/CoolProp/issues/783>`_ : Failure in PsychScript
* `#777 <https://github.com/CoolProp/CoolProp/issues/777>`_ : No input passed with PT_INPUTS and tabular backed
* `#776 <https://github.com/CoolProp/CoolProp/issues/776>`_ : Fix docs for IF97 backend
* `#773 <https://github.com/CoolProp/CoolProp/issues/773>`_ : Missing files in LabVIEW wrapper folder or documentation needed
* `#772 <https://github.com/CoolProp/CoolProp/issues/772>`_ : Acentric factor of air
* `#770 <https://github.com/CoolProp/CoolProp/issues/770>`_ : Make clear() overridable / clear Helmholtz cache
* `#769 <https://github.com/CoolProp/CoolProp/issues/769>`_ : Improve docs for second partial derivatives
* `#768 <https://github.com/CoolProp/CoolProp/issues/768>`_ : Fix solver for first criticality contour crossing
* `#767 <https://github.com/CoolProp/CoolProp/issues/767>`_ : When tracing criticality contour, make sure that delta is always increasing
* `#764 <https://github.com/CoolProp/CoolProp/issues/764>`_ : Add `calc_speed_sound` to tabular backend
* `#763 <https://github.com/CoolProp/CoolProp/issues/763>`_ : Add and implement all phase functions to tabular backends
* `#762 <https://github.com/CoolProp/CoolProp/issues/762>`_ : Temperature with `HmassP_INPUTS` with twophase fluid and tabular
* `#761 <https://github.com/CoolProp/CoolProp/issues/761>`_ : Add auto-generated docs for configuration variables
* `#760 <https://github.com/CoolProp/CoolProp/issues/760>`_ : Add `surface tension` to tabular backend
* `#759 <https://github.com/CoolProp/CoolProp/issues/759>`_ : Add comprehensive docs for REFPROP interface
* `#757 <https://github.com/CoolProp/CoolProp/issues/757>`_ : Cannot evaluate PT (or PH?) below p_triple
* `#756 <https://github.com/CoolProp/CoolProp/issues/756>`_ : HAPropsSI does not converge for T= 299.8 K
* `#754 <https://github.com/CoolProp/CoolProp/issues/754>`_ : Failure with sat derivative with QT and tables
* `#753 <https://github.com/CoolProp/CoolProp/issues/753>`_ : Relative humidity calculation error
* `#751 <https://github.com/CoolProp/CoolProp/issues/751>`_ : D-P is far slower than it should be
* `#750 <https://github.com/CoolProp/CoolProp/issues/750>`_ : Invalid index to calc_first_saturation_deriv in TabularBackends
* `#747 <https://github.com/CoolProp/CoolProp/issues/747>`_ : Plotting example on coolprop.org does not work - potentially related to issue #351
* `#746 <https://github.com/CoolProp/CoolProp/issues/746>`_ : Implement viscosity models for HFO (ECS?)
* `#745 <https://github.com/CoolProp/CoolProp/issues/745>`_ : Undocumented high level interface for saturation derivatives
* `#742 <https://github.com/CoolProp/CoolProp/issues/742>`_ : Expedite the D+Y flash routines
* `#741 <https://github.com/CoolProp/CoolProp/issues/741>`_ : Expedite the single-phase T+Y flash routines
* `#740 <https://github.com/CoolProp/CoolProp/issues/740>`_ : HapropsSI("T", "B", 299.15, "R", 0, "P", 101325) lead to an error
* `#739 <https://github.com/CoolProp/CoolProp/issues/739>`_ : Quality-related updates with tabular backend
* `#738 <https://github.com/CoolProp/CoolProp/issues/738>`_ : TTSE ranges
* `#737 <https://github.com/CoolProp/CoolProp/issues/737>`_ : Missing bib entry IAPWS-SurfaceTension-1994
* `#735 <https://github.com/CoolProp/CoolProp/issues/735>`_ : phase is wrong for water at STP
* `#734 <https://github.com/CoolProp/CoolProp/issues/734>`_ : F is missing from mixture interaction parameters on the web
* `#733 <https://github.com/CoolProp/CoolProp/issues/733>`_ : Typo in excess term in mixture docs
* `#731 <https://github.com/CoolProp/CoolProp/issues/731>`_ : Add EOS for Novec 649 from McLinden
* `#730 <https://github.com/CoolProp/CoolProp/issues/730>`_ : Merge references from paper about CoolProp into main bib file
* `#727 <https://github.com/CoolProp/CoolProp/issues/727>`_ : HapropsSI("T", "B", 299.15, "R", 0, "P", 101325) lead to an error
* `#726 <https://github.com/CoolProp/CoolProp/issues/726>`_ : Improve caching of derivative terms when using mixtures
* `#725 <https://github.com/CoolProp/CoolProp/issues/725>`_ : Implement dipole moment

5.1.1
-----

New features:

* A wrapper for the R language
* Tabular integration with tables from REFPROP only for now
* The Python wrapper is now also available on binstar: https://anaconda.org/CoolProp/coolprop
* ... and a lot of little bugfixes (see issues)

Issues Closed:

* `#724 <https://github.com/CoolProp/CoolProp/issues/724>`_ : Gibbs not working as output (mass or molar)
* `#722 <https://github.com/CoolProp/CoolProp/issues/722>`_ : Predefined mixtures crash python
* `#721 <https://github.com/CoolProp/CoolProp/issues/721>`_ : v5.1.1
* `#714 <https://github.com/CoolProp/CoolProp/issues/714>`_ : Possible error in isobaric thermal expansion coefficient
* `#713 <https://github.com/CoolProp/CoolProp/issues/713>`_ : Bicubic backend and first_saturation_deriv
* `#712 <https://github.com/CoolProp/CoolProp/issues/712>`_ : Expose saturation derivatives from PropsSI [wishlist]
* `#708 <https://github.com/CoolProp/CoolProp/issues/708>`_ : CoolPropsetup.m needs to be installed
* `#707 <https://github.com/CoolProp/CoolProp/issues/707>`_ : conda builds
* `#703 <https://github.com/CoolProp/CoolProp/issues/703>`_ : 2/ HapropsSI ( "T" , "B" , ValueB, "W" , 0 , "P" , 101325) lead to an error
* `#702 <https://github.com/CoolProp/CoolProp/issues/702>`_ : 1 : HapropsSI ( "T" , "H" , ValueH, "W" , 0 , "P" , 101325) lead to an error
* `#700 <https://github.com/CoolProp/CoolProp/issues/700>`_ : If git is not found, still compile properly
* `#699 <https://github.com/CoolProp/CoolProp/issues/699>`_ : Fugacity using Python wrapper
* `#697 <https://github.com/CoolProp/CoolProp/issues/697>`_ : Get State (old-style) class working with predefined mixtures
* `#696 <https://github.com/CoolProp/CoolProp/issues/696>`_ : cp0 broken for tabular backends
* `#695 <https://github.com/CoolProp/CoolProp/issues/695>`_ : Problem with reference state
* `#691 <https://github.com/CoolProp/CoolProp/issues/691>`_ : variable names for second_partial_deriv
* `#688 <https://github.com/CoolProp/CoolProp/issues/688>`_ : PropsSI in saturation region
* `#685 <https://github.com/CoolProp/CoolProp/issues/685>`_ : Problem with Hazard output
* `#684 <https://github.com/CoolProp/CoolProp/issues/684>`_ : some problem and questions for calc in Excel
* `#681 <https://github.com/CoolProp/CoolProp/issues/681>`_ : Mix call failure after release update
* `#680 <https://github.com/CoolProp/CoolProp/issues/680>`_ : Tabular backend data range too small for (P,H) inputs and R245fa
* `#675 <https://github.com/CoolProp/CoolProp/issues/675>`_ : Get consistency plots working with Tabular backends
* `#674 <https://github.com/CoolProp/CoolProp/issues/674>`_ : QT inputs do not work for Tabular backends
* `#673 <https://github.com/CoolProp/CoolProp/issues/673>`_ : Mass-based saturation derivatives not supported
* `#672 <https://github.com/CoolProp/CoolProp/issues/672>`_ : Tabular methods returns hmolar for smolar for saturation
* `#671 <https://github.com/CoolProp/CoolProp/issues/671>`_ : MATLAB on OSX cannot load REFPROP
* `#670 <https://github.com/CoolProp/CoolProp/issues/670>`_ : Low-Level interfacing with MATLAB
* `#668 <https://github.com/CoolProp/CoolProp/issues/668>`_ : R wrapper
* `#664 <https://github.com/CoolProp/CoolProp/issues/664>`_ : Re-enable triple point for REFPROP backend for mixtures
* `#663 <https://github.com/CoolProp/CoolProp/issues/663>`_ : Vapor mass quality = 1 generates error for pseudo-pures
* `#662 <https://github.com/CoolProp/CoolProp/issues/662>`_ : Write function to determine phase after an update with PT and a guess for rho
* `#661 <https://github.com/CoolProp/CoolProp/issues/661>`_ : Predefined mixtures not working properly with Tabular backends
* `#660 <https://github.com/CoolProp/CoolProp/issues/660>`_ : T,X and PS, PD, PU not working with BICUBIC, but does with TTSE
* `#659 <https://github.com/CoolProp/CoolProp/issues/659>`_ : Add "PIP" as parameter
* `#658 <https://github.com/CoolProp/CoolProp/issues/658>`_ : Implement PIP for REFPROP
* `#657 <https://github.com/CoolProp/CoolProp/issues/657>`_ : Describe how to call REFPROP
* `#654 <https://github.com/CoolProp/CoolProp/issues/654>`_ : Add ability to calculate Ideal curves
* `#653 <https://github.com/CoolProp/CoolProp/issues/653>`_ : Implement update_with_guesses for P,T for REFPROP backend
* `#652 <https://github.com/CoolProp/CoolProp/issues/652>`_ : Implement solver for "true" critical point using REFPROP
* `#650 <https://github.com/CoolProp/CoolProp/issues/650>`_ : MATLAB examples not on website
* `#648 <https://github.com/CoolProp/CoolProp/issues/648>`_ : Link to examples broken
* `#647 <https://github.com/CoolProp/CoolProp/issues/647>`_ : Implement the new REFPROP header file and make necessary changes
* `#646 <https://github.com/CoolProp/CoolProp/issues/646>`_ : Add B,C virial coefficients for REFPROP backend
* `#645 <https://github.com/CoolProp/CoolProp/issues/645>`_ : PQ_INPUTS don't work with TTSE backend
* `#644 <https://github.com/CoolProp/CoolProp/issues/644>`_ : Get first_two_phase_deriv working with Tabular backends
* `#641 <https://github.com/CoolProp/CoolProp/issues/641>`_ : Install psyrc file
* `#640 <https://github.com/CoolProp/CoolProp/issues/640>`_ : Expose saturation_ancillary equation through python
* `#639 <https://github.com/CoolProp/CoolProp/issues/639>`_ : Incorrect error when non two-phase inputs to two-phase deriv
* `#638 <https://github.com/CoolProp/CoolProp/issues/638>`_ : Heavy Water Viscosity Unavailable
* `#636 <https://github.com/CoolProp/CoolProp/issues/636>`_ : Error surface tension in CoolProp v5.1.0
* `#635 <https://github.com/CoolProp/CoolProp/issues/635>`_ : Implement first_saturation_deriv for TTSE/BICUBIC
* `#631 <https://github.com/CoolProp/CoolProp/issues/631>`_ : Methane conductivity
* `#630 <https://github.com/CoolProp/CoolProp/issues/630>`_ : Make HS use DH rather than PH
* `#629 <https://github.com/CoolProp/CoolProp/issues/629>`_ : Handle PT inputs around saturation in a better way with BICUBIC
* `#628 <https://github.com/CoolProp/CoolProp/issues/628>`_ : Dry air enthalpy
* `#627 <https://github.com/CoolProp/CoolProp/issues/627>`_ : Test that H and S are the same for all the state points
* `#626 <https://github.com/CoolProp/CoolProp/issues/626>`_ : Improve docs for low-level interface
* `#622 <https://github.com/CoolProp/CoolProp/issues/622>`_ : TTSE fails around saturated liquid
* `#617 <https://github.com/CoolProp/CoolProp/issues/617>`_ : Block Tabular backend use with PropsSI somehow

5.1.0
-----

New features:

* Tabular interpolation using TTSE or Bicubic interpolation (https://coolprop.org/coolprop/Tabular.html)
* Equation of state for heavy water
* Added IF97 backend for industrial formulation for properties of pure water
* Lots of little bugfixes (see issues)

Issues Closed:

* `#624 <https://github.com/CoolProp/CoolProp/issues/624>`_ : Stability in two-phase region
* `#621 <https://github.com/CoolProp/CoolProp/issues/621>`_ : TTSE Input Param (Water)
* `#620 <https://github.com/CoolProp/CoolProp/issues/620>`_ : TTSE Problem (Water)
* `#618 <https://github.com/CoolProp/CoolProp/issues/618>`_ : H,S not working for pseudo-pure
* `#615 <https://github.com/CoolProp/CoolProp/issues/615>`_ : Ammonia T-P saturation calculation deviation
* `#614 <https://github.com/CoolProp/CoolProp/issues/614>`_ : Typos in parameter descriptions.
* `#612 <https://github.com/CoolProp/CoolProp/issues/612>`_ : Added missing cell "Input/Output" for enthalpy row.
* `#611 <https://github.com/CoolProp/CoolProp/issues/611>`_ : Splined Output Doubt
* `#609 <https://github.com/CoolProp/CoolProp/issues/609>`_ : Some Windows builds fail (error removing non-existent directory)
* `#608 <https://github.com/CoolProp/CoolProp/issues/608>`_ : MinGW builds fail
* `#605 <https://github.com/CoolProp/CoolProp/issues/605>`_ : CMake changes
* `#602 <https://github.com/CoolProp/CoolProp/issues/602>`_ : TTSE fails for two-phase H,P with heavy water
* `#601 <https://github.com/CoolProp/CoolProp/issues/601>`_ : Benzene conductivity bibtex is wrong
* `#599 <https://github.com/CoolProp/CoolProp/issues/599>`_ : Something is messed up with water properties
* `#595 <https://github.com/CoolProp/CoolProp/issues/595>`_ : add DOIs to bibliography
* `#591 <https://github.com/CoolProp/CoolProp/issues/591>`_ : Request for extension: table of quantities in the documentation for HAPropsSI like for PropsSI
* `#588 <https://github.com/CoolProp/CoolProp/issues/588>`_ : matplotlib and numpy should not be explicit dependencies
* `#586 <https://github.com/CoolProp/CoolProp/issues/586>`_ : HAProps humidity ratio calculation issue
* `#585 <https://github.com/CoolProp/CoolProp/issues/585>`_ : HAProps at low humidity ratio
* `#584 <https://github.com/CoolProp/CoolProp/issues/584>`_ : [Tabular] pure fluid AbstractState returns the wrong mole fractions
* `#583 <https://github.com/CoolProp/CoolProp/issues/583>`_ : Development docs only available on dreamhosters
* `#579 <https://github.com/CoolProp/CoolProp/issues/579>`_ : Issue with Excel Wrapper for Coolprop for OS X Excel 2011
* `#578 <https://github.com/CoolProp/CoolProp/issues/578>`_ : Update examples to show how to call TTSE and BICUBIC backends
* `#577 <https://github.com/CoolProp/CoolProp/issues/577>`_ : Unicode characters in bibtex not being escaped properly
* `#575 <https://github.com/CoolProp/CoolProp/issues/575>`_ : Phase envelopes should be able to be constructed for pure fluids too
* `#574 <https://github.com/CoolProp/CoolProp/issues/574>`_ : Methane (and pentane) transport properties
* `#573 <https://github.com/CoolProp/CoolProp/issues/573>`_ : Bug in derivatives from Matlab
* `#570 <https://github.com/CoolProp/CoolProp/issues/570>`_ : Implement EOS for heavy water
* `#569 <https://github.com/CoolProp/CoolProp/issues/569>`_ : REFPROP SPLNval for rhomolar_vap wrong
* `#568 <https://github.com/CoolProp/CoolProp/issues/568>`_ : Reference of state not working for Refprop backend
* `#567 <https://github.com/CoolProp/CoolProp/issues/567>`_ : Add IF97 Backend
* `#566 <https://github.com/CoolProp/CoolProp/issues/566>`_ : Retrieve phase envelopes from REFPROP using SPLNVAL function
* `#564 <https://github.com/CoolProp/CoolProp/issues/564>`_ : Molecular Formulas as Trivial Property
* `#562 <https://github.com/CoolProp/CoolProp/issues/562>`_ : Add docs about how to set the reference state
* `#556 <https://github.com/CoolProp/CoolProp/issues/556>`_ : [Tabular] Saturation curves for mixtures
* `#555 <https://github.com/CoolProp/CoolProp/issues/555>`_ : [Tabular] Re-enable the PHI0dll function for REFPROP
* `#552 <https://github.com/CoolProp/CoolProp/issues/552>`_ : IsFluidType function
* `#549 <https://github.com/CoolProp/CoolProp/issues/549>`_ : Implement up to 4th order derivatives of all Helmholtz terms (except SAFT)
* `#548 <https://github.com/CoolProp/CoolProp/issues/548>`_ : Problem with HAPropsSI
* `#546 <https://github.com/CoolProp/CoolProp/issues/546>`_ : Small speed enhancement for Julia wrapper
* `#541 <https://github.com/CoolProp/CoolProp/issues/541>`_ : Update CoolProp.jl
* `#540 <https://github.com/CoolProp/CoolProp/issues/540>`_ : Update CoolProp.jl
* `#539 <https://github.com/CoolProp/CoolProp/issues/539>`_ : Add SATTP to REFPROP wrapper
* `#537 <https://github.com/CoolProp/CoolProp/issues/537>`_ : [Tabular] rebuild tables if limits (especially enthalpies) have shifted
* `#536 <https://github.com/CoolProp/CoolProp/issues/536>`_ : Add low level interface to Julia wrapper as discussed in #534 + Fixes #497
* `#535 <https://github.com/CoolProp/CoolProp/issues/535>`_ : When using high-level wrapper of low-level interface, errors don't bubble properly
* `#534 <https://github.com/CoolProp/CoolProp/issues/534>`_ : Add error handling to Julia's wrapper
* `#532 <https://github.com/CoolProp/CoolProp/issues/532>`_ : More Coverity cleanups
* `#530 <https://github.com/CoolProp/CoolProp/issues/530>`_ : When reference state is changed, reducing/critical and hs_anchor states need to be changed
* `#529 <https://github.com/CoolProp/CoolProp/issues/529>`_ : First bunch of Coverity Scan static analysis warning fixes
* `#528 <https://github.com/CoolProp/CoolProp/issues/528>`_ : PQ Flash Failure for CO2+Water
* `#527 <https://github.com/CoolProp/CoolProp/issues/527>`_ : Silence all output to screen when building phase envelopes
* `#526 <https://github.com/CoolProp/CoolProp/issues/526>`_ : When building phase envelopes, stop when the composition is almost pure
* `#524 <https://github.com/CoolProp/CoolProp/issues/524>`_ : set_reference_state does not create expected output
* `#523 <https://github.com/CoolProp/CoolProp/issues/523>`_ : error: thermal conductivity R32:  _phase is unknown
* `#522 <https://github.com/CoolProp/CoolProp/issues/522>`_ : [Tabular] Implement solver when one of the inputs is not a native input
* `#521 <https://github.com/CoolProp/CoolProp/issues/521>`_ : [Tabular] Fix derivatives, and c_p
* `#520 <https://github.com/CoolProp/CoolProp/issues/520>`_ : [Tabular] Fix transport properties
* `#519 <https://github.com/CoolProp/CoolProp/issues/519>`_ : [Tabular] Fix cells close to the saturation curves
* `#518 <https://github.com/CoolProp/CoolProp/issues/518>`_ : Tabular methods implemented
* `#517 <https://github.com/CoolProp/CoolProp/issues/517>`_ : Isobaric expansion coefficient is not implemented
* `#516 <https://github.com/CoolProp/CoolProp/issues/516>`_ : [Tabular] Actually zip up the tables using zlib
* `#515 <https://github.com/CoolProp/CoolProp/issues/515>`_ : Kill off the CRT deprecate warning (#512)
* `#513 <https://github.com/CoolProp/CoolProp/issues/513>`_ : Primitive structures simplification attempt 2
* `#512 <https://github.com/CoolProp/CoolProp/issues/512>`_ : Kill off the CRT deprecate warning
* `#511 <https://github.com/CoolProp/CoolProp/issues/511>`_ : Python version should be 5.1.0dev, not just 5.1.0
* `#508 <https://github.com/CoolProp/CoolProp/issues/508>`_ : Add a ways of using the shared_ptr directly through shared library
* `#507 <https://github.com/CoolProp/CoolProp/issues/507>`_ : Add possibility to disable a backend at compile-time
* `#506 <https://github.com/CoolProp/CoolProp/issues/506>`_ : [Tabular] Add docs for TTSE and bicubic usage
* `#497 <https://github.com/CoolProp/CoolProp/issues/497>`_ : Julia and C++ Low Level Interface for faster Computation
* `#490 <https://github.com/CoolProp/CoolProp/issues/490>`_ : Add partial pressure of water as an output in HAPropsSI
* `#481 <https://github.com/CoolProp/CoolProp/issues/481>`_ : A bug is found when pressure approximates Critical Pressure for Air
* `#455 <https://github.com/CoolProp/CoolProp/issues/455>`_ : HS Inputs in PropsSI function working in two-phase region?
* `#297 <https://github.com/CoolProp/CoolProp/issues/297>`_ : Call matlab script from command line, with no window, catching errors, and never going interactive
* `#296 <https://github.com/CoolProp/CoolProp/issues/296>`_ : Update examples for v5
* `#262 <https://github.com/CoolProp/CoolProp/issues/262>`_ : Re-implement tabular methods
* `#43 <https://github.com/CoolProp/CoolProp/issues/43>`_ : [Tabular] Warn about tabular folder size

5.0.8
-----

New features:

* Added a Smath Studio native wrapper (thanks to Mike Kaganski for all his help)
* Lots of little cleanups to the code (thanks to Mike Kaganski)

Issues Closed:

* `#510 <https://github.com/CoolProp/CoolProp/issues/510>`_ : const, ref and iterator optimization
* `#509 <https://github.com/CoolProp/CoolProp/issues/509>`_ : Exceptions restructured
* `#505 <https://github.com/CoolProp/CoolProp/issues/505>`_ : AbstractState in python should implement phase() function
* `#504 <https://github.com/CoolProp/CoolProp/issues/504>`_ : More ref args
* `#503 <https://github.com/CoolProp/CoolProp/issues/503>`_ : Add compressibility factor for humid air
* `#502 <https://github.com/CoolProp/CoolProp/issues/502>`_ : thread_local broken on OSX
* `#501 <https://github.com/CoolProp/CoolProp/issues/501>`_ : thread_local: one more (hopefully portable) attempt
* `#500 <https://github.com/CoolProp/CoolProp/issues/500>`_ : Fix directory size calculations
* `#499 <https://github.com/CoolProp/CoolProp/issues/499>`_ : Longdouble remap
* `#498 <https://github.com/CoolProp/CoolProp/issues/498>`_ : HAProp - Conductivity & Viscosity
* `#496 <https://github.com/CoolProp/CoolProp/issues/496>`_ : Implement checking of directory size
* `#495 <https://github.com/CoolProp/CoolProp/issues/495>`_ : CoolPropDbl
* `#493 <https://github.com/CoolProp/CoolProp/issues/493>`_ : Avoid copying of parameters; some fixes for _HAPropsSI_inputs
* `#492 <https://github.com/CoolProp/CoolProp/issues/492>`_ : Add docs for Low-Level Interface
* `#488 <https://github.com/CoolProp/CoolProp/issues/488>`_ : Some more static analyser warning fixes
* `#487 <https://github.com/CoolProp/CoolProp/issues/487>`_ : Cannot use REFPROP to get reducing state variables
* `#485 <https://github.com/CoolProp/CoolProp/issues/485>`_ : Rewrite HAPropsSI to call _HAPropsSI
* `#484 <https://github.com/CoolProp/CoolProp/issues/484>`_ : Kill off all warnings in 64-bit compilation
* `#483 <https://github.com/CoolProp/CoolProp/issues/483>`_ : Problems noted by VS2013 static analysis
* `#479 <https://github.com/CoolProp/CoolProp/issues/479>`_ : RelativeHumidity simplification
* `#478 <https://github.com/CoolProp/CoolProp/issues/478>`_ : Julia 0.3 wrapper
* `#476 <https://github.com/CoolProp/CoolProp/issues/476>`_ : buildbot failure messages don't have the correct URL
* `#473 <https://github.com/CoolProp/CoolProp/issues/473>`_ : Wrapper for Julia 0.3
* `#472 <https://github.com/CoolProp/CoolProp/issues/472>`_ : Fix potential buffer overflow with get_parameter_information_string
* `#471 <https://github.com/CoolProp/CoolProp/issues/471>`_ : Document which inputs are possible in Props1SI
* `#470 <https://github.com/CoolProp/CoolProp/issues/470>`_ : Consider evaluating water at Tdb,p for transport properties in humid air
* `#469 <https://github.com/CoolProp/CoolProp/issues/469>`_ : Initialize fluids in HAProps_Aux
* `#468 <https://github.com/CoolProp/CoolProp/issues/468>`_ : Sanitize internal code in HAPropsSI
* `#467 <https://github.com/CoolProp/CoolProp/issues/467>`_ : Cp in HAPropsSI cannot be calculated in 5.0.7
* `#466 <https://github.com/CoolProp/CoolProp/issues/466>`_ : Prandtl number cannot be returned directly


5.0.7
-----

New Features:

* Added a Lua wrapper

Issues Closed:

* `#460 <https://github.com/CoolProp/CoolProp/issues/460>`_ : PropsSI ("Q", "P", valueP, "H", valueH, "REFPROP-R410A") only return 0
* `#459 <https://github.com/CoolProp/CoolProp/issues/459>`_ : PropsSI ("D", "P", valueP, "T", valueT, "R407C") return bad result in L+V Phasis
* `#456 <https://github.com/CoolProp/CoolProp/issues/456>`_ : Slave alert
* `#454 <https://github.com/CoolProp/CoolProp/issues/454>`_ : Add density dependency to entropy and enthalpy of incomprerssible fluids
* `#452 <https://github.com/CoolProp/CoolProp/issues/452>`_ : Allow mixtures to have zero mole fractions
* `#450 <https://github.com/CoolProp/CoolProp/issues/450>`_ : Calling PropsSI to get thermal conductivity throws an exception
* `#448 <https://github.com/CoolProp/CoolProp/issues/448>`_ : Retrieving acentric factor through Props1SI fails
* `#443 <https://github.com/CoolProp/CoolProp/issues/443>`_ : Javascript index.html is missing
* `#437 <https://github.com/CoolProp/CoolProp/issues/437>`_ : REFPROP predefined mixtures no longer work
* `#434 <https://github.com/CoolProp/CoolProp/issues/434>`_ : R404A Refprop value differs from Refprop Value in Excel
* `#432 <https://github.com/CoolProp/CoolProp/issues/432>`_ : All the mixture interaction parameters of Gernert are wrong
* `#431 <https://github.com/CoolProp/CoolProp/issues/431>`_ : REFPROP should not be reloaded after every call to PropsSI
* `#430 <https://github.com/CoolProp/CoolProp/issues/430>`_ : HAPropsSI is missing from the SWIG wrapper
* `#429 <https://github.com/CoolProp/CoolProp/issues/429>`_ : Entropy of Melinder fluids giving wrong results
* `#428 <https://github.com/CoolProp/CoolProp/issues/428>`_ : On windows, do not error out if REFPROP fluid files are not found in c:\Program Files\REFPROP
* `#427 <https://github.com/CoolProp/CoolProp/issues/427>`_ : HapropsSi("W","B", 279.15, "T", 293.15, "P", 101325) lead to a "-1.#IND" value
* `#425 <https://github.com/CoolProp/CoolProp/issues/425>`_ : Incompressible viscosity
* `#419 <https://github.com/CoolProp/CoolProp/issues/419>`_ : HapropSI ("T","B",273.15+37,"D",273.15+36.44,"P",101325) lead to an error ...
* `#416 <https://github.com/CoolProp/CoolProp/issues/416>`_ : Sphinx docs
* `#413 <https://github.com/CoolProp/CoolProp/issues/413>`_ : Incompressible entropy
* `#410 <https://github.com/CoolProp/CoolProp/issues/410>`_ : Phase envelope segfaults for pure fluids
* `#409 <https://github.com/CoolProp/CoolProp/issues/409>`_ : Trivial outputs
* `#408 <https://github.com/CoolProp/CoolProp/issues/408>`_ : HapropsSI function issues
* `#403 <https://github.com/CoolProp/CoolProp/issues/403>`_ : Error in new CoolProp version in the function HAPropsSI (variable combination 'PH' and 'W')
* `#401 <https://github.com/CoolProp/CoolProp/issues/401>`_ : Linux/OSX error with refprop 9.1* and mixtures containing  R1234YF
* `#400 <https://github.com/CoolProp/CoolProp/issues/400>`_ : HAPropsSI(Output, "B",valueB, "R", 1, "P", 101325) lead to an error
* `#398 <https://github.com/CoolProp/CoolProp/issues/398>`_ : HAPropsSI(Output, "B",252.84, "D";250.85, "P", 101325) lead to an infinite value
* `#387 <https://github.com/CoolProp/CoolProp/issues/387>`_ : Vectorised PropSI breaks plotting functions
* `#386 <https://github.com/CoolProp/CoolProp/issues/386>`_ : Bibtex numbering
* `#307 <https://github.com/CoolProp/CoolProp/issues/307>`_ : Transport Properties for Mixtures


5.0.6
-----

New Features:

* Mathematica wrapper finished

Issues Closed:

* `#396 <https://github.com/CoolProp/CoolProp/issues/396>`_ : Initialize fail for HEOS in mixture with Argon and CarbonDioxide (in Matlab)
* `#395 <https://github.com/CoolProp/CoolProp/issues/395>`_ : keyed_output and incompressibles
* `#394 <https://github.com/CoolProp/CoolProp/issues/394>`_ : Python list inputs
* `#391 <https://github.com/CoolProp/CoolProp/issues/391>`_ : release.bsh and source file
* `#390 <https://github.com/CoolProp/CoolProp/issues/390>`_ : Transport properties of water
* `#389 <https://github.com/CoolProp/CoolProp/issues/389>`_ : HAPropsSI("D", "T",273.15+20, "R", 0.8, "P", 101325) lead to an error
* `#384 <https://github.com/CoolProp/CoolProp/issues/384>`_ : Put the example.nb Mathematica file in the main folder
* `#383 <https://github.com/CoolProp/CoolProp/issues/383>`_ : When doing release, force a full build of the docs
* `#382 <https://github.com/CoolProp/CoolProp/issues/382>`_ : Fix up the mathematica docs
* `#379 <https://github.com/CoolProp/CoolProp/issues/379>`_ : After a release is done, delete the release folder
* `#378 <https://github.com/CoolProp/CoolProp/issues/378>`_ : Also integrate the sphinx docs into the binaries/release/unstable folder output
* `#377 <https://github.com/CoolProp/CoolProp/issues/377>`_ : Remove old mathematica files
* `#376 <https://github.com/CoolProp/CoolProp/issues/376>`_ : Add python to list of prerequisites for self-compilation in the docs
* `#329 <https://github.com/CoolProp/CoolProp/issues/329>`_ : Configure buildbot to send emails when we break things

5.0.5
-----

New Features:

* Added Mathematica wrapper
* Added ``Prandtl()`` function to ``AbstractState``
* Added vectorized ``PropsSImulti`` function that can return a matrix of outputs for vectors of state inputs and desired outputs

Removed Features:

* All the ``PropsSI`` overloads.  For all other types of inputs, the ``PropsSImulti`` function is now used

Issues Closed:

* `#375 <https://github.com/CoolProp/CoolProp/issues/375>`_ : If one input and one output to PropsSI, bubble error cleanly
* `#373 <https://github.com/CoolProp/CoolProp/issues/373>`_ : Move predefined mixture parsing to HelmholtzEOS initializer function
* `#372 <https://github.com/CoolProp/CoolProp/issues/372>`_ : Prandtl number is missing from AbstractState
* `#371 <https://github.com/CoolProp/CoolProp/issues/371>`_ : Parse inputs to PropsSI/PropsSI(vectorized) and turn into a vector of inputs
* `#370 <https://github.com/CoolProp/CoolProp/issues/370>`_ : Docs are missing all the fluid files
* `#368 <https://github.com/CoolProp/CoolProp/issues/368>`_ : CoolProp on iOS
* `#367 <https://github.com/CoolProp/CoolProp/issues/367>`_ : Python module architecture
* `#366 <https://github.com/CoolProp/CoolProp/issues/366>`_ : Get value of universal gas constant
* `#365 <https://github.com/CoolProp/CoolProp/issues/365>`_ : REFPROP_lib.h is missed in 5.0.4 source code zip
* `#364 <https://github.com/CoolProp/CoolProp/issues/364>`_ : Liquid and vapor saturation pressures are not the same for some fluids
* `#363 <https://github.com/CoolProp/CoolProp/issues/363>`_ : Revision synchronisation
* `#359 <https://github.com/CoolProp/CoolProp/issues/359>`_ : Add high-level function that allows for multiple outputs
* `#357 <https://github.com/CoolProp/CoolProp/issues/357>`_ : Vector functions and state class
* `#349 <https://github.com/CoolProp/CoolProp/issues/349>`_ : Host v4 docs

5.0.4
-----

BUGFIX: Lots of bugs squashed. 

New features: 

* Julia wrapper added
* Derivatives along the saturation line for pure fluids implemented
* Exposed the configuration getter/setter through SWIG (except for MATLAB)
* Added transport properties for xylenes and Ethylbenzene
* Surface tension for HFC pseudo-pures added

Issues Closed:

* `#355 <https://github.com/CoolProp/CoolProp/issues/355>`_ : In MSVC, too many symbols are exported in SWIG+MATLAB
* `#354 <https://github.com/CoolProp/CoolProp/issues/354>`_ : REFPROP headers
* `#353 <https://github.com/CoolProp/CoolProp/issues/353>`_ : Using HAPropsSI within circular reference on Mac Excel 2011 causes div/0 error!
* `#350 <https://github.com/CoolProp/CoolProp/issues/350>`_ : Python module docs
* `#347 <https://github.com/CoolProp/CoolProp/issues/347>`_ : Implement calc_melting_line for incompressibles
* `#346 <https://github.com/CoolProp/CoolProp/issues/346>`_ : Memory sanitizer is reporting errors with RPVersion function call
* `#344 <https://github.com/CoolProp/CoolProp/issues/344>`_ : skip typeerror in Excel to make 32-bit xlam work in 64-bit excel
* `#342 <https://github.com/CoolProp/CoolProp/issues/342>`_ : Refprop mixture with 4 components error
* `#339 <https://github.com/CoolProp/CoolProp/issues/339>`_ : Some SWIG tests fail due to the inclusion of rapidjson header
* `#337 <https://github.com/CoolProp/CoolProp/issues/337>`_ : ECS not yielding the proper values for eta and lambda
* `#332 <https://github.com/CoolProp/CoolProp/issues/332>`_ : Make the REFPROP wrapper code 1% more sane
* `#331 <https://github.com/CoolProp/CoolProp/issues/331>`_ : Excel wapper shouts errors (in Excel 2013)
* `#330 <https://github.com/CoolProp/CoolProp/issues/330>`_ : Implement ECS model for viscosity of xylenes and ethylbenzene
* `#326 <https://github.com/CoolProp/CoolProp/issues/326>`_ : expose configuration through SWIG
* `#325 <https://github.com/CoolProp/CoolProp/issues/325>`_ : Implement the generalized derivatives for REFPROP as well
* `#324 <https://github.com/CoolProp/CoolProp/issues/324>`_ : SetPath for Refprop
* `#322 <https://github.com/CoolProp/CoolProp/issues/322>`_ : Add method to AbstractState to return mixture component names
* `#321 <https://github.com/CoolProp/CoolProp/issues/321>`_ : Add more R-number aliases
* `#320 <https://github.com/CoolProp/CoolProp/issues/320>`_ : HAPropsSI("T", "V", 0.83, "R", 1, "P", 101325) & lead to infinite value
* `#319 <https://github.com/CoolProp/CoolProp/issues/319>`_ : Error in entropy calculation with TH inputs
* `#314 <https://github.com/CoolProp/CoolProp/issues/314>`_ : Add surface tension reference information to docs
* `#312 <https://github.com/CoolProp/CoolProp/issues/312>`_ : Small examples of the use of derivatives should be in docs
* `#309 <https://github.com/CoolProp/CoolProp/issues/309>`_ : MEG properties
* `#308 <https://github.com/CoolProp/CoolProp/issues/308>`_ : Set maximum states for saturation curves for pseudo-pures properly
* `#306 <https://github.com/CoolProp/CoolProp/issues/306>`_ : Surface Tension for HFC Pseudo-Pure is missing
* `#304 <https://github.com/CoolProp/CoolProp/issues/304>`_ : Develop some docs about hooking up with Julia code
* `#294 <https://github.com/CoolProp/CoolProp/issues/294>`_ : Add the clang sanitize tests to buildbot
* `#247 <https://github.com/CoolProp/CoolProp/issues/247>`_ : Implement thermal conductivity for o-Xylene, m-Xylene, p-Xylene, and Ethylbenzene
* `#238 <https://github.com/CoolProp/CoolProp/issues/238>`_ : add a function to retrieve derivatives along the saturation curve


5.0.3
-----
Bugfix, with some new functionality

The most important fix is for users of Microsoft Excel on windows. It is imperative to download a new CoolProp.dll, there was a serious bug in how Excel and CoolProp interact that has been fixed.

Issues Closed:

* `#293 <https://github.com/CoolProp/CoolProp/issues/293>`_ : Requirement for zipped source code file
* `#292 <https://github.com/CoolProp/CoolProp/issues/292>`_ : Update CycloHexane EOS
* `#289 <https://github.com/CoolProp/CoolProp/issues/289>`_ : Two-phase states don't work for DY flash
* `#288 <https://github.com/CoolProp/CoolProp/issues/288>`_ : Some calls in Excel throw FPU exceptions which throw error messages
* `#287 <https://github.com/CoolProp/CoolProp/issues/287>`_ : Predefined mixtures cannot be used in PropsSI
* `#285 <https://github.com/CoolProp/CoolProp/issues/285>`_ : Cannot solve for conductivity and viscosity
* `#284 <https://github.com/CoolProp/CoolProp/issues/284>`_ : Create build steps on the master that allow us to automate the releasing even more
* `#283 <https://github.com/CoolProp/CoolProp/issues/283>`_ : Change fullclean logic to use git pull to wipe the folder completely
* `#282 <https://github.com/CoolProp/CoolProp/issues/282>`_ : SWIG wrappers not converting errors in PropsSI to exceptions
* `#280 <https://github.com/CoolProp/CoolProp/issues/280>`_ : Describe the predefined mixtures with examples on website

5.0.2
-----
Bugfix.

Issues Closed:

* `#281 <https://github.com/CoolProp/CoolProp/issues/281>`_ : Surface Tension Errors
* `#278 <https://github.com/CoolProp/CoolProp/issues/278>`_ : Add script to generate milestone text automatically
* `#277 <https://github.com/CoolProp/CoolProp/issues/277>`_ : Fix doxygen docs for generalized residual helmholtz term
* `#275 <https://github.com/CoolProp/CoolProp/issues/275>`_ : Logscale densities for consistency plots
* `#274 <https://github.com/CoolProp/CoolProp/issues/274>`_ : P and D as inputs produces some errors
* `#273 <https://github.com/CoolProp/CoolProp/issues/273>`_ : hmolar, smolar etc. are incorrect for HEOS backend with PD inputs
* `#272 <https://github.com/CoolProp/CoolProp/issues/272>`_ : 32bit Pre-compiled Binary for C#
* `#254 <https://github.com/CoolProp/CoolProp/issues/254>`_ : Error : hapropsSI("R";"T";253;"B";252;"P";101325) lead to an error

5.0.1
-----
The first release with the automated release script. No major code changes.

5.0.0
-----
**MAJOR** The new version of CoolProp implementing the new structure based on AbstractState
**MAJOR** Some features have been temporarily (or permanently) deprecated
**MAJOR** CoolProp now supports mixtures
**MAJOR** Buildbot system powered by CMake set up to run builds after every commit