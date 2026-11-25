# Agent notes

Durable, hard-won knowledge about working on CoolProp: gotchas, design
decisions and verification techniques that are not obvious from the code.
Read the section for the area you are about to touch.  Add to this file (in
the PR that taught you the lesson) rather than keeping private notes; open
work belongs in Linear, not here.

Migrated from the Beads memory store on 2026-09-26; entries that only
recorded history of finished or abandoned work were dropped (they remain in
`.beads/issues.jsonl`).

## Build, CI and tooling

- **clang-format CI diffs against the live PR base.** `dev/ci/clang-format.sh`
  formats every file differing between `pull_request.base.sha` and HEAD, so a
  branch cut from a stale master inherits failures from misformatted files
  that landed on master since.  Fix by rebasing onto the current master tip.
  Check locally: `CLANG_FORMAT='uvx clang-format@18.1.8' bash dev/ci/clang-format.sh HEAD origin/master`.
- **clang-tidy in CI is line-anchored; preflight is whole-file.** CI runs
  `clang-tidy-diff.py` on `git diff -U0`, so touching one line of a legacy
  file does not surface its backlog.  `preflight.sh` passes whole files to
  `run-clang-tidy-staged.sh`, so it does surface the untouched backlog of any
  file you edit — despite its step being labelled "diff-only" (Linear
  COO-14).  `.clang-tidy` is in whitelist mode (PR #2802) because the broad
  set was ~60% style noise; clang-tidy CI is informational by design.
- **clang-tidy in Docker:** delete `dev/.fluiddepcache` and
  `dev/.incompdepcache` first; the host's pickle is incompatible with the
  container's Python.
- **clang-tidy `-fix` sweeps can corrupt code.** Anchor `--header-filter` to
  the repo root (`^<repo>/(include|src)/`) or it rewrites vendored headers
  under `build_catch/_deps`.  `run-clang-tidy`'s export+apply double-applies
  shared-header fixes (`: a(a) : a(a) {`), even at `-j1`; run `clang-tidy -fix`
  one TU at a time instead.  `cppcoreguidelines-prefer-member-initializer`
  can hoist a read of a member that is only filled later in the constructor
  body, silently zeroing results; always review a member-initializer sweep
  adversarially.  `modernize-use-auto` / `modernize-use-using` are low-risk
  and welcome as standalone sweeps.
- **Catch2 v3 in this repo:** use `Catch::Approx`, and `~[slow]` (not
  `[!slow]`) to exclude slow tests.
- **cppcheck `syntaxError` means the rest of the file went unchecked.**
  cppcheck runs without Catch2's headers, so `dev/ci/cppcheck/catch2.cfg`
  defines `TEST_CASE` / `TEST_CASE_METHOD`; without it, either macro inside a
  namespace drops the whole file (COO-137).  cppcheck 2.x also cannot parse
  unary-plus lambdas in a braced list (`for (auto f : {+[](...){...}})`);
  use a named function-pointer array.  cppcheck reports one syntaxError per
  file and the passes find different ones, so re-run after each fix.  To
  locate one, cut the file at top-level boundaries with namespaces still
  balanced, or the cut itself is reported.  Never suppress `syntaxError`.
- **CMake per-target properties must follow every target-creation path.**
  `Main` is created under both `COOLPROP_MY_MAIN` (CodeQL/Coverity/IWYU CI)
  and `COOLPROP_MAIN_MODULE`.  A `set_property(TARGET ...)` nested in one of
  them silently no-ops on the other; put instrumentation blocks at the bottom
  of `CMakeLists.txt` with `if(TARGET ...)`.  Symptom: "Manually-specified
  variables were not used by the project".
- **Moving headers: forwarding shims don't cover tooling that reads headers
  by path.** `wrappers/Python/generate_constants_module.py` scrapes enums from
  `include/*.h` by file I/O; after a move it reads a shim, matches nothing and
  the wheel build dies.  Grep `wrappers/ dev/ setup.py CMakeLists` for
  `open()` / `os.path.join('include', ...)` and `DEPENDS` on header paths, and
  repoint them in the same PR.
- **CodeQL vendored-code noise:** `paths-ignore` is ignored for C/C++ in a
  build-based analysis.  `build-mode: none` would honour it but was rejected
  (PR #3013 closed): the buildless extractor is less accurate on a
  macro/template-heavy codebase.  Dismiss vendored alerts in the Security UI
  instead.
- **Docs build installs only the Python wheel.** `docs_docker-run.yml`
  executes notebooks with `jupyter nbconvert --execute` after `pip wheel .`;
  it does not build the C++ test runner.  Figures in doc notebooks must come
  from the Python API alone, never from Catch2 test output.
- **nanobind stub drift gate is reproducible locally**, but only
  with the pinned pair, not the ambient interpreter:

  ```bash
  uv venv --python 3.12 /tmp/stubvenv
  VIRTUAL_ENV=/tmp/stubvenv uv pip install 'nanobind==2.12.0' cython scikit-build-core ninja
  COOLPROP_NANOBIND=ON VIRTUAL_ENV=/tmp/stubvenv uv build --wheel --no-build-isolation \
    --config-setting=cmake.define.CMAKE_BUILD_TYPE=Release \
    --config-setting=cmake.define.COOLPROP_NANOBIND_REGEN_STUB=ON --out-dir /tmp/stubdist .
  unzip -o /tmp/stubdist/*.whl 'CoolProp/CoolProp/*.pyi' -d /tmp/stubcheck
  diff -u wrappers/Python/_nanobind/CoolProp.pyi /tmp/stubcheck/CoolProp/CoolProp/__init__.pyi
  diff -u wrappers/Python/_nanobind/Plot.pyi /tmp/stubcheck/CoolProp/CoolProp/Plot.pyi
  ```

  Take the generated stubs verbatim; never hand-edit them.  A default value
  that stubgen writes with a name the stub does not define (an enum of the
  core in `Plot.pyi`) is fixed in the binding with `nb::arg(...).sig("...")`.
- **TestPyPI pruning is manual** (since 2025-11-14, PyPI's new-device email
  verification blocks unattended TOTP logins).  Stay logged into TestPyPI in
  a browser, then `python3 dev/testpypi_delete.py --keep 10` (dry run) and
  `--do-it`.  See PR #3202.
- **RapidJSON was removed** (#3112); use nlohmann/json (fetched via
  `cmake/dependencies.cmake`) and do not reintroduce it.

- **GitHub's ubuntu-26.04 image builds with FMA by default.**  Its GCC 15.2
  reports `-march=x86-64-v3` / `-mfma [enabled]` with no flags passed, while
  the `gcc` in a stock `ubuntu:26.04` Docker container (same `15.2.0-16ubuntu1`
  version string) defaults to plain `x86-64`.  With GCC's default
  `-ffp-contract=fast` (GNU dialects, which CoolProp builds with), any test that
  compares two separately compiled code paths bit-for-bit can fail on CI and
  pass in the container (COO-182: `all_deltaonly()` vs `all()`).  To
  reproduce, add `-march=x86-64-v3` in the container; to diagnose, check
  `gcc -Q --help=target | grep -E 'march|mfma'` and count `vfmadd` in the
  object.  Compare such paths to a scaled tolerance, not bit-for-bit.
- **Catch2 test names must not contain commas.**  A comma in a test spec
  is OR, so "Cubic T,P root choice" cannot be selected by name.  Write
  "TP", "P-input", etc.
- **Rewriting a stack of PR branches:** commit the fix on the owning
  branch, then from the top branch `git rebase --update-refs --onto
  <owner-new> <owner-old>`; it moves every intermediate branch ref in one
  pass.  Use the old SHA the upper branch actually sits on (after an
  amend it differs from the branch's own old tip).  Enable `rerere`.

## REFPROP

- **Running `[refprop]` Catch2 tests locally:** without configuration they
  SKIP silently (and `preflight.sh` reports green).  Point CoolProp at a
  REFPROP install via the `COOLPROP_REFPROP_ROOT` environment variable, or
  via the config keys `ALTERNATIVE_REFPROP_PATH` and
  `ALTERNATIVE_REFPROP_LIBRARY_PATH` set as `COOLPROP_ALTERNATIVE_REFPROP_PATH`
  etc. (config keys are read from `COOLPROP_<KEY>` environment variables).
- **Fork PRs never run REFPROP tests.** `test_catch2.yml` builds REFPROP only
  when the PR head is in CoolProp/CoolProp (it needs
  `secrets.REFPROP_GPG_PASSPHRASE`).  A green fork PR is not evidence that
  REFPROP tests pass.
- **The REFPROP backend is not thread-safe**, and never was: REFPROP itself
  is not reentrant.  Callers must serialize REFPROP calls.
- **Composition is a per-call argument** to every REFPROP routine; the only
  global state is which fluids are loaded.  Saturated-phase states can be
  modelled as lightweight backend instances that reuse the host's loaded
  component string and carry their own composition — use
  `mole_fractions_liq` / `mole_fractions_vap`, not the bulk composition.
- **REFPROP 10.1 beta** ships Windows DLLs only, and its FLD files carry
  materially newer transport models than 10.0 (argon, nitrogen, methane,
  ethanol, R134a, R32, D2O viscosity all differ).  Audit fluid content
  against the 10.1 FLD set, not 10.0.  FLD files also carry corrected
  coefficients without citing the erratum.

## Thermodynamics and flash routines

- **HS single-phase non-uniqueness is spurious.** A second (h,s)-matching
  root can be mechanically and adiabatically stable yet thermally unstable
  (cv < 0) — seen on dense supercritical hydrogen.  Accept a root only if it
  reproduces (h,s), lies in [Tmin, Tmax], and has dp/drho|_T > 0 AND cv > 0.
- **Generate flash test points in (p,T), not density bands.** A (T,rho) grid
  lands states inside the spinodal (cv < 0), producing false failures and
  meaningless comparisons; (p,T) always lands on the stable root.
- **The backend's `get_superanc()` throws for mixtures and pseudo-pure
  fluids** (it returns nullptr only for a pure fluid that has none); guard it
  with `is_pure()`.  H,S flashes without a superancillary (pseudo-pure, or
  pure without one) use the dome-free legs of `hs_cascade` (PR #3182).
- **Cubic backends have no classic ancillaries.** `components[0].ancillaries`
  (pL, pV, hL, ...) is empty on SRK/PR, so any path that reaches it throws
  "type not set" (#3441).  `get_superanc()` is null there too.  Reach the
  cubic Chebyshev superancillary through the virtuals
  `superanc_saturation_at_p` / `superanc_saturation_at_T` (or
  `calc_saturation_ancillary`); inside the cubic backend everything goes
  through the private `cubic_superanc_eval`.  Q,T, P,Q and the T,P root
  choice use it too (COO-136); the old secant saturation solve is no
  fallback (~1e-7 in p near Tc, and silently the trivial or wrong root
  outside the table).  When adding a superancillary, check that every input
  pair needing saturation reaches it: the cubic one (#2739) sat unused for
  Q,T and P,Q because `AbstractCubicBackend::update` intercepts both.
  Direct `get_superanc()` callers that remain are listed in COO-146.
- **Cubic mixtures: never solve a density with neither a phase nor a guess.**
  With three roots the cubic `solver_rho_Tp` then calls `p_critical()`, a full
  critical-point search (about 7 ms for 2 components, steeply more with N).  #3451's
  recovery pass hit it on every density solve (#3597); a guess is now followed.
  Tests should count `calc_all_critical_points` calls
  (`CoolProp-Tests-CubicMixtureGuess.cpp`), not time them.
- **Cubic `solver_rho_Tp_global` obeys an imposed phase; the HEOS one does
  not.**  `SatL` always has liquid imposed, so the stability test evaluates a
  gas feed with two roots at its metastable liquid root, calls it unstable, and
  only the split solver's fallback recovers the gas (#3597; still open).  The
  obvious fixes lose real splits.  Of 8800 PQ-to-PT dome states master misses
  273; solving the feed unimposed misses 1633 (the liquid-like trial then stays
  on the gas branch), and ignoring the imposed phase altogether misses 284
  (CO2/ethane and R410A: `SatV`'s imposed gas phase stands in for the
  opposite-branch trials, which need Chebyshev roots and skip cubics).  Feed unimposed plus liquid imposed on the liquid-like
  Wilson trial misses 111.  Judge such a change by a dome sweep, not timings.
- **Cubic gaps a Helmholtz-fluid path assumes away.**  No triple point
  (`triple_liquid.rhomolar` is +inf: seed with `calc_rhomolar_max_bound()`);
  `rhomolar_critical()` must be the EOS's own (Z_c R Tc/pc), not a
  real-fluid fit, because it decides gas vs liquid; the Soave alpha passes
  through zero at T* = Tc(1+1/m)^2, where T,U is degenerate, Ttilde turns
  over and the closed-form cubic loses ~10 digits (COO-145).  Within ~1e-9
  of the critical point the expanded cubic's coefficients cancel: polish
  roots against p(rho) in EOS form, never the polynomial.
- **Flash non-uniqueness is physical, not a bug.**  H,T is non-unique below
  the max Joule-Thomson inversion temperature (~5 Tc): a cold compressed
  liquid's (T, h) is also a two-phase state.  Pure cubics refuse H,T; HEOS
  and cubic mixtures still return the two-phase root (COO-158).  T,S is
  unique for pure cubics but not for cold liquid water (alpha_p < 0 below
  the density maximum): HEOS and REFPROP both return the two-phase state.
  A round-trip test across such a region must not demand the reference
  state back; exempt it only with an independent check that the input is
  inside the saturated [L, V] interval at that T.
- **Compute long constants, never transcribe them.**  #2745's "exact"
  SRK/PR Omegas had doubled and dropped digits (fixed in #3515), which put
  SRK's own critical point ~1e-9 below the input Tc/pc and surfaced as a
  phantom "end gap" in the superancillary.  Evaluate the closed form in
  high precision in a committed script and store the correctly rounded
  double as a hex literal (`dev/scripts/cubic_constants.py` writes
  `src/Backends/Cubics/CubicConstants.h`; evaluating in double at run time
  is 1-2 ulp off and platform-dependent; the C++ build has no
  multiprecision).  Test the defining property (critical conditions), not
  digits.
- **`recalculate_singlephase_phase()` reads `_rhomolar`.**  Assign the new
  density before calling it (cubic `solver_rho_Tp` did not; fixed with #3441).
- **HS cascade rules:** every cascade leg must be exception-guarded (a leg
  leaving the EOS domain otherwise aborts the cascade), and phase
  determination for two-phase (h,s) inputs is the dispatcher's job, done
  before the single-phase cascade — the cascade only finds metastable
  in-dome roots there.
- **`solver_rho_Tp` is fragile on its own.** Called standalone with no
  guess it fails on 25–65% of mixture (T,p) cells where the full PT flash
  succeeds, clustered near phase boundaries; the flash succeeds only because
  `update_TP` falls through to the bracketed `solver_rho_Tp_global`.  Don't
  call it directly without a bracket/fallback.
- **P+{H,S,U} is a 1-D problem.** A 2-D homotopy/continuation flash (good for
  HS) regressed P+X in the bulk and was reverted; speed P+X by improving the
  1-D solver (TOMS748, warm starts).  Validate flash speed on a whole
  consistency-plot grid, never on hand-picked points — a cherry-picked bench
  hid that regression.
- **HS flash scope:** the shipped Helmholtz HS happy path (PR #2999) is the
  intended end state for the EOS backend.  Fast 2-D (h,s) inversion and
  call-to-call caching belong to SVDSBTL, not the Helmholtz path.
- **Domain clamps need slack sized between two scales.** Water's (h,s)
  homotopy folds ~0.03% below Tmin near the density anomaly, while a
  spurious ortho-hydrogen basin sits ~30% below; a 2% relative slack admits
  the first and blocks the second.
- **XN_DEPENDENT const-(T,rho) composition derivatives are derivable from the
  XN_INDEPENDENT (mole-number) branch.** `ndln_fugacity_i_dnj__constT_V_xi` is
  valid for both flags; build composition Jacobians from it plus an explicit
  rho-coupling term rather than via simplex projection.
- **Thread-safety contract is one AbstractState per thread.** TSan reports on
  `IdealHelmholtzContainer` caches (#2844) are false positives: the ideal
  cache path is dead (`cache_values` defaults to false everywhere) and the
  residual cache is per-instance.  Sharing one backend across threads would
  genuinely race.
- **Windows is where NaN/_HUGE FP flags hurt** (#3012): Delphi unmasks FP
  exceptions and Excel/VBA poll the status word.  `CoolProp::fpu_guard`
  (`include/CoolProp/FPUGuard.h`) keys on `_WIN32` (any compiler, incl.
  MinGW), not `_MSC_VER`.
- **Mixture critical points: filter on `stable`, nothing else** (#3528).
  `calc_all_critical_points()` also returns spurious points; most sit at
  liquid-like density and negative pressure, which `calc_critical_point`
  already marks unstable (negative pressure is physical, so don't filter on
  p).  The TPD verdict on the ~100 K / 100s-of-MPa points is *marginal*: it
  flips between builds (R470A), so don't assert per-blend stable counts there
  (COO-155).  The N×N L*/M* formulation is general-N, but Air (ternary) yields
  four stable points (COO-154).
- **`MixtureDerivatives` scalar chains nest O(N) sums** (#3560).  The
  third-order composition derivatives (`n2Aijk` -> ... ->
  `d2_ndalphardni_dxj_dxk`) each re-run a `sum_m x_m ...` that doesn't
  depend on the outer index, so filling an N×N×N tensor from the scalars is
  O(N^7): minutes for a 10-component natural gas.  Build whole arrays with
  `nAij_matrix` / `n2Aijk_tensor` (each sub-term once); keep the scalars as
  the test reference.  What remains is per-leaf: the GERG reducing
  derivatives recompute rho_r(x), T_r(x) on every call.  Profile with a
  RelWithDebInfo shared lib and `sample <pid>`: the pip wheel is stripped.
- **Stability-test shortcuts must be swept on near-azeotropes** (COO-193).
  #3427 took the SS proximity heuristic ("trial converging to the trivial
  solution") as proof of stability and skipped the TPD minimizer.  For
  R410A or R507A the real incipient phase is also close to z with a shallow
  tm, so the PT flash called most dome states single-phase.  Gas-mixture
  benchmarks (Amarillo, air) never showed it.  Check any change to the
  stability test with a QT-to-PT dome sweep over refrigerant blends
  (`CoolProp-Tests-NearAzeotrope.cpp`), not only with timings.  The same
  test decides critical-point stability: with the minimizer back, two of
  R452A's three "stable" critical points (COO-154) became unstable.
- **Stability trials can sit on the feed's density branch** (COO-193).
  Trials take the lowest-Gibbs root of their composition.  A Wilson trial
  that starts close to the feed (all K near or below 1) lands on the feed's
  own branch, and SS collapses to the trivial solution while the incipient
  phase, on the other branch, has tpd < 0 (CO2/ethane 20/80 at 207 K:
  -0.036, reported "stable").  So "close to z in composition" means trivial
  only on the same branch: the proximity shortcut is gated on the trial's
  root being within 1.5x of the feed's.  A trial stuck on the feed's branch
  is re-run on the other one (`opposite_branch_trial`).  Its roots must come
  from the Chebyshev solver (`solver_rho_Tp_cheb_branch`), never a
  phase-imposed Newton solve: that finds alpha^r-well roots that pass the
  pressure and mechanical-stability checks (R454B, 20 MPa: tm = -1e12, then
  a pure-R32/pure-R1234yf "split").  Measure with a sweep over many
  mixtures, not the six blends a fix was tuned on: a threshold rule fitted
  on six blends recovered only half the misses over 203.

## Fluid data, superancillaries and transport

- **Validate E/Z SMILES against the stored InChI.** Slash directions depend
  on traversal and branch order; matching slash characters alone do not
  identify E versus Z. An `@` on an alkene carbon does not encode double-bond
  stereochemistry. For metadata corrections, parse the SMILES with a chemistry
  toolkit and compare its InChI and InChIKey with `INFO.INCHI_STRING` and
  `INFO.INCHI_KEY`.
- **Superancillary provenance.** fastchebpure fits superancillaries from the
  CoolProp EOS via its submodule pin, and the docs deviation plots download a
  pinned fastchebpure release.  Changing a fluid's EOS in `dev/fluids` trips the
  `superanc-pin` gate (`dev/scripts/check_superanc_release_pin.py`) until
  fastchebpure regenerates, is re-tagged, and the docs pin is bumped — or the
  fluid is added to `PENDING_UPSTREAM`.  Two independent freshness axes:
  EOS vs embedded superancillary (`source_eos_hash` + check points) and EOS vs
  docs reference.  Suspect provenance lag before the superancillary math.
- **Don't differentiate a fitted Chebyshev expansion** when the derivative
  has a closed form: sample the derivative directly at the nodes.
  Differentiation truncates the top coefficient (a 180x accuracy loss was
  measured at loose build tolerance).
- **Regenerating a superancillary:** never hand-edit `source_eos_hash`.  Use
  fastchebpure's `fitcheb` CLI: `fitcheb -f FLUID -d <coolprop checkout>` to
  fit and check, then `fitcheb inject -f FLUID -d <coolprop checkout>` to
  replace the SUPERANCILLARY block, add check points and stamp the hash.
- **Incompressible caloric fits come in two kinds, both loaded.**  A fluid
  may ship `density`/`specific_heat` (classic polynomial) and
  `density_cheb`/`specific_heat_cheb`; `IncompressibleFluid` keeps both and
  every caloric function takes a `chebyshev` flag (`useChebyshevFit` picks, with
  fallback to whichever exists).  `IncompressibleBackend` snapshots
  `INCOMP_PREFER_CHEBYSHEV_CALORIC` in its constructor and must pass it to
  every caloric call, including residuals inside `T_rho`/`T_c`: a call that
  falls back to the default argument reads the live config and can mix fits
  with the cached `_hmass_ref`/`_smass_ref`.  Derived matrices (`ddT`,
  `int_dT`, `int_dTdivT`, `int_dTdivT_log`) exist only after `validate()`;
  evaluation throws on an empty one rather than returning 0.  `setTbase()`
  drops the polynomial entropy matrices (they embed the value of Tbase), so
  call `validate()` again after it.  A Chebyshev fit may sit in a classic slot
  only when the fluid has no separate `_cheb` entry for that property.
- **Incompressible fits are evaluated in three places that must agree to the
  last bit.**  The full evaluation (`polyEvalMatrix`, `chebEvalMatrix`) and
  the per-state `CompositionFits` (`FitAtComposition::eval`, built by
  `fixComposition` when the composition changes) share `rowAtComposition`,
  `hornerSum` and `clenshawSum`.  Change the arithmetic only there.  To
  check, dump rho, cp, h, s, mu, k as `%a` hex floats for every fluid before
  and after and `diff` them (26485 values per fit kind, zero differences
  expected on x86-64 with the default flags).  A compiler that contracts
  multiply-adds (aarch64, `-ffp-contract=fast`) can round the two paths
  differently by a few ulp; there use the backend test's 8-epsilon check.  A backend rebuilds its `CompositionFits` when
  `IncompressibleFluid::data_version` changes, which `validate()` and every
  setter of cached data (fits, `setTbase`, `setxbase`) bump through
  `nextDataVersion()`; a new setter for cached data must bump it too.
  `add_fluids_as_JSON` replaces a fluid in place under existing backends.
- **Incompressible inversions use `BoundedNewton` (src/Solvers.cpp).**
  `T_rho`, `T_c` and the `HmassP`/`PSmass` flashes search the liquid range
  through `solveInLiquidRange`: `[TminAt(p, x), Tmax]` first, where
  `TminAt` is the freezing point when it is above `Tmin` (density fits
  extrapolated below freezing turn over, and on the full range a valid
  density was reached twice); and `[TminAt, TmaxAt]` only if that lands on a
  boiling state or finds no root (fits extrapolated above the boiling point
  turn over too, COO-190).  `TmaxAt` solves psat(T) = p with Brent, about
  300 ns: computing it on every flash doubled the cost when p changes with
  every call, hence the lazy second pass.  `BoundedNewton` needs a
  sign change over the range and throws otherwise, like Brent; it must never
  orient its bracket from the sign of f' alone (that silently returned the
  other root, 77 K off).  The flash slopes are exact (`c - p T d2v/dT2` for h, `c/T - p d2v/dT2`
  for s): cp
  alone needs 8 iterations at 500 bar instead of 4.  Round-trip every shipped
  fluid (T -> D, h, s -> T) after touching any of this; `INCOMP::Air` is a
  1-atm gas fit and fails at 10 bar in every version.
- **Where the liquid ends is one rule, shared by three places.**
  `isLiquidAt(T, p, x)` decides liquid vs boiling for both `checkP` and the
  first pass of `solveInLiquidRange`; it accepts a state up to
  `saturationMargin()` (2e-9 Tmax, about 1 uK) above the boiling point, because
  a flash from a saturated state lands a few ulp above T_sat.  `TmaxAt` returns
  0.3 to 0.7 margin above where the liquid ends (T_sat, or TminPsat when the
  liquid boils just above it), and the second pass searches 0.75 margin
  beyond that and checks its result with `isLiquidAt`.  Change one of these
  and the others drift apart: every step of PR #3567 that moved only one of
  them broke saturated round trips or left a root unreachable.  Two traps:
  CoolProp's `Brent` stops once half its bracket is below `t`, so its result
  is good to 2 t, not t; and it also returns an end of its interval when
  |f| there is below `t`, so scale a residual (here to ulps of p) if `t` is
  meant for x only.  Test with a slowly rising psat (0.001 T Pa) and a
  density fit that turns over, and with the 8-fluid saturation sweep
  (QT_INPUTS Q = 0, then H,P and P,S back, about 2000 temperatures each).
- **Undoing the reference shift rounds.**  A state reports
  h = h_ref + raw - hmass_ref, and the flash rebuilds raw_target = h - h_ref
  + hmass_ref; these cancel only to about 1 ulp.  At an end of the search
  range that ulp removes the sign change, or (on a fit that is not monotonic,
  TX22 and DowJ2 entropy at 100 MPa) keeps one and leads to another root.
  `CaloricResidual` returns exactly 0 where the reported value (the same
  expression as `calc_hmass`/`calc_smass`) equals the input bit for bit.  Do
  not use a tolerance band instead: one sized from |hmass_ref| grows with an
  extreme reference state and snaps a bracketed root to an end, and recovering
  only after the solve fails misses the case that keeps a sign change.  Any new
  flash that rebuilds a target from a shifted value needs the same.
- **Incompressible cp is `dh/dT|p`, not the fitted `c` (#3577).**  h and s
  carry a pressure term from the density fit, taken from `p_fit_Pa` (1 atm,
  the pressure of the data), not from 0 Pa: from 0 Pa, `INCOMP::Air`'s fit
  artefacts in d2v/dT2 made cp negative near Tmax at 1 atm.  `cpmass`,
  `dhdTatPx`, `dsdTatPx` and the flash slope all go through
  `dhdT_from`/`dsdT_from` in `IncompressibleBackend.h`; `cmass()` is the raw
  fit and `cvmass` returns it, because v(T) alone has zero compressibility
  and leaves cv undefined.  `C` equals the fit only at 1 atm: a test that
  checks a fitted heat capacity at another pressure must read `O` (Cvmass).
- **The default incompressible reference temperature depends on the range
  (COO-203).**  `default_T_ref_K()` clamps 293.15 K into [Tmin, Tmax]: fits
  extrapolated outside the range gave IceNA a reference entropy of 2.4e11
  J/kg/K.  A default reference (`T_ref_is_default`) is re-resolved in
  `drop_stale_reference_values()` when the fluid's `data_version` changes, so
  a fluid replaced with another range moves it; an explicit
  `set_reference_state` keeps its T0.  `setTmin`/`setTmax` bump
  `data_version` for this reason; a new setter of anything the reference or
  the fits depend on must do the same.  The same check drops every property
  the backend caches (`clear_property_caches`), and every cached accessor
  calls it first.  Viscosity and conductivity are cached in `AbstractState`
  itself and are not covered: after a replacement they stay stale until
  `update()` (COO-204).
- **Profile incompressibles from C++, not Python.**  The Python wrapper adds
  about 0.12 us per call, as much as a whole property evaluation.  Build a
  static `RelWithDebInfo` library, time `AbstractState` in a small C++
  program, and use `valgrind --tool=callgrind` (there is no `perf` in the
  cloud container).  A thrown and caught exception costs about 1.3 us: one in
  `checkP` (psat below `TminPsat`) used to dominate every update of 19 fluids.
- **Config side effects belong in `on_config_changed`** (`src/Configuration.cpp`),
  which every setter calls (the JSON setter only for changed values).  A side
  effect wired into one setter is skipped when the key is set another way.
- **Changing one fluid's viscosity can move other fluids.** Two couplings:
  (1) ECS — `viscosity_ECS` evaluates the reference fluid's
  `calc_viscosity_background()` (R134a moves seven refrigerants by 1–2.6%);
  (2) conductivity — `conductivity_dilute` type `eta0_and_poly`,
  `conductivity_critical_simplified_Olchowy_Sengers` (divides by viscosity)
  and `conductivity_dilute_hardcoded_ethane`.  Before changing a viscosity
  model, grep `dev/fluids` for `reference_fluid` pointing at it and for those
  three conductivity types; measure by swapping only that fluid's JSON.
- **rhosr-CS viscosity:** `rhosr_critical` is not a fitted parameter — it is
  `rho_c*R*(tau*dalphar_dtau - alphar)` at the EOS critical point and must be
  recomputed whenever the EOS changes.  `C` and `rhosr_critical` are coupled
  and must be refit together (see `dev/scripts/fit_R1233zdE_viscosity.py`).
- **Validate a transport correlation against its paper's own verification
  points.** That is what catches typeset-equation defects (missing exp, sign
  in a denominator, wrong exponent — including in a numerator).  Do not
  prioritise by "signs hide in denominators".

## Working with papers

- **Never diagnose an equation from PDF text extraction.** Display equations
  are images; extraction loses fraction bars and superscripts.  Read the
  rendered page.
- **Recovering an equation from a .docx** (NIST `get_pdf.cfm` sometimes serves
  one): it is a ZIP; `word/document.xml` maps the paragraph with the equation
  number to `word/media/imageN.wmf`.  In the WMF, `ExtTextOut` (0x0A32)
  records carry glyph + x/y and `CreateFontIndirect` (0x02FB) the font
  height; sort by x and read y against the baseline to separate sub- and
  superscripts.
- **Before reporting an error in a published paper, check for a correction:**
  `curl -s https://api.crossref.org/works/<doi> | jq '.message."updated-by"'`.
- **Bollengier, Brown & Shaw 2019 water EOS** (J. Chem. Phys. 151:054501):
  a B-spline Gibbs surface for liquid water to 2300 MPa, covering the
  high-pressure domain IAPWS-95 in CoolProp refuses.  Its reference code
  (SeaFreeze) is GPL-3 and cannot be vendored into MIT CoolProp; the
  coefficients are in the paper's supplementary material, so a clean-room
  implementation is the path.  Its reference state differs from IAPWS
  (constant h0, s0 shift needed).

## Tabular backends and verification

- **Adversarial mid-cell probe for tabular accuracy claims:** evaluate at the
  centres of training-grid cells and compare with the at-node error.  Same
  order of magnitude means the interpolation is real; orders of magnitude
  worse means the benchmark was accidentally near-node.  (The original
  probe, used on SVDSBTL Water in May 2026, never reached master.)
- **BICUBIC inversion failures on (P,T) tables** (#1301): the fix is the
  SVDSBTL DT-indexed surface, not patching `invert_single_phase_y`.  The
  `find_nearest_neighbor` bisection improvement from PR #2892 is worth
  reviving on its own (in-table p99 error 33% → 8.8% on a CO2 grid).
- **Phase-envelope tracer baseline** (`PhaseEnvelopeRoutines::build`, Sep
  2026): fast but unreliable — natural-gas mixtures never set `built=true`
  (bubble curve hits T_min before p drops below the 100 Pa start), some
  mixtures follow the trivial solution past the critical point to GPa
  pressures, and heavy traces underflow to NaN.  Use these cases when
  judging a replacement tracer.

## Property plots (CoolProp::Plot)

- **`src/CoolPropPlot.cpp` is a port of CoolPlot's numerics** (the Python
  package: `thermo/isolines.py`, `tracer.py`, `stepping.py`, `fluid.py` of
  CoolPlot 0.7+).  Names and the order of steps follow the Python files on
  purpose; fix a numerics bug in both.  CoolPlot uses the port through
  `CoolProp.CoolProp.Plot` when it is present and falls back to its Python
  copy, which its tests treat as the reference.
- **The nanobind submodule lives at `CoolProp.CoolProp.Plot`**, not
  `CoolProp.Plot`: the package `__init__` is kept identical to the legacy
  one (import-tree parity gate), and `CoolProp.Plots` is the old Python
  plotting package.  Its stubs are a stub package,
  `CoolProp/CoolProp/__init__.pyi` and `CoolProp/CoolProp/Plot.pyi`: a stub
  module `CoolProp/CoolProp.pyi` cannot have submodules, and a `Plot.pyi`
  next to it types the non-existent `CoolProp.Plot`.  The regeneration
  imports the built module as `CoolProp.CoolProp` from a staging package;
  imported as a top-level `CoolProp`, stubgen names the types in `Plot.pyi`
  `CoolProp.<type>`, which do not resolve.
- **Releasing the GIL** in `PropertyPlot.calc_isolines`, `isoline_range`
  and `value_at` is safe for Python objects: the arguments are converted
  to C++ before the call (the settings are copied, see the binding) and
  the results after it.  `AbstractState.update` keeps the GIL: existing
  code that shares one state between Python threads relies on the GIL to
  serialise the calls.  `AbstractState.update_nogil` is the opt-in that
  releases it, for code that gives each thread its own state.  CoolProp's
  lazily built shared data (the fluid and mixture libraries, the caloric
  superancillaries) is built under `std::call_once` or a mutex, so distinct
  HEOS, cubic and GERG states are safe in parallel.  Not safe: a thread
  that changes the configuration (`set_config_*` writes a global without a
  lock), and a thread that builds tabular tables (creating a TTSE or
  BICUBIC state, or setting its fractions, when its tables are not on
  disk yet): the build toggles `DONT_CHECK_PROPERTY_LIMITS` globally
  (`TabularBackends.cpp`, COO-144).
  REFPROP states are serialised by `REFPROP_mutex()`.  Every other binding keeps the GIL; add
  `nb::call_guard<nb::gil_scoped_release>()` only to calls that never
  re-enter Python.  A REFPROP plot is safe without the GIL too: the
  backend serialises its calls through `REFPROP_mutex()`.  But another
  thread's REFPROP state can now run between two calls of a plot, and
  REFPROP keeps one set of loaded fluids: interaction parameters set
  with `set_binary_interaction_*`, a reference state or a custom `.MIX`
  composition are lost when that other state loads different fluids
  (see `REFPROPMixtureBackend.h`).
- **Every isoline tracer builds its own phase envelope**, about a quarter
  of the time of a traced line for a blend.  Do not share or cache that
  saturation state between tracers: the mixture saturation flashes take
  their guesses from it, so a line would depend on the lines computed
  before it, and lines are cached and computed one by one in any order
  (tried in COO-142 and reverted; CoolPlot's `tracer.py` has the same note).
- **No `==` or `!=` between doubles in the plot numerics**, and no
  tolerance in their place where an exact answer is meant.  The envelope's
  Q flags (exactly 0 or 1) and the markers of `generate_update_pair` are
  compared by side (`Q > 0.5`); an exact root or a vanished denominator is
  `is_zero` (a classification); sign tests use `sign_of` (integers); a
  cache key or a copied sample is `same_bits`; the start and end of a
  stepped range are passed as flags, not found by comparing `u`; tests
  that mean "the same double" use `WithinULP(x, 0)`.  `-Wfloat-equal` on
  `src/CoolPropPlot.cpp` is clean (the warnings left come from the bundled
  fmt).  The rewrite was checked bit-identical on about 1100 isolines of
  pure fluids, blends and an incompressible.  CoolPlot has the same
  helpers in `coolplot/thermo/floats.py`; its `is_zero` classifies the
  bits too, since `abs(v) < math.ulp(0.0)` is false for 0.0 when an
  extension built with fast-math sets denormals-are-zero.
- **A `PropertyPlot` serialises its own calls**: every line of a plot is
  computed with the plot's one `AbstractState`, updated in place, so
  `calc_isolines`, `isoline_range` and `value_at` hold a mutex (they
  release the GIL, and two Python threads really do call them at once).
  Parallel work needs one plot per thread (CoolPlot's adapter keeps one
  per thread, after a review of its thread workers found a shared plot).
- **The mixture H-p flash on master refuses two-phase states 8.0.0 solved**
  (R32/R125 70/30 at h = 4e5 J/kg, p = 9.65 bar: the residual gate after the
  temperature search throws "the (T,p) flash is misclassifying the phase").
  Tests of the tracer against the flash compare only where the flash
  converges; the isoline tracer itself never needs that flash.
- **Building the Python wheel offline:** the catch build fills
  `.cpm_cache/` in the repo; `pip wheel .` configures a fresh tree, so pass
  `CMAKE_ARGS="-DCPM_SOURCE_CACHE=$PWD/.cpm_cache"` (plus the valijson
  source if downloads are blocked) or Eigen is fetched again (gitlab
  returns 429 quickly).
