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
  file does not surface its backlog.  `preflight.sh` runs clang-tidy on the
  whole file and will (tracked in Linear, "preflight.sh and CI gates
  overhaul").  `.clang-tidy` is in whitelist mode (PR #2802) because the
  broad set was ~60% style noise; clang-tidy CI is informational by design.
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
- **nanobind `CoolProp.pyi` drift gate is reproducible locally**, but only
  with the pinned pair, not the ambient interpreter:

  ```bash
  uv venv --python 3.12 /tmp/stubvenv
  VIRTUAL_ENV=/tmp/stubvenv uv pip install 'nanobind==2.12.0' cython scikit-build-core ninja
  COOLPROP_NANOBIND=ON VIRTUAL_ENV=/tmp/stubvenv uv build --wheel --no-build-isolation \
    --config-setting=cmake.define.CMAKE_BUILD_TYPE=Release \
    --config-setting=cmake.define.COOLPROP_NANOBIND_REGEN_STUB=ON --out-dir /tmp/stubdist .
  unzip -o /tmp/stubdist/*.whl 'CoolProp/CoolProp.pyi' -d /tmp/stubcheck
  diff -u wrappers/Python/_nanobind/CoolProp.pyi /tmp/stubcheck/CoolProp/CoolProp.pyi
  ```

  Take the generated stub verbatim; never hand-edit it.
- **TestPyPI pruning is manual** (since 2025-11-14, PyPI's new-device email
  verification blocks unattended TOTP logins).  Stay logged into TestPyPI in
  a browser, then `python3 dev/testpypi_delete.py --keep 10` (dry run) and
  `--do-it`.  See PR #3202.
- **rapidjson is deprecated** and being removed; use nlohmann/json
  (fetched via `cmake/dependencies.cmake`) for new JSON code.

## REFPROP

- **Running `[refprop]` Catch2 tests locally:** without configuration they
  SKIP silently (and `preflight.sh` reports green).  Point CoolProp at a
  REFPROP install via `COOLPROP_REFPROP_ROOT`, or
  `COOLPROP_ALTERNATIVE_REFPROP_PATH` + `COOLPROP_ALTERNATIVE_REFPROP_LIBRARY_PATH`
  (CoolProp reads any config key from `COOLPROP_<KEY>`).
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
- **`get_superanc()` throws for pseudo-pure fluids** — it is not a
  null-returning probe.  Guard with `is_pure()`.  Pseudo-pure H,S / D+X flashes
  use the dome-free legs of `hs_cascade` without a superancillary (PR #3182).
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

## Fluid data, superancillaries and transport

- **Superancillary provenance.** fastchebpure fits superancillaries from the
  CoolProp EOS via its submodule pin, and the docs deviation plots download a
  pinned fastchebpure release.  Changing any EOS in `dev/fluids` trips the
  `superanc-pin` gate (`dev/scripts/check_superanc_release_pin.py`) until
  fastchebpure regenerates, is re-tagged, and the docs pin is bumped — or the
  fluid is added to `PENDING_UPSTREAM`.  Two independent freshness axes:
  EOS vs embedded superancillary (`source_eos_hash` + check points) and EOS vs
  docs reference.  Suspect provenance lag before the superancillary math.
- **Regenerating a superancillary:** never hand-edit `source_eos_hash`.  Use
  fastchebpure's `fitcheb` CLI: `fitcheb -f FLUID -d <coolprop checkout>` to
  fit and check, then `fitcheb inject -f FLUID -d <coolprop checkout>` to
  replace the SUPERANCILLARY block, add check points and stamp the hash.
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
