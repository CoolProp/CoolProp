# nanobind v8 package files

Files installed into the `CoolProp` package for the **nanobind** (`COOLPROP_NANOBIND=ON`)
build — the v8 core, its capsule `State` shim, and the committed type stub.

## `CoolProp.pyi` + `py.typed` — the committed PEP 561 stub

`CoolProp.pyi` is the type stub for the compiled nanobind core. It is
**committed** (so editors/type-checkers see the typed API with no build step) and
**generated** by `nanobind`'s `stubgen` from the built `CoolProp.abi3.so`. The
wheel installs this committed copy **verbatim** — it is *not* regenerated at
wheel-build time.

### Why committed instead of generated per wheel

`stubgen` **imports** the freshly-built `CoolProp.abi3.so` to introspect it. On a
**cross-compiled** wheel the host interpreter cannot load a target-arch `.so`
(macOS x86_64 on an arm64 runner, Windows ARM64 on x64) and dies with
`ImportError: incompatible architecture`, which fails the entire wheel build
(bd CoolProp-1gas). The generated stub is architecture-independent (verified
byte-identical across macOS arm64 / Linux x86_64 with the pinned nanobind), so a
single committed copy correctly serves every wheel.

> ⚠️ **Do not edit `CoolProp.pyi` by hand.** It must stay byte-for-byte equal to
> `stubgen` output or the CI drift gate fails. Fix the bindings (or
> `stub_patterns.txt`) and regenerate.

### Regenerate

Build a **native** abi3 wheel with the regenerate flag, then copy the stub out:

```bash
# nanobind is pinned so the byte-exact drift gate is reproducible — keep this in
# sync with the pin in .github/workflows/dev_checks.yml (nanobind-wheels job).
uv pip install "nanobind==2.12.0" cython scikit-build-core ninja
COOLPROP_NANOBIND=ON uv build --wheel --no-build-isolation \
  --config-setting=cmake.define.CMAKE_BUILD_TYPE=Release \
  --config-setting=cmake.define.COOLPROP_NANOBIND_REGEN_STUB=ON \
  --out-dir dist .
unzip -o -j dist/coolprop-*-cp312-abi3-*.whl \
  CoolProp/CoolProp.pyi CoolProp/py.typed -d wrappers/Python/_nanobind/
```

`COOLPROP_NANOBIND_REGEN_STUB` runs `stubgen` and is for native builds only;
published/cross builds leave it `OFF` and ship the committed copy.

### How sync is enforced

The `nanobind-wheels` CI job (`.github/workflows/dev_checks.yml`, PR-only) builds
the cp312-abi3 wheel with `COOLPROP_NANOBIND_REGEN_STUB=ON` and diffs the
regenerated stub against this committed copy. A diff means a binding changed
without regenerating — fix by regenerating as above.

`stub_patterns.txt` is the `stubgen` pattern file applied during regeneration
(re-adds numpy imports, replaces the auto `PropsSI`/`HAPropsSI` overloads with
typed ones, etc.). The legacy **Cython** wrapper's stub is a separate pipeline —
see `dev/stubs/README.md`.
