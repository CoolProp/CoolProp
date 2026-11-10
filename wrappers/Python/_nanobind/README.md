# nanobind v8 package files

Files installed into the `CoolProp` package for the **nanobind** (`COOLPROP_NANOBIND=ON`)
build — the v8 core, its capsule `State` shim, and the committed type stub.

## `CoolProp.pyi`, `Plot.pyi` + `py.typed` — the committed PEP 561 stubs

`CoolProp.pyi` is the type stub for the compiled nanobind core
(`CoolProp.CoolProp`) and `Plot.pyi` the stub for its submodule
`CoolProp.CoolProp.Plot`. They are **committed** (so editors/type-checkers see
the typed API with no build step) and **generated** by `nanobind`'s `stubgen`
from the built `CoolProp.abi3.so`. The wheel installs these committed copies
**verbatim** — they are *not* regenerated at wheel-build time.

Because the module has a submodule, the wheel installs the stubs as a **stub
package** next to the extension module:

| committed file | installed as |
|---|---|
| `CoolProp.pyi` | `CoolProp/CoolProp/__init__.pyi` |
| `Plot.pyi` | `CoolProp/CoolProp/Plot.pyi` |
| `py.typed` | `CoolProp/py.typed` |

A type checker then finds `CoolProp.CoolProp.Plot` where the runtime has it. (A
stub module `CoolProp/CoolProp.pyi` cannot have submodules: its
`from . import Plot` would mean `CoolProp.Plot`, which does not exist.) The
folder holds only `.pyi` files, so it does not shadow the extension module at
runtime.

### Why committed instead of generated per wheel

`stubgen` **imports** the freshly-built `CoolProp.abi3.so` to introspect it. On a
**cross-compiled** wheel the host interpreter cannot load a target-arch `.so`
(macOS x86_64 on an arm64 runner, Windows ARM64 on x64) and dies with
`ImportError: incompatible architecture`, which fails the entire wheel build
(bd CoolProp-1gas). The generated stub is architecture-independent (verified
byte-identical across macOS arm64 / Linux x86_64 with the pinned nanobind), so a
single committed copy correctly serves every wheel.

> ⚠️ **Do not edit `CoolProp.pyi` or `Plot.pyi` by hand.** They must stay
> byte-for-byte equal to `stubgen` output or the CI drift gate fails. Fix the
> bindings (or `stub_patterns.txt`) and regenerate.

### Regenerate

Build a **native** abi3 wheel with the regenerate flag, then copy the stubs out:

```bash
# nanobind is pinned so the byte-exact drift gate is reproducible — keep this in
# sync with the pin in .github/workflows/dev_checks.yml (nanobind-wheels job).
uv pip install "nanobind==2.12.0" cython scikit-build-core ninja
COOLPROP_NANOBIND=ON uv build --wheel --no-build-isolation \
  --config-setting=cmake.define.CMAKE_BUILD_TYPE=Release \
  --config-setting=cmake.define.COOLPROP_NANOBIND_REGEN_STUB=ON \
  --out-dir dist .
WHL=$(ls dist/coolprop-*-cp312-abi3-*.whl)
unzip -p "$WHL" CoolProp/CoolProp/__init__.pyi > wrappers/Python/_nanobind/CoolProp.pyi
unzip -p "$WHL" CoolProp/CoolProp/Plot.pyi > wrappers/Python/_nanobind/Plot.pyi
unzip -p "$WHL" CoolProp/py.typed > wrappers/Python/_nanobind/py.typed
```

`COOLPROP_NANOBIND_REGEN_STUB` runs `stubgen` and is for native builds only;
published/cross builds leave it `OFF` and ship the committed copies. It imports
the built module as `CoolProp.CoolProp` from a staging package (a copy of the
module and an empty `__init__.py` in the build folder), so the stub of the
submodule names the core's types `CoolProp.CoolProp.<type>`, and writes both
stubs in one recursive run.

### How sync is enforced

The `nanobind-wheels` CI job (`.github/workflows/dev_checks.yml`, PR-only) builds
the cp312-abi3 wheel with `COOLPROP_NANOBIND_REGEN_STUB=ON` and diffs both
regenerated stubs against these committed copies. A diff means a binding changed
without regenerating — fix by regenerating as above.
`wrappers/Python/pytest/test_stub.py` checks the installed layout, the symbols
against the runtime, and with mypy that `CoolProp.CoolProp.Plot` resolves.

`stub_patterns.txt` is the `stubgen` pattern file applied during regeneration
(re-adds numpy imports, replaces the auto `PropsSI`/`HAPropsSI` overloads with
typed ones, etc.). The legacy **Cython** wrapper's stub is a separate pipeline —
see `dev/stubs/README.md`.
