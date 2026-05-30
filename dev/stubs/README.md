# Type stubs (`.pyi`) for the CoolProp Cython wrapper

These tools generate and verify the checked-in `.pyi` stubs that give editors
and type checkers a typed view of the compiled CoolProp extension modules.

## Why checked-in + generated

The stubs are **committed** to the repo (so any IDE sees them with no build
step) but **generated** from the Cython sources (so they can't silently drift).
`stubgen-pyx` reads the `.pyx`/`.pxd` declarations; a deterministic
post-process turns the leftover Cython/C++ types into valid Python typing.

Scope: one `CoolProp.pyi` for the whole `CoolProp.CoolProp` extension module.
`CoolProp.pyx` `include`s `HumidAirProp.pyx` and `AbstractState.pyx`, so the
module is built as a single `.so` and a single stub covers it: `AbstractState`
(~130 typed methods), `State`, the Chebyshev/SuperAncillary helpers, and the
module-level functions. The scalar/array `PropsSI`/`PhaseSI`/`HAPropsSI` family
gets hand-written `@overload`s (see `overloads_coolprop.pyi`).

## Regenerate

```bash
dev/stubs/gen_stubs.sh
```

This creates an isolated, version-pinned venv (Cython + stubgen-pyx — unrelated
to the Cython that *compiles* CoolProp), runs the generator, and writes
`wrappers/Python/CoolProp/*.pyi`. The output is deterministic.

## How sync is enforced (3 layers)

Layers 1–2 are wired as a PR gate in
`.github/workflows/python_typestubs.yml` (build-free, runs on every PR so it can
be a required check without deadlocking path-filtered PRs).  Layer 3 (runtime
parity) needs a compiled CoolProp and runs in the wheel-build workflow instead.

1. **Drift gate (exact).** CoolProp.pyi is fully generated, so CI can:
   ```bash
   dev/stubs/gen_stubs.sh
   git diff --exit-code -- 'wrappers/Python/CoolProp/*.pyi'
   ```
   A nonzero diff means a `.pyx` changed without regenerating — fail the build.
   Needs Cython 3.x (the pinned venv), which `gen_stubs.sh` provisions itself.

2. **Type correctness (static, build-free).** `dev/stubs/typecheck_smoke.py` is
   checked by both `pyright` and `mypy` (not executed) and asserts e.g.
   `AbstractState("HEOS", "Water").T()` is `float` and that `PropsSI` resolves
   2-arg->float / 6-arg-scalar->float / 6-arg-array->ndarray via the overloads.
   A type regression in the stub fails here. The checkers must resolve
   `CoolProp` from the source tree (where `CoolProp.pyi` lives) and need numpy
   stubs and a Python version with `typing.assert_type`:

   ```bash
   # run from the repo root (config paths are repo-root-relative):
   uvx --with numpy pyright --project dev/stubs/pyrightconfig.smoke.json
   uvx --with numpy mypy --ignore-missing-imports wrappers/Python/CoolProp/CoolProp.pyi
   ```

   Verified locally: pyright reads `T() -> float` from the stub (a deliberately
   wrong `assert_type(..., str)` is correctly rejected), so the pass is real,
   not vacuous.

3. **Symbol parity (needs a built CoolProp).** `dev/stubs/test_stub_parity.py`
   is a normal pytest that imports the compiled module and asserts the stub and
   the runtime class expose the same public methods. It needs a source-built
   CoolProp (a PyPI wheel could have a different surface), so it runs in the
   wheel-build workflow, not the build-free gate. (It lives in `dev/stubs/`, not
   under the CoolProp package, so pytest imports it standalone rather than
   against the installed package.)

## Post-process rules

See `postprocess.py`. In short: drop the non-importable `libcpp`/relative
imports; `string`→`str`; `vector`→`list[float]`; the enum modules
(`parameters`/`phases`/`input_pairs`)→`int` (they are integer constants at the
Python boundary); `CoolPropDbl`→`float`; unresolved `-> ...`→`-> float`;
`__cinit__`→`__init__`; drop `__dealloc__` and `@cython.*` decorators.

Hand-written content (the scalar/array `@overload`s) lives in
`overloads_coolprop.pyi`; `postprocess.py` drops the generated untyped defs that
the fragment redefines and appends it, so the result stays deterministic and the
drift gate holds.
