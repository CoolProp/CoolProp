# nanobind interface pytest suite

Modern pytest tests for the v8 (nanobind) CoolProp Python interface. The legacy
`CoolProp/tests/` suite is mostly `nose`-based (dead on Python 3.12+), so this is
the test-verification vehicle for the nanobind migration (bd CoolProp-r9sq.5),
and it backfills the q2sh "existing suite green" acceptance check.

## What it covers

- High-level API: `PropsSI`, `HAPropsSI`
- `AbstractState` low-level surface, including the q2sh parity additions
  (idealgas/residual decompositions, `get_fluid_constant`, `neff`, the
  `set_reference_state(FluidName, *args)` dispatcher)
- module-level convenience wrappers (`FluidsList`, `get_aliases`,
  `get_REFPROPname`, `get_BibTeXKey`, `get_config_int`)
- the capsule `State` shim unit conventions (kPa/kJ getters, g/mol `get_MM`,
  SI `.pAS`) — the split the original PDSim regression hinged on
- the approved v8 removals (`Props`/`HAProps` absent; `HAProps` raises)

The v8-specific removal tests `skipif` the legacy build, so the file is safe to
point at either wheel.

## Running

Build and install a nanobind wheel, then run pytest against it:

```bash
# from the repo root, in a venv with: nanobind, cython, scikit-build-core, uv, pytest
uv build --wheel --no-build-isolation \
    --config-setting=cmake.define.COOLPROP_NANOBIND=ON --out-dir dist .
pip install --force-reinstall dist/coolprop-*.whl
pytest wrappers/Python/pytest -q
```
