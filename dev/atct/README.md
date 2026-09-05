# ATcT standard enthalpies of formation

Regenerates `INFO.STANDARD_STATE` in `dev/fluids/*.json` from the Active
Thermochemical Tables.

Source: B. Ruscic and D. H. Bross, *Active Thermochemical Tables (ATcT) values
based on ver. 1.220 of the Thermochemical Network*, DOI: 10.17038/CSE/2568691.

## Re-running

    python dev/atct/fetch_atct_formation.py --version 1.220 --cache /tmp/atct.html --write

The run **fails** if any fluid's coverage differs from `expected_coverage.json`.
That is deliberate: a future ATcT version that renames or drops a species must
surface as a failure, not a silent gap. After reviewing the change, re-run with
`--update-ledger` and commit the new ledger alongside the new values.

## What is committed

- `expected_coverage.json` — every fluid and its expected state
- `atct_report.csv` — per-fluid audit trail including the source page SHA-256
- the regenerated `dev/fluids/*.json`

The 4.3 MB source page is **not** committed; use `--cache` to keep a local copy.

## Rebuild after regenerating

Fluid JSON is embedded in the library as CBOR, so values only take effect after:

    python dev/generate_headers.py
    cmake --build build_catch --target CatchTestRunner -j8

## Scope

Enthalpy only. ATcT publishes no standard entropies, so S°(298.15 K) needs a
different source and is a separate tier — see the spec, section 8.
