#!/usr/bin/env python3
"""Validate CoolProp's source JSON data files against their JSON schemas.

This is the build-time correctness gate for embedded fluid data
(RapidJSON->nlohmann migration spec, section 5). Run from the repo root.
Exits non-zero on the first validation failure.

The committed schemas (``*_schema.json`` under ``dev/``) all declare a
top-level ``type: array`` and describe the WHOLE data document (an array of
fluids / departure functions), mirroring how the C++ loaders use them:
e.g. ``cpjson::validate_schema(pcsaft_fluids_schema_JSON, JSON, errstr)`` in
``src/Backends/PCSAFT/PCSAFTLibrary.cpp`` validates the entire blob, not each
item. We therefore validate each data document as a whole against its schema.
"""
import json
import sys
from pathlib import Path

try:
    import jsonschema
except ImportError:
    sys.exit("jsonschema is required: pip install jsonschema")

REPO = Path(__file__).resolve().parent.parent

# (data file, schema file) pairs. Each schema declares a top-level
# ``type: array`` and validates the whole corresponding data document.
PAIRS = [
    (REPO / "dev/pcsaft/all_pcsaft_fluids.json", REPO / "dev/pcsaft/pcsaft_fluids_schema.json"),
    (REPO / "dev/cubics/all_cubic_fluids.json", REPO / "dev/cubics/cubic_fluids_schema.json"),
    (REPO / "dev/mixtures/mixture_departure_functions.json",
     REPO / "dev/mixtures/mixture_departure_functions_schema.json"),
]


def validate_pair(data_path: Path, schema_path: Path) -> int:
    if not data_path.exists() or not schema_path.exists():
        print(f"SKIP (missing): {data_path.name} / {schema_path.name}")
        return 0
    schema = json.loads(schema_path.read_text())
    data = json.loads(data_path.read_text())
    n_items = len(data) if isinstance(data, (list, dict)) else 1
    try:
        jsonschema.validate(instance=data, schema=schema)
    except jsonschema.ValidationError as e:
        loc = "/".join(str(p) for p in e.absolute_path) or "<root>"
        print(f"FAIL {data_path.name} at {loc}: {e.message}")
        return 1
    print(f"OK  {data_path.name} ({n_items} items) vs {schema_path.name}")
    return 0


def main() -> int:
    total = sum(validate_pair(d, s) for d, s in PAIRS)
    if total:
        print(f"\n{total} schema validation failure(s)")
        return 1
    print("\nAll fluid data files validate against their schemas.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
