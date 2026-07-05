"""Add the optional ``*_cheb`` caloric entries to every committed json/ file.

One-shot migration tool for the Chebyshev caloric fits (see
CPIncomp/ChebyshevFits.py for the schema). Unlike a full pipeline
regeneration, this deliberately does NOT refit the polynomial entries --
numpy/scipy version drift would churn every committed coefficient. Instead:

- the polynomial entries (and every other pre-existing key) are taken from
  the committed file and re-emitted VALUE-IDENTICAL (asserted below);
- ``density_cheb``/``specific_heat_cheb`` are added, fitted from the raw
  tabular data when the fluid carries any, otherwise exact-basis-converted
  from the committed polynomial;
- files are rewritten with the writer's canonical json.dumps(indent=2,
  sort_keys=True) format, which also normalises the handful of files still
  in the old Python-2 writer format (trailing spaces, %e floats) -- a
  formatting-only change, verified value-identical.

Future full regenerations get the same entries automatically through
SolutionDataWriter.toJSON; this script exists so the *_cheb rollout does
not entangle with a full refit.

Fluids sampled from CoolProp's own EOS (Air, Acetone, Ethanol, Hexane) are
basis-converted rather than refit when the CoolProp package is absent --
identical information either way, since their committed polynomial was
fitted from the same cached grids.
"""

import glob
import json
import os
import sys

import numpy as np

from CPIncomp import getPureFluids, getSolutionFluids, getSecCoolFluids, getDigitalFluids, getExampleNames
from CPIncomp.WriterObjects import SolutionDataWriter
from CPIncomp import ChebyshevFits

JSON_DIR = os.path.join(os.path.dirname(__file__), "json")


def collect_fluid_objects():
    objects = {}
    for getter in (getPureFluids, getSolutionFluids, getDigitalFluids, getSecCoolFluids):
        for obj in getter():
            objects[obj.name] = obj
    for obj in getExampleNames(obj=True):
        objects.setdefault(obj.name, obj)
    return objects


def raw_grids(obj, prop):
    if obj is None:
        return None, None, None
    if hasattr(obj, "fitFluid") and getattr(obj, prop).data is None:
        # SecCool fluids load their grids during fitFluid(); the refit result
        # is discarded -- only the loaded raw arrays are used here.
        obj.fitFluid()
    return SolutionDataWriter.getRawGrid(obj, prop)


def deep_equal(a, b):
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(deep_equal(a[k], b[k]) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(deep_equal(x, y) for x, y in zip(a, b))
    return a == b


def main():
    objects = collect_fluid_objects()
    added, converted, skipped = [], [], []
    for path in sorted(glob.glob(os.path.join(JSON_DIR, "*.json"))):
        with open(path) as fh:
            fluid = json.load(fh)
        before = json.loads(json.dumps(fluid))  # deep copy of pre-existing state
        name = fluid["name"]
        obj = objects.get(name)

        for prop in ChebyshevFits.CALORIC_PROPERTIES:
            rawT, rawX, rawGrid = raw_grids(obj, prop)
            entry = ChebyshevFits.build_entry(fluid, prop, rawT, rawX, rawGrid)
            if entry is None:
                skipped.append((name, prop))
                continue
            fluid[prop + "_cheb"] = entry
            (added if rawGrid is not None else converted).append((name, prop))

        # every pre-existing key must be value-identical, except the *_cheb
        # entries this script owns (re-running it refreshes them)
        for key, value in before.items():
            if key.endswith("_cheb"):
                continue
            assert deep_equal(fluid[key], value), (name, key)

        with open(path, "w") as fh:
            fh.write(json.dumps(fluid, indent=2, sort_keys=True))

    print("fitted from raw data : {0} entries".format(len(added)))
    print("basis-converted      : {0} entries".format(len(converted)))
    if skipped:
        print("skipped (no usable fit): {0}".format(sorted(set(skipped))))
    return 0


if __name__ == "__main__":
    sys.exit(main())
