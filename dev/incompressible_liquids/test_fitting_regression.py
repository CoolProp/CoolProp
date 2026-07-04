"""Golden-master regression test for the incompressible-fluid fitting pipeline.

Issue #2488 ("Incompressible fitting code is no longer working") happened
because a numpy upgrade silently changed how the least-squares fits behaved,
and nobody noticed until a user reported wrong vapor pressure/viscosity
output (#2447). There was no automated test that re-ran the fit generator and
compared its output to what's actually checked in under json/*.json -- this
is that test.

It re-fits a small, representative set of fluids (one CoefficientFluids-style
pure fluid, one SolutionFluids mixture, one SecCoolFluids mixture) and checks
the regenerated coefficients are close to the committed JSON. A *tiny*
numeric drift (~1e-4 relative) is expected and tolerated -- different
numpy/scipy versions solve the same least-squares problem with slightly
different rounding. A *large* drift, a different fit type, or a crash is
exactly the class of regression #2488 needs caught before release.

Prerequisites: numpy and scipy (see requirements.txt). matplotlib and a
built CoolProp Python package are only needed for the optional report/plot
paths, not for the fitting exercised here. If numpy/scipy aren't available,
the test collects as skipped rather than failing.
"""

import json
import os

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("scipy")

from CPIncomp.WriterObjects import SolutionDataWriter
from CPIncomp import getPureFluids, getSolutionFluids, getSecCoolFluids

JSON_DIR = os.path.join(os.path.dirname(__file__), "json")

# Tolerance is numpy.allclose-style (absolute + relative): a fixed global
# absolute floor doesn't work here because a single fit's coefficients span
# many orders of magnitude (leading terms vs. high-order ones), so the
# absolute floor is tied to that property's own dominant coefficient instead
# of a universal constant -- otherwise a tiny higher-order term (e.g. DowJ's
# conductivity has terms ~1e-11) gets judged against a floor of the same
# order as itself and flakes on ordinary solver-version noise. This
# comfortably covers numpy/scipy solver-version noise (observed ~1e-4
# relative on the fluids below) while still catching an order-of-magnitude
# or wrong-fit-type regression like #2488/#2447.
RELATIVE_TOLERANCE = 1e-2
ABSOLUTE_TOLERANCE_FRACTION = 1e-6  # of the property's own largest coefficient
PROPERTIES = ["density", "specific_heat", "viscosity", "conductivity"]


def _load_json_coeffs(name, key):
    with open(os.path.join(JSON_DIR, f"{name}.json")) as fh:
        data = json.load(fh)
    entry = data.get(key, {})
    coeffs = entry.get("coeffs")
    if coeffs in (None, "null"):
        return None
    return np.array(coeffs, dtype=float)


def _assert_close_to_disk(fluidObject):
    name = fluidObject.name
    for prop in PROPERTIES:
        onDisk = _load_json_coeffs(name, prop)
        fitted = getattr(fluidObject, prop, None)
        fittedCoeffs = None if fitted is None else getattr(fitted, "coeffs", None)
        if onDisk is None or fittedCoeffs is None:
            continue  # property not defined/fitted for this fluid -- nothing to compare
        fittedCoeffs = np.array(fittedCoeffs, dtype=float)
        assert fittedCoeffs.shape == onDisk.shape, f"{name}.{prop}: shape changed, {fittedCoeffs.shape} vs {onDisk.shape}"
        magnitude = np.max(np.abs(onDisk)) if onDisk.size else 0.0
        absoluteFloor = max(magnitude * ABSOLUTE_TOLERANCE_FRACTION, 1e-12)
        tolerance = absoluteFloor + RELATIVE_TOLERANCE * np.abs(onDisk)
        worstExcess = np.max(np.abs(fittedCoeffs - onDisk) - tolerance)
        assert worstExcess <= 0, f"{name}.{prop}: refit drifted beyond tolerance (worst excess {worstExcess:.3e}) from json/{name}.json"


def test_pure_fluid_refit_matches_disk():
    fluids = {f.name: f for f in getPureFluids()}
    target = fluids["DowJ"]
    SolutionDataWriter().fitFluidList([target])
    _assert_close_to_disk(target)


def test_solution_fluid_refit_matches_disk():
    fluids = {f.name: f for f in getSolutionFluids()}
    target = fluids["LiBr"]
    SolutionDataWriter().fitFluidList([target])
    _assert_close_to_disk(target)


def test_seccool_fluid_refit_matches_disk():
    fluids = {f.name: f for f in getSecCoolFluids()}
    target = fluids["AKF"]
    SolutionDataWriter().fitSecCoolList([target])
    _assert_close_to_disk(target)
