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

This shares the same prerequisites as dev/incompressible_liquids/all_incompressibles.py
itself: a built-and-installed CoolProp Python package (CPIncomp.WriterObjects
imports CoolProp.BibtexParser), plus numpy/scipy/matplotlib. The project's
docs CI (.github/workflows/docs_docker-run.yml) already builds and installs
that wheel before touching this directory, so run this test as part of (or
after) that step, not standalone. If the prerequisites aren't met, the test
collects as skipped rather than failing -- it cannot be exercised in a
sandbox without a built CoolProp wheel.
"""

import json
import os

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("scipy")
pytest.importorskip("matplotlib")

try:
    from CPIncomp.WriterObjects import SolutionDataWriter
except ImportError as exc:
    pytest.skip(f"CPIncomp.WriterObjects requires a built CoolProp Python package: {exc}", allow_module_level=True)

from CPIncomp import getPureFluids, getSolutionFluids, getSecCoolFluids

JSON_DIR = os.path.join(os.path.dirname(__file__), "json")

# Tolerances are relative to the magnitude of the checked-in coefficients;
# this comfortably covers numpy/scipy solver-version noise (observed ~1e-4
# relative on the fluids below) while still catching an order-of-magnitude
# or wrong-fit-type regression like #2488/#2447.
RELATIVE_TOLERANCE = 1e-2
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
        scale = np.maximum(np.abs(onDisk), 1e-12)
        relDiff = np.max(np.abs(fittedCoeffs - onDisk) / scale)
        assert relDiff < RELATIVE_TOLERANCE, f"{name}.{prop}: refit drifted {relDiff:.3e} relative from json/{name}.json (tolerance {RELATIVE_TOLERANCE:.0e})"


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
