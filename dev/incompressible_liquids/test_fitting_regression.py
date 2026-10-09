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
# All six fitted properties, not just the polynomial ones: saturation_pressure
# and T_freeze are the exponential-family fits (logexponential / exppolynomial)
# implicated in #2488, so omitting them left the golden master blind to exactly
# the fit type that motivated it. Properties a fluid does not define are skipped
# below, so the list can be exhaustive.
PROPERTIES = ["density", "specific_heat", "viscosity", "conductivity", "saturation_pressure", "T_freeze"]


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
        if onDisk is None:
            continue  # property not defined for this fluid -- nothing to compare
        # fitFluidList swallows TypeError/ValueError from a failed fit and only
        # prints, leaving coeffs unset. Skipping that case here would report a
        # green golden-master run for exactly the crash this test exists to
        # catch, so a property committed to json/ must also come back fitted.
        assert fittedCoeffs is not None, (
            f"{name}.{prop}: committed in json/{name}.json but the refit produced no coefficients "
            f"(a swallowed fit failure -- check the fitter output above)")
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


@pytest.mark.parametrize("name", ["IceEA", "IceNA", "IcePG"])
def test_seccool_ice_refit_matches_disk(name):
    """Refit each ice slurry and compare against the committed json."""
    # Issue #3303: the conductivity and viscosity committed for these three
    # fluids could not be regenerated, because their source csv tables were
    # latin-1 encoded and the read that failed on them sat inside a bare
    # try/except. A refit cleared both properties to "notdefined" and nothing
    # noticed -- the golden master above pins DowJ, LiBr and AKF only.
    #
    # Re-encoding the csvs to ASCII recovered the data: the refit reproduces
    # the committed coefficients to ~3e-7 relative, well inside the tolerance
    # above, which confirms the shipped values came from exactly this data.
    # Pin all three here so the next silent data-loading failure fails a test.
    fluids = {f.name: f for f in getSecCoolFluids()}
    target = fluids[name]
    SolutionDataWriter().fitSecCoolList([target])
    for prop in ["conductivity", "viscosity"]:
        assert _load_json_coeffs(name, prop) is not None, (
            "json/{0}.json no longer ships {1} -- if that is deliberate, drop this "
            "test rather than the data".format(name, prop))
    _assert_close_to_disk(target)


# ---------------------------------------------------------------------------
# "Does the committed fit still describe the data it was fitted to?"
#
# The golden master above pins refit-vs-json. It cannot catch a fit that was
# committed long ago from data that has since been rescaled or replaced,
# because both sides of that comparison come from the same fitter run.
# PCL is the case that motivated this: json/PCL.json shipped a viscosity a
# factor of 100 above the Paracryol table times its declared
# viscosityFactor=1e-5, i.e. 0.26 Pa s at 0 degC for a light hydrocarbon whose
# own density (792) and conductivity (0.149) put it in Therminol D-12
# territory, where 2.6 mPa s is right. The committed coefficients corresponded
# to an older factor of 1e-3.
# ---------------------------------------------------------------------------

# A fit smooths its data, so some deviation is expected and healthy. Across
# the whole corpus the worst legitimate deviation is 33% (VMA viscosity, a
# steep low-temperature curve fitted by a low-order polynomial), and the
# median is 0.08%. A committed fit that is more than 100% off its own data is
# not smoothing, it is describing something else. The threshold sits three
# times above the worst real fluid and a hundred times below the PCL failure,
# so it is not a hair trigger.
MAX_DEVIATION_FROM_SOURCE_DATA = 1.0

DATA_BACKED_PROPERTIES = ["density", "specific_heat", "conductivity", "viscosity"]


def _evaluate(kind, coeffs, T, x, Tbase, xbase):
    """Mirror of IncompressibleFluid's evaluators for the shipped fit types.

    Kept deliberately literal rather than importing the fitter, so that a
    change in the fitter cannot quietly change what this test measures.
    """
    import math

    c = np.atleast_2d(np.array(coeffs, dtype=float))
    if kind in ("polynomial", "exppolynomial"):
        dT, dx = T - Tbase, x - xbase
        value = sum(c[i, j] * dT ** i * dx ** j
                    for i in range(c.shape[0]) for j in range(c.shape[1]))
        return math.exp(value) if kind == "exppolynomial" else value
    flat = c.ravel()
    if kind == "exponential":
        # exp(c0 / (T + c1) - c2), with a pole at T = -c1
        denominator = T + flat[1]
        if abs(denominator) < 1e-8:
            return float("nan")
        return math.exp(flat[0] / denominator - flat[2])
    if kind == "logexponential":
        denominator = T + flat[0]
        if denominator < 1e-8:
            return float("nan")
        return math.exp(math.log(1.0 / denominator + 1.0 / denominator / denominator) * flat[1] + flat[2])
    return float("nan")  # notdefined, or a type this test does not model


def _worst_deviation_from_data(fluidObject, prop, entry, Tbase, xbase):
    """Largest relative gap between the committed fit and the loaded grid.

    Returns None when there is nothing to compare -- no loaded data, no
    committed coefficients, or a fit type this test does not model.
    """
    data = getattr(fluidObject, prop, None)
    if data is None or data.data is None or data.xData is None or data.yData is None:
        return None
    coeffs = entry.get("coeffs")
    if coeffs in (None, "null"):
        return None
    grid = np.atleast_2d(np.array(data.data, dtype=float))
    Ts = np.atleast_1d(np.array(data.xData, dtype=float))
    xs = np.atleast_1d(np.array(data.yData, dtype=float))
    if grid.shape != (len(Ts), len(xs)):
        return None  # a reshaped or transposed grid is not this test's business

    worst = None
    for i, T in enumerate(Ts):
        for j, x in enumerate(xs):
            measured = grid[i, j]
            if not np.isfinite(measured) or measured == 0.0:
                continue  # NaN sentinels mark points the source does not cover
            fitted = _evaluate(entry.get("type"), coeffs, T, x, Tbase, xbase)
            if not np.isfinite(fitted):
                continue
            deviation = abs(fitted - measured) / abs(measured)
            if worst is None or deviation > worst:
                worst = deviation
    return worst


def test_committed_fits_still_describe_their_source_data():
    """Every committed fit must still match the grid its fluid loads.

    See the block comment above for why the golden master cannot catch this
    and what PCL looked like when it did not.
    """
    fluids = getSecCoolFluids() + getSolutionFluids() + getPureFluids()
    offenders = []
    compared = 0
    for fluidObject in fluids:
        path = os.path.join(JSON_DIR, fluidObject.name + ".json")
        if not os.path.isfile(path):
            continue
        with open(path) as fh:
            disk = json.load(fh)
        for prop in DATA_BACKED_PROPERTIES:
            entry = disk.get(prop, {})
            worst = _worst_deviation_from_data(fluidObject, prop, entry, disk["Tbase"], disk["xbase"])
            if worst is None:
                continue
            compared += 1
            if worst > MAX_DEVIATION_FROM_SOURCE_DATA:
                offenders.append("{0}.{1}: committed fit is off its own data by {2:.0f}x".format(
                    fluidObject.name, prop, worst))

    # Without this the test passes vacuously if the loaders ever stop
    # populating .data -- which is exactly the failure mode of #3303.
    assert compared > 150, "only {0} property blocks could be compared; the data loaders look broken".format(compared)
    assert not offenders, "committed coefficients no longer match their source data:\n" + "\n".join(offenders)
