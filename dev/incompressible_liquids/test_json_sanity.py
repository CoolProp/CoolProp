"""Sanity checks on the committed json/*.json fluid files.

These catch the "leaked placeholder" class of bug (GitHub #1331, #2567): a
property whose fit failed or never ran kept the optimizer's starting guess
(or the all-zero polynomial template) in its coefficient array, and that got
written to JSON as if it were a real fit. The C++ backend then evaluated
those coefficients to plausible-looking garbage instead of throwing -- LiBr
conductivity = 0 W/m/K, LiBr viscosity = exp(0) = 1 Pa.s, LiBr/MITSW
T_freeze ~ 0 K. A property without a usable fit must be typed "notdefined"
so the backend raises a clear error.

Unlike test_fitting_regression.py this needs no third-party packages at all;
it only reads the committed JSON files.
"""

import glob
import json
import os

JSON_DIR = os.path.join(os.path.dirname(__file__), "json")

# Optimizer starting guesses from SolutionDataWriter.fitAll / fitCoeffs.
# If a committed file carries exactly these values, the fit never happened.
KNOWN_GUESSES = [
    [500.0, -60.0, 10.0],  # viscosity, exponential
    [-5000.0, 60.0, -10.0],  # saturation pressure, exponential
    [700.0, -60.0, 10.0],  # freezing temperature, exponential
    [-250.0, 1.5, 10.0],  # log-exponential retry in IncompressibleData.fitCoeffs
]
PROPERTIES = ["density", "specific_heat", "conductivity", "viscosity", "saturation_pressure", "T_freeze"]


def _flatten(coeffs):
    flat = []
    for value in coeffs:
        flat.extend(_flatten(value) if isinstance(value, list) else [value])
    return flat


def _iter_property_coeffs():
    files = sorted(glob.glob(os.path.join(JSON_DIR, "*.json")))
    assert files, "no fluid JSON files found in {0}".format(JSON_DIR)
    for path in files:
        with open(path) as fh:
            fluid = json.load(fh)
        for prop in PROPERTIES:
            entry = fluid.get(prop, {})
            coeffs = entry.get("coeffs")
            if coeffs in (None, "null"):
                assert entry.get("type", "notdefined") == "notdefined", \
                    "{0}: {1} has no coefficients but type {2!r}".format(os.path.basename(path), prop, entry.get("type"))
                continue
            yield os.path.basename(path), prop, _flatten(coeffs)


def test_no_placeholder_guesses_shipped():
    offenders = [(name, prop) for name, prop, flat in _iter_property_coeffs() if flat in KNOWN_GUESSES]
    assert not offenders, "unfitted optimizer starting guesses shipped as coefficients: {0}".format(offenders)


def test_no_all_zero_fits_shipped():
    offenders = [(name, prop) for name, prop, flat in _iter_property_coeffs() if all(v == 0.0 for v in flat)]
    assert not offenders, "all-zero coefficient matrices shipped as fits: {0}".format(offenders)


def test_all_coefficients_finite():
    for name, prop, flat in _iter_property_coeffs():
        for value in flat:
            assert isinstance(value, (int, float)), "{0}: {1} has non-numeric coefficient {2!r}".format(name, prop, value)
            assert value == value and abs(value) != float("inf"), "{0}: {1} has non-finite coefficient".format(name, prop)
