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
import math
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

# The three composition-conversion blocks. add_one parses these through the
# same parse_coefficients as the properties above, so the key requirements
# below apply to them identically -- but they are not part of the property fit
# loop, so their coefficient VALUES are not fit output and are not checked by
# the placeholder/all-zero tests.
CONVERSIONS = ["mass2input", "volume2input", "mole2input"]

# Properties that JSONIncompressibleLibrary::add_one parses with vital=true.
# For these, parse_coefficients THROWS on an unrecognised "type" (including
# "notdefined") rather than returning an empty block, and add_one rethrows --
# which aborts add_many's loop, so every fluid after the offender silently
# vanishes from the library. A cleared vital property is therefore not a
# "property you cannot query", it is a broken file.
VITAL_PROPERTIES = ("density", "specific_heat")

# The only "type" strings parse_coefficients dispatches on. Anything else --
# "notdefined", a typo, a future type this checkout predates -- means "no
# usable fit": tolerated for a non-vital property, fatal for a vital one.
# Keyed on the type string rather than on whether coeffs happen to be cleared,
# because that is what the C++ actually branches on: a vital property with a
# real-looking coeffs array and an unrecognised type throws just the same.
RECOGNISED_TYPES = ("polynomial", "exponential", "logexponential", "exppolynomial", "polyoffset")


def _assert_loader_contract(name, prop, entry):
    """Assert one block satisfies what JSONIncompressibleLibrary demands of it.

    parse_coefficients throws *unconditionally* -- irrespective of the vital
    flag -- when either the "type" or the "coeffs" key is missing from a block
    that is present, and throws for a vital property whose type it does not
    recognise. Every such throw propagates out of add_one and aborts
    add_many's loop, so the library silently truncates from that fluid onward.
    Checking presence with ``.get(..., default)`` would make an absent key
    indistinguishable from a cleared one and fail open on exactly that.
    """
    assert isinstance(entry, dict), "{0}: {1} block is {2}, not an object".format(name, prop, type(entry).__name__)
    assert "type" in entry, \
        "{0}: {1} has no \"type\" key; the C++ loader throws on that".format(name, prop)
    assert "coeffs" in entry, \
        "{0}: {1} has no \"coeffs\" key; the C++ loader throws on that".format(name, prop)
    if entry["type"] not in RECOGNISED_TYPES:
        assert prop not in VITAL_PROPERTIES, \
            "{0}: {1} is vital and has unrecognised type {2!r}; the C++ loader throws and truncates the fluid list".format(
                name, prop, entry["type"])


def _flatten(coeffs):
    """Flatten an arbitrarily nested coefficient list into a flat list of scalars.

    Block shape follows the fit *type*, not pure-vs-solution: polynomial and
    exppolynomial fits are 2D (temperature x composition) while exponential,
    logexponential and polyoffset fits are 1D, and both shapes occur in pure
    fluids and in solutions. The checks below only care about the scalar
    values, so the shape is irrelevant here.
    """
    flat = []
    for value in coeffs:
        flat.extend(_flatten(value) if isinstance(value, list) else [value])
    return flat


def _iter_property_coeffs():
    """Yield ``(filename, property, flat_coeffs)`` for every fitted property shipped.

    Walks every fluid JSON in ``JSON_DIR``. Properties that carry no
    coefficients are asserted to be typed ``notdefined`` and are then skipped,
    so the three checks below see only blocks that claim to be real fits.
    """
    files = sorted(glob.glob(os.path.join(JSON_DIR, "*.json")))
    assert files, "no fluid JSON files found in {0}".format(JSON_DIR)
    for path in files:
        with open(path) as fh:
            fluid = json.load(fh)
        for prop in PROPERTIES:
            # Require the block to exist rather than defaulting to {}: toJSON
            # always emits all of PROPERTIES, so an omitted one means data was
            # lost, and .get(prop, {}) would silently let that pass as
            # "notdefined" -- the guard failing open on exactly the case it is
            # meant to catch.
            assert prop in fluid, "{0}: missing property block {1}".format(os.path.basename(path), prop)
            entry = fluid[prop]
            _assert_loader_contract(os.path.basename(path), prop, entry)
            coeffs = entry["coeffs"]
            # Both spellings are accepted on purpose. The shipped files use the
            # STRING "null" as the cleared-coefficients placeholder -- that is
            # the pre-existing convention, not a serialization bug: every
            # entry.at("coeffs") in parse_coefficients sits inside a branch that
            # matches on "type", so for an unrecognised type the value is never
            # parsed and only its presence matters. Real JSON null behaves
            # identically; accept either so this test does not depend on which
            # is used. A cleared value with a RECOGNISED type is still a bug --
            # the loader would try to parse "null" as a coefficient array.
            if coeffs in (None, "null"):
                assert entry["type"] not in RECOGNISED_TYPES, \
                    "{0}: {1} has no coefficients but claims fitted type {2!r}".format(
                        os.path.basename(path), prop, entry["type"])
                continue
            yield os.path.basename(path), prop, _flatten(coeffs)


def test_loader_contract_satisfied():
    """Every block the C++ loader parses is loadable by it.

    Covers all nine blocks add_one feeds to parse_coefficients -- the six
    properties plus the three composition conversions -- not just the fitted
    properties. A missing "type"/"coeffs" key, or an unrecognised type on a
    vital property, makes the loader throw and silently truncate the fluid
    list; none of the coefficient-value checks below would notice.
    """
    files = sorted(glob.glob(os.path.join(JSON_DIR, "*.json")))
    assert files, "no fluid JSON files found in {0}".format(JSON_DIR)
    checked = 0
    for path in files:
        with open(path) as fh:
            fluid = json.load(fh)
        name = os.path.basename(path)
        assert isinstance(fluid, dict), "{0}: top level is {1}, not an object".format(name, type(fluid).__name__)
        for prop in PROPERTIES:
            assert prop in fluid, "{0}: missing property block {1}".format(name, prop)
            _assert_loader_contract(name, prop, fluid[prop])
            checked += 1
        for conv in CONVERSIONS:
            # An absent block is legal for these: parse_coefficients only
            # throws on a missing block when vital, and none of the three is.
            # Present-but-malformed is not legal, hence the contract check.
            if conv in fluid:
                _assert_loader_contract(name, conv, fluid[conv])
                checked += 1
    # Guard against the checks above becoming vacuous through a renamed key or
    # an empty JSON_DIR: 127 files x 9 blocks today.
    assert checked >= 9 * len(files), "only checked {0} blocks across {1} files".format(checked, len(files))


def test_no_placeholder_guesses_shipped():
    """No committed property carries an optimizer starting guess as its fit.

    An exact match against KNOWN_GUESSES means the fit never converged (or
    never ran) and the guess was serialized as though it were a result.
    """
    offenders = [(name, prop) for name, prop, flat in _iter_property_coeffs() if flat in KNOWN_GUESSES]
    assert not offenders, "unfitted optimizer starting guesses shipped as coefficients: {0}".format(offenders)


def test_no_all_zero_fits_shipped():
    """No committed property ships the all-zero polynomial template as a fit.

    An all-zero block is the untouched template: it evaluates without error to
    a constant 0 (or exp(0) = 1), which is the "plausible-looking garbage"
    failure mode this suite exists to prevent.
    """
    offenders = [(name, prop) for name, prop, flat in _iter_property_coeffs() if all(v == 0.0 for v in flat)]
    assert not offenders, "all-zero coefficient matrices shipped as fits: {0}".format(offenders)


def test_all_coefficients_finite():
    """Every shipped coefficient is a real, finite number.

    Rejects NaN and both infinities, and non-numeric JSON values -- including
    booleans, which subclass int and would otherwise pass both checks.
    """
    for name, prop, flat in _iter_property_coeffs():
        for value in flat:
            # `not isinstance(value, bool)` is load-bearing: bool subclasses
            # int, so JSON true/false would otherwise satisfy both this check
            # and the finiteness one below (math.isnan(True) is False).
            assert isinstance(value, (int, float)) and not isinstance(value, bool), \
                "{0}: {1} has non-numeric coefficient {2!r}".format(name, prop, value)
            assert not math.isnan(value) and abs(value) != float("inf"), "{0}: {1} has non-finite coefficient".format(name, prop)
