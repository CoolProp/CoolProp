"""Sanity checks on the committed json/*.json fluid files.

These catch the "leaked placeholder" class of bug (GitHub #1331, #2567): a
property whose fit failed or never ran kept the optimizer's starting guess
(or the all-zero polynomial template) in its coefficient array, and that got
written to JSON as if it were a real fit. The C++ backend then evaluated
those coefficients to plausible-looking garbage instead of throwing -- LiBr
conductivity = 0 W/m/K, LiBr viscosity = exp(0) = 1 Pa.s, LiBr/MITSW
T_freeze ~ 0 K. A property without a usable fit must be typed "notdefined"
so the backend raises a clear error.

Scope limit worth knowing about: the coefficient-VALUE checks cover the six
fitted properties only, not the three composition-conversion blocks. 24 of
those currently ship as an all-zero "polynomial" -- the very shape this module
rejects elsewhere -- because the SecCool writer fits them with the standard
template, swallows the failure, and clearUnfittedCoefficients does not cover
them. That data is inert today (the inputFromMass/inputFromVolume evaluation
switch is commented out behind a NotImplementedError), so cleaning it up is a
separate data change; the loader-contract check below does apply to those
blocks. See the CONVERSIONS comment for detail.

Unlike test_fitting_regression.py this needs no third-party packages at all;
it only reads the committed JSON files.
"""

import glob
import json
import math
import os
import re

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
# same parse_coefficients as the properties above, so the loader-contract
# requirements below apply to them identically.
#
# Their coefficient VALUES are deliberately NOT checked by the
# placeholder/all-zero tests, and not because they are hand-written constants
# -- they are fit output like everything else (SecCoolSolutionData fits
# mass2input/volume2input with the same all-zero template and swallows the
# failure), and clearUnfittedCoefficients does not cover them, so 24 blocks
# currently ship as type "polynomial" with an all-zero 4x6 matrix and
# NRMS null: exactly the leaked-template shape this module exists to reject.
# Extending the value checks here would fail on that pre-existing data, which
# needs regenerating or clearing first. It is currently inert -- the
# inputFromMass/inputFromVolume evaluation switch is commented out behind a
# NotImplementedError -- but it turns 24 fluids into "any composition maps to
# 0" the moment that is uncommented. Tracked separately; see the module
# docstring.
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

# Coefficient rank each recognised type is parsed with: get_double_array2D for
# the 2-D ones, get_double_array for the 1-D ones. Both throw
# *unconditionally* -- not gated on vital -- when the value is not an array of
# the right rank containing only numbers, so a recognised type with malformed
# coeffs is loader-fatal for any property, vital or not.
TYPE_NDIM = {
    "polynomial": 2,
    "exppolynomial": 2,
    "exponential": 1,
    "logexponential": 1,
    "polyoffset": 1,
}


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
    # cpjson::get_string throws on a non-string "type" before any dispatch, so
    # this is fatal for every property, not just the vital ones.
    assert isinstance(entry["type"], str), \
        "{0}: {1} type is {2}, not a string; cpjson::get_string throws".format(
            name, prop, type(entry["type"]).__name__)
    ndim = TYPE_NDIM.get(entry["type"])
    if ndim is None:
        # Unrecognised type: no branch matches, so coeffs is never parsed and
        # its value is free. Fatal only when the property is vital.
        assert prop not in VITAL_PROPERTIES, \
            "{0}: {1} is vital and has unrecognised type {2!r}; the C++ loader throws and truncates the fluid list".format(
                name, prop, entry["type"])
        return
    # Recognised type => coeffs really is parsed, at a rank fixed by the type.
    coeffs = entry["coeffs"]
    assert isinstance(coeffs, list) and coeffs, \
        "{0}: {1} claims fitted type {2!r} but coeffs is {3!r}; the loader needs a non-empty array".format(
            name, prop, entry["type"], coeffs)
    is_nested = isinstance(coeffs[0], list)
    assert is_nested == (ndim == 2), \
        "{0}: {1} type {2!r} is parsed as {3}-D but coeffs is {4}-D; the loader throws".format(
            name, prop, entry["type"], ndim, 2 if is_nested else 1)
    for value in _flatten(coeffs):
        assert isinstance(value, (int, float)) and not isinstance(value, bool), \
            "{0}: {1} has non-numeric coefficient {2!r}; the loader throws".format(name, prop, value)


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


def test_recognised_types_match_the_loader():
    """RECOGNISED_TYPES/TYPE_NDIM still match the C++ dispatch they duplicate.

    Staleness here fails *closed*, not open: a type C++ gained but this module
    lacks is treated as "no usable fit", which for a vital property asserts and
    blames the data for what is really an out-of-date list. PR 4/5 of this stack
    adds a ``chebyshev`` type, so the divergence is scheduled, not theoretical.
    Cross-checking turns that into one obvious failure telling you to update the
    tuple.
    """
    loader = os.path.join(os.path.dirname(__file__), "..", "..", "src", "Backends", "Incompressible",
                          "IncompressibleLibrary.cpp")
    assert os.path.exists(loader), "loader not found at {0}".format(loader)
    with open(loader) as fh:
        source = fh.read()
    start = source.find("IncompressibleData JSONIncompressibleLibrary::parse_coefficients")
    assert start != -1, "could not locate parse_coefficients in {0}".format(loader)
    body = source[start:source.index("\n}", start)]
    found = tuple(re.findall(r'type\.compare\("([^"]+)"\)', body))
    assert found, "could not parse the loader's type dispatch; update this test alongside it"
    assert set(found) == set(RECOGNISED_TYPES), \
        "loader dispatches on {0} but RECOGNISED_TYPES has {1}; update RECOGNISED_TYPES and TYPE_NDIM".format(
            sorted(found), sorted(RECOGNISED_TYPES))
    assert set(TYPE_NDIM) == set(RECOGNISED_TYPES), \
        "TYPE_NDIM covers {0} but RECOGNISED_TYPES has {1}".format(sorted(TYPE_NDIM), sorted(RECOGNISED_TYPES))


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
    checked_props = 0
    checked_convs = 0
    for path in files:
        with open(path) as fh:
            fluid = json.load(fh)
        name = os.path.basename(path)
        assert isinstance(fluid, dict), "{0}: top level is {1}, not an object".format(name, type(fluid).__name__)
        for prop in PROPERTIES:
            assert prop in fluid, "{0}: missing property block {1}".format(name, prop)
            _assert_loader_contract(name, prop, fluid[prop])
            checked_props += 1
        for conv in CONVERSIONS:
            # An absent block is legal for these: parse_coefficients only
            # throws on a missing block when vital, and none of the three is.
            # Present-but-malformed is not legal, hence the contract check.
            if conv in fluid:
                _assert_loader_contract(name, conv, fluid[conv])
                checked_convs += 1
    # Guard against the loop above going vacuous through a renamed key. Counted
    # in two halves on purpose: the properties are mandatory so their count is
    # exact, while conversion blocks are optional, so folding them into one
    # total would demand all three in every file and REJECT LEGAL DATA.
    assert checked_props == len(PROPERTIES) * len(files), \
        "checked {0} property blocks, expected {1}".format(checked_props, len(PROPERTIES) * len(files))
    assert checked_convs > 0, "no conversion blocks checked across {0} files".format(len(files))


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
