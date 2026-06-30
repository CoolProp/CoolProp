"""Sanity checks on the committed json/*.json fluid files.

These catch the "leaked placeholder" class of bug (GitHub #1331, #2567): a
property whose fit failed or never ran kept the optimizer's starting guess
(or the all-zero polynomial template) in its coefficient array, and that got
written to JSON as if it were a real fit. The C++ backend then evaluated
those coefficients to plausible-looking garbage instead of throwing -- LiBr
conductivity = 0 W/m/K, LiBr viscosity = exp(0) = 1 Pa.s, LiBr/MITSW
T_freeze ~ 0 K. A property without a usable fit must be typed "notdefined"
so the backend raises a clear error.

Needs no third-party packages; it only reads the committed JSON files.
"""

import glob
import json
import math
import os

JSON_DIR = os.path.join(os.path.dirname(__file__), "json")

# Optimizer starting guesses from SolutionDataWriter.fitAll / fitCoeffs.
# If a committed file carries exactly these values, the fit never happened.
#
# Read from the fitter itself rather than restated here: a copy that drifted
# would silently stop recognising a failed fit, which is the whole failure
# mode (#1331 / #2567) this check exists to prevent. Falls back to the literal
# values when numpy/CPIncomp are unavailable, so the stdlib-only checks in this
# module still run.
try:
    from CPIncomp.BaseObjects import IncompressibleData

    KNOWN_GUESSES = [
        list(IncompressibleData.VISCOSITY_GUESS),
        list(IncompressibleData.PSAT_GUESS),
        list(IncompressibleData.TFREEZE_GUESS),
        list(IncompressibleData.LOGEXP_GUESS),
    ]
except ImportError:  # pragma: no cover - exercised only without numpy installed
    KNOWN_GUESSES = [
        [500.0, -60.0, 10.0],  # viscosity, exponential
        [-5000.0, 60.0, -10.0],  # saturation pressure, exponential
        [700.0, -60.0, 10.0],  # freezing temperature, exponential
        [-250.0, 1.5, 10.0],  # log-exponential retry in IncompressibleData.fitCoeffs
    ]
PROPERTIES = ["density", "specific_heat", "conductivity", "viscosity", "saturation_pressure", "T_freeze"]

# The three composition-conversion blocks. add_one parses these through the
# same parse_coefficients as the properties, so the loader-contract check
# applies to them identically. Their coefficient VALUES are not checked: 24
# blocks already ship as an all-zero "polynomial", which the value tests would
# reject, and clearing that data is a separate change. It is inert today --
# inputFromMass/inputFromVolume throw NotImplementedError wherever a conversion
# is actually needed (the cross-basis case, which is every one of the 24).
CONVERSIONS = ["mass2input", "volume2input", "mole2input"]

# Parsed by add_one with vital=true: for these, an unrecognised "type"
# (including "notdefined") throws instead of yielding an empty block, which
# truncates the fluid list. A cleared vital property is a broken file, not an
# unqueryable property.
VITAL_PROPERTIES = ("density", "specific_heat")

# Types parse_coefficients dispatches on, and the coeffs rank each is parsed
# with (get_double_array2D vs get_double_array). Anything else means "no usable
# fit": tolerated for a non-vital property, fatal for a vital one. Note the rank
# parsers throw *unconditionally*, not gated on vital, so a recognised type with
# malformed coeffs is loader-fatal for any property -- which is why the
# structural checks below run regardless of VITAL_PROPERTIES.
RECOGNISED_TYPES = ("polynomial", "exponential", "logexponential", "exppolynomial", "polyoffset")

TYPE_NDIM = {
    "polynomial": 2,
    "exppolynomial": 2,
    "exponential": 1,
    "logexponential": 1,
    "polyoffset": 1,
}


def _assert_loader_contract(name, prop, entry):
    """Assert one block is loadable by JSONIncompressibleLibrary.

    A throw here propagates out of add_one and aborts add_many's loop, so the
    library truncates from that fluid onward with only a message on stdout.
    Presence must be tested with ``in``, not ``.get(..., default)``, which
    cannot tell an absent key from a cleared one.
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
    assert isinstance(coeffs, list), \
        "{0}: {1} claims fitted type {2!r} but coeffs is {3!r}; the loader needs an array".format(
            name, prop, entry["type"], coeffs)
    # Not loader-fatal (it yields an empty matrix), but a fit with no
    # coefficients is meaningless data.
    assert coeffs, "{0}: {1} claims fitted type {2!r} with an empty coeffs array".format(
        name, prop, entry["type"])
    # Check structure BEFORE flattening -- _flatten erases the nesting the
    # loader validates. Rectangularity matters most: num_cols() falls back to
    # max_cols() for a non-squared input, so vec_to_eigen reads past every short
    # row via unchecked operator[], fabricating a coefficient instead of
    # throwing.
    if ndim == 2:
        for index, row in enumerate(coeffs):
            assert isinstance(row, list), \
                "{0}: {1} type {2!r} is parsed as 2-D but row {3} is {4!r}; get_double_array2D throws".format(
                    name, prop, entry["type"], index, row)
        widths = {len(row) for row in coeffs}
        assert len(widths) == 1, \
            "{0}: {1} has non-rectangular coeffs (row widths {2}); vec_to_eigen reads out of bounds".format(
                name, prop, sorted(widths))
        scalars = [value for row in coeffs for value in row]
    else:
        scalars = coeffs
    for value in scalars:
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
    # Pin the key sets, not the counts: a count compared against
    # len(PROPERTIES)*len(files) is tautological, since the loop increments once
    # per name it iterates and so passes even for an empty list.
    assert set(PROPERTIES) == {"density", "specific_heat", "conductivity", "viscosity",
                               "saturation_pressure", "T_freeze"}, \
        "PROPERTIES no longer lists the six fitted properties: {0}".format(sorted(PROPERTIES))
    assert set(CONVERSIONS) == {"mass2input", "volume2input", "mole2input"}, \
        "CONVERSIONS no longer lists the three conversion blocks: {0}".format(sorted(CONVERSIONS))
    # Conversion blocks are optional per file, so demanding three per file would
    # reject legal data; a per-file floor still catches a wholesale key rename.
    assert checked_convs >= len(files), \
        "only {0} conversion blocks across {1} files; a key rename would look like this".format(
            checked_convs, len(files))


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


def test_known_guesses_match_the_fitter():
    # The stdlib fallback above must not drift from the fitter's real starting
    # guesses: if it did, this module would stop recognising an unfitted
    # placeholder and quietly pass the very check it exists to enforce.
    pytest = __import__("pytest")
    IncompressibleData = pytest.importorskip("CPIncomp.BaseObjects").IncompressibleData
    canonical = [
        list(IncompressibleData.VISCOSITY_GUESS),
        list(IncompressibleData.PSAT_GUESS),
        list(IncompressibleData.TFREEZE_GUESS),
        list(IncompressibleData.LOGEXP_GUESS),
    ]
    assert KNOWN_GUESSES == canonical, (
        "KNOWN_GUESSES has drifted from IncompressibleData's starting guesses; "
        "update the stdlib fallback in this module")
