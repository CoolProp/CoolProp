"""Checks on what SolutionDataWriter actually emits.

Two separate outputs, two separate hazards:

1. json/*.json -- the fitted coefficients. json.dumps writes all 17 decimal
   digits of a double, which is far more than the source data supports (the
   SecCool tables carry four). Without rounding, re-running the pipeline on a
   different numpy/scipy rewrites every coefficient in all 116 files with
   last-digit noise, and a real change cannot be told apart from that churn.
   That churn is what let issue #3303 sit unnoticed.

2. The composition limits published in the online fluid tables. Users read
   xmin/xmax off those tables and feed them back into set_mass_fractions(),
   so a limit printed outside the real range turns into a ValueError from
   the backend (issue #2567).

Needs numpy and scipy (WriterObjects imports both); skips cleanly without.
"""

import glob
import json
import math
import os

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("scipy")

from CPIncomp.WriterObjects import SolutionDataWriter, roundNestedNumbers, roundToSignificantDigits, SIGNIFICANT_DIGITS

JSON_DIR = os.path.join(os.path.dirname(__file__), "json")

# IncompressibleFluid::checkX in the C++ backend, which is what a user hits
# after reading a limit off the published table:
#   x < xmin * (1 - eps) or x > xmax * (1 + eps)  ->  ValueError
# Note the tolerance is relative to the bound, so it buys nothing at all when
# the bound is 0.0 and very little when it is small.
INCOMP_EPSILON = 1e-6


def _significant_digits(value):
    """Count the significant decimal digits in value's shortest repr."""
    text = repr(float(value))
    if "e" in text or "E" in text:
        text = text.split("e")[0].split("E")[0]
    digits = text.lstrip("-").replace(".", "").lstrip("0")
    return len(digits.rstrip("0")) or 1


def test_round_to_significant_digits_handles_every_magnitude():
    assert roundToSignificantDigits(0.5582565205627504) == 0.5582565
    assert roundToSignificantDigits(1.3048060293070506e-08) == 1.304806e-08
    assert roundToSignificantDigits(278.27819999999997) == 278.2782
    assert roundToSignificantDigits(-67.98607219465411) == -67.98607
    # Zero has no exponent to take out, and log10(0) would raise.
    assert roundToSignificantDigits(0.0) == 0.0
    # Non-finite values must pass through rather than blow up: an unfittable
    # property can carry NaN, and the caller decides what to do with it.
    assert math.isnan(roundToSignificantDigits(float("nan")))
    assert math.isinf(roundToSignificantDigits(float("inf")))


def test_round_nested_numbers_leaves_non_numbers_alone():
    """The null placeholders must survive: json.dumps(NaN) is not valid JSON.

    toJSON writes the string "null" for absent coefficients and None for an
    absent NRMS. Rounding must not coerce either into a number -- an earlier
    attempt turned "null" coefficients into NaN, which json.dumps happily
    wrote out as a bare NaN token that no strict JSON parser accepts.
    """
    payload = {
        "type": "notdefined",
        "coeffs": "null",
        "NRMS": None,
        "nested": [{"a": 1.23456789012}, ["null", None, 2.98765432109]],
        "flag": True,  # bool is a subclass of int -- must not become 1.0
    }
    out = roundNestedNumbers(payload)
    assert out["type"] == "notdefined"
    assert out["coeffs"] == "null"
    assert out["NRMS"] is None
    assert out["nested"][0]["a"] == 1.234568
    assert out["nested"][1][0] == "null"
    assert out["nested"][1][1] is None
    assert out["nested"][1][2] == 2.987654
    assert out["flag"] is True
    # The whole point: the result must serialise as strict JSON.
    json.loads(json.dumps(out))


def test_round_nested_numbers_caps_the_digits_it_emits():
    noisy = {"coeffs": [[0.5582565205627504, 1.3048060293070506e-08],
                        [-67.98607219465411, 948.7603675742353]]}
    for value in np.array(roundNestedNumbers(noisy)["coeffs"]).ravel():
        assert _significant_digits(value) <= SIGNIFICANT_DIGITS, value


@pytest.mark.parametrize("path", sorted(glob.glob(os.path.join(JSON_DIR, "*.json"))),
                         ids=lambda p: os.path.splitext(os.path.basename(p))[0])
def test_published_composition_range_is_a_valid_range(path):
    """Sweeping the documented limits must not trip the backend's own check.

    Issue #2567 reported T_freeze raising near the composition limits for
    MAM2, VCA, VKC and VMG. Those fluids are fine -- the fits evaluate
    sensibly across their whole real range. What was wrong is the published
    range: MAM2 stops at 0.236 and was printed as 0.24, so the documented
    upper limit is a composition the backend rejects.
    """
    with open(path) as fh:
        fluid = json.load(fh)
    if fluid.get("xid") not in ("mass", "volume"):
        return  # pure fluids carry no composition axis
    xmin, xmax = fluid["xmin"], fluid["xmax"]
    if xmin == xmax:
        return

    writer = SolutionDataWriter()
    shownMin = float(writer.x(xmin))
    shownMax = float(writer.x(xmax))

    name = os.path.basename(path)
    assert shownMin >= xmin * (1 - INCOMP_EPSILON), (
        "{0}: published lower limit {1} is below the real {2}".format(name, shownMin, xmin))
    assert shownMax <= xmax * (1 + INCOMP_EPSILON), (
        "{0}: published upper limit {1} is above the real {2}".format(name, shownMax, xmax))
    # And every point in between, the way a user would sweep it.
    for x in np.linspace(shownMin, shownMax, 50):
        assert xmin * (1 - INCOMP_EPSILON) <= x <= xmax * (1 + INCOMP_EPSILON), (
            "{0}: sweeping the published range reaches {1}, outside [{2}, {3}]".format(
                name, x, xmin, xmax))
