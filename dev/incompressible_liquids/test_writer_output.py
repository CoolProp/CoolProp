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
import sys

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("scipy")

from CPIncomp.WriterObjects import SolutionDataWriter, roundNestedNumbers, roundToSignificantDigits, SIGNIFICANT_DIGITS

JSON_DIR = os.path.join(os.path.dirname(__file__), "json")

# IncompressibleFluid::checkX in the C++ backend, which is what a user hits
# after reading a limit off the published table:
#   x < xmin * (1 - eps) or x > xmax * (1 + eps)  ->  ValueError
# This must be the backend's own value, from IncompressibleFluid.cpp:
#   constexpr double INCOMP_EPSILON = DBL_EPSILON * 100.0;
# An earlier draft used 1e-6 here, which is around 4.5e7 times looser, so a
# published limit up to 1e-6 outside the enforced range would have passed this
# test and still been rejected at runtime. Note the tolerance is relative to
# the bound, so it buys nothing at all when the bound is 0.0.
INCOMP_EPSILON = sys.float_info.epsilon * 100.0


def _significant_digits(value):
    """Count the significant decimal digits in value's shortest repr."""
    text = repr(float(value))
    if "e" in text or "E" in text:
        text = text.split("e")[0].split("E")[0]
    digits = text.lstrip("-").replace(".", "").lstrip("0")
    return len(digits.rstrip("0")) or 1


def test_round_to_significant_digits_handles_every_magnitude():
    """Rounding must work across the exponent range, and pass non-finites through."""
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
    """No value that survives rounding may carry more than 7 significant digits."""
    noisy = {"coeffs": [[0.5582565205627504, 1.3048060293070506e-08],
                        [-67.98607219465411, 948.7603675742353]]}
    for value in np.array(roundNestedNumbers(noisy)["coeffs"]).ravel():
        assert _significant_digits(value) <= SIGNIFICANT_DIGITS, value


JSON_PATHS = sorted(glob.glob(os.path.join(JSON_DIR, "*.json")))

# pytest's default empty_parameter_set_mark is "skip", so an empty JSON_PATHS
# would collapse the test below into a single skipped item and the run would
# still be green. Fail at import instead.
assert JSON_PATHS, "no files found in {0} -- the composition-range test would silently skip".format(JSON_DIR)


@pytest.mark.parametrize("path", JSON_PATHS,
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
    # skip, not return: a bare return reports these as passes, so if the "xid"
    # key were ever renamed all 127 cases would go vacuous with no signal.
    if fluid.get("xid") not in ("mass", "volume"):
        pytest.skip("pure fluid, no composition axis")
    xmin, xmax = fluid["xmin"], fluid["xmax"]
    if xmin == xmax:
        pytest.skip("single-composition fluid")

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


def test_write_fluid_list_raises_when_a_fluid_cannot_be_serialised(monkeypatch, tmp_path):
    """A fluid that fails to serialise must not be reported as written.

    toJSON can now raise (allow_nan=False), and writeFluidList used to log
    that and carry on, printing "done" and returning normally while
    json/<name>.json stayed missing or stale. Shipping a stale coefficient
    file while the pipeline claims success is exactly the failure this branch
    exists to remove, so the batch writer has to surface it.
    """
    writer = SolutionDataWriter()

    class FakeFluid(object):
        def __init__(self, name):
            self.name = name

    # BadOne comes FIRST on purpose. With the healthy fluid first, a fail-fast
    # writeFluidList that re-raised immediately would still have written
    # GoodOne before raising, and the assertions below would pass -- the test
    # would then be pinning "it raises" without pinning "it keeps going".
    # Putting the failure first means only an implementation that continues
    # past it can reach GoodOne at all.
    fluids = [FakeFluid("BadOne"), FakeFluid("GoodOne")]
    written = []

    def fakeToJSON(data, quiet=False):
        if data.name == "BadOne":
            raise ValueError("Out of range float values are not JSON compliant")
        written.append(data.name)

    monkeypatch.setattr(writer, "toJSON", fakeToJSON)
    monkeypatch.setattr(writer, "printStatusID", lambda objs, obj: None)

    with pytest.raises(ValueError, match="BadOne"):
        writer.writeFluidList(fluids)

    # This is the assertion that distinguishes aggregation from fail-fast, and
    # it must not rely on the raised message happening to name the fluid.
    assert written == ["GoodOne"], (
        "the healthy fluid after the failure was not written, so writeFluidList "
        "stopped at the first error instead of attempting every fluid")


# IAPWS-08 seawater freezing temperature at 0.1013 MPa, across MITSW's own
# salinity range. Generated once with iapws.iapws08._Tf and committed here so
# the check needs no runtime dependency and cannot silently skip when a
# package is absent. Regenerate with:
#     from iapws.iapws08 import _Tf; _Tf(0.1013, x)
IAPWS08_SEAWATER_FREEZE = [
    (0.00, 273.1525), (0.01, 272.6139), (0.02, 272.0748), (0.03, 271.5225),
    (0.04, 270.9536), (0.05, 270.3652), (0.06, 269.7542), (0.07, 269.1178),
    (0.08, 268.4532), (0.09, 267.7583), (0.10, 267.0318), (0.11, 266.2730),
    (0.12, 265.4820),
]

# The shipped curve is a cubic, deliberately smoothed so it never crosses Tmin
# (see the long comment on MITSeaWater in CPIncomp/SolutionFluids.py). Its
# worst departure from IAPWS-08 is ~10 mK, at zero salinity. 20 mK leaves room
# for a coefficient re-fit without being loose enough to miss a real error: a
# wrong sign, a dropped term or a stale xbase all move it by whole kelvin.
SEAWATER_FREEZE_TOLERANCE_K = 0.02


def _evaluate_T_freeze(entry, x, xbase):
    """Mirror of IncompressibleFluid::Tfreeze for a single-row coefficient set.

    The backend evaluates poly.evaluate(coeffs, p, x, 0, 0, 0.0, xbase), so
    rows are powers of pressure and columns are powers of (x - xbase). One row
    means no pressure dependence, which is what this fluid ships.
    """
    coeffs = np.atleast_2d(np.array(entry["coeffs"], dtype=float))
    assert coeffs.shape[0] == 1, "expected no pressure dependence, got {0}".format(coeffs.shape)
    return sum(coeffs[0, j] * (x - xbase) ** j for j in range(coeffs.shape[1]))


def test_mitsw_freezing_curve_matches_iapws08():
    """MITSW's committed T_freeze must be the IAPWS-08 seawater freeze curve.

    Issue #2567 asked for this. Sharqawy 2010 supplies every other MITSW
    property but publishes no freezing temperature, and neither does the MIT
    seawater library, so IAPWS-08 is the source. Checking the committed
    coefficients against IAPWS values directly keeps this honest: refitting
    and comparing to the same coefficients would prove nothing.
    """
    with open(os.path.join(JSON_DIR, "MITSW.json")) as fh:
        fluid = json.load(fh)
    entry = fluid["T_freeze"]

    assert entry["type"] == "polynomial", (
        "MITSW T_freeze is {0}; if it was deliberately removed, delete this test "
        "rather than loosening it".format(entry["type"]))

    worst = 0.0
    for x, expected in IAPWS08_SEAWATER_FREEZE:
        got = _evaluate_T_freeze(entry, x, fluid["xbase"])
        worst = max(worst, abs(got - expected))
    assert worst <= SEAWATER_FREEZE_TOLERANCE_K, (
        "MITSW T_freeze departs from IAPWS-08 by {0:.4f} K, tolerance {1} K".format(
            worst, SEAWATER_FREEZE_TOLERANCE_K))


def test_mitsw_freezing_curve_never_crosses_tmin():
    """The freeze curve must stay at or below Tmin, or checkT starts rejecting.

    IncompressibleFluid::checkT throws when T < T_freeze(x), on top of the
    Tmin bound. The true freezing point of pure water at 1 atm is 273.1525 K,
    slightly ABOVE MITSW's Tmin of 273.15 K, because 273.15 K is the
    air-saturated ice point. So a curve fitted too closely to IAPWS-08 at
    x = 0 would make PropsSI reject T = Tmin at zero salinity, which works
    today. The shipped cubic is smoothed to avoid exactly that, and this is
    the test that stops someone "improving" the fit and breaking it.
    """
    with open(os.path.join(JSON_DIR, "MITSW.json")) as fh:
        fluid = json.load(fh)
    entry = fluid["T_freeze"]
    if entry.get("coeffs") in (None, "null"):
        pytest.skip("MITSW ships no T_freeze")

    # Endpoints plus the real roots of the derivative, NOT a sampling grid. A
    # polynomial attains its extremes on an interval only at an endpoint or a
    # critical point, so this is exhaustive; a grid is not. An earlier version
    # sampled 241 points, and a parabola peaked exactly between two samples
    # passes that at 273.0875 K while actually reaching 274.15 K, a whole
    # kelvin above Tmin. The point of this test is to constrain coefficients
    # nobody has written yet, so it has to hold for any polynomial, not just
    # the well-behaved cubic shipped today.
    coeffs = np.atleast_2d(np.array(entry["coeffs"], dtype=float))[0]
    derivativeRoots = np.polynomial.Polynomial(coeffs).deriv().roots()
    interior = [root.real + fluid["xbase"] for root in np.atleast_1d(derivativeRoots)
                if np.isclose(root.imag, 0.0)
                and fluid["xmin"] <= root.real + fluid["xbase"] <= fluid["xmax"]]
    xs = np.array([fluid["xmin"], fluid["xmax"]] + interior)

    values = np.array([_evaluate_T_freeze(entry, x, fluid["xbase"]) for x in xs])
    worstIndex = int(np.argmax(values))
    assert values[worstIndex] <= fluid["Tmin"], (
        "T_freeze reaches {0:.6f} K at x = {1:.4f}, above Tmin = {2}. checkT would "
        "reject T = Tmin there.".format(values[worstIndex], xs[worstIndex], fluid["Tmin"]))


# ---------------------------------------------------------------------------
# Two guards for the rounding contract itself.
#
# Both exist because an adversarial review of the regeneration found that the
# contract was asserted in comments and demonstrated once by hand, with
# nothing mechanical holding it. Reverting the writer fix entirely, while
# keeping the regenerated json, left the whole suite green.
# ---------------------------------------------------------------------------

def _significant_digit_count(value):
    text = repr(float(value))
    if "e" in text or "E" in text:
        text = text.split("e")[0].split("E")[0]
    digits = text.lstrip("-").replace(".", "").lstrip("0")
    return len(digits.rstrip("0")) or 1


def _walk_numbers(obj, path=""):
    if isinstance(obj, dict):
        for k, v in obj.items():
            for r in _walk_numbers(v, path + "/" + str(k)):
                yield r
    elif isinstance(obj, (list, tuple)):
        for v in obj:
            for r in _walk_numbers(v, path + "[]"):
                yield r
    elif isinstance(obj, float) or (isinstance(obj, int) and not isinstance(obj, bool)):
        yield path, float(obj)


# Acetone, Air, Ethanol and Hexane are DigitalFluids generated from the
# CoolProp HEOS backend, so regenerating them needs a built CoolProp Python
# package. They predate the rounding and are still committed at full
# precision. They are excluded here rather than hand-edited, because rounding
# a polynomial without re-deriving its Chebyshev conversion is exactly the bug
# this module's other tests exist to catch. The set is pinned below so it can
# neither grow silently nor go stale once someone regenerates them.
UNROUNDED_LEGACY_FLUIDS = frozenset({"Acetone", "Air", "Ethanol", "Hexane"})


def _digit_cap_offenders(path):
    """Numbers in one file that carry more digits than their rule allows.

    Returns (how many numbers were examined, the offending ones). The count
    is returned so a caller can tell "nothing was over the cap" apart from
    "nothing was looked at".
    """
    from CPIncomp import ChebyshevFits

    with open(path) as fh:
        fluid = json.load(fh)

    exactCoeffPaths = set()
    for prop in ChebyshevFits.CALORIC_PROPERTIES:
        entry = fluid.get(prop + "_cheb")
        if entry and entry.get("fit_source") in ChebyshevFits.EXACT_FIT_SOURCES:
            exactCoeffPaths.add("/" + prop + "_cheb/coeffs")

    offenders, checked = [], 0
    for where, value in _walk_numbers(fluid):
        checked += 1
        cap = (ChebyshevFits.EXACT_CONVERSION_DIGITS
               if where.replace("[]", "") in exactCoeffPaths else SIGNIFICANT_DIGITS)
        if _significant_digit_count(value) > cap:
            offenders.append("{0}: {1!r} has {2} significant digits, cap {3}".format(
                where, value, _significant_digit_count(value), cap))
    return checked, offenders


def test_unrounded_legacy_files_are_exactly_the_known_set():
    """The exemption list must match reality, in both directions.

    If a fifth file starts shipping full-precision numbers this fails, and if
    one of the four is finally regenerated this also fails, prompting its
    removal from the list. An exemption nobody revisits is how a temporary
    carve-out becomes permanent.
    """
    offending = set()
    for path in JSON_PATHS:
        _, offenders = _digit_cap_offenders(path)
        if offenders:
            offending.add(os.path.splitext(os.path.basename(path))[0])
    assert offending == set(UNROUNDED_LEGACY_FLUIDS), (
        "files exceeding the digit cap are {0}, expected exactly {1}".format(
            sorted(offending), sorted(UNROUNDED_LEGACY_FLUIDS)))


@pytest.mark.parametrize("path", JSON_PATHS,
                         ids=lambda p: os.path.splitext(os.path.basename(p))[0])
def test_committed_numbers_respect_the_digit_cap(path):
    """No committed number may carry more digits than its rule allows.

    The rounding exists so that regenerating json/ on a different numpy or
    scipy does not rewrite every coefficient in every file. That property is
    only real if it holds for the FILES, not just for the helper: a value
    added to the writer's dict after the rounding pass would ship at 17
    digits and nothing would notice. This is the check that notices.

    One documented exception: an exact Chebyshev conversion is held to
    ChebyshevFits.EXACT_CONVERSION_TOLERANCE against the polynomial it comes
    from, which 7 digits cannot deliver, so it gets EXACT_CONVERSION_DIGITS.
    """
    if os.path.splitext(os.path.basename(path))[0] in UNROUNDED_LEGACY_FLUIDS:
        pytest.skip("legacy full-precision file; pinned by "
                    "test_unrounded_legacy_files_are_exactly_the_known_set")

    checked, offenders = _digit_cap_offenders(path)

    # A fluid file always carries numbers; zero would mean this walked nothing.
    assert checked > 0, "no numbers found in {0}".format(path)
    assert not offenders, "{0}:\n  {1}".format(
        os.path.basename(path), "\n  ".join(offenders[:10]))


def _cheb_minus_poly(entry, poly_coeffs, Tbase, xmin, xmax):
    """Worst relative gap between a Chebyshev entry and a centered polynomial.

    Both describe the same property, so on the fit domain they must agree.
    Sampled on a lattice, which is enough: the two are polynomials of the
    same degree, so they cannot diverge only between the samples.
    """
    from numpy.polynomial import chebyshev as ncheb

    Trange, xbase = entry["Trange"], entry["xbase"]
    M = np.array(entry["coeffs"], dtype=float)
    P = np.atleast_2d(np.array(poly_coeffs, dtype=float))
    Ts = np.linspace(Trange[0], Trange[1], 25)
    u = (2.0 * Ts - (Trange[1] + Trange[0])) / (Trange[1] - Trange[0])

    worst = 0.0
    for x in np.linspace(xmin, xmax, 5):
        a = np.array([sum(M[i, j] * (x - xbase) ** j for j in range(M.shape[1]))
                      for i in range(M.shape[0])])
        cheb = ncheb.chebval(u, a)
        poly = sum(P[i, j] * (Ts - Tbase) ** i * (x - xbase) ** j
                   for i in range(P.shape[0]) for j in range(P.shape[1]))
        worst = max(worst, float(np.max(np.abs(cheb - poly)
                                        / np.maximum(np.abs(poly), 1e-30))))
    return worst


# Fluids used by the writer test below. They must be ones whose caloric
# polynomials are least-squares FITS, because only then does the 7-digit
# rounding actually change the coefficients. The Melinder fluids (MEG and
# friends) carry book coefficients that already fit in 7 digits, so rounding
# them is a no-op and the ordering bug is invisible on them -- an earlier
# version of this test used MEG and passed with the bug reinstated.
WRITER_TEST_FLUIDS = ("TVP1869", "ZMC")


def test_writer_keeps_basis_conversions_exact_against_what_it_serialises(monkeypatch, tmp_path):
    """Drive the writer and check the SERIALISED polynomial, not the file.

    test_chebyshev_entries.py checks the same invariant, but only by reading
    json/. That cannot see a writer that has stopped producing it: reverting
    this fix and keeping the regenerated files leaves the suite green. The
    bug it protects against is an ordering one, so it only appears when the
    writer actually runs.

    A basis conversion is exact with respect to the polynomial it is derived
    from. If the writer derives it from the full-precision polynomial in
    memory and rounds that polynomial afterwards, the two committed
    representations disagree by ~1e-7 even though each is individually fine.
    """
    from CPIncomp import ChebyshevFits
    from add_chebyshev_entries import collect_fluid_objects

    objects = collect_fluid_objects()
    writer = SolutionDataWriter()

    # Fit first: pure numerics, but it reads the SecCool data files by a path
    # relative to this directory.
    targets = []
    for name in WRITER_TEST_FLUIDS:
        target = objects[name]
        writer.fitFluidList([target])
        # toJSON ends by reloading the written file into the object, which
        # replaces these with their rounded selves, so snapshot them now.
        unrounded = {prop: np.atleast_2d(np.array(getattr(target, prop).coeffs,
                                                  dtype=float))
                     for prop in ChebyshevFits.CALORIC_PROPERTIES}
        targets.append((target, unrounded))

    # toJSON writes json/<name>.json and the hash cache relative to the
    # current directory, so running this from the repo root used to drop a
    # stray json/ into the working tree.
    monkeypatch.chdir(tmp_path)

    compared, sensitive = 0, 0
    for target, unrounded in targets:
        captured = {}
        realGetHash = writer.get_hash

        def captureDump(dump, _captured=captured, _real=realGetHash):
            _captured["dump"] = dump
            return _real(dump)

        writer.get_hash = captureDump
        try:
            writer.toJSON(target, quiet=True)
        finally:
            writer.get_hash = realGetHash

        assert "dump" in captured, "the writer did not serialise anything"
        fluid = json.loads(captured["dump"])
        Tbase = float(fluid.get("Tbase", 0.0) or 0.0)

        for prop in ChebyshevFits.CALORIC_PROPERTIES:
            entry = fluid.get(prop + "_cheb")
            committed = fluid.get(prop, {})
            if (entry is None
                    or entry.get("fit_source") not in ChebyshevFits.EXACT_FIT_SOURCES
                    or committed.get("type") != "polynomial"):
                continue

            gap = _cheb_minus_poly(entry, committed["coeffs"], Tbase,
                                   fluid["xmin"], fluid["xmax"])
            assert gap < ChebyshevFits.EXACT_CONVERSION_TOLERANCE, (
                "{0} {1}: conversion disagrees with the committed polynomial "
                "by {2:.2e}".format(fluid["name"], prop, gap))
            compared += 1

            # Non-vacuity, measured rather than assumed: would the WRONG
            # order (convert from the unrounded fit, round the polynomial
            # afterwards) actually be caught here? If not, this property
            # proves nothing and must not be counted.
            wrongOrder = ChebyshevFits.convert_polynomial(
                unrounded[prop], Tbase, entry["xbase"], entry["Trange"])
            if _cheb_minus_poly(dict(entry, coeffs=wrongOrder.tolist()),
                                committed["coeffs"], Tbase,
                                fluid["xmin"], fluid["xmax"]) \
                    >= ChebyshevFits.EXACT_CONVERSION_TOLERANCE:
                sensitive += 1

    assert compared > 0, ("no exact Chebyshev conversion was produced, so this "
                          "test checked nothing; pick fluids that produce one")
    assert sensitive > 0, (
        "none of {0} has a caloric polynomial that the 7-digit rounding "
        "actually changes, so this test would pass with the ordering bug "
        "reinstated; pick fluids whose polynomials are fitted, not "
        "tabulated".format(list(WRITER_TEST_FLUIDS)))
