#!/usr/bin/env python3
"""Generate GERG-2004/GERG-2008 reference values from teqp for Tasks 8-9.

teqp is the reference implementation the CoolProp GERG backend (Task 9) is
validated against; see
.superpowers/sdd/2026-07-25-gerg-strict-backend/task-7-brief.md.  This
script is a DEVELOPER TOOL, not part of CoolProp's build: nothing in
CMakeLists.txt invokes it, and CoolProp's build has no dependency on teqp.
The committed output, src/Backends/GERG/GERGReferenceValues.h, is what
Tasks 8/9 actually compile against.

Run (exact commands used to produce the committed header; teqp 0.23.2):

    python3 -m venv /tmp/gergenv
    /tmp/gergenv/bin/pip install "teqp>=0.18" numpy
    /tmp/gergenv/bin/python dev/gerg/generate_reference_values.py \
        > src/Backends/GERG/GERGReferenceValues.h
    uvx clang-format@18.1.8 -i src/Backends/GERG/GERGReferenceValues.h

The clang-format pass is required, not cosmetic: this script's own emission
(4-space indent, one point per source line) is NOT what
`clang-format --dry-run -Werror` accepts, and the pre-PR gate runs exactly
that check. clang-format re-indents to 2 spaces and, for the AGA8 mixture
rows (21-element names/z vectors per row), explodes each row across ~20
lines (one element per line) to respect the column limit -- this makes
those rows harder to eyeball as single-record rows, but every value is
still complete and correct, and there is no `// clang-format off` precedent
in this repo to suppress it. Regenerating and re-formatting is expected to
reproduce byte-identical output unless teqp's own GERG implementation
changes; see dev/gerg/README.md.

The system python3's teqp (0.15.3 as of this writing) predates the GERG
factory registration and fails with "Unknown kind:GERG2008resid" -- use the
venv above, not the system interpreter.

DERIVATIVE-ACCESSOR NAMING (load-bearing, verify against a fresh teqp
checkout before trusting this comment on a version bump): teqp's Python
accessors are get_Ar<TAU order><DELTA order>, i.e. get_Ar20 is the
tau-tau second derivative, get_Ar02 is delta-delta. An earlier draft of
this generator had these reversed, which still produces a finite plausible
number silently. This was verified empirically on the ideal-gas model:
get_Ar20 gives -cv_ideal/R (3.303 for methane at 300 K, which is right,
since cv_ideal/R -> ~3.3 for a polyatomic like methane well above its
Debye-like vibrational thresholds), while get_Ar02 gives exactly -1.0 (the
ideal gas has no delta-dependence in alpha^o beyond the ln(delta) term, so
its OWN first delta-derivative is 1 and everything past that vanishes --
-1.0 is a tell that the wrong derivative order was read, not a real cv).

The p/cv/w expressions below mirror teqp's own test file exactly
(~/Code/teqp/src/tests/catch_test_GERG.cxx:675-694, the "Validate all
GERG2008 models" test case) rather than being re-derived; any future
disagreement between this script and that test is a bug in this script,
not in teqp.  That test reproduces the published AGA8 gas-2 validation
point (T = 190.68 K, D = 11.0 mol/L) to 3.96e-13 relative on pressure when
(and only when) the mixture is built with the AGA8 component ordering
below, not GERG2008::component_names order (see AGA8_NAMES docstring).
"""
import math
import sys
from pathlib import Path

import numpy as np
import teqp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _mixture_comps import MIXTURE_COMPS  # noqa: E402  (verified byte-identical to
from _validation_data import VALIDATION_DATA  # noqa: E402   teqp's own tables, see dev/gerg/README.md)

R = 8.314472  # J/mol/K -- the GERG gas constant, teqp GERG.hpp:470, confirmed via m.get_R(z).

# GERG2004::component_names / GERG2008::component_names order (teqp GERG.hpp:431,973).
# This is the order src/Backends/GERG/GERGData.h::component_names() uses, and the
# order pure-fluid and same-model binary-pair points below are generated in.
GERG2004_NAMES = [
    "methane", "nitrogen", "carbondioxide", "ethane", "propane", "n-butane",
    "isobutane", "n-pentane", "isopentane", "n-hexane", "n-heptane", "n-octane",
    "hydrogen", "oxygen", "carbonmonoxide", "water", "helium", "argon",
]  # fmt: skip
GERG2008_NAMES = GERG2004_NAMES + ["hydrogensulfide", "n-nonane", "n-decane"]

# CRITICAL: the AGA8 test-code component ordering, DISTINCT from GERG2008_NAMES
# above (isobutane precedes n-butane, isopentane precedes n-pentane, helium and
# argon are last).  mixture_comps/validation_data columns follow THIS order --
# see catch_test_GERG.cxx:512, transcribed verbatim.  Using GERG2008_NAMES to
# build the AGA8 mixture models instead silently misassigns every composition:
# still sums to 1, still evaluates, wrong by a few tenths of a percent (verified:
# gas 2 gives p=4.6142620561 MPa with the wrong order vs 4.62270367011 MPa,
# matching the published value to 4e-13, with this one).
AGA8_NAMES = [
    "methane", "nitrogen", "carbondioxide", "ethane", "propane", "isobutane", "n-butane",
    "isopentane", "n-pentane", "n-hexane", "n-heptane", "n-octane", "n-nonane", "n-decane",
    "hydrogen", "oxygen", "carbonmonoxide", "water", "hydrogensulfide", "helium", "argon",
]  # fmt: skip

# Reducing parameters and molar masses, transcribed in Task 2 into
# src/Backends/GERG/GERGData.h (pure_info_2004()/pure_info_2008_overrides()) from
# the same GERG-2004 monograph Table A3.5 / GERG-2008 Table A2 teqp itself uses.
# This is a SECOND on-disk copy of the same numbers, not an independent source --
# deliberately so: see the get_Tcvec/get_vcvec cross-check in pure_points() below,
# which verifies this copy against teqp's OWN internal reducing parameters at
# generation time rather than trusting the transcription silently.
# Units as tabulated: Tc_K, rhoc in mol/dm^3, M in kg/kmol.
PURE_INFO_2004 = {
    "methane": (190.564000000, 10.139342719, 16.042460),
    "nitrogen": (126.192000000, 11.183900000, 28.013400),
    "carbondioxide": (304.128200000, 10.624978698, 44.009500),
    "ethane": (305.322000000, 6.870854540, 30.069040),
    "propane": (369.825000000, 5.000043088, 44.095620),
    "n-butane": (425.125000000, 3.920016792, 58.122200),
    "isobutane": (407.817000000, 3.860142940, 58.122200),
    "n-pentane": (469.700000000, 3.215577588, 72.148780),
    "isopentane": (460.350000000, 3.271018581, 72.148780),
    "n-hexane": (507.820000000, 2.705877875, 86.175360),
    "n-heptane": (540.130000000, 2.315324434, 100.201940),
    "n-octane": (569.320000000, 2.056404127, 114.228520),
    "hydrogen": (33.190000000, 14.940000000, 2.015880),
    "oxygen": (154.595000000, 13.630000000, 31.998800),
    "carbonmonoxide": (132.800000000, 10.850000000, 28.010100),
    "water": (647.096000000, 17.873716090, 18.015280),
    "helium": (5.195300000, 17.399000000, 4.002602),
    "argon": (150.687, 13.407429659, 39.948000),
}
PURE_INFO_2008_OVERRIDES = {
    "carbonmonoxide": (132.86, 10.85, 28.010100),
    "isopentane": (460.35, 3.271, 72.148780),
    "n-nonane": (594.55, 1.81, 128.2551),
    "n-decane": (617.7, 1.64, 142.28168),
    "hydrogensulfide": (373.1, 10.19, 34.08088),
}


def pure_info_si(year, name):
    """(Tc_K, rhoc_molm3, M_kgmol) for one component under the given model year."""
    Tc, rhoc, M = PURE_INFO_2008_OVERRIDES.get(name) if year == 2008 and name in PURE_INFO_2008_OVERRIDES else PURE_INFO_2004[name]
    return Tc, rhoc * 1000.0, M / 1000.0


def molar_mass(name):
    """kg/mol for any of the 21 GERG-2008 names; molar mass never differs
    between the GERG-2004 and GERG-2008 tables for the 18 shared components."""
    return pure_info_si(2008, name)[2]


def models(year, names):
    resid = teqp.make_model({"kind": f"GERG{year}resid", "model": {"names": names}})
    ideal = teqp.make_model({"kind": f"GERG{year}idealgas", "model": {"names": names}})
    return resid, ideal


def _safe(fn):
    """Evaluate a zero-arg thunk, returning NaN (never raising) if it either
    raises or returns a non-finite value.  Used so that ONE quantity failing
    at a given (T, rho) -- e.g. w blowing up on the mechanically-unstable
    branch -- never prevents the OTHER quantities computed at the exact same
    state from being reported."""
    try:
        v = fn()
    except Exception:
        return float("nan")
    return v if math.isfinite(v) else float("nan")


def point(resid, ideal, z, T, rho, M):
    """(alphar, alphaig, p_Pa, cv_JmolK, w_ms) at (T, rho) for composition z
    and mixture molar mass M (kg/mol).  Mirrors teqp's own test file exactly
    (catch_test_GERG.cxx:675-694).

    Each of the 5 returned values is computed and null'd to NaN
    INDEPENDENTLY -- never by discarding the whole point.  w legitimately
    goes NaN on the mechanically-unstable branch (w^2 < 0) while p and cv
    stay perfectly finite at the exact same (T, rho); the pure-fluid grid
    used to drop the WHOLE row whenever any single one of the 5 values was
    non-finite, which silently threw away ~50 perfectly good p/cv fixture
    points clustered right at the near-phase-boundary states where a wiring
    bug is most likely to show up (see dev/gerg/README.md). Callers must
    NOT drop a row because one field is NaN; keep the row, null just that
    field, and let downstream consumers skip per-field."""
    z = np.asarray(z, dtype=float)
    Ar00 = _safe(lambda: resid.get_Ar00(T, rho, z))
    Ar01 = _safe(lambda: resid.get_Ar01(T, rho, z))  # delta * dalphar/ddelta
    Ar02 = _safe(lambda: resid.get_Ar02(T, rho, z))  # delta^2 * d2alphar/ddelta^2
    Ar11 = _safe(lambda: resid.get_Ar11(T, rho, z))  # delta*tau * d2alphar/(ddelta dtau)
    Ar20 = _safe(lambda: resid.get_Ar20(T, rho, z))  # tau^2 * d2alphar/dtau^2
    Aig00 = _safe(lambda: ideal.get_Ar00(T, rho, z))
    Aig20 = _safe(lambda: ideal.get_Ar20(T, rho, z))

    p = rho * R * T * (1.0 + Ar01) if math.isfinite(Ar01) else float("nan")
    cv = -(Aig20 + Ar20) * R if all_finite(Aig20, Ar20) else float("nan")
    w = float("nan")
    if all_finite(Ar01, Ar02, Ar11, Aig20, Ar20):
        with np.errstate(invalid="ignore", divide="ignore"):
            w2 = (R * T / M) * (1 + 2 * Ar01 + Ar02 - (1 + Ar01 - Ar11) ** 2 / (Aig20 + Ar20))
        if math.isfinite(w2) and w2 > 0:
            w = math.sqrt(w2)
    return Ar00, Aig00, p, cv, w


def all_finite(*vals):
    return all(math.isfinite(v) for v in vals)


# ---------------------------------------------------------------------------
# Pure-fluid grid
# ---------------------------------------------------------------------------

# Every fluid is expected to retain a finite alphar at (nearly) all 16 grid
# points -- alphar (teqp's get_Ar00) is the most basic quantity computed (a
# single function evaluation, no ratio/subtraction that can blow up the way
# w's does near the mechanically-unstable branch) and should only go
# non-finite when teqp's underlying solver itself gives up, not from benign
# near-phase-boundary cancellation. This is a PER-FLUID floor, not a global
# total, precisely so one fluid failing broadly cannot hide behind other
# fluids that happen to have full coverage (see dev/gerg/README.md).
MIN_FINITE_ALPHAR_PER_FLUID = 14  # of 16; see README for the empirical basis


def pure_points(year, names):
    rows = []
    # Per-fluid, per-field counts of NaN'd (non-finite) values, for the
    # coverage summary -- NOT for dropping rows; see point()'s docstring.
    field_nulls = {name: {"alphar": 0, "alphaig": 0, "p": 0, "cv": 0, "w": 0} for name in names}
    for name in names:
        Tc, rhoc, M = pure_info_si(year, name)
        resid, ideal = models(year, [name])

        # Cross-check the transcribed reducing parameters against teqp's OWN
        # internal values for this exact model, rather than trusting the
        # second on-disk copy above silently.
        Tc_teqp = resid.get_Tcvec()[0]
        rhoc_teqp = 1.0 / resid.get_vcvec()[0]
        assert abs(Tc_teqp - Tc) <= 1e-9 * Tc, f"{year}/{name}: Tc mismatch {Tc_teqp} vs {Tc}"
        assert abs(rhoc_teqp - rhoc) <= 1e-9 * rhoc, f"{year}/{name}: rhoc mismatch {rhoc_teqp} vs {rhoc}"

        n_alphar_finite = 0
        for Tfac in (0.7, 1.0, 1.3, 2.0):
            for rfac in (0.01, 0.5, 1.5, 2.5):
                T = Tfac * Tc
                rho = rfac * rhoc
                Ar00, Aig00, p, cv, w = point(resid, ideal, [1.0], T, rho, M)
                for label, val in (("alphar", Ar00), ("alphaig", Aig00), ("p", p), ("cv", cv), ("w", w)):
                    if not math.isfinite(val):
                        field_nulls[name][label] += 1
                if math.isfinite(Ar00):
                    n_alphar_finite += 1
                # ALWAYS keep the row -- individual fields may be NaN (see
                # point()'s docstring); do not drop it because one field
                # failed while the fixed grid still deterministically
                # produces 16 rows/fluid.
                rows.append((name, T, rho, Ar00, Aig00, p, cv, w))

        assert n_alphar_finite >= MIN_FINITE_ALPHAR_PER_FLUID, (
            f"{year}/{name}: only {n_alphar_finite}/16 grid points have a finite alphar -- "
            "this fluid is failing broadly (not just isolated near-phase-boundary points); "
            "investigate before shipping, do not silently ship a near-empty fluid."
        )
    return rows, field_nulls


# ---------------------------------------------------------------------------
# Same-model binary pairs: every (f_i, f_j), i<j, at one mid-range state.
# This is what makes Task 9 exercise every reducing-parameter row (153 for
# GERG-2004, 210 for GERG-2008) numerically -- see task-7-brief.md.
# ---------------------------------------------------------------------------


def binary_pair_points(year, names):
    T, rho, z = 250.0, 5000.0, [0.4, 0.6]
    rows = []
    skipped = []
    nan_w = []
    nan_alpha = []  # rows whose alphar or alphaig came back non-finite; see below
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            f1, f2 = names[i], names[j]
            resid, ideal = models(year, [f1, f2])
            M = z[0] * molar_mass(f1) + z[1] * molar_mass(f2)
            try:
                Ar00, Aig00, p, cv, w = point(resid, ideal, z, T, rho, M)
            except Exception:
                skipped.append((f1, f2))
                continue
            # This fixed (T=250K, rho=5000 mol/m^3) state falls inside the
            # mechanically-unstable (spinodal) branch of the single-phase EOS
            # for some heavy/light or highly non-ideal pairs -- w^2 legitimately
            # goes negative there (see dev/gerg/README.md for the CO2/ethane
            # example and the physical argument). p and cv are still analytic
            # and finite on that branch, so KEEP the row (this is the guard the
            # brief requires for every reducing-parameter row) and only NaN out
            # w, rather than dropping the whole pair.
            if not all_finite(p, cv):
                skipped.append((f1, f2))
                continue
            if not math.isfinite(w):
                nan_w.append((f1, f2))
            # alphar/alphaig are the two halves of the EOS in isolation, and
            # alphaig is the ONLY column that can see an ideal-gas
            # integration-constant or Tc-vs-T_red error (Task 8 measured this:
            # a 1e-9 shift in n0[1] failed 624 alphaig assertions with alphar,
            # p, cv and w all clean).  A row silently arriving here with a NaN
            # in either column would be skipped per-field downstream and would
            # remove that guard without failing anything, so it is counted and
            # asserted on in main() rather than tolerated.
            if not all_finite(Ar00, Aig00):
                nan_alpha.append((f1, f2))
            rows.append(([f1, f2], list(z), T, rho, Ar00, Aig00, p, cv, w))
    return rows, skipped, nan_w, nan_alpha


# ---------------------------------------------------------------------------
# AGA8 gas mixtures (GERG-2008 only: these compositions routinely carry
# nonzero hydrogen sulfide / n-nonane / n-decane, which GERG-2004 has no
# coefficients for at all).
# ---------------------------------------------------------------------------


def aga8_points():
    resid, ideal = models(2008, AGA8_NAMES)
    mw = np.array([molar_mass(n) for n in AGA8_NAMES])
    rows = []
    skipped = []
    nan_w = []
    nan_alpha = []  # same rationale as binary_pair_points'
    p_reldiffs = []  # (reldiff, gasno) for every emitted row -- used to report both
    # the worst-case and the median disagreement against validation_data.P_MPa.
    for gasno, T, D_molL, P_MPa, cv_ref, cp_ref, w_ref in VALIDATION_DATA:
        z = np.asarray(MIXTURE_COMPS[gasno - 2], dtype=float) / 100.0
        assert abs(z.sum() - 1.0) < 1e-6, f"gas {gasno}: mole fractions sum to {z.sum()}"
        rho = D_molL * 1000.0
        M = float(np.dot(z, mw))
        try:
            Ar00, Aig00, p, cv, w = point(resid, ideal, z, T, rho, M)
        except Exception:
            skipped.append(gasno)
            continue
        # Same rationale as binary_pair_points: keep p/cv even if this
        # particular (T, rho, z) happens to sit on the mechanically-unstable
        # branch for this composition (w^2 < 0); none of the 187 AGA8 rows
        # actually hit that in practice, but the guard is symmetric with
        # binary_pair_points for robustness against future re-generation.
        if not all_finite(p, cv):
            skipped.append(gasno)
            continue
        if not math.isfinite(w):
            nan_w.append(gasno)
        if not all_finite(Ar00, Aig00):
            nan_alpha.append(gasno)
        p_reldiffs.append((abs(p / 1e6 - P_MPa) / P_MPa, gasno))
        # Emit all 21 AGA8_NAMES/mole-fractions (not a trimmed nonzero subset) so
        # component order is unambiguous -- see the module docstring on AGA8_NAMES.
        rows.append((list(AGA8_NAMES), list(z), T, rho, Ar00, Aig00, p, cv, w))
    return rows, skipped, nan_w, nan_alpha, p_reldiffs


# ---------------------------------------------------------------------------
# C++ emission
# ---------------------------------------------------------------------------


def fnum(x):
    """Emit `x` as a C++ literal/expression, full round-trip precision via
    %.17g for finite values. A handful of mix_points_2008/mix_points_2004
    rows carry a legitimately NaN `w` (see binary_pair_points): "nan" is not
    a valid C++ floating-point literal, so those need the <limits> spelling
    instead -- plain %.17g would silently emit invalid, non-compiling source."""
    if math.isnan(x):
        return "std::numeric_limits<double>::quiet_NaN()"
    if math.isinf(x):
        return ("-" if x < 0 else "") + "std::numeric_limits<double>::infinity()"
    return "%.17g" % x


def emit_pure(rows, varname, out):
    out.append(f"inline const std::vector<PureRefPoint> {varname} = {{")
    for name, T, rho, alphar, alphaig, p, cv, w in rows:
        out.append(
            f'    {{"{name}", {fnum(T)}, {fnum(rho)}, {fnum(alphar)}, {fnum(alphaig)}, {fnum(p)}, {fnum(cv)}, {fnum(w)}}},'
        )
    out.append("};")
    out.append("")


def emit_mix(rows, varname, out):
    out.append(f"inline const std::vector<MixRefPoint> {varname} = {{")
    for names, z, T, rho, alphar, alphaig, p, cv, w in rows:
        names_str = ", ".join(f'"{n}"' for n in names)
        z_str = ", ".join(fnum(v) for v in z)
        out.append(
            f"    {{{{{names_str}}}, {{{z_str}}}, {fnum(T)}, {fnum(rho)}, "
            f"{fnum(alphar)}, {fnum(alphaig)}, {fnum(p)}, {fnum(cv)}, {fnum(w)}}},"
        )
    out.append("};")
    out.append("")


def main():
    pure_2004, nulls_2004 = pure_points(2004, GERG2004_NAMES)
    pure_2008, nulls_2008 = pure_points(2008, GERG2008_NAMES)
    bin_2004, binskip_2004, nanw_2004, nanalpha_2004 = binary_pair_points(2004, GERG2004_NAMES)
    bin_2008, binskip_2008, nanw_2008, nanalpha_2008 = binary_pair_points(2008, GERG2008_NAMES)
    aga8_rows, aga8_skip, aga8_nanw, aga8_nanalpha, p_reldiffs = aga8_points()

    mix_2008 = bin_2008 + aga8_rows

    # --- Coverage assertions: fail loudly instead of silently shipping a
    # header with a missing reducing-parameter row or AGA8 gas. Task 5's
    # review established that a single mistyped digit in one of the 225
    # reducing-parameter rows (153 GERG-2004 + 72 GERG-2008-override) is
    # guarded by NOTHING except these binary-pair points existing and being
    # numerically checked in Task 8/9 -- a pair silently vanishing from the
    # header (e.g. a future teqp bump making it evaluate non-finite instead
    # of finite-but-wrong) would remove that row's only guard while leaving
    # every test green. See dev/gerg/README.md.
    n_pairs_2004 = len(GERG2004_NAMES) * (len(GERG2004_NAMES) - 1) // 2
    n_pairs_2008 = len(GERG2008_NAMES) * (len(GERG2008_NAMES) - 1) // 2
    assert n_pairs_2004 == 153, f"GERG2004_NAMES has {len(GERG2004_NAMES)} names -> C(n,2)={n_pairs_2004}, expected 153"
    assert n_pairs_2008 == 210, f"GERG2008_NAMES has {len(GERG2008_NAMES)} names -> C(n,2)={n_pairs_2008}, expected 210"
    assert len(bin_2004) == n_pairs_2004, (
        f"mix_points_2004 has {len(bin_2004)} binary-pair rows, expected exactly {n_pairs_2004} "
        f"-- {len(binskip_2004)} pair(s) dropped: {binskip_2004}. A missing pair removes the ONLY "
        "numeric guard on that reducing-parameter row; investigate, do not ship."
    )
    assert len(bin_2008) == n_pairs_2008, (
        f"mix_points_2008 binary pairs: {len(bin_2008)} rows, expected exactly {n_pairs_2008} "
        f"-- {len(binskip_2008)} pair(s) dropped: {binskip_2008}"
    )
    assert len(aga8_rows) == len(VALIDATION_DATA) == 187, (
        f"AGA8 rows: {len(aga8_rows)} emitted from {len(VALIDATION_DATA)} validation_data entries, "
        f"expected 187 of both -- {len(aga8_skip)} gas(es) dropped: {aga8_skip}"
    )
    # Every emitted mixture row must carry a finite alphar AND a finite alphaig.
    # Task 9's gate skips a NaN reference field rather than failing on it (a
    # NaN w is legitimate on the spinodal branch), so a NaN that crept into the
    # alpha columns would quietly delete the only instrument that can see an
    # ideal-gas error -- exactly the failure mode these columns were added for.
    # Assert here, at generation time, instead of shipping it.
    assert not nanalpha_2004, f"mix_points_2004: non-finite alphar/alphaig on {len(nanalpha_2004)} pair(s): {nanalpha_2004}"
    assert not nanalpha_2008, f"mix_points_2008 binary pairs: non-finite alphar/alphaig on {len(nanalpha_2008)} pair(s): {nanalpha_2008}"
    assert not aga8_nanalpha, f"mix_points_2008 AGA8: non-finite alphar/alphaig on {len(aga8_nanalpha)} gas(es): {aga8_nanalpha}"

    reldiffs_sorted = sorted(p_reldiffs)
    worst_reldiff, worst_gasno = reldiffs_sorted[-1]
    median_reldiff = reldiffs_sorted[len(reldiffs_sorted) // 2][0]
    n_over_1e6 = sum(1 for d, _ in p_reldiffs if d > 1e-6)

    # --- AGA8 component-ordering gate. ASSERTED, not merely reported.
    #
    # This median is the discriminator for the ordering bug the module
    # docstring calls CRITICAL: MIXTURE_COMPS is indexed in AGA8_NAMES order,
    # and pairing it with GERG2008_NAMES instead silently misassigns every
    # composition -- it does not raise, it does not produce NaN, it produces
    # a full, plausible-looking 187-row fixture of wrong numbers.  Measured
    # both ways in this tree (see dev/gerg/README.md):
    #
    #   correct AGA8_NAMES ordering  : median 3.5e-12, worst 6.8e-04
    #   wrong  GERG2008_NAMES ordering: median 7.1e-04, worst 2.7e-01
    #
    # 1e-9 sits ~285x above the correct-ordering median and ~7e5x below the
    # wrong-ordering median, so it discriminates with enormous margin while
    # leaving room for ordinary teqp-version / libm drift in the last digits.
    # Deliberately the MEDIAN and not the worst case: the worst rows are a
    # known table-vs-reference-code provenance issue on a handful of rich
    # CO2/H2S/He blends (documented in the emitted header below and in
    # README.md), which is why worst_reldiff and n_over_1e6 stay reported
    # rather than asserted -- but nothing in that provenance issue can move
    # the MEDIAN off ~1e-12.
    #
    # Before this assert existed, a regeneration that reintroduced the
    # ordering bug emitted the full 187-row fixture and exited 0.
    assert median_reldiff < 1e-9, (
        f"AGA8 p_teqp vs validation_data.P_MPa: median relative difference is {median_reldiff:.3e}, "
        "expected < 1e-9. The overwhelmingly likely cause is a component-ORDERING bug -- "
        "MIXTURE_COMPS rows are ordered per AGA8_NAMES and must not be paired with GERG2008_NAMES "
        "(that mistake gives a median of ~7.1e-4). Do NOT relax this threshold; fix the ordering."
    )

    def null_summary(nulls, names):
        totals = {"alphar": 0, "alphaig": 0, "p": 0, "cv": 0, "w": 0}
        for name in names:
            for k in totals:
                totals[k] += nulls[name][k]
        return ", ".join(f"{k}={v}" for k, v in totals.items())

    report = [
        f"# pure_points_2004: {len(pure_2004)} emitted ({len(GERG2004_NAMES)} fluids x 16 grid pts, "
        f"floor {MIN_FINITE_ALPHAR_PER_FLUID}/16 finite alphar enforced per fluid), NaN'd-field counts: "
        f"{null_summary(nulls_2004, GERG2004_NAMES)}",
        f"# pure_points_2008: {len(pure_2008)} emitted ({len(GERG2008_NAMES)} fluids x 16 grid pts, "
        f"floor {MIN_FINITE_ALPHAR_PER_FLUID}/16 finite alphar enforced per fluid), NaN'd-field counts: "
        f"{null_summary(nulls_2008, GERG2008_NAMES)}",
        "# No pure-fluid row is ever dropped: every field (alphar/alphaig/p/cv/w) is nulled "
        "independently when non-finite, per point()'s docstring -- Tasks 8/9 must skip per-FIELD "
        "(std::isnan check on that field alone), not per-row.",
        f"# mix_points_2004 (binary pairs only): {len(bin_2004)} of {n_pairs_2004} expected emitted "
        f"(assert enforced), {len(nanw_2004)} with w = NaN (mechanically-unstable branch; p, cv still valid)",
        f"# mix_points_2008 (binary pairs + AGA8): {len(mix_2008)} emitted "
        f"({len(bin_2008)} of {n_pairs_2008} expected pairs + {len(aga8_rows)} of 187 expected AGA8 gases, "
        f"both counts assert-enforced), {len(nanw_2008)} pairs with w = NaN, {len(aga8_nanw)} AGA8 with w = NaN",
        "# Mixture alphar/alphaig: 0 NaN on every emitted row of both vectors (assert-enforced). w is the "
        "ONLY mixture column that is ever NaN, so Task 9's per-field skip may fire on w and on nothing else.",
        f"# AGA8 p_teqp vs validation_data.P_MPa: worst {worst_reldiff:.3e} (gas {worst_gasno}), "
        f"median {median_reldiff:.3e} across {len(p_reldiffs)} rows, {n_over_1e6} rows > 1e-6 -- see",
        "# dev/gerg/README.md for why this is a table-vs-reference-code provenance issue (confirmed against",
        "# NIST's own bundled AGA8 reference code, not just teqp), not a composition/ordering bug: the median",
        "# matches the ~1e-12 verified-correct-ordering baseline, teqp's OWN test disables this exact",
        "# comparison (catch_test_GERG.cxx:696-698), and the worst rows are all rich, complex multi-component",
        "# natural-gas blends (CO2/H2S/He-rich), not the uniform ~2e-3 offset a swapped-column bug would produce.",
    ]
    print("\n".join(f"// {line[2:]}" for line in report), file=sys.stderr)

    out = []
    out.append("// GENERATED FILE -- do not edit by hand.")
    out.append("//")
    out.append("// Produced by dev/gerg/generate_reference_values.py from teqp 0.23.2")
    out.append("// (https://github.com/usnistgov/teqp), the reference implementation the")
    out.append("// CoolProp GERG-2004/GERG-2008 backend is validated against. Regenerating")
    out.append("// this file is expected to produce a byte-identical result unless teqp's")
    out.append("// own GERG implementation changes. See dev/gerg/README.md for the exact")
    out.append("// venv commands and why this header is committed rather than generated")
    out.append("// at CoolProp build time (CoolProp's build has no dependency on teqp).")
    out.append("//")
    for line in report:
        # `report` entries are plain text prefixed with "# " (Python-comment
        # style, used verbatim for the human-readable stderr summary above).
        # Re-prefix with "//" for the header itself so these are valid C++
        # comments -- a bare line starting with '#' is a preprocessor
        # directive, and an earlier version of this script emitted `report`
        # here unmodified, which produced "invalid preprocessing directive"
        # compile errors for every one of these lines.
        out.append(f"// {line[2:]}")
    out.append("#ifndef GERGREFERENCEVALUES_H_")
    out.append("#define GERGREFERENCEVALUES_H_")
    out.append("")
    out.append("#include <limits>")
    out.append("#include <vector>")
    out.append("")
    out.append("namespace CoolProp {")
    out.append("namespace GERG {")
    out.append("namespace reference {")
    out.append("")
    out.append("/// One (T, rho) evaluation of a single pure-fluid GERG residual/ideal-gas")
    out.append("/// EOS, generated directly from teqp -- NOT from CoolProp's own tables, so")
    out.append("/// comparing against these values tests the ASSEMBLED equation of state,")
    out.append("/// not just the coefficient tables Tasks 1-6 already cross-checked.")
    out.append("/// `name` is the GERG component name (matches GERGData.h component_names()).")
    out.append("/// `alphar`/`alphaig` are the dimensionless residual/ideal-gas Helmholtz")
    out.append("/// energies (teqp's get_Ar00 on the resid/ideal-gas models respectively).")
    out.append("///")
    out.append("/// Any of alphar/alphaig/p_Pa/cvmolar/w may independently be NaN")
    out.append("/// (std::isnan(...) true) at a given row -- e.g. w legitimately goes NaN on")
    out.append("/// the mechanically-unstable branch while p_Pa/cvmolar stay perfectly finite")
    out.append("/// at the exact same (T, rho). Rows are NEVER dropped for having one")
    out.append("/// non-finite field: Task 8/9 must check std::isnan(...) per FIELD being")
    out.append("/// compared, not skip the whole row when any single field is NaN.")
    out.append("struct PureRefPoint")
    out.append("{")
    out.append("    const char* name;")
    out.append("    double T_K, rhomolar, alphar, alphaig, p_Pa, cvmolar, w;")
    out.append("};")
    out.append("")
    out.append("/// One (T, rho, z) evaluation of a mixture GERG EOS. `names`/`z` are")
    out.append("/// parallel vectors of equal length so the composition is self-describing")
    out.append("/// -- component ORDER matters for GERG reducing parameters (see")
    out.append("/// GERGData.h get_betasgammas' orientation note), and mixing up the order")
    out.append("/// of `names` and `z` independently in Task 9 must not be possible to do")
    out.append("/// silently the way a bare std::vector<double> composition would allow.")
    out.append("///")
    out.append("/// `alphar`/`alphaig` are the dimensionless residual/ideal-gas Helmholtz")
    out.append("/// energies of the MIXTURE (teqp's get_Ar00 on the GERG200Xresid and")
    out.append("/// GERG200Xidealgas models respectively), so `alphaig` includes teqp's own")
    out.append("/// sum_i x_i*(alpha0_i + ln x_i) mixing term.  `alphaig` is the only column")
    out.append("/// in this struct that can detect an ideal-gas error: the mixture branch of")
    out.append("/// CoolProp's calc_alpha0_deriv_nocache evaluates each component's alpha0 at")
    out.append("/// the component's own Tc_i/T, and getting the Tc-vs-T_red convention wrong")
    out.append("/// there perturbs alphaig long before it is visible in p (which is blind to")
    out.append("/// the ideal-gas part entirely).  Both are ALWAYS finite -- the generator")
    out.append("/// asserts it -- unlike `w` below.")
    out.append("///")
    out.append("/// `w` is NaN (std::isnan(w) true) for a small number of rows in")
    out.append("/// mix_points_2004/mix_points_2008: this fixed (T=250K, rho=5000 mol/m^3)")
    out.append("/// state falls on the mechanically-unstable (spinodal, d p/d rho|T < 0)")
    out.append("/// branch of the single-phase EOS for some pairs of dissimilar components")
    out.append("/// (e.g. carbondioxide/ethane), where w^2 is negative. `p_Pa` and `cvmolar`")
    out.append("/// remain finite and meaningful there (alphar and its low-order T,rho")
    out.append("/// derivatives are still analytic on that branch) and are NOT dropped --")
    out.append("/// Task 9 must compare p_Pa/cvmolar unconditionally but SKIP the w")
    out.append("/// comparison whenever std::isnan(w) is true for that row.")
    out.append("struct MixRefPoint")
    out.append("{")
    out.append("    std::vector<const char*> names;")
    out.append("    std::vector<double> z;")
    out.append("    double T_K, rhomolar, alphar, alphaig, p_Pa, cvmolar, w;")
    out.append("};")
    out.append("")
    out.append("/// Pure fluids: each of the 18 GERG-2004 / 21 GERG-2008 component EOS at")
    out.append("/// T in {0.7, 1.0, 1.3, 2.0} x Tc crossed with rho in {0.01, 0.5, 1.5, 2.5}")
    out.append("/// x rhoc -- EXACTLY 16 points per fluid, always: no row is ever dropped for")
    out.append("/// a non-finite value (see PureRefPoint's doc comment above and the")
    out.append("/// generator's per-field NaN-count report). Together these two vectors cover")
    out.append("/// all 23 distinct pure EOS (see dev/gerg/verify_transcription.py's Family-1")
    out.append("/// check): the 21 GERG-2008 component names include isopentane and")
    out.append("/// carbonmonoxide with DIFFERENT coefficients than their GERG-2004 entries.")
    emit_pure(pure_2004, "pure_points_2004", out)
    emit_pure(pure_2008, "pure_points_2008", out)
    out.append("/// Mixtures, GERG-2004: every one of the C(18,2) = 153 binary pairs among")
    out.append("/// GERG-2004's own component set, at T = 250 K, rho = 5000 mol/m^3,")
    out.append("/// z = [0.4, 0.6] -- one numeric point per reducing-parameter row so a")
    out.append("/// single mistyped digit anywhere in that table cannot pass silently.")
    out.append("/// GERG-2004 has no coefficients for hydrogen sulfide/n-nonane/n-decane,")
    out.append("/// so it carries no AGA8 gas-mixture points (see mix_points_2008).")
    emit_mix(bin_2004, "mix_points_2004", out)
    out.append("/// Mixtures, GERG-2008: every one of the C(21,2) = 210 binary pairs among")
    out.append("/// GERG-2008's component set (same rationale as mix_points_2004), PLUS")
    out.append("/// every AGA8/GERG-2008 published natural-gas composition and state point")
    out.append("/// from teqp's own validation_data/mixture_comps tables (dev/gerg/_mixture_comps.py,")
    out.append("/// dev/gerg/_validation_data.py -- verified byte-identical to teqp's C++")
    out.append("/// source at generation time). Every AGA8 row carries all 21 component")
    out.append("/// names/mole fractions (many exactly 0.0) rather than a trimmed subset,")
    out.append("/// so Task 9 cannot reintroduce the AGA8-vs-component_names ordering bug")
    out.append("/// (see AGA8_NAMES above) by re-deriving which components are 'present'.")
    emit_mix(mix_2008, "mix_points_2008", out)
    out.append("}  // namespace reference")
    out.append("}  // namespace GERG")
    out.append("}  // namespace CoolProp")
    out.append("#endif /* GERGREFERENCEVALUES_H_ */")

    print("\n".join(out))


if __name__ == "__main__":
    main()
