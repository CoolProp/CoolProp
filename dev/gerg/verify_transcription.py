#!/usr/bin/env python3
"""Mechanical re-verification of the GERG-2004/GERG-2008 coefficient tables.

This is a DEVELOPER TOOL, not part of CoolProp's build or CI. Nothing in
CMakeLists.txt or CTest invokes it. Run it by hand after touching any table
in src/Backends/GERG/GERGData.h or src/Backends/GERG/GERGBackend.cpp:

    python3 dev/gerg/verify_transcription.py

It requires a local checkout of teqp (https://github.com/usnistgov/teqp) --
the reference implementation these tables were transcribed from -- at
~/Code/teqp by default; override with --teqp-root if yours lives elsewhere.

WHY THIS EXISTS: Tasks 3-5 of the GERG backend each did an equivalent
mechanical check (parse teqp's GERG.hpp with a throwaway script, dump the
same tables from a temporary Catch2 [.gergdump] test case, diff the two) but
the scaffold tests and scripts were deleted before committing, per the
project's "no scaffold left behind" convention. That means none of those
checks can be re-run today, and any future edit to these tables (a typo fix,
a new fluid, a value correction) has no mechanical check to fall back on.
This script re-implements all four of those checks -- ONE script, covering
ALL FOUR coefficient-table families -- and commits it so it stays runnable.
A fifth family, pure_info (Table A3.5's tabulated Tc/rhoc/M), was added
later: it was checked once by hand (against teqp GERG.hpp:431-465/973-997)
but that check only ever existed in a review transcript, so it is re-done
here mechanically like the other four.

It re-extracts each table family directly from teqp's GERG.hpp with regexes
(not by building/linking against teqp, which would pull in Eigen/Boost/
autodiff as build dependencies for a lookup-table sanity check) and diffs
every extracted value against the equivalent regex-extracted values from
this repository's own GERGData.h/GERGBackend.cpp. Where the two sources
structure the same data differently (e.g. teqp inlines the pure-fluid
exponent sets into each if/else branch; CoolProp factors them out into
main12_exponents()/mne24_exponents() so multiple fluids share one literal)
the script resolves both down to the same canonical (model, fluid) ->
(n, t, d, c, l) shape before comparing, mirroring each source's own
fallthrough logic (GERG-2008 overrides checked first, else GERG-2004).

Families checked:
  1. Pure-fluid residual coefficients (Table A3.2) -- 23 distinct pure EOS.
  2. Ideal-gas coefficients (Table A3.1) -- 39 (model, component) pairs.
  3. Binary reducing parameters (Table A3.8 / Table A8) -- 153 GERG-2004
     rows, 72 GERG-2008 override/addition rows.
  4. Departure functions and F_ij (Table A3.6) -- 15 pairs with a departure
     function, 7 of them scaled by F_ij != 1.
  5. Pure-fluid reducing/critical parameters (Table A3.5) -- 18 GERG-2004
     rows, 5 GERG-2008 override/addition rows. Compares the raw pre-conversion
     table literals (mol/dm^3, K, kg/kmol) in both teqp's data_map and
     CoolProp's pure_info_2004()/pure_info_2008_overrides(), by parsed float
     value AND by normalised raw literal token (a tolerance-based float
     compare alone can miss a transcription error that a raw-token compare
     catches).

A SIXTH family, the acentric factors in GERGAcentric.h, was added with the
fix for bd CoolProp-deut. It is different in kind from the other five: GERG
publishes no acentric factor, so those values are DERIVED from the GERG pure
equations of state (omega = -1 - log10(p_sat(0.7 T_c)/p_c)) rather than
transcribed from anything, and there is no table in teqp's GERG.hpp to diff
them against. Verifying them means RECOMPUTING them, which is the one check
here that needs teqp's python MODULE at run time rather than just a source
checkout. A missing teqp module is therefore a FAILURE, not a silent skip;
pass --skip-acentric to verify only the five transcribed families, and the
summary line says so.

Exits 0 and prints "OK" for every family on success; exits 1 and prints a
diff-shaped mismatch report on any discrepancy.
"""
import argparse
import os
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
GERG_DATA_H = REPO_ROOT / "src" / "Backends" / "GERG" / "GERGData.h"
GERG_BACKEND_CPP = REPO_ROOT / "src" / "Backends" / "GERG" / "GERGBackend.cpp"
GERG_ACENTRIC_H = REPO_ROOT / "src" / "Backends" / "GERG" / "GERGAcentric.h"

FLOAT_RE = r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?"


def fail(msg):
    print(f"FAIL: {msg}", file=sys.stderr)
    sys.exit(1)


def parse_floats(blob):
    """Parse a comma-separated blob of C++ numeric literals into floats."""
    vals = [v.strip() for v in blob.split(",")]
    vals = [v for v in vals if v != ""]
    out = []
    for v in vals:
        try:
            out.append(float(v))
        except ValueError:
            fail(f"could not parse numeric literal {v!r} from blob {blob!r}")
    return out


def close(a, b, tol=1e-13):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def vec_close(a, b, tol=1e-13):
    return len(a) == len(b) and all(close(x, y, tol) for x, y in zip(a, b))


def require(cond, msg):
    if not cond:
        fail(msg)


def require_new_key(rows, key, label):
    """Refuse a duplicate key while building a parsed table.

    THE WHOLE POINT OF THIS SCRIPT depends on this guard being present in
    every parser.  C++'s `std::map` initializer-list constructor keeps the
    FIRST of two entries with equal keys; a Python `dict` built by
    successive assignment keeps the LAST.  So a table that accidentally
    carries a row twice -- a bad first copy followed by a correct second
    copy -- is what C++ actually evaluates from the BAD row, while this
    verifier silently compares the GOOD one and reports the table clean.
    Row-count asserts do not catch it either: the C++ map and the Python
    dict both collapse to the same (short) length, so the counts still
    agree with each other and with the expected total.

    Verified live: before this guard was added to every parser, injecting a
    bogus `{{"methane", "nitrogen"}, {1.5, 1.5, 1.5, 1.5}}` row above the
    real one in `betasgammas_2004()` left this script at exit 0 with all of
    families 1-5 reported PASS.

    `check_departure`'s specific-table parser and `parse_acentric_rows`
    already had this guard; `parse_bip_rows`, `parse_pure_info_rows`,
    `parse_alphaig_rows`, `parse_fij_rows` and `parse_named_vec_block` did
    not, which is why it lives here as one shared helper rather than five
    hand-written `require`s that can drift apart.
    """
    require(key not in rows, f"{label}: duplicate row for key {key!r} -- a C++ std::map initializer list keeps the FIRST "
                             "of two equal keys while this script would otherwise keep the LAST, so the row that is "
                             "actually compiled would go unverified. Remove the duplicate.")


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def extract_component_names(coolprop_data_h_text):
    m = re.search(r'names_2004\s*=\s*\{(.*?)\};', coolprop_data_h_text, re.DOTALL)
    require(m, "could not find names_2004 in GERGData.h")
    names_2004 = re.findall(r'"([a-z0-9\-]+)"', m.group(1))
    m2 = re.search(r'v\.insert\(v\.end\(\),\s*\{(.*?)\}\)', coolprop_data_h_text, re.DOTALL)
    require(m2, "could not find the GERG-2008 name-extension insert() in GERGData.h")
    extra = re.findall(r'"([a-z0-9\-]+)"', m2.group(1))
    names_2008 = names_2004 + extra
    return names_2004, names_2008


def bounded(text, start_anchor, end_anchor, label):
    """Return the substring of `text` starting at the literal `start_anchor`
    up to (not including) the first `end_anchor` occurring after it."""
    i = text.find(start_anchor)
    require(i >= 0, f"could not find start anchor {start_anchor!r} for {label}")
    j = text.find(end_anchor, i + len(start_anchor))
    require(j >= 0, f"could not find end anchor {end_anchor!r} after start for {label}")
    return text[i + len(start_anchor):j]


# ---------------------------------------------------------------------------
# Family 3: binary reducing parameters (betaV, gammaV, betaT, gammaT)
# ---------------------------------------------------------------------------

# Per the brief: the naive '[a-z0-9.\-]' value character class undercounts
# (152/153) because it cannot span some rows' whitespace forms. This pattern
# (value group is "everything that isn't a brace") reproduces 153/72 exactly.
# Note this also allows optional whitespace between the two quoted names and
# their separating comma: teqp's own GERG.hpp has none ("methane","nitrogen"),
# but clang-format inserts one in CoolProp's copy ("methane", "nitrogen").
BIP_ROW_RE = re.compile(r'\{\{"([a-z0-9\-]+)"\s*,\s*"([a-z0-9\-]+)"\}\s*,\s*\{([^{}]*)\}\s*\}')


def parse_bip_rows(blob, label):
    rows = {}
    for f1, f2, valstr in BIP_ROW_RE.findall(blob):
        vals = parse_floats(valstr)
        require(len(vals) == 4, f"BIP row for ({f1},{f2}) does not have 4 values: {valstr!r}")
        require_new_key(rows, (f1, f2), label)
        rows[(f1, f2)] = tuple(vals)
    return rows


def check_reducing_parameters(teqp_text, coolprop_text):
    print("== Family 3: binary reducing parameters (Table A3.8 / Table A8) ==")

    teqp_2004_blob = bounded(teqp_text, "BIP_data = {", "\n    };\n    auto pair_normal", "teqp GERG2004 BIP_data")
    # The GERG2008 BIP_data block is the SECOND occurrence of this anchor.
    second_start = teqp_text.find("BIP_data = {", teqp_text.find("BIP_data = {") + 1)
    require(second_start >= 0, "could not find second (GERG2008) BIP_data block in teqp")
    teqp_2008_blob = bounded(teqp_text[second_start:], "BIP_data = {", "\n    };\n    auto pair_normal", "teqp GERG2008 BIP_data")

    cp_2004_blob = bounded(coolprop_text, "betasgammas_2004() {\n    static const std::map<BIPKey, BetasGammas> data = {", "\n    };\n    return data;",
                            "CoolProp betasgammas_2004")
    cp_2008_blob = bounded(coolprop_text, "betasgammas_2008_overrides() {\n    static const std::map<BIPKey, BetasGammas> data = {",
                            "\n    };\n    return data;", "CoolProp betasgammas_2008_overrides")

    teqp_2004 = parse_bip_rows(teqp_2004_blob, "teqp GERG2004 BIP_data")
    teqp_2008 = parse_bip_rows(teqp_2008_blob, "teqp GERG2008 BIP_data")
    cp_2004 = parse_bip_rows(cp_2004_blob, "CoolProp betasgammas_2004")
    cp_2008 = parse_bip_rows(cp_2008_blob, "CoolProp betasgammas_2008_overrides")

    require(len(teqp_2004) == 153, f"teqp GERG-2004 BIP row count is {len(teqp_2004)}, expected 153")
    require(len(teqp_2008) == 72, f"teqp GERG-2008 BIP override row count is {len(teqp_2008)}, expected 72")
    require(len(cp_2004) == 153, f"CoolProp GERG-2004 BIP row count is {len(cp_2004)}, expected 153")
    require(len(cp_2008) == 72, f"CoolProp GERG-2008 BIP override row count is {len(cp_2008)}, expected 72")

    mismatches = []
    for key, tvals in teqp_2004.items():
        cvals = cp_2004.get(key)
        if cvals is None:
            mismatches.append(f"2004 pair {key} present in teqp, missing in CoolProp")
        elif not vec_close(tvals, cvals):
            mismatches.append(f"2004 pair {key}: teqp {tvals} != CoolProp {cvals}")
    for key in cp_2004:
        if key not in teqp_2004:
            mismatches.append(f"2004 pair {key} present in CoolProp, missing in teqp")
    for key, tvals in teqp_2008.items():
        cvals = cp_2008.get(key)
        if cvals is None:
            mismatches.append(f"2008-override pair {key} present in teqp, missing in CoolProp")
        elif not vec_close(tvals, cvals):
            mismatches.append(f"2008-override pair {key}: teqp {tvals} != CoolProp {cvals}")
    for key in cp_2008:
        if key not in teqp_2008:
            mismatches.append(f"2008-override pair {key} present in CoolProp, missing in teqp")

    require(not mismatches, "reducing-parameter mismatches:\n  " + "\n  ".join(mismatches))
    print(f"  2004 rows: {len(teqp_2004)}, 2008-override rows: {len(teqp_2008)} -- all match. OK")


# ---------------------------------------------------------------------------
# Family 2: ideal-gas coefficients (n0[1..7], theta0[4..7])
# ---------------------------------------------------------------------------

ALPHAIG_ROW_RE = re.compile(r'\{"([a-z0-9\-]+)",\s*\{\{([^{}]*)\},\s*\{([^{}]*)\}\}\}', re.DOTALL)


def parse_alphaig_rows(blob, label):
    rows = {}
    for name, n0str, theta0str in ALPHAIG_ROW_RE.findall(blob):
        n0 = parse_floats(n0str)
        theta0 = parse_floats(theta0str)
        require(len(n0) == 7, f"{name}: alphaig n0 does not have 7 values: {n0str!r}")
        require(len(theta0) == 4, f"{name}: alphaig theta0 does not have 4 values: {theta0str!r}")
        require_new_key(rows, name, label)
        rows[name] = (tuple(n0), tuple(theta0))
    return rows


def check_ideal_gas(teqp_text, coolprop_text, names_2004, names_2008):
    print("== Family 2: ideal-gas coefficients (Table A3.1) ==")

    teqp_2004_blob = bounded(teqp_text, "inline AlphaigCoeffs get_alphaig_coeffs(const std::string& fluid){", "\n    if (dict.find(fluid)",
                              "teqp GERG2004 alphaig dict")
    gerg2008_ns = teqp_text[teqp_text.find("namespace GERG2008{"):]
    teqp_2008_blob = bounded(gerg2008_ns, "inline AlphaigCoeffs get_alphaig_coeffs(const std::string& fluid){", "\n    if (dict.find(fluid)",
                              "teqp GERG2008 alphaig dict")

    cp_2004_blob = bounded(coolprop_text, "alphaig_2004() {\n    static const std::map<std::string, RawAlphaig> data = {", "};\n    return data;",
                            "CoolProp alphaig_2004")
    cp_2008_blob = bounded(coolprop_text, "alphaig_2008_overrides() {\n    static const std::map<std::string, RawAlphaig> data = {",
                            "};\n    return data;", "CoolProp alphaig_2008_overrides")

    teqp_2004 = parse_alphaig_rows(teqp_2004_blob, "teqp GERG2004 alphaig dict")
    teqp_2008 = parse_alphaig_rows(teqp_2008_blob, "teqp GERG2008 alphaig dict")
    cp_2004 = parse_alphaig_rows(cp_2004_blob, "CoolProp alphaig_2004")
    cp_2008 = parse_alphaig_rows(cp_2008_blob, "CoolProp alphaig_2008_overrides")

    require(len(teqp_2004) == 18, f"teqp GERG-2004 alphaig row count is {len(teqp_2004)}, expected 18")
    require(len(teqp_2008) == 5, f"teqp GERG-2008 alphaig override row count is {len(teqp_2008)}, expected 5")
    require(len(cp_2004) == 18, f"CoolProp GERG-2004 alphaig row count is {len(cp_2004)}, expected 18")
    require(len(cp_2008) == 5, f"CoolProp GERG-2008 alphaig override row count is {len(cp_2008)}, expected 5")

    def resolve(name, teqp_base, teqp_ov, cp_base, cp_ov):
        t = teqp_ov.get(name, teqp_base.get(name))
        c = cp_ov.get(name, cp_base.get(name))
        return t, c

    mismatches = []
    pairs_checked = 0
    for model, names in (("GERG_2004", names_2004), ("GERG_2008", names_2008)):
        base_t, ov_t = (teqp_2004, {}) if model == "GERG_2004" else (teqp_2004, teqp_2008)
        base_c, ov_c = (cp_2004, {}) if model == "GERG_2004" else (cp_2004, cp_2008)
        for name in names:
            tval, cval = resolve(name, base_t, ov_t, base_c, ov_c)
            if tval is None:
                mismatches.append(f"{model}/{name}: no teqp alphaig row found")
                continue
            if cval is None:
                mismatches.append(f"{model}/{name}: no CoolProp alphaig row found")
                continue
            tn0, tth0 = tval
            cn0, cth0 = cval
            if not vec_close(tn0, cn0):
                mismatches.append(f"{model}/{name}: n0 teqp {tn0} != CoolProp {cn0}")
            if not vec_close(tth0, cth0):
                mismatches.append(f"{model}/{name}: theta0 teqp {tth0} != CoolProp {cth0}")
            pairs_checked += 1

    require(pairs_checked == 39, f"checked {pairs_checked} (model, component) ideal-gas pairs, expected 39")
    require(not mismatches, "ideal-gas mismatches:\n  " + "\n  ".join(mismatches))
    print(f"  (model, component) pairs checked: {pairs_checked} -- all match. OK")


# ---------------------------------------------------------------------------
# Family 1: pure-fluid residual coefficients (n, t, d, c, l)
# ---------------------------------------------------------------------------

NAMED_VEC_RE = re.compile(r'\{"([a-z0-9\-]+)",\s*\{([^{}]*)\}\}', re.DOTALL)
QUAD_AFTER_RE = lambda anchor: re.compile(
    re.escape(anchor) + r"\s*;\s*pc\.t\s*=\s*\{([^{}]*)\};\s*pc\.d\s*=\s*\{([^{}]*)\};\s*pc\.c\s*=\s*\{([^{}]*)\};\s*pc\.l\s*=\s*\{([^{}]*)\};"
)
BESPOKE_RE = lambda name: re.compile(
    r'fluid\s*==\s*"' + re.escape(name) + r'"\)\{\s*PureCoeffs pc;\s*pc\.n\s*=\s*\{([^{}]*)\};\s*pc\.t\s*=\s*\{([^{}]*)\};'
    r'\s*pc\.d\s*=\s*\{([^{}]*)\};\s*pc\.c\s*=\s*\{([^{}]*)\};\s*pc\.l\s*=\s*\{([^{}]*)\};',
    re.DOTALL,
)


def parse_named_vec_block(blob, label):
    rows = {}
    for name, v in NAMED_VEC_RE.findall(blob):
        require_new_key(rows, name, label)
        rows[name] = tuple(parse_floats(v))
    return rows


def teqp_shared_quad(text, anchor):
    m = QUAD_AFTER_RE(anchor).search(text)
    require(m, f"could not find shared t/d/c/l quad after {anchor!r} in teqp")
    return tuple(tuple(parse_floats(g)) for g in m.groups())


def teqp_bespoke(text, name):
    m = BESPOKE_RE(name).search(text)
    require(m, f"could not find bespoke pure-fluid block for {name!r} in teqp")
    n, t, d, c, l = (tuple(parse_floats(g)) for g in m.groups())
    return n, t, d, c, l


def cp_shared_exponents(text, fn_name):
    blob = bounded(text, f"{fn_name}() {{\n    static const SharedExponents e = {{", "};\n    return e;", f"CoolProp {fn_name}")
    # Four sibling bracketed groups: t, d, c, l (see SharedExponents{t,d,c,l}).
    m = re.match(r'\s*\{([^{}]*)\}\s*,\s*\{([^{}]*)\}\s*,\s*\{([^{}]*)\}\s*,\s*\{([^{}]*)\}\s*', blob, re.DOTALL)
    require(m, f"could not split {fn_name} into 4 exponent groups")
    return tuple(tuple(parse_floats(g)) for g in m.groups())


def cp_bespoke(text, fn_name):
    blob = bounded(text, f"PureCoeffs {fn_name}() {{\n    return {{", "};\n}", f"CoolProp {fn_name}")
    m = re.match(r'\s*\{([^{}]*)\}\s*,\s*\{([^{}]*)\}\s*,\s*\{([^{}]*)\}\s*,\s*\{([^{}]*)\}\s*,\s*\{([^{}]*)\}\s*', blob, re.DOTALL)
    require(m, f"could not split {fn_name} into 5 fields (n,t,d,c,l)")
    return tuple(tuple(parse_floats(g)) for g in m.groups())


BESPOKE_FLUIDS = ["carbondioxide", "hydrogen", "water", "helium"]
BESPOKE_CP_FN = {
    "carbondioxide": "carbondioxide_2004",
    "hydrogen": "hydrogen_2004",
    "water": "water_2004",
    "helium": "helium_2004",
}


def check_pure_residual(teqp_text, coolprop_text, names_2004, names_2008):
    print("== Family 1: pure-fluid residual coefficients (Table A3.2) ==")

    teqp_2004_fn = bounded(teqp_text, "inline PureCoeffs get_pure_coeffs(const std::string& fluid){", "\n\n\ninline BetasGammas",
                            "teqp GERG2004 get_pure_coeffs")
    gerg2008_ns = teqp_text[teqp_text.find("namespace GERG2008{"):]
    teqp_2008_fn = bounded(gerg2008_ns, "inline PureCoeffs get_pure_coeffs(const std::string& fluid){", "\ninline AlphaigCoeffs",
                            "teqp GERG2008 get_pure_coeffs")

    # main12/mne24 n-only maps (2004) and shared exponents.
    teqp_main_n = parse_named_vec_block(bounded(teqp_2004_fn, "n_dict_main = {", "\n    };", "teqp n_dict_main (2004)"), "teqp n_dict_main (2004)")
    teqp_mne_n = parse_named_vec_block(bounded(teqp_2004_fn, "n_dict_mne = {", "\n    };", "teqp n_dict_mne (2004)"), "teqp n_dict_mne (2004)")
    require(len(teqp_main_n) == 11, f"teqp n_dict_main (2004) has {len(teqp_main_n)} fluids, expected 11")
    require(len(teqp_mne_n) == 3, f"teqp n_dict_mne (2004) has {len(teqp_mne_n)} fluids, expected 3")
    teqp_main_tdcl = teqp_shared_quad(teqp_2004_fn, "pc.n = n_dict_main[fluid]")
    teqp_mne_tdcl = teqp_shared_quad(teqp_2004_fn, "pc.n = n_dict_mne.at(fluid)")

    teqp_bespoke_2004 = {name: teqp_bespoke(teqp_2004_fn, name) for name in BESPOKE_FLUIDS}

    teqp_main_n_2008 = parse_named_vec_block(bounded(teqp_2008_fn, "n_dict_main = {", "\n    };", "teqp n_dict_main (2008 override)"),
                                              "teqp n_dict_main (2008 override)")
    require(len(teqp_main_n_2008) == 5, f"teqp n_dict_main (2008 override) has {len(teqp_main_n_2008)} fluids, expected 5")

    # CoolProp side.
    cp_main_n = parse_named_vec_block(bounded(coolprop_text, "n_main_2004() {\n    static const std::map<std::string, std::vector<double>> data = {",
                                               "};\n    return data;", "CoolProp n_main_2004"), "CoolProp n_main_2004")
    cp_mne_n = parse_named_vec_block(bounded(coolprop_text, "n_mne_2004() {\n    static const std::map<std::string, std::vector<double>> data = {",
                                              "};\n    return data;", "CoolProp n_mne_2004"), "CoolProp n_mne_2004")
    require(len(cp_main_n) == 11, f"CoolProp n_main_2004 has {len(cp_main_n)} fluids, expected 11")
    require(len(cp_mne_n) == 3, f"CoolProp n_mne_2004 has {len(cp_mne_n)} fluids, expected 3")

    cp_main_tdcl = cp_shared_exponents(coolprop_text, "main12_exponents")
    cp_mne_tdcl = cp_shared_exponents(coolprop_text, "mne24_exponents")

    cp_bespoke_2004 = {name: cp_bespoke(coolprop_text, BESPOKE_CP_FN[name]) for name in BESPOKE_FLUIDS}

    cp_main_n_2008 = parse_named_vec_block(
        bounded(coolprop_text, "n_main_2008_overrides() {\n    static const std::map<std::string, std::vector<double>> data = {",
                "};\n    return data;", "CoolProp n_main_2008_overrides"), "CoolProp n_main_2008_overrides")
    require(len(cp_main_n_2008) == 5, f"CoolProp n_main_2008_overrides has {len(cp_main_n_2008)} fluids, expected 5")

    def resolve_teqp(name, model):
        if model == "GERG_2008" and name in teqp_main_n_2008:
            return (teqp_main_n_2008[name],) + teqp_main_tdcl
        if name in BESPOKE_FLUIDS:
            n, t, d, c, l = teqp_bespoke_2004[name]
            return n, t, d, c, l
        if name in teqp_main_n:
            return (teqp_main_n[name],) + teqp_main_tdcl
        if name in teqp_mne_n:
            return (teqp_mne_n[name],) + teqp_mne_tdcl
        return None

    def resolve_cp(name, model):
        if model == "GERG_2008" and name in cp_main_n_2008:
            return (cp_main_n_2008[name],) + cp_main_tdcl
        if name in BESPOKE_FLUIDS:
            n, t, d, c, l = cp_bespoke_2004[name]
            return n, t, d, c, l
        if name in cp_main_n:
            return (cp_main_n[name],) + cp_main_tdcl
        if name in cp_mne_n:
            return (cp_mne_n[name],) + cp_mne_tdcl
        return None

    mismatches = []
    distinct_eos = set()
    for model, names in (("GERG_2004", names_2004), ("GERG_2008", names_2008)):
        for name in names:
            tres = resolve_teqp(name, model)
            cres = resolve_cp(name, model)
            require(tres is not None, f"{model}/{name}: no teqp pure-residual row resolved")
            require(cres is not None, f"{model}/{name}: no CoolProp pure-residual row resolved")
            for field, tv, cv in zip("ntdcl", tres, cres):
                if not vec_close(tv, cv):
                    mismatches.append(f"{model}/{name}: field {field} teqp {tv} != CoolProp {cv}")
            n, t, d, c, l = cres
            for i, (ci, li) in enumerate(zip(c, l)):
                expect = 1.0 if li > 0 else 0.0
                if ci != expect:
                    mismatches.append(f"{model}/{name}: c[{i}]={ci} but l[{i}]={li} (c/l invariant violated)")
            distinct_eos.add(cres)

    require(not mismatches, "pure-residual mismatches:\n  " + "\n  ".join(mismatches))
    require(len(distinct_eos) == 23, f"found {len(distinct_eos)} distinct pure EOS tuples, expected 23")
    print(f"  distinct pure EOS: {len(distinct_eos)} -- all match, c/l invariant holds. OK")


# ---------------------------------------------------------------------------
# Family 4: departure functions and F_ij (Table A3.6)
# ---------------------------------------------------------------------------

FIJ_ROW_RE = re.compile(r'\{\{"([a-z0-9\-]+)"\s*,\s*"([a-z0-9\-]+)"\}\s*,\s*(' + FLOAT_RE + r')\s*\}')


def parse_fij_rows(blob, label):
    """A named function rather than the dict comprehension this used to be,
    solely so require_new_key has somewhere to live -- see its docstring."""
    rows = {}
    for f1, f2, v in FIJ_ROW_RE.findall(blob):
        require_new_key(rows, (f1, f2), label)
        rows[(f1, f2)] = float(v)
    return rows

# Departure-coefficient rows are multi-line dc.FIELD = {...}; assignments.
# This pattern captures one field assignment; a full row is 7 consecutive
# assignments (n,d,t,eta,epsilon,beta,gamma for teqp; n,t,d,eta,beta,gamma,
# epsilon for CoolProp's DepartureCoeffs -- field NAMES match, only the
# on-disk order differs, so this is parsed by name, not position).
DC_FIELD_RE = re.compile(r'dc\.(n|d|t|eta|epsilon|beta|gamma)\s*=\s*\{([^{}]*)\}\s*;', re.DOTALL)


def row_key(row):
    """A hashable, order-independent fingerprint for one departure row."""
    return tuple(row[f] for f in ("n", "d", "t", "eta", "epsilon", "beta", "gamma"))


def canon_pair(f1, f2):
    return tuple(sorted((f1, f2)))


# One `if (sortedpair == sortpair("f1","f2")){ ... return dc; }` branch per
# fluid-pair-specific departure function -- captures which PAIR each row
# belongs to, not just the row's content.
SPECIFIC_IF_RE = re.compile(r'if \(sortedpair == sortpair\("([a-z0-9\-]+)"\s*,\s*"([a-z0-9\-]+)"\)\)\{(.*?)return dc\s*;', re.DOTALL)

# The `generalized` std::set<pair<string,string>> that selects which pairs
# use the shared generalized departure function, and the one branch that
# returns it.
GENERALIZED_SET_RE = re.compile(r'sortpair\("([a-z0-9\-]+)"\s*,\s*"([a-z0-9\-]+)"\)')
GENERALIZED_ROW_RE = re.compile(r'if \(generalized\.find\(sortedpair\) != generalized\.end\(\)\)\{(.*?)return dc\s*;', re.DOTALL)

# CoolProp: departure_specific_table()'s map keys, each pointing at a named
# departure_<pair>() function; and generalized_departure_pairs()'s contents.
CP_SPECIFIC_TABLE_ROW_RE = re.compile(r'\{\{"([a-z0-9\-]+)"\s*,\s*"([a-z0-9\-]+)"\}\s*,\s*(departure_\w+)\(\)\}')
CP_NAMED_FN_RE = re.compile(r'DepartureCoeffs (departure_\w+|generalized_departure)\(\) \{(.*?)return dc\s*;', re.DOTALL)
CP_GENERALIZED_PAIR_RE = re.compile(r'\{"([a-z0-9\-]+)"\s*,\s*"([a-z0-9\-]+)"\}')


def parse_row_body(body):
    row = {}
    for name, valstr in DC_FIELD_RE.findall(body):
        row[name] = tuple(parse_floats(valstr))
    require(len(row) == 7 and len({len(v) for v in row.values()}) == 1, f"malformed departure row body: {body[:200]!r}")
    return row


def build_teqp_pair_rows(teqp_dep_blob):
    """PAIR-KEYED (not just content-keyed) map of every one of the 15
    departure pairs teqp defines to its 7-field row. This is what actually
    proves a pair is wired to the RIGHT row, not merely that 8 distinct row
    shapes exist somewhere in the source."""
    generalized_set_blob = bounded(teqp_dep_blob, "const std::set<std::pair<std::string, std::string>> generalized = {", "\n    };",
                                    "teqp generalized set")
    generalized_pairs = {canon_pair(f1, f2) for f1, f2 in GENERALIZED_SET_RE.findall(generalized_set_blob)}
    require(len(generalized_pairs) == 8, f"teqp `generalized` set has {len(generalized_pairs)} pairs, expected 8")

    specific = {}
    for f1, f2, body in SPECIFIC_IF_RE.findall(teqp_dep_blob):
        pair = canon_pair(f1, f2)
        require(pair not in specific, f"teqp: pair {pair} appears in more than one specific if-branch")
        specific[pair] = parse_row_body(body)
    require(len(specific) == 7, f"teqp has {len(specific)} fluid-pair-specific departure branches, expected 7")

    gen_m = GENERALIZED_ROW_RE.search(teqp_dep_blob)
    require(gen_m, "could not find teqp's generalized-departure-function branch")
    gen_row = parse_row_body(gen_m.group(1))

    pair_rows = dict(specific)
    for pair in generalized_pairs:
        require(pair not in pair_rows, f"teqp: pair {pair} is claimed by both a specific branch and the generalized set")
        pair_rows[pair] = gen_row
    require(len(pair_rows) == 15, f"teqp resolves {len(pair_rows)} total departure pairs, expected 15")
    return pair_rows


def build_cp_pair_rows(cp_dep_blob):
    """Same PAIR-KEYED map as build_teqp_pair_rows, built from CoolProp's
    departure_specific_table() (pair -> named function) plus
    generalized_departure_pairs() (pair -> the shared generalized row),
    resolving each named function to its own parsed row. This is what
    actually catches a miswired pair (e.g. "methane","nitrogen" pointing at
    departure_methane_ethane()) or a pair wrongly added to/missing from the
    generalized set -- comparing unkeyed row content alone cannot."""
    fn_rows = {}
    for fn_name, body in CP_NAMED_FN_RE.findall(cp_dep_blob):
        fn_rows[fn_name] = parse_row_body(body)
    require("generalized_departure" in fn_rows, "could not find CoolProp's generalized_departure() body")

    specific_table_blob = bounded(cp_dep_blob, "departure_specific_table() {\n    static const std::map<BIPKey, DepartureCoeffs> data = {",
                                   "\n    };\n    return data;", "CoolProp departure_specific_table")
    specific = {}
    for f1, f2, fn_name in CP_SPECIFIC_TABLE_ROW_RE.findall(specific_table_blob):
        pair = canon_pair(f1, f2)
        require(fn_name in fn_rows, f"departure_specific_table() references undefined function {fn_name}() for pair {pair}")
        require(pair not in specific, f"CoolProp: pair {pair} appears twice in departure_specific_table()")
        specific[pair] = fn_rows[fn_name]
    require(len(specific) == 7, f"CoolProp departure_specific_table() has {len(specific)} rows, expected 7")

    gen_pairs_blob = bounded(cp_dep_blob, "generalized_departure_pairs() {\n    static const std::vector<BIPKey> pairs = {", "\n    };\n    return pairs;",
                              "CoolProp generalized_departure_pairs")
    generalized_pairs = {canon_pair(f1, f2) for f1, f2 in CP_GENERALIZED_PAIR_RE.findall(gen_pairs_blob)}
    require(len(generalized_pairs) == 8, f"CoolProp generalized_departure_pairs() has {len(generalized_pairs)} pairs, expected 8")

    pair_rows = dict(specific)
    for pair in generalized_pairs:
        require(pair not in pair_rows, f"CoolProp: pair {pair} is claimed by both departure_specific_table() and "
                                        "generalized_departure_pairs()")
        pair_rows[pair] = fn_rows["generalized_departure"]
    require(len(pair_rows) == 15, f"CoolProp resolves {len(pair_rows)} total departure pairs, expected 15")
    return pair_rows


def check_departure(teqp_text, coolprop_text):
    print("== Family 4: departure functions and F_ij (Table A3.6) ==")

    teqp_fij_blob = bounded(teqp_text, "static std::map<std::pair<std::string, std::string>, double> Fij_dict = {", "\n    };",
                             "teqp Fij_dict")
    cp_fij_blob = bounded(coolprop_text, "fij_table() {\n    static const std::map<BIPKey, double> data = {", "\n    };\n    return data;",
                           "CoolProp fij_table")

    teqp_fij = parse_fij_rows(teqp_fij_blob, "teqp Fij_dict")
    cp_fij = parse_fij_rows(cp_fij_blob, "CoolProp fij_table")

    require(len(teqp_fij) == 15, f"teqp Fij_dict has {len(teqp_fij)} rows, expected 15")
    require(len(cp_fij) == 15, f"CoolProp fij_table has {len(cp_fij)} rows, expected 15")

    mismatches = []
    for key, tval in teqp_fij.items():
        cval = cp_fij.get(key)
        if cval is None:
            mismatches.append(f"F_ij pair {key} present in teqp, missing in CoolProp")
        elif not close(tval, cval, 1e-11):
            mismatches.append(f"F_ij pair {key}: teqp {tval} != CoolProp {cval}")
    for key in cp_fij:
        if key not in teqp_fij:
            mismatches.append(f"F_ij pair {key} present in CoolProp, missing in teqp")
    with_scaled_f = sum(1 for v in teqp_fij.values() if v != 1.0)
    require(with_scaled_f == 7, f"teqp Fij_dict has {with_scaled_f} entries with F_ij != 1, expected 7")

    # Departure coefficients: get_departurecoeffs (teqp) vs the specific-pair
    # functions + generalized_departure() (CoolProp). PAIR-KEYED: this
    # verifies which PAIR each row belongs to, not just that the same 8
    # distinct row shapes exist somewhere on both sides. An earlier version
    # of this check compared unkeyed sets of row fingerprints, which cannot
    # detect a miswired pair (e.g. methane/nitrogen accidentally pointing at
    # methane/ethane's row) or a pair wrongly added to/missing from the
    # generalized set -- the set of 8 distinct blobs would be unchanged.
    # build_teqp_pair_rows/build_cp_pair_rows resolve each side's own
    # pair -> function/branch -> row wiring from the source text itself,
    # so a miswiring shows up as a value mismatch (or a missing/extra pair)
    # for that specific pair, exactly like the F_ij check above.
    teqp_dep_blob = bounded(teqp_text, "inline DepartureCoeffs get_departurecoeffs(const std::string&fluid1, const std::string &fluid2){",
                            '\n    throw std::invalid_argument("could not get departure coeffs', "teqp get_departurecoeffs")
    cp_dep_blob = bounded(coolprop_text, "// Departure-function tables: F_ij scaling factors and departure coefficients.",
                          "bool is_generalized_departure_pair", "CoolProp departure tables")

    teqp_pairs = build_teqp_pair_rows(teqp_dep_blob)
    cp_pairs = build_cp_pair_rows(cp_dep_blob)

    dep_mismatches = []
    for pair, trow in teqp_pairs.items():
        crow = cp_pairs.get(pair)
        if crow is None:
            dep_mismatches.append(f"departure pair {pair} present in teqp, missing in CoolProp")
            continue
        for field in ("n", "d", "t", "eta", "epsilon", "beta", "gamma"):
            if not vec_close(trow[field], crow[field]):
                dep_mismatches.append(f"departure pair {pair} field {field}: teqp {trow[field]} != CoolProp {crow[field]}")
    for pair in cp_pairs:
        if pair not in teqp_pairs:
            dep_mismatches.append(f"departure pair {pair} present in CoolProp, missing in teqp")

    # Report the pair-keyed mismatches (most specific: names the exact pair
    # and field) before the secondary distinct-shape count below, so a
    # miswired pair is diagnosed by WHICH pair, not just a raw count.
    require(not dep_mismatches, "departure-coefficient mismatches:\n  " + "\n  ".join(dep_mismatches))
    require(not mismatches, "F_ij mismatches:\n  " + "\n  ".join(mismatches))

    # Secondary sanity check retained from the original version: the number
    # of DISTINCT row shapes should still be 8 (7 specific + 1 generalized)
    # on both sides -- catches accidental duplication/typo'd near-duplicate
    # rows that the pair-keyed check above wouldn't flag as a "missing pair"
    # (e.g. two genuinely different pairs both accidentally wired to the
    # same specific-function name, keeping every pair "present" but wrong).
    teqp_distinct = {row_key(r) for r in teqp_pairs.values()}
    cp_distinct = {row_key(r) for r in cp_pairs.values()}
    require(len(teqp_distinct) == 8, f"teqp departure rows resolve to {len(teqp_distinct)} distinct shapes, expected 8")
    require(len(cp_distinct) == 8, f"CoolProp departure rows resolve to {len(cp_distinct)} distinct shapes, expected 8")
    print(f"  F_ij rows: {len(teqp_fij)} (scaled: {with_scaled_f}), departure pairs: {len(teqp_pairs)} "
          f"(distinct row shapes: {len(teqp_distinct)}) -- all match, pair-keyed. OK")


# ---------------------------------------------------------------------------
# Family 5: pure-fluid reducing/critical parameters (Table A3.5)
# ---------------------------------------------------------------------------

PURE_INFO_ROW_RE = re.compile(r'\{"([a-z0-9\-]+)"\s*,\s*\{([^{}]*)\}\}')


def parse_pure_info_rows(blob, label):
    """name -> (raw literal tokens, parsed floats), both length-3
    (rhoc_molm3, Tc_K, M_kgmol) in that order -- teqp's PureInfo brace-init
    order and CoolProp's PureInfo field order match, so no reordering is
    needed here (unlike BetasGammas' betaV/gammaV/betaT/gammaT)."""
    rows = {}
    for name, valstr in PURE_INFO_ROW_RE.findall(blob):
        raw = [v.strip() for v in valstr.split(",") if v.strip() != ""]
        require(len(raw) == 3, f"{name}: pure_info row does not have 3 values: {valstr!r}")
        require_new_key(rows, name, label)
        rows[name] = (tuple(raw), tuple(parse_floats(valstr)))
    return rows


def _norm_literal(tok):
    """Formatting-insensitive normalisation of a numeric literal TOKEN, for a
    raw-token compare that still tolerates pure formatting (unary '+',
    trailing zeros) without being a second float parse. Deliberately
    narrower than parse_floats/close(): it exists because a tolerance-based
    float compare can silently accept a transcription slip that changes the
    literal but not the parsed value enough to trip `close()`'s tolerance --
    e.g. a truncated digit inside a long mantissa. Comparing the normalised
    tokens as well as the parsed floats catches that class of error; two
    literals that mean the same number after this normalisation (e.g.
    "1.00000000" and "1.000000000") are treated as equivalent, not as a
    mismatch, since only the trailing-zero formatting differs.
    """
    t = tok.strip()
    if t.startswith("+"):
        t = t[1:]
    if "." in t and "e" not in t and "E" not in t:
        t = t.rstrip("0")
        if t.endswith("."):
            t += "0"
    return t


def check_pure_info(teqp_text, coolprop_data_h_text):
    print("== Family 5: pure-fluid reducing/critical parameters (Table A3.5) ==")

    # UNITS: this compares the RAW TABLE LITERALS as they appear in source --
    # mol/dm^3 for rhoc, K for Tc, kg/kmol for M -- i.e. BEFORE the
    # mol/m^3 / kg/mol conversion that get_pure_info() applies on the way out
    # (GERGData.h's pure_info_2004()/pure_info_2008_overrides() are the
    # pre-conversion tables; GERGData.h::get_pure_info() does `*= 1000` /
    # `/= 1000` afterwards, and teqp's get_pure_info does the identical
    # conversion). Comparing pre-conversion literals on both sides is "like
    # with like"; running either side's post-conversion accessor would only
    # re-verify the same *1000/1000 arithmetic against itself, not the
    # transcription of the tabulated values.
    teqp_start_anchor = "static std::map<std::string, PureInfo> data_map = {"
    teqp_2004_blob = bounded(teqp_text, teqp_start_anchor, "\n    };\n    if (data_map.find(name) == data_map.end())",
                              "teqp GERG2004 pure_info data_map")
    second_start = teqp_text.find(teqp_start_anchor, teqp_text.find(teqp_start_anchor) + 1)
    require(second_start >= 0, "could not find second (GERG2008) pure_info data_map block in teqp")
    teqp_2008_blob = bounded(teqp_text[second_start:], teqp_start_anchor, "\n    };\n    if (data_map.find(name) != data_map.end())",
                              "teqp GERG2008 pure_info data_map")

    cp_2004_blob = bounded(coolprop_data_h_text, "pure_info_2004() {\n    static const std::map<std::string, PureInfo> data = {",
                            "};\n    return data;", "CoolProp pure_info_2004")
    cp_2008_blob = bounded(coolprop_data_h_text, "pure_info_2008_overrides() {\n    static const std::map<std::string, PureInfo> data = {",
                            "};\n    return data;", "CoolProp pure_info_2008_overrides")

    teqp_2004 = parse_pure_info_rows(teqp_2004_blob, "teqp GERG2004 pure_info data_map")
    teqp_2008 = parse_pure_info_rows(teqp_2008_blob, "teqp GERG2008 pure_info data_map")
    cp_2004 = parse_pure_info_rows(cp_2004_blob, "CoolProp pure_info_2004")
    cp_2008 = parse_pure_info_rows(cp_2008_blob, "CoolProp pure_info_2008_overrides")

    require(len(teqp_2004) == 18, f"teqp GERG-2004 pure_info row count is {len(teqp_2004)}, expected 18")
    require(len(teqp_2008) == 5, f"teqp GERG-2008 pure_info override row count is {len(teqp_2008)}, expected 5")
    require(len(cp_2004) == 18, f"CoolProp GERG-2004 pure_info row count is {len(cp_2004)}, expected 18")
    require(len(cp_2008) == 5, f"CoolProp GERG-2008 pure_info override row count is {len(cp_2008)}, expected 5")

    mismatches = []

    def compare(model_label, teqp_rows, cp_rows):
        for name, (traw, tvals) in teqp_rows.items():
            entry = cp_rows.get(name)
            if entry is None:
                mismatches.append(f"{model_label} {name}: present in teqp, missing in CoolProp")
                continue
            craw, cvals = entry
            if not vec_close(tvals, cvals):
                mismatches.append(f"{model_label} {name}: teqp {tvals} != CoolProp {cvals}")
            if tuple(_norm_literal(v) for v in traw) != tuple(_norm_literal(v) for v in craw):
                mismatches.append(f"{model_label} {name}: raw literal mismatch teqp {traw} != CoolProp {craw}")
        for name in cp_rows:
            if name not in teqp_rows:
                mismatches.append(f"{model_label} {name}: present in CoolProp, missing in teqp")

    compare("2004", teqp_2004, cp_2004)
    compare("2008-override", teqp_2008, cp_2008)

    require(not mismatches, "pure_info mismatches:\n  " + "\n  ".join(mismatches))
    print(f"  2004 rows: {len(teqp_2004)}, 2008-override rows: {len(teqp_2008)} -- all match (float + raw-literal). OK")


# ---------------------------------------------------------------------------
# Family 6: acentric factors (DERIVED, not transcribed)
# ---------------------------------------------------------------------------

ACENTRIC_ROW_RE = re.compile(r'\{"([a-z0-9\-]+)"\s*,\s*(' + FLOAT_RE + r')\}')

#: Largest tolerated |omega_committed - omega_rederived|.  The solve is a
#: bisection plus a Newton polish to a fixed residual, so it reproduces to
#: ~1e-12 across numpy/teqp builds -- unlike the ancillary fits, whose last
#: digits move with the least-squares implementation.  1e-9 absorbs that
#: without absorbing a real change: the smallest meaningful edit to a
#: coefficient table moves omega by ~1e-6 at worst.
ACENTRIC_TOL = 1e-9


def parse_acentric_rows(blob):
    """name -> omega, from one acentric_*() function body."""
    rows = {}
    for name, val in ACENTRIC_ROW_RE.findall(blob):
        require(name not in rows, f"{name}: duplicate acentric row")
        rows[name] = float(val)
    return rows


def check_acentric():
    """Re-derive every acentric factor from teqp's GERG implementation and diff
    it against the literals committed in GERGAcentric.h.

    UNLIKE THE OTHER FIVE FAMILIES this is not a transcription check, and it
    cannot be done by regex against teqp's source: GERG publishes no acentric
    factor, so there is no table in GERG.hpp to compare against.  The value is
    DERIVED -- omega = -1 - log10(p_sat(0.7 T_c)/p_c) evaluated on the GERG
    pure equation -- so verifying it means recomputing it.  That needs teqp's
    PYTHON MODULE at run time, which the other five checks deliberately avoid.

    What it catches, and what it does not.  It catches the failure mode a
    generated table actually has: someone edits a pure-fluid coefficient row or
    a Table A3.5 reducing parameter and forgets to rerun
    dev/gerg/compute_acentric.py, leaving a table that no longer matches the
    equation it claims to describe.  It does NOT independently confirm the VLE
    solver, because it reuses compute_acentric.py's own machinery -- that
    confirmation lives in the Catch2 case "GERG acentric factors satisfy their
    own definition", which re-derives the same 23 values through CoolProp's
    entirely separate saturation solver.
    """
    print("== Family 6: acentric factors (derived from the GERG EOS) ==")

    # Imported here rather than at module scope so that the five regex-only
    # families keep working in an environment without teqp installed --
    # main() decides whether reaching this point without teqp is an error.
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from compute_acentric import compute_all  # noqa: PLC0415

    header = GERG_ACENTRIC_H.read_text()
    cp_2004 = parse_acentric_rows(bounded(header, "acentric_2004()", "return data;", "acentric_2004"))
    cp_2008 = parse_acentric_rows(bounded(header, "acentric_2008_overrides()", "return data;", "acentric_2008_overrides"))
    require(len(cp_2004) == 18, f"parsed {len(cp_2004)} GERG-2004 acentric rows, expected 18")
    require(len(cp_2008) == 5, f"parsed {len(cp_2008)} GERG-2008 acentric override rows, expected 5")

    committed = {(2004, n): v for n, v in cp_2004.items()}
    committed.update({(2008, n): v for n, v in cp_2008.items()})

    mismatches = []
    rederived = compute_all()
    require(len(rederived) == 23, f"re-derived {len(rederived)} acentric factors, expected 23")
    for r in rederived:
        key = (r.year, r.name)
        if key not in committed:
            mismatches.append(f"{r.year}/{r.name}: re-derived but missing from GERGAcentric.h")
            continue
        delta = abs(committed[key] - r.omega)
        if not (delta <= ACENTRIC_TOL):
            mismatches.append(
                f"{r.year}/{r.name}: committed {committed[key]!r} != re-derived {r.omega!r} (|delta| = {delta:.3e} > {ACENTRIC_TOL:.1e})"
            )
    for key in committed:
        if not any((r.year, r.name) == key for r in rederived):
            mismatches.append(f"{key[0]}/{key[1]}: present in GERGAcentric.h, not re-derived")

    require(not mismatches, "acentric mismatches (rerun dev/gerg/compute_acentric.py?):\n  " + "\n  ".join(mismatches))
    worst = max(abs(committed[(r.year, r.name)] - r.omega) for r in rederived)
    print(f"  23 rows re-derived from the GERG EOS -- all match (worst |delta| = {worst:.3e}). OK")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--teqp-root", default=os.path.expanduser("~/Code/teqp"), help="Path to a teqp checkout (default: ~/Code/teqp)")
    ap.add_argument(
        "--skip-acentric",
        action="store_true",
        help="Skip family 6 (acentric factors). Family 6 is the one check that needs teqp's PYTHON MODULE "
        "at run time, not just a source checkout; without this flag a missing teqp module is a FAILURE, "
        "not a silent skip.",
    )
    args = ap.parse_args()

    teqp_hpp = Path(args.teqp_root) / "include" / "teqp" / "models" / "GERG" / "GERG.hpp"
    if not teqp_hpp.is_file():
        fail(f"teqp reference source not found at {teqp_hpp} -- pass --teqp-root, or clone "
             "https://github.com/usnistgov/teqp there. This script is a developer tool: it needs "
             "a local teqp checkout and is not part of CoolProp's build or CI.")

    teqp_text = teqp_hpp.read_text()
    coolprop_text = GERG_BACKEND_CPP.read_text()
    coolprop_data_h_text = GERG_DATA_H.read_text()

    names_2004, names_2008 = extract_component_names(coolprop_data_h_text)
    require(len(names_2004) == 18, f"parsed {len(names_2004)} GERG-2004 component names, expected 18")
    require(len(names_2008) == 21, f"parsed {len(names_2008)} GERG-2008 component names, expected 21")

    check_pure_residual(teqp_text, coolprop_text, names_2004, names_2008)
    check_ideal_gas(teqp_text, coolprop_text, names_2004, names_2008)
    check_reducing_parameters(teqp_text, coolprop_text)
    check_departure(teqp_text, coolprop_text)
    check_pure_info(teqp_text, coolprop_data_h_text)

    if args.skip_acentric:
        # An EXPLICIT opt-out, printed loudly and reflected in the summary
        # line, so that "OK" never silently means "five of six".
        print("== Family 6: acentric factors -- SKIPPED at your request (--skip-acentric) ==")
        print("\nFive of six table families verified against teqp; acentric factors SKIPPED. OK")
        return

    try:
        import teqp  # noqa: PLC0415

        # CAPABILITY probe, not a version or attribute probe. An importable
        # teqp is not necessarily a USABLE one, and the two obvious cheap
        # checks are both wrong here:
        #
        #   * a version compare would have to be maintained by hand;
        #   * `hasattr(teqp, "make_model")` does NOT discriminate -- teqp
        #     0.15.3, the exact version that has cost this project several
        #     round trips, HAS make_model and fails only later with
        #     "Unknown kind:GERG2004resid".
        #
        # So ask teqp the question that actually matters: build the smallest
        # possible GERG model. Anything at all going wrong there means this
        # teqp cannot do GERG, whatever its version string says. This is a
        # construction of a trivial model, not a numerical evaluation, so it
        # cannot mask a real disagreement in the values -- those are computed
        # after this point and are never caught here.
        try:
            teqp.make_model({"kind": "GERG2004resid", "model": {"names": ["methane"]}})
        except Exception as probe_exc:  # noqa: BLE001 -- teqp raises a bare RuntimeError for an unknown kind
            raise ImportError(
                f"teqp imported from {teqp.__file__!r} cannot build a GERG model "
                f"(teqp.make_model(kind='GERG2004resid') raised {type(probe_exc).__name__}: {probe_exc}). "
                "teqp before 0.18 predates the GERG factory registration."
            ) from probe_exc
    except ImportError as exc:
        # NOT a skip.  Family 6 is the only check that cannot be done from
        # teqp's source alone, and letting a missing module quietly drop it
        # would mean this script prints OK while the one derived table in the
        # backend goes unverified.
        fail(
            f"family 6 (acentric factors) needs a working teqp python module, and this one is unusable ({exc}). "
            "Install it into a venv -- python3 -m venv /tmp/gergenv && "
            '/tmp/gergenv/bin/pip install "teqp>=0.18" numpy -- and rerun with that interpreter, '
            "or pass --skip-acentric to verify only the five transcribed families."
        )
    check_acentric()

    print("\nAll six table families verified against teqp. OK")


if __name__ == "__main__":
    main()
