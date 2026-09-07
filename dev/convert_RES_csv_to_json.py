"""Build dev/res_transport_parameters.json from the published RES parameter tables.

Every input is read straight from the vendored supporting information under dev/RES_reference/,
which is the authors' own copy in each case -- there is no second, re-headed copy of any table in
this repository, deliberately.  See dev/RES_reference/README.md for provenance and licensing.

    property      block                     source file (under dev/RES_reference/)
    ------------  ------------------------  --------------------------------------------------
    viscosity     dilute                    martinek2025_viscosity/Dilute_gas_viscosity.txt
    viscosity     HEOS  (REFPROP-fitted)    martinek2025_viscosity/RES_Parameter.txt
    viscosity     PR, SRK                   yang2025/RES_V_Parameters_{PR,SRK}.txt
    conductivity  dilute                    li2024_conductivity/Dilute_gas_TC.txt
    conductivity  HEOS  (REFPROP-fitted)    li2024_conductivity/RES_Parameter.txt
    conductivity  PR, SRK                   yang2025/RES_TC_Parameters_{PR,SRK}.txt
    conductivity  critical_enhancement      li2024_conductivity/Fluid_Constants.txt

Run from the CoolProp root:

    python dev/convert_RES_csv_to_json.py
"""

import functools
import glob
import json
import os

import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
REFERENCE_DIR = os.path.join(SCRIPT_DIR, "RES_reference")
VIS_REF_DIR = os.path.join(REFERENCE_DIR, "martinek2025_viscosity")
TC_REF_DIR = os.path.join(REFERENCE_DIR, "li2024_conductivity")
YANG_REF_DIR = os.path.join(REFERENCE_DIR, "yang2025")
FLUIDS_DIR = os.path.join(SCRIPT_DIR, "fluids")
CUBIC_FLUIDS_JSON = os.path.join(SCRIPT_DIR, "cubics", "all_cubic_fluids.json")

# Literature sources per parameter group, embedded in the generated JSON so a consumer of the
# shipped table can cite it without reaching back into this repository.  Matches the file table
# in the module docstring above.
SOURCES = {
    "viscosity": {
        "dilute": {"citation": "Martinek et al. 2025", "doi": "10.1021/acs.jced.4c00451"},
        "HEOS": {"citation": "Martinek et al. 2025", "doi": "10.1021/acs.jced.4c00451"},
        "PR": {"citation": "Yang 2025", "doi": "10.1021/acsomega.4c10815"},
        "SRK": {"citation": "Yang 2025", "doi": "10.1021/acsomega.4c10815"},
    },
    "conductivity": {
        "dilute": {"citation": "Li et al. 2024", "doi": "10.1021/acs.iecr.4c02946"},
        "HEOS": {"citation": "Li et al. 2024", "doi": "10.1021/acs.iecr.4c02946"},
        "PR": {"citation": "Yang 2025", "doi": "10.1021/acsomega.4c10815"},
        "SRK": {"citation": "Yang 2025", "doi": "10.1021/acsomega.4c10815"},
        "critical_enhancement": {"citation": "Li et al. 2024", "doi": "10.1021/acs.iecr.4c02946"},
    },
}

# Fluids whose REFPROP-fitted RES coefficients are not trustworthy on the HEOS backend, because
# HEOS's residual entropy differs enough from REFPROP's to move the transport property.  Their
# "HEOS" entry is omitted below so the C++ backend raises a clear exception instead of silently
# returning an inaccurate value.  The PR / SRK entries are unaffected -- those coefficients were
# regressed against those equations of state directly.
#
# DERIVED FROM MEASUREMENT, not from the published sample points.  Reproduce with:
#
#     python dev/convert_RES_csv_to_json.py --keep-all-heos   # so excluded fluids stay measurable
#     python dev/RES_grid_build.py                            # reduced-coordinate grid
#     python dev/RES_reference_run.py --eos REFPROP
#     python dev/RES_reference_run.py --eos HEOS
#     python dev/RES_grid_report.py
#
# Both columns come from the AUTHORS' OWN reference code, vendored under dev/RES_reference/ and
# run over the same grid with the same model choices and only the equation of state swapped.  One
# implementation, two equations of state, so what is left is the equation of state and nothing
# else -- and, usefully, no CoolProp RES implementation is involved, so this list can be derived
# and shipped before any C++ model exists.
#
# Three earlier attempts at this list were wrong in instructive ways, and the harness is built
# around avoiding all three:
#
#   - One published sample point per fluid is not enough.  Several of those points are almost
#     entirely dilute gas, where the property barely depends on s_res at all, so a large s_res
#     error cannot show up.  The grid therefore judges each fluid only at states where the
#     residual term supplies at least 20% of the property.
#   - Leaving the two sides on their native dilute-gas and enhancement-viscosity sources measured
#     that source choice, not the parameters.  RES_reference_run.py pins both.
#   - A grid point at exactly rho/rhoc = 2 sits on Li's critical-enhancement gate.  The two
#     equations of state have slightly different rhoc, so such a point has the enhancement on for
#     one and off for the other; that step alone inflated an earlier conductivity list from 12
#     fluids to 30.
#
# The criterion is a >1% transport deviation at residual-dominated states away from the critical
# enhancement, OR a >5% deviation in s_res itself at those states.  The second clause matters: a
# fluid can agree on the transport property at every sampled state while its s_res is badly off,
# which means the states happen to be insensitive rather than that the parameters transfer.
#
# Measured over 3570 grid points on 122 fluids.  For the 108 fluids NOT listed below the
# parameters transfer essentially exactly: the median deviation across all comparable points is
# 0.00004 % (viscosity) and 0.00007 % (conductivity).
#
# NARROW vs BROAD is a real choice, not a formality.  Judging everywhere instead of only where the
# critical enhancement is inactive would exclude 28 fluids for conductivity rather than 12 -- but
# the extra 16 move because the enhancement consumes each equation of state's OWN critical point
# and derivatives, which says nothing about the RES coefficients.  Viscosity has no enhancement
# and its two criteria select the same 14 fluids either way.
#
# One caveat the report prints and this list does not encode: a few fluids per property are
# judged on 5 or fewer points, because not every sampled state puts enough of the property in the
# residual term to be informative about it.  Those are kept rather than withheld -- every one
# agrees to better than 0.01%, four orders below the threshold, so the thin evidence is thin
# agreement rather than a near miss.  Re-run dev/RES_grid_report.py for the per-fluid counts
# before trusting any single entry here.

HEOS_TRANSFER_EXCLUDE = {
    "BENZENE", "D5", "HEPTANE", "MD3M", "MD4M", "R1123", "R1224YDZ",
    "R1233ZDE", "R1234YF", "R1243ZF", "R13", "R161", "R41", "VINYLCHLORIDE",
}
# Kept as separate names so a future property-specific refit can diverge again.
HEOS_VISCOSITY_EXCLUDE = set(HEOS_TRANSFER_EXCLUDE)
HEOS_CONDUCTIVITY_EXCLUDE = set(HEOS_TRANSFER_EXCLUDE)


def check_table_alignment(reference, tables: dict) -> None:
    """Verify every table lists the same fluids in the same order as the dilute table.

    The dilute tables are semicolon-delimited and the residual / critical tables are
    whitespace-delimited, so they are parsed independently and nothing but row ORDER ties a
    coefficient set to its fluid.  Reorder or add a row in any one file and this script would
    otherwise emit one fluid's coefficients under another fluid's name, silently, into the
    shipped header.

    Checked here for the WHOLE column rather than per emitted row: rows are dropped further down
    for fluids no CoolProp backend can build, and a misalignment that starts inside a dropped
    run would go unseen by a per-row check while still corrupting every row after it.
    """

    def names(df):
        return [str(v).strip().upper() for v in df["Material"]]

    want = names(reference)
    for label, df in tables.items():
        got = names(df)
        if got == want:
            continue
        if len(got) != len(want):
            raise SystemExit(f"{label}: has {len(got)} rows, dilute table has {len(want)}")
        first = next(i for i, (a, b) in enumerate(zip(got, want)) if a != b)
        raise SystemExit(
            f"{label}: row order mismatch at row {first} -- dilute table has {want[first]!r}, "
            f"this table has {got[first]!r}"
        )


@functools.cache
def backend_fluid_names():
    """{backend: frozenset(every name that backend answers to, uppercased)}.

    Read from the REPOSITORY's own fluid data rather than by constructing states through an
    installed CoolProp: the JSON produced here is compiled into THIS tree's binary, so the fluid
    set has to be this tree's.  Resolving it through `import CoolProp` instead makes the shipped
    table depend on whichever build happens to be on sys.path -- which is how an earlier run
    emitted 115 fluids where the tree supports 122.

    Mirrors the lookup in JSONFluidLibrary::load_RES_transport_parameters, which tries the
    fluid's name, its REFPROP name and each alias, in both the original and uppercased form.
    Everything is uppercased on both sides here, since the RES tables are uppercased already.
    """
    heos = set()
    cubic = set()

    def add(into, value):
        if isinstance(value, str) and value:
            into.add(value.strip().upper())

    # HEOS: dev/fluids/*.json is the source dev/all_fluids.cbor is built from.  `.json_disabled`
    # files are deliberately excluded -- those fluids do not reach the library either.
    for path in glob.glob(os.path.join(FLUIDS_DIR, "*.json")):
        with open(path, encoding="utf-8") as fh:
            info = json.load(fh).get("INFO", {})
        add(heos, info.get("NAME"))
        add(heos, info.get("REFPROP_NAME"))
        for alias in info.get("ALIASES", []):
            add(heos, alias)

    # PR and SRK share one fluid list.
    with open(CUBIC_FLUIDS_JSON, encoding="utf-8") as fh:
        for entry in json.load(fh):
            add(cubic, entry.get("name"))
            for alias in entry.get("aliases", []):
                add(cubic, alias)

    # A renamed directory or partial checkout would leave these empty, and every fluid would then
    # be judged unavailable -- producing an EMPTY table and exiting 0.
    if not heos or not cubic:
        raise SystemExit(
            f"no fluid names found (HEOS: {len(heos)}, cubic: {len(cubic)}); expected dev/fluids/*.json "
            "and dev/cubics/all_cubic_fluids.json to be present"
        )
    return {"HEOS": frozenset(heos), "PR": frozenset(cubic), "SRK": frozenset(cubic)}


def backend_has_fluid(eos_key: str, name: str) -> bool:
    """True if THAT backend knows this fluid.

    Gated per EOS key, not per fluid.  10 of the fitted fluids are HEOS-only -- C4F10, C5F12,
    C6F14, CHLORINE, R1123, R1224YDZ, R1234ZEZ, R1243ZF, R1336MZZZ and VINYLCHLORIDE are in
    dev/fluids/ but not in dev/cubics/all_cubic_fluids.json -- so this withholds their PR and SRK
    blocks, which no cubic backend could ever reach.  The gate runs in the other direction too;
    it just happens that no fitted fluid is cubic-only today.
    """
    return name.strip().upper() in backend_fluid_names()[eos_key]


def coolprop_has_fluid(name: str) -> bool:
    """True if ANY native CoolProp backend knows this fluid.

    Gates the whole entry: the RES fits cover fluids from REFPROP's much larger database, and an
    entry no backend at all can reach is pure dead weight.
    """
    return any(backend_has_fluid(eos, name) for eos in ("HEOS", "PR", "SRK"))


def load_vis_params(eos: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    dilute = pd.read_csv(
        os.path.join(VIS_REF_DIR, "Dilute_gas_viscosity.txt"),
        skipinitialspace=True,
        header=1,
        delimiter=";",
    )
    res_path = (
        os.path.join(VIS_REF_DIR, "RES_Parameter.txt")
        if eos == "REFPROP"
        else os.path.join(YANG_REF_DIR, f"RES_V_Parameters_{eos}.txt")
    )
    res = pd.read_csv(res_path, sep=r"\s+", header=0, skipinitialspace=True)
    return dilute, res


def load_tc_params(eos: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    dilute = pd.read_csv(
        os.path.join(TC_REF_DIR, "Dilute_gas_TC.txt"),
        skipinitialspace=True,
        header=1,
        delimiter=";",
    )
    res_path = (
        os.path.join(TC_REF_DIR, "RES_Parameter.txt")
        if eos == "REFPROP"
        else os.path.join(YANG_REF_DIR, f"RES_TC_Parameters_{eos}.txt")
    )
    res = pd.read_csv(res_path, sep=r"\s+", header=0, skipinitialspace=True)
    crit = pd.read_csv(os.path.join(TC_REF_DIR, "Fluid_Constants.txt"), sep=r"\s+", header=0, skipinitialspace=True)
    return dilute, res, crit


def vis_table_name(json_eos: str) -> str:
    """The file a viscosity residual table was actually read from, for error messages."""
    return "martinek2025_viscosity/RES_Parameter.txt" if json_eos == "HEOS" else f"yang2025/RES_V_Parameters_{json_eos}.txt"


def tc_table_name(json_eos: str) -> str:
    """The file a conductivity residual table was actually read from, for error messages."""
    return "li2024_conductivity/RES_Parameter.txt" if json_eos == "HEOS" else f"yang2025/RES_TC_Parameters_{json_eos}.txt"


def dilute_ascending(row) -> list:
    """The five dilute-gas coefficients in ASCENDING powers of T, as CoolProp evaluates them.

    Read n4..n0, not n0..n4.  Both papers' dilute tables head the T^4 column `n0` and the
    constant column `n4`, i.e. DESCENDING, and their reference code multiplies index 0 by T^4 to
    match.  CoolProp's loader evaluates n[0] + T*(n[1] + ...), so the order has to be flipped
    exactly here, once.

    An earlier revision of this repository kept its own copy of these tables with the column
    names reversed to make this reindexing unnecessary.  Two copies of one table under two
    conventions is a trap -- read either one with the other's convention and the dilute viscosity
    comes out ~1e9 times too large -- so the copies are gone and the flip is explicit instead.
    Guarded by dilute_order_is_descending() below.
    """
    return [float(row["n4"]), float(row["n3"]), float(row["n2"]), float(row["n1"]), float(row["n0"])]


def check_dilute_order(df, table: str) -> None:
    """Fail if a dilute table is not in the descending order dilute_ascending() assumes.

    Compares COLUMN MAXIMA, not row by row.  Deliberately: the per-row property does not hold --
    PROPANE in the viscosity table and RE347MCC in the conductivity table both have a constant
    term around 1e-15, so |n0| >= |n4| for those rows -- and a per-row check would reject the
    real tables.  Over the whole column the separation is unambiguous (~1e1 vs ~1e-9).

    This catches a whole-table convention flip, which is the error that would otherwise pass
    silently and scale the dilute term by ~1e9.  It does NOT catch a single corrupted row or an
    n1<->n3 swap; check_table_alignment covers row identity, nothing covers the latter.
    """
    small = df["n0"].abs().max()   # claimed T^4 coefficient
    large = df["n4"].abs().max()   # claimed constant term
    if not (small < large):
        raise SystemExit(
            f"{table}: column n0 (max |{small:g}|) is not the small T^4 coefficient it should be "
            f"relative to n4 (max |{large:g}|); the table's column convention has changed and "
            "dilute_ascending() would emit the polynomial backwards."
        )


# Fluids whose individual-fit row was all zeros and therefore fell back to the group row.
# Collected rather than printed at the point of use so the report is one line per property per
# equation of state instead of one line per occurrence.
ZERO_IND_FALLBACK: list[str] = []


def _select_res_coefficients(res_row, n_terms: int, label: str) -> tuple[list[float], float, bool]:
    """Pick the individual or the group coefficient row, per the table's own ind_fit flag.

    With one correction.  The SI data files carry exactly one row whose ind_fit flag is set while
    all four of its individual coefficients are zero -- NEOPENTN, conductivity, in Li 2024 and in
    both of Yang 2025's tables.  Li's published paper lists that fluid with ind_fit = 0, so the
    flag is a transcription error in the SI and the group row is the intended fit.  Followed
    literally it deletes the residual term instead: -89% against REFPROP in the liquid.

    An all-zero individual row is therefore treated as absent.  dev/RES_reference/ carries the
    matching `zero_ind_fit` option so both sides of the comparison make the same choice, and every
    occurrence is printed rather than silently repaired -- a second one would not have the paper
    behind it.
    """
    ind_keys = ["n{}_ind".format(k) for k in range(1, n_terms + 1)]
    glb_keys = ["n{}_glb".format(k) for k in range(1, n_terms + 1)]
    ind_fit = int(res_row["ind_fit"]) == 1
    if ind_fit:
        n_res = [float(res_row[k]) for k in ind_keys]
        if any(n_res):
            return n_res, 1.0, True
        ZERO_IND_FALLBACK.append(label)
        ind_fit = False
    return [float(res_row[k]) for k in glb_keys], float(res_row["xita"]), ind_fit


def build_vis_entry(res_row, label: str = "?") -> dict:
    n_res, xita, ind_fit = _select_res_coefficients(res_row, 3, label)
    return {
        "n": n_res,
        "xita": xita,
        "group": int(res_row["GroupN"]),
        "ind_fit": ind_fit,
    }


def build_tc_entry(res_row, label: str = "?") -> dict:
    n_res, xita, ind_fit = _select_res_coefficients(res_row, 4, label)
    return {
        "n": n_res,
        "xita": xita,
        "group": int(res_row["GroupN"]),
        "ind_fit": ind_fit,
    }


def main() -> None:
    # EOS label in source files -> JSON key
    eos_map = {
        "REFPROP": "HEOS",  # REFPROP-fitted params are used as HEOS approximation
        "PR": "PR",
        "SRK": "SRK",
    }

    out: dict = {"viscosity": {}, "conductivity": {}, "sources": SOURCES}

    # ── viscosity ──────────────────────────────────────────────────────────────
    # load dilute-gas params once (same for all EOS)
    dilute_v, _ = load_vis_params("REFPROP")

    # load per-EOS residual params
    res_by_eos: dict[str, pd.DataFrame] = {}
    for csv_eos, json_eos in eos_map.items():
        _, res_df = load_vis_params(csv_eos)
        res_by_eos[json_eos] = res_df

    check_dilute_order(dilute_v, "Dilute_gas_viscosity.txt")
    check_table_alignment(dilute_v, {vis_table_name(e): d for e, d in res_by_eos.items()})

    n_skipped_unavailable_vis = 0
    skipped_vis: list = []
    n_excluded_heos_vis = 0
    for idx, (_label, drow) in enumerate(dilute_v.iterrows()):
        fluid = str(drow["Material"]).strip().upper()
        if not coolprop_has_fluid(fluid):
            n_skipped_unavailable_vis += 1
            skipped_vis.append(fluid)
            continue
        n_dilute = dilute_ascending(drow)
        entry: dict = {"dilute": {"n": n_dilute}}
        for json_eos, res_df in res_by_eos.items():
            if not backend_has_fluid(json_eos, fluid):
                continue
            entry[json_eos] = build_vis_entry(res_df.iloc[idx], "{} viscosity/{}".format(fluid, json_eos))
        if fluid in HEOS_VISCOSITY_EXCLUDE:
            entry.pop("HEOS", None)
            n_excluded_heos_vis += 1
        out["viscosity"][fluid] = entry

    # ── thermal conductivity ───────────────────────────────────────────────────
    dilute_tc, _, _ = load_tc_params("REFPROP")

    res_tc_by_eos: dict[str, pd.DataFrame] = {}
    for csv_eos, json_eos in eos_map.items():
        _, res_df, _ = load_tc_params(csv_eos)
        res_tc_by_eos[json_eos] = res_df

    # critical enhancement params are EOS-independent (universal constants)
    _, _, crit_df = load_tc_params("REFPROP")

    check_dilute_order(dilute_tc, "Dilute_gas_TC.txt")
    check_table_alignment(
        dilute_tc,
        {**{tc_table_name(e): d for e, d in res_tc_by_eos.items()},
         "Fluid_Constants.txt": crit_df},
    )

    n_skipped_unavailable_tc = 0
    skipped_tc: list = []
    n_excluded_heos_tc = 0
    for idx, (_label, drow) in enumerate(dilute_tc.iterrows()):
        fluid = str(drow["Material"]).strip().upper()
        if not coolprop_has_fluid(fluid):
            n_skipped_unavailable_tc += 1
            skipped_tc.append(fluid)
            continue
        n_dilute = dilute_ascending(drow)
        entry: dict = {"dilute": {"n": n_dilute}}
        for json_eos, res_df in res_tc_by_eos.items():
            if not backend_has_fluid(json_eos, fluid):
                continue
            entry[json_eos] = build_tc_entry(res_df.iloc[idx], "{} conductivity/{}".format(fluid, json_eos))

        crow = crit_df.iloc[idx]
        q_D_inv = float(crow["qDinv"])
        # An all-zero row in Fluid_Constants.txt means "no critical-enhancement parameters
        # fitted for this fluid" (e.g. D2O, HELIUM, ORTHOHYD), not "zero enhancement".  Li 2024
        # gates the enhancement on Tref > 0, so omit the block entirely rather than shipping zeros
        # that would divide by Tref/Gamma downstream.
        if float(crow["Tref"]) > 0 and float(crow["Gamma"]) > 0 and float(crow["xi0"]) > 0 and q_D_inv > 0:
            entry["critical_enhancement"] = {
                # Shipped as published, and DELIBERATELY NOT USED.  Li's reference code overwrites
                # this column with a flat 1.02 for every fluid (conductivity_RES.py, get_paramters)
                # and never reads the table value; CoolProp's model does the same.  The column is
                # not constant -- 43 of the 119 records here are 1.01 or 1.03 -- so anyone who
                # wires it in will diverge from both papers by ~1% on the enhancement term for a
                # third of the table.  Kept rather than normalised to 1.02 because the JSON is a
                # faithful transcription of the source table; the place to enforce 1.02 is the
                # model, where both reference implementations enforce it.
                "R_D": float(crow["R_D"]),
                "gamma_uni": float(crow["gamma_uni"]),
                "Gamma": float(crow["Gamma"]),
                "phi0": float(crow["xi0"]),
                "t_ref": float(crow["Tref"]),
                "q_D": 1.0 / q_D_inv,
            }
        if fluid in HEOS_CONDUCTIVITY_EXCLUDE:
            entry.pop("HEOS", None)
            n_excluded_heos_tc += 1
        out["conductivity"][fluid] = entry

    out_path = os.path.join(SCRIPT_DIR, "res_transport_parameters.json")
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    print(f"Wrote {out_path}")
    print(f"  viscosity entries : {len(out['viscosity'])}"
          f" ({n_skipped_unavailable_vis} fluids skipped: unavailable in CoolProp;"
          f" {n_excluded_heos_vis} kept without a HEOS entry: HEOS/REFPROP s_res mismatch)")
    print(f"  conductivity entries : {len(out['conductivity'])}"
          f" ({n_skipped_unavailable_tc} fluids skipped: unavailable in CoolProp;"
          f" {n_excluded_heos_tc} kept without a HEOS entry: HEOS/REFPROP s_res mismatch)")
    # Named, not counted, for the same reason the skip list below is: one occurrence is a known
    # defect in the published tables, a second would be news.
    if ZERO_IND_FALLBACK:
        print(f"  individual-fit row was all zeros, used the group row instead: "
              f"{', '.join(sorted(ZERO_IND_FALLBACK))}")
    # Named, not merely counted.  A count hides the case where CoolProp DOES have the fluid under
    # a name the papers do not use: R150 is 1,2-dichloroethane, shipped as `Dichloroethane` with
    # no alias spelled R150, so it is skipped here despite being fully supported.  Anything in
    # this list that looks like a fluid CoolProp knows is an alias gap worth closing.
    if skipped_vis or skipped_tc:
        both = sorted(set(skipped_vis) | set(skipped_tc))
        print(f"  skipped (no native backend knows the name): {', '.join(both)}")


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description="Build dev/res_transport_parameters.json from dev/RES_reference/.")
    ap.add_argument("--keep-all-heos", action="store_true",
                    help="ignore HEOS_*_EXCLUDE and emit every HEOS entry. Required before re-deriving "
                         "the exclusion lists, since an excluded fluid has no HEOS entry and so cannot "
                         "be measured at all; see the comment on HEOS_TRANSFER_EXCLUDE for the full "
                         "procedure. Do not ship the result.")
    args = ap.parse_args()
    if args.keep_all_heos:
        HEOS_VISCOSITY_EXCLUDE.clear()
        HEOS_CONDUCTIVITY_EXCLUDE.clear()
        print("WARNING: --keep-all-heos set; emitting unvetted HEOS entries (measurement only).")
    main()
