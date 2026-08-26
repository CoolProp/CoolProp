"""Convert RES transport property parameter CSV files to res_transport_parameters.json.

Source text files live under dev/RES_params/, taken from the supporting information of the
papers cited in dev/RES_params/source_documentation.txt, which records the origin of each file.
Run from the CoolProp root:
    python dev/convert_RES_csv_to_json.py [RES_params_dir]

RES_params_dir defaults to dev/RES_params next to this script.
"""

import functools
import json
import os
import sys

import pandas as pd

try:
    import CoolProp.CoolProp as CP
except ImportError:
    sys.exit(
        "CoolProp (Python package) is required to filter out fluids that no CoolProp "
        "backend can construct. Activate the environment that has it installed."
    )

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
DEFAULT_PARAMS_DIR = os.path.join(SCRIPT_DIR, "RES_params")

# Literature sources per parameter group, transcribed from
# dev/RES_params/source_documentation.txt.
SOURCES = {
    "viscosity": {
        "dilute": {"citation": "Martinek et al. 2025", "doi": "10.1021/acs.jced.4c00451"},
        "REFPROP": {"citation": "Martinek et al. 2025", "doi": "10.1021/acs.jced.4c00451"},
        "HEOS": {"citation": "Martinek et al. 2025", "doi": "10.1021/acs.jced.4c00451"},
        "PR": {"citation": "Yang et al. 2025", "doi": "10.1021/acsomega.4c10815"},
        "SRK": {"citation": "Yang et al. 2025", "doi": "10.1021/acsomega.4c10815"},
    },
    "conductivity": {
        "dilute": {"citation": "Li et al. 2024", "doi": "10.1021/acs.iecr.4c02946"},
        "REFPROP": {"citation": "Li et al. 2024", "doi": "10.1021/acs.iecr.4c02946"},
        "HEOS": {"citation": "Li et al. 2024", "doi": "10.1021/acs.iecr.4c02946"},
        "PR": {"citation": "Yang et al. 2025", "doi": "10.1021/acsomega.4c10815"},
        "SRK": {"citation": "Yang et al. 2025", "doi": "10.1021/acsomega.4c10815"},
        "critical_enhancement": {"citation": "Li et al. 2024", "doi": "10.1021/acs.iecr.4c02946"},
    },
}

# Fluids whose REFPROP-fitted RES coefficients are not trustworthy on the HEOS backend, because
# HEOS's residual entropy differs enough from REFPROP's to move the transport property.
# Their "HEOS" entry is omitted below so the C++ backend raises a clear exception instead of
# silently returning an inaccurate value.  The "REFPROP" entry is unaffected -- there the
# coefficients are exact, since that is the EOS they were regressed against.
#
# DERIVED FROM MEASUREMENT, not from the published sample points.  Reproduce with:
#
#     python dev/convert_RES_csv_to_json.py --keep-all-heos   # so excluded fluids can be measured
#     python dev/generate_headers.py && <rebuild>
#     python dev/RES_grid_build.py                            # 3828 points over 147 fluids
#     ./build_tests/Release/CatchTestRunner.exe "[RES_grid]"  # evaluate both backends
#     python dev/RES_grid_report.py
#
# The grid compares the SAME implementation on both backends at the same (T, p), with both pinned
# to the same dilute-gas and enhancement-viscosity source, so the equation of state is the only
# thing that differs.  Three earlier attempts at this list were wrong in instructive ways:
#
#   - One published sample point per fluid is not enough.  Several of those points are almost
#     entirely dilute gas, where the property barely depends on s_res at all, so a large s_res
#     error cannot show up.  The grid therefore judges each fluid only at states where the
#     residual term supplies at least 20% of the property.
#   - Leaving the backends on their default source policies measured that policy, not the
#     parameters: REFPROP takes the dilute term and the enhancement viscosity from its own
#     transport model and HEOS cannot, which alone moved ~40 fluids past the threshold.
#   - A grid point at exactly rho/rhoc = 2 sits on Li's critical-enhancement gate.  The two
#     backends have slightly different rhoc, so such a point has the enhancement on for one and
#     off for the other; that step alone inflated the conductivity list from 12 fluids to 30.
#
# The criterion is a >1% transport deviation at residual-dominated states away from the critical
# enhancement, OR a >5% deviation in s_res itself at those states.  The second clause matters:
# MD4M's conductivity agrees to 0.27% while its s_res is off by 27%, which means the sampled
# states happen to be insensitive rather than that the parameters transfer.  Both clauses select
# the SAME 14 fluids for viscosity and for conductivity -- unsurprising in hindsight, since s_res
# is the one input both models share.
#
# Median deviation across all 3243 comparable grid points is 0.00006% (viscosity) and 0.00005%
# (conductivity): for the other 133 fluids the parameters transfer essentially exactly.
#
# NOTE BENZENE.  An earlier round kept it, on the strength of a single sample point that agreed to
# -0.008% once the phase was imposed.  On the grid it deviates by 1.68% at residual-dominated
# states and its s_res is off by 5.6%, so that one point was simply not diagnostic.
HEOS_TRANSFER_EXCLUDE = {
    "BENZENE", "D5", "HEPTANE", "MD3M", "MD4M", "R1123", "R1224YDZ",
    "R1233ZDE", "R1234YF", "R1243ZF", "R13", "R161", "R41", "VINYLCHLORIDE",
}
# Kept as separate names so a future property-specific refit can diverge again.
HEOS_VISCOSITY_EXCLUDE = set(HEOS_TRANSFER_EXCLUDE)
HEOS_CONDUCTIVITY_EXCLUDE = set(HEOS_TRANSFER_EXCLUDE)


@functools.cache
def coolprop_has_fluid(name: str) -> bool:
    """True if any native CoolProp backend (HEOS, PR, SRK) can construct this fluid.

    This gates the *native-backend keys*, not the fluid entry itself.  29 of the 151 fitted
    fluids exist only in REFPROP's much larger database; they still get a "REFPROP" block
    (the REFPROP backend can load them), but emitting HEOS/PR/SRK blocks for them would ship
    parameters no backend can ever reach.
    """
    for eos in ("HEOS", "PR", "SRK"):
        try:
            CP.AbstractState(eos, name)
            return True
        except Exception:
            pass
    return False


def load_vis_params(data_dir: str, eos: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    dilute = pd.read_csv(
        os.path.join(data_dir, "viscosity", "dilute_gas_V_Parameters.txt"),
        skipinitialspace=True,
        header=1,
        delimiter=";",
    )
    res = pd.read_csv(
        os.path.join(data_dir, "viscosity", f"RES_V_Parameters_{eos}.txt"),
        sep=r"\s+",
        header=0,
        skipinitialspace=True,
    )
    return dilute, res


def load_tc_params(data_dir: str, eos: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    dilute = pd.read_csv(
        os.path.join(data_dir, "thermal_conductivity", "dilute_gas_TC_Parameters.txt"),
        skipinitialspace=True,
        header=1,
        delimiter=";",
    )
    res = pd.read_csv(
        os.path.join(data_dir, "thermal_conductivity", f"RES_TC_Parameters_{eos}.txt"),
        sep=r"\s+",
        header=0,
        skipinitialspace=True,
    )
    crit = pd.read_csv(
        os.path.join(data_dir, "thermal_conductivity", "critical_TC_Parameters.txt"),
        sep=r"\s+",
        header=0,
        skipinitialspace=True,
    )
    return dilute, res, crit


def check_same_fluid(fluid: str, row, table: str) -> None:
    """Fail loudly when a positionally-paired row belongs to a different fluid.

    The dilute tables are semicolon-delimited and the residual / critical tables are
    whitespace-delimited, so they are parsed independently and nothing but row ORDER ties a
    coefficient set to its fluid.  Reorder or add a row in any one file and this script would
    otherwise emit one fluid's coefficients under another fluid's name, silently, into the
    shipped header.  Every table carries a `Material` column, so the check is cheap.
    """
    other = str(row["Material"]).strip().upper()
    if other != fluid:
        raise SystemExit(f"{table}: row order mismatch -- dilute table has {fluid!r}, this table has {other!r}")


def build_vis_entry(res_row) -> dict:
    ind_fit = int(res_row["ind_fit"]) == 1
    if ind_fit:
        n_res = [float(res_row["n1_ind"]), float(res_row["n2_ind"]), float(res_row["n3_ind"])]
        xita = 1.0
    else:
        n_res = [float(res_row["n1_glb"]), float(res_row["n2_glb"]), float(res_row["n3_glb"])]
        xita = float(res_row["xita"])
    return {
        "n": n_res,
        "xita": xita,
        "group": int(res_row["GroupN"]),
        "ind_fit": ind_fit,
    }


def build_tc_entry(res_row) -> dict:
    ind_fit = int(res_row["ind_fit"]) == 1
    if ind_fit:
        n_res = [
            float(res_row["n1_ind"]),
            float(res_row["n2_ind"]),
            float(res_row["n3_ind"]),
            float(res_row["n4_ind"]),
        ]
        xita = 1.0
    else:
        n_res = [
            float(res_row["n1_glb"]),
            float(res_row["n2_glb"]),
            float(res_row["n3_glb"]),
            float(res_row["n4_glb"]),
        ]
        xita = float(res_row["xita"])
    return {
        "n": n_res,
        "xita": xita,
        "group": int(res_row["GroupN"]),
        "ind_fit": ind_fit,
    }


def main(params_dir: str) -> None:
    # EOS label in the source files -> JSON key(s).
    #
    # The REFPROP-fitted coefficients are emitted TWICE, under two keys that mean different
    # things: as "REFPROP" they are exact (that is the EOS they were regressed against), and
    # as "HEOS" they are an approximation whose validity is decided by HEOS_*_EXCLUDE below.
    # Shipping the duplicate costs ~370 kB of header and keeps the two meanings separable.
    eos_map = {
        "REFPROP": ["REFPROP", "HEOS"],
        "PR": ["PR"],
        "SRK": ["SRK"],
    }

    out: dict = {"viscosity": {}, "conductivity": {}, "sources": SOURCES}

    # ── viscosity ──────────────────────────────────────────────────────────────
    # load dilute-gas params once (same for all EOS)
    dilute_v, _ = load_vis_params(params_dir, "REFPROP")

    # load per-EOS residual params
    res_by_eos: dict[str, pd.DataFrame] = {}
    for csv_eos, json_keys in eos_map.items():
        _, res_df = load_vis_params(params_dir, csv_eos)
        for json_eos in json_keys:
            res_by_eos[json_eos] = res_df

    n_refprop_only_vis = 0
    n_excluded_heos_vis = 0
    for idx, (_label, drow) in enumerate(dilute_v.iterrows()):
        fluid = str(drow["Material"]).strip().upper()
        n_dilute = [
            float(drow["n0"]),
            float(drow["n1"]),
            float(drow["n2"]),
            float(drow["n3"]),
            float(drow["n4"]),
        ]
        entry: dict = {"dilute": {"n": n_dilute}}
        for json_eos, res_df in res_by_eos.items():
            res_row = res_df.iloc[idx]
            check_same_fluid(fluid, res_row, f"RES_V_Parameters_{json_eos}")
            entry[json_eos] = build_vis_entry(res_row)
        # Filters apply PER KEY, not per fluid: "dilute" and "REFPROP" are always valid.
        if not coolprop_has_fluid(fluid):
            for dead_key in ("HEOS", "PR", "SRK"):
                entry.pop(dead_key, None)
            n_refprop_only_vis += 1
        if fluid in HEOS_VISCOSITY_EXCLUDE:
            entry.pop("HEOS", None)
            n_excluded_heos_vis += 1
        out["viscosity"][fluid] = entry

    # ── thermal conductivity ───────────────────────────────────────────────────
    dilute_tc, _, _ = load_tc_params(params_dir, "REFPROP")

    res_tc_by_eos: dict[str, pd.DataFrame] = {}
    for csv_eos, json_keys in eos_map.items():
        _, res_df, _ = load_tc_params(params_dir, csv_eos)
        for json_eos in json_keys:
            res_tc_by_eos[json_eos] = res_df

    # critical enhancement params are EOS-independent (universal constants)
    _, _, crit_df = load_tc_params(params_dir, "REFPROP")

    n_refprop_only_tc = 0
    n_excluded_heos_tc = 0
    for idx, (_label, drow) in enumerate(dilute_tc.iterrows()):
        fluid = str(drow["Material"]).strip().upper()
        n_dilute = [
            float(drow["n0"]),
            float(drow["n1"]),
            float(drow["n2"]),
            float(drow["n3"]),
            float(drow["n4"]),
        ]
        entry: dict = {"dilute": {"n": n_dilute}}
        for json_eos, res_df in res_tc_by_eos.items():
            res_row = res_df.iloc[idx]
            check_same_fluid(fluid, res_row, f"RES_TC_Parameters_{json_eos}")
            entry[json_eos] = build_tc_entry(res_row)

        crow = crit_df.iloc[idx]
        check_same_fluid(fluid, crow, "critical_TC_Parameters")
        q_D_inv = float(crow["qDinv"])
        # An all-zero row in critical_TC_Parameters.txt means "no critical-enhancement parameters
        # fitted for this fluid" (e.g. D2O, HELIUM, ORTHOHYD), not "zero enhancement".  Li 2024
        # gates the enhancement on Tref > 0, so omit the block entirely rather than shipping zeros
        # that would divide by Tref/Gamma downstream.
        if float(crow["Tref"]) > 0 and float(crow["Gamma"]) > 0 and float(crow["xi0"]) > 0 and q_D_inv > 0:
            entry["critical_enhancement"] = {
                "R_D": float(crow["R_D"]),
                "gamma_uni": float(crow["gamma_uni"]),
                "Gamma": float(crow["Gamma"]),
                "phi0": float(crow["xi0"]),
                "t_ref": float(crow["Tref"]),
                "q_D": 1.0 / q_D_inv,
            }
        # Filters apply PER KEY, not per fluid (see the viscosity loop above).
        if not coolprop_has_fluid(fluid):
            for dead_key in ("HEOS", "PR", "SRK"):
                entry.pop(dead_key, None)
            n_refprop_only_tc += 1
        if fluid in HEOS_CONDUCTIVITY_EXCLUDE:
            entry.pop("HEOS", None)
            n_excluded_heos_tc += 1
        out["conductivity"][fluid] = entry

    out_path = os.path.join(SCRIPT_DIR, "res_transport_parameters.json")
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    print(f"Wrote {out_path}")
    print(f"  viscosity entries : {len(out['viscosity'])}"
          f" ({n_refprop_only_vis} REFPROP-only: no native CoolProp backend can build them;"
          f" {n_excluded_heos_vis} without a HEOS entry: HEOS/REFPROP s_res mismatch)")
    print(f"  conductivity entries : {len(out['conductivity'])}"
          f" ({n_refprop_only_tc} REFPROP-only: no native CoolProp backend can build them;"
          f" {n_excluded_heos_tc} without a HEOS entry: HEOS/REFPROP s_res mismatch)")


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description="Build dev/res_transport_parameters.json from dev/RES_params/.")
    ap.add_argument("params_dir", nargs="?", default=DEFAULT_PARAMS_DIR,
                    help="directory holding the RES source tables (default: dev/RES_params)")
    ap.add_argument("--keep-all-heos", action="store_true",
                    help="ignore HEOS_*_EXCLUDE and emit every HEOS entry. Required before re-deriving "
                         "the exclusion lists, since an excluded fluid cannot be evaluated on HEOS at "
                         "all; see the comment on HEOS_TRANSFER_EXCLUDE for the full procedure. "
                         "Do not ship the result.")
    args = ap.parse_args()
    if args.keep_all_heos:
        HEOS_VISCOSITY_EXCLUDE.clear()
        HEOS_CONDUCTIVITY_EXCLUDE.clear()
        print("WARNING: --keep-all-heos set; emitting unvetted HEOS entries (measurement only).")
    main(args.params_dir)
