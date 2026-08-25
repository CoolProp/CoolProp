"""Convert RES transport property parameter CSV files to res_transport_parameters.json.

Source text files live under dev/RES_params/ (copied from kktoolbox/property-calculation;
see dev/RES_params/source_documentation.txt for the literature source of each file).
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
        "HEOS": {"citation": "Martinek et al. 2025", "doi": "10.1021/acs.jced.4c00451"},
        "PR": {"citation": "Yang et al. 2025", "doi": "10.1021/acsomega.4c10815"},
        "SRK": {"citation": "Yang et al. 2025", "doi": "10.1021/acsomega.4c10815"},
    },
    "conductivity": {
        "dilute": {"citation": "Li et al. 2024", "doi": "10.1021/acs.iecr.4c02946"},
        "HEOS": {"citation": "Li et al. 2024", "doi": "10.1021/acs.iecr.4c02946"},
        "PR": {"citation": "Yang et al. 2025", "doi": "10.1021/acsomega.4c10815"},
        "SRK": {"citation": "Yang et al. 2025", "doi": "10.1021/acsomega.4c10815"},
        "critical_enhancement": {"citation": "Li et al. 2024", "doi": "10.1021/acs.iecr.4c02946"},
    },
}

# Fluids whose REFPROP-fitted RES coefficients are not trustworthy on the HEOS backend,
# because HEOS's smolar_residual differs enough from REFPROP's to move the transport
# property. Their "HEOS" entry is omitted below so the C++ backend raises a clear exception
# instead of silently returning an inaccurate value.
#
# Two measurements decide an exclusion, and BOTH matter:
#
#   1. deviation of the transport property from the published reference (Martinek vis_res /
#      Li TC_RES), via dev/compare_HEOS_vs_REFPROP_RES.py at a 1% threshold;
#   2. deviation of s_res itself over a range of states.
#
# Criterion (1) alone is NOT sufficient, because there is only ONE published sample point
# per fluid and some of those points are almost entirely dilute gas -- where the transport
# property barely depends on s_res, so a large s_res error cannot show up. Measured residual
# share at the sample point: R161 4.3% (viscosity) and -0.2% (conductivity), VINYLCHLORIDE
# 0.3%, R13 0.9%, R1234YF 0.2%. Those samples are blind to the error, and each of those
# fluids has a systematic s_res deviation elsewhere on the (T,p) surface (R161 12%,
# VINYLCHLORIDE 32%, R13 17%, R1234YF 6%), so they stay excluded despite a <1% transport
# deviation at their sample point.
#
# By contrast the sample points for D5 (99.9%), R1233ZDE (98.4%) and R1224YDZ (96.1%) are
# residual-dominated, so their transport deviations are meaningful on their own.
#
# Two cases need care, because HEOS and REFPROP can disagree about the *phase* at a sample
# point that straddles their (nearly identical) saturation curves. That inflates the
# deviation without saying anything about the parameters:
#   BENZENE  -- deviation collapses from -41.6% to -0.008% once the phase is imposed, so its
#               parameters are fine and it is NOT excluded.
#   R41      -- deviation collapses from +266.8% to -1.12%, which still exceeds 1%, so it
#               stays excluded on the strength of that genuine residual error.
#
# Re-run the comparison script (against a REFPROP-enabled CoolProp build, `pip install .`)
# to refresh these lists if the fluid library, the EOS, or the parameter fits change.
HEOS_VISCOSITY_EXCLUDE = {"D5", "R1224YDZ", "R1233ZDE", "R13", "R161", "R41"}
HEOS_CONDUCTIVITY_EXCLUDE = {"D5", "R1233ZDE", "R1234YF", "R13", "R161", "VINYLCHLORIDE"}


@functools.cache
def coolprop_has_fluid(name: str) -> bool:
    """True if any CoolProp backend (HEOS, PR, SRK) can construct this fluid.

    Fluids the local CoolProp build has no data for at all (only fitted against REFPROP's
    much larger fluid database) can never be looked up by FluidLibrary::load_RES_transport_parameters,
    so keeping their entry around in the JSON is dead weight.
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
    # EOS label in source files -> JSON key
    eos_map = {
        "REFPROP": "HEOS",  # REFPROP-fitted params are used as HEOS approximation
        "PR": "PR",
        "SRK": "SRK",
    }

    out: dict = {"viscosity": {}, "conductivity": {}, "sources": SOURCES}

    # ── viscosity ──────────────────────────────────────────────────────────────
    # load dilute-gas params once (same for all EOS)
    dilute_v, _ = load_vis_params(params_dir, "REFPROP")

    # load per-EOS residual params
    res_by_eos: dict[str, pd.DataFrame] = {}
    for csv_eos, json_eos in eos_map.items():
        _, res_df = load_vis_params(params_dir, csv_eos)
        res_by_eos[json_eos] = res_df

    n_skipped_unavailable_vis = 0
    n_excluded_heos_vis = 0
    for idx, drow in dilute_v.iterrows():
        fluid = str(drow["Material"]).strip().upper()
        if not coolprop_has_fluid(fluid):
            n_skipped_unavailable_vis += 1
            continue
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
            entry[json_eos] = build_vis_entry(res_row)
        if fluid in HEOS_VISCOSITY_EXCLUDE:
            entry.pop("HEOS", None)
            n_excluded_heos_vis += 1
        out["viscosity"][fluid] = entry

    # ── thermal conductivity ───────────────────────────────────────────────────
    dilute_tc, _, _ = load_tc_params(params_dir, "REFPROP")

    res_tc_by_eos: dict[str, pd.DataFrame] = {}
    for csv_eos, json_eos in eos_map.items():
        _, res_df, _ = load_tc_params(params_dir, csv_eos)
        res_tc_by_eos[json_eos] = res_df

    # critical enhancement params are EOS-independent (universal constants)
    _, _, crit_df = load_tc_params(params_dir, "REFPROP")

    n_skipped_unavailable_tc = 0
    n_excluded_heos_tc = 0
    for idx, drow in dilute_tc.iterrows():
        fluid = str(drow["Material"]).strip().upper()
        if not coolprop_has_fluid(fluid):
            n_skipped_unavailable_tc += 1
            continue
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
            entry[json_eos] = build_tc_entry(res_row)

        crow = crit_df.iloc[idx]
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


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description="Build dev/res_transport_parameters.json from dev/RES_params/.")
    ap.add_argument("params_dir", nargs="?", default=DEFAULT_PARAMS_DIR,
                    help="directory holding the RES source tables (default: dev/RES_params)")
    ap.add_argument("--keep-all-heos", action="store_true",
                    help="ignore HEOS_*_EXCLUDE and emit every HEOS entry. Use this to re-measure the "
                         "per-fluid deviation with dev/compare_HEOS_vs_REFPROP_RES.py before deciding "
                         "which fluids actually warrant exclusion; do not ship the result.")
    args = ap.parse_args()
    if args.keep_all_heos:
        HEOS_VISCOSITY_EXCLUDE.clear()
        HEOS_CONDUCTIVITY_EXCLUDE.clear()
        print("WARNING: --keep-all-heos set; emitting unvetted HEOS entries (measurement only).")
    main(args.params_dir)
