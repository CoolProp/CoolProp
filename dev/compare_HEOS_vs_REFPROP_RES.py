"""
Quantify how far CoolProp's HEOS RES transport properties drift from the
REFPROP-fitted reference, per fluid.

Purpose
-------
The RES parameters stored under the "HEOS" key in res_transport_parameters.json were
fitted against REFPROP's smolar_residual (see dev/RES_transport_plan.md).  The HEOS
backend uses a different multiparameter EOS, so its s_res differs, and the fitted
coefficients are therefore only approximately valid there.  This script measures the
resulting error so that fluids which drift too far can have their HEOS entry removed
from res_transport_parameters.json -- making the C++ backend raise a clear exception
instead of silently returning an inaccurate number.

What is actually being tested
-----------------------------
The HEOS column is produced by the **real C++ implementation**
(TransportRoutines::viscosity_RES / conductivity_RES) reached through the Python
wrapper via AbstractState.use_viscosity_RES(True) / use_conductivity_RES(True).
Nothing about the RES model is reimplemented here.

The reference column is the published model output shipped with the papers --
`vis_res` in Martinek 2025's Samples_pure_fluids_output.txt and `TC_RES` in Li 2024's
Table_S5_SI.txt -- both computed with REFPROP as the thermodynamic backend.  So the
comparison is exactly "same RES model, same fitted coefficients, REFPROP s_res vs HEOS
s_res".

REFPROP is still opened per fluid, but only to read `smolar_residual()`, `rhomass()` and
`phase()` as diagnostics -- these attribute a deviation to its cause rather than merely
flagging it:

  * a large `s_res_dev%` with a matching large `rho_dev%` and `phase_mismatch=True`
    means the sample point straddles the two backends' saturation curves.  That is an
    artefact of the individual state point, NOT evidence that the fluid's HEOS
    parameters are unusable -- check other state points before excluding it.
  * a `s_res_dev%` that is small while `transport_dev%` is large means the RES formula
    is amplifying the s_res error (the viscosity residual is exponential in s_plus), a
    genuine reason to exclude the fluid.

Requires CoolProp built from this working tree (`pip install .`) so that the RES API is
present, and REFPROP available for the diagnostic columns.

Usage
-----
    python dev/compare_HEOS_vs_REFPROP_RES.py [--threshold PCT]

Output
------
    dev/RES_comparison/compare_vis_pure.csv
    dev/RES_comparison/compare_tc_pure.csv
    Printed summary of fluids exceeding the threshold.
"""

import argparse
import os

import CoolProp.CoolProp as CP
import numpy as np
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
SAMPLES_DIR = os.path.join(SCRIPT_DIR, "RES_samples")
OUT_DIR = os.path.join(SCRIPT_DIR, "RES_comparison")
os.makedirs(OUT_DIR, exist_ok=True)

VIS_PURE_OUT = os.path.join(SAMPLES_DIR, "Martinek_2025_viscosity", "Samples_pure_fluids_output.txt")
TC_PURE_REF = os.path.join(SAMPLES_DIR, "Li_2024_conductivity", "Table_S5_SI.txt")

DEFAULT_THRESHOLD_PCT = 1.0


def pct(a, b):
    """Signed percentage deviation of a relative to b."""
    if b == 0 or not np.isfinite(a) or not np.isfinite(b):
        return np.nan
    return 100.0 * (a - b) / abs(b)


def load_vis_samples(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", header=0, skipinitialspace=True)
    # Material T/K p/MPa den/kg/m3 s_resi/JmolK vis_exp/uPas vis_res vis_REF
    df.columns = ["material", "T", "p_MPa", "den", "s_res_ref", "vis_exp", "vis_res_ref", "vis_REF"]
    df["p_Pa"] = df["p_MPa"] * 1e6
    return df


def load_tc_samples(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", header=0, skipinitialspace=True)
    # Material T/K p/kPa den s TC_EXP TC_RES TC_REF
    df.columns = ["material", "T", "p_kPa", "den", "s_res_ref", "tc_exp", "tc_res_ref", "tc_REF"]
    df["p_Pa"] = df["p_kPa"] * 1e3
    return df


def state(backend: str, fluid: str, p_Pa: float, T: float):
    """Return an updated AbstractState, or None if unavailable/unsolvable."""
    try:
        AS = CP.AbstractState(backend, fluid)
        AS.update(CP.PT_INPUTS, p_Pa, T)
        return AS
    except Exception:
        return None


def compare(samples: pd.DataFrame, prop: str) -> pd.DataFrame:
    """prop is 'viscosity' or 'conductivity'."""
    is_vis = prop == "viscosity"
    ref_col = "vis_res_ref" if is_vis else "tc_res_ref"
    exp_col = "vis_exp" if is_vis else "tc_exp"
    # Martinek tabulates viscosity in uPa*s; CoolProp returns Pa*s.
    scale = 1e6 if is_vis else 1.0

    rows = []
    for _, r in samples.iterrows():
        name, T, p_Pa = r["material"], r["T"], r["p_Pa"]
        rec = {"material": name, "T": T, "p_Pa": p_Pa}

        heos = state("HEOS", name, p_Pa, T)
        if heos is None:
            rows.append({**rec, "status": "HEOS unavailable"})
            continue

        # Drive the actual C++ RES implementation.
        try:
            if is_vis:
                heos.use_viscosity_RES(True)
                value = heos.viscosity() * scale
            else:
                heos.use_conductivity_RES(True)
                value = heos.conductivity() * scale
        except Exception as e:
            # Fluids whose HEOS entry was deliberately removed land here, by design.
            rows.append({**rec, "status": f"HEOS RES error: {str(e)[:90]}"})
            continue

        rec["s_res_HEOS"] = heos.smolar_residual()
        rec["rho_HEOS"] = heos.rhomass()
        rec["phase_HEOS"] = str(heos.phase())

        # REFPROP is diagnostic only -- it cannot run RES (the API lives on the
        # Helmholtz/cubic backends), and it does not need to: the reference column
        # already is the REFPROP-based RES value published with the paper.
        ref = state("REFPROP", name, p_Pa, T)
        if ref is not None:
            rec["s_res_REFPROP"] = ref.smolar_residual()
            rec["rho_REFPROP"] = ref.rhomass()
            rec["phase_REFPROP"] = str(ref.phase())
            rec["s_res_dev%"] = pct(rec["s_res_HEOS"], rec["s_res_REFPROP"])
            rec["rho_dev%"] = pct(rec["rho_HEOS"], rec["rho_REFPROP"])
            rec["phase_mismatch"] = rec["phase_HEOS"] != rec["phase_REFPROP"]

        rec[f"{prop}_HEOS"] = value
        rec[f"{prop}_reference"] = r[ref_col]
        rec["transport_dev%"] = pct(value, r[ref_col])
        rec["experiment"] = r[exp_col]
        rec["dev_vs_experiment%"] = pct(value, r[exp_col])
        rec["status"] = "ok"
        rows.append(rec)

    return pd.DataFrame(rows)


def report(df: pd.DataFrame, prop: str, threshold: float, csv_path: str):
    label = prop.upper()
    print(f"\n{'-' * 78}")
    ok = df[df["status"] == "ok"].copy()
    print(f"{label}: {len(ok)} points evaluated through the C++ RES implementation")

    if ok.empty:
        print("  nothing evaluated")
    else:
        flagged = ok[ok["transport_dev%"].abs() > threshold]
        if flagged.empty:
            print(f"  no fluid deviates from the reference by more than {threshold}%")
        else:
            print(f"\n  deviating from REFPROP-fitted reference by more than {threshold}%:")
            cols = ["material", "T", "s_res_dev%", "rho_dev%", "phase_mismatch", "transport_dev%"]
            cols = [c for c in cols if c in flagged.columns]
            print(flagged[cols].to_string(index=False))

            if "phase_mismatch" in flagged.columns:
                # fillna guards rows where REFPROP was unavailable and the flag is NaN.
                straddle = flagged[flagged["phase_mismatch"].fillna(False).astype(bool)]
                if not straddle.empty:
                    print(
                        "\n  NOTE: the following are phase-selection artefacts -- HEOS and REFPROP\n"
                        "  disagree on the phase at this specific state point, so the deviation says\n"
                        "  nothing about the fluid's parameters. Verify at other state points before\n"
                        "  excluding them:"
                    )
                    print("    " + ", ".join(sorted(set(straddle["material"]))))

    skipped = df[df["status"] != "ok"]
    if not skipped.empty:
        print(f"\n  {len(skipped)} points not evaluated:")
        for status, grp in skipped.groupby("status"):
            names = sorted(set(grp["material"]))
            shown = ", ".join(names[:8]) + (f", ... (+{len(names) - 8})" if len(names) > 8 else "")
            print(f"    {len(grp):3d}  {status}")
            print(f"         {shown}")

    print(f"\n  -> {csv_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD_PCT,
                    help=f"flag deviations above this %% (default {DEFAULT_THRESHOLD_PCT})")
    args = ap.parse_args()

    if not hasattr(CP.AbstractState("HEOS", "Water"), "use_viscosity_RES"):
        raise SystemExit(
            "The installed CoolProp has no RES API. Rebuild from this working tree:\n"
            "    python -m pip install ."
        )

    vis_df = compare(load_vis_samples(VIS_PURE_OUT), "viscosity")
    vis_csv = os.path.join(OUT_DIR, "compare_vis_pure.csv")
    vis_df.to_csv(vis_csv, index=False, float_format="%.6g")
    report(vis_df, "viscosity", args.threshold, vis_csv)

    tc_df = compare(load_tc_samples(TC_PURE_REF), "conductivity")
    tc_csv = os.path.join(OUT_DIR, "compare_tc_pure.csv")
    tc_df.to_csv(tc_csv, index=False, float_format="%.6g")
    report(tc_df, "conductivity", args.threshold, tc_csv)

    print(f"\n{'-' * 78}")
    print("Fluids flagged above (excluding phase-selection artefacts) are candidates for")
    print("removal from the HEOS entries in dev/res_transport_parameters.json, via the")
    print("exclusion sets in dev/convert_RES_csv_to_json.py.")


if __name__ == "__main__":
    main()
