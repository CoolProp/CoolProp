"""Evaluate the RES model over dev/RES_comparison/grid_points.csv using the AUTHORS' OWN code.

This is the reference side of the RES validation harness.  It drives the two vendored reference
implementations in dev/RES_reference/ -- Martinek 2025 for viscosity, Li 2024 for conductivity --
and writes one CSV per equation of state:

    --eos REFPROP  ->  dev/RES_comparison/grid_ref_refprop.csv
    --eos HEOS     ->  dev/RES_comparison/grid_ref_heos.csv

Those two answer different questions:

  * REFPROP vs HEOS, both from THIS script, measures whether the REFPROP-fitted coefficients
    transfer to CoolProp's own equation of state.  One implementation, two equations of state, so
    what is left is the equation of state and nothing else.  dev/RES_grid_report.py turns that
    into the HEOS_TRANSFER_EXCLUDE list in dev/convert_RES_csv_to_json.py.  No CoolProp RES
    implementation is involved, which is what lets the exclusion list be derived before the C++
    model exists.

  * HEOS from this script vs CoolProp's own RES implementation measures whether CoolProp
    implements the model correctly -- same model, same equation of state, so any deviation is a
    defect.  That comparison arrives with the model; nothing in this repository consumes the
    HEOS column yet.

Both runs pin the model to CoolProp's own choices -- the dilute-gas term to the fitted polynomial,
and the conductivity enhancement viscosity to the RES viscosity -- because the two papers made
opposite choices there per code path, and leaving them alone would measure that choice rather
than the equation of state or the implementation.  Those are options on the vendored modules,
whose defaults remain the published behaviour; see dev/RES_reference/README.md.

    python dev/RES_reference_run.py --eos REFPROP
    python dev/RES_reference_run.py --self-test      # vendored code still reproduces the papers
"""

import argparse
import csv
import filecmp
import math
import os
import sys

import CoolProp.CoolProp as CP

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
OUT_DIR = os.path.join(SCRIPT_DIR, "RES_comparison")
GRID_CSV = os.path.join(OUT_DIR, "grid_points.csv")
REFERENCE_DIR = os.path.join(SCRIPT_DIR, "RES_reference")

sys.path.insert(0, SCRIPT_DIR)
from RES_reference.li2024_conductivity import conductivity_RES as tc_ref  # noqa: E402
from RES_reference.martinek2025_viscosity import viscosity_RES as vis_ref  # noqa: E402

# What CoolProp does, and therefore what both reference runs are pinned to.  See the module
# docstring and dev/RES_reference/README.md.
PINNED_VIS = dict(dilute_source="polynomial")
# zero_ind_fit: matches the shipped table, which follows Li's paper rather than the SI's flag for
# NEOPENTN (conductivity only).  See dev/RES_reference/README.md.
PINNED_TC = dict(dilute_source="polynomial", enhancement_viscosity="res", zero_ind_fit="global")


def scalar(v):
    """Both modules return a 1-element numpy array for pure fluids and a float for mixtures."""
    try:
        return float(v[0]) if hasattr(v, "__len__") else float(v)
    except (TypeError, IndexError):
        return float(v)


def as_reference_args(fluid, mole_fractions):
    """Both modules take a bare name for a pure fluid and a LIST for a mixture.

    That is how they pick their pure vs mixture branch, so a binary has to be passed as a list
    even where it could be spelled as one string.
    """
    if "&" not in fluid:
        return fluid, [1.0]
    names = fluid.split("&")
    return names, [float(x) for x in mole_fractions.split(";")]


def state_columns(eos, names, z, T, p):
    """Bulk EOS quantities at the grid point, for the report's weighting and gate checks."""
    AS = CP.AbstractState(eos, "&".join(names))
    if len(names) > 1:
        AS.set_mole_fractions(z)
    AS.update(CP.PT_INPUTS, p, T)
    out = {"rho": AS.rhomass(), "s_res": AS.smolar_residual()}
    # Informational only -- the report uses it to COUNT points where the two sides disagree about
    # the phase, never to decide anything.  REFPROP implements calc_phase for pure fluids but not
    # for mixtures, so letting this raise would discard every mixture point over a column nothing
    # depends on.
    try:
        out["phase"] = AS.phase()
    except Exception:
        out["phase"] = ""
    # A mixture critical point costs a solve that may not converge; the report only needs the
    # enhancement gate for pure fluids, where the exclusions are decided.
    if len(names) == 1:
        out["Tc"] = AS.T_critical()
        out["rhoc"] = AS.rhomass_critical()
    else:
        out["Tc"] = ""
        out["rhoc"] = ""
    return out


def tref_deriv_ok(eos, names, T_ref, rho):
    """Whether the reference code can reach (d p / d rho)_T at its reference temperature.

    Li's enhancement needs that derivative at T_ref (typically 1.5*Tc, i.e. at or just past some
    fluids' upper EOS limit).  It reaches it through PropsSI, and where PropsSI raises, TC_RES
    catches and silently returns zero enhancement.  CoolProp evaluates alpha^r there directly and
    keeps the enhancement -- a deliberate divergence -- so the report has to be able to tell the
    two populations apart rather than reading the difference as an implementation defect.
    """
    if not (T_ref and T_ref > 0) or len(names) > 1:
        return ""
    try:
        CP.PropsSI("d(P)/d(Dmass)|T", "T", T_ref, "Dmass", rho, "{}::{}".format(eos, names[0]))
        return 1
    except Exception:
        return 0


def self_test():
    """Check the vendored copies still reproduce the papers' own published sample tables.

    Every adaptation in dev/RES_reference/ is a keyword argument defaulting to the published
    behaviour, so calling them with no options must regenerate the supporting information's own
    output files byte-for-byte.  That is a direct check that adapting the code did not change the
    model -- something the earlier monkey-patching version of this script could not do; it could
    only assert that its text substitutions had matched something.

    Runs on REFPROP, because that is what the papers ran on.
    """
    import shutil
    import tempfile

    cases = [
        ("Martinek 2025 viscosity", vis_ref.write_sample_table,
         os.path.join(REFERENCE_DIR, "martinek2025_viscosity", "Samples_pure_fluids_output.txt")),
        ("Li 2024 conductivity", tc_ref.write_sample_table,
         os.path.join(REFERENCE_DIR, "li2024_conductivity", "Table_S5_SI.txt")),
        # The binary tables matter separately: they are the only thing that exercises the Wilke
        # rule, the mole-fraction averaging of the residual coefficients, and (for conductivity)
        # the mixture branch of the critical enhancement.
        ("Martinek 2025 viscosity, binaries", vis_ref.write_binary_sample_table,
         os.path.join(REFERENCE_DIR, "martinek2025_viscosity", "Samples_binaries_output.txt")),
        ("Li 2024 conductivity, binaries", tc_ref.write_binary_sample_table,
         os.path.join(REFERENCE_DIR, "li2024_conductivity", "Table_S6_SI.txt")),
    ]
    failures = 0
    # Not TemporaryDirectory(): the vendored writers close their output file on the normal path
    # only, so a writer that raises leaves a handle open, and on Windows the automatic cleanup
    # then fails with a PermissionError that buries the actual error underneath it.
    tmp = tempfile.mkdtemp(prefix="res_reference_selftest_")
    try:
        for label, writer, published in cases:
            produced = os.path.join(tmp, os.path.basename(published))
            writer(produced)
            if filecmp.cmp(produced, published, shallow=False):
                print("  OK   {}: regenerates {} byte-for-byte".format(label, os.path.basename(published)))
                continue
            failures += 1
            print("  FAIL {}: output differs from the published {}".format(label, os.path.basename(published)))
            with open(produced, encoding="utf-8") as a, open(published, encoding="utf-8") as b:
                for i, (la, lb) in enumerate(zip(a, b), 1):
                    if la != lb:
                        print("       first difference at line {}:".format(i))
                        print("         ours:      {}".format(la.rstrip()))
                        print("         published: {}".format(lb.rstrip()))
                        break
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return failures


def main(eos, out_path, limit):
    import json

    with open(os.path.join(SCRIPT_DIR, "res_transport_parameters.json"), encoding="utf-8") as fh:
        res_json = json.load(fh)
    t_ref_of = {
        f: e["critical_enhancement"]["t_ref"] for f, e in res_json["conductivity"].items() if "critical_enhancement" in e
    }

    with open(GRID_CSV, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if limit:
        rows = rows[:limit]

    fields = ["fluid", "T_K", "p_Pa", "region", "mole_fractions", "eos", "ok", "err",
              "eta", "tc", "eta0", "tc0", "rho", "s_res", "Tc", "rhoc", "phase", "tref_deriv_ok"]
    n_ok = n_fail = 0
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            fluid, T, p = r["fluid"], float(r["T_K"]), float(r["p_Pa"])
            ref_fluid, z = as_reference_args(fluid, r["mole_fractions"])
            names = fluid.split("&")
            out = dict(r, eos=eos, ok=0, err="", eta="", tc="", eta0="", tc0="",
                       rho="", s_res="", Tc="", rhoc="", phase="", tref_deriv_ok="")
            try:
                eta, _rho, _s, _g, eta0 = vis_ref.viscosity_RES(ref_fluid, z, p, T, eos=eos, **PINNED_VIS)
                tc, _rho, _s, _g, tc0 = tc_ref.TC_RES(ref_fluid, z, p, T, eos=eos, **PINNED_TC)
                eta, tc, eta0, tc0 = (scalar(v) for v in (eta, tc, eta0, tc0))
                out.update(state_columns(eos, names, z, T, p))
                out.update(eta=repr(eta), tc=repr(tc), eta0=repr(eta0), tc0=repr(tc0), ok=1,
                           tref_deriv_ok=tref_deriv_ok(eos, names, t_ref_of.get(fluid), out["rho"]))
                if not (math.isfinite(eta) and math.isfinite(tc)):
                    out["ok"] = 0
                    out["err"] = "non-finite"
            # SystemExit is expected -- both vendored modules call sys.exit() on an unknown
            # fluid -- but NOT KeyboardInterrupt, which BaseException would swallow into a
            # CSV cell and make a 3570-point run impossible to stop.
            except (Exception, SystemExit) as exc:
                out["err"] = "{}: {}".format(type(exc).__name__, str(exc).splitlines()[0][:80])
            n_ok += out["ok"] == 1
            n_fail += out["ok"] != 1
            w.writerow(out)

    print("Wrote {}".format(out_path))
    print("  {} points evaluated on {}, {} failed".format(n_ok, eos, n_fail))
    if n_ok == 0:
        raise SystemExit("no grid point evaluated successfully; the report would read this as "
                         "'no evidence' rather than as a failure")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--eos", choices=["REFPROP", "HEOS"], help="equation of state to run the reference code on")
    ap.add_argument("--self-test", action="store_true",
                    help="check the vendored reference code still reproduces the papers' published sample tables")
    ap.add_argument("--out", default=None, help="output CSV (default: dev/RES_comparison/grid_ref_<eos>.csv)")
    ap.add_argument("--limit", type=int, default=0, help="evaluate only the first N grid points (smoke test)")
    a = ap.parse_args()
    if a.self_test:
        print("Checking the vendored reference code against the published sample tables:")
        raise SystemExit(1 if self_test() else 0)
    if not a.eos:
        ap.error("--eos is required unless --self-test is given")
    out = a.out or os.path.join(OUT_DIR, "grid_ref_{}.csv".format(a.eos.lower()))
    main(a.eos, out, a.limit)
