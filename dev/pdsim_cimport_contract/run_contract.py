#!/usr/bin/env python3
"""
PDSim cimport-contract runner (standalone + the worker the pytest wrapper calls).

Builds pdsim_surface.pyx against the *installed* CoolProp and verifies the
values it produces through the legacy Cython surface.  Run directly:

    python dev/pdsim_cimport_contract/run_contract.py

Exits 0 on success, 1 on any failure.  See SURFACE.md for the contract and
test_pdsim_contract.py for why this is a subprocess (sys.path hygiene).
"""
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FLUID, T, RHO = "Water", 320.0, 1000.0   # compressed liquid: every property defined


def build():
    proc = subprocess.run(
        [sys.executable, "setup_contract.py", "build_ext", "--inplace"],
        cwd=HERE, capture_output=True, text=True,
    )
    if proc.returncode != 0:
        sys.stderr.write(
            "\n*** PDSim cimport contract FAILED TO COMPILE ***\n"
            "The cimportable State / AbstractState / constants_header surface "
            "drifted from what PDSim depends on:\n\n"
            + proc.stdout[-4000:] + "\n" + proc.stderr[-4000:] + "\n")
        raise SystemExit(1)


def check():
    sys.path.insert(0, HERE)
    import pdsim_surface
    from CoolProp.CoolProp import PropsSI

    SI = lambda name: PropsSI(name, "T", T, "Dmass", RHO, FLUID)
    out = pdsim_surface.exercise(FLUID, T, RHO)
    g = out["getters"]

    fails = []

    def expect(label, got, want, rel=1e-6):
        if not math.isclose(got, want, rel_tol=rel):
            fails.append(f"{label}: got {got:.8g}, want {want:.8g}")

    # --- THE contract: legacy State is kPa/kJ, pAS/Props are SI ---
    expect("get_p [kPa]", g["p"], SI("P") / 1000.0)
    expect("get_h [kJ/kg]", g["h"], SI("Hmass") / 1000.0)
    expect("get_s [kJ/kg/K]", g["s"], SI("Smass") / 1000.0)
    expect("get_u [kJ/kg]", g["u"], SI("Umass") / 1000.0)
    expect("get_cp [kJ/kg/K]", g["cp"], SI("Cpmass") / 1000.0)
    expect("get_cv [kJ/kg/K]", g["cv"], SI("Cvmass") / 1000.0)
    expect("get_speed_sound [m/s, SI]", g["speed_sound"], SI("speed_of_sound"))
    expect("get_rho [kg/m^3, SI]", g["rho"], RHO, rel=1e-9)

    expect("cached T_", out["T_"], T, rel=1e-9)
    expect("cached rho_", out["rho_"], RHO, rel=1e-9)
    expect("cached p_ (kPa)", out["p_"], g["p"], rel=1e-9)

    expect("Props() [SI]", out["h_via_Props"], SI("Hmass"))
    expect("pAS.keyed_output [SI]", out["ko_h"], SI("Hmass"))
    expect("pAS.rhomass [SI]", out["rho_pAS"], RHO, rel=1e-9)
    expect("pAS.T [SI]", out["T_pAS"], T, rel=1e-9)
    expect("pAS.cpmass [SI]", out["cp_pAS"], SI("Cpmass"))
    expect("pAS.cvmass [SI]", out["cv_pAS"], SI("Cvmass"))

    expect("get_dpdT == fpd(iP,iT,iDmolar)/1000", out["dpdT_pAS"] / 1000.0, g["dpdT"])
    expect("copy().get_h [kJ/kg]", out["h_copy"], SI("Hmass") / 1000.0)

    # two-phase get_Q + construction
    tp = pdsim_surface.exercise_twophase(FLUID, 373.0, 0.4)
    expect("two-phase get_Q", tp["Q"], 0.4, rel=1e-9)
    expect("two-phase p [kPa]", tp["p"],
           PropsSI("P", "T", 373.0, "Q", 0.4, FLUID) / 1000.0)

    # hot loop (PDSim's inner pattern) just has to run and produce a finite sum
    acc = pdsim_surface.hot_loop(FLUID, T, RHO, 1000)
    if not (math.isfinite(acc) and acc != 0.0):
        fails.append(f"hot_loop returned {acc}")

    if not (out["names"] and FLUID.lower() in out["names"][0].lower()):
        fails.append(f"fluid_names mismatch: {out['names']}")

    if fails:
        sys.stderr.write("\n*** PDSim contract VALUE MISMATCHES ***\n  "
                         + "\n  ".join(fails) + "\n")
        raise SystemExit(1)
    print(f"PDSim cimport contract OK: {FLUID} @ T={T} K, rho={RHO} kg/m^3 "
          f"-- full State/AbstractState/constants surface compiles, links "
          f"(no -lCoolProp), and matches PropsSI under the legacy kPa/kJ "
          f"unit convention.")


if __name__ == "__main__":
    build()
    check()
