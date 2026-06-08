"""
Increment 3: build the UNCHANGED #3079 contract module (pdsim_surface.pyx, with
`from CoolProp.State cimport State`) against a `CoolProp` *package* whose `State`
is the capsule shim -- proving PDSim's exact cimport surface works with no
Cython AbstractState.

Package layout assembled here (the v8 target):
    CoolProp/__init__.py            ->  from .CoolProp import *   (nanobind core)
    CoolProp/CoolProp.*.so          ->  the nanobind core (exports _capi, enums, PropsSI)
    CoolProp/State.{pxd,*.so}       ->  the frozen shim (cimport target)
    CoolProp/constants_header.pxd   ->  generated, real enum values
Then `pdsim_surface.pyx` (copied verbatim from #3079) cimports CoolProp.State
(the shim) + CoolProp.constants_header (real), and run_contract's checks run.
"""
import glob
import math
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WT = os.path.abspath(os.path.join(HERE, "..", ".."))            # nanobind worktree root
MAIN = "/Users/ianbell/Code/CoolProp"                            # for build*/_deps
INCLUDE = os.path.join(WT, "include")
NB_BUILD = os.path.join(WT, "wrappers", "Python", "nanobind", "build")
PKG_ROOT = os.path.join(HERE, "_pkg")
PKG = os.path.join(PKG_ROOT, "CoolProp")
PDSIM_SURFACE = ("/Users/ianbell/Code/CoolProp/.claude/worktrees/pdsim-contract/"
                 "dev/pdsim_cimport_contract/pdsim_surface.pyx")
FLUID, T, RHO = "Water", 320.0, 1000.0


def find_dep(header_rel, fetch_glob):
    cands = [INCLUDE] + glob.glob(os.path.join(MAIN, fetch_glob)) + ["/opt/homebrew/include"]
    for c in cands:
        if os.path.isfile(os.path.join(c, header_rel)):
            return c
    return None


def build_package():
    shutil.rmtree(PKG_ROOT, ignore_errors=True)
    os.makedirs(PKG)
    core = glob.glob(os.path.join(NB_BUILD, "CoolProp.cpython-*.so"))[0]
    shutil.copy(core, os.path.join(PKG, os.path.basename(core)))         # -> CoolProp.CoolProp
    with open(os.path.join(PKG, "__init__.py"), "w") as fp:
        fp.write("from .CoolProp import *\n")
    # generate constants_header.pxd from this worktree, copy into the package
    subprocess.run([sys.executable, "generate_constants_module.py"],
                   cwd=os.path.join(WT, "wrappers", "Python"), check=True, capture_output=True)
    shutil.copy(os.path.join(WT, "wrappers", "Python", "CoolProp", "constants_header.pxd"),
                os.path.join(PKG, "constants_header.pxd"))
    for f in ("State.pyx", "State.pxd"):
        shutil.copy(os.path.join(HERE, f), os.path.join(PKG, f))
    shutil.copy(PDSIM_SURFACE, os.path.join(PKG_ROOT, "pdsim_surface.pyx"))  # VERBATIM from #3079


def build_extensions():
    from setuptools import Extension, setup
    from Cython.Build import cythonize

    incs = [INCLUDE, os.path.join(sys.prefix, "include", f"python{sys.version_info.major}.{sys.version_info.minor}")]
    for hdr, g in [("fmt/format.h", "build*/_deps/fmt-src/include")]:
        d = find_dep(hdr, g)
        if d:
            incs.append(d)
    exts = [
        Extension("CoolProp.State", ["CoolProp/State.pyx"], language="c++",
                  include_dirs=incs, extra_compile_args=["-std=c++17"]),
        Extension("pdsim_surface", ["pdsim_surface.pyx"], language="c++",
                  include_dirs=incs, extra_compile_args=["-std=c++17"]),
    ]
    # run setup with argv pointing at PKG_ROOT
    old = sys.argv
    sys.argv = ["setup.py", "build_ext", "--inplace"]
    cwd = os.getcwd()
    os.chdir(PKG_ROOT)
    try:
        setup(name="state_capsule_inc3", ext_modules=cythonize(exts, language_level=3, include_path=[PKG_ROOT]),
              script_args=["build_ext", "--inplace"])
    finally:
        os.chdir(cwd)
        sys.argv = old


def run_checks():
    sys.path.insert(0, PKG_ROOT)
    import CoolProp as CP          # the shim package (core re-exported)
    import pdsim_surface           # built against CoolProp.State (shim) + constants

    def SI(name):
        return CP.PropsSI(name, "T", T, "Dmass", RHO, FLUID)

    fails = []

    def expect(label, got, want, rel=1e-6):
        if not math.isclose(got, want, rel_tol=rel):
            fails.append(f"{label}: got {got:.8g}, want {want:.8g}")

    out = pdsim_surface.exercise(FLUID, T, RHO)
    g = out["getters"]
    expect("get_p [kPa]", g["p"], SI("P") / 1000.0)
    expect("get_h [kJ/kg]", g["h"], SI("Hmass") / 1000.0)
    expect("get_cp0 [kJ/kg/K]", g["cp0"], SI("Cp0mass") / 1000.0)
    expect("get_speed_sound [SI]", g["speed_sound"], SI("speed_of_sound"))
    expect("cached p_ [kPa]", out["p_"], g["p"], rel=1e-9)
    expect("Props() [SI]", out["h_via_Props"], SI("Hmass"))
    expect("pAS.keyed_output [SI]", out["ko_h"], SI("Hmass"))
    expect("pAS.cpmass [SI]", out["cp_pAS"], SI("Cpmass"))
    expect("get_dpdT==fpd/1000", out["dpdT_pAS"] / 1000.0, g["dpdT"])
    expect("copy().get_h [kJ/kg]", out["h_copy"], SI("Hmass") / 1000.0)
    if not (out["names"] and FLUID.lower() in out["names"][0].lower()):
        fails.append(f"fluid_names: {out['names']}")
    tp = pdsim_surface.exercise_twophase(FLUID, 373.0, 0.4)
    expect("two-phase Q", tp["Q"], 0.4, rel=1e-9)
    acc = pdsim_surface.hot_loop(FLUID, T, RHO, 1000)
    if not (math.isfinite(acc) and acc != 0.0):
        fails.append(f"hot_loop {acc}")

    if fails:
        sys.stderr.write("\n*** INCREMENT 3 MISMATCHES ***\n  " + "\n  ".join(fails) + "\n")
        return 1
    print("STATE CAPSULE (increment 3): PASS")
    print("  The UNCHANGED #3079 pdsim_surface.pyx (`from CoolProp.State cimport State`)")
    print("  compiled + ran against the capsule shim package -- PDSim's cimport surface")
    print("  is satisfied with no Cython AbstractState.")
    return 0


if __name__ == "__main__":
    build_package()
    build_extensions()
    raise SystemExit(run_checks())
