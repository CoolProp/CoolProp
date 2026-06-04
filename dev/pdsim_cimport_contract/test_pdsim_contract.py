"""
Pytest wrapper for the PDSim cimport contract.

Why a subprocess?  This test must exercise the *installed* CoolProp, but it
lives in a repo that also contains the uncompiled *source* CoolProp package.
If the pytest process imported CoolProp directly, sys.path could shadow the
installed build with the unbuilt source.  So the pytest process imports
nothing from CoolProp -- it just runs run_contract.py (build + value checks)
in a subprocess whose cwd makes the installed CoolProp authoritative, and
asserts it exits 0.

The real guarantees live in run_contract.py / pdsim_surface.pyx:
  * compile failure  => the cimport surface drifted (PDSim would break)
  * value mismatch   => the kPa/kJ unit contract drifted (PDSim silently wrong)
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def test_pdsim_cimport_contract():
    proc = subprocess.run(
        [sys.executable, "run_contract.py"],
        cwd=HERE, capture_output=True, text=True,
    )
    assert proc.returncode == 0, (
        "PDSim cimport contract failed -- the v8 Cython surface or its unit "
        "convention drifted from what PDSim depends on:\n\n"
        + proc.stdout + "\n" + proc.stderr
    )
