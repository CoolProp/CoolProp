"""HAProps consumer migration (bd CoolProp-r9sq.27).

The non-SI ``HAProps`` was removed in v8; the live in-repo consumers
(``Web/.../NelsonValidation.py`` and ``CoolProp/GUI/PsychScript.py``) now call
``HAPropsSI`` directly in SI units.  (``CoolProp/Plots/psy.py`` is legacy
Python-2 / PyQt4 code that does not import on Python 3 and is left untouched.)

This runs the standalone Nelson validation script end-to-end (it needs only
CoolProp, no matplotlib) and checks it produces the expected moist-air values --
proving the SI->display conversion direction is right (notably that the enthalpy
column comes back in kJ/kg, i.e. the SI J/kg value divided by 1000).
"""
import os
import sys
import subprocess

import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_NELSON = os.path.normpath(
    os.path.join(_HERE, "..", "..", "..", "Web", "fluid_properties", "Validation", "NelsonValidation.py")
)


@pytest.mark.skipif(not os.path.exists(_NELSON), reason="NelsonValidation.py not found next to the repo")
def test_nelson_validation_runs_and_is_sane():
    proc = subprocess.run([sys.executable, _NELSON], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    out = proc.stdout
    assert "kJ/kg_da" in out  # the enthalpy column header

    # The Tdb=25, Twb=20 row carries known psychrometric values at 1 atm:
    # W ~ 0.0127 kg/kg, h ~ 57.4 kJ/kg, v ~ 0.861 m^3/kg.  Crucially h is ~57 (kJ/kg),
    # NOT ~57000 (J/kg) -- i.e. the SI enthalpy was divided by 1000 for the column.
    # (The Twb==Tdb "saturation" rows are legitimately omitted: the script's
    # tdb+273.13 vs twb+273.15 typo puts Twb 0.02 K above Tdb -> RH slightly > 1.)
    for line in out.splitlines():
        cols = line.split()
        if len(cols) == 7 and cols[0] == "25.00" and cols[1] == "20.00":
            W, h, v = float(cols[4]), float(cols[5]), float(cols[6])
            assert 0.011 < W < 0.014, "W (humidity ratio) out of range: %r" % W
            assert 55.0 < h < 60.0, "h not in kJ/kg (kSI conversion wrong?): %r" % h
            assert 0.85 < v < 0.87, "v (specific volume) out of range: %r" % v
            break
    else:
        pytest.fail("expected Tdb=25,Twb=20 data row not found in NelsonValidation output")
