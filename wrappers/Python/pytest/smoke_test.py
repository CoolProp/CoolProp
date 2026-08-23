"""
CI Smoke Test Suite for CoolProp Wheel Verification.

Validates that compiled C++ extensions (nanobind) load successfully, executes 
basic thermodynamic property lookups (PropsSI), and verifies core exception 
handling. 

Designed to run under both native desktop wheels and WebAssembly (Pyodide) 
runtimes with minimal external dependencies.

Usage:
    pytest wrappers/Python/pytest/smoke_test.py
"""
import pytest

# Guard against broken binary wheel installations in CI
import numpy 

import CoolProp
from CoolProp.CoolProp import PropsSI, PhaseSI
from CoolProp.HumidAirProp import HAPropsSI

# Verify the Cython sidecar capsule shim loads alongside the nanobind core
from CoolProp.State import State  

def test_build_version_injection():
    """Ensure CMake properly resolved and injected the git revision during build."""
    assert '?' not in CoolProp.__gitrevision__, "Git revision was not properly injected during build"


def test_simple_propssi():
    assert round(PropsSI("T", "P", 101325, "Q", 0, "Water"), 3) == 373.124

    with pytest.raises(ValueError):
        PropsSI("T", "P", 101325, "Q", 0, "Walter")


def test_simple_phasesi():
    assert PhaseSI("P", 101325, "Q", 0, "Water") == "twophase"


def test_simple_HAPropsSI():
    assert round(HAPropsSI("H", "T", 298.15, "P", 101325, "R", 0.5), 3) == 50423.450

    with pytest.raises(ValueError):
        HAPropsSI("H", "T", 298.15, "P", -101325, "R", 0.5)