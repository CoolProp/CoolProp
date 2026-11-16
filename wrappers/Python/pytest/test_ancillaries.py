"""Saturation ancillaries against the equation of state (#2266).

The pS, rhoL and rhoV ancillaries give the initial guesses for saturation
solves and, with superancillaries disabled, decide the phase of P,T inputs.
#2266 came from rhoV fits of n-Hexane and Ammonia that carried ``using_tau_r``
by mistake, which made TTSE/BICUBIC saturation return inf.
"""
import json
import os

import numpy as np
import pytest

from CoolProp import CoolProp as CP

# Transcription errors are off by tens of percent.  The worst genuine fit in
# the library is about 5 % between 0.6 Tc and 0.95 Tc; closer to the triple
# point some heavy fluids' fits are much worse, so that range is not checked.
TOLERANCE_PCT = 10.0

FLUIDS = sorted(CP.get_global_param_string("FluidsList").split(","))


@pytest.fixture
def maxwell_saturation():
    """Solve saturation from the EOS, not from the superancillary."""
    if not hasattr(CP, "ENABLE_SUPERANCILLARIES"):
        yield
        return
    prev = CP.get_config_bool(CP.ENABLE_SUPERANCILLARIES)
    CP.set_config_bool(CP.ENABLE_SUPERANCILLARIES, False)
    try:
        yield
    finally:
        CP.set_config_bool(CP.ENABLE_SUPERANCILLARIES, prev)


@pytest.mark.parametrize("fluid", FLUIDS)
def test_ancillaries_match_eos(fluid, maxwell_saturation):
    jfluid = json.loads(CP.get_fluid_param_string(fluid, "JSON"))[0]
    if jfluid["EOS"][0].get("pseudo_pure", False):
        pytest.skip("pseudo-pure fluids have no Maxwell saturation state")
    AS = CP.AbstractState("HEOS", fluid)
    Tc = AS.T_critical()
    bad = []
    for key, param, Q in (("pS", CP.iP, 1), ("rhoL", CP.iDmolar, 0), ("rhoV", CP.iDmolar, 1)):
        anc = jfluid["ANCILLARIES"][key]
        Tlo = max(0.6 * Tc, anc["Tmin"], AS.Ttriple())
        Thi = min(0.95 * Tc, anc["Tmax"])
        for T in np.linspace(Tlo, Thi, 8):
            AS.update(CP.QT_INPUTS, float(Q), T)
            err = 100 * abs(AS.saturation_ancillary(param, Q, CP.iT, T) / AS.keyed_output(param) - 1)
            if not err <= TOLERANCE_PCT:
                bad.append(f"{key} at T = {T:.2f} K: {err:.3g} %")
    assert not bad, f"{fluid} ancillaries disagree with the EOS: " + "; ".join(bad)


@pytest.fixture(scope="module")
def tables_dir(tmp_path_factory):
    prev = CP.get_config_string(CP.ALTERNATIVE_TABLES_DIRECTORY)
    # path_to_tables() appends the backend name without a separator
    CP.set_config_string(CP.ALTERNATIVE_TABLES_DIRECTORY, str(tmp_path_factory.mktemp("tables")) + os.sep)
    yield
    CP.set_config_string(CP.ALTERNATIVE_TABLES_DIRECTORY, prev)


@pytest.mark.parametrize("backend", ["TTSE", "BICUBIC"])
@pytest.mark.parametrize("fluid", ["R717", "Hexane"])
def test_tabular_saturation_2266(fluid, backend, tables_dir, maxwell_saturation):
    heos = CP.AbstractState("HEOS", fluid)
    tab = CP.AbstractState(backend + "&HEOS", fluid)
    for Q, T in ((0.3, 283.15), (0.5, 293.15), (0.7, 303.15)):
        heos.update(CP.QT_INPUTS, Q, T)
        tab.update(CP.QT_INPUTS, Q, T)
        assert tab.p() == pytest.approx(heos.p(), rel=1e-6)
