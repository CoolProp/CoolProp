"""
AbstractState.update_nogil: update() with the GIL released (COO-141).

update() keeps the GIL, so code that shares one state between Python threads
stays safe. update_nogil is the opt-in for code that gives every thread its
own state: other Python threads run while CoolProp computes. These tests
check that it computes the same states as update(), that distinct states,
also of the same fluid, can be updated from several threads at once, and
that a failed update still raises ValueError.
"""
from concurrent.futures import ThreadPoolExecutor

import pytest

import CoolProp
import CoolProp.CoolProp as CP


def _state(fluid):
    """A HEOS state of a pure fluid, or of a mixture given as (names, mole fractions)."""
    if isinstance(fluid, tuple):
        names, z = fluid
        state = CP.AbstractState("HEOS", "&".join(names))
        state.set_mole_fractions(z)
        return state
    return CP.AbstractState("HEOS", fluid)


# Pressure and temperature points on either side of the dome, for a pure fluid and a blend
CASES = [("R134a", [(1e5, T) for T in (250.0, 300.0, 350.0)] + [(3e6, T) for T in (300.0, 400.0)]),
         ((["R32", "R125"], [0.7, 0.3]), [(1e5, T) for T in (250.0, 300.0)] + [(2e6, 350.0)])]


def _densities(fluid, points, update_name):
    state = _state(fluid)
    update = getattr(state, update_name)
    out = []
    for p_Pa, T_K in points:
        update(CoolProp.PT_INPUTS, p_Pa, T_K)
        out.append(state.rhomolar())
    return out


@pytest.mark.parametrize("fluid,points", CASES)
def test_update_nogil_gives_the_states_of_update(fluid, points):
    assert _densities(fluid, points, "update_nogil") == _densities(fluid, points, "update")


@pytest.mark.parametrize("fluid,points", CASES)
def test_distinct_states_update_in_parallel(fluid, points):
    """Eight threads, each with its own state of the same fluid, get the serial results."""
    expected = _densities(fluid, points, "update")
    with ThreadPoolExecutor(max_workers=8) as pool:
        results = list(pool.map(lambda _: _densities(fluid, points, "update_nogil"), range(16)))
    assert all(r == expected for r in results)


def test_a_failed_update_nogil_raises_value_error():
    state = _state("R134a")
    with pytest.raises(ValueError):
        state.update_nogil(CoolProp.PT_INPUTS, -1.0, 300.0)
