import unittest

import numpy as np
import pytest

from CoolProp.CoolProp import PropsSI

_FLUID = 'Water'
# scalar, list, and ndarray temperature inputs (exercises the scalar + vectorized
# PropsSI paths). Computed at import via the 2-arg trivial PropsSI form.
_TVALS = [
    0.5 * PropsSI(_FLUID, 'Tmin') + 0.5 * PropsSI(_FLUID, 'Tcrit'),
    [PropsSI(_FLUID, 'Tmin') + 1e-5, PropsSI(_FLUID, 'Tcrit') - 1e-5],
    np.linspace(PropsSI(_FLUID, 'Tmin') + 1e-5, PropsSI(_FLUID, 'Tcrit') - 1e-5, 30),
]


@pytest.mark.parametrize("Tvals", _TVALS, ids=["scalar", "list", "ndarray"])
def test_input_types(Tvals):
    """Scalar / list / ndarray temperature inputs are all accepted."""
    PropsSI('P', 'T', Tvals, 'Q', 0, _FLUID)


class PropsFailures(unittest.TestCase):
    def testUnmatchedLengths(self):
        # mismatched array lengths -> TypeError (matches the legacy wrapper)
        self.assertRaises(TypeError, PropsSI, 'P', 'T', [280, 290, 300], 'Q', [0, 1], 'R134a')

    def testMatrix(self):
        # multi-dimensional input -> ValueError "not one-dimensional"
        # (the legacy wrapper rejected this; the original assertion of TypeError
        # here never matched legacy and is corrected to ValueError).
        self.assertRaises(
            ValueError,
            PropsSI,
            'P',
            'T',
            np.array([280, 290, 300, 280, 290, 300]).reshape(2, 3),
            'Q',
            np.array([0, 0.5, 1, 0.0, 0.5, 1]).reshape(2, 3),
            'R134a',
        )
