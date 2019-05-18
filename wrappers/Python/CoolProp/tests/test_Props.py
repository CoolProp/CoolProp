import unittest
from CoolProp.CoolProp import PropsSI
import CoolProp
import numpy as np


def test_input_types():
    for Fluid in ['Water']:
        for Tvals in [0.5 * PropsSI(Fluid, 'Tmin') + 0.5 * PropsSI(Fluid, 'Tcrit'),
                      [PropsSI(Fluid, 'Tmin') + 1e-5, PropsSI(Fluid, 'Tcrit') - 1e-5],
                      np.linspace(PropsSI(Fluid, 'Tmin') + 1e-5, PropsSI(Fluid, 'Tcrit') - 1e-5, 30)
                      ]:
            yield check_type, Fluid, Tvals


def check_type(fluid, Tvals):
    PropsSI('P', 'T', Tvals, 'Q', 0, fluid)


class PropsFailures(unittest.TestCase):
    def testUnmatchedLengths(self):
        self.assertRaises(TypeError, PropsSI, 'P', 'T', [280, 290, 300], 'Q', [0, 1], 'R134a')

    def testMatrix(self):
        self.assertRaises(TypeError, PropsSI, 'P', 'T', np.array([280, 290, 300, 280, 290, 300]).reshape(2, 3), 'Q', np.array([0, 0.5, 1, 0.0, 0.5, 1]).reshape(2, 3), 'R134a')


if __name__ == '__main__':
    import nose
    nose.runmodule()
