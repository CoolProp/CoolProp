from __future__ import division, print_function
import numpy as np
from .DataObjects import CoefficientData, PureData


class NitrateSalt(PureData, CoefficientData):
    """
    Heat transfer fluid based on 60% NaNO3 and 40% KNO3
    """

    def __init__(self):
        CoefficientData.__init__(self)
        PureData.__init__(self)
        self.name = "NaK"
        self.description = "Nitrate salt, 0.6 NaNO3 and 0.4 KNO3"
        self.reference = "Zavoico2001"

        self.Tmin = 300 + 273.15
        self.Tmax = 600 + 273.15
        self.TminPsat = self.Tmax

        self.Tbase = 273.15

        #self.temperature.data         = self.getTrange()
        # self.concentration.data       = np.array([     0 ]) # mass fraction

        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.source = self.density.SOURCE_COEFFS
        self.density.coeffs = np.array([[2090], [-0.636]])

        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.specific_heat.coeffs = np.array([[1443], [+0.172]])

        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.conductivity.coeffs = np.array([[0.443], [+1.9e-4]])

        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_POLYNOMIAL
        self.viscosity.source = self.viscosity.SOURCE_COEFFS
        self.viscosity.coeffs = np.array([[22.714], [-0.120], [2.281 * 1e-4], [-1.474 * 1e-7]]) / 1e3


class FoodProtein(PureData, CoefficientData):

    def __init__(self):
        CoefficientData.__init__(self)
        PureData.__init__(self)
        self.name = "FoodProtein"
        self.description = "Food protein model from the 2006 ASHRAE Handbook based on data from Choi and Okos (1986)"
        self.reference = "ASHRAE2006"

        self.Tmin = -40 + 273.15
        self.Tmax = 150 + 273.15
        self.TminPsat = self.Tmax
        self.Tbase = 0.00 + 273.15

        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = np.array([[1.7881e-1], [1.1958e-3], [-2.7178e-6]])

        self.density.source = self.density.SOURCE_COEFFS
        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = np.array([[1.3299e3], [-5.1840e-1]])

        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = np.array([[2.0082], [1.2089e-3], [-1.3129e-6]]) * 1e3
        

class FoodFat(PureData, CoefficientData):

    def __init__(self):
        CoefficientData.__init__(self)
        PureData.__init__(self)
        self.name = "FoodFat"
        self.description = "Food fat model from the 2006 ASHRAE Handbook based on data from Choi and Okos (1986)"
        self.reference = "ASHRAE2006"

        self.Tmin = -40 + 273.15
        self.Tmax = 150 + 273.15
        self.TminPsat = self.Tmax
        self.Tbase = 0.00 + 273.15

        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = np.array([[1.8071e-1], [-2.7604e-4], [-1.7749e-7]])
                                                                         
        self.density.source = self.density.SOURCE_COEFFS
        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = np.array([[9.2559e2], [-4.1757e-1]])

        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = np.array([[1.9842], [1.4733e-3], [-4.8008e-6]]) * 1e3


class FoodCarbohydrate(PureData, CoefficientData):

    def __init__(self):
        CoefficientData.__init__(self)
        PureData.__init__(self)
        self.name = "FoodCarbohydrate"
        self.description = "Food carbohydrate model from the 2006 ASHRAE Handbook based on data from Choi and Okos (1986)"
        self.reference = "ASHRAE2006"

        self.Tmin = -40 + 273.15
        self.Tmax = 150 + 273.15
        self.TminPsat = self.Tmax
        self.Tbase = 0.00 + 273.15

        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = np.array([[2.0141e-1], [1.3874e-3], [-4.3312e-6]])

        self.density.source = self.density.SOURCE_COEFFS
        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = np.array([[1.5991e3], [-3.1046e-1]])

        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = np.array([[1.5488], [1.9625e-3], [-5.9399e-6]]) * 1e3


class FoodFiber(PureData, CoefficientData):

    def __init__(self):
        CoefficientData.__init__(self)
        PureData.__init__(self)
        self.name = "FoodFiber"
        self.description = "Food fiber model from the 2006 ASHRAE Handbook based on data from Choi and Okos (1986)"
        self.reference = "ASHRAE2006"

        self.Tmin = -40 + 273.15
        self.Tmax = 150 + 273.15
        self.TminPsat = self.Tmax
        self.Tbase = 0.00 + 273.15

        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = np.array([[1.8331e-1], [1.2497e-3], [-3.1683e-6]])

        self.density.source = self.density.SOURCE_COEFFS
        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = np.array([[1.3115e3], [-3.6589e-1]])

        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = np.array([[1.8459], [1.8306e-3], [-4.6509e-6]]) * 1e3
                                                                      

class FoodAsh(PureData, CoefficientData):

    def __init__(self):
        CoefficientData.__init__(self)
        PureData.__init__(self)
        self.name = "FoodAsh"
        self.description = "Food ash model from the 2006 ASHRAE Handbook based on data from Choi and Okos (1986)"
        self.reference = "ASHRAE2006"

        self.Tmin = -40 + 273.15
        self.Tmax = 150 + 273.15
        self.TminPsat = self.Tmax
        self.Tbase = 0.00 + 273.15

        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = np.array([[3.2962e-1], [1.4011e-3], [-2.9069e-6]])

        self.density.source = self.density.SOURCE_COEFFS
        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = np.array([[2.4238e3], [-2.8063e-1]])

        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = np.array([[1.0926], [1.8896e-3], [-3.6817e-6]]) * 1e3


class FoodWater(PureData, CoefficientData):
    
    def __init__(self):
        CoefficientData.__init__(self)
        PureData.__init__(self)
        self.name = "FoodWater"
        self.description = "Food water model from the 2006 ASHRAE Handbook based on data from Choi and Okos (1986)"
        self.reference = "ASHRAE2006"

        self.Tmin = -40 + 273.15
        self.Tmax = 150 + 273.15
        self.TminPsat = self.Tmax
        self.Tbase = 0.00 + 273.15

        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = np.array([[5.7109e-1], [1.7625e-3], [-6.7036e-6]])

        self.density.source = self.density.SOURCE_COEFFS
        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = np.array([[9.9718e2], [3.1439e-3], [-3.7574e-3]])

        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = np.array([[4.1289], [-9.0864e-5], [5.4731e-6]]) * 1e3


class FoodIce(PureData, CoefficientData):
    
    def __init__(self):
        CoefficientData.__init__(self)
        PureData.__init__(self)
        self.name = "FoodIce"
        self.description = "Food ice model from the 2006 ASHRAE Handbook based on data from Choi and Okos (1986)"
        self.reference = "ASHRAE2006"

        self.Tmin = -40 + 273.15
        self.Tmax = 150 + 273.15
        self.TminPsat = self.Tmax
        self.Tbase = 0.00 + 273.15

        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = np.array([[2.2196], [-6.2489e-3], [1.0154e-4]])

        self.density.source = self.density.SOURCE_COEFFS
        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = np.array([[9.1689e2], [-1.3071e-1]])

        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = np.array([[2.0623], [6.0769e-3]]) * 1e3
