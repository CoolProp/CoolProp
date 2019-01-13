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
