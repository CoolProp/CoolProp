"""
Some fluids that are part of CoolProp 4 or another
software and that have to be reimplemented for
CoolProp 5.

The most famous one is probably the aqueous Lithium
bromide solution. This one is going to be transferred
from the equations from the publication to the standard
parameter form.

"""
from __future__ import division, print_function
import numpy as np
from .DataObjects import DigitalData, PureData


class HyCool20(PureData, DigitalData):
    def __init__(self):
        DigitalData.__init__(self)
        PureData.__init__(self)

        self.name = "HY20"
        self.description = "HYCOOL 20, Potassium formate"
        self.reference = "Hydro2000"

        self.Tmax = 50 + 273.15
        self.Tmin = -20 + 273.15
        self.TminPsat = self.Tmax
        self.Tbase = 0.00 + 273.15
        self.temperature.data = self.getTrange()

        self.density.source = self.density.SOURCE_COEFFS
        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = np.array([[1202.2], [-0.42918]])

        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = np.array([[2.955], [0.0023]]) * 1e3

        key = 'Cond'

        def funcCond(T, x):
            T = (T - self.Tbase)
            if T <= 20: return 0.001978 * T + 0.5172
            else: return 0.001005 * T + 0.5368
        self.conductivity.xData, self.conductivity.yData, self.conductivity.data = self.getArray(dataID=key, func=funcCond, x_in=self.temperature.data, y_in=self.concentration.data, DEBUG=self.conductivity.DEBUG)
        self.conductivity.source = self.conductivity.SOURCE_EQUATION
        funcCond = None

        key = 'Mu'

        def funcMu(T, x):
            T = (T - self.Tbase)
            if T <= 20: mPas = 0.07190 * np.exp(524.75 / (T + 142.05))
            else: mPas = T * (0.0005524 * T - 0.06281) + 2.8536
            return mPas / 1e3
        self.viscosity.xData, self.viscosity.yData, self.viscosity.data = self.getArray(dataID=key, func=funcMu, x_in=self.temperature.data, y_in=self.concentration.data, DEBUG=self.viscosity.DEBUG)
        self.viscosity.source = self.viscosity.SOURCE_EQUATION
        funcMu = None


class HyCool30(PureData, DigitalData):
    def __init__(self):
        DigitalData.__init__(self)
        PureData.__init__(self)

        self.name = "HY30"
        self.description = "HyCool 30, Potassium formate"
        self.reference = "Hydro2000"

        self.Tmax = 50 + 273.15
        self.Tmin = -30 + 273.15
        self.TminPsat = self.Tmax
        self.Tbase = 0.00 + 273.15
        self.temperature.data = self.getTrange()

        self.density.source = self.density.SOURCE_COEFFS
        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = np.array([[1257.5], [-0.475350]])

        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = np.array([[2.783], [0.0023]]) * 1e3

        key = 'Cond'

        def funcCond(T, x):
            T = (T - self.Tbase)
            if T <= 20: return 0.001840 * T + 0.4980
            else: return 0.001000 * T + 0.5140
        self.conductivity.xData, self.conductivity.yData, self.conductivity.data = self.getArray(dataID=key, func=funcCond, x_in=self.temperature.data, y_in=self.concentration.data, DEBUG=self.conductivity.DEBUG)
        self.conductivity.source = self.conductivity.SOURCE_EQUATION
        funcCond = None

        key = 'Mu'

        def funcMu(T, x):
            T = (T - self.Tbase)
            if T <= 20: mPas = 0.11100 * np.exp(408.17 / (T + 123.15))
            else: mPas = T * (0.000295 * T - 0.0441) + 2.6836
            return mPas / 1e3
        self.viscosity.xData, self.viscosity.yData, self.viscosity.data = self.getArray(dataID=key, func=funcMu, x_in=self.temperature.data, y_in=self.concentration.data, DEBUG=self.viscosity.DEBUG)
        self.viscosity.source = self.viscosity.SOURCE_EQUATION
        funcMu = None


class HyCool40(PureData, DigitalData):
    def __init__(self):
        DigitalData.__init__(self)
        PureData.__init__(self)

        self.name = "HY40"
        self.description = "HyCool 40, Potassium formate"
        self.reference = "Hydro2000"

        self.Tmax = 20 + 273.15
        self.Tmin = -40 + 273.15
        self.TminPsat = self.Tmax
        self.Tbase = 0.00 + 273.15
        self.temperature.data = self.getTrange()

        self.density.source = self.density.SOURCE_COEFFS
        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = np.array([[1304.5], [-0.512290]])

        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = np.array([[2.646], [0.0023]]) * 1e3

        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = np.array([[0.4826], [0.001730]])

        key = 'Mu'

        def funcMu(T, x):
            T = (T - self.Tbase)
            mPas = 0.07830 * np.exp(498.13 / (T + 130.25))
            return mPas / 1e3
        self.viscosity.xData, self.viscosity.yData, self.viscosity.data = self.getArray(dataID=key, func=funcMu, x_in=self.temperature.data, y_in=self.concentration.data, DEBUG=self.viscosity.DEBUG)
        self.viscosity.source = self.viscosity.SOURCE_EQUATION
        funcMu = None


class HyCool45(PureData, DigitalData):
    def __init__(self):
        DigitalData.__init__(self)
        PureData.__init__(self)

        self.name = "HY45"
        self.description = "HyCool 45, Potassium formate"
        self.reference = "Hydro2000"

        self.Tmax = 20 + 273.15
        self.Tmin = -45 + 273.15
        self.TminPsat = self.Tmax
        self.Tbase = 0.00 + 273.15
        self.temperature.data = self.getTrange()

        self.density.source = self.density.SOURCE_COEFFS
        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = np.array([[1328.7], [-0.530754]])

        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = np.array([[2.578], [0.0023]]) * 1e3

        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = np.array([[0.4750], [0.001674]])

        key = 'Mu'

        def funcMu(T, x):
            T = (T - self.Tbase)
            mPas = 0.08990 * np.exp(479.09 / (T + 126.55))
            return mPas / 1e3
        self.viscosity.xData, self.viscosity.yData, self.viscosity.data = self.getArray(dataID=key, func=funcMu, x_in=self.temperature.data, y_in=self.concentration.data, DEBUG=self.viscosity.DEBUG)
        self.viscosity.source = self.viscosity.SOURCE_EQUATION
        funcMu = None


class HyCool50(PureData, DigitalData):
    def __init__(self):
        DigitalData.__init__(self)
        PureData.__init__(self)

        self.name = "HY50"
        self.description = "HyCool 50, Potassium formate"
        self.reference = "Hydro2000"

        self.Tmax = 20 + 273.15
        self.Tmin = -50 + 273.15
        self.TminPsat = self.Tmax
        self.Tbase = 0.00 + 273.15
        self.temperature.data = self.getTrange()

        self.density.source = self.density.SOURCE_COEFFS
        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = np.array([[1359.0], [-0.552300]])

        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = np.array([[2.498], [0.0023]]) * 1e3

        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = np.array([[0.4660], [0.001610]])

        key = 'Mu'

        def funcMu(T, x):
            T = (T - self.Tbase)
            res = 0.0491 * np.exp(581.12 / (T + 129.05))
            if T > -10: mPas = res + 0.2
            else: mPas = res
            return mPas / 1e3
        self.viscosity.xData, self.viscosity.yData, self.viscosity.data = self.getArray(dataID=key, func=funcMu, x_in=self.temperature.data, y_in=self.concentration.data, DEBUG=self.viscosity.DEBUG)
        self.viscosity.source = self.viscosity.SOURCE_EQUATION
        funcMu = None
