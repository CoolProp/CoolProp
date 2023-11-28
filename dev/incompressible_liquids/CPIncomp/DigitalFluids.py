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

try:
    import CoolProp.CoolProp as CP

    def get_coolprop_class_Tp(fluid_name: str, fluid_desc: str, backend_str: str, fluid_str: str, p_ref: float, T_min: float, T_max: float):

        class NewClass(PureData, DigitalData):
            def __init__(self):
                DigitalData.__init__(self)
                PureData.__init__(self)

                self.name = fluid_name
                self.description = fluid_desc

                cp_state = CP.AbstractState(backend_str, fluid_str)

                references = []
                references.append(CP.get_BibTeXKey(fluid_str, "EOS"))
                references.append(CP.get_BibTeXKey(fluid_str, "CP0"))
                references.append(CP.get_BibTeXKey(fluid_str, "VISCOSITY"))
                references.append(CP.get_BibTeXKey(fluid_str, "CONDUCTIVITY"))
                # references.append(CP.get_BibTeXKey(fluid_str, "ECS_LENNARD_JONES"))
                # references.append(CP.get_BibTeXKey(fluid_str, "ECS_FITS"))
                # references.append(CP.get_BibTeXKey(fluid_str, "SURFACE_TENSION"))
                # references.append(CP.get_BibTeXKey(fluid_str, "MELTING_LINE"))

                self.reference = "; ".join(references)

                self.Tmax = min(cp_state.Tmax(), T_max)
                self.Tmin = max(cp_state.Tmin(), T_min)

                if p_ref < cp_state.p_critical():# and p_ref > cp_state.p_triple():
                    cp_state.update(CP.PQ_INPUTS, p_ref, 0.0)
                    # prefer the liquid phase
                    if self.Tmin < cp_state.T():
                        self.Tmax = min(self.Tmax, cp_state.T() - 1e-3)
                    elif self.Tmax > cp_state.T():
                        self.Tmin = max(self.Tmin, cp_state.T() + 1e-3)

                self.TminPsat = self.Tmax
                self.Tbase = (self.Tmax + self.Tmin) / 2.0
                self.temperature.data = self.getTrange()
                self.concentration.data = self.getxrange()

                def funcRho(T, x):
                    cp_state.update(CP.PT_INPUTS, p_ref, T)
                    return cp_state.rhomass()
                
                self.density.xData, self.density.yData, self.density.data = self.getArray(dataID='Rho', func=funcRho, x_in=self.temperature.data, y_in=self.concentration.data, DEBUG=self.density.DEBUG)
                self.density.source = self.density.SOURCE_EQUATION


                def funcCp(T, x):
                    cp_state.update(CP.PT_INPUTS, p_ref, T)
                    return cp_state.cpmass()
                
                self.specific_heat.xData, self.specific_heat.yData, self.specific_heat.data = self.getArray(dataID='Cp', func=funcCp, x_in=self.temperature.data, y_in=self.concentration.data, DEBUG=self.specific_heat.DEBUG)
                self.specific_heat.source = self.specific_heat.SOURCE_EQUATION

                try:
                    def funcMu(T, x):
                        cp_state.update(CP.PT_INPUTS, p_ref, T)
                        return cp_state.viscosity()
                    
                    self.viscosity.xData, self.viscosity.yData, self.viscosity.data = self.getArray(dataID='Mu', func=funcMu, x_in=self.temperature.data, y_in=self.concentration.data, DEBUG=self.viscosity.DEBUG)
                    self.viscosity.source = self.viscosity.SOURCE_EQUATION
                except:
                    pass

                
                try:
                    def funcCond(T, x):
                        cp_state.update(CP.PT_INPUTS, p_ref, T)
                        return cp_state.conductivity()
                    
                    self.conductivity.xData, self.conductivity.yData, self.conductivity.data = self.getArray(dataID='Cond', func=funcCond, x_in=self.temperature.data, y_in=self.concentration.data, DEBUG=self.conductivity.DEBUG)
                    self.conductivity.source = self.conductivity.SOURCE_EQUATION
                except:
                    pass

        return NewClass

    class Air(get_coolprop_class_Tp(fluid_name="Air", fluid_desc="Air, gaseous phase at 1 atm (101325 Pa)", backend_str="HEOS", fluid_str="Air", p_ref=101325.0, T_min=-75 + 273.15, T_max=250 + 273.15)):
        pass

    class Ethanol(get_coolprop_class_Tp(fluid_name="Ethanol", fluid_desc="Ethanol, liquid phase at 10 bar", backend_str="HEOS", fluid_str="Ethanol", p_ref=10e5, T_min=-75 + 273.15, T_max=250 + 273.15)):
        pass

    class Acetone(get_coolprop_class_Tp(fluid_name="Acetone", fluid_desc="Acetone, liquid phase at 10 bar", backend_str="HEOS", fluid_str="Acetone", p_ref=10e5, T_min=-75 + 273.15, T_max=250 + 273.15)):
        pass

    class Hexane(get_coolprop_class_Tp(fluid_name="Hexane",   fluid_desc="Hexane, liquid phase at 10 bar",  backend_str="HEOS", fluid_str="Hexane", p_ref=10e5, T_min=-75 + 273.15, T_max=250 + 273.15)):
        pass




            
except ImportError:
    # Do not handle cases where CoolProp is not installed
    pass