import numpy as np
from CPIncomp.CoefficientObjects import CoefficientData

class DEBLiquidClass(CoefficientData):
    """ 
    Pure fluid according to Melinder's book
    """ 
    def __init__(self):
        CoefficientData.__init__(self) 
        self.name = "DEB"
        self.description = "Diethylbenzene mixture - Dowtherm J Dow Chemical Co."
        self.reference = "Melinder-BOOK-2010"

        self.Tmin = -80.0 + 273.15;
        self.Tmax = 100.0 + 273.15;
        self.TminPsat =  self.Tmax

        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        _,_,self.density.coeffs = self.density.shapeArray(np.array([1076.5,-0.731182]))

        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        _,_,self.specific_heat.coeffs = self.specific_heat.shapeArray(np.array([999.729,2.87576]))

        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
        _,_,self.viscosity.coeffs = self.viscosity.shapeArray(np.array([3.5503,-0.0566396,7.03331e-05]))

        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        _,_,self.conductivity.coeffs = self.conductivity.shapeArray(np.array([0.000189132,-2.06364e-07]))


class HCMLiquidClass(CoefficientData):
    """ 
    Pure fluid according to Melinder's book
    """ 
    def __init__(self):
        CoefficientData.__init__(self) 
        self.name = "HCM"
        self.description = "Hydrocarbon mixture (synthetic) - Therminol D12 (Gilotherm D12) Solutia"
        self.reference = "Melinder-BOOK-2010"

        self.Tmin = -80.0 + 273.15;
        self.Tmax = 100.0 + 273.15;
        self.TminPsat =  self.Tmax

        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        _,_,self.density.coeffs = self.density.shapeArray(np.array([971.725,-0.718788]))

        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        _,_,self.specific_heat.coeffs = self.specific_heat.shapeArray(np.array([844.023,4.31212]))

        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
        _,_,self.viscosity.coeffs = self.viscosity.shapeArray(np.array([18.3237,-0.14706,0.000209096]))

        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        _,_,self.conductivity.coeffs = self.conductivity.shapeArray(np.array([0.000153716,-1.51212e-07]))


class HFELiquidClass(CoefficientData):
    """ 
    Pure fluid according to Melinder's book
    """ 
    def __init__(self):
        CoefficientData.__init__(self) 
        self.name = "HFE"
        self.description = "Hydrofluoroether - HFE-7100 3M Novec"
        self.reference = "Melinder-BOOK-2010"

        self.Tmin = -80.0 + 273.15;
        self.Tmax = 100.0 + 273.15;
        self.TminPsat =  self.Tmax

        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        _,_,self.density.coeffs = self.density.shapeArray(np.array([1822.37,-0.918485]))

        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        _,_,self.specific_heat.coeffs = self.specific_heat.shapeArray(np.array([871.834,858788]))

        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
        _,_,self.viscosity.coeffs = self.viscosity.shapeArray(np.array([-4.22878,-0.0114765,7.39823e-06]))

        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        _,_,self.conductivity.coeffs = self.conductivity.shapeArray(np.array([9.92958e-05,-8.33333e-08]))


class PMS1LiquidClass(CoefficientData):
    """ 
    Pure fluid according to Melinder's book
    """ 
    def __init__(self):
        CoefficientData.__init__(self) 
        self.name = "PMS1"
        self.description = "Polydimethylsiloxan 1. - Baysilone KT3"
        self.reference = "Melinder-BOOK-2010"

        self.Tmin = -80.0 + 273.15;
        self.Tmax = 100.0 + 273.15;
        self.TminPsat =  self.Tmax

        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        _,_,self.density.coeffs = self.density.shapeArray(np.array([1172.35,-0.9025]))

        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        _,_,self.specific_heat.coeffs = self.specific_heat.shapeArray(np.array([1223.69,1.48417]))

        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
        _,_,self.viscosity.coeffs = self.viscosity.shapeArray(np.array([6.36183,-0.0636352,7.51428e-05]))

        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        _,_,self.conductivity.coeffs = self.conductivity.shapeArray(np.array([0.000207526,-2.84167e-07]))


class PMS2LiquidClass(CoefficientData):
    """ 
    Pure fluid according to Melinder's book
    """ 
    def __init__(self):
        CoefficientData.__init__(self) 
        self.name = "PMS2"
        self.description = "Polydimethylsiloxan 2. - Syltherm XLT Dow Corning Co."
        self.reference = "Melinder-BOOK-2010"

        self.Tmin = -80.0 + 273.15;
        self.Tmax = 100.0 + 273.15;
        self.TminPsat =  self.Tmax

        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        _,_,self.density.coeffs = self.density.shapeArray(np.array([1155.94,-1.02576]))

        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        _,_,self.specific_heat.coeffs = self.specific_heat.shapeArray(np.array([1153.55,2.10788]))

        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
        _,_,self.viscosity.coeffs = self.viscosity.shapeArray(np.array([5.66926,-0.065582,8.09988e-05]))

        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        _,_,self.conductivity.coeffs = self.conductivity.shapeArray(np.array([0.000172305,-2.11212e-07]))


class SABLiquidClass(CoefficientData):
    """ 
    Pure fluid according to Melinder's book
    """ 
    def __init__(self):
        CoefficientData.__init__(self) 
        self.name = "SAB"
        self.description = "Synthetic alkyl benzene - Marlotherm X"
        self.reference = "Melinder-BOOK-2010"

        self.Tmin = -80.0 + 273.15;
        self.Tmax = 100.0 + 273.15;
        self.TminPsat =  self.Tmax

        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        _,_,self.density.coeffs = self.density.shapeArray(np.array([1102.34,-0.801667]))

        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        _,_,self.specific_heat.coeffs = self.specific_heat.shapeArray(np.array([1360.94,1.51667]))

        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
        _,_,self.viscosity.coeffs = self.viscosity.shapeArray(np.array([5.21288,-0.0665792,8.5066e-05]))

        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        _,_,self.conductivity.coeffs = self.conductivity.shapeArray(np.array([0.000208374,-2.61667e-07]))


class HCBLiquidClass(CoefficientData):
    """ 
    Pure fluid according to Melinder's book
    """ 
    def __init__(self):
        CoefficientData.__init__(self) 
        self.name = "HCB"
        self.description = "Hydrocarbon blend - Dynalene MV"
        self.reference = "Melinder-BOOK-2010"

        self.Tmin = -80.0 + 273.15;
        self.Tmax = 100.0 + 273.15;
        self.TminPsat =  self.Tmax

        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        _,_,self.density.coeffs = self.density.shapeArray(np.array([1071.78,-0.772024]))

        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        _,_,self.specific_heat.coeffs = self.specific_heat.shapeArray(np.array([761.393,3.52976]))

        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
        _,_,self.viscosity.coeffs = self.viscosity.shapeArray(np.array([7.16819,-0.0863212,0.000130604]))

        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        _,_,self.conductivity.coeffs = self.conductivity.shapeArray(np.array([0.000203186,-2.3869e-07]))


class TCOLiquidClass(CoefficientData):
    """ 
    Pure fluid according to Melinder's book
    """ 
    def __init__(self):
        CoefficientData.__init__(self) 
        self.name = "TCO"
        self.description = "Terpene from citrus oils - d-Limonene"
        self.reference = "Melinder-BOOK-2010"

        self.Tmin = -80.0 + 273.15;
        self.Tmax = 100.0 + 273.15;
        self.TminPsat =  self.Tmax

        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        _,_,self.density.coeffs = self.density.shapeArray(np.array([1071.02,-0.778166]))

        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        _,_,self.specific_heat.coeffs = self.specific_heat.shapeArray(np.array([223.775,5.2159]))

        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
        _,_,self.viscosity.coeffs = self.viscosity.shapeArray(np.array([-3.47971,-0.0107031,1.14086e-06]))

        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        _,_,self.conductivity.coeffs = self.conductivity.shapeArray(np.array([0.000174156,-1.85052e-07]))
