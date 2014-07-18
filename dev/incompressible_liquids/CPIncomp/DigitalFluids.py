"""
Some fluids that are part of CoolProp 4 or another 
software and that have to be reimplemented for 
CoolProp 5. 

The most famous one is probably the aqueous Lithium 
bromide solution. This one is going to be transferred 
from the equations from the publication to the standard
parameter form.

"""
import numpy as np
from CPIncomp.DataObjects import DigitalData
from CPIncomp.DataObjects import PureData


class NitrateSalt(PureData,DigitalData):
    """ 
    Heat transfer fluid based on 60% NaNO3 and 40% KNO3
    """
    def __init__(self):
        DigitalData.__init__(self) 
        PureData.__init__(self)
        self.name        = "NaK" 
        self.description = "NitrateSalt"
        self.reference   = "Solar Power Tower Design Basis Document,  Alexis B. Zavoico, Sandia Labs, USA"
        
        self.Tmin        = 300 + 273.15
        self.Tmax        = 600 + 273.15
        self.TminPsat    = self.Tmax 
        
        self.temperature.data         = self.getTrange()
        self.concentration.data       = np.array([     0 ]) # mass fraction
        
        def f_rho( T,x):
            return 2090 - 0.636 * (T-273.15)
        def f_cp(  T,x):
            return 1443 + 0.172 * (T-273.15)
        def f_mu(  T,x):
            return ( 22.714 - 0.120 * (T-273.15) + 2.281 * 1e-4 * (T-273.15)*(T-273.15) - 1.474 * 1e-7 * (T-273.15)*(T-273.15)*(T-273.15) )/1e3
        def f_lam( T,x):
            return 0.443 + 1.9e-4 * (T-273.15)

        self.density.data       = self.getArray(f_rho,'D')
        self.specific_heat.data = self.getArray(f_cp ,'C')
        self.viscosity.data     = self.getArray(f_mu ,'V')
        self.conductivity.data  = self.getArray(f_lam,'L')
        

class LiBrData(DigitalData):
    """ 
    Lithium Bromide solution from CoolProp 4
    Source: Patek et al.
    """ 
    def __init__(self):
        DigitalData.__init__(self) 
        
        import CoolProp.CoolProp as CP
        self.name = "LiBr"
        self.description = "Lithium-Bromide solution from Patek2006"
        self.reference = "Patek2006"
        
        self.Tmin = 273.00;
        self.Tmax = 500.00;
        self.xmax = 0.75
        self.xmin = 0.00
        self.TminPsat = self.Tmin;
        
        self.temperature.data         = self.getTrange()
        self.concentration.data       = self.getxrange()
        
#        data = [self.density.data,self.specific_heat.data,self.saturation_pressure.data]
#        keys = ["D",              "C",                    "Psat"]
#        
#        import os
#        for i in range(len(keys)):
#            def func(T,x):
#                return CP.PropsSI(keys[i],'T',T,'P',1e8,self.name+"-{0:.4f}%".format(x*100.0))
#            #if os.path.isfile(self.getFile(key)): os.remove(self.getFile(key))
#            data[i] = self.getArray(func,keys[i])

        key = 'D'
        def funcD(T,x):
            return CP.PropsSI(key,'T',T,'P',1e8,self.name+"-{0:.4f}%".format(x*100.0))
        self.density.data = self.getArray(funcD,key)
        
        key = 'C'
        def funcC(T,x):
            return CP.PropsSI(key,'T',T,'P',1e8,self.name+"-{0:.4f}%".format(x*100.0))
        self.specific_heat.data = self.getArray(funcC,key)
        
        key = 'Psat'
        def funcP(T,x):
            return CP.PropsSI(key,'T',T,'P',1e8,self.name+"-{0:.4f}%".format(x*100.0))
        self.saturation_pressure.data = self.getArray(funcP,key)

        

