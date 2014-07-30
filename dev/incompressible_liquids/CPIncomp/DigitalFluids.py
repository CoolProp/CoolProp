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
from CPIncomp.DataObjects import DigitalData  

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
        
        self.Tmin = 273.00
        self.Tmax = 500.00
        self.xmax = 0.75
        self.xmin = 0.00
        self.xid  = self.ifrac_mass
        self.TminPsat = self.Tmin
        
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
        self.density.source           = self.density.SOURCE_EQUATION
        
        key = 'C'
        def funcC(T,x):
            return CP.PropsSI(key,'T',T,'P',1e8,self.name+"-{0:.4f}%".format(x*100.0))
        self.specific_heat.data = self.getArray(funcC,key)
        self.specific_heat.source     = self.specific_heat.SOURCE_EQUATION
        
        key = 'Psat'
        def funcP(T,x):
            return CP.PropsSI(key,'T',T,'P',1e8,self.name+"-{0:.4f}%".format(x*100.0))
        self.saturation_pressure.data   = self.getArray(funcP,key)
        self.saturation_pressure.source = self.saturation_pressure.SOURCE_EQUATION

