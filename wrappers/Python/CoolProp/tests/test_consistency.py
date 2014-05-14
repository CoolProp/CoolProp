from __future__ import division, print_function
import CoolProp
from CoolProp import unit_systems_constants
from CoolProp.CoolProp import Props, get_REFPROPname, IsFluidType, set_standard_unit_system
import CoolProp.CoolProp as CP
import numpy as np

modes = []
modes.append('pure')
#modes.append('pseudo-pure')

# Check if REFPROP is supported, the Props call should work without throwing exception if it is supported
## try:
##     Props('D','T',300,'Q',1,'REFPROP-Propane')
##     modes.append('REFPROP')
## except ValueError:
##     pass
            
twophase_inputs = [('T','D'),('T','Q'),('P','Q'),('P','H'),('P','S'),('P','D'),('T','S'),('H','S')] #
singlephase_inputs = [('T','D'),('T','P'),('P','H'),('P','S'),('P','D'),('H','S'),('T','S')]
singlephase_outputs = ['T','P','H','S','A','O','C','G','V','L','C0','U']

## def test_subcrit_singlephase_consistency():
##     for Fluid in sorted(CoolProp.__fluids__):
##         T = (Props(Fluid,'Tmin')+Props(Fluid,'Tcrit'))/2.0
##         for mode in modes:
##             rhoL = Props('D','T',T,'Q',0,Fluid)
##             rhoV = Props('D','T',T,'Q',1,Fluid)
##             for rho in [rhoL+0.1, rhoV*0.9]:
##                 for inputs in singlephase_inputs:
##                     for unit_system in ['SI','kSI']:
##                         yield check_consistency,Fluid,mode,unit_system,T,rho,inputs
                    
def test_subcrit_twophase_consistency():
    
    for Fluid in reversed(sorted(CoolProp.__fluids__)):
        Tmin = Props(Fluid,'Tmin')
        Tcrit = Props(Fluid,'Tcrit')
        for T in [Tmin + 1, (Tmin+Tcrit)/2.0, 0.95*Tcrit]:
            for mode in modes:
                rhoL = Props('D','T',T,'Q',0,Fluid)
                rhoV = Props('D','T',T,'Q',1,Fluid)
                for Q in [0.0, 0.5, 1.0]:
                    rho = 1/((1-Q)/rhoL+Q/rhoV)
                    for inputs in twophase_inputs:
                        for unit_system in ['kSI','SI']:
                            yield check_consistency,Fluid,mode,unit_system, T,rho,inputs

def check_consistency(Fluid,mode,unit_system,T,rho,inputs):
        
    if unit_system == 'SI':
        set_standard_unit_system(unit_systems_constants.UNIT_SYSTEM_SI)
    elif unit_system == 'kSI':
        set_standard_unit_system(unit_systems_constants.UNIT_SYSTEM_KSI)
    else:
        raise ValueError('invalid unit_system:'+str(unit_system))
    
    if get_REFPROPname(Fluid) == 'N/A':
        return
        
    if mode == 'REFPROP':
        Fluid = 'REFPROP-' + get_REFPROPname(Fluid)
    
    if mode == 'pure' and not IsFluidType(Fluid,'PureFluid'):
        return
        
    #  Evaluate the inputs; if inputs is ('T','P'), calculate the temperature and the pressure
    Input1 = Props(inputs[0],'T',T,'D',rho,Fluid)
    Input2 = Props(inputs[1],'T',T,'D',rho,Fluid)
    
    #  Evaluate using the inputs given --> back to T,rho
    TEOS = Props('T',inputs[0],Input1,inputs[1],Input2,Fluid)
    DEOS = Props('D',inputs[0],Input1,inputs[1],Input2,Fluid)
    
    print('T',inputs[0],Input1,inputs[1],Input2,Fluid)
    #  Check they are consistent
    if abs(TEOS -T) > 1e-1 or abs(DEOS/rho-1) > 0.05:
        raise AssertionError("{T:g} K {D:g} kg/m^3 inputs: \"D\",'{ins1:s}',{in1:.12g},'{ins2:s}',{in2:.12g},\"{fluid:s}\"  || T: {TEOS:g} D: {DEOS:g}".format(T = T, 
                                                                                                                       D = rho, 
                                                                                                                       TEOS = TEOS, 
                                                                                                                       DEOS = DEOS, 
                                                                                                                       inputs = str(inputs), 
                                                                                                                       in1 = Input1, 
                                                                                                                       in2 = Input2,
                                                                                                                       ins1 = inputs[0],
                                                                                                                       ins2 = inputs[1],
                                                                                                                       fluid = Fluid)
                            )
    
if __name__=='__main__':
    import nose
    nose.runmodule()
