from __future__ import division, print_function
import CoolProp.CoolProp as CP
import CoolProp.unit_systems_constants
from CoolProp import param_constants
from CoolProp.State import State

S = State('R134a',dict(T=300,D=1))

# factor is SI / kSI value
State_listing = [('T',1),
                 ('get_speed_sound',1),
                 ('rho',1),
                 ('p',1000),
                 ('h',1000),
                 ('s',1000),
                 ('cv',1000),
                 ('cp',1000),
                 ('visc',1),
                 ('k',1000),
                ]
    
def test_State():    
    for parameter,SI_over_kSI in State_listing:
        
        CP.set_standard_unit_system(CoolProp.unit_systems_constants.UNIT_SYSTEM_SI)
        val_SI = getattr(S, parameter)
        
        CP.set_standard_unit_system(CoolProp.unit_systems_constants.UNIT_SYSTEM_KSI)
        val_kSI = getattr(S, parameter)
        
        yield check, val_SI, val_kSI, SI_over_kSI
    
Props_listing = [('T',1),
                 ('A',1),
           ('D',1),
           ('P',1000),
           ('H',1000),
           ('S',1000),
           ('C',1000),
           ('C0',1000),
           ('O',1000),
           ('V',1),
           ('L',1000),
           ]
def test_PROPS():
    for parameter, SI_over_kSI in Props_listing:
        yield check_Props, parameter, SI_over_kSI
        
def check_Props(parameter, SI_over_kSI):
    
    CP.set_standard_unit_system(CoolProp.unit_systems_constants.UNIT_SYSTEM_SI)
    val_SI = CP.Props(parameter,'T',300.0,'D',1.0,'R134a')
    
    CP.set_standard_unit_system(CoolProp.unit_systems_constants.UNIT_SYSTEM_KSI)
    val_kSI = CP.Props(parameter,'T',300.0,'D',1.0,'R134a')
    
    try:
        val_SI = val_SI()
        val_kSI = val_kSI()
    except:
        pass
            
    print(val_SI,val_kSI, val_SI/val_kSI - SI_over_kSI)
    if abs(val_SI/val_kSI - SI_over_kSI) > 1e-12:
        raise ValueError(val_SI/val_kSI-SI_over_kSI)
        
State_Props_listing = [(param_constants.iT,1),
                       (param_constants.iA,1),
                       (param_constants.iD,1),
                       (param_constants.iP,1000),
                       (param_constants.iH,1000),
                       (param_constants.iS,1000),
                       (param_constants.iC,1000),
                       (param_constants.iC0,1000),
                       (param_constants.iO,1000),
                       (param_constants.iV,1),
                       (param_constants.iL,1000),
                       ]
def test_State_PROPS():
    for parameter, SI_over_kSI in State_Props_listing:
        
        CP.set_standard_unit_system(CoolProp.unit_systems_constants.UNIT_SYSTEM_SI)
        val_SI = S.Props(parameter)
        
        CP.set_standard_unit_system(CoolProp.unit_systems_constants.UNIT_SYSTEM_KSI)
        val_kSI = S.Props(parameter)
        
        yield check, val_SI, val_kSI, SI_over_kSI
    
def check(val_SI, val_kSI, SI_over_kSI):
    
    try:
        val_SI = val_SI()
        val_kSI = val_kSI()
    except:
        pass
            
    if abs(val_SI/val_kSI - SI_over_kSI) > 1e-12:
        raise ValueError(val_SI/val_kSI-SI_over_kSI)
        
if __name__=='__main__':
    import nose
    nose.runmodule()