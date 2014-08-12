import CoolProp.CoolProp as CP
import CoolProp.unit_systems_constants
from CoolProp import param_constants
from CoolProp.State import State

def test_global_param_strings():
    for p in ['version','errstring','gitrevision','FluidsList']:
        yield check_global, p

def check_global(p):
    CP.get_global_param_string(p)

def test_fluid_param_strings():
    for fluid in CoolProp.__fluids__:
        for p in ['aliases','CAS','ASHRAE34','REFPROPName','TTSE_mode']:
            yield check_fluid, fluid, p

def check_fluid(fluid, p):
    CP.get_fluid_param_string(fluid,p)
        
if __name__=='__main__':
    import nose
    nose.runmodule()