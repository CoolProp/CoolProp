
from __future__ import print_function,division

import CoolProp
#from CoolProp.CoolProp import PropsSI

import timeit

def get_speed_data():
    
    H_TP = 350e3
    P_TP = 400e3

    H_SP = 250e3
    P_SP = 1000e3

    P_PT = 101325
    T_PT = 300

    fluid = 'R245fa'
    number = 10000
    repeat = 3
    version = CoolProp.__version__
    
    if int(CoolProp.__version__[0])>4:
        loaded = 5
        print("Loaded CoolProp version 5")
        from CoolProp.CoolProp import generate_update_pair,get_parameter_index,set_debug_level
        TTSE = CoolProp.AbstractState('TTSE&HEOS',fluid)
        BICUBIC = CoolProp.AbstractState('BICUBIC&HEOS',fluid)
        HEOS = CoolProp.AbstractState('HEOS',fluid)
        def two_phase_TTSE():
            TTSE.update(CoolProp.HmassP_INPUTS, H_TP, P_TP)
            TTSE.rhomolar()
        def single_phase_TTSE():
            TTSE.update(CoolProp.HmassP_INPUTS, H_SP, P_SP)
            TTSE.rhomolar()
        def single_phase_pT_TTSE():
            TTSE.update(CoolProp.PT_INPUTS, P_PT, T_PT)
            TTSE.rhomolar()
            
        def two_phase_BICUBIC():
            BICUBIC.update(CoolProp.HmassP_INPUTS, H_TP, P_TP)
            BICUBIC.rhomolar()
        def single_phase_BICUBIC():
            BICUBIC.update(CoolProp.HmassP_INPUTS, H_SP, P_SP)
            BICUBIC.rhomolar()
        def single_phase_pT_BICUBIC():
            BICUBIC.update(CoolProp.PT_INPUTS, P_PT, T_PT)
            BICUBIC.rhomolar()
            
        def two_phase_HEOS():
            HEOS.update(CoolProp.HmassP_INPUTS, H_TP, P_TP)
            HEOS.rhomolar()
        def single_phase_HEOS():
            HEOS.update(CoolProp.HmassP_INPUTS, H_SP, P_SP)
            HEOS.rhomolar()
        def single_phase_pT_HEOS():
            HEOS.update(CoolProp.PT_INPUTS, P_PT, T_PT)
            HEOS.rhomolar()
        
    else:
        loaded = 4
        print("Loaded CoolProp version 4")
        #from CoolProp.CoolProp import set_debug_level,set_standard_unit_system,enable_TTSE_LUT,disable_TTSE_LUT
        CoolProp.CoolProp.set_standard_unit_system(CoolProp.UNIT_SYSTEM_SI)
        state = CoolProp.State.State(fluid,{"H":H_TP*2,"P":P_TP})
        def two_phase_HP():
            state.update({"H":H_TP,"P":P_TP})
            state.get_rho()
        def single_phase_HP():
            state.update({"H":H_SP,"P":P_SP})
            state.get_rho()
        def single_phase_PT():
            state.update({"P":P_PT,"T":T_PT})
            state.get_rho()

    if loaded==4:
        CoolProp.CoolProp.disable_TTSE_LUT(fluid)
        two_phase_hp_heos    = min(timeit.Timer(two_phase_HP).repeat(repeat=repeat, number=number))/number
        single_phase_hp_heos = min(timeit.Timer(single_phase_HP).repeat(repeat=repeat, number=number))/number
        single_phase_pt_heos = min(timeit.Timer(single_phase_PT).repeat(repeat=repeat, number=number))/number
        CoolProp.CoolProp.enable_TTSE_LUT(fluid)
        two_phase_hp_ttse    = min(timeit.Timer(two_phase_HP).repeat(repeat=repeat, number=number))/number
        single_phase_hp_ttse = min(timeit.Timer(single_phase_HP).repeat(repeat=repeat, number=number))/number
        single_phase_pt_ttse = min(timeit.Timer(single_phase_PT).repeat(repeat=repeat, number=number))/number
        CoolProp.CoolProp.disable_TTSE_LUT(fluid)
    elif loaded==5:
        two_phase_hp_heos    = min(timeit.Timer(two_phase_HEOS).repeat(repeat=repeat, number=number))/number
        single_phase_hp_heos = min(timeit.Timer(single_phase_HEOS).repeat(repeat=repeat, number=number))/number
        single_phase_pt_heos = min(timeit.Timer(single_phase_pT_HEOS).repeat(repeat=repeat, number=number))/number
        two_phase_hp_ttse    = min(timeit.Timer(two_phase_TTSE).repeat(repeat=repeat, number=number))/number
        single_phase_hp_ttse = min(timeit.Timer(single_phase_TTSE).repeat(repeat=repeat, number=number))/number
        single_phase_pt_ttse = min(timeit.Timer(single_phase_pT_TTSE).repeat(repeat=repeat, number=number))/number
        two_phase_hp_bicubic    = min(timeit.Timer(two_phase_BICUBIC).repeat(repeat=repeat, number=number))/number
        single_phase_hp_bicubic = min(timeit.Timer(single_phase_BICUBIC).repeat(repeat=repeat, number=number))/number
        single_phase_pt_bicubic = min(timeit.Timer(single_phase_pT_BICUBIC).repeat(repeat=repeat, number=number))/number
    else:
        raise ValueError("Unknown CoolProp version.")
        
    return locals()
    
if __name__ == '__main__':

    print("{0:15s}: {1:6d} calls, {2:2d} repetitions, CoolProp version {3:1d}".format(fluid,number,repeat,loaded))
    print("{0:15s}: 2P HEOS: {1:6.2f} us,  1P HEOS: {2:6.2f} us, PT HEOS: {3:6.2f} us".format(fluid,two_phase_hp_heos*1e6,single_phase_hp_heos*1e6,single_phase_pt_heos*1e6))
    print("{0:15s}: 2P TTSE: {1:6.2f} us,  1P TTSE: {2:6.2f} us, PT TTSE: {3:6.2f} us".format(fluid,two_phase_hp_ttse*1e6,single_phase_hp_ttse*1e6,single_phase_pt_ttse*1e6))
    print("{0:15s}: 2P BICU: {1:6.2f} us,  1P BICU: {2:6.2f} us, PT BICU: {3:6.2f} us".format(fluid,two_phase_hp_bicubic*1e6,single_phase_hp_bicubic*1e6,single_phase_pt_bicubic*1e6))
