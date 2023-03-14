
# Humid air example from Sphinx
from CoolProp.HumidAirProp import HAPropsSI
h = HAPropsSI('H','T',298.15,'P',101325,'R',0.5); print(h)
T = HAPropsSI('T','P',101325,'H',h,'R',1.0); print(T)
T = HAPropsSI('T','H',h,'R',1.0,'P',101325); print(T)

# Verification script
import CoolProp.CoolProp as CP
import numpy as np
import itertools
from multiprocessing import Pool

def generate_values(TR,P=101325):
    """ Starting with T,R as inputs, generate all other values """
    T,R = TR
    psi_w = CP.HAPropsSI('psi_w','T',T,'R',R,'P',P)
    other_output_keys = ['T_wb','T_dp','Hda','Sda','Vda','Omega']
    outputs = {'psi_w':psi_w,'T':T,'P':P,'R':R}
    for k in other_output_keys:
        outputs[k] = CP.HAPropsSI(k,'T',T,'R',R,'P',P)
    return outputs

def get_supported_input_pairs():
    """ Determine which input pairs are supported """
    good_ones = []
    inputs = generate_values((300, 0.5))
    for k1, k2 in itertools.product(inputs.keys(), inputs.keys()):
        if 'P' in [k1,k2] or k1==k2:
            continue
        args = ('psi_w', k1, inputs[k1], k2, inputs[k2], 'P', inputs['P'])
        try:
            psi_w_new = CP.HAPropsSI(*args)
            good_ones.append((k1,k2))
        except BaseException as BE:
            pass
            if 'currently at least one of' in str(BE) or 'cannot provide two inputs' in str(BE):
                pass
            else:
                print(BE)
                good_ones.append((k1,k2))
    return good_ones

def calculate(inputs):
    """ For a given input, try all possible input pairs """
    errors = []
    supported_pairs = get_supported_input_pairs()
    for k1, k2 in supported_pairs:
        psi_w_input = inputs['psi_w']
        args = 'psi_w',k1,inputs[k1],k2,inputs[k2],'P',inputs['P']
        try:
            psi_w_new = CP.HAPropsSI(*args)
        except BaseException as BE:
            errors.append((str(BE),args, inputs))
    return errors

if __name__ == '__main__':
    TR = itertools.product(np.linspace(240, 360, 31), np.linspace(0, 1, 31))
    with Pool(processes=2) as pool:
        input_values = pool.map(generate_values, TR)
        errors = pool.map(calculate, input_values)
        for err in itertools.chain.from_iterable(errors):
            print(err)
