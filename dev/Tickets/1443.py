from __future__ import print_function, division

import sys, threading

from timeit import default_timer as timer
from tqdm import tqdm

import numpy as np

import CoolProp
from CoolProp.CoolProp import AbstractState as AS
from CoolProp.CoolProp import PropsSI as PS
from CoolProp.CoolProp import get_fluid_param_string, get_config_bool, set_config_bool, generate_update_pair, get_parameter_index
from CoolProp.Plots.Common import process_fluid_state

CONF_CPP_DIR = "Call PropsSI_multi"
CONF_CPP_GUE = "Call PropsSI_multi, guesses enabled"
CONF_PYT_DIR = "Call AbstractState.update"
CONF_PYT_GUE = "Call AbstractState.update_with_guesses"


class ResClass(object):
    def __init__(self):
        self.time = None
        self.err_count = None
        self.size = None
        self.comment = None
        self.fluid = None

    def __str__(self):
        return "{0:40s} - {1:8.2f} us, {2} errors, {3} elements ({4})".format(self.fluid, self.time, self.err_count, self.size, self.comment)


def calc_isolines(fluid="HEOS::water", kind_a="P", range_a=[1e6, 10e6], kind_b="T", range_b=[50 + 273.15], kind_o=["Hmass"], CUR_CONF=CONF_PYT_DIR):

    state = process_fluid_state(fluid)

    kind_a_index = get_parameter_index(kind_a)
    kind_b_index = get_parameter_index(kind_b)
    kind_o_index = [get_parameter_index(kind_o_single) for kind_o_single in kind_o]

    (input_pair, one, two) = generate_update_pair(kind_a_index, 0, kind_b_index, 1)
    if one == 0:
        swap = False
    else:
        swap = True

    range_one = np.asanyarray(range_b)
    range_two = np.asanyarray(range_a)

    result = []

    # for i in range(10):
    #    pbar.update(10)
    # pbar.close()

    progress_bar = tqdm(total=range_two.size * range_one.size, desc="{0:40s}".format(fluid))
    progress_bar.clear()

    for cnt, one in np.ndenumerate(range_one):
        if not swap:
            kind_in1 = kind_b
            kind_in2 = kind_a
            vector_in1 = np.zeros_like(range_two) + one
            vector_in2 = range_two
        else:
            kind_in2 = kind_b
            kind_in1 = kind_a
            vector_in2 = np.zeros_like(range_two) + one
            vector_in1 = range_two

        # At this point, everything is prepared for the loop
        # vector_in1 contains the first set of inputs and
        # vector_in2 contains the second set of inputs

        #CONF_CPP_DIR = "Call PropsSI_multi"
        #CONF_CPP_GUE = "Call PropsSI_multi, guesses enabled"
        #CONF_PYT_DIR = "Call AbstractState.update"
        # CONF_PYT_GUE

        tmp_result = 0.0
        err_count = 0
        single_result = ResClass()
        guesses = CoolProp.CoolProp.PyGuessesStructure()

        if CUR_CONF == CONF_CPP_DIR:
            set_config_bool(CoolProp.USE_GUESSES_IN_PROPSSI, False)
            start = timer()
            try:
                tmp_result = PS(kind_o, kind_in1, vector_in1, kind_in2, vector_in2, fluid)
            except:
                err_count += 1
                pass
            end = timer()

        elif CUR_CONF == CONF_CPP_GUE:
            set_config_bool(CoolProp.USE_GUESSES_IN_PROPSSI, True)
            start = timer()
            try:
                tmp_result = PS(kind_o, kind_in1, vector_in1, kind_in2, vector_in2, fluid)
            except:
                err_count += 1
                pass
            end = timer()

        elif CUR_CONF == CONF_PYT_DIR:
            start = timer()
            for i, two in np.ndenumerate(vector_in1):
                try:
                    state.update(input_pair, two, one)
                    for kind_o_single in kind_o_index:
                        tmp_result = state.keyed_output(kind_o_single)
                except:
                    err_count += 1
                    pass
            end = timer()

        elif CUR_CONF == CONF_PYT_GUE:
            start = timer()
            for i, two in np.ndenumerate(vector_in1):
                try:
                    if i < 1: state.update(input_pair, two, one)
                    else: state.update_with_guesses(input_pair, two, one, guess)
                    guesses.rhomolar = state.rhomolar()
                    guesses.T = state.T()
                    for kind_o_single in kind_o_index:
                        tmp_result = state.keyed_output(kind_o_single)
                except:
                    err_count += 1
                    guesses.rhomolar = np.NaN
                    guesses.T = np.NaN
                    pass
            end = timer()

        single_result.time = (end - start) * 1e6 / (1.0 * range_two.size)
        single_result.err_count = err_count
        single_result.size = range_two.size
        single_result.fluid = fluid
        single_result.comment = CUR_CONF
        result.append(single_result)

        progress_bar.update(range_two.size)

    progress_bar.close()

    return result


# if __name__ == "__main__":
#    print("two.py is being run directly")
# else:
#    sys.exit(1)

# def __main__():
#    CoolProp.CoolProp.set_config_bool(configuration_keys key, bool value)


import numpy as np


def get_fluid_strings(mix=False):
    if mix:
        fluids_in = ["Propane[0.5]&Ethane[0.5]", "ISOBUTAN[0.8]&PROPANE[0.2]", "R32[0.697615]&R125[0.302385]"]
    else:
        fluids_in = ["Water", "R134a", "Air"]
    fluids = []
    for fld in fluids_in:
        for bac in ["HEOS", "REFPROP"]:
            if bac == "REFPROP" and not mix:
                fluids.append(bac + "::" + get_fluid_param_string(fld, "REFPROP_name"))
            else:
                fluids.append(bac + "::" + fld)
    return fluids


def get_state_objects(fluids=None):
    if fluids is None:
        fluids = get_fluid_strings()
    states = []
    for fld in fluids:
        states.append(process_fluid_state(fld))
    return states


steps = 250


def get_p_range(state=AS("HEOS", "Water")):
    p_max = 100.0e5
    p_min = 0.001e5
    p_range = np.logspace(np.log10(p_min), np.log10(p_max), steps, base=10)
    return p_range


def get_T_range(state=AS("HEOS", "Water")):
    T_max = 100.0 + 273.15
    T_min = -50.0 + 273.15
    T_range = np.linspace(T_min, T_max, steps)
    return T_range


def get_T_iso(state=AS("HEOS", "Water")):
    try:
        state.build_phase_envelope("dummy")
    except:
        pass
    try:
        T_c = state.T_critical()
    except:
        try:
            T_c, r_c = state.true_critical_point()
        except:
            T_c = (state.Tmin() + state.Tmax()) / 2.0
    T_min = state.Tmin()
    T_max = state.Tmax()
    return np.array([(T_c + T_min) / 2.0, T_c, (T_c + T_max) / 2.0])


def get_h_iso(state=AS("HEOS", "Water")):
    res = []
    try:
        T_c = state.T_critical()
        r_c = state.rhomolar_critical()
        state.update(CoolProp.DmolarT_INPUTS, r_c, T_c)
        res.append(state.hmass())
        state.update(CoolProp.QT_INPUTS, 0, 0.75 * T_c)
        res.append(state.hmass())
        state.update(CoolProp.QT_INPUTS, 1, 0.75 * T_c)
        res.append(state.hmass())
        return np.array(res)
    except:
        state.build_phase_envelope("dummy")
        PE = state.get_phase_envelope_data()
        T_c = PE.T[int(len(PE.T) / 2)]
        p_c = PE.p[int(len(PE.p) / 2)]
        state.update(CoolProp.PQ_INPUTS, p_c, 0.0)
        res.append(state.hmass())
        state.update(CoolProp.PQ_INPUTS, p_c * 0.75, 0.0)
        res.append(state.hmass())
        state.update(CoolProp.PQ_INPUTS, p_c * 0.75, 1.0)
        res.append(state.hmass())
        return np.array(res)


results = {}
mix = True
for fluid, state in zip(get_fluid_strings(mix), get_state_objects(get_fluid_strings(mix))):
    p_range = get_p_range(state)
    T_range = get_T_iso(state)
    h_range = get_h_iso(state)

    results[fluid] = {}

    for CNF in [CONF_CPP_DIR, CONF_CPP_GUE, CONF_PYT_DIR, CONF_PYT_GUE]:
        results[fluid][CNF] = {}
        results[fluid][CNF]["PT"] = calc_isolines(fluid=fluid, kind_a="P", range_a=p_range, kind_b="T", range_b=T_range, kind_o=["Hmass", "Dmass", "T", "P"], CUR_CONF=CNF)
        # for result in results:
        #    print(str(result) + "p,T-inputs")
        results[fluid][CNF]["HP"] = calc_isolines(fluid=fluid, kind_a="P", range_a=p_range, kind_b="Hmass", range_b=h_range, kind_o=["Hmass", "Dmass", "T", "P"], CUR_CONF=CNF)
        # for result in results:
        #    print(str(result) + "h,p-inputs")

for FLD in sorted(results.keys()):
    for CNF in sorted(results[FLD].keys()):
        for INP in sorted(results[FLD][CNF].keys()):
            cur_res = ResClass()
            cur_res.err_count = 0
            cur_res.fluid = FLD
            cur_res.comment = CNF
            cur_res.size = 0
            cur_res.time = 0
            for res in results[FLD][CNF][INP]:
                cur_res.err_count += res.err_count
                if cur_res.time == 0:
                    cur_res.time = res.time
                else:
                    cur_res.time = (cur_res.time * cur_res.size + res.time * res.size) / (cur_res.size + res.size)
                cur_res.size += res.size
            print(str(cur_res) + " " + INP + "-inputs")
