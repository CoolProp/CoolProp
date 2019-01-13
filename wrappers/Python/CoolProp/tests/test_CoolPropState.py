from __future__ import division, print_function
import CoolProp.CoolProp as CP
from CoolProp.State import State


def first_derivative(S, func, iVal, Val, iConstant, Constant, epsilon=1e-3):

    S.update({iVal: Val, iConstant: Constant})
    val1 = func()

    S.update({iVal: Val + epsilon, iConstant: Constant})
    val2 = func()

    S.update({iVal: Val, iConstant: Constant})

    return (val2 - val1) / epsilon


def second_derivative(S, func, iVal, Val, iConstant, Constant, epsilon=2):

    S.update({iVal: Val - epsilon, iConstant: Constant})
    val1 = func()

    S.update({iVal: Val, iConstant: Constant})
    val2 = func()

    S.update({iVal: Val + epsilon, iConstant: Constant})
    val3 = func()

    S.update({iVal: Val, iConstant: Constant})

    print(val1, val2, val3, S.T, S.p, S.rho, (val1 - 2 * val2 + val3))

    return (val1 - 2 * val2 + val3) / (epsilon * epsilon)


def teest_1phase_first_derivatives():

    for US in [CoolProp.UNIT_SYSTEM_SI, CoolProp.UNIT_SYSTEM_KSI]:
        CP.set_standard_unit_system(US)

        S = State('R134a', dict(T=300, D=1))

        l = [(S.get_rho, 'T', S.T, 'P', S.p, S.PFC.drhodT_constp),
             (S.get_rho, 'P', S.p, 'T', S.T, S.PFC.drhodp_constT),
             (S.get_p, 'D', S.rho, 'T', S.T, S.PFC.dpdrho_constT),
             # (S.get_p,'D',S.rho,'H',S.h,S.PFC.dpdrho_consth), #(these inputs not supported)
             (S.get_p, 'T', S.T, 'D', S.rho, S.PFC.dpdT_constrho),
             # (S.get_p,'T',S.T,'H',S.h,S.PFC.dpdT_consth),     #(these inputs not supported)
             (S.get_h, 'D', S.rho, 'T', S.T, S.PFC.dhdrho_constT),
             (S.get_h, 'D', S.rho, 'P', S.p, S.PFC.dhdrho_constp),
             (S.get_h, 'T', S.T, 'D', S.rho, S.PFC.dhdT_constrho),
             (S.get_h, 'T', S.T, 'P', S.p, S.PFC.dhdT_constp),
             (S.get_h, 'P', S.p, 'T', S.T, S.PFC.dhdp_constT),
             (S.get_s, 'D', S.rho, 'T', S.T, S.PFC.dsdrho_constT),
             (S.get_s, 'T', S.T, 'D', S.rho, S.PFC.dsdT_constrho),
             (S.get_s, 'D', S.rho, 'P', S.p, S.PFC.dsdrho_constp),
             (S.get_s, 'T', S.T, 'P', S.p, S.PFC.dsdT_constp),
             (S.get_s, 'P', S.p, 'T', S.T, S.PFC.dsdp_constT),

            ]
        for args in l:
            yield (check_1phase_first_derivatives,) + (S,) + args


def check_1phase_first_derivatives(S, func, iVal, Val, iConstant, Constant, deriv_func):

    Deriv_val = first_derivative(S, func, iVal, Val, iConstant, Constant)
    EOS_val = deriv_func()
    if abs(EOS_val / Deriv_val - 1) > 1e-2:
        raise ValueError('Finite Diff: ' + str(Deriv_val) + ' EOS: ' + str(EOS_val))


def teest_sat_first_derivatives():

    for US in [CoolProp.UNIT_SYSTEM_SI, CoolProp.UNIT_SYSTEM_KSI]:
        CP.set_standard_unit_system(US)

        S = State('R134a', dict(T=300, Q=1))

        l = [(S.get_T, 'P', S.p, 'Q', 0, S.PFC.dTdp_along_sat),
             (S.get_rho, 'P', S.p, 'Q', 0, S.PFC.drhodp_along_sat_liquid),
             (S.get_rho, 'P', S.p, 'Q', 1, S.PFC.drhodp_along_sat_vapor),
             (S.get_rho, 'T', S.T, 'Q', 0, S.PFC.drhodT_along_sat_liquid),
             (S.get_rho, 'T', S.T, 'Q', 1, S.PFC.drhodT_along_sat_vapor),
             (S.get_h, 'P', S.p, 'Q', 0, S.PFC.dhdp_along_sat_liquid),
             (S.get_h, 'P', S.p, 'Q', 1, S.PFC.dhdp_along_sat_vapor),
             (S.get_s, 'P', S.p, 'Q', 0, S.PFC.dsdp_along_sat_liquid),
             (S.get_s, 'P', S.p, 'Q', 1, S.PFC.dsdp_along_sat_vapor),
            ]
        for args in l:
            yield (check_sat_first_derivatives,) + (S,) + args


def check_sat_first_derivatives(S, func, iVal, Val, iConstant, Constant, deriv_func):

    Deriv_val = first_derivative(S, func, iVal, Val, iConstant, Constant)
    EOS_val = deriv_func()
    if abs(EOS_val / Deriv_val - 1) > 1e-2:
        raise ValueError('Finite Diff: ' + str(Deriv_val) + ' EOS: ' + str(EOS_val))


def teest_sat_second_derivatives():
    for US in [CoolProp.UNIT_SYSTEM_SI, CoolProp.UNIT_SYSTEM_KSI]:
        CP.set_standard_unit_system(US)

        S = State('R134a', dict(T=290, Q=1))

        l = [(S.get_T, 'P', S.p, 'Q', 0, S.PFC.d2Tdp2_along_sat),
             (S.get_rho, 'P', S.p, 'Q', 0, S.PFC.d2rhodp2_along_sat_liquid),
             (S.get_rho, 'P', S.p, 'Q', 1, S.PFC.d2rhodp2_along_sat_vapor),
             (S.get_h, 'P', S.p, 'Q', 0, S.PFC.d2hdp2_along_sat_liquid),
             (S.get_h, 'P', S.p, 'Q', 1, S.PFC.d2hdp2_along_sat_vapor),
             (S.get_s, 'P', S.p, 'Q', 0, S.PFC.d2sdp2_along_sat_liquid),
             (S.get_s, 'P', S.p, 'Q', 1, S.PFC.d2sdp2_along_sat_vapor),
            ]
        for args in l:
            yield (check_sat_second_derivatives,) + (S,) + args


def check_sat_second_derivatives(S, func, iVal, Val, iConstant, Constant, deriv_func):

    Deriv_val = second_derivative(S, func, iVal, Val, iConstant, Constant)
    EOS_val = deriv_func()
    if abs(EOS_val / Deriv_val - 1) > 1e-2:
        raise ValueError('Finite Diff: ' + str(Deriv_val) + ' EOS: ' + str(EOS_val))


if __name__ == '__main__':
    import nose
    nose.runmodule()
