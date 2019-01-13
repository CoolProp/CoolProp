'''
Created on 20 Sep 2013

@author: jowr
'''

# New example with R407F mixture
from pyrp.refpropClasses import RefpropSI
import CoolProp.CoolProp as cp

p = 30000
T = 273.15

ref = False

if ref:
    xkg = [0.473194694453358, 0.205109095413331, 0.321696210133311]
    names = "R32|R125|R134a"
    RP = RefpropSI()
    RP.SETUPFLEX(xkg=xkg, FluidNames=names)
    T_A, p_A, D_A, Dl_A, Dv_A, q_A, e_A, h_A, s_A, cv_A, cp_A, w_A = RP.PQFLSH(p, 0)
    T_B, p_B, D_B, Dl_B, Dv_B, q_B, e_B, h_B, s_B, cv_B, cp_B, w_B = RP.PQFLSH(p, 1)
    T_C, p_C, D_C, Dl_C, Dv_C, q_C, e_C, h_C, s_C, cv_C, cp_C, w_C = RP.TQFLSH(T, 0)
    hlb = h_A / 1000.
    hrb = h_B / 1000.
    h200 = h_C / 1000.
    print("Refprop: %s %s %s" % (hlb, hrb, h200))
else:
    R407F = 'REFPROP-MIX:R32[0.473194694453358]&R125[0.205109095413331]&R134a[0.321696210133311]'
    # R407F='REFPROP-MIX:R32[0.651669604033581]&R125[0.122438378639971]&R134a[0.225892017326446]'
    hlb = cp.Props('H', 'P', 30, 'Q', 0, R407F)  # 30 kPa saturated liquid
    hrb = cp.Props('H', 'P', 30, 'Q', 1, R407F)  # 30 kPa saturated vapour
    h200 = cp.Props('H', 'T', 273.15, 'Q', 0, R407F)  # saturated liquid at 0C IIR
    print("CoolProp: %s %s %s" % (hlb, hrb, h200))
