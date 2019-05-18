from __future__ import division
from CoolProp.CoolProp import Props, DerivTerms


def finite_diff(Output, Input1, Val1, Input2, Val2, Fluid, index, deriv, order, type='centered', dx=1e-3):
    def f(index, dx):
        if index == 1:
            return Props(Output, Input1, Val1 + dx, Input2, Val2, Fluid)
        else:
            return Props(Output, Input1, Val1, Input2, Val2 + dx, Fluid)

    if deriv == 1:
        if order == 1:
            return (f(index, dx) - f(index, 0)) / (dx)
        elif order == 2:
            return (f(index, dx) - f(index, -dx)) / (2 * dx)
        elif order == 4:
            return (f(index, -2 * dx) - 8 * f(index, -dx) + 8 * f(index, dx) - f(index, 2 * dx)) / (12 * dx)
    elif deriv == 2:
        if order == 2:
            return (f(index, -dx) - 2 * f(index, 0) + f(index, dx)) / (dx**2)
        elif order == 4:
            return (-f(index, -2 * dx) + 16 * f(index, -dx) - 30 * f(index, 0) + 16 * f(index, dx) - f(index, 2 * dx)) / (12 * dx**2)
        else:
            raise ValueError
    else:
        raise ValueError

## fluid = 'Water'
# T = 647.74374374374372#647#Props(fluid,'Tcrit')+1
# rho = 322.32199999#358#Props(fluid,'rhocrit')+1
## H = Props('H','T',T,'D',rho,fluid)
## P = Props('P','T',T,'D',rho,fluid)
# print T,rho,H,P


fluid = 'R125'
T = 300
rho = 1.5
H = Props('H', 'T', T, 'D', rho, fluid)
P = Props('P', 'T', T, 'D', rho, fluid)
print("%s %s %s %s" % (T, rho, H, P))

dpdT__rho = DerivTerms('dpdT|rho', T, rho, fluid)
dpdrho__T = DerivTerms('dpdrho|T', T, rho, fluid)
dhdT__rho = DerivTerms('dhdT|rho', T, rho, fluid)
dhdrho__T = DerivTerms('dhdrho|T', T, rho, fluid)
print('*******************************************')
print('CHECKING DERIVATIVES FROM EOS')
print('*******************************************')
dpdT__rho_num = finite_diff('P', 'T', T, 'D', rho, fluid, 1, 1, 4)
print('dpdT|rho')
print(dpdT__rho)
print(dpdT__rho_num)

dpdrho__T_num = finite_diff('P', 'T', T, 'D', rho, fluid, 2, 1, 4)
print('dpdrho|T')
print(dpdrho__T)
print(dpdrho__T_num)

dhdT__rho_num = finite_diff('H', 'T', T, 'D', rho, fluid, 1, 1, 4)
print('dhdT|rho')
print(dhdT__rho)
print(dhdT__rho_num)

dhdrho__T_num = finite_diff('H', 'T', T, 'D', rho, fluid, 2, 1, 4)
print('dhdrho|T')
print(dhdrho__T)
print(dhdrho__T_num)
print("")
print('*******************************************')
print('CHECKING HELMHOLTZ DERIVATIVES')
print('*******************************************')
dx = 1e-10
Tc = Props(fluid, 'Tcrit')
rhoc = Props(fluid, 'rhocrit')
delta = rho / rhoc
tau = Tc / T


def f(dx):
    rho = rhoc * (delta + dx)
    return DerivTerms('dphir_dDelta', T, rho, fluid)


print('d2phir_dDelta2')
print((f(-2 * dx) - 8 * f(-dx) + 8 * f(dx) - f(2 * dx)) / (12 * dx))
print(DerivTerms('d2phir_dDelta2', T, rho, fluid))


def f(dx):
    rho = rhoc * (delta + dx)
    return DerivTerms('d2phir_dDelta2', T, rho, fluid)


print('d3phir_dDelta3')
print((f(-2 * dx) - 8 * f(-dx) + 8 * f(dx) - f(2 * dx)) / (12 * dx))
print(DerivTerms('d3phir_dDelta3', T, rho, fluid))


def f(dx):
    rho = rhoc * (delta + dx)
    return DerivTerms('dphir_dTau', T, rho, fluid)


print('d2phir_dDelta_dTau')
print((f(-2 * dx) - 8 * f(-dx) + 8 * f(dx) - f(2 * dx)) / (12 * dx))
print(DerivTerms('d2phir_dDelta_dTau', T, rho, fluid))


def f(dx):
    rho = rhoc * (delta + dx)
    return DerivTerms('d2phir_dTau2', T, rho, fluid)


print('d3phir_dDelta_dTau2')
print((f(-2 * dx) - 8 * f(-dx) + 8 * f(dx) - f(2 * dx)) / (12 * dx))
print(DerivTerms('d3phir_dDelta_dTau2', T, rho, fluid))


def f(dx):
    T = Tc / (tau + dx)
    return DerivTerms('d2phir_dDelta2', T, rho, fluid)


print('d3phir_dDelta2_dTau')
print((f(-2 * dx) - 8 * f(-dx) + 8 * f(dx) - f(2 * dx)) / (12 * dx))
print(DerivTerms('d3phir_dDelta2_dTau', T, rho, fluid))


def f(dx):
    T = Tc / (tau + dx)
    return DerivTerms('dphir_dTau', T, rho, fluid)


print('d2phir_dTau2')
print((f(-2 * dx) - 8 * f(-dx) + 8 * f(dx) - f(2 * dx)) / (12 * dx))
print(DerivTerms('d2phir_dTau2', T, rho, fluid))


def f(dx):
    T = Tc / (tau + dx)
    return DerivTerms('d2phir_dTau2', T, rho, fluid)


print('d3phir_dTau3')
print((f(-2 * dx) - 8 * f(-dx) + 8 * f(dx) - f(2 * dx)) / (12 * dx))
print(DerivTerms('d3phir_dTau3', T, rho, fluid))


def f(dx):
    T = Tc / (tau + dx)
    return DerivTerms('dphi0_dTau', T, rho, fluid)


print('d2phi0_dTau2')
print((f(-2 * dx) - 8 * f(-dx) + 8 * f(dx) - f(2 * dx)) / (12 * dx))
print(DerivTerms('d2phi0_dTau2', T, rho, fluid))


def f(dx):
    T = Tc / (tau + dx)
    return DerivTerms('d2phi0_dTau2', T, rho, fluid)


print('d3phi0_dTau3')
print((f(-2 * dx) - 8 * f(-dx) + 8 * f(dx) - f(2 * dx)) / (12 * dx))
print(DerivTerms('d3phi0_dTau3', T, rho, fluid))
print("")
print('*******************************************')
print('CHECKING FIRST DERIVATIVES OF PROPERTIES')
print('*******************************************')

A = dpdT__rho * dhdrho__T - dpdrho__T * dhdT__rho

dT_dh_num = finite_diff('T', 'H', H, 'P', P, fluid, 1, 1, 4)
dT_dh = 1 / Props('C', 'T', T, 'D', rho, fluid)
print('dTdh|p')
print(dT_dh)
print(dT_dh_num)

dT_dp_num = finite_diff('T', 'H', H, 'P', P, fluid, 2, 1, 4)
dT_dp = 1 / A * dhdrho__T
print('dTdp|h')
print(dT_dp)
print(dT_dp_num)

drho_dh_num = finite_diff('D', 'H', H, 'P', P, fluid, 1, 1, 4)
drho_dh = 1 / A * dpdT__rho
print('drhodh|p')
print(drho_dh_num)
print(drho_dh_num)

drho_dp_num = finite_diff('D', 'H', H, 'P', P, fluid, 2, 1, 4)
drho_dp = -1 / A * dhdT__rho
print('drhodp|h')
print(drho_dp)
print(drho_dp_num)

ds_dh_num = finite_diff('S', 'H', H, 'P', P, fluid, 1, 1, 4)
ds_dh = 1 / T
print('dsdh|p')
print(ds_dh)
print(ds_dh_num)

ds_dp_num = finite_diff('S', 'H', H, 'P', P, fluid, 2, 1, 4)
ds_dp = -1 / (T * rho)
print('dsdp|h')
print(ds_dp)
print(ds_dp_num)


print("")
print('*******************************************')
print('CHECKING SECOND DERIVATIVES OF PROPERTIES')
print('*******************************************')

## d2T_dh2_num = finite_diff('T','H',H,'P',P,fluid,1,2,4)
## d2T_dh2 = 1/Props('C','T',T,'D',rho,fluid)
# print 'dTdh|p'
# print dT_dh
# print dT_dh_num

d2s_dh2 = -(1 / T / T) * dT_dh
d2s_dh2_num = finite_diff('S', 'H', H, 'P', P, fluid, 1, 2, 4, dx=10)
print('d2sdh2|p')
print(d2s_dh2)
print(d2s_dh2_num)

d2s_dp2 = 1 / (T**2 * rho) * dT_dp + 1 / (T * rho**2) * drho_dp
d2s_dp2_num = finite_diff('S', 'H', H, 'P', P, fluid, 2, 2, 4, dx=10)
print('d2sdp2|h')
print(d2s_dp2)
print(d2s_dp2_num)

d2s_dhdp = -(1 / T / T) * dT_dp
d2s_dhdp_num = finite_diff('S', 'H', H, 'P', P, fluid, 1, 2, 4)
print('d2s_dhdp')
print(d2s_dhdp)
print(d2s_dhdp_num)

## d2s_dhdp = -(1/T/T)*dT_dp
# print d2s_dhdp
## d2s_dhdp = 1/(T**2*rho)*dT_dh+1/(T*rho**2)*drho_dh
# print d2s_dhdp
