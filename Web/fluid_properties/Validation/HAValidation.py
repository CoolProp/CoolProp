from __future__ import print_function

from CoolProp.HumidAirProp import HAPropsSI
import numpy as np

print(' Replicating the tables from ASHRAE RP-1485')
print('  ')
print('A.6.1 Psychrometric Properties of Moist Air at 0C and Below')
print('Saturated air at 101.325 kPa')
s5 = ' ' * 5
print('====================================================')
print("{T:8s}{W:10s}{v:10s}{h:10s}{s:10s}".format(W=s5 + ' Ws', v=s5 + '  v', h=s5 + 'h', s=s5 + ' s', T='   T'))
print("{T:8s}{W:10s}{v:10s}{h:10s}{s:10s}".format(W='  kgw/kg_da', v='   m3/kgda', h='  kJ/kgda', s=' kJ/kgda/K', T='   C'))
print('----------------------------------------------------')
for T in np.linspace(-60, 0, 13) + 273.15:
    h = HAPropsSI('H', 'T', T, 'R', 1.0, 'P', 101325) / 1000
    Twb = HAPropsSI('Twb', 'T', T, 'R', 1.0, 'P', 101325) - 273.15
    W = HAPropsSI('W', 'T', T, 'R', 1.0, 'P', 101325)
    v = HAPropsSI('V', 'T', T, 'R', 1.0, 'P', 101325)
    s = HAPropsSI('S', 'T', T, 'R', 1.0, 'P', 101325) / 1000
    print("{T:8.0f}{W:10.7f}{v:10.4f}{h:10.3f}{s:10.4f}".format(W=W, T=T - 273.15, v=v, h=h, s=s))
print('====================================================')
print(' ')
print('A.6.2 Psychrometric Properties of Moist Air at 0C and Above')
print('Saturated air at 101.325 kPa')
s5 = ' ' * 5
print('====================================================')
print("{T:8s}{W:10s}{v:10s}{h:10s}{s:10s}".format(W=s5 + ' Ws', v=s5 + '  v', h=s5 + 'h', s=s5 + ' s', T='   T'))
print("{T:8s}{W:10s}{v:10s}{h:10s}{s:10s}".format(W='  kgw/kg_da', v='   m3/kgda', h='  kJ/kgda', s=' kJ/kgda/K', T='   C'))
print('----------------------------------------------------')
for T in np.linspace(0, 90, 19) + 273.15:
    h = HAPropsSI('H', 'T', T, 'R', 1.0, 'P', 101325) / 1000
    Twb = HAPropsSI('Twb', 'T', T, 'R', 1.0, 'P', 101325) - 273.15
    W = HAPropsSI('W', 'T', T, 'R', 1.0, 'P', 101325)
    v = HAPropsSI('V', 'T', T, 'R', 1.0, 'P', 101325)
    s = HAPropsSI('S', 'T', T, 'R', 1.0, 'P', 101325) / 1000
    print("{T:8.0f}{W:10.7f}{v:10.3f}{h:10.2f}{s:10.4f}".format(W=W, T=T - 273.15, v=v, h=h, s=s))
print('====================================================')
print(' ')


def HotAir(num):
    from CoolProp.HumidAirProp import HAPropsSI
    if num == '8':
        Temp = str(200)
        T = 200 + 273.15
    elif num == '9':
        Temp = str(320)
        T = 320 + 273.15
    print('A.' + num + '.1 Psychrometric Properties of Moist Air at 101.325 kPa ')
    print('Dry Bulb temperature of ' + Temp + 'C')
    s5 = ' ' * 5
    print('================================================================')
    print("{W:10s}{Twb:10s}{v:10s}{h:10s}{s:10s}{R:10s}".format(W=s5 + ' W', Twb=s5 + 'Twb', v=s5 + '  v', h=s5 + 'h', s=s5 + ' s', R=s5 + 'RH'))
    print("{W:10s}{Twb:10s}{v:10s}{h:10s}{s:10s}{R:10s}".format(W='  kgw/kg_da', Twb='      C', v='   m3/kgda', h='  kJ/kgda', s=' kJ/kgda/K', R='    %'))
    print("----------------------------------------------------------------")
    for W in [0.0, 0.05, 0.1, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.0]:
        h = HAPropsSI('H', 'T', T, 'W', W, 'P', 101325) / 1000
        Twb = HAPropsSI('Twb', 'T', T, 'W', W, 'P', 101325) - 273.15
        R = HAPropsSI('R', 'T', T, 'W', W, 'P', 101325) * 100
        v = HAPropsSI('V', 'T', T, 'W', W, 'P', 101325)
        s = HAPropsSI('S', 'T', T, 'W', W, 'P', 101325) / 1000
        print("{W:10.2f}{Twb:10.2f}{v:10.3f}{h:10.2f}{s:10.4f}{R:10.4f}".format(W=W, Twb=Twb, v=v, h=h, s=s, R=R))
    print('================================================================')
    print(' ')
    print('A.' + num + '.2 Psychrometric Properties of Moist Air at 1000 kPa ')
    print('Dry Bulb temperature of ' + Temp + 'C')
    print('================================================================')
    print("{W:10s}{Twb:10s}{v:10s}{h:10s}{s:10s}{R:10s}".format(W=s5 + ' W', Twb=s5 + 'Twb', v=s5 + '  v', h=s5 + 'h', s=s5 + ' s', R=s5 + 'RH'))
    print("{W:10s}{Twb:10s}{v:10s}{h:10s}{s:10s}{R:10s}".format(W='  kgw/kg_da', Twb='      C', v='   m3/kgda', h='  kJ/kgda', s=' kJ/kgda/K', R='    %'))
    print('----------------------------------------------------------------')
    for W in [0.0, 0.05, 0.1, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.0]:
        h = HAPropsSI('H', 'T', T, 'W', W, 'P', 1000e3) / 1000
        Twb = HAPropsSI('Twb', 'T', T, 'W', W, 'P', 1000e3) - 273.15
        R = HAPropsSI('R', 'T', T, 'W', W, 'P', 1000e3) * 100
        v = HAPropsSI('V', 'T', T, 'W', W, 'P', 1000e3)
        s = HAPropsSI('S', 'T', T, 'W', W, 'P', 1000e3) / 1000
        print("{W:10.2f}{Twb:10.2f}{v:10.3f}{h:10.2f}{s:10.4f}{R:10.4f}".format(W=W, Twb=Twb, v=v, h=h, s=s, R=R))
    print('================================================================')
    print(' ')
    s5 = ' ' * 5
    print('A.' + num + '.3 Psychrometric Properties of Moist Air at 2000 kPa ')
    print('Dry Bulb temperature of ' + Temp + 'C')
    print('================================================================')
    print("{W:10s}{Twb:10s}{v:10s}{h:10s}{s:10s}{R:10s}".format(W=s5 + ' W', Twb=s5 + 'Twb', v=s5 + '  v', h=s5 + 'h', s=s5 + ' s', R=s5 + 'RH'))
    print("{W:10s}{Twb:10s}{v:10s}{h:10s}{s:10s}{R:10s}".format(W='  kgw/kg_da', Twb='      C', v='   m3/kgda', h='  kJ/kgda', s=' kJ/kgda/K', R='    %'))
    print('----------------------------------------------------------------')
    for W in [0.0, 0.05, 0.1, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.0]:
        h = HAPropsSI('H', 'T', T, 'W', W, 'P', 2000e3) / 1000
        Twb = HAPropsSI('Twb', 'T', T, 'W', W, 'P', 2000e3) - 273.15
        R = HAPropsSI('R', 'T', T, 'W', W, 'P', 2000e3) * 100
        v = HAPropsSI('V', 'T', T, 'W', W, 'P', 2000e3)
        s = HAPropsSI('S', 'T', T, 'W', W, 'P', 2000e3) / 1000
        print("{W:10.2f}{Twb:10.2f}{v:10.3f}{h:10.2f}{s:10.4f}{R:10.4f}".format(W=W, Twb=Twb, v=v, h=h, s=s, R=R))
    print('================================================================')
    print(' ')
    s5 = ' ' * 5
    print('A.' + num + '.4 Psychrometric Properties of Moist Air at 5000 kPa ')
    print('Dry Bulb temperature of ' + Temp + 'C')
    print('================================================================')
    print("{W:10s}{Twb:10s}{v:10s}{h:10s}{s:10s}{R:10s}".format(W=s5 + ' W', Twb=s5 + 'Twb', v=s5 + '  v', h=s5 + 'h', s=s5 + ' s', R=s5 + 'RH'))
    print("{W:10s}{Twb:10s}{v:10s}{h:10s}{s:10s}{R:10s}".format(W='  kgw/kg_da', Twb='      C', v='   m3/kgda', h='  kJ/kgda', s=' kJ/kgda/K', R='    %'))
    print('----------------------------------------------------------------')
    if Temp == '200':
        Wrange = [0.0, 0.05, 0.1, 0.15, 0.20, 0.25, 0.30]
    else:
        Wrange = [0.0, 0.05, 0.1, 0.15, 0.20, 0.25, 0.30, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    for W in Wrange:
        h = HAPropsSI('H', 'T', T, 'W', W, 'P', 5000e3) / 1000
        Twb = HAPropsSI('Twb', 'T', T, 'W', W, 'P', 5000e3) - 273.15
        R = HAPropsSI('R', 'T', T, 'W', W, 'P', 5000e3) * 100
        v = HAPropsSI('V', 'T', T, 'W', W, 'P', 5000e3)
        s = HAPropsSI('S', 'T', T, 'W', W, 'P', 5000e3) / 1000
        print("{W:10.2f}{Twb:10.2f}{v:10.3f}{h:10.2f}{s:10.4f}{R:10.4f}".format(W=W, Twb=Twb, v=v, h=h, s=s, R=R))
    print('================================================================')
    print(' ')
    s5 = ' ' * 5
    print('A.' + num + '.5 Psychrometric Properties of Moist Air at 10,000 kPa ')
    print('Dry Bulb temperature of ' + Temp + 'C')
    print('================================================================')
    print("{W:10s}{Twb:10s}{v:10s}{h:10s}{s:10s}{R:10s}".format(W=s5 + ' W', Twb=s5 + 'Twb', v=s5 + '  v', h=s5 + 'h', s=s5 + ' s', R=s5 + 'RH'))
    print("{W:10s}{Twb:10s}{v:10s}{h:10s}{s:10s}{R:10s}".format(W='  kgw/kg_da', Twb='      C', v='   m3/kgda', h='  kJ/kgda', s=' kJ/kgda/K', R='    %'))
    print('----------------------------------------------------------------')

    if Temp == '200':
        Wrange = [0.0, 0.05, 0.1]
    else:
        Wrange = [0.0, 0.05, 0.1, 0.15, 0.20, 0.25, 0.30, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    for W in Wrange:
        h = HAPropsSI('H', 'T', T, 'W', W, 'P', 10000e3) / 1000
        Twb = HAPropsSI('Twb', 'T', T, 'W', W, 'P', 10000e3) - 273.15
        R = HAPropsSI('R', 'T', T, 'W', W, 'P', 10000e3) * 100
        v = HAPropsSI('V', 'T', T, 'W', W, 'P', 10000e3)
        s = HAPropsSI('S', 'T', T, 'W', W, 'P', 10000e3) / 1000
        print("{W:10.2f}{Twb:10.2f}{v:10.3f}{h:10.2f}{s:10.4f}{R:10.4f}".format(W=W, Twb=Twb, v=v, h=h, s=s, R=R))
    print('================================================================')


HotAir('8')
print(' ')
HotAir('9')
##############################
#### Virial Coefficients #####
##############################


def Virials(variables):
    from CoolProp.HumidAirProp import HAProps_Aux
    import numpy as np

    varString = "%-10s" % ('T')
    units = "%-10s" % ('C')
    # Build the header
    for var in variables:
        varString += "%-20s" % (var)
        units += "%-20s" % (HAProps_Aux(var, 300, 100, 0.0)[1])
    print(varString)
    print(units)

    # Build the table
    for T in np.linspace(-60, 200, 27) + 273.15:
        values = "%-10.1f" % (T - 273.15)
        for var in variables:
            values += "%-20.10e" % (HAProps_Aux(var, T, 100, 0.0)[0])
        print(values)


print("")
print("Pure fluid Virial Coefficients")
print("------------------------------")
Virials(['Baa', 'Caaa', 'Bww', 'Cwww'])
Virials(['Baw', 'Caaw', 'Caww'])

print("")
print("Pure fluid Virial Coefficients Derivatives")
print("------------------------------------------")
Virials(['dBaa', 'dCaaa', 'dBww', 'dCwww'])
Virials(['dBaw', 'dCaaw', 'dCaww'])


##############################
####### Water Saturation #####
##############################

print("")
print("Water saturation pressure p_ws [kPa]")
from CoolProp.HumidAirProp import HAProps_Aux
import numpy as np
Tv = np.linspace(-60, 300, 13) + 273.15
print("%-10s %-20s" % ('T', 'p_ws'))
print("%-10s %-20s" % ('C', HAProps_Aux('p_ws', Tv[-1], 100, 0.0)[1]))
# Build the table
for T in Tv:
    values = "%-10.2f" % (T - 273.15)
    values += "%-20.10e" % (HAProps_Aux('p_ws', T, 100, 0.0)[0])
    print(values)

##############################
####### Henry Constant #######
##############################

print("")
print("Henry Constant (zero for T < 273.15 K)")
from CoolProp.HumidAirProp import HAProps_Aux
import numpy as np
Tv = np.linspace(0, 300, 11) + 273.16
print("%-10s %-20s" % ('T', 'beta_H'))
print("%-10s %-20s" % ('C', HAProps_Aux('beta_H', Tv[-1], 100, 0.0)[1]))
# Build the table
for T in Tv:
    values = "%-10.2f" % (T - 273.15)
    values += "%-20.10e" % (HAProps_Aux('beta_H', T, 100, 0.0)[0])
    print(values)

##########################################
####### Isothermal Compressibility #######
##########################################

print("")
print("Isothermal Compressibility of water (kT) [1/Pa]")
from CoolProp.HumidAirProp import HAProps_Aux
import numpy as np
Tv = np.linspace(-60, 300, 13) + 273.15
Pv = [101325, 200000, 500000, 1000000]
variables = "%-10s" % ('T')
for p in Pv:
    variables += "%-20s" % ("p = %-0.3f Pa " % (p))
print(variables)
# Build the actual table
for T in Tv:
    values = "%-10.2f" % (T - 273.15)
    for p in Pv:
        values += "%-20.10e" % (HAProps_Aux('kT', T, p, 0.0)[0])
    print(values)

##########################################
####### Saturated Molar Volume Water #####
##########################################

print("")
print("Molar volume of saturated liquid water or ice (vbar_ws) [m^3/mol_H2O]")
from CoolProp.HumidAirProp import HAProps_Aux
import numpy as np
Tv = np.linspace(-60, 300, 13) + 273.15
Pv = [101325, 200000, 500000, 1000000]
variables = "%-10s" % ('T')
for p in Pv:
    variables += "%-20s" % ("p = %-0.3f Pa " % (p))
print(variables)
# Build the actual table
for T in Tv:
    values = "%-10.2f" % (T - 273.15)
    for p in Pv:
        values += "%-20.10e" % (HAProps_Aux('vbar_ws', T, p, 0.0)[0])
    print(values)

##########################################
########### Enhancement Factor ###########
##########################################

print("")
print("Enhancement factor (f) [no units]")
from CoolProp.HumidAirProp import HAProps_Aux
import numpy as np
Tv = np.array([-60, -40, -20, 0, 40, 80, 120, 160, 200, 250, 300, 350]) + 273.15
Pv = [101325, 200000, 500000, 1000000, 10000000]
variables = "%-10s" % (u'T')
for p in Pv:
    variables += "%-20s" % ("p = %-0.3f Pa " % (p))
print(variables)
# Build the actual table
for T in Tv:
    values = "%-10.2f" % (T - 273.15)
    for p in Pv:
        values += "%-20.10e" % (HAProps_Aux('f', T, p, 0.0)[0])
    print(values)
