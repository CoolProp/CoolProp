from CoolProp.HumidAirProp import HAProps


print("Validation against H.F. Nelson and H.J. Sauer,\"Formulation for High-Temperature Properties for Moist Air\", HVAC&R Research v.8 #3, 2002")
print("Note: More accurate formulation employed than in Nelson.  Just for sanity checking")
print("Yields a negative relative humidity for Tdb=5C,Twb=-3C, point omitted")
tdb = [5, 5, 5, 25, 25, 25, 25, 50, 50, 50, 50, 50, 50, 50]
twb = [5, 2, -1, 25, 20, 15, 10, 50, 40, 30, 25, 22, 20, 19]
print(" ")
print("Table 6: Adiabatic Saturation")
print("P=101325 Pa, Altitude = 0 m")
print("========================================================================")
print("{Tdb:10s}{Twb:10s}{Tdp:10s}{R:10s}{W:10s}{h:10s}{v:10s}".format(W='W', Twb='Twb', Tdp='Tdp', Tdb='Tdb', v='v', h='h', s='s', R='RH'))
print("{Tdb:10s}{Twb:10s}{Tdp:10s}{R:10s}{W:10s}{h:10s}{v:10s}".format(W='-', Twb='C', Tdp='C', Tdb='C', v='m^3/kg_da', h='kJ/kg_da', s='kJ/kg_da/K', R='%'))
print("------------------------------------------------------------------------")
for (tdb_, twb_) in zip(tdb, twb):
    h = HAProps('H', 'T', tdb_ + 273.13, 'Twb', twb_ + 273.15, 'P', 101.325)
    tdp = HAProps('Tdp', 'T', tdb_ + 273.13, 'Twb', twb_ + 273.15, 'P', 101.325) - 273.15
    W = HAProps('W', 'T', tdb_ + 273.13, 'Twb', twb_ + 273.15, 'P', 101.325)
    R = HAProps('R', 'T', tdb_ + 273.13, 'Twb', twb_ + 273.15, 'P', 101.325) * 100
    v = HAProps('V', 'T', tdb_ + 273.13, 'Twb', twb_ + 273.15, 'P', 101.325)
    s = 0
    print("{Tdb:10.2f}{Twb:10.2f}{Tdp:10.2f}{R:10.1f}{W:10.5f}{h:10.2f}{v:10.3f}".format(W=W, Twb=twb_, Tdp=tdp, Tdb=tdb_, v=v, h=h, s=s, R=R))
print("------------------------------------------------------------------------")
print(" ")
print("Table 7: Adiabatic Saturation")
print("P=84,556 Pa, Altitude = 1500 m")
print("========================================================================")
print("{Tdb:10s}{Twb:10s}{Tdp:10s}{R:10s}{W:10s}{h:10s}{v:10s}".format(W='W', Twb='Twb', Tdp='Tdp', Tdb='Tdb', v='v', h='h', s='s', R='RH'))
print("{Tdb:10s}{Twb:10s}{Tdp:10s}{R:10s}{W:10s}{h:10s}{v:10s}".format(W='-', Twb='C', Tdp='C', Tdb='C', v='m^3/kg_da', h='kJ/kg_da', s='kJ/kg_da/K', R='%'))
print("------------------------------------------------------------------------")
for (tdb_, twb_) in zip(tdb, twb):
    h = HAProps('H', 'T', tdb_ + 273.13, 'Twb', twb_ + 273.15, 'P', 84.556)
    tdp = HAProps('Tdp', 'T', tdb_ + 273.13, 'Twb', twb_ + 273.15, 'P', 84.556) - 273.15
    W = HAProps('W', 'T', tdb_ + 273.13, 'Twb', twb_ + 273.15, 'P', 84.556)
    R = HAProps('R', 'T', tdb_ + 273.13, 'Twb', twb_ + 273.15, 'P', 84.556) * 100
    v = HAProps('V', 'T', tdb_ + 273.13, 'Twb', twb_ + 273.15, 'P', 84.556)
    s = 0
    print("{Tdb:10.2f}{Twb:10.2f}{Tdp:10.2f}{R:10.1f}{W:10.5f}{h:10.2f}{v:10.3f}".format(W=W, Twb=twb_, Tdp=tdp, Tdb=tdb_, v=v, h=h, s=s, R=R))
print("------------------------------------------------------------------------")
