from CoolProp.HumidAirProp import HAPropsSI

# All humid-air calls are SI (the v8 interface): pressure in Pa, enthalpy in
# J/kg of dry air, temperatures in K.  (The non-SI ``HAProps`` was removed in
# CoolProp v8; bd CoolProp-r9sq.27.)  The ``h`` column below is printed in
# kJ/kg_da, so the SI enthalpy is divided by 1000 on the way out.

print("Validation against H.F. Nelson and H.J. Sauer,\"Formulation for High-Temperature Properties for Moist Air\", HVAC&R Research v.8 #3, 2002")
print("Note: More accurate formulation employed than in Nelson.  Just for sanity checking")
print("Yields a negative relative humidity for Tdb=5C,Twb=-3C, point omitted")
tdb = [5, 5, 5, 25, 25, 25, 25, 50, 50, 50, 50, 50, 50, 50]
twb = [5, 2, -1, 25, 20, 15, 10, 50, 40, 30, 25, 22, 20, 19]


def _table(P):
    print("{Tdb:10s}{Twb:10s}{Tdp:10s}{R:10s}{W:10s}{h:10s}{v:10s}".format(W='W', Twb='Twb', Tdp='Tdp', Tdb='Tdb', v='v', h='h', R='RH'))
    print("{Tdb:10s}{Twb:10s}{Tdp:10s}{R:10s}{W:10s}{h:10s}{v:10s}".format(W='-', Twb='C', Tdp='C', Tdb='C', v='m^3/kg_da', h='kJ/kg_da', R='%'))
    print("------------------------------------------------------------------------")
    for (tdb_, twb_) in zip(tdb, twb):
        try:
            h = HAPropsSI('H', 'T', tdb_ + 273.13, 'Twb', twb_ + 273.15, 'P', P) / 1000.0
            tdp = HAPropsSI('Tdp', 'T', tdb_ + 273.13, 'Twb', twb_ + 273.15, 'P', P) - 273.15
            W = HAPropsSI('W', 'T', tdb_ + 273.13, 'Twb', twb_ + 273.15, 'P', P)
            R = HAPropsSI('R', 'T', tdb_ + 273.13, 'Twb', twb_ + 273.15, 'P', P) * 100
            v = HAPropsSI('V', 'T', tdb_ + 273.13, 'Twb', twb_ + 273.15, 'P', P)
        except ValueError:
            # Some near-saturation points fall outside CoolProp's RH validity range
            # (the header notes one such omission); skip rather than abort the table.
            print("  (Tdb={Tdb:.0f}C, Twb={Twb:.0f}C omitted: outside validity range)".format(Tdb=tdb_, Twb=twb_))
            continue
        print("{Tdb:10.2f}{Twb:10.2f}{Tdp:10.2f}{R:10.1f}{W:10.5f}{h:10.2f}{v:10.3f}".format(W=W, Twb=twb_, Tdp=tdp, Tdb=tdb_, v=v, h=h, R=R))
    print("------------------------------------------------------------------------")


print(" ")
print("Table 6: Adiabatic Saturation")
print("P=101325 Pa, Altitude = 0 m")
print("========================================================================")
_table(101325.0)
print(" ")
print("Table 7: Adiabatic Saturation")
print("P=84,556 Pa, Altitude = 1500 m")
print("========================================================================")
_table(84556.0)
