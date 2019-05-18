import CoolProp
from CoolProp.CoolProp import PropsSI
from CoolProp.CoolProp import set_reference_state

print("CoolProp version %s" % CoolProp.__version__)
print("CoolProp revision %s" % CoolProp.__gitrevision__)

REF = 'R134a'

T0 = 273.15
RefState = 'IIR'
set_reference_state(REF, RefState)
print(REF, RefState)
print("HMass at %s C %s" % (T0 - 273.15, "C", PropsSI('HMASS', 'T', T0, 'Q', 0, REF) / 1000))
print("SMass at %s C %s" % (T0 - 273.15, "C", PropsSI('SMASS', 'T', T0, 'Q', 0, REF) / 1000))

T0 = 273.15 - 40
RefState = 'ASHRAE'
set_reference_state(REF, RefState)
print(REF, RefState)
print("HMass at %s C %s" % (T0 - 273.15, PropsSI('HMASS', 'T', T0, 'Q', 0, REF) / 1000))
print("SMass at %s C %s" % (T0 - 273.15, PropsSI('SMASS', 'T', T0, 'Q', 0, REF) / 1000))

P0 = 101325
RefState = 'NBP'
set_reference_state(REF, RefState)
print(REF, RefState)
print("HMass at %s Pa %s" % (P0, PropsSI('HMASS', 'P', P0, 'Q', 0, REF) / 1000))
print("SMass at %s Pa %s" % (P0, PropsSI('SMASS', 'P', P0, 'Q', 0, REF) / 1000))

T0 = 273.15
RefState = 'IIR'
set_reference_state(REF, RefState)
print(REF, RefState)
print("HMass at %s C %s" % (T0 - 273.15, PropsSI('HMASS', 'T', T0, 'Q', 0, REF) / 1000))
print("SMass at %s C %s" % (T0 - 273.15, PropsSI('SMASS', 'T', T0, 'Q', 0, REF) / 1000))

T0 = 273.15
RefState = 'DEF'
set_reference_state(REF, RefState)
print(REF, RefState)
print("HMass at %s C %s" % (T0 - 273.15, "C", PropsSI('HMASS', 'T', T0, 'Q', 0, REF) / 1000))
print("SMass at %s C %s" % (T0 - 273.15, "C", PropsSI('SMASS', 'T', T0, 'Q', 0, REF) / 1000))
