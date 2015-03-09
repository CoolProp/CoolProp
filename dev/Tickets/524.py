import CoolProp
from CoolProp.CoolProp import PropsSI
from CoolProp.CoolProp import set_reference_state

print "CoolProp version",CoolProp.__version__
print "CoolProp revision",CoolProp.__gitrevision__

REF = 'R134a'

T0 = 273.15
RefState = 'IIR'
set_reference_state(REF,RefState)
print REF, RefState
print "HMass at",T0-273.15,"C",PropsSI('HMASS','T',T0,'Q',0,REF)/1000
print "SMass at",T0-273.15,"C",PropsSI('SMASS','T',T0,'Q',0,REF)/1000

T0 = 273.15-40
RefState = 'ASHRAE'
set_reference_state(REF,RefState)
print REF, RefState
print "HMass at",T0-273.15,"C",PropsSI('HMASS','T',T0,'Q',0,REF)/1000
print "SMass at",T0-273.15,"C",PropsSI('SMASS','T',T0,'Q',0,REF)/1000

P0 = 101325
RefState = 'NBP'
set_reference_state(REF,RefState)
print REF, RefState
print "HMass at",P0,"Pa",PropsSI('HMASS','P',P0,'Q',0,REF)/1000
print "SMass at",P0,"Pa",PropsSI('SMASS','P',P0,'Q',0,REF)/1000

T0 = 273.15
RefState = 'IIR'
set_reference_state(REF,RefState)
print REF, RefState
print "HMass at",T0-273.15,"C",PropsSI('HMASS','T',T0,'Q',0,REF)/1000
print "SMass at",T0-273.15,"C",PropsSI('SMASS','T',T0,'Q',0,REF)/1000

T0 = 273.15
RefState = 'DEF'
set_reference_state(REF,RefState)
print REF, RefState
print "HMass at",T0-273.15,"C",PropsSI('HMASS','T',T0,'Q',0,REF)/1000
print "SMass at",T0-273.15,"C",PropsSI('SMASS','T',T0,'Q',0,REF)/1000
