#!/usr/bin/python
import sys
import CoolProp.CoolProp as CP
# print "{0:14.8f}".format(CP.Props('V','D',13,'P',500,'n-Pentane'))
# print "{0:14.8f}".format(CP.Props('V','H',158,'P',1000,'TX22'))
#T = 300
T = float(sys.argv[1])
P = float(sys.argv[2])
print("Temperature: " + str(T))
print("Pressure:    " + str(P))
print("")
print("Viscosity: ")
print("{0:14.8f}".format(CP.Props('V', 'T', T, 'P', P, 'SecCoolSolution-20%')))
print("{0:14.8f}".format(CP.Props('V', 'T', T, 'P', P, 'MEG-20%')))
print("{0:14.8f}".format(CP.Props('V', 'T', T, 'P', P, 'EG-20%')))
print("{0:14.8f}".format(CP.Props('V', 'T', T, 'P', P, 'water')))
# print
# print "{0:14.8f}".format(CP.Props('L','D',13,'P',500,'n-Pentane'))
# print "{0:14.8f}".format(CP.Props('L','H',158,'P',P,'TX22'))
print("")
print("Conductivity: ")
print("{0:14.8f}".format(CP.Props('L', 'T', T, 'P', P, 'SecCoolSolution-20%')))
print("{0:14.8f}".format(CP.Props('L', 'T', T, 'P', P, 'MEG-20%')))
print("{0:14.8f}".format(CP.Props('L', 'T', T, 'P', P, 'EG-20%')))
print("{0:14.8f}".format(CP.Props('L', 'T', T, 'P', P, 'water')))
print("")
print("Density: ")
print("{0:14.8f}".format(CP.Props('D', 'T', T, 'P', P, 'SecCoolSolution-20%')))
print("{0:14.8f}".format(CP.Props('D', 'T', T, 'P', P, 'MEG-20%')))
print("{0:14.8f}".format(CP.Props('D', 'T', T, 'P', P, 'EG-20%')))
print("{0:14.8f}".format(CP.Props('D', 'T', T, 'P', P, 'water')))
print("")
print("Capacity: ")
print("{0:14.8f}".format(CP.Props('C', 'T', T, 'P', P, 'SecCoolSolution-20%')))
print("{0:14.8f}".format(CP.Props('C', 'T', T, 'P', P, 'MEG-20%')))
print("{0:14.8f}".format(CP.Props('C', 'T', T, 'P', P, 'EG-20%')))
print("{0:14.8f}".format(CP.Props('C', 'T', T, 'P', P, 'water')))
print("")
print("Enthalpy: ")
print("{0:14.8f}".format(CP.Props('H', 'T', T, 'P', P, 'SecCoolSolution-20%')))
print("{0:14.8f}".format(CP.Props('H', 'T', T, 'P', P, 'MEG-20%')))
print("{0:14.8f}".format(CP.Props('H', 'T', T, 'P', P, 'EG-20%')))
print("{0:14.8f}".format(CP.Props('H', 'T', T, 'P', P, 'water')))
print("")
print("Internal energy: ")
print("{0:14.8f}".format(CP.Props('U', 'T', T, 'P', P, 'SecCoolSolution-20%')))
print("{0:14.8f}".format(CP.Props('U', 'T', T, 'P', P, 'MEG-20%')))
print("-")  # "{0:14.8f}".format(CP.Props('U','T',T,'P',P,'EG-20%'))
print("{0:14.8f}".format(CP.Props('U', 'T', T, 'P', P, 'water')))
print("")
print("Entropy: ")
print("{0:14.8f}".format(CP.Props('S', 'T', T, 'P', P, 'SecCoolSolution-20%')))
print("{0:14.8f}".format(CP.Props('S', 'T', T, 'P', P, 'MEG-20%')))
print("{0:14.8f}".format(CP.Props('S', 'T', T, 'P', P, 'EG-20%')))
print("{0:14.8f}".format(CP.Props('S', 'T', T, 'P', P, 'water')))
print("")
print("Freezing point: ")
print("{0:14.8f}".format(CP.Props('Tfreeze', 'T', T, 'P', P, 'SecCoolSolution-20%')))
print("{0:14.8f}".format(CP.Props('Tfreeze', 'T', T, 'P', P, 'MEG-20%')))
print("{0:14.8f}".format(CP.Props('F', 'T', T, 'P', P, 'EG-20%')))
print("")
print("TX22: ")
print("{0:14.8f}".format(CP.Props('H', 'T', T + 50, 'P', P, 'TX22')))
print("{0:14.8f}".format(CP.Props('U', 'T', T + 50, 'P', P, 'TX22')))
print("{0:14.8f}".format(CP.Props('S', 'T', T + 50, 'P', P, 'TX22')))
print("")
# print "HCB: "
# print "{0:14.8f}".format(CP.Props('H','T',T+50,'P',P,'HCB'))
# print "{0:14.8f}".format(CP.Props('U','T',T+50,'P',P,'HCB'))
# print "{0:14.8f}".format(CP.Props('S','T',T+50,'P',P,'HCB'))
# print
# print CP.Props('T','D',13,'P',500,'n-Pentane')
# print 378.29349366868433
# print
# print CP.Props('D','T',373,'P',500,'n-Pentane')
# print  13.28572300574231
# print
# print CP.Props('C','T',375,'P',1000,'TX22')
# print 2.1829549872783205
# print
# print CP.Props('H','T',375,'P',1000,'TX22')
# print 158.1306573393763
# print
# print CP.Props('C','H',158,'P',1000,'TX22')
# print 2.18273570194
# print
# print CP.Props('C','T',275,'P',1000,'TestSolution-20%')
# print 4.0917018571
# print
# print CP.Props('H','T',275,'P',1000,'TestSolution-20%')
#print -93.759593807
# print
# print CP.Props('C','H',-94,'P',1000,'TestSolution-0.2')
# print 4.09169432982
# print
# print CP.Props('C','T',275,'P',1000,'EG-20%')
# print 4.07399544409
# print
# print CP.Props('H','T',275,'P',1000,'EG-20%')
#print -5.8439608317
# print
# print CP.Props('C','H',-5.8,'P',1000,'EG-20%')
# print 4.07402880961
