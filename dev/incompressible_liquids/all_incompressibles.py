from __future__ import division, absolute_import, print_function
import inspect
import numpy as np
import itertools,scipy.interpolate

import CoolProp.CoolProp as CP

import CPIncomp.PureFluids
import CPIncomp.DataObjects

from CPIncomp.WriterObjects import SolutionDataWriter
from CPIncomp.DataObjects import PureExample, SolutionExample
from CPIncomp.CoefficientObjects import SecCoolExample, MelinderExample

def getExampleData():
    return [PureExample(), SolutionExample()]

def getExampleCoef():
    return [SecCoolExample(), MelinderExample()]

def getExampleObjects():
    return getExampleData() + getExampleCoef()

def getBaseClassNames():
    ignList = []
    for i in inspect.getmembers(CPIncomp.DataObjects):
        ignList.append(i[0])
    return ignList

def getPureDataObjects():
    classes = []
    ignList = getBaseClassNames()

    for name, obj in inspect.getmembers(CPIncomp.PureFluids):
        if inspect.isclass(obj):
            #print(name)
            if not name in ignList: # Ignore the base classes
                classes += [obj()]
    return classes 

def getSolutionDataObjects():
    return []
    classes = []
    ignList = getBaseClassNames()

    for name, obj in inspect.getmembers(CPIncomp.SolutionFluids):
        if inspect.isclass(obj):
            #print(name)
            if not name in ignList: # Ignore the base classes
                classes += [obj()]
    return classes 

def getCoefficientObjects():
    return []
    classes = []
    ignList = getBaseClassNames()

    for name, obj in inspect.getmembers(CPIncomp.CoefficientFluids):
        if inspect.isclass(obj):
            #print(name)
            if not name in ignList: # Ignore the base classes
                classes += [obj()]
    return classes 


def fitFluidList(fluidObjs):
    for obj in fluidObjs:
        if obj==fluidObjs[0]:
            print(" {0}".format(obj.name), end="")
        elif obj==fluidObjs[-1]:
            print(", {0}".format(obj.name), end="")
        else:
            print(", {0}".format(obj.name), end="")
            
        try: 
            writer.fitAll(obj)
        except (TypeError, ValueError) as e:
            print("An error occurred for fluid: {0}".format(obj.name))
            print(obj)
            print(e)
            pass 
    return 

if __name__ == '__main__':   
    
    writer = SolutionDataWriter()
    
    dataObjs = getExampleData()
    for obj in dataObjs:
        writer.fitAll(obj)
        
    dataObjs += getExampleCoef()
    for obj in dataObjs:
        writer.toJSON(obj)
        
    # If the examples did not cause any errors, 
    # we can proceed to the real data.
    dataObjs = getPureDataObjects()
    print("Fitting pure fluids:".format(obj.name), end="")
    fitFluidList(dataObjs)
    print(" ... done")
    for obj in dataObjs:
        writer.toJSON(obj)
    
    dataObjs = getSolutionDataObjects()
    print("Fitting solutions:".format(obj.name), end="")
    fitFluidList(dataObjs)
    print(" ... done")
    for obj in dataObjs:
        writer.toJSON(obj)
    
#    data = SecCoolExample()    
#    writer.toJSON(data)

#    data = PureExample()
#    data.density.coeffs = np.zeros((4,1))
#    data.density.type = data.density.INCOMPRESSIBLE_POLYNOMIAL
#    data.density.data = data.density.data[0:4]
#    data.temperature.data = data.temperature.data[0:4]
#    data.density.fit(data.temperature.data)
#    #writer.toJSON(data)
#    print data.density.data[0][0]
#    print np.polynomial.polynomial.polyval2d(data.temperature.data[0], 0, data.density.coeffs)
#    print data.density.data[1][0]
#    print np.polynomial.polynomial.polyval2d(data.temperature.data[1], 0, data.density.coeffs)
    
    
#    data = SolutionExample()
#    data.density.coeffs = np.zeros((3,2))
#    data.density.type = data.density.INCOMPRESSIBLE_POLYNOMIAL
#    data.density.data = data.density.data[0:4][:,0:2]
#    data.temperature.data = data.temperature.data[0:4]
#    data.density.fit(data.temperature.data,data.concentration.data[0:2])
#    #writer.toJSON(data)
#    print data.density.data[0][0]
#    print np.polynomial.polynomial.polyval2d(data.temperature.data[0], data.concentration.data[0], data.density.coeffs)
#    print data.density.data[1][1]
#    print np.polynomial.polynomial.polyval2d(data.temperature.data[1], data.concentration.data[1], data.density.coeffs)

#    test = True  
#    #if test: import CoolProp.CoolProp as CP
#    if test: from scipy import interpolate
#    if test: p = 10e5
#    
#    def printInfo(data):
#        print("{0:s} : {1:.4e}, {2:.4e}".format(data.name, data.Tbase, data.xbase))
#        
#    def printValue(data, T, p, x, fluid='', f=None, dataFunc=None, dataLetter=''):
#        if f!=None:
#            try:
#                print("{0:s} : {7:s} : {1:.4e}, {2:.4e}, {3:.4e}, inputs: {4:.4e}, {5:.4e}, {6:.4e} ".format(data.name, dataFunc(T, p, x), CP.PropsSI(dataLetter,'T',T,'P',p,fluid), float(f(T,x)), T, p, x, dataLetter))
#            except:
#                print("{0:s} : {7:s} : {1:.4e}, {2:.4e}, {3:.4e}, inputs: {4:.4e}, {5:.4e}, {6:.4e} ".format(data.name, dataFunc(T, p, x), CP.PropsSI(dataLetter,'T',T,'P',p,fluid), float(f(T)), T, p, x, dataLetter))
#        else: 
#            print    ("{0:s} : {7:s} : {1:.4e}, {2:.4e}, {3:.4e}, inputs: {4:.4e}, {5:.4e}, {6:.4e} ".format(data.name, dataFunc(T, p, x), CP.PropsSI(dataLetter,'T',T,'P',p,fluid), 0.0, T, p, x, dataLetter))
#    
#    def printDens(data, T, p, x, fluid='', f=None):
#        printValue(data, T, p, x, fluid=fluid, f=f, dataFunc=data.rho, dataLetter='D')
#    
#    def printHeat(data, T, p, x, fluid='', f=None):
#        printValue(data, T, p, x, fluid=fluid, f=f, dataFunc=data.c, dataLetter='C')
#        
#    def printEnergy(data, T, p, x, fluid='', f=None):
#        printValue(data, T, p, x, fluid=fluid, f=f, dataFunc=data.u, dataLetter='U')
#
#    def printEnthalpy(data, T, p, x, fluid='', f=None):
#        printValue(data, T, p, x, fluid=fluid, f=f, dataFunc=data.h, dataLetter='H')
#    
#    def printVisc(data, T, p, x, fluid='', f=None):
#        printValue(data, T, p, x, fluid=fluid, f=f, dataFunc=data.visc, dataLetter='V')
#    
#    def printCond(data, T, p, x, fluid='', f=None):
#        printValue(data, T, p, x, fluid=fluid, f=f, dataFunc=data.cond, dataLetter='L')
#    
#    def printPsat(data, T, p, x, fluid='', f=None):
#        printValue(data, T, p, x, fluid=fluid, f=f, dataFunc=data.psat, dataLetter='Psat')
#        
#    def printTfreeze(data, T, p, x, fluid='', f=None):
#        printValue(data, T, p, x, fluid=fluid, f=f, dataFunc=data.Tfreeze, dataLetter='Tfreeze')
#        
#    def printAll(data, T, p, x, fluid=''):
#        printDens(data, T, p, x, fluid, None)
#        printHeat(data, T, p, x, fluid, None)
#        printEnergy(data, T, p, x, fluid, None)
#        #printEnthalpy(data, T, p, x, fluid, None)
#        printVisc(data, T, p, x, fluid, None)
#        printCond(data, T, p, x, fluid, None)
#        #printPsat(data, T, p, x, fluid, None)
#        #printTfreeze(data, T, p, x, fluid, None)
#    
# 
#    
#    
#    data = PureExample()
#    writer.fitAll(data)
#    writer.toJSON(data)
#    #printInfo(data)
#    if test: T = 55+273.15
#    if test: x = 0.0 
#    #if test: f = interpolate.interp1d(data.temperature.data, data.density.data.T[0])
#    #if test: printDens(data, T, p, x, fluid='TD12', f=f)
#    #if test: f = interpolate.interp1d(data.temperature.data, data.specific_heat.data.T[0])
#    #if test: printHeat(data, T, p, x, fluid='TD12', f=f)
#    if test: printAll(data, T, p, x, fluid='TD12')
#   
#    
#    data = SolutionExample()
#    writer.fitAll(data)
#    writer.toJSON(data)
#    #printInfo(data)
#    if test: T = -15+273.15
#    if test: x = 0.10
#    #if test: f = interpolate.interp2d(data.temperature.data, data.concentration.data, data.density.data.T)
#    #if test: printDens(data, T, p, x, fluid='IceEA-{0:.4f}%'.format(x*100.0), f=f)
#    #if test: f = None
#    if test: printAll(data, T, p, x, fluid='IceEA-{0:.4f}%'.format(x*100.0))
#
#        
#    data = SecCoolExample()
#    writer.toJSON(data)
#    #printInfo(data)
#    if test: T = -5+273.15
#    if test: x = 0.40
#    if test: f = None #interpolate.interp2d(data.temperature.data, data.concentration.data, data.density.data.T)
#    if test: printAll(data, T, p, x, fluid='SecCoolSolution-{0:.4f}%'.format(x*100.0))    
#    
#    data = MelinderExample()
#    writer.toJSON(data)
#    #printInfo(data)
#    if test: T = -5+273.15
#    if test: x = 0.3
#    if test: f = None #interpolate.interp2d(data.temperature.data, data.concentration.data, data.density.data.T)
#    if test: printAll(data, T, p, x, fluid='MMA-{0:.4f}%'.format(x*100.0))   
#    
    