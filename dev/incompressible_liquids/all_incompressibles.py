from __future__ import division, absolute_import, print_function
from CPIncomp.WriterObjects import SolutionDataWriter
from CPIncomp import getExampleNames, getPureFluids, getCoefficientFluids,\
    getDigitalFluids, getSecCoolFluids, getMelinderFluids
import sys


if __name__ == '__main__':   
    
    writer = SolutionDataWriter()
    doneObjs = []
    
    # To debug single fluids
    from CPIncomp.SecCoolFluids import SecCoolSolutionData,SecCoolIceData
    from CPIncomp.PureFluids import Therminol72
    #solObjs  = [SecCoolSolutionData(sFile='Melinder, Ammonia'            ,sFolder='xMass',name='MAM2',desc='Melinder, Ammonia'            ,ref='Melinder-BOOK-2010, SecCool software')]
    #solObjs += [SecCoolIceData(sFile='IceNA'   ,sFolder='xMass',name='IceNA',desc='Ice slurry with NaCl' ,ref='Danish Technological Institute, SecCool software')]
    #solObjs = [Freezium()]
    #solObjs[0].density.DEBUG = True
    #solObjs[0].specific_heat.DEBUG = True
    #solObjs[0].conductivity.DEBUG = True
    #solObjs[0].viscosity.DEBUG = True
    #solObjs[0].T_freeze.DEBUG = True
    #writer.fitSecCoolList(solObjs)
    solObjs = [Therminol72()]
    solObjs[0].viscosity.DEBUG=True
    #solObjs[0].saturation_pressure.DEBUG=True
    writer.fitFluidList([solObjs[-1]])
    #
    ##from CPIncomp.ExampleObjects import SecCoolExample
    ##solObjs = [SecCoolExample()]
    writer.writeFluidList(solObjs)
    writer.writeReportList(solObjs)
    sys.exit(0)
    
    fluidObjs = getExampleNames(obj=True)
    examplesToFit = ["ExamplePure","ExampleSolution","ExampleDigital","ExampleDigitalPure"]
    for obj in fluidObjs:
        if obj.name in examplesToFit:
            writer.fitAll(obj)
    doneObjs += fluidObjs[:]
    
    print("\nProcessing example fluids")
    writer.writeFluidList(doneObjs)
    
    writer.writeReportList(doneObjs)
    #sys.exit(0)
        
    # If the examples did not cause any errors, 
    # we can proceed to the real data.
    doneObjs = []
    
    print("\nProcessing fluids with given coefficients")
    fluidObjs = getCoefficientFluids()
    doneObjs += fluidObjs[:]
    
    print("\nProcessing digital fluids")
    fluidObjs = getDigitalFluids()
    writer.fitFluidList(fluidObjs)
    doneObjs += fluidObjs[:]
    
    print("\nProcessing Melinder fluids")
    fluidObjs = getMelinderFluids()
    doneObjs += fluidObjs[:]
    
    print("\nProcessing pure fluids")
    fluidObjs = getPureFluids()
    writer.fitFluidList(fluidObjs)
    doneObjs += fluidObjs[:]
    
    print("\nProcessing SecCool fluids")
    fluidObjs = getSecCoolFluids()
    writer.fitSecCoolList(fluidObjs)
    doneObjs += fluidObjs[:]
    
    print("\nAll {0} fluids processed, all coefficients should be set.".format(len(doneObjs)))
    print("Checking the list of fluid objects.")
    #doneObjs += getCoefficientObjects()[:]
    doneObjs = sorted(doneObjs, key=lambda x: x.name)
    for i in range(len(doneObjs)-1):
        if doneObjs[i].name==doneObjs[i+1].name:
            print("Conflict between {0} and {1}, aborting".format(doneObjs[i],doneObjs[i+1]))
            raise ValueError("Two elements have the same name, that does not work: {0}".format(doneObjs[i].name))
    
    
    print("All checks passed, going to write to disk.")
    
    writer.writeFluidList(doneObjs)
    writer.writeReportList(doneObjs)

    print("All done, bye")
    sys.exit(0)
    
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
    