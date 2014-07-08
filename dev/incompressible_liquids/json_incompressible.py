import numpy as np
from data_incompressible import *
    
class SolutionDataWriter(object):
    """ 
    A base class that defines all the variables needed 
    in order to make a proper fit. You can copy this code
    put in your data and add some documentation for where the
    information came from. 
    """
    def __init__(self):
        self.verbose = False
    
    def fitAll(self, data):
        T = data.temperature.data
        x = data.concentration.data
        
        #if data.Tbase==0.0:
        #    data.Tbase = (np.min(T) + np.max(T)) / 2.0
        #if data.xbase==0.0:
        #    data.xbase = (np.min(x) + np.max(x)) / 2.0
            
        # Set the standard order for polynomials
        std_xorder = 3
        std_yorder = 5
        std_coeffs = np.zeros((std_xorder,std_yorder))
        
        errList = (ValueError, AttributeError, TypeError)
        
        try:
            data.density.coeffs = std_coeffs[:]
            data.density.type   = data.density.INCOMPRESSIBLE_POLYNOMIAL
            data.density.fit(T,x,data.Tbase,data.xbase)
        except errList as ve:
            if self.verbose: print "Could not fit density coefficients:", ve
            pass

        try:        
            data.specific_heat.coeffs = std_coeffs[:]
            data.specific_heat.type   = data.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
            data.specific_heat.fit(T,x,data.Tbase,data.xbase)
        except errList as ve:
            if self.verbose: print "Could not fit specific heat coefficients:", ve
            pass
        
        try:
            data.viscosity.coeffs = std_coeffs[:]
            data.viscosity.type   = data.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
            data.viscosity.fit(T,x,data.Tbase,data.xbase)
        except errList as ve:
            if self.verbose: print "Could not fit viscosity coefficients:", ve
            pass

        try:        
            data.conductivity.coeffs = std_coeffs[:]
            data.conductivity.type   = data.conductivity.INCOMPRESSIBLE_POLYNOMIAL
            data.conductivity.fit(T,x,data.Tbase,data.xbase)
        except errList as ve:
            if self.verbose: print "Could not fit conductivity coefficients:", ve
            pass
        
        try:
            data.saturation_pressure.coeffs = std_coeffs[:]
            data.saturation_pressure.type   = data.saturation_pressure.INCOMPRESSIBLE_EXPPOLYNOMIAL
            data.saturation_pressure.fit(T,x,data.Tbase,data.xbase)
        except errList as ve:
            if self.verbose: print "Could not fit saturation pressure coefficients:", ve
            pass

        try:        
            data.T_freeze.coeffs = std_coeffs[:]
            data.T_freeze.type   = data.T_freeze.INCOMPRESSIBLE_POLYNOMIAL
            data.T_freeze.fit(T,x,data.Tbase,data.xbase)
        except errList as ve:
            if self.verbose: print "Could not fit TFreeze coefficients:", ve
            pass

        try:        
            data.volume2mass.coeffs = std_coeffs[:]
            data.volume2mass.type   = data.volume2mass.INCOMPRESSIBLE_POLYNOMIAL
            data.volume2mass.fit(T,x,data.Tbase,data.xbase)
        except errList as ve:
            if self.verbose: print "Could not fit V2M coefficients:", ve
            pass

        try:        
            data.mass2mole.coeffs = std_coeffs[:]
            data.mass2mole.type   = data.mass2mole.INCOMPRESSIBLE_POLYNOMIAL
            data.mass2mole.fit(T,x,data.Tbase,data.xbase)
        except errList as ve:
            if self.verbose: print "Could not fit M2M coefficients:", ve
            pass
        
    
    def toJSON(self,data):
        jobj = {}
        
        jobj['name'] = data.name # Name of the current fluid
        jobj['description'] = data.description # Description of the current fluid
        jobj['reference'] = data.reference # Reference data for the current fluid
        
        jobj['Tmax'] = data.Tmax         # Maximum temperature in K
        jobj['Tmin'] = data.Tmin         # Minimum temperature in K
        jobj['xmax'] = data.xmax         # Maximum concentration
        jobj['xmin'] = data.xmin         # Minimum concentration
        jobj['TminPsat'] = data.TminPsat     # Minimum saturation temperature in K
        jobj['Tbase'] = data.Tbase       # Base value for temperature fits
        jobj['xbase'] = data.xbase       # Base value for concentration fits
    
        #data.temperature  # Temperature for data points in K
        #data.concentration  # Concentration data points in weight fraction
        jobj['density'] = data.density.toJSON() # Density in kg/m3
        jobj['specific_heat'] = data.specific_heat.toJSON() # Heat capacity in J/(kg.K)
        jobj['viscosity'] = data.viscosity.toJSON() # Dynamic viscosity in Pa.s
        jobj['conductivity'] = data.conductivity.toJSON() # Thermal conductivity in W/(m.K)
        jobj['saturation_pressure'] = data.saturation_pressure.toJSON() # Saturation pressure in Pa
        jobj['T_freeze'] = data.T_freeze.toJSON() # Freezing temperature in K
        jobj['volume2mass'] = data.volume2mass.toJSON() # dd
        jobj['mass2mole'] = data.mass2mole.toJSON() # dd
        
        import json
        
        
        original_float_repr = json.encoder.FLOAT_REPR 
        json.encoder.FLOAT_REPR = lambda o: format(o, ' .8e') 
        #print json.dumps(1.0001) 
        dump = json.dumps(jobj, indent = 2, sort_keys = True)
        json.encoder.FLOAT_REPR = original_float_repr
        
        #print dump
                
        fp = open(jobj['name']+'.json', 'w')
        fp.write(dump)
        fp.close()
            
    
    
if __name__ == '__main__':
    writer = SolutionDataWriter()
    
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

    test = True  
    if test: import CoolProp.CoolProp as CP
    if test: from scipy import interpolate
    if test: p = 10e5
    
    def printInfo(data):
        print "{0:s} : {1:.4e}, {2:.4e}".format(data.name, data.Tbase, data.xbase)
        
    def printValue(data, T, p, x, fluid='', f=None, dataFunc=None, dataLetter=''):
        if f!=None:
            try:
                print "{0:s} : {1:.4e}, {2:.4e}, {3:.4e}, inputs: {4:.4e}, {5:.4e}, {6:.4e} ".format(data.name, dataFunc(T, p, x), CP.Props(dataLetter,'T',T,'P',p,fluid), float(f(T,x)), T, p, x)
            except:
                print "{0:s} : {1:.4e}, {2:.4e}, {3:.4e}, inputs: {4:.4e}, {5:.4e}, {6:.4e} ".format(data.name, dataFunc(T, p, x), CP.Props(dataLetter,'T',T,'P',p,fluid), float(f(T)), T, p, x)
        else: 
            print "{0:s} : {1:.4e}, {2:.4e}, {3:.4e}, inputs: {4:.4e}, {5:.4e}, {6:.4e} ".format(data.name, dataFunc(T, p, x), CP.Props(dataLetter,'T',T,'P',p,fluid), 0.0, T, p, x)
    
    def printDens(data, T, p, x, fluid='', f=None):
        printValue(data, T, p, x, fluid=fluid, f=f, dataFunc=data.rho, dataLetter='D')
    
    def printHeat(data, T, p, x, fluid='', f=None):
        printValue(data, T, p, x, fluid=fluid, f=f, dataFunc=data.c, dataLetter='C')
        
    
    data = PureExample()
    writer.fitAll(data)
    writer.toJSON(data)
    printInfo(data)
    if test: T = data.Tbase+55+273.15
    if test: x = 0.0 
    if test: f = interpolate.interp1d(data.temperature.data, data.density.data.T[0])
    if test: printDens(data, T, p, x, fluid='TD12', f=f)
   
    
    data = SolutionExample()
    writer.fitAll(data)
    writer.toJSON(data)
    printInfo(data)
    if test: T = data.Tbase-15+273.15
    if test: x = 0.22 
    if test: f = interpolate.interp2d(data.temperature.data, data.concentration.data, data.density.data.T)
    if test: printDens(data, T, p, x, fluid='IceEA-{0:.4f}%'.format(x*100.0), f=f)

        
    data = SecCoolExample()
    writer.toJSON(data)
    printInfo(data)
    if test: T = data.Tbase+0
    if test: x = 0.3157
    if test: f = None #interpolate.interp2d(data.temperature.data, data.concentration.data, data.density.data.T)
    if test: printDens(data, T, p, x, fluid='SecCoolSolution-{0:.4f}%'.format(x*100.0), f=f)
    
    print 
    print "It works for temperature changes:"
    printDens(data, T-10, p, x, fluid='SecCoolSolution-{0:.4f}%'.format(x*100.0), f=f)
    printDens(data, T+10, p, x, fluid='SecCoolSolution-{0:.4f}%'.format(x*100.0), f=f)
    printHeat(data, T-10, p, x, fluid='SecCoolSolution-{0:.4f}%'.format(x*100.0), f=f)
    printHeat(data, T+10, p, x, fluid='SecCoolSolution-{0:.4f}%'.format(x*100.0), f=f)
    print "but not for concentration changes:"
    x = x-0.1
    printDens(data, T, p, x, fluid='SecCoolSolution-{0:.4f}%'.format(x*100.0), f=f)
    x = x+0.2
    printDens(data, T, p, x, fluid='SecCoolSolution-{0:.4f}%'.format(x*100.0), f=f)
    x = x-0.2
    printHeat(data, T, p, x, fluid='SecCoolSolution-{0:.4f}%'.format(x*100.0), f=f)
    x = x+0.2
    printHeat(data, T, p, x, fluid='SecCoolSolution-{0:.4f}%'.format(x*100.0), f=f)
    print 
    
    
    data = MelinderExample()
    writer.toJSON(data)
    printInfo(data)
    if test: T = data.Tbase+10
    if test: x = 0.22 
    if test: f = None #interpolate.interp2d(data.temperature.data, data.concentration.data, data.density.data.T)
    if test: printDens(data, T, p, x, fluid='MMA-22%', f=f)
    
    
    
    

    
    
    
    
    
    
    
    
    
