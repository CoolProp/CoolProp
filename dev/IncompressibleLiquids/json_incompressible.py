import numpy
from data_incompressible import SolutionData
    
class SolutionDataWriter(object):
    """ 
    A base class that defines all the variables needed 
    in order to make a proper fit. You can copy this code
    put in your data and add some documentation for where the
    information came from. 
    """
    #def __init__(self):
    #    SolutionData.__init__(self)
    
    def fitAll(self):
    
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
        
        dump = json.dumps(jobj, indent = 2, sort_keys = True)
        
        print dump
                
        fp = open(jobj['name']+'.json', 'w')
        fp.write(dump)
        fp.close()
            
    
    
if __name__ == '__main__':
    writer = SolutionDataWriter()
    
#    data = SecCoolExample()    
#    writer.toJSON(data)

    data = PureExample()
    data.density.coeffs = numpy.zeros((4,1))
    data.density.type = data.density.INCOMPRESSIBLE_POLYNOMIAL
    
    data.density.data = data.density.data[0:4]
    data.temperature.data = data.temperature.data[0:4]
    
    data.density.fit(data.temperature.data)
    #writer.toJSON(data)
    
    print data.density.data[0][0]
    print numpy.polynomial.polynomial.polyval2d(data.temperature.data[0], 0, data.density.coeffs)
    
    print data.density.data[1][0]
    print numpy.polynomial.polynomial.polyval2d(data.temperature.data[1], 0, data.density.coeffs)
    
    
    data = SolutionExample()
    data.density.coeffs = numpy.zeros((3,2))
    data.density.type = data.density.INCOMPRESSIBLE_POLYNOMIAL
    
    data.density.data = data.density.data[0:4][:,0:2]
    data.temperature.data = data.temperature.data[0:4]
    
    data.density.fit(data.temperature.data,data.concentration.data[0:2])
    #writer.toJSON(data)
    
    print data.density.data[0][0]
    print numpy.polynomial.polynomial.polyval2d(data.temperature.data[0], data.concentration.data[0], data.density.coeffs)
    
    print data.density.data[1][1]
    print numpy.polynomial.polynomial.polyval2d(data.temperature.data[1], data.concentration.data[1], data.density.coeffs)
    
    
    

    
    
    
    
    
    
    
    
    
