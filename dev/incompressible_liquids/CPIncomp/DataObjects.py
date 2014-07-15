from __future__ import division, absolute_import, print_function
import numpy as np
from CPIncomp.BaseObjects import IncompressibleData


class SolutionData(object):
    """ 
    A base class that defines all the variables needed 
    in order to make a proper fit. You can copy this code
    put in your data and add some documentation for where the
    information came from. 
    """
    def __init__(self):
        self.name        = None # Name of the current fluid
        self.description = None # Description of the current fluid
        self.reference   = None # Reference data for the current fluid
        
        self.Tmax        = None # Maximum temperature in K
        self.Tmin        = None # Minimum temperature in K
        self.xmax        = 1.0 # Maximum concentration
        self.xmin        = 0.0 # Minimum concentration
        self.TminPsat    = None # Minimum saturation temperature in K
        self.Tbase       = 0.0  # Base value for temperature fits
        self.xbase       = 0.0  # Base value for concentration fits
    
        self.temperature   = IncompressibleData() # Temperature for data points in K
        self.concentration = IncompressibleData() # Concentration data points in weight fraction
        self.density       = IncompressibleData() # Density in kg/m3
        self.specific_heat = IncompressibleData() # Heat capacity in J/(kg.K)
        self.viscosity     = IncompressibleData() # Dynamic viscosity in Pa.s
        self.conductivity  = IncompressibleData() # Thermal conductivity in W/(m.K)
        self.saturation_pressure = IncompressibleData() # Saturation pressure in Pa
        self.T_freeze      = IncompressibleData() # Freezing temperature in K
        self.volume2mass   = IncompressibleData() # dd
        self.mass2mole     = IncompressibleData() # dd
        
        # Some of the functions might need a guess array
        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPONENTIAL
        self.viscosity.coeffs = np.array([+7e+2, -6e+1, +1e+1])
        self.saturation_pressure.type = self.saturation_pressure.INCOMPRESSIBLE_EXPONENTIAL
        self.saturation_pressure.coeffs = np.array([-5e+3, +3e+1, -1e+1])
        
        self.xref = 0.0
        self.Tref = 0.0
        self.pref = 0.0
        self.href = 0.0
        self.sref = 0.0
        self.uref = 0.0
        self.rhoref = 0.0
        
    
    def rho (self, T, p, x=0.0, c=None):
        if c==None: 
            c=self.density.coeffs
        if self.density.type==self.density.INCOMPRESSIBLE_POLYNOMIAL:
            return np.polynomial.polynomial.polyval2d(T-self.Tbase, x-self.xbase, c)
        else:  raise ValueError("Unknown function.")
    
    def c   (self, T, p, x=0.0, c=None):
        if c==None: 
            c = self.specific_heat.coeffs
        if self.specific_heat.type==self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL:
            return np.polynomial.polynomial.polyval2d(T-self.Tbase, x-self.xbase, c)
        else:  raise ValueError("Unknown function.")
    
    def cp  (self, T, p, x=0.0, c=None):
        return self.c(T,p,x,c)
    
    def cv  (self, T, p, x=0.0, c=None):
        return self.c(T,p,x,c)
    
    def u   (self, T, p, x=0.0, c=None):
        if c==None: 
            c = self.specific_heat.coeffs
        if self.specific_heat.type==self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL:
            c_tmp = np.polynomial.polynomial.polyint(c)
            return np.polynomial.polynomial.polyval2d(T-self.Tbase, x-self.xbase, c_tmp)
        else:  raise ValueError("Unknown function.")
    #def u   (T, p, x);
    
    def h   (self, T, p, x=0.0):
        return self.h_u(T,p,x)
    
    def visc(self, T, p, x=0.0, c=None):
        return self.viscosity.baseFunction(T, x, self.Tbase, self.xbase, c=c)
    
    def cond(self, T, p, x=0.0, c=None):
        return self.conductivity.baseFunction(T, x, self.Tbase, self.xbase, c=c)
        
    def psat(self, T, x=0.0, c=None):
        return self.saturation_pressure.baseFunction(T, x, self.Tbase, self.xbase, c=c)
        
    def Tfreeze(self, p, x=0.0, c=None):
        return self.T_freeze.baseFunction(x, 0.0, self.xbase, 0.0, c=c)
    
    #def V2M (T,           y);
    #def M2M (T,           x);
    
    
    def h_u(self, T, p, x):
        return self.u(T,p,x)+p/self.rho(T,p,x)-self.href

    def u_h(self, T, p, x):
        return self.h(T,p,x)-p/self.rho(T,p,x)+self.href
    
    def set_reference_state(self, T0, p0, x0=0.0, h0=0.0, s0=0.0):
        self.rhoref = self.rho(T0,p0,x0)
        self.pref = p0
        self.uref = h0 - p0/self.rhoref
        self.uref = self.u(T0,p0,x0)
        self.href = h0
        self.sref = s0
        self.href = self.h(T0,p0,x0)
        self.sref = self.s(T0,p0,x0)
        
    
class PureExample(SolutionData):
    def __init__(self):
        SolutionData.__init__(self) 
        self.name = "ExamplePure"
        self.description = "Heat transfer fluid TherminolD12 by Solutia"
        self.reference = "Solutia data sheet"
        self.Tmax = 150 + 273.15
        self.Tmin =  50 + 273.15
        self.TminPsat =  self.Tmax
        
        self.temperature.data         = np.array([    50 ,     60 ,     70 ,     80 ,     90 ,    100 ,    110 ,    120 ,    130 ,    140 ,    150 ])+273.15 # Kelvin
        self.density.data             = np.array([[  740],[   733],[   726],[   717],[   710],[   702],[   695],[   687],[   679],[   670],[   662]])        # kg/m3
        self.specific_heat.data       = np.array([[ 2235],[  2280],[  2326],[  2361],[  2406],[  2445],[  2485],[  2528],[  2571],[  2607],[  2645]])        # J/kg-K
        self.viscosity.data           = np.array([[0.804],[ 0.704],[ 0.623],[ 0.556],[ 0.498],[ 0.451],[ 0.410],[ 0.374],[ 0.346],[ 0.317],[ 0.289]])        # Pa-s
        self.conductivity.data        = np.array([[0.105],[ 0.104],[ 0.102],[ 0.100],[ 0.098],[ 0.096],[ 0.095],[ 0.093],[ 0.091],[ 0.089],[ 0.087]])        # W/m-K
        self.saturation_pressure.data = np.array([[  0.5],[   0.9],[   1.4],[   2.3],[   3.9],[   6.0],[   8.7],[  12.4],[  17.6],[  24.4],[  33.2]])        # Pa


class SolutionExample(SolutionData):
    def __init__(self):
        SolutionData.__init__(self) 
        self.name = "ExampleSolution"
        self.description = "Ethanol ice slurry"
        self.reference = "SecCool software"
        self.Tmax = -10 + 273.15
        self.Tmin = -45 + 273.15
        self.TminPsat =  self.Tmax
        
        self.temperature.data         = np.array([   -45 ,    -40 ,    -35 ,    -30 ,    -25 ,    -20 ,    -15 ,    -10])+273.15 # Kelvin
        self.concentration.data       = np.array([     5 ,     10 ,     15 ,     20 ,     25 ,     30 ,     35 ])/100.0 # mass fraction
        
        self.density.data             = np.array([
          [1064.0,    1054.6,    1045.3,    1036.3,    1027.4,    1018.6,    1010.0],
          [1061.3,    1052.1,    1043.1,    1034.3,    1025.6,    1017.0,    1008.6],
          [1057.6,    1048.8,    1040.1,    1031.5,    1023.1,    1014.8,    1006.7],
          [1053.1,    1044.6,    1036.2,    1028.0,    1019.9,    1012.0,    1004.1],
          [1047.5,    1039.4,    1031.5,    1023.7,    1016.0,    1008.4,    1000.9],
          [1040.7,    1033.2,    1025.7,    1018.4,    1011.2,    1004.0,     997.0],
          [1032.3,    1025.3,    1018.5,    1011.7,    1005.1,     998.5,     992.0],
          [1021.5,    1015.3,    1009.2,    1003.1,     997.1,     991.2,     985.4]]) # kg/m3
        
#        self.density.data             = np.array([
#          [np.nan,    np.nan,    np.nan,    np.nan,    np.nan,    np.nan,    np.nan],
#          [np.nan,    np.nan,    np.nan,    np.nan,    np.nan,    np.nan,    np.nan],
#          [np.nan,    np.nan,    np.nan,    np.nan,    np.nan,    np.nan,    np.nan],
#          [np.nan,    np.nan,    np.nan,    np.nan,    np.nan,    np.nan,    np.nan],
#          [np.nan,    np.nan,    np.nan,    np.nan,    1016.0,    1008.4,    np.nan],
#          [np.nan,    1033.2,    np.nan,    np.nan,    np.nan,    np.nan,    np.nan],
#          [np.nan,    1025.3,    1018.5,    np.nan,    np.nan,     998.5,     992.0],
#          [np.nan,    np.nan,    1009.2,    np.nan,    np.nan,    np.nan,    np.nan]]) # kg/m3


if __name__ == '__main__':
    obj = PureExample()
    obj.saturation_pressure.type = obj.saturation_pressure.INCOMPRESSIBLE_EXPONENTIAL
    print(obj.saturation_pressure.coeffs)
    xData = obj.temperature.data
    coeffs = obj.saturation_pressure.getCoeffsIterative1D(xData)
    obj.saturation_pressure.coeffs = coeffs
    print(obj.saturation_pressure.coeffs)
    
    obj = PureExample()
    print(obj.saturation_pressure.coeffs)
    obj.saturation_pressure.fit(obj.temperature.data)
    print(obj.saturation_pressure.coeffs)
    