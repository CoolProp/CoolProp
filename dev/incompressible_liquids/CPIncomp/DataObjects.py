from __future__ import division, absolute_import, print_function
import numpy as np
from scipy import interpolate
from CPIncomp.BaseObjects import IncompressibleData
import os, CPIncomp, math

class SolutionData(object):
    """ 
    A base class that defines all the variables needed 
    in order to make a proper fit. You can copy this code
    put in your data and add some documentation for where the
    information came from. 
    """
    def __init__(self):
        self.ifrac_mass      = 0
        self.ifrac_mole      = 1
        self.ifrac_volume    = 2
        self.ifrac_undefined = 3
        self.ifrac_pure      = 4
        
        self.significantDigits = 6
        
        self.name        = None # Name of the current fluid
        self.description = None # Description of the current fluid
        self.reference   = None # Reference data for the current fluid
        
        self.Tmax        = None # Maximum temperature in K
        self.Tmin        = None # Minimum temperature in K
        self.xmax        = 1.0 # Maximum concentration
        self.xmin        = 0.0 # Minimum concentration
        self.xid         = self.ifrac_undefined # Concentration is mole, mass or volume-based
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
        self.mass2input    = IncompressibleData() # dd
        self.volume2input  = IncompressibleData() # dd
        self.mole2input    = IncompressibleData() # dd
        
        ## Some of the functions might need a guess array
        #self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPONENTIAL
        #self.viscosity.coeffs = np.array([+7e+2, -6e+1, +1e+1])
        #self.saturation_pressure.type = self.saturation_pressure.INCOMPRESSIBLE_EXPONENTIAL
        #self.saturation_pressure.coeffs = np.array([-5e+3, +3e+1, -1e+1])

        
        self.xref = 0.0
        self.Tref = 0.0
        self.pref = 0.0
        self.href = 0.0
        self.sref = 0.0
        self.uref = 0.0
        self.rhoref = 0.0
    
#    def getDataObjects(self):
#        objList  = {}
#        objList["temperature"] = self.temperature
#        objList["concentration"] = self.concentration
#        objList["density"] = self.density
#        objList["specific_heat"] = self.specific_heat
#        objList["viscosity"] = self.viscosity
#        objList["conductivity"] = self.conductivity
#        objList["saturation_pressure"] = self.saturation_pressure
#        objList["T_freeze"] = self.T_freeze
#        objList["volume2mass"] = self.volume2mass
#        objList["mass2mole"] = self.mass2mole
#        return objList

    def roundSingle(self,x):
        if x==0.0: return 0.0
        return round(x, self.significantDigits-int(math.floor(math.log10(abs(x))))-1) 

    def round(self, x):
        r,c,res = IncompressibleData.shapeArray(x)
        #digits = -1*np.floor(np.log10(res))+self.significantDigits-1 
        for i in range(r):
            for j in range(c):
                if np.isfinite(res[i,j]):
                    res[i,j] = self.roundSingle(res[i,j])
        return res

#    def getPolyObjects(self):
#        objList  = {}
#        objList["density"] = self.density
#        objList["specific heat"] = self.specific_heat
##        objList["viscosity"] = self.viscosity
#        objList["conductivity"] = self.conductivity
##        objList["saturation_pressure"] = self.saturation_pressure
##        objList["T_freeze"] = self.T_freeze
##        objList["volume2mass"] = self.volume2mass
##        objList["mass2mole"] = self.mass2mole
#        return objList
#    
#    def getExpPolyObjects(self):
#        objList  = {}
#        objList["viscosity"] = self.viscosity
#        objList["saturation pressure"] = self.saturation_pressure
#        return objList
        
    
    def rho (self, T, p=0.0, x=0.0, c=None):
        if c==None: 
            c=self.density.coeffs
        if self.density.type==self.density.INCOMPRESSIBLE_POLYNOMIAL:
            return np.polynomial.polynomial.polyval2d(T-self.Tbase, x-self.xbase, c)
        else:  raise ValueError("Unknown function.")
    
    def c   (self, T, p=0.0, x=0.0, c=None):
        if c==None: 
            c = self.specific_heat.coeffs
        if self.specific_heat.type==self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL:
            return np.polynomial.polynomial.polyval2d(T-self.Tbase, x-self.xbase, c)
        else:  raise ValueError("Unknown function.")
    
    def cp  (self, T, p=0.0, x=0.0, c=None):
        return self.c(T,p,x,c)
    
    def cv  (self, T, p=0.0, x=0.0, c=None):
        return self.c(T,p,x,c)
    
    def u   (self, T, p=0.0, x=0.0, c=None):
        if c==None: 
            c = self.specific_heat.coeffs
        if self.specific_heat.type==self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL:
            c_tmp = np.polynomial.polynomial.polyint(c)
            return np.polynomial.polynomial.polyval2d(T-self.Tbase, x-self.xbase, c_tmp)
        else:  raise ValueError("Unknown function.")
    #def u   (T, p, x);
    
    def h   (self, T, p=0.0, x=0.0):
        return self.h_u(T,p,x)
    
    def visc(self, T, p=0.0, x=0.0, c=None):
        return self.viscosity.baseFunction(T, x, self.Tbase, self.xbase, c=c)
    
    def cond(self, T, p=0.0, x=0.0, c=None):
        return self.conductivity.baseFunction(T, x, self.Tbase, self.xbase, c=c)
        
    def psat(self, T,        x=0.0, c=None):
        return self.saturation_pressure.baseFunction(T, x, self.Tbase, self.xbase, c=c)
        
    def Tfreeze(self, p=0.0, x=0.0, c=None):
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


class PureData(SolutionData):
    """ 
    An extension of the solution data that makes it 
    easier to gather data for pure fluids. 
    """
    def __init__(self):
        SolutionData.__init__(self)
        self.xid         = self.ifrac_pure
        self.concentration.data       = np.array([     0 ]) # mass fraction
        
    def reshapeData(self, dataArray, length):
        """
        Makes any array 1-dimensional and implicitly 
        checks the length.
        """
        if dataArray!=None:
            return np.reshape(dataArray, (length,1))
        else:
            return None
        
        
    def reshapeAll(self):
        len_T = len(self.temperature.data)
        #len_x = len(self.concentration.data)
        self.density.data       = self.reshapeData(self.density.data, len_T)
        self.specific_heat.data = self.reshapeData(self.specific_heat.data, len_T)
        self.viscosity.data     = self.reshapeData(self.viscosity.data, len_T)
        self.conductivity.data  = self.reshapeData(self.conductivity.data, len_T)
        self.saturation_pressure.data = self.reshapeData(self.saturation_pressure.data, len_T)
        #self.T_freeze.data      = self.reshapeData(self.T_freeze.data, len_T)
        #self.volume2mass.data   = self.reshapeData(self.volume2mass.data, len_T)
        #self.mass2mole.data     = self.reshapeData(self.mass2mole.data, len_T)        


class DigitalData(SolutionData):
    """ 
    An extension of the solution data that makes it 
    easier to generate fitting data from fluids available
    as Python packages.
    """
    def __init__(self):
        SolutionData.__init__(self)
        
    def getFile(self, data):
        return os.path.join(CPIncomp.__path__[0], 'data', self.name+"_"+data+".txt")
    
    def getFromFile(self, data):
        fullPath = self.getFile(data)
        _,_,res = IncompressibleData.shapeArray(np.loadtxt(fullPath))
        return res 
    
    def writeToFile(self, data, array):
        fullPath = self.getFile(data)
        stdFmt = "%1.{0}e".format(int(self.significantDigits-1))
        return np.savetxt(fullPath, array, fmt=stdFmt)
    
    def getTrange(self):
        if self.Tmin<self.Tmax:
            return np.linspace(self.Tmin, self.Tmax, 20)
        else:
            return np.array([0.0])
    
    def getxrange(self):
        if self.xmin<self.xmax:
            return np.linspace(self.xmin, self.xmax, 20)
        else:
            return np.array([0.0])
    
    def getArray(self, func=None, data=None):
        """
        func is a callable object that takes T,x as inputs
        and data is the file name for the data. 
        We try to read the file and if unsuccessful, we 
        generate the data and write it.
        """
        
        forceUpdate = False
        readFromFile = False
        fileArray = None
        
        if self.temperature.data==None or self.concentration.data==None: # no data set, try to get it from file
            if self.temperature.data!=None: raise ValueError("Temperature is not None, but concentration is.")
            if self.concentration.data!=None: raise ValueError("Concentration is not None, but temperature is.")
            if (os.path.isfile(self.getFile(data))): # File found
                fileArray = self.getFromFile(data)
                self.temperature.data = np.copy(fileArray[1:,0])
                self.concentration.data = np.copy(fileArray[0,1:])
                readFromFile = True
            else:
                raise ValueError("No temperature and concentration data given and no readable file found.")
            
        tData = self.round(self.temperature.data)[:,0]
        xData = self.round(self.concentration.data)[:,0]

        baseArray = np.zeros( (len(tData)+1,len(xData)+1) )
        
        if (os.path.isfile(self.getFile(data)) and not forceUpdate): # File found and no update wanted
            if fileArray==None: fileArray = self.getFromFile(data)
            
#            tFile = fileArray[1:,0]
#            xFile = fileArray[0,1:]
#            curDataLim = np.array([np.min(tData),np.max(tData),np.min(xData),np.max(xData)])
#            curFileLim = np.array([np.min(tFile),np.max(tFile),np.min(xFile),np.max(xFile)])
#            if np.allclose(curDataLim, curFileLim, rtol=1e-2) and fileArray.shape!=baseArray.shape: # We might have to interpolate
#                if len(tData)<len(tFile) or len(xData)<len(xFile): # OK, we can interpolate
#                    data = fileArray[1:,1:]
#                    if len(tFile)==1: # 1d in concentration
#                        f = interpolate.interp1d(xFile, data.flat)#, kind='cubic')
#                        dataNew = f(xData).reshape((1,len(xData)))
#                    elif len(xFile)==1: # 1d in temperature
#                        f = interpolate.interp1d(tFile, data.flat)#, kind='cubic')
#                        dataNew = f(tData).reshape((len(tData),1))
#                    else: # 2d
#                        f = interpolate.interp2d(xFile, tFile, data)#, kind='cubic')
#                        dataNew = f(xData,tData)
#                    fileArray = np.copy(baseArray)
#                    fileArray[1:,1:] = dataNew
#                    fileArray[1:,0] = tData
#                    fileArray[0,1:] = xData
#                else:
#                    raise ValueError("Less points in file, not enough data for interpolation.")
#            elif fileArray.shape==baseArray.shape:
#                pass
#            else:
#                raise ValueError("The array shapes do not match. Check {0}".format(self.getFile(data)))

            if readFromFile or fileArray.shape==baseArray.shape: # Shapes match
                if readFromFile or np.allclose(tData, fileArray[1:,0]): # Temperature data matches
                    if readFromFile or np.allclose(xData, fileArray[0,1:]): # Concentration data matches
                        baseArray = fileArray 
                    else:
                        raise ValueError("Concentration arrays do not match. Check {0} \n and have a look at \n {1} \n vs \n {2} \n yielding \n {3}".format(self.getFile(data),xData,fileArray[0,1:],(xData==fileArray[0,1:])))
                else:
                    raise ValueError("Temperature arrays do not match. Check {0} \n and have a look at \n {1} \n vs \n {2} \n yielding \n {3}".format(self.getFile(data),tData,fileArray[1:,0],(tData==fileArray[1:,0])))
            else:
                raise ValueError("The array shapes do not match. Check {0}".format(self.getFile(data)))
        else: # file not found or update forced
            if func==None and forceUpdate:
                raise ValueError("Need a function to update the data file.")
            for cT,T in enumerate(self.temperature.data):
                for cx,x in enumerate(self.concentration.data):
                    baseArray[cT+1][cx+1] = func(T,x)
            baseArray[0,0] = np.NaN
            baseArray[1:,0] = np.copy(self.temperature.data)
            baseArray[0,1:] = np.copy(self.concentration.data)
            self.writeToFile(data, baseArray)

        return np.copy(baseArray.T[1:].T[1:]) # Remove the first column and row and return


class PureExample(PureData):
    def __init__(self):
        PureData.__init__(self) 
        self.name = "ExamplePure"
        self.description = "Heat transfer fluid TherminolD12 by Solutia"
        self.reference = "Solutia data sheet"
        self.Tmax = 150 + 273.15
        self.Tmin =  50 + 273.15
        self.TminPsat =  self.Tmax
        
        self.temperature.data         = np.array([   50,   60,    70,     80,    90,   100,   110,   120,   130,   140,   150])+273.15 # Kelvin
        self.density.data             = np.array([  740,   733,   726,   717,   710,   702,   695,   687,   679,   670,   662])        # kg/m3
        self.specific_heat.data       = np.array([ 2235,  2280,  2326,  2361,  2406,  2445,  2485,  2528,  2571,  2607,  2645])        # J/kg-K
        self.viscosity.data           = np.array([0.804, 0.704, 0.623, 0.556, 0.498, 0.451, 0.410, 0.374, 0.346, 0.317, 0.289])        # Pa-s
        self.conductivity.data        = np.array([0.105, 0.104, 0.102, 0.100, 0.098, 0.096, 0.095, 0.093, 0.091, 0.089, 0.087])        # W/m-K
        self.saturation_pressure.data = np.array([  0.5,   0.9,   1.4,   2.3,   3.9,   6.0,   8.7,  12.4,  17.6,  24.4,  33.2])        # Pa
        self.reshapeAll()


class SolutionExample(SolutionData):
    def __init__(self):
        SolutionData.__init__(self) 
        self.name = "ExampleSolution"
        self.description = "Ethanol ice slurry"
        self.reference = "SecCool software"
        
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
        
        self.specific_heat.data = np.copy(self.density.data)
        
        self.Tmax = np.max(self.temperature.data)
        self.Tmin = np.min(self.temperature.data)
        self.xmax = np.max(self.concentration.data)
        self.xmin = np.min(self.concentration.data)
        self.xid         = self.ifrac_mass
        self.TminPsat =  self.Tmax
        
#        self.density.data             = np.array([
#          [np.nan,    np.nan,    np.nan,    np.nan,    np.nan,    np.nan,    np.nan],
#          [np.nan,    np.nan,    np.nan,    np.nan,    np.nan,    np.nan,    np.nan],
#          [np.nan,    np.nan,    np.nan,    np.nan,    np.nan,    np.nan,    np.nan],
#          [np.nan,    np.nan,    np.nan,    np.nan,    np.nan,    np.nan,    np.nan],
#          [np.nan,    np.nan,    np.nan,    np.nan,    1016.0,    1008.4,    np.nan],
#          [np.nan,    1033.2,    np.nan,    np.nan,    np.nan,    np.nan,    np.nan],
#          [np.nan,    1025.3,    1018.5,    np.nan,    np.nan,     998.5,     992.0],
#          [np.nan,    np.nan,    1009.2,    np.nan,    np.nan,    np.nan,    np.nan]]) # kg/m3

class DigitalExample(DigitalData):

    def __init__(self):
        DigitalData.__init__(self) 

        self.name = "ExampleDigital"
        self.description = "some fluid"
        self.reference = "none"
        
        self.Tmin = 273.00;
        self.Tmax = 500.00;
        self.xmax = 1.0
        self.xmin = 0.0
        self.xid         = self.ifrac_mass
        self.TminPsat = self.Tmin;
        
        self.temperature.data         = self.getTrange()
        self.concentration.data       = self.getxrange()
        
        def funcRho(T,x):
            return T + x*100.0 + T*(x+0.5)
        self.density.data             = self.getArray(funcRho,"rho")
        
        def funcCp(T,x):
            return T + x*50.0 + T*(x+0.6)
        self.specific_heat.data             = self.getArray(funcCp,"cp")
        

if __name__ == '__main__':
    pass 
##   An example with a pure fluid
#    obj = PureExample()
#    obj.density.type   = obj.density.INCOMPRESSIBLE_POLYNOMIAL
#    obj.density.coeffs = np.zeros((4,6))
#    print(obj.density.coeffs)
#    obj.density.fitCoeffs(obj.temperature.data)
#    print(obj.density.coeffs)
#    obj.saturation_pressure.type   = obj.density.INCOMPRESSIBLE_EXPPOLYNOMIAL
#    obj.saturation_pressure.coeffs = np.zeros((4,6))
#    print(obj.saturation_pressure.coeffs)
#    obj.saturation_pressure.fitCoeffs(obj.temperature.data)
#    print(obj.saturation_pressure.coeffs)
#    print(obj.density.data[2][0],obj.rho(obj.temperature.data[2],10e5,obj.concentration.data[0]))
#    print(obj.density.data[5][0],obj.rho(obj.temperature.data[5],10e5,obj.concentration.data[0]))
#    print(obj.saturation_pressure.data[2][0],obj.psat(obj.temperature.data[2],obj.concentration.data[0]))
#    print(obj.saturation_pressure.data[5][0],obj.psat(obj.temperature.data[5],obj.concentration.data[0]))

##   An example with a solution
#    obj = SolutionExample()
#    obj.density.type   = obj.density.INCOMPRESSIBLE_POLYNOMIAL
#    obj.density.coeffs = np.ones((3,5))
#    print(obj.density.coeffs)
#    obj.density.fitCoeffs(obj.temperature.data,obj.concentration.data)
#    print(obj.density.coeffs)
#    print(obj.density.data[2][0],obj.rho(obj.temperature.data[2],10e5,obj.concentration.data[0]))
#    print(obj.density.data[2][2],obj.rho(obj.temperature.data[2],10e5,obj.concentration.data[2]))

##   An example with a digital fluid
#    obj = DigitalExample()
#    obj.density.type   = obj.density.INCOMPRESSIBLE_POLYNOMIAL
#    obj.density.coeffs = np.ones((3,5))
#    print(obj.density.coeffs)
#    obj.density.fitCoeffs(obj.temperature.data,obj.concentration.data)
#    print(obj.density.coeffs)
#    print(obj.density.data[2][0],obj.rho(obj.temperature.data[2],10e5,obj.concentration.data[0]))
#    print(obj.density.data[2][2],obj.rho(obj.temperature.data[2],10e5,obj.concentration.data[2]))