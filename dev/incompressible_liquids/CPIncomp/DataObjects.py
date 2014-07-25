from __future__ import division, print_function
import numpy as np
import os, CPIncomp, math
from CPIncomp.BaseObjects import IncompressibleData
from abc import ABCMeta

class SolutionData(object):
    """ 
    A base class that defines all the variables needed 
    in order to make a proper fit. You can copy this code
    put in your data and add some documentation for where the
    information came from. 
    """
    __metaclass__ = ABCMeta
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
    __metaclass__ = ABCMeta
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
    __metaclass__ = ABCMeta
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
            if (data!=None and os.path.isfile(self.getFile(data))): # File found
                fileArray = self.getFromFile(data)
                self.temperature.data = np.copy(fileArray[1:,0])
                self.concentration.data = np.copy(fileArray[0,1:])
                readFromFile = True
            else:
                raise ValueError("No temperature and concentration data given and no readable file found for {0}".format(data))
            
        tData = self.round(self.temperature.data)[:,0]
        xData = self.round(self.concentration.data)[:,0]

        baseArray = np.zeros( (len(tData)+1,len(xData)+1) )
        
        if (data!=None and os.path.isfile(self.getFile(data)) and not forceUpdate): # File found and no update wanted
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
            if data!=None: self.writeToFile(data, baseArray)

        return np.copy(baseArray.T[1:].T[1:]) # Remove the first column and row and return


class CoefficientData(SolutionData):
    """ 
    A class to convert parameter arrays from different other sources
    """ 
    __metaclass__ = ABCMeta
    def __init__(self):
        SolutionData.__init__(self) 
        self.reference = "Some other software"
        
    def convertSecCoolArray(self, array):
        if len(array)!=18:
            raise ValueError("The lenght is not equal to 18!")
        
        self.reference = "SecCool software"
        array = np.array(array)        
        tmp = np.zeros((4,6))
        
        tmp[0][0] = array[0]
        tmp[0][1] = array[1]
        tmp[0][2] = array[2]
        tmp[0][3] = array[3]
        tmp[0][4] = array[4]
        tmp[0][5] = array[5]
        
        tmp[1][0] = array[6]
        tmp[1][1] = array[7]
        tmp[1][2] = array[8]
        tmp[1][3] = array[9]
        tmp[1][4] = array[10]
        #tmp[1][5] = array[11]
        
        tmp[2][0] = array[11]
        tmp[2][1] = array[12]
        tmp[2][2] = array[13]
        tmp[2][3] = array[14]
        #tmp[2][4] = array[4]
        #tmp[2][5] = array[5]
        
        tmp[3][0] = array[15]
        tmp[3][1] = array[16]
        tmp[3][2] = array[17]
        #tmp[3][3] = array[3]
        #tmp[3][4] = array[4]
        #tmp[3][5] = array[5]
        
        # Concentration is no longer handled in per cent!
        for i in range(6):
            tmp.T[i] *= 100.0**i 
                
        return tmp
    
    
    def convertSecCoolTfreeze(self, array):
        expo = -1.0
        for i in range(len(array)):
            array[i] = array[i]*np.power(100.0,expo+i)
        array[1] = array[1] + 273.15
        #self.T_freeze.type = self.T_freeze.INCOMPRESSIBLE_POLYOFFSET
        return array



    def convertMelinderArray(self, array):
        """The same function as the SecCool converter, 
        the original source code is slightly different though.
        That is why the implementation is in a transposed form..."""
        
        if len(array)!=18:
            raise ValueError("The lenght is not equal to 18!")
        
        self.reference = "Melinder Book"
        array = np.array(array)        
        tmp = np.zeros((6,4))
        
        tmp[0][0] = array[0] 
        tmp[0][1] = array[6] 
        tmp[0][2] = array[11]
        tmp[0][3] = array[15]
        
        tmp[1][0] = array[1] 
        tmp[1][1] = array[7] 
        tmp[1][2] = array[12]
        tmp[1][3] = array[16]
        
        tmp[2][0] = array[2] 
        tmp[2][1] = array[8] 
        tmp[2][2] = array[13]
        tmp[2][3] = array[17]
        
        tmp[3][0] = array[3] 
        tmp[3][1] = array[9] 
        tmp[3][2] = array[14]
        
        tmp[4][0] = array[4] 
        tmp[4][1] = array[10]
        
        tmp[5][0] = array[5] 
        
        # Concentration is no longer handled in per cent!
        for i in range(6):
            tmp[i] *= 100.0**i 
                
        return tmp.T
    
    def convertMelinderMatrix(self, array):
        """Function to convert the full coefficient array
        from the very first CoolProp implementation
        based on the book by Melinder"""
        if len(array)!=18:
            raise ValueError("The lenght is not equal to 18!")
        if len(array[0])!=5:
            raise ValueError("The lenght is not equal to 5!")
        array = np.array(array)
        tmp = np.zeros((18,5))
        
        for j in range(5):
            tmp[ 0][j] = array[ 0][j]
            tmp[ 1][j] = array[ 4][j]
            tmp[ 2][j] = array[ 8][j]
            tmp[ 3][j] = array[12][j]
            tmp[ 4][j] = array[15][j]
            tmp[ 5][j] = array[17][j]
            tmp[ 6][j] = array[ 1][j]
            tmp[ 7][j] = array[ 5][j]
            tmp[ 8][j] = array[ 9][j]
            tmp[ 9][j] = array[13][j]
            tmp[10][j] = array[16][j]
            tmp[11][j] = array[ 2][j]
            tmp[12][j] = array[ 6][j]
            tmp[13][j] = array[10][j]
            tmp[14][j] = array[14][j]
            tmp[15][j] = array[ 3][j]
            tmp[16][j] = array[ 7][j]
            tmp[17][j] = array[11][j]
            
        return tmp
    
    
    
    def setMelinderMatrix(self, matrix):
#        matrix = np.array([
#        [-26.29           , 958.1           ,3887           ,   0.4175            ,   1.153           ],
#        [ -0.000002575    ,  -0.4151        ,   7.201       ,   0.0007271         ,  -0.03866         ],
#        [ -0.000006732    ,  -0.002261      ,  -0.08979     ,   0.0000002823      ,   0.0002779       ],
#        [  0.000000163    ,   0.0000002998  ,  -0.000439    ,   0.000000009718    ,  -0.000001543     ],
#        [ -1.187          ,  -1.391         , -18.5         ,  -0.004421          ,   0.005448        ],
#        [ -0.00001609     ,  -0.0151        ,   0.2984      ,  -0.00002952        ,   0.0001008       ],
#        [  0.000000342    ,   0.0001113     ,  -0.001865    ,   0.00000007336     ,  -0.000002809     ],
#        [  0.0000000005687,  -0.0000003264  ,  -0.00001718  ,   0.0000000004328   ,   0.000000009811  ],
#        [ -0.01218        ,  -0.01105       ,  -0.03769     ,   0.00002044        ,  -0.0005552       ],
#        [  0.0000003865   ,   0.0001828     ,  -0.01196     ,   0.0000003413      ,   0.000008384     ],
#        [  0.000000008768 ,  -0.000001641   ,   0.00009801  ,  -0.000000003665    ,  -0.00000003997   ],
#        [ -0.0000000002095,   0.0000000151  ,   0.000000666 ,  -0.00000000002791  ,  -0.0000000003466 ],
#        [ -0.00006823     ,  -0.0001208     ,  -0.003776    ,   0.0000002943      ,   0.000003038     ],
#        [  0.00000002137  ,   0.000002992   ,  -0.00005611  ,  -0.0000000009646   ,  -0.00000007435   ],
#        [ -0.0000000004271,   0.000000001455,  -0.0000007811,   0.00000000003174  ,   0.0000000007442 ],
#        [  0.0000001297   ,   0.000004927   ,  -0.0001504   ,  -0.0000000008666   ,   0.00000006669   ],
#        [ -0.0000000005407,  -0.0000001325  ,   0.000007373 ,  -0.0000000000004573,  -0.0000000009105 ],
#        [  0.00000002363  ,  -0.00000007727 ,   0.000006433 ,  -0.0000000002033   ,  -0.0000000008472 ]
#        ])
        
        coeffs = self.convertMelinderMatrix(matrix).T
        
        self.T_freeze.type = self.T_freeze.INCOMPRESSIBLE_POLYNOMIAL
        self.T_freeze.coeffs = self.convertMelinderArray(coeffs[0])
        
        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = self.convertMelinderArray(coeffs[1])
        
        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = self.convertMelinderArray(coeffs[2])
        
        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = self.convertMelinderArray(coeffs[3])
        
        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_POLYNOMIAL
        self.viscosity.coeffs = self.convertMelinderArray(coeffs[4])
    
    
    
    
    
