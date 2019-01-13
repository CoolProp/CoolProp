from __future__ import division, print_function
import numpy as np
import os, math
from .BaseObjects import IncompressibleData, IncompressibleFitter
from abc import ABCMeta


class SolutionData(object):
    """
    A base class that defines all the variables needed
    in order to make a proper fit. You can copy this code
    put in your data and add some documentation for where the
    information came from.
    """
    ifrac_mass = "mass"
    ifrac_mole = "mole"
    ifrac_volume = "volume"
    ifrac_undefined = "not defined"
    ifrac_pure = "pure"

    def __init__(self):

        self.significantDigits = 7

        self.name = None  # Name of the current fluid
        self.description = None  # Description of the current fluid
        self.reference = None  # Reference data for the current fluid

        self.Tmax = None  # Maximum temperature in K
        self.Tmin = None  # Minimum temperature in K
        self.xmax = 1.0  # Maximum concentration
        self.xmin = 0.0  # Minimum concentration
        self.xid = self.ifrac_undefined  # Concentration is mole, mass or volume-based
        self.TminPsat = None  # Minimum saturation temperature in K
        self.Tbase = None  # Base value for temperature fits
        self.xbase = None  # Base value for concentration fits

        self.temperature = IncompressibleData()  # Temperature for data points in K
        self.concentration = IncompressibleData()  # Concentration data points in weight fraction
        self.density = IncompressibleData()  # Density in kg/m3
        self.specific_heat = IncompressibleData()  # Heat capacity in J/(kg.K)
        self.viscosity = IncompressibleData()  # Dynamic viscosity in Pa.s
        self.conductivity = IncompressibleData()  # Thermal conductivity in W/(m.K)
        self.saturation_pressure = IncompressibleData()  # Saturation pressure in Pa
        self.T_freeze = IncompressibleData()  # Freezing temperature in K
        self.mass2input = IncompressibleData()  # dd
        self.volume2input = IncompressibleData()  # dd
        self.mole2input = IncompressibleData()  # dd

        # Some of the functions might need a guess array
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

    def roundSingle(self, x):
        if x == 0.0: return 0.0
        return round(x, self.significantDigits - int(math.floor(math.log10(abs(x)))) - 1)

    def round(self, x):
        r, c, res = IncompressibleFitter.shapeArray(x)
        #digits = -1*np.floor(np.log10(res))+self.significantDigits-1
        for i in range(r):
            for j in range(c):
                if np.isfinite(res[i, j]):
                    res[i, j] = self.roundSingle(res[i, j])
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

    def checkT(self, T, p, x):
        if self.Tmin <= 0.: raise ValueError("Please specify the minimum temperature.")
        if self.Tmax <= 0.: raise ValueError("Please specify the maximum temperature.");
        if ((self.Tmin > T) or (T > self.Tmax)): raise ValueError("Your temperature {0} is not between {1} and {2}.".format(T, self.Tmin, self.Tmax))
        TF = 0.0
        if (self.T_freeze.type != IncompressibleData.INCOMPRESSIBLE_NOT_SET): TF = self.Tfreeze(T, p, x)
        if (T < TF): raise ValueError("Your temperature {0} is below the freezing point of {1}.".format(T, TF))
        else: return True
        return False

    def checkP(self, T, p, x):
        ps = 0.0
        if (self.saturation_pressure.type != IncompressibleData.INCOMPRESSIBLE_NOT_SET): ps = self.psat(T, p, x)
        if (p < 0.0): raise ValueError("You cannot use negative pressures: {0} < {1}. ".format(p, 0.0))
        if (p < ps): raise ValueError("Equations are valid for liquid phase only: {0} < {1}. ".format(p, ps))
        else: return True
        return False

    def checkX(self, x):
        if (self.xmin < 0.0 or self.xmin > 1.0): raise ValueError("Please specify the minimum concentration between 0 and 1.");
        if (self.xmax < 0.0 or self.xmax > 1.0): raise ValueError("Please specify the maximum concentration between 0 and 1.");
        if ((self.xmin > x) or (x > self.xmax)): raise ValueError("Your composition {0} is not between {1} and {2}.".format(x, self.xmin, self.xmax))
        else: return True
        return False

    def checkTPX(self, T, p, x):
        try:
            return (self.checkT(T, p, x) and self.checkP(T, p, x) and self.checkX(x))
        except ValueError as ve:
            #print("Check failed: {0}".format(ve))
            pass
        return False

    def rho(self, T, p=0.0, x=0.0, c=None):
        if not self.checkTPX(T, p, x): return np.NAN
        if c is None:
            c = self.density.coeffs
        if self.density.type == self.density.INCOMPRESSIBLE_POLYNOMIAL:
            return np.polynomial.polynomial.polyval2d(T - self.Tbase, x - self.xbase, c)
        else: raise ValueError("Unknown function.")

    def c(self, T, p=0.0, x=0.0, c=None):
        if not self.checkTPX(T, p, x): return np.NAN
        if c is None:
            c = self.specific_heat.coeffs
        if self.specific_heat.type == self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL:
            return np.polynomial.polynomial.polyval2d(T - self.Tbase, x - self.xbase, c)
        else: raise ValueError("Unknown function.")

    def cp(self, T, p=0.0, x=0.0, c=None):
        return self.c(T, p, x, c)

    def cv(self, T, p=0.0, x=0.0, c=None):
        return self.c(T, p, x, c)

    def u(self, T, p=0.0, x=0.0, c=None):
        if not self.checkTPX(T, p, x): return np.NAN
        if c is None:
            c = self.specific_heat.coeffs
        if self.specific_heat.type == self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL:
            c_tmp = np.polynomial.polynomial.polyint(c)
            return np.polynomial.polynomial.polyval2d(T - self.Tbase, x - self.xbase, c_tmp)
        else: raise ValueError("Unknown function.")
    # def u   (T, p, x);

    def h(self, T, p=0.0, x=0.0):
        return self.h_u(T, p, x)

    def visc(self, T, p=0.0, x=0.0, c=None):
        if not self.checkTPX(T, p, x): return np.NAN
        return self.viscosity.baseFunction(T, x, self.Tbase, self.xbase, c=c)

    def cond(self, T, p=0.0, x=0.0, c=None):
        if not self.checkTPX(T, p, x): return np.NAN
        return self.conductivity.baseFunction(T, x, self.Tbase, self.xbase, c=c)

    def psat(self, T, p=0.0, x=0.0, c=None):
        if (T <= self.TminPsat): return 0.0
        return self.saturation_pressure.baseFunction(T, x, self.Tbase, self.xbase, c=c)

    def Tfreeze(self, T, p=0.0, x=0.0, c=None):
        if c is None:
            c = self.T_freeze.coeffs

        if self.T_freeze.type == self.T_freeze.INCOMPRESSIBLE_POLYNOMIAL:
            return np.polynomial.polynomial.polyval2d(p - 0.0, x - self.xbase, c)

        elif self.T_freeze.type == self.T_freeze.INCOMPRESSIBLE_POLYOFFSET:
            #if y!=0.0: raise ValueError("This is 1D only, use x not y.")
            return self.T_freeze.basePolyOffset(c, x)  # offset included in coeffs

        elif self.T_freeze.type == self.T_freeze.INCOMPRESSIBLE_EXPONENTIAL:
            #if y!=0.0: raise ValueError("This is 1D only, use x not y.")
            return self.T_freeze.baseExponential(c, x)

        elif self.T_freeze.type == IncompressibleData.INCOMPRESSIBLE_LOGEXPONENTIAL:
            #if y!=0.0: raise ValueError("This is 1D only, use x not y.")
            return self.T_freeze.baseLogexponential(c, x)

        elif self.T_freeze.type == self.T_freeze.INCOMPRESSIBLE_EXPPOLYNOMIAL:
            return np.exp(np.polynomial.polynomial.polyval2d(p - 0.0, x - self.xbase, c))

        else:
            raise ValueError("Unknown function: {0}.".format(self.T_freeze.type))

    # def V2M (T,           y);
    # def M2M (T,           x);

    def h_u(self, T, p, x):
        return self.u(T, p, x) + p / self.rho(T, p, x) - self.href

    def u_h(self, T, p, x):
        return self.h(T, p, x) - p / self.rho(T, p, x) + self.href

    def set_reference_state(self, T0, p0, x0=0.0, h0=0.0, s0=0.0):
        self.rhoref = self.rho(T0, p0, x0)
        self.pref = p0
        self.uref = h0 - p0 / self.rhoref
        self.uref = self.u(T0, p0, x0)
        self.href = h0
        self.sref = s0
        self.href = self.h(T0, p0, x0)
        self.sref = self.s(T0, p0, x0)


class PureData(SolutionData):
    """
    An extension of the solution data that makes it
    easier to gather data for pure fluids.
    """
    __metaclass__ = ABCMeta

    def __init__(self):
        SolutionData.__init__(self)
        self.xbase = 0.0
        self.xid = self.ifrac_pure
        self.concentration.data = np.array([0])  # mass fraction

    def reshapeData(self, dataArray, length):
        """
        Makes any array 1-dimensional and implicitly
        checks the length.
        """
        if not dataArray is None:
            return np.reshape(dataArray, (length, 1))
        else:
            return None

    def reshapeAll(self):
        len_T = len(self.temperature.data)
        #len_x = len(self.concentration.data)
        self.density.data = self.reshapeData(self.density.data, len_T)
        self.specific_heat.data = self.reshapeData(self.specific_heat.data, len_T)
        self.viscosity.data = self.reshapeData(self.viscosity.data, len_T)
        self.conductivity.data = self.reshapeData(self.conductivity.data, len_T)
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
        return os.path.join(os.path.dirname(__file__), 'data', self.name + "_" + data + ".txt")

    def getFromFile(self, data):
        fullPath = self.getFile(data)
        _, _, res = IncompressibleFitter.shapeArray(np.loadtxt(fullPath))
        return res

    def writeToFile(self, data, array):
        fullPath = self.getFile(data)
        if not os.path.exists(os.path.dirname(fullPath)):
            os.makedirs(os.path.dirname(fullPath))
        stdFmt = "%1.{0}e".format(int(self.significantDigits - 1))
        return np.savetxt(fullPath, array, fmt=stdFmt)

    def getTrange(self):
        if self.Tmin < self.Tmax:
            return np.linspace(self.Tmin, self.Tmax, 20)
        else:
            return np.array([0.0])

    def getxrange(self):
        if self.xmin < self.xmax and self.xid != self.ifrac_pure:
            return np.linspace(self.xmin, self.xmax, 20)
        else:
            return np.array([0.0])

    def getArray(self, dataID=None, func=None, x_in=None, y_in=None, DEBUG=False):
        """ Tries to read a data file, overwrites it if x or y do not match

        :param dataID  : ID to construct the path to the data file
        :param func    : Callable object that can take x_in and y_in
        :param x_in    : a 1D array in x direction or 2D with one column, most likely temperature
        :param y_in    : a 1D array in y direction or 2D with one row, most likely cocentration
        :param DEBUG   : a boolean that controls verbosity

        :returns       : Returns a tuple with three entries: x(1D),y(1D),data(2D)
        """

        x = None
        y = None
        z = None

        # First we try to read the file
        if (not dataID is None and os.path.isfile(self.getFile(dataID))):  # File found
            fileArray = self.getFromFile(dataID)
            x = np.copy(fileArray[1:, 0])
            y = np.copy(fileArray[0, 1:])
            z = np.copy(fileArray[1:, 1:])
        else:
            if DEBUG: print("No readable file found for {0}: {1}".format(dataID, self.getFile(dataID)))

        updateFile = DEBUG

        if not x_in is None:  # Might need update
            if not x is None:  # Both given, check if different
                mask = np.isfinite(x)
                if IncompressibleFitter.allClose(x[mask], x_in[mask]):
                    if DEBUG: print("Both x-arrays are the same, no action required.")
                    updateFile = (updateFile or False)  # Do not change a True value to False
                else:
                    updateFile = True
                    if DEBUG: print("x-arrays do not match. {0} contains \n {1} \n and will be updated with \n {2}".format(self.getFile(dataID), x, x_in))
            else: updateFile = True
        elif x is None:
            if DEBUG: print("Could not load x from file {0} and no x_in provided, aborting.".format(self.getFile(dataID)))
            return None, None, None
        else: updateFile = (updateFile or False)  # Do not change a True value to False

        if not y_in is None:  # Might need update
            if not y is None:  # Both given, check if different
                mask = np.isfinite(y)
                if IncompressibleFitter.allClose(y[mask], y_in[mask]):
                    if DEBUG: print("Both y-arrays are the same, no action required.")
                    updateFile = (updateFile or False)  # Do not change a True value to False
                else:
                    updateFile = True
                    if DEBUG: print("y-arrays do not match. {0} contains \n {1} \n and will be updated with \n {2}".format(self.getFile(dataID), y, y_in))
            else: updateFile = True
        elif y is None:
            if DEBUG: print("Could not load y from file {0} and no y_in provided, aborting.".format(self.getFile(dataID)))
            return None, None, None
        else: updateFile = (updateFile or False)  # Do not change a True value to False

        if DEBUG: print("Updating data file {0}".format(updateFile))

        if not updateFile: return x, y, z  # Done, data read from file

        # Overwrite inputs
        x = x_in
        y = y_in
        z = np.zeros((len(x) + 1, len(y) + 1))
        r, c = z.shape

        if func is None: raise ValueError("Need a function to update the data file.")

        for i in range(r - 1):
            for j in range(c - 1):
                z[i + 1, j + 1] = func(x[i], y[j])
        z[0, 0] = np.NaN
        z[1:, 0] = x
        z[0, 1:] = y

        if not dataID is None:
            self.writeToFile(dataID, z)
        else:
            if DEBUG: print("Not updating data file, dataID is missing.")

        return x, y, z[1:, 1:]


class CoefficientData(SolutionData):
    """
    A class to convert parameter arrays from different other sources
    """
    __metaclass__ = ABCMeta

    def __init__(self):
        SolutionData.__init__(self)
        self.reference = "Some other software"

    def convertSecCoolArray(self, array):
        if len(array) != 18:
            raise ValueError("The length is not equal to 18!")

        self.reference = "SecCool software"
        array = np.array(array)
        tmp = np.zeros((4, 6))

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
            array[i] = array[i] * np.power(100.0, expo + i)
        array[1] = array[1] + 273.15
        #self.T_freeze.type = self.T_freeze.INCOMPRESSIBLE_POLYOFFSET
        return array

    def convertMelinderArray(self, array):
        """The same function as the SecCool converter,
        the original source code is slightly different though.
        That is why the implementation is in a transposed form..."""

        if len(array) != 18:
            raise ValueError("The length is not equal to 18!")

        #self.reference = "Melinder Book"
        array = np.array(array)
        tmp = np.zeros((6, 4))

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
        if len(array) != 18:
            raise ValueError("The length is not equal to 18!")
        if len(array[0]) != 5:
            raise ValueError("The length is not equal to 5!")
        array = np.array(array)
        tmp = np.zeros((18, 5))

        for j in range(5):
            tmp[0][j] = array[0][j]
            tmp[1][j] = array[4][j]
            tmp[2][j] = array[8][j]
            tmp[3][j] = array[12][j]
            tmp[4][j] = array[15][j]
            tmp[5][j] = array[17][j]
            tmp[6][j] = array[1][j]
            tmp[7][j] = array[5][j]
            tmp[8][j] = array[9][j]
            tmp[9][j] = array[13][j]
            tmp[10][j] = array[16][j]
            tmp[11][j] = array[2][j]
            tmp[12][j] = array[6][j]
            tmp[13][j] = array[10][j]
            tmp[14][j] = array[14][j]
            tmp[15][j] = array[3][j]
            tmp[16][j] = array[7][j]
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

        self.T_freeze.source = self.T_freeze.SOURCE_COEFFS
        self.T_freeze.type = self.T_freeze.INCOMPRESSIBLE_POLYNOMIAL
        self.T_freeze.coeffs = self.convertMelinderArray(coeffs[0])
        self.T_freeze.coeffs[0, 0] += 273.15
        self.T_freeze.coeffs = np.array([self.T_freeze.coeffs[0]])
        # print(self.T_freeze.coeffs)

        self.density.source = self.density.SOURCE_COEFFS
        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = self.convertMelinderArray(coeffs[1])

        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = self.convertMelinderArray(coeffs[2])

        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = self.convertMelinderArray(coeffs[3])

        self.viscosity.source = self.viscosity.SOURCE_COEFFS
        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
        self.viscosity.coeffs = self.convertMelinderArray(coeffs[4])
        self.viscosity.coeffs[0, 0] -= math.log(1000)  # Fixes the units mPa s -> Pa s
