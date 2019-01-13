from __future__ import division, print_function
import numpy as np
from .BaseObjects import IncompressibleData, IncompressibleFitter
from .DataObjects import DigitalData, PureData
import os


class SecCoolSolutionData(DigitalData):
    """
    A base class that can be fed with a fluid ID from SecCool
    to read data files sitting in data/SecCool/xMass and
    data/SecCool/xVolume.
    """

    def __init__(self, sFile=None, sFolder=None, name=None, desc=None, ref='SecCool software', densityFactor=None, heatFactor=None, conductivityFactor=None, viscosityFactor=None):
        DigitalData.__init__(self)
        self.allowNegativeData = False

        if sFile == None: raise ValueError("File name cannot be empty.")
        if sFolder == None: raise ValueError("Folder name cannot be empty.")
        if name == None: raise ValueError("Fluid name cannot be empty.")
        if desc == None: raise ValueError("Description cannot be empty.")

        self.sFile = sFile
        self.sFolder = sFolder

        self.name = name
        self.description = desc
        self.reference = ref

        self.temperature.data, self.concentration.data, self.density.data = self.getArray(dataID='Rho')

        self.Tmax = np.max(self.temperature.data)
        self.Tmin = np.min(self.temperature.data)
        self.xmax = np.max(self.concentration.data)
        self.xmin = np.min(self.concentration.data)
        if self.sFolder == 'xVolume':
            self.xid = self.ifrac_volume
        elif self.sFolder == 'xMass':
            self.xid = self.ifrac_mass
        elif self.sFolder == 'xPure':
            self.xid = self.ifrac_pure
            self.xbase = 0.0  # Disables the reset of xmax and xmin
            self.xmax = 1.0
            self.xmin = 0.0
        else:
            raise ValueError("Unknown folder type specified.")
        self.TminPsat = self.Tmax

        try:
            self.density.xData, self.density.yData, self.density.data = self.getArray(dataID="Rho")
            # while np.max(self.density.data[np.isfinite(self.density.data)])<500: # Expect values around 1e3
            #    self.density.data *= 10.0
            if densityFactor is not None: self.density.data *= densityFactor
            self.density.source = self.density.SOURCE_DATA
        except:
            if self.density.DEBUG: print("Could not load {}".format(self.getFile("Rho")))
            pass

        try:
            self.specific_heat.xData, self.specific_heat.yData, self.specific_heat.data = self.getArray(dataID='Cp')
            # while np.max(self.specific_heat.data[np.isfinite(self.specific_heat.data)])<1000: # Expect values around 2e3
            #    self.specific_heat.data *= 10.0
            if heatFactor is not None: self.specific_heat.data *= heatFactor
            self.specific_heat.source = self.specific_heat.SOURCE_DATA
        except:
            if self.specific_heat.DEBUG: print("Could not load {}".format(self.getFile("Cp")))
            pass

        try:
            self.conductivity.xData, self.conductivity.yData, self.conductivity.data = self.getArray(dataID='Cond')
            # while np.max(self.conductivity.data[np.isfinite(self.conductivity.data)])>2: # Expect value below 1
            #    self.conductivity.data *=  0.1
            if conductivityFactor is not None: self.conductivity.data *= conductivityFactor
            self.conductivity.source = self.conductivity.SOURCE_DATA
        except:
            if self.conductivity.DEBUG: print("Could not load {}".format(self.getFile("Cond")))
            pass

        try:
            self.viscosity.xData, self.viscosity.yData, self.viscosity.data = self.getArray(dataID='Mu')
            # while np.max(self.viscosity.data[np.isfinite(self.viscosity.data)])>0.2: # Expect value below 0.1
            #    self.viscosity.data *=  0.1
            if viscosityFactor is not None: self.viscosity.data *= viscosityFactor
            self.viscosity.source = self.viscosity.SOURCE_DATA
        except:
            if self.viscosity.DEBUG: print("Could not load {}".format(self.getFile("Mu")))
            pass


#    def interp(self, tOld, xOld, dataNew):
#        tCur = self.temperature.data
#        xCur = self.concentration.data
#
#        data = None
#        if len(tCur)==1: # 1d in concentration
#            f = interpolate.interp1d(xCur, dataNew.flat)#, kind='cubic')
#            data = f(xOld).reshape((1,len(xOld)))
#        elif len(xCur)==1: # 1d in temperature
#            f = interpolate.interp1d(tCur, dataNew.flat)#, kind='cubic')
#            data = f(tOld).reshape((len(tOld),1))
#        else: # 2d
#            f = interpolate.interp2d(xCur, tCur, dataNew)#, kind='cubic')
#            data = f(xOld,tOld)
#        return data

#    def interp(self, tOld, xOld, dataNew):
#        tCur = self.temperature.data
#        xCur = self.concentration.data
#
#        data = None
#        if len(tCur)==1: # 1d in concentration
#            f = interpolate.InterpolatedUnivariateSpline(xCur, dataNew.flat, k=1)
#            data = f(xOld).reshape((1,len(xOld)))
#        elif len(xCur)==1: # 1d in temperature
#            f = interpolate.InterpolatedUnivariateSpline(tCur, dataNew.flat, k=1)
#            data = f(tOld).reshape((len(tOld),1))
#        else: # 2d
#            #f = interpolate.SmoothBivariateSpline(xCur, tCur, dataNew, kx=1, ky=1)
#            f = interpolate.interp2d(xCur, tCur, dataNew)#, kind='cubic')
#            data = f(xOld,tOld)
#        return data


    def fitFluid(self):

        if self.Tbase == None:
            self.Tbase = (self.Tmin + self.Tmax) / 2.0
        if self.xbase == None:
            self.xbase = (self.xmin + self.xmax) / 2.0

        std_coeffs = np.zeros((4, 6))
        errList = (ValueError, AttributeError, TypeError, RuntimeError)

        try:
            self.density.coeffs = np.copy(std_coeffs)
            self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
            self.density.fitCoeffs(self.Tbase, self.xbase)
        except errList as ve:
            if self.density.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(self.name, 'density', ve))
            pass

        try:
            self.specific_heat.coeffs = np.copy(std_coeffs)
            self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
            self.specific_heat.fitCoeffs(self.Tbase, self.xbase)
        except errList as ve:
            if self.specific_heat.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(self.name, 'specific heat', ve))
            pass

        try:
            self.conductivity.coeffs = np.copy(std_coeffs)
            self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
            self.conductivity.fitCoeffs(self.Tbase, self.xbase)
        except errList as ve:
            if self.conductivity.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(self.name, 'conductivity', ve))
            pass

#        try:
#            self.viscosity.coeffs = np.copy(std_coeffs)
#            self.viscosity.type   = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
#            self.viscosity.fitCoeffs(self.Tbase,self.xbase)
#        except errList as ve:
#            if self.viscosity.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(self.name,'viscosity',ve))
#            pass

        try:
            tried = False
            if len(self.viscosity.yData) == 1:  # and np.isfinite(fluidObject.viscosity.data).sum()<10:
                self.viscosity.coeffs = np.array([+7e+2, -6e+1, +1e+1])
                self.viscosity.type = IncompressibleData.INCOMPRESSIBLE_EXPONENTIAL
                self.viscosity.fitCoeffs(self.Tbase, self.xbase)
                if self.viscosity.coeffs is None or IncompressibleFitter.allClose(self.viscosity.coeffs, np.array([+7e+2, -6e+1, +1e+1])):  # Fit failed
                    tried = True
            if len(self.viscosity.yData) > 1 or tried:
                self.viscosity.coeffs = np.copy(std_coeffs)  # np.zeros(np.round(np.array(std_coeffs.shape) * 1.5))
                self.viscosity.type = IncompressibleData.INCOMPRESSIBLE_EXPPOLYNOMIAL
                self.viscosity.fitCoeffs(self.Tbase, self.xbase)
        except errList as ve:
            if self.viscosity.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(self.name, 'viscosity', ve))
            pass

        # reset data for getArray and read special files
        if self.xid != self.ifrac_pure and self.xid != self.ifrac_undefined:
            allowNegativeData_org = self.allowNegativeData
            self.allowNegativeData = True

            try:
                x, _, z = self.getArray(dataID='TFreeze')
                self.T_freeze.yData = (x - 273.15) / 100.0
                self.T_freeze.xData = [0.0]
                if np.min(z) < 150: z += 273.15
                self.T_freeze.data = z.T
                try:
                    self.T_freeze.source = self.T_freeze.SOURCE_DATA
                    #self.T_freeze.type   = self.T_freeze.INCOMPRESSIBLE_EXPONENTIAL
                    #self.T_freeze.coeffs = np.array([+7e+6, +6e+4, +1e+1])
                    self.T_freeze.coeffs = np.zeros(std_coeffs.shape)
                    self.T_freeze.type = self.T_freeze.INCOMPRESSIBLE_EXPPOLYNOMIAL
                    # if np.isfinite(self.T_freeze.data).sum()<10:
                    #    self.T_freeze.coeffs = np.array([+7e+6, +6e+4, +1e+1])
                    #    self.T_freeze.type   = self.T_freeze.INCOMPRESSIBLE_EXPONENTIAL
                    # else:
                    #    self.T_freeze.coeffs = np.zeros(np.round(np.array(std_coeffs.shape) * 2))
                    #    self.T_freeze.type   = self.T_freeze.INCOMPRESSIBLE_EXPPOLYNOMIAL
                    self.T_freeze.fitCoeffs(self.Tbase, self.xbase)
                except errList as ve:
                    if self.T_freeze.DEBUG: print("{0}: Could not fit {1} coefficients: {2}".format(self.name, "T_freeze", ve))
                    pass
            except errList as ve:
                if self.T_freeze.DEBUG: print("{0}: Could not load {1} data: {2}".format(self.name, "T_freeze", ve))
                pass

            # reset data for getArray again
            self.allowNegativeData = allowNegativeData_org
            try:
                x, _, z = self.getArray(dataID='Vol2Mass')
                massData = z[:, 0] / 100.0
                volData = (x - 273.15) / 100.0

                if self.xid == self.ifrac_volume:
                    _, _, self.mass2input.data = IncompressibleFitter.shapeArray(volData, axs=1)
                    self.mass2input.xData = [0.0]
                    self.mass2input.yData = massData
                    try:
                        self.mass2input.coeffs = np.copy(std_coeffs)
                        self.mass2input.source = self.mass2input.SOURCE_DATA
                        self.mass2input.type = self.mass2input.INCOMPRESSIBLE_POLYNOMIAL
                        self.mass2input.fitCoeffs(self.Tbase, self.xbase)
                    except errList as ve:
                        if self.mass2input.DEBUG: print("{0}: Could not fit {1} coefficients: {2}".format(self.name, "mass2input", ve))
                        pass
                elif self.xid == self.ifrac_mass:
                    _, _, self.volume2input.data = IncompressibleFitter.shapeArray(massData, axs=1)
                    self.volume2input.xData = [0.0]
                    self.volume2input.yData = volData
                    try:
                        self.volume2input.coeffs = np.copy(std_coeffs)
                        self.volume2input.source = self.volume2input.SOURCE_DATA
                        self.volume2input.type = self.volume2input.INCOMPRESSIBLE_POLYNOMIAL
                        self.volume2input.fitCoeffs(self.Tbase, self.xbase)
                    except errList as ve:
                        if self.volume2input.DEBUG: print("{0}: Could not fit {1} coefficients: {2}".format(self.name, "volume2input", ve))
                        pass
                else:
                    raise ValueError("Unknown xid specified.")
            except errList as ve:
                if self.mass2input.DEBUG or self.volume2input.DEBUG: print("{0}: Could not load {1} data: {2}".format(self.name, "Vol2Mass", ve))
                pass

    # Redefine some functions to avoid data loss
    def getFile(self, data):
        return os.path.join(os.path.dirname(__file__), 'data', 'SecCool', self.sFolder, self.sFile + "_" + data + ".txt")

    def getFromFile(self, data):
        fullPath = self.getFile(data)
        # TODO: improve this case
        if not os.path.isfile(fullPath):
            arr = np.empty((1, 1))
            arr[:] = np.NAN
            return arr
        r, c, res = IncompressibleFitter.shapeArray(np.loadtxt(fullPath, dtype=type('string')))
#        numbers = res.astype(np.float)
        numbers = np.zeros((r, c))
        for i in range(r):
            for j in range(c):
                nu = np.NAN
                try:
                    nu = np.float(res[i, j])
                    if i == 0: nu *= 1e-2  # Percent to fraction
                    if j == 0: nu += 273.15  # Celsius to Kelvin
                    if not self.allowNegativeData and nu < 0:
                        nu = np.NAN  # invalid entries
                except (ValueError, TypeError) as ve:
                    if False: print("Could not convert entry: {0}".format(ve))
                    if i == 0: nu = 0.0  # Dummy for tables without concentration (TFreeze and Vol2Mass)
                    pass
                numbers[i, j] = nu
        if numbers[1, 0] > numbers[-1, 0]:
            numbers[1:, :] = numbers[1:, :][::-1, :]
        if numbers[0, 1] > numbers[0, -1]:
            numbers[:, 1:] = numbers[:, 1:][:, ::-1]
        return numbers

    def writeToFile(self, data, array):
        raise ValueError("You should not overwrite the SecCool data.")

    @staticmethod
    def factory():
        """
        A static method that returns a list of SecCool
        solutions. Coefficients are fitted by calling
        fitFluid and then the objects should be ready
        to be written to JSON.
        """
        print("Loading SecCool fluids: ", end="")
        sec = []
        sec += [SecCoolSolutionData(sFile='Antifrogen KF', sFolder='xVolume', name='AKF', desc='Antifrogen KF, Potassium Formate', ref='Clariant2000,Skovrup2013', densityFactor=None, heatFactor=1e3, conductivityFactor=None, viscosityFactor=1e-3)]
        print("{0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Antifrogen L', sFolder='xVolume', name='AL', desc='Antifrogen L, Propylene Glycol', ref='Clariant2000,Skovrup2013', densityFactor=None, heatFactor=1e3, conductivityFactor=None, viscosityFactor=1e-3)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Antifrogen N', sFolder='xVolume', name='AN', desc='Antifrogen N, Ethylene Glycol', ref='Clariant2000,Skovrup2013', densityFactor=None, heatFactor=1e3, conductivityFactor=None, viscosityFactor=1e-3)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='ASHRAE, Ethylene Glycol', sFolder='xVolume', name='AEG', desc='ASHRAE, Ethylene Glycol', ref='ASHRAE2001,Skovrup2013', densityFactor=None, heatFactor=None, conductivityFactor=1e-3, viscosityFactor=1e-5)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='ASHRAE, Propylene Glycol', sFolder='xVolume', name='APG', desc='ASHRAE, Propylene Glycol', ref='ASHRAE2001,Skovrup2013', densityFactor=None, heatFactor=None, conductivityFactor=1e-3, viscosityFactor=1e-5)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Glykosol N', sFolder='xVolume', name='GKN', desc='Glykosol N, Ethylene Glycol', ref='PKS2005,Skovrup2013', densityFactor=1e3, heatFactor=1e3, conductivityFactor=None, viscosityFactor=1e-3)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Pekasol 2000', sFolder='xVolume', name='PK2', desc='Pekasol 2000, K acetate/formate', ref='PKS2005,Skovrup2013', densityFactor=1e3, heatFactor=1e3, conductivityFactor=None, viscosityFactor=1e-3)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Pekasol L', sFolder='xVolume', name='PKL', desc='Pekasol L, Propylene Glycol', ref='PKS2005,Skovrup2013', densityFactor=1e3, heatFactor=1e3, conductivityFactor=None, viscosityFactor=1e-3)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Zitrec AC', sFolder='xVolume', name='ZAC', desc='Zitrec AC, Corrosion Inhibitor', ref='Arteco2010,Skovrup2013', densityFactor=None, heatFactor=1e3, conductivityFactor=None, viscosityFactor=None)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Zitrec FC', sFolder='xVolume', name='ZFC', desc='Zitrec FC, Propylene Glycol', ref='Arteco2010,Skovrup2013', densityFactor=None, heatFactor=1e3, conductivityFactor=None, viscosityFactor=None)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Zitrec LC', sFolder='xVolume', name='ZLC', desc='Zitrec LC, Propylene Glycol', ref='Arteco2010,Skovrup2013', densityFactor=None, heatFactor=1e3, conductivityFactor=None, viscosityFactor=None)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Zitrec MC', sFolder='xVolume', name='ZMC', desc='Zitrec MC, Ethylene Glycol', ref='Arteco2010,Skovrup2013', densityFactor=None, heatFactor=1e3, conductivityFactor=None, viscosityFactor=None)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Zitrec M', sFolder='xVolume', name='ZM', desc='Zitrec M, Ethylene Glycol', ref='Arteco2010,Skovrup2013', densityFactor=None, heatFactor=1e3, conductivityFactor=None, viscosityFactor=None)]
        print(", {0}".format(sec[-1].name), end="")

        sec += [SecCoolSolutionData(sFile='Melinder, Ammonia', sFolder='xMass', name='MAM2', desc='Melinder, Ammonia', ref='Melinder2010,Skovrup2013', densityFactor=None, heatFactor=None, conductivityFactor=None, viscosityFactor=1e-5)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Melinder, Calcium Cloride', sFolder='xMass', name='MCA2', desc='Melinder, Calcium Chloride', ref='Melinder2010,Skovrup2013', densityFactor=None, heatFactor=None, conductivityFactor=None, viscosityFactor=1e-5)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Melinder, Ethanol', sFolder='xMass', name='MEA2', desc='Melinder, Ethanol', ref='Melinder2010,Skovrup2013', densityFactor=None, heatFactor=None, conductivityFactor=None, viscosityFactor=1e-5)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Melinder, Ethylene glycol', sFolder='xMass', name='MEG2', desc='Melinder, Ethylene Glycol', ref='Melinder2010,Skovrup2013', densityFactor=None, heatFactor=None, conductivityFactor=None, viscosityFactor=1e-5)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Melinder, Glycerol', sFolder='xMass', name='MGL2', desc='Melinder, Glycerol', ref='Melinder2010,Skovrup2013', densityFactor=None, heatFactor=None, conductivityFactor=None, viscosityFactor=1e-5)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Melinder, Magnesium Chloride', sFolder='xMass', name='MMG2', desc='Melinder, Magnesium Chloride', ref='Melinder2010,Skovrup2013', densityFactor=None, heatFactor=None, conductivityFactor=None, viscosityFactor=1e-5)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Melinder, Methanol', sFolder='xMass', name='MMA2', desc='Melinder, Methanol', ref='Melinder2010,Skovrup2013', densityFactor=None, heatFactor=None, conductivityFactor=None, viscosityFactor=1e-5)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Melinder, Potassium Acetate', sFolder='xMass', name='MKA2', desc='Melinder, Potassium Acetate', ref='Melinder2010,Skovrup2013', densityFactor=None, heatFactor=None, conductivityFactor=None, viscosityFactor=1e-5)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Melinder, Potassium Carbonate', sFolder='xMass', name='MKC2', desc='Melinder, Potassium Carbonate', ref='Melinder2010,Skovrup2013', densityFactor=None, heatFactor=None, conductivityFactor=None, viscosityFactor=1e-5)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Melinder, Propylene Glycol', sFolder='xMass', name='MPG2', desc='Melinder, Propylene Glycol', ref='Melinder2010,Skovrup2013', densityFactor=None, heatFactor=None, conductivityFactor=None, viscosityFactor=1e-5)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Melinder, Sodium Chloride', sFolder='xMass', name='MNA2', desc='Melinder, Sodium Chloride', ref='Melinder2010,Skovrup2013', densityFactor=None, heatFactor=None, conductivityFactor=None, viscosityFactor=1e-5)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='VDI, Calcium Cloride', sFolder='xMass', name='VCA', desc='VDI, Calcium Cloride', ref='Preisegger2010,Skovrup2013', densityFactor=None, heatFactor=None, conductivityFactor=1e-3, viscosityFactor=1e-5)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='VDI, Magnesium Chloride', sFolder='xMass', name='VMG', desc='VDI, Magnesium Chloride', ref='Preisegger2010,Skovrup2013', densityFactor=None, heatFactor=None, conductivityFactor=1e-3, viscosityFactor=1e-5)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='VDI, Methanol', sFolder='xMass', name='VMA', desc='VDI, Methanol', ref='Preisegger2010,Skovrup2013', densityFactor=None, heatFactor=None, conductivityFactor=1e-3, viscosityFactor=1e-5)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='VDI, Potassium Carbonate', sFolder='xMass', name='VKC', desc='VDI, Potassium Carbonate', ref='Preisegger2010,Skovrup2013', densityFactor=None, heatFactor=None, conductivityFactor=1e-3, viscosityFactor=1e-5)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='VDI, Sodium Chloride', sFolder='xMass', name='VNA', desc='VDI, Sodium Chloride', ref='Preisegger2010,Skovrup2013', densityFactor=None, heatFactor=None, conductivityFactor=1e-3, viscosityFactor=1e-5)]
        print(", {0}".format(sec[-1].name), end="")

        sec += [SecCoolSolutionData(sFile='HFE-7100', sFolder='xPure', name='HFE2', desc='HFE-7100, Hydrofluoroether', ref='3M2007,Skovrup2013', densityFactor=None, heatFactor=None, conductivityFactor=None, viscosityFactor=1e-3)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='NBS, Water', sFolder='xPure', name='NBS', desc='NBS, Water', ref='Schmidt1979,Skovrup2013', densityFactor=None, heatFactor=None, conductivityFactor=None, viscosityFactor=1e-6)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Paracryol', sFolder='xPure', name='PCL', desc='Paracryol, Aliphatic Hydrocarbon', ref='Sulzer1999,Skovrup2013', densityFactor=None, heatFactor=1e3, conductivityFactor=None, viscosityFactor=1e-5)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Paratherm NF', sFolder='xPure', name='PNF2', desc='Paratherm NF, Hydrotreated mineral oil', ref='Paratherm2013,Skovrup2013', densityFactor=None, heatFactor=1e3, conductivityFactor=None, viscosityFactor=1e-3)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Tyfoxit 1.10', sFolder='xPure', name='TY10', desc='Tyfoxit 1.10, Potassium Acetate', ref='Tyfoprop1999,Skovrup2013', densityFactor=None, heatFactor=1e3, conductivityFactor=None, viscosityFactor=1e-5)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Tyfoxit 1.15', sFolder='xPure', name='TY15', desc='Tyfoxit 1.15, Potassium Acetate', ref='Tyfoprop1999,Skovrup2013', densityFactor=None, heatFactor=None, conductivityFactor=1e-3, viscosityFactor=1e-5)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Tyfoxit 1.20', sFolder='xPure', name='TY20', desc='Tyfoxit 1.20, Potassium Acetate', ref='Tyfoprop1999,Skovrup2013', densityFactor=None, heatFactor=None, conductivityFactor=1e-3, viscosityFactor=1e-5)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Tyfoxit 1.24', sFolder='xPure', name='TY24', desc='Tyfoxit 1.24, Potassium Acetate', ref='Tyfoprop1999,Skovrup2013', densityFactor=None, heatFactor=None, conductivityFactor=1e-3, viscosityFactor=1e-5)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Zitrec S10', sFolder='xPure', name='ZS10', desc='Zitrec S10, Potassium formate/Sodium propionate', ref='Arteco2010,Skovrup2013', densityFactor=None, heatFactor=1e3, conductivityFactor=None, viscosityFactor=None)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Zitrec S25', sFolder='xPure', name='ZS25', desc='Zitrec S25, Potassium formate/Sodium propionate', ref='Arteco2010,Skovrup2013', densityFactor=None, heatFactor=1e3, conductivityFactor=None, viscosityFactor=None)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Zitrec S40', sFolder='xPure', name='ZS40', desc='Zitrec S40, Potassium formate/Sodium propionate', ref='Arteco2010,Skovrup2013', densityFactor=None, heatFactor=1e3, conductivityFactor=None, viscosityFactor=None)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Zitrec S45', sFolder='xPure', name='ZS45', desc='Zitrec S45, Potassium formate/Sodium propionate', ref='Arteco2010,Skovrup2013', densityFactor=None, heatFactor=1e3, conductivityFactor=None, viscosityFactor=None)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Zitrec S55', sFolder='xPure', name='ZS55', desc='Zitrec S55, Potassium formate/Sodium propionate', ref='Arteco2010,Skovrup2013', densityFactor=None, heatFactor=1e3, conductivityFactor=None, viscosityFactor=None)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Syltherm XLT', sFolder='xPure', name='XLT2', desc='Syltherm XLT, Polydimethylsiloxan', ref='Dow1997,Skovrup2013', densityFactor=None, heatFactor=1e3, conductivityFactor=None, viscosityFactor=1e-3)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Dowtherm J', sFolder='xPure', name='DowJ2', desc='Dowtherm J, Diethylbenzene mixture', ref='Dow1997,Skovrup2013', densityFactor=None, heatFactor=None, conductivityFactor=None, viscosityFactor=1e-5)]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Dowtherm Q', sFolder='xPure', name='DowQ2', desc='Dowtherm Q, Diphenylethane/alkylated aromatics', ref='Dow1997,Skovrup2013', densityFactor=None, heatFactor=None, conductivityFactor=None, viscosityFactor=1e-5)]
        print(", {0}".format(sec[-1].name), end="")

        sec += [SecCoolIceData(sFile='IceEA', sFolder='xMass', name='IceEA', desc='Ice slurry with Ethanol', ref='Kauffeld2001,Skovrup2013')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolIceData(sFile='IceNA', sFolder='xMass', name='IceNA', desc='Ice slurry with NaCl', ref='Kauffeld2001,Skovrup2013')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolIceData(sFile='IcePG', sFolder='xMass', name='IcePG', desc='Ice slurry with Propylene Glycol', ref='Kauffeld2001,Skovrup2013')]
        print(", {0}".format(sec[-1].name), end="")

        sec += [ThermogenVP1869()]
        print(", {0}".format(sec[-1].name), end="")
        sec += [Freezium()]
        print(", {0}".format(sec[-1].name), end="")
        sec += [AS10()]
        print(", {0}".format(sec[-1].name), end="")
        sec += [AS20()]
        print(", {0}".format(sec[-1].name), end="")
        sec += [AS30()]
        print(", {0}".format(sec[-1].name), end="")
        sec += [AS40()]
        print(", {0}".format(sec[-1].name), end="")
        sec += [AS55()]
        print(", {0}".format(sec[-1].name), end="")

        print(" ... done")

        return sec


class SecCoolIceData(SecCoolSolutionData):
    """
    A base class that can be fed with a fluid ID from SecCool
    to read data files sitting in data/SecCool/xTables.
    """

    def __init__(self, sFile=None, sFolder=None, name=None, desc=None, ref='Danish Technological Institute,Skovrup2013'):
        SecCoolSolutionData.__init__(self, sFile=sFile, sFolder=sFolder, name=name, desc=desc, ref=ref)

        #self.density.xData,self.density.yData,self.density.data = self.getArray(dataID="Rho")
        #self.density.source = self.density.SOURCE_DATA

        self.specific_heat.xData, self.specific_heat.yData, self.specific_heat.data = self.getArray(dataID='Hfusion')
        self.specific_heat.source = self.specific_heat.SOURCE_DATA

        #self.conductivity.xData,self.conductivity.yData,self.conductivity.data   = self.getArray(dataID='Cond')
        #self.conductivity.source = self.conductivity.SOURCE_DATA

        #self.viscosity.xData,self.viscosity.yData,self.viscosity.data   = self.getArray(dataID='Mu')
        #self.viscosity.source = self.viscosity.SOURCE_DATA


#    def fitFluid(self):
#        if self.Tbase==None:
#            self.Tbase = (self.Tmin + self.Tmax) / 2.0
#        if self.xbase==None:
#            self.xbase = (self.xmin + self.xmax) / 2.0
#
#        std_coeffs = np.zeros((4,6))
#        errList    = (ValueError, AttributeError, TypeError, RuntimeError)
#
#        try:
#            self.density.coeffs = np.copy(std_coeffs)
#            self.density.type   = self.density.INCOMPRESSIBLE_POLYNOMIAL
#            self.density.fitCoeffs(self.Tbase,self.xbase)
#        except errList as ve:
#            if self.density.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(self.name,'density',ve))
#            pass
#
#        try:
#            self.specific_heat.coeffs = np.copy(std_coeffs)
#            self.specific_heat.type   = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
#            self.specific_heat.fitCoeffs(self.Tbase,self.xbase)
#        except errList as ve:
#            if self.specific_heat.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(self.name,'specific heat',ve))
#            pass
#
#        try:
#            self.conductivity.coeffs = np.copy(std_coeffs)
#            self.conductivity.type   = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
#            self.conductivity.fitCoeffs(self.Tbase,self.xbase)
#        except errList as ve:
#            if self.conductivity.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(self.name,'conductivity',ve))
#            pass
#
#        try:
#            self.viscosity.coeffs = np.copy(std_coeffs)
#            self.viscosity.type   = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
#            self.viscosity.fitCoeffs(self.Tbase,self.xbase)
#        except errList as ve:
#            if self.viscosity.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(self.name,'viscosity',ve))
#            pass

    # Redefine some functions to avoid data loss

    def getFile(self, data):
        return os.path.join(os.path.dirname(__file__), 'data', 'SecCool', 'xTables', self.sFolder, self.sFile + "_" + data + ".csv")

    def getFromFile(self, data):
        fullPath = self.getFile(data)
        content = np.loadtxt(fullPath, dtype=type('string'), delimiter=",", skiprows=3)
        r, c, res = IncompressibleFitter.shapeArray(content)
#        numbers = res.astype(np.float)
        numbers = np.zeros((r, c))
        for i in range(r):
            for j in range(c):
                nu = np.NAN
                try:
                    nu = np.float(res[i, j])
                    if i == 0: nu *= 1e-2  # Percent to fraction
                    if not self.allowNegativeData and nu < 0:
                        nu = np.NAN  # invalid entries
                except (ValueError, TypeError) as ve:
                    if False: print("Could not convert entry: {0}".format(ve))
                    if i == 0: nu = 0.0  # Dummy for tables without concentration (TFreeze and Vol2Mass)
                    pass
                numbers[i, j] = nu
        if numbers[1, 0] > numbers[-1, 0]:
            numbers[1:, :] = numbers[1:, :][::-1, :]
        if numbers[0, 1] > numbers[0, -1]:
            numbers[:, 1:] = numbers[:, 1:][:, ::-1]
        return numbers


class ThermogenVP1869(PureData, DigitalData):
    """
    Source: SecCool Software
    """

    def __init__(self):
        DigitalData.__init__(self)
        PureData.__init__(self)
        self.name = "TVP1869"
        self.description = "Thermogen VP 1869"
        self.reference = "Hoechst1995,Skovrup2013"

        self.Tmax = 20 + 273.15
        self.Tmin = -80 + 273.15
        self.TminPsat = self.Tmax

        self.Tbase = 0.00 + 273.15

        self.density.source = self.density.SOURCE_COEFFS
        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = np.array([[945.5454545], [-1.054545455]])

        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = np.array([[2.322218182], [0.003843636]]) * 1000

        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = np.array([[0.15], [-0.000154545]])

        self.temperature.data = self.getTrange()
        self.concentration.data = np.array([0.0])  # mass fraction

    def fitFluid(self):
        key = 'Mu'

        def funcMu(T, x):
            T = T - self.Tbase
            return (341.3688975 + T * (-0.713408301 + 0.017723992 * T)) / (1.0 + T * (0.034502393 + T * (0.000401319 + 1.57288E-06 * T))) * 1e-2 * 1e-3

        self.viscosity.xData, self.viscosity.yData, self.viscosity.data = self.getArray(dataID=key, func=funcMu, x_in=self.temperature.data, y_in=self.concentration.data)

        try:
            self.viscosity.source = self.viscosity.SOURCE_EQUATION
            self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
            self.viscosity.coeffs = np.zeros((4, 6))
            self.viscosity.fitCoeffs(self.Tbase, self.xbase)
        except (ValueError, AttributeError, TypeError, RuntimeError) as e:
            if self.viscosity.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(self.name, 'viscosity', e))
            pass


class Freezium(DigitalData):
    def __init__(self):
        DigitalData.__init__(self)

        self.name = "FRE"
        self.description = "Freezium, Potassium Formate"
        self.reference = "Kemira1998,Skovrup2013"

        self.Tmin = -40 + 273.00
        self.Tmax = +40 + 273.00
        self.xmax = 0.50
        self.xmin = 0.19
        self.xid = self.ifrac_mass
        self.TminPsat = self.Tmax

        self.Tbase = 273.15
        self.xbase = 0.0

        self.temperature.data = self.getTrange()
        self.concentration.data = self.getxrange()

        self.density.source = self.density.SOURCE_COEFFS
        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = np.array([[1015, 462, 406], [-40 / 100.0, 0.0, 0.0]])

        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = np.array([[0.55, -0.15], [0.18 / 100.0, -0.16 / 100.0]])

    def fitFluid(self):

        key = 'Cp'

        def funcCp(T, x):
            T = (T - self.Tbase) / 100.0
            return (4.15 * np.exp(-0.9 * x) + 0.63 * T * x) * 1000.0

        self.specific_heat.xData, self.specific_heat.yData, self.specific_heat.data = self.getArray(dataID=key, func=funcCp, x_in=self.temperature.data, y_in=self.concentration.data, DEBUG=self.specific_heat.DEBUG)

        try:
            self.specific_heat.source = IncompressibleData.SOURCE_EQUATION
            self.specific_heat.type = IncompressibleData.INCOMPRESSIBLE_POLYNOMIAL
            self.specific_heat.coeffs = np.zeros((4, 6))
            self.specific_heat.fitCoeffs(self.Tbase, self.xbase)
        except (ValueError, AttributeError, TypeError, RuntimeError) as e:
            if self.specific_heat.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(self.name, 'specific heat', e))
            pass
        funcCp = None

        key = 'Mu'

        def funcMu(T, x):
            Tr = (T - self.Tbase) / 100.0
            result = 0.32 + x * (-0.70 + x * 2.26) + Tr * (-1.26 + Tr * (1.12 - Tr * 0.894))
            return self.rho(T, 1e6, x) * np.power(10, result) * 1E-6;

        self.viscosity.xData, self.viscosity.yData, self.viscosity.data = self.getArray(dataID=key, func=funcMu, x_in=self.temperature.data, y_in=self.concentration.data, DEBUG=self.viscosity.DEBUG)

        try:
            self.viscosity.source = IncompressibleData.SOURCE_EQUATION
            self.viscosity.type = IncompressibleData.INCOMPRESSIBLE_EXPPOLYNOMIAL
            self.viscosity.coeffs = np.zeros((4, 6))
            self.viscosity.fitCoeffs(self.Tbase, self.xbase)
        except (ValueError, AttributeError, TypeError, RuntimeError) as e:
            if self.viscosity.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(self.name, 'viscosity', e))
            pass
        funcMu = None

        # Changed the coefficient order for TFreeze
        key = 'Tfreeze'

        def funcTf(T, x):
            x = x * 100.0
            a = 0.03422039835160944
            b = -0.05425629002714395
            c = -0.007991085555390726
            d = 0.001036937163846157
            e = 0.0003268582531827402
            f = -7.721483884155584E-06
            g = -4.841293026057464E-06
            h = 1.216929917247388E-08
            i = 2.41704396074865E-08
            j = 4.314861246570078E-11
            return ((a + x * (c + x * (e + x * (g + i * x)))) / (1 + x * (b + x * (d + x * (f + x * (h + j * x))))) * 100) + 273.15

        self.T_freeze.xData, self.T_freeze.yData, self.T_freeze.data = self.getArray(dataID=key, func=funcTf, x_in=np.array([0.0]), y_in=self.concentration.data, DEBUG=self.T_freeze.DEBUG)
        try:
            self.T_freeze.coeffs = np.zeros((4, 6))
            self.T_freeze.source = IncompressibleData.SOURCE_EQUATION
            self.T_freeze.type = IncompressibleData.INCOMPRESSIBLE_POLYNOMIAL
            self.T_freeze.fitCoeffs(self.Tbase, self.xbase)
        except (ValueError, AttributeError, TypeError, RuntimeError) as e:
            if self.T_freeze.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(self.name, 'T freeze', e))
            pass
        funcTf = None


class AS10(PureData, DigitalData):
    def __init__(self):
        DigitalData.__init__(self)
        PureData.__init__(self)
        self.name = "AS10"
        self.description = "Aspen Temper -10, Potassium acetate/formate"
        self.reference = "Aspen2001,Skovrup2013"

        self.Tmax = 30 + 273.15
        self.Tmin = -10 + 273.15
        self.TminPsat = self.Tmax
        self.Tbase = 0.00 + 273.15
        self.temperature.data = self.getTrange()

        self.density.source = self.density.SOURCE_COEFFS
        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = np.array([[-1.89735e-18], [+ 1.66533e-16], [- 0.2], [+ 1090]])[::-1]

        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = np.array([[-0.0132589], [+ 2.01], [+ 3541.83]])[::-1]

        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = np.array([[0.001483], [+ 0.514259]])[::-1]

        self.viscosity.source = self.viscosity.SOURCE_COEFFS
        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_POLYNOMIAL
        self.viscosity.coeffs = np.array([[-2.11481e-5], [+ 0.00235381], [- 0.10631376], [+ 2.80154921]])[::-1]
        self.viscosity.coeffs *= 1e-3

    def fitFluid(self):
        pass


class AS20(PureData, DigitalData):
    def __init__(self):
        DigitalData.__init__(self)
        PureData.__init__(self)
        self.name = "AS20"
        self.description = "Aspen Temper -20, Potassium acetate/formate"
        self.reference = "Aspen2001,Skovrup2013"

        self.Tmax = 30 + 273.15
        self.Tmin = -20 + 273.15
        self.TminPsat = self.Tmax
        self.Tbase = 0.00 + 273.15
        self.temperature.data = self.getTrange()

        self.density.source = self.density.SOURCE_COEFFS
        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = np.array([[2.98156e-19], [- 0.00142857], [- 0.22142857], [+ 1147]])[::-1]

        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = np.array([[-0.0122673], [+ 2.86], [+ 3262.52]])[::-1]

        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = np.array([[0.001342], [+ 0.480766]])[::-1]

    def fitFluid(self):
        key = 'Mu'

        def funcMu(T, x):
            T = (T - self.Tbase)
            mPas = 2.43708721027941 * np.exp(-0.0537593944541809 * T) + 0.97244
            return mPas / 1e3

        self.viscosity.xData, self.viscosity.yData, self.viscosity.data = self.getArray(dataID=key, func=funcMu, x_in=self.temperature.data, y_in=self.concentration.data, DEBUG=self.viscosity.DEBUG)

        try:
            self.viscosity.source = self.viscosity.SOURCE_EQUATION
            self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
            self.viscosity.coeffs = np.zeros((4, 6))
            self.viscosity.fitCoeffs(self.Tbase, self.xbase)
        except (ValueError, AttributeError, TypeError, RuntimeError) as e:
            if self.viscosity.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(self.name, 'viscosity', e))
            pass
        funcMu = None


class AS30(PureData, DigitalData):
    def __init__(self):
        DigitalData.__init__(self)
        PureData.__init__(self)
        self.name = "AS30"
        self.description = "Aspen Temper -30, Potassium acetate/formate"
        self.reference = "Aspen2001,Skovrup2013"

        self.Tmax = 30 + 273.15
        self.Tmin = -30 + 273.15
        self.TminPsat = self.Tmax
        self.Tbase = 0.00 + 273.15
        self.temperature.data = self.getTrange()

        self.density.source = self.density.SOURCE_COEFFS
        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = np.array([[1.76768e-5], [- 0.00103896], [- 0.31085859], [+ 1183.85930736]])[::-1]

        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = np.array([[-0.0270232], [+ 2.99], [+ 3075.04]])[::-1]

        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = np.array([[0.001256], [+ 0.460388]])[::-1]

    def fitFluid(self):
        key = 'Mu'

        def funcMu(T, x):
            T = (T - self.Tbase)
            mPas = 2.65653950695888 * np.exp(-0.0598806339442954 * T) + 1.30143
            return mPas / 1e3

        self.viscosity.xData, self.viscosity.yData, self.viscosity.data = self.getArray(dataID=key, func=funcMu, x_in=self.temperature.data, y_in=self.concentration.data, DEBUG=self.viscosity.DEBUG)

        try:
            self.viscosity.source = self.viscosity.SOURCE_EQUATION
            self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
            self.viscosity.coeffs = np.zeros((4, 6))
            self.viscosity.fitCoeffs(self.Tbase, self.xbase)
        except (ValueError, AttributeError, TypeError, RuntimeError) as e:
            if self.viscosity.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(self.name, 'viscosity', e))
            pass
        funcMu = None


class AS40(PureData, DigitalData):
    def __init__(self):
        DigitalData.__init__(self)
        PureData.__init__(self)
        self.name = "AS40"
        self.description = "Aspen Temper -40, Potassium acetate/formate"
        self.reference = "Aspen2001,Skovrup2013"

        self.Tmax = 30 + 273.15
        self.Tmin = -40 + 273.15
        self.TminPsat = self.Tmax
        self.Tbase = 0.00 + 273.15
        self.temperature.data = self.getTrange()

        self.density.source = self.density.SOURCE_COEFFS
        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = np.array([[1.6835e-5], [- 0.00109307], [- 0.37819865], [+ 1214.83982684]])[::-1]

        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = np.array([[-0.0387227], [+ 2.28], [+ 2977.88]])[::-1]

        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = np.array([[0.001099], [+ 0.443327]])[::-1]

    def fitFluid(self):
        key = 'Mu'

        def funcMu(T, x):
            T = (T - self.Tbase)
            mPas = 0.714976365635003 * np.exp(-0.100050525515385 * T) + 4.38768154440393 * np.exp(-0.0260039000649317 * T)
            return mPas / 1e3

        self.viscosity.xData, self.viscosity.yData, self.viscosity.data = self.getArray(dataID=key, func=funcMu, x_in=self.temperature.data, y_in=self.concentration.data, DEBUG=self.viscosity.DEBUG)

        try:
            self.viscosity.source = self.viscosity.SOURCE_EQUATION
            self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
            self.viscosity.coeffs = np.zeros((4, 6))
            self.viscosity.fitCoeffs(self.Tbase, self.xbase)
        except (ValueError, AttributeError, TypeError, RuntimeError) as e:
            if self.viscosity.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(self.name, 'viscosity', e))
            pass
        funcMu = None


class AS55(PureData, DigitalData):
    def __init__(self):
        DigitalData.__init__(self)
        PureData.__init__(self)
        self.name = "AS55"
        self.description = "Aspen Temper -55, Potassium acetate/formate"
        self.reference = "Aspen2001,Skovrup2013"

        self.Tmax = 30 + 273.15
        self.Tmin = -55 + 273.15
        self.TminPsat = self.Tmax
        self.Tbase = 0.00 + 273.15
        self.temperature.data = self.getTrange()

        self.density.source = self.density.SOURCE_COEFFS
        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = np.array([[1.98824e-5], [- 0.00117189], [- 0.47629615], [+ 1249.7534665]])[::-1]

        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = np.array([[-0.0248618], [+ 2.29], [+ 2839.85]])[::-1]

        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = np.array([[2.287e-6], [+ 0.000937902], [+ 0.425799423]])[::-1]

    def fitFluid(self):
        key = 'Mu'

        def funcMu(T, x):
            T = (T - self.Tbase)
            mPas = 0.159583223482554 * np.exp(-0.138097704125669 * T) + 6.3176967296442 * np.exp(-0.0380509974688477 * T)
            return mPas / 1e3

        self.viscosity.xData, self.viscosity.yData, self.viscosity.data = self.getArray(dataID=key, func=funcMu, x_in=self.temperature.data, y_in=self.concentration.data, DEBUG=self.viscosity.DEBUG)

        try:
            self.viscosity.source = self.viscosity.SOURCE_EQUATION
            self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
            self.viscosity.coeffs = np.zeros((4, 6))
            self.viscosity.fitCoeffs(self.Tbase, self.xbase)
        except (ValueError, AttributeError, TypeError, RuntimeError) as e:
            if self.viscosity.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(self.name, 'viscosity', e))
            pass
        funcMu = None
