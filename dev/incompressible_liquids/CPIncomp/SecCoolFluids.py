from __future__ import division, print_function
import numpy as np
from BaseObjects import IncompressibleData
from DataObjects import DigitalData, PureData 
import os, sys
from scipy import interpolate 
from CPIncomp.BaseObjects import IncompressibleFitter

class SecCoolSolutionData(DigitalData):
    """ 
    A base class that can be fed with a fluid ID from SecCool
    to read data files sitting in data/SecCool/xMass and 
    data/SecCool/xVolume. 
    """ 
    def __init__(self,sFile=None,sFolder=None,name=None,desc=None,ref='SecCool software'):
        DigitalData.__init__(self)
        self.allowNegativeData = False
        
        if sFile   == None: raise ValueError("File name cannot be empty.")
        if sFolder == None: raise ValueError("Folder name cannot be empty.")
        if name    == None: raise ValueError("Fluid name cannot be empty.")
        if desc    == None: raise ValueError("Description cannot be empty.")
        
        self.sFile   = sFile
        self.sFolder = sFolder
        
        self.name = name
        self.description = desc
        self.reference = ref
                
        self.temperature.data,self.concentration.data,self.density.data = self.getArray(dataID='Rho')
                
        self.Tmax = np.max(self.temperature.data)
        self.Tmin = np.min(self.temperature.data)
        self.xmax = np.max(self.concentration.data)
        self.xmin = np.min(self.concentration.data)
        if self.sFolder=='xVolume':
            self.xid  = self.ifrac_volume
        elif self.sFolder=='xMass':
            self.xid  = self.ifrac_mass
        elif self.sFolder=='xPure':
            self.xid  = self.ifrac_pure
        else:
            raise ValueError("Unknown folder type specified.")        
        self.TminPsat = self.Tmax
        
        
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
        
        if self.Tbase==None:
            self.Tbase = (self.Tmin + self.Tmax) / 2.0
        if self.xbase==None:
            self.xbase = (self.xmin + self.xmax) / 2.0
        
        std_coeffs = np.zeros((4,6))
        errList    = (ValueError, AttributeError, TypeError, RuntimeError)
        
        try:
            self.density.xData,self.density.yData,self.density.data = self.getArray(dataID="Rho")
            self.density.coeffs = np.copy(std_coeffs)
            self.density.source = self.density.SOURCE_DATA
            self.density.type   = self.density.INCOMPRESSIBLE_POLYNOMIAL
            self.density.fitCoeffs(self.Tbase,self.xbase)
        except errList as ve:
            if self.density.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(self.name,'density',ve))
            pass
        
        try:
            self.specific_heat.xData,self.specific_heat.yData,self.specific_heat.data = self.getArray(dataID='Cp')
            while np.max(self.specific_heat.data[np.isfinite(self.specific_heat.data)])<500: # Expect values around 1e3
                self.specific_heat.data *= 1e3 
            self.specific_heat.coeffs = np.copy(std_coeffs)
            self.specific_heat.source = self.specific_heat.SOURCE_DATA
            self.specific_heat.type   = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
            self.specific_heat.fitCoeffs(self.Tbase,self.xbase)          
        except errList as ve:
            if self.specific_heat.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(self.name,'specific heat',ve))
            pass 
        
        try:
            self.conductivity.xData,self.conductivity.yData,self.conductivity.data   = self.getArray(data='Cond')
            while np.max(self.conductivity.data[np.isfinite(self.conductivity.data)])>10: # Expect value below 1
                self.conductivity.data *= 1e-3 
            self.conductivity.coeffs = np.copy(std_coeffs)
            self.conductivity.source = self.conductivity.SOURCE_DATA
            self.conductivity.type   = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
            self.conductivity.fitCoeffs(self.Tbase,self.xbase)
        except errList as ve:
            if self.conductivity.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(self.name,'conductivity',ve))
            pass
        
        try:
            self.viscosity.xData,self.viscosity.yData,self.viscosity.data   = self.getArray(data='Mu')
            while np.max(self.viscosity.data[np.isfinite(self.viscosity.data)])>100: # Expect value below 10
                self.viscosity.data *= 1e-3 
            self.viscosity.coeffs = np.copy(std_coeffs)
            self.viscosity.source = self.viscosity.SOURCE_DATA
            self.viscosity.type   = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
            self.viscosity.fitCoeffs(self.Tbase,self.xbase)
        except errList as ve:
            if self.viscosity.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(self.name,'viscosity',ve))
            pass
        
        # reset data for getArray and read special files
        if self.xid!=self.ifrac_pure and self.xid!=self.ifrac_undefined:
            allowNegativeData_org = self.allowNegativeData
            self.allowNegativeData = True
            
            x,_,z = self.getArray(dataID='TFreeze')
            self.T_freeze.yData = (x - 273.15) / 100.0
            self.T_freeze.xData = [0.0]
            if np.min(z)<150: z += 273.15 
            self.T_freeze.data = z.T
            
            try:
                self.T_freeze.source = self.T_freeze.SOURCE_DATA
                if np.isfinite(self.T_freeze.data).sum()<2:
                    self.T_freeze.coeffs = np.array([+7e+6, +6e+4, +1e+1])
                    self.T_freeze.type   = self.T_freeze.INCOMPRESSIBLE_EXPONENTIAL
                else:   
                    self.T_freeze.coeffs = np.zeros(np.round(np.array(std_coeffs.shape) * 2))
                    self.T_freeze.type   = self.T_freeze.INCOMPRESSIBLE_EXPPOLYNOMIAL
                self.T_freeze.fitCoeffs(self.Tbase,self.xbase)
            except errList as ve:
                if self.T_freeze.DEBUG: print("{0}: Could not fit {1} coefficients: {2}".format(self.name,"T_freeze",ve))
                pass

            # reset data for getArray again
            self.allowNegativeData = allowNegativeData_org
            x,_,z = self.getArray(dataID='Vol2Mass')            
            massData =    z[:,0]    /100.0
            volData  = (x - 273.15) /100.0
            if self.xid==self.ifrac_volume:
                _,_,self.mass2input.data = IncompressibleFitter.shapeArray(volData,axs=1)
                self.mass2input.xData = [0.0]
                self.mass2input.yData = massData
                try:
                    self.mass2input.coeffs = np.copy(std_coeffs)
                    self.mass2input.source = self.mass2input.SOURCE_DATA
                    self.mass2input.type   = self.mass2input.INCOMPRESSIBLE_POLYNOMIAL
                    self.mass2input.fitCoeffs(self.Tbase,self.xbase)
                except errList as ve:
                    if self.mass2input.DEBUG: print("{0}: Could not fit {1} coefficients: {2}".format(self.name,"mass2input",ve))
                    pass
            elif self.xid==self.ifrac_mass:
                _,_,self.volume2input.data = IncompressibleFitter.shapeArray(massData,axs=1)
                self.volume2input.xData = [0.0]
                self.volume2input.yData = volData
                try:
                    self.volume2input.coeffs = np.copy(std_coeffs)
                    self.volume2input.source = self.volume2input.SOURCE_DATA
                    self.volume2input.type   = self.volume2input.INCOMPRESSIBLE_POLYNOMIAL
                    self.volume2input.fitCoeffs(self.Tbase,self.xbase)
                except errList as ve:
                    if self.volume2input.DEBUG: print("{0}: Could not fit {1} coefficients: {2}".format(self.name,"volume2input",ve))
                    pass
            else:
                raise ValueError("Unknown xid specified.")

        
    # Redefine some functions to avoid data loss
    def getFile(self, data):
        return os.path.join(os.path.dirname(__file__), 'data','SecCool', self.sFolder, self.sFile+"_"+data+".txt")
        
    
    def getFromFile(self, data):
        fullPath = self.getFile(data)
        r,c,res = IncompressibleFitter.shapeArray(np.loadtxt(fullPath,dtype=type('string')))
#        numbers = res.astype(np.float)
        numbers = np.zeros((r,c))
        for i in range(r):
            for j in range(c):
                nu = np.NAN
                try:
                    nu = np.float(res[i,j])
                    if i==0: nu *= 1e-2 # Percent to fraction
                    if j==0: nu += 273.15 # Celsius to Kelvin
                    if not self.allowNegativeData and nu<0:
                        nu = np.NAN # invalid entries
                except (ValueError, TypeError) as ve:
                    if False: print("Could not convert entry: {0}".format(ve))
                    if i==0: nu = 0.0 # Dummy for tables without concentration (TFreeze and Vol2Mass)
                    pass
                numbers[i,j] = nu
        if numbers[1,0]>numbers[-1,0]: 
            numbers[1:,:] = numbers[1:,:][::-1,:]
        if numbers[0,1]>numbers[0,-1]: 
            numbers[:,1:] = numbers[:,1:][:,::-1]
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
        sec += [SecCoolSolutionData(sFile='Antifrogen KF'           ,sFolder='xVolume',name='AKF',desc='Antifrogen KF, Potassium Formate'       ,ref='Clariant GmbH Jan. 2000, SecCool software')]
        print("{0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Antifrogen L'            ,sFolder='xVolume',name='AL' ,desc='Antifrogen L, Propylene Glycol'         ,ref='Clariant GmbH Jan. 2000, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Antifrogen N'            ,sFolder='xVolume',name='AN' ,desc='Antifrogen N, Ethylene Glycol'          ,ref='Clariant GmbH Jan. 2000, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='ASHRAE, Ethylene Glycol' ,sFolder='xVolume',name='AEG',desc='ASHRAE, Ethylene Glycol'                ,ref='ASHRAE Fundamentals Handbook 2001, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='ASHRAE, Propylene Glycol',sFolder='xVolume',name='APG',desc='ASHRAE, Propylene Glycol'               ,ref='ASHRAE Fundamentals Handbook 2001, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Glykosol N'              ,sFolder='xVolume',name='GKN',desc='Glykosol N, Ethylene Glycol'            ,ref='pro KUEHLSOLE GmbH, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Pekasol 2000'            ,sFolder='xVolume',name='PK2',desc='Pekasol 2000, Potassium acetate/formate',ref='pro KUEHLSOLE GmbH, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Pekasol L'               ,sFolder='xVolume',name='PKL',desc='Pekasol L, Propylene Glycol'            ,ref='pro KUEHLSOLE GmbH, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Zitrec AC'               ,sFolder='xVolume',name='ZAC',desc='Zitrec AC, Corrosion Inhibitor'         ,ref='Arteco, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Zitrec FC'               ,sFolder='xVolume',name='ZFC',desc='Zitrec FC, Propylene Glycol'            ,ref='Arteco, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Zitrec LC'               ,sFolder='xVolume',name='ZLC',desc='Zitrec LC, Propylene Glycol'            ,ref='Arteco, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Zitrec MC'               ,sFolder='xVolume',name='ZMC',desc='Zitrec MC, Ethylene Glycol'             ,ref='Arteco, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Zitrec M'                ,sFolder='xVolume',name='ZM' ,desc='Zitrec M, Ethylene Glycol'              ,ref='Arteco, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        
        sec += [SecCoolSolutionData(sFile='Melinder, Ammonia'            ,sFolder='xMass',name='MAM2',desc='Melinder, Ammonia'            ,ref='Melinder-BOOK-2010, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Melinder, Calcium Cloride'    ,sFolder='xMass',name='MCA2',desc='Melinder, Calcium Chloride'   ,ref='Melinder-BOOK-2010, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Melinder, Ethanol'            ,sFolder='xMass',name='MEA2',desc='Melinder, Ethanol'            ,ref='Melinder-BOOK-2010, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Melinder, Ethylene glycol'    ,sFolder='xMass',name='MEG2',desc='Melinder, Ethylene Glycol'    ,ref='Melinder-BOOK-2010, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Melinder, Glycerol'           ,sFolder='xMass',name='MGL2',desc='Melinder, Glycerol'           ,ref='Melinder-BOOK-2010, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Melinder, Magnesium Chloride' ,sFolder='xMass',name='MMG2',desc='Melinder, Magnesium Chloride' ,ref='Melinder-BOOK-2010, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Melinder, Methanol'           ,sFolder='xMass',name='MMA2',desc='Melinder, Methanol'           ,ref='Melinder-BOOK-2010, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Melinder, Potassium Acetate'  ,sFolder='xMass',name='MKA2',desc='Melinder, Potassium Acetate'  ,ref='Melinder-BOOK-2010, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Melinder, Potassium Carbonate',sFolder='xMass',name='MKC2',desc='Melinder, Potassium Carbonate',ref='Melinder-BOOK-2010, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Melinder, Propylene Glycol'   ,sFolder='xMass',name='MPG2',desc='Melinder, Propylene Glycol'   ,ref='Melinder-BOOK-2010, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Melinder, Sodium Chloride'    ,sFolder='xMass',name='MNA2',desc='Melinder, Sodium Chloride'    ,ref='Melinder-BOOK-2010, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='VDI, Calcium Cloride'         ,sFolder='xMass',name='VCA' ,desc='VDI, Calcium Cloride'         ,ref='VDI Waermeatlas 9th Edition 2002, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='VDI, Magnesium Chloride'      ,sFolder='xMass',name='VMG' ,desc='VDI, Magnesium Chloride'      ,ref='VDI Waermeatlas 9th Edition 2002, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='VDI, Methanol'                ,sFolder='xMass',name='VMA' ,desc='VDI, Methanol'                ,ref='VDI Waermeatlas 9th Edition 2002, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='VDI, Potassium Carbonate'     ,sFolder='xMass',name='VKC' ,desc='VDI, Potassium Carbonate'     ,ref='VDI Waermeatlas 9th Edition 2002, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='VDI, Sodium Chloride'         ,sFolder='xMass',name='VNA' ,desc='VDI, Sodium Chloride'         ,ref='VDI Waermeatlas 9th Edition 2002, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")

        sec += [SecCoolSolutionData(sFile='HFE-7100'     ,sFolder='xPure',name='HFE2' ,desc='HFE-7100, Hydrofluoroether'                      ,ref='3M Novec, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='NBS, Water'   ,sFolder='xPure',name='NBS' ,desc='NBS, Water'                                      ,ref='Properties of Water and Steam in SI-Units, 2nd Revised and Updated Printing, Springer 1979, pp. 175 ff., SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Paracryol'    ,sFolder='xPure',name='PCL' ,desc='Paracryol, Aliphatic Hydrocarbon'                ,ref='Sulzer Chemtech AG, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Paratherm NF' ,sFolder='xPure',name='PNF' ,desc='Paratherm NF, Hydrotreated mineral oil'          ,ref='Paratherm Ltd, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Tyfoxit 1.10' ,sFolder='xPure',name='TY10',desc='Tyfoxit 1.10, Potassium Acetate'                 ,ref='Tyforop Chemie Gmbh - Technical information 09/99, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Tyfoxit 1.15' ,sFolder='xPure',name='TY15',desc='Tyfoxit 1.15, Potassium Acetate'                 ,ref='Tyforop Chemie Gmbh - Technical information 09/99, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Tyfoxit 1.20' ,sFolder='xPure',name='TY20',desc='Tyfoxit 1.20, Potassium Acetate'                 ,ref='Tyforop Chemie Gmbh - Technical information 09/99, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Tyfoxit 1.24' ,sFolder='xPure',name='TY24',desc='Tyfoxit 1.24, Potassium Acetate'                 ,ref='Tyforop Chemie Gmbh - Technical information 09/99, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Zitrec S10'   ,sFolder='xPure',name='ZS10',desc='Zitrec S10, Potassium formate/Sodium propionate' ,ref='Arteco, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Zitrec S25'   ,sFolder='xPure',name='ZS25',desc='Zitrec S25, Potassium formate/Sodium propionate' ,ref='Arteco, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Zitrec S40'   ,sFolder='xPure',name='ZS40',desc='Zitrec S40, Potassium formate/Sodium propionate' ,ref='Arteco, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Zitrec S45'   ,sFolder='xPure',name='ZS45',desc='Zitrec S45, Potassium formate/Sodium propionate' ,ref='Arteco, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Zitrec S55'   ,sFolder='xPure',name='ZS55',desc='Zitrec S55, Potassium formate/Sodium propionate' ,ref='Arteco, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Syltherm XLT'   ,sFolder='xPure',name='XLT2',desc='Syltherm XLT, Polydimethylsiloxan' ,ref='Dow Chemicals, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Dowtherm J'   ,sFolder='xPure',name='DowJ2',desc='Dowtherm J, Diethylbenzene mixture' ,ref='Dow Chemicals, SecCool software')]
        print(", {0}".format(sec[-1].name), end="")
        sec += [SecCoolSolutionData(sFile='Dowtherm Q'   ,sFolder='xPure',name='DowQ2',desc='Dowtherm Q, Diphenylethane/alkylated aromatics' ,ref='Dow Chemicals, SecCool software')]
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


class ThermogenVP1869(PureData,DigitalData):
    """ 
    Source: SecCool Software
    """ 
    def __init__(self):        
        DigitalData.__init__(self)
        PureData.__init__(self) 
        self.name = "TVP1869"
        self.description = "Thermogen VP 1869"
        self.reference = "Hoechst, SecCool software"
        
        self.Tmax =  20 + 273.15
        self.Tmin = -80 + 273.15
        self.TminPsat =  self.Tmax
    
        self.Tbase =   0.00 + 273.15
       
        self.density.source = self.density.SOURCE_COEFFS
        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = np.array([[945.5454545],[-1.054545455]])
        
        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = np.array([[2.322218182],[0.003843636]])*1000

        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = np.array([[0.15],[-0.000154545]])
        
        self.temperature.data         = self.getTrange()
        self.concentration.data       = np.array([0]) # mass fraction
        
        
    def fitFluid(self):        
        key = 'Mu'
        def funcMu(T,x):
            T = T-self.Tbase
            return (341.3688975+T*(-0.713408301+0.017723992*T))/ \
              (1+T*(0.034502393+T*(0.000401319+1.57288E-06*T)))*1e-2

        self.viscosity.xData,self.viscosity.yData,self.viscosity.data = self.getArray(dataID=key,func=funcMu,x_in=self.temperature.data,y_in=self.concentration.data)
        
        try:
            self.viscosity.source = self.viscosity.SOURCE_EQUATION
            self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
            self.viscosity.coeffs = np.zeros((4,6))
            self.viscosity.fitCoeffs(self.temperature.data,self.concentration.data,self.Tbase,self.xbase)
        except (ValueError, AttributeError, TypeError, RuntimeError) as e:
            if self.viscosity.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(self.name,'viscosity',e))
            pass 


        
class Freezium(DigitalData):
    def __init__(self):
        DigitalData.__init__(self) 
        
        self.name = "FRE"
        self.description = "Freezium, Potassium Formate"
        self.reference = "Kemira Chemicals OY, SecCool software"
        
        self.Tmin = -40 + 273.00
        self.Tmax = +40 + 273.00
        self.xmax = 0.50
        self.xmin = 0.19
        self.xid  = self.ifrac_mass
        self.TminPsat = self.Tmax
        
        self.Tbase = 273.15
        self.xbase = 0.0
        
        self.temperature.data         = self.getTrange()
        self.concentration.data       = self.getxrange()
        
        self.density.source = self.density.SOURCE_COEFFS
        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = np.array([[1015,462,406],[-40/100.0,0.0,0.0]])

        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = np.array([[0.55,-0.15],[0.18/100.0,-0.16/100.0]])
        
        
    def fitFluid(self):
            
        key = 'Cp'
        def funcCp(T,x):
            T = (T-self.Tbase)/100.0
            return (4.15*np.exp(-0.9*x)+0.63*T*x)*1000.0

        self.specific_heat.xData,self.specific_heat.yData,self.specific_heat.data = self.getArray(dataID=key,func=funcCp,x_in=self.temperature.data,y_in=self.concentration.data,DEBUG=self.specific_heat.DEBUG)
        
        try:
            self.specific_heat.source = IncompressibleData.SOURCE_EQUATION
            self.specific_heat.type = IncompressibleData.INCOMPRESSIBLE_POLYNOMIAL
            self.specific_heat.coeffs = np.zeros((4,6))
            self.specific_heat.fitCoeffs(self.Tbase,self.xbase)
        except (ValueError, AttributeError, TypeError, RuntimeError) as e:
            if self.specific_heat.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(self.name,'specific heat',e))
            pass 
        funcCp = None
        

        key = 'Mu'
        def funcMu(T,x):
            Tr = (T-self.Tbase)/100.0
            result = 0.32+x*(-0.70+x*2.26)+Tr*(-1.26+Tr*(1.12-Tr*0.894))
            return self.rho(T, 1e6, x)*np.power(10,result)*1E-3;

        self.viscosity.xData,self.viscosity.yData,self.viscosity.data = self.getArray(dataID=key,func=funcMu,x_in=self.temperature.data,y_in=self.concentration.data,DEBUG=self.viscosity.DEBUG)
        
        try:
            self.viscosity.source = IncompressibleData.SOURCE_EQUATION
            self.viscosity.type = IncompressibleData.INCOMPRESSIBLE_EXPPOLYNOMIAL
            self.viscosity.coeffs = np.zeros((4,6))
            self.viscosity.fitCoeffs(self.Tbase,self.xbase)
        except (ValueError, AttributeError, TypeError, RuntimeError) as e:
            if self.viscosity.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(self.name,'viscosity',e))
            pass
        funcMu = None
        
        
        # Changed the coefficient order for TFreeze
        key = 'Tfreeze'
        def funcTf(T,x):
            x = x * 100.0
            a =  0.03422039835160944
            b = -0.05425629002714395
            c = -0.007991085555390726
            d =  0.001036937163846157
            e =  0.0003268582531827402
            f = -7.721483884155584E-06
            g = -4.841293026057464E-06
            h =  1.216929917247388E-08
            i =  2.41704396074865E-08
            j =  4.314861246570078E-11
            return ((a+x*(c+x*(e+x*(g+i*x))))/(1+x*(b+x*(d+x*(f+x*(h+j*x)))))*100)+273.15
        
        self.T_freeze.xData,self.T_freeze.yData,self.T_freeze.data = self.getArray(dataID=key,func=funcTf,x_in=np.array([0.0]),y_in=self.concentration.data,DEBUG=self.T_freeze.DEBUG)
        try:        
            self.T_freeze.coeffs = np.zeros((4,6))
            self.T_freeze.source = IncompressibleData.SOURCE_EQUATION
            self.T_freeze.type   = IncompressibleData.INCOMPRESSIBLE_POLYNOMIAL
            self.T_freeze.fitCoeffs(self.Tbase,self.xbase)
        except (ValueError, AttributeError, TypeError, RuntimeError) as e:
            if self.T_freeze.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(self.name,'T freeze',e))
            pass 
        funcTf = None
        
        
        

class AS10(PureData,DigitalData):
    def __init__(self):
        DigitalData.__init__(self)
        PureData.__init__(self)
        self.name = "AS10"
        self.description = "Aspen Temper -10, Potassium acetate/formate"
        self.reference = "Aspen Petroleum AB, SecCool software"

        self.Tmax = 30 + 273.15
        self.Tmin = -10 + 273.15
        self.TminPsat = self.Tmax
        self.Tbase = 0.00 + 273.15
        self.temperature.data         = self.getTrange()

        self.density.source = self.density.SOURCE_COEFFS
        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = np.array([[-1.89735e-18],[+ 1.66533e-16],[- 0.2],[+ 1090]])[::-1]

        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = np.array([[-0.0132589],[+ 2.01],[+ 3541.83]])[::-1]

        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = np.array([[0.001483],[+ 0.514259]])[::-1]

        self.viscosity.source = self.viscosity.SOURCE_COEFFS
        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_POLYNOMIAL
        self.viscosity.coeffs = np.array([[-2.11481e-5],[+ 0.00235381],[- 0.10631376],[+ 2.80154921]])[::-1]

    def fitFluid(self):
        pass




class AS20(PureData,DigitalData):
    def __init__(self):
        DigitalData.__init__(self)
        PureData.__init__(self)
        self.name = "AS20"
        self.description = "Aspen Temper -20, Potassium acetate/formate"
        self.reference = "Aspen Petroleum AB, SecCool software"

        self.Tmax = 30 + 273.15
        self.Tmin = -20 + 273.15
        self.TminPsat = self.Tmax
        self.Tbase = 0.00 + 273.15
        self.temperature.data         = self.getTrange()

        self.density.source = self.density.SOURCE_COEFFS
        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = np.array([[2.98156e-19],[- 0.00142857],[- 0.22142857],[+ 1147]])[::-1]

        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = np.array([[-0.0122673],[+ 2.86],[+ 3262.52]])[::-1]

        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = np.array([[0.001342],[+ 0.480766]])[::-1]
        

    def fitFluid(self):
        key = 'Mu'
        def funcMu(T,x):
            T = (T-self.Tbase)
            return 2.43708721027941*np.exp(-0.0537593944541809*T) + 0.97244

        self.viscosity.xData,self.viscosity.yData,self.viscosity.data = self.getArray(dataID=key,func=funcMu,x_in=self.temperature.data,y_in=self.concentration.data,DEBUG=self.viscosity.DEBUG)
        
        try:
            self.viscosity.source = self.viscosity.SOURCE_EQUATION
            self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
            self.viscosity.coeffs = np.zeros((4,6))
            self.viscosity.fitCoeffs(self.Tbase,self.xbase)
        except (ValueError, AttributeError, TypeError, RuntimeError) as e:
            if self.viscosity.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(self.name,'viscosity',e))
            pass
        funcMu = None




class AS30(PureData,DigitalData):
    def __init__(self):
        DigitalData.__init__(self)
        PureData.__init__(self)
        self.name = "AS30"
        self.description = "Aspen Temper -30, Potassium acetate/formate"
        self.reference = "Aspen Petroleum AB, SecCool software"

        self.Tmax = 30 + 273.15
        self.Tmin = -30 + 273.15
        self.TminPsat = self.Tmax
        self.Tbase = 0.00 + 273.15
        self.temperature.data         = self.getTrange()

        self.density.source = self.density.SOURCE_COEFFS
        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = np.array([[1.76768e-5],[- 0.00103896],[- 0.31085859],[+ 1183.85930736]])[::-1]

        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = np.array([[-0.0270232],[+ 2.99],[+ 3075.04]])[::-1]

        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = np.array([[0.001256],[+ 0.460388]])[::-1]


    def fitFluid(self):
        key = 'Mu'
        def funcMu(T,x):
            T = (T-self.Tbase)
            return 2.65653950695888*np.exp(-0.0598806339442954*T) + 1.30143

        self.viscosity.xData,self.viscosity.yData,self.viscosity.data = self.getArray(dataID=key,func=funcMu,x_in=self.temperature.data,y_in=self.concentration.data,DEBUG=self.viscosity.DEBUG)
        
        try:
            self.viscosity.source = self.viscosity.SOURCE_EQUATION
            self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
            self.viscosity.coeffs = np.zeros((4,6))
            self.viscosity.fitCoeffs(self.Tbase,self.xbase)
        except (ValueError, AttributeError, TypeError, RuntimeError) as e:
            if self.viscosity.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(self.name,'viscosity',e))
            pass
        funcMu = None




class AS40(PureData,DigitalData):
    def __init__(self):
        DigitalData.__init__(self)
        PureData.__init__(self)
        self.name = "AS40"
        self.description = "Aspen Temper -40, Potassium acetate/formate"
        self.reference = "Aspen Petroleum AB, SecCool software"

        self.Tmax = 30 + 273.15
        self.Tmin = -40 + 273.15
        self.TminPsat = self.Tmax
        self.Tbase = 0.00 + 273.15
        self.temperature.data         = self.getTrange()

        self.density.source = self.density.SOURCE_COEFFS
        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = np.array([[1.6835e-5],[- 0.00109307],[- 0.37819865],[+ 1214.83982684]])[::-1]

        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = np.array([[-0.0387227],[+ 2.28],[+ 2977.88]])[::-1]

        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = np.array([[0.001099],[+ 0.443327]])[::-1]

    def fitFluid(self):
        key = 'Mu'
        def funcMu(T,x):
            T = (T-self.Tbase)
            return 0.714976365635003*np.exp(-0.100050525515385*T) + 4.38768154440393*np.exp(-0.0260039000649317*T)

        self.viscosity.xData,self.viscosity.yData,self.viscosity.data = self.getArray(dataID=key,func=funcMu,x_in=self.temperature.data,y_in=self.concentration.data,DEBUG=self.viscosity.DEBUG)
        
        try:
            self.viscosity.source = self.viscosity.SOURCE_EQUATION
            self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
            self.viscosity.coeffs = np.zeros((4,6))
            self.viscosity.fitCoeffs(self.Tbase,self.xbase)
        except (ValueError, AttributeError, TypeError, RuntimeError) as e:
            if self.viscosity.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(self.name,'viscosity',e))
            pass
        funcMu = None




class AS55(PureData,DigitalData):
    def __init__(self):
        DigitalData.__init__(self)
        PureData.__init__(self)
        self.name = "AS55"
        self.description = "Aspen Temper -55, Potassium acetate/formate"
        self.reference = "Aspen Petroleum AB, SecCool software"

        self.Tmax = 30 + 273.15
        self.Tmin = -55 + 273.15
        self.TminPsat = self.Tmax
        self.Tbase = 0.00 + 273.15
        self.temperature.data         = self.getTrange()

        self.density.source = self.density.SOURCE_COEFFS
        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = np.array([[1.98824e-5],[- 0.00117189],[- 0.47629615],[+ 1249.7534665]])[::-1]

        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = np.array([[-0.0248618],[+ 2.29],[+ 2839.85]])[::-1]

        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = np.array([[2.287e-6],[+ 0.000937902],[+ 0.425799423]])[::-1]


    def fitFluid(self):
        key = 'Mu'
        def funcMu(T,x):
            T = (T-self.Tbase)
            return 0.159583223482554*np.exp(-0.138097704125669*T) + 6.3176967296442*np.exp(-0.0380509974688477*T)

        self.viscosity.xData,self.viscosity.yData,self.viscosity.data = self.getArray(dataID=key,func=funcMu,x_in=self.temperature.data,y_in=self.concentration.data,DEBUG=self.viscosity.DEBUG)
        
        try:
            self.viscosity.source = self.viscosity.SOURCE_EQUATION
            self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
            self.viscosity.coeffs = np.zeros((4,6))
            self.viscosity.fitCoeffs(self.Tbase,self.xbase)
        except (ValueError, AttributeError, TypeError, RuntimeError) as e:
            if self.viscosity.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(self.name,'viscosity',e))
            pass
        funcMu = None



    
    
if __name__ == '__main__':   
    
    print("The main method converts the AspenTemper coefficients from SecCool")
    names = []
    descs = []
    refs  = []
    
    Tmax = []
    Tmin = []
    TminPsat = []
    Tbase = []
    
    cond = []
    cp   = []
    mu   = []
    rho  = []
    
    ID = 10
    names.append("AS{0}".format(ID))
    descs.append("Aspen Temper -{0}, Potassium acetate/formate".format(ID))
    refs.append ("Aspen Petroleum AB, SecCool software")
    
    Tmax.append("30 + 273.15")
    Tmin.append("-{0}".format(ID))
    TminPsat.append("self.Tmax")
    Tbase.append("0.00 + 273.15")
    
    cond.append ("0.001483*(273+T)+0.1094")
    cp.append   ("1000*(3.54183+T*(0.00201-0.0000132589*T))")
    mu.append   ("2.80154921+T*(-0.10631376+T*(0.00235381-0.0000211481*T))")
    rho.append  ("1090+T*(-0.2+T*(1.66533E-16-1.89735E-18*T))")
    
    ID = 20
    names.append("AS{0}".format(ID))
    descs.append("Aspen Temper -{0}, Potassium acetate/formate".format(ID))
    refs.append ("Aspen Petroleum AB, SecCool software")
        
    Tmax.append("30 + 273.15")
    Tmin.append("-{0}".format(ID))
    TminPsat.append("self.Tmax")
    Tbase.append("0.00 + 273.15")
    
    cond.append ("0.001342*(273+T)+0.1144")
    cp.append   ("1000*(3.26252+T*(0.00286-0.0000122673*T))")
    mu.append   ("(0.97244+(7.14199*exp((-20-T)/18.6014)))")
    rho.append  ("1147+T*(-0.22142857+T*(-0.00142857+2.98156E-19*T))")


    ID = 30
    names.append("AS{0}".format(ID))
    descs.append("Aspen Temper -{0}, Potassium acetate/formate".format(ID))
    refs.append ("Aspen Petroleum AB, SecCool software")
        
    Tmax.append("30 + 273.15")
    Tmin.append("-{0}".format(ID))
    TminPsat.append("self.Tmax")
    Tbase.append("0.00 + 273.15")
    
    cond.append ("0.001256*(273+T)+0.1175")
    cp.append   ("1000*(3.07504+T*(0.00299-0.0000270232*T))")
    mu.append   ("(1.30143+(16.01368*exp((-30-T)/16.69989)))")
    rho.append  ("1183.85930736+T*(-0.31085859+T*(-0.00103896+0.0000176768*T))")


    ID = 40
    names.append("AS{0}".format(ID))
    descs.append("Aspen Temper -{0}, Potassium acetate/formate".format(ID))
    refs.append ("Aspen Petroleum AB, SecCool software")
        
    Tmax.append("30 + 273.15")
    Tmin.append("-{0}".format(ID))
    
    cond.append ("0.001099*(273+T)+0.1433")
    cp.append   ("1000*(2.97788+T*(0.00228-0.0000387227*T))")
    mu.append   ("(39.11536*exp((-40-T)/9.99495)+(12.41564*exp((-40-T)/38.45577)))")
    rho.append  ("1214.83982684+T*(-0.37819865+T*(-0.00109307+0.000016835*T))")

    ID = 55
    names.append("AS{0}".format(ID))
    descs.append("Aspen Temper -{0}, Potassium acetate/formate".format(ID))
    refs.append ("Aspen Petroleum AB, SecCool software")
        
    Tmax.append("30 + 273.15")
    Tmin.append("-{0}".format(ID))
    TminPsat.append("self.Tmax")
    Tbase.append("0.00 + 273.15")
    
    cond.append ("(273+T)*(0.000002287*(273+T)-0.0003108)+0.3402")
    cp.append   ("1000*(2.83985+T*(0.00229-0.0000248618*T))")
    mu.append   ("(317.40673*exp((-55-T)/7.24125)+(51.22151*exp((-55-T)/26.28052)))")
    rho.append  ("1249.7534665+T*(-0.47629615+T*(-0.00117189+0.0000198824*T))")
    
    spa = "        "
    
    from sympy import symbols, expand
    x, T = symbols('x T') 
    
    for i in range(len(names)):
        print("class {0}(PureData):\n    def __init__(self):".format(names[i]))
        print("{0}PureData.__init__(self)".format(spa))
        print("{0}self.name = \"{1}\"".format(spa,names[i]))
        print("{0}self.description = \"{1}\"".format(spa,descs[i]))
        print("{0}self.reference = \"{1}\"".format(spa,refs[i]))
        print("")
        print("{0}self.Tmax = {1}".format(spa,Tmax[i]))
        print("{0}self.Tmin = {1}".format(spa,Tmin[i]))
        print("{0}self.TminPsat = self.Tmax".format(spa))
        print("{0}self.Tbase = 0.00 + 273.15".format(spa))
        print("")
        print("{0}self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL".format(spa))
        print("{0}self.density.coeffs = np.array([{1}])[::-1]".format(spa,expand(rho[i])))
        print("")
        print("{0}self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL".format(spa))
        print("{0}self.specific_heat.coeffs = np.array([{1}])[::-1]".format(spa,expand(cp[i])))
        print("")
        print("{0}self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL".format(spa))
        print("{0}self.conductivity.coeffs = np.array([{1}])[::-1]".format(spa,expand(cond[i])))
        print("")
        print("{0}self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL".format(spa))
        print("{0}self.viscosity.coeffs = np.array([{1}])[::-1]".format(spa,expand(mu[i])))
        print("")
        print("{0}def fitFluid(self):".format("    "))
        print("{0}key = \"Mu\"".format(spa))
        print("{0}def funcMu(T,x):".format(spa))
        print("")
        print("")
        print("{0}funcMu=None".format(spa))
        print("")
        print("")
        print("")
        print("")
        
#        PureData.__init__(self)
#        self.name = "AS10"
#        self.description = "Aspen Temper -10, Potassium acetate/formate"
#        self.reference = "Aspen Petroleum AB, SecCool software"
#        
#        self.Tmax =  30 + 273.15
#        self.Tmin = -10 + 273.15
#        self.TminPsat =  self.Tmax
#    
#        self.Tbase =   0.00 + 273.15
#       
#        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
#        self.density.coeffs = np.array([[945.5454545],[-1.054545455]])
#        
#        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
#        self.specific_heat.coeffs = np.array([[2.322218182],[0.003843636]])*1000
#
#        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
#        self.conductivity.coeffs = np.array([[0.001483*273.15+0.1094],[0.001483]])
#        
#        self.temperature.data         = self.getTrange()
#        self.concentration.data       = np.array([0]) # mass fraction


#
#
#{ THyCool_20 }
#
#constructor THyCool_20.Create;
#begin
#  inherited;
#  TFreeze   := -20;
#  FluidType := 'Potassium Formate';
#  TradeName := 'HYCOOL 20';
#  Reference := 'Hydro Chemicals';
#  TMin      := -20;
#  TMax      := 50;
#end;
#
#function THyCool_20.GetCondT(T: Double): Double;
#begin
#  if T <= 20 then
#    Result := 0.001978*T+0.5172
#  else
#    Result := 0.001005*T+0.5368;
#end;
#
#function THyCool_20.GetCpT(T: Double): Double;
#begin
#  Result := 1000*(0.0023*T+2.955);
#end;
#
#function THyCool_20.GetMuT(T: Double): Double;
#begin
#  if T <= 20 then
#    Result := 0.07190*exp(524.75/(T+142.05))
#  else
#    Result := T*(0.0005524*T - 0.06281)+2.8536;
#end;
#
#function THyCool_20.GetRhoT(T: Double): Double;
#begin
#  Result := -0.429180*T+1202.2;
#end;
#
#{ THyCool_30 }
#
#constructor THyCool_30.Create;
#begin
#  inherited;
#  TFreeze   := -30;
#  FluidType := 'Potassium Formate';
#  TradeName := 'HYCOOL 30';
#  Reference := 'Hydro Chemicals';
#  TMin      := -30;
#  TMax      := 50;
#end;
#
#function THyCool_30.GetCondT(T: Double): Double;
#begin
#  if T <= 20 then
#    Result := 0.001840*T+0.4980
#  else
#    Result := 0.001000*T+0.514;
#end;
#
#function THyCool_30.GetCpT(T: Double): Double;
#begin
#  Result := 1000*(0.0023*T+2.783);
#end;
#
#function THyCool_30.GetMuT(T: Double): Double;
#begin
#  if T <= 20 then
#    Result := 0.11100*exp(408.17/(T+123.15))
#  else
#    Result := T*(0.000295*T - 0.0441)+2.6836;
#end;
#
#function THyCool_30.GetRhoT(T: Double): Double;
#begin
#  Result := -0.475350*T+1257.5;
#end;
#
#{ THyCool_40 }
#
#constructor THyCool_40.Create;
#begin
#  inherited;
#  TFreeze   := -40;
#  FluidType := 'Potassium Formate';
#  TradeName := 'HYCOOL 40';
#  Reference := 'Hydro Chemicals';
#  TMin      := -40;
#  TMax      := 20;
#end;
#
#function THyCool_40.GetCondT(T: Double): Double;
#begin
#  Result := 0.001730*T+0.4826;
#end;
#
#function THyCool_40.GetCpT(T: Double): Double;
#begin
#  Result := 1000*(0.0023*T+2.646);
#end;
#
#function THyCool_40.GetMuT(T: Double): Double;
#begin
#  Result := 0.07830*exp(498.13/(T+130.25));
#end;
#
#function THyCool_40.GetRhoT(T: Double): Double;
#begin
#  Result := -0.512290*T+1304.5;
#end;
#
#{ THyCool_45 }
#
#constructor THyCool_45.Create;
#begin
#  inherited;
#  TFreeze   := -45;
#  FluidType := 'Potassium Formate';
#  TradeName := 'HYCOOL 45';
#  Reference := 'Hydro Chemicals';
#  TMin      := -40;
#  TMax      := 20;
#end;
#
#function THyCool_45.GetCondT(T: Double): Double;
#begin
#  Result := 0.001674*T+0.4750;
#end;
#
#function THyCool_45.GetCpT(T: Double): Double;
#begin
#  Result := 1000*(0.0023*T+2.578);
#end;
#
#function THyCool_45.GetMuT(T: Double): Double;
#begin
#  Result := 0.08990*exp(479.09/(T+126.55));
#end;
#
#function THyCool_45.GetRhoT(T: Double): Double;
#begin
#  Result := -0.530754*T+1328.7;
#end;
#
#{ THyCool_50 }
#
#constructor THyCool_50.Create;
#begin
#  inherited;
#  TFreeze   := -50;
#  FluidType := 'Potassium Formate';
#  TradeName := 'HYCOOL 50';
#  Reference := 'Hydro Chemicals';
#  TMin      := -50;
#  TMax      := 20;
#end;
#
#function THyCool_50.GetCondT(T: Double): Double;
#begin
#  Result := 0.001610*T+0.4660;
#end;
#
#function THyCool_50.GetCpT(T: Double): Double;
#begin
#  Result := 1000*(0.0023*T+2.498);
#end;
#
#function THyCool_50.GetMuT(T: Double): Double;
#begin
#  Result := 0.0491*exp(581.12/(T+129.05));
#  if T > -10 then
#    Result := Result+0.2;
#end;
#
#function THyCool_50.GetRhoT(T: Double): Double;
#begin
#  Result := -0.552300*T+1359.0;
#end;
#
#        
#        
#        
#        
#        
#        
#        
#        
#        
#        
#        
#        
#        
        
        
        
        
        
        
        
        
        
        
        
        

