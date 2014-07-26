from __future__ import division, print_function
import numpy as np

import hashlib, os, json

class SolutionDataWriter(object):
    """ 
    A base class that defines all the variables needed 
    in order to make a proper fit. You can copy this code
    put in your data and add some documentation for where the
    information came from. 
    """
    def __init__(self):
        pass 
    
    def fitAll(self, fluidObject):
        
        tempData = fluidObject.temperature.data
        concData = fluidObject.concentration.data
        
        if fluidObject.Tbase==0.0 or fluidObject.Tbase==None:
            fluidObject.Tbase = (np.min(tempData) + np.max(tempData)) / 2.0
        if fluidObject.xbase==0.0 or fluidObject.xbase==None:
            fluidObject.xbase = (np.min(concData) + np.max(concData)) / 2.0
            
        # Set the standard order for polynomials
        std_xorder = 3+1
        std_yorder = 5+1
        std_coeffs = np.zeros((std_xorder,std_yorder))
        
        errList = (ValueError, AttributeError, TypeError, RuntimeError)
        
        try:
            fluidObject.density.coeffs = np.copy(std_coeffs)
            fluidObject.density.type   = fluidObject.density.INCOMPRESSIBLE_POLYNOMIAL
            fluidObject.density.fitCoeffs(tempData,concData,fluidObject.Tbase,fluidObject.xbase)
        except errList as ve:
            if fluidObject.density.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(fluidObject.name,'density',ve))
            pass
        
        try:
            fluidObject.specific_heat.coeffs = np.copy(std_coeffs)
            fluidObject.specific_heat.type   = fluidObject.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
            fluidObject.specific_heat.fitCoeffs(tempData,concData,fluidObject.Tbase,fluidObject.xbase)
        except errList as ve:
            if fluidObject.specific_heat.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(fluidObject.name,'specific heat',ve))
            pass 
        
        try:
            fluidObject.conductivity.coeffs = np.copy(std_coeffs)
            fluidObject.conductivity.type   = fluidObject.conductivity.INCOMPRESSIBLE_POLYNOMIAL
            fluidObject.conductivity.fitCoeffs(tempData,concData,fluidObject.Tbase,fluidObject.xbase)
        except errList as ve:
            if fluidObject.conductivity.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(fluidObject.name,'conductivity',ve))
            pass
        
        try:
            fluidObject.viscosity.coeffs = np.copy(std_coeffs)
            fluidObject.viscosity.type   = fluidObject.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
            fluidObject.viscosity.fitCoeffs(tempData,concData,fluidObject.Tbase,fluidObject.xbase)
        except errList as ve:
            if fluidObject.viscosity.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(fluidObject.name,'viscosity',ve))
            pass
        
        # reset data for getArray and read special files
        if fluidObject.xid!=fluidObject.ifrac_pure and fluidObject.xid!=fluidObject.ifrac_undefined:
            try:        
                fluidObject.T_freeze.coeffs = np.copy(std_coeffs)
                fluidObject.T_freeze.type   = fluidObject.T_freeze.INCOMPRESSIBLE_POLYNOMIAL
                fluidObject.T_freeze.fitCoeffs(concData,0.0,fluidObject.xbase,0.0)
            except errList as ve:
                if fluidObject.T_freeze.DEBUG: print("{0}: Could not fit {1} coefficients: {2}".format(fluidObject.name,"T_freeze",ve))
                pass
# 
#            # reset data for getArray again
#            if fluidObject.xid==fluidObject.ifrac_volume:
#                try:
#                    fluidObject.mass2input.coeffs = np.copy(std_coeffs)
#                    fluidObject.mass2input.type   = fluidObject.mass2input.INCOMPRESSIBLE_POLYNOMIAL
#                    fluidObject.mass2input.fitCoeffs([fluidObject.Tbase],massData,fluidObject.Tbase,fluidObject.xbase)
#                except errList as ve:
#                    if fluidObject.mass2input.DEBUG: print("{0}: Could not fit {1} coefficients: {2}".format(fluidObject.name,"mass2input",ve))
#                    pass
#            elif fluidObject.xid==fluidObject.ifrac_mass:
#                _,_,fluidObject.volume2input.data = IncompressibleData.shfluidObject.ray(massData,axs=1)
#                #_,_,volData = IncompressibleData.shapeArray(volData,axs=1)
#                try:
#                    fluidObject.volume2input.coeffs = np.copy(std_coeffs)
#                    fluidObject.volume2input.type   = fluidObject.volume2input.INCOMPRESSIBLE_POLYNOMIAL
#                    fluidObject.volume2input.fitCoeffs([fluidObject.Tbase],volData,fluidObject.Tbase,fluidObject.xbase)
#                except errList as ve:
#                    if fluidObject.volume2input.DEBUG: print("{0}: Could not fit {1} coefficients: {2}".format(fluidObject.name,"volume2input",ve))
#                    pass
#            else:
#                raise ValueError("Unknown xid specified.")


    def get_hash(self,data):
        return hashlib.sha224(data).hexdigest()
    
    def get_hash_file(self):
        return os.path.join(os.path.dirname(__file__), 'data', "hashes.json")
    
    def load_hashes(self):
        hashes_fname = self.get_hash_file()
        if os.path.exists(hashes_fname):
            hashes = json.load(open(hashes_fname,'r'))
        else:
            hashes = dict()
        return hashes
    
    def write_hashes(self, hashes):
        hashes_fname = self.get_hash_file()
        fp = open(hashes_fname,'w')
        fp.write(json.dumps(hashes))
        fp.close()
        return True
            
    
    def toJSON(self,data,quiet=False):
        jobj = {}
        
        jobj['name'] = data.name # Name of the current fluid
        jobj['description'] = data.description # Description of the current fluid
        jobj['reference'] = data.reference # Reference data for the current fluid
        
        jobj['Tmax'] = data.Tmax         # Maximum temperature in K
        jobj['Tmin'] = data.Tmin         # Minimum temperature in K
        jobj['xmax'] = data.xmax         # Maximum concentration
        jobj['xmin'] = data.xmin         # Minimum concentration
        jobj['xid'] = data.xid           # Concentration is mole, mass or volume-based
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
        jobj['mass2input'] = data.mass2input.toJSON() # dd
        jobj['volume2input'] = data.volume2input.toJSON() # dd
        jobj['mole2input'] = data.mole2input.toJSON() # dd
               
        original_float_repr = json.encoder.FLOAT_REPR 
        #print json.dumps(1.0001) 
        stdFmt = " .{0}e".format(int(data.significantDigits-1))
        #pr = np.finfo(float).eps * 10.0
        #pr = np.finfo(float).precision - 2 # stay away from numerical precision
        #json.encoder.FLOAT_REPR = lambda o: format(np.around(o,decimals=pr), stdFmt) 
        json.encoder.FLOAT_REPR = lambda o: format(o, stdFmt)
        dump = json.dumps(jobj, indent = 2, sort_keys = True)
        json.encoder.FLOAT_REPR = original_float_repr
        
        #print dump
        hashes = self.load_hashes()
        hash   = self.get_hash(dump)
        
        name = jobj['name']
               
        if name not in hashes or \
          hashes[name] != hash: # update hashes and write file
            
            hashes[name] = hash
            self.write_hashes(hashes)           
            
            fp = open(name+'.json', 'w')
            fp.write(dump)
            fp.close()
            
            if not quiet: print(" ({0})".format("w"), end="")
        else:
            if not quiet: print(" ({0})".format("i"), end="")
                
        
    def printStatusID(self, fluidObjs, obj): 
        #obj = fluidObjs[num]
        if obj==fluidObjs[0]:
            print(" {0}".format(obj.name), end="")
        elif obj==fluidObjs[-1]:
            print(", {0}".format(obj.name), end="")
        else:
            print(", {0}".format(obj.name), end="")

        return


    def fitFluidList(self, fluidObjs):
        print("Fitting normal fluids:", end="")
        for obj in fluidObjs:
            self.printStatusID(fluidObjs, obj) 
            try: 
                self.fitAll(obj)
            except (TypeError, ValueError) as e:
                print("An error occurred for fluid: {0}".format(obj.name))
                print(obj)
                print(e)
                pass 
        print(" ... done")
        return
    
    
    def fitSecCoolList(self, fluidObjs):
        print("Fitting SecCool fluids:", end="")
        for obj in fluidObjs:
            self.printStatusID(fluidObjs, obj) 
            try: 
                obj.fitFluid()
            except (TypeError, ValueError) as e:
                print("An error occurred for fluid: {0}".format(obj.name))
                print(obj)
                print(e)
                pass 
        print(" ... done")
        return
    
    
    def writeFluidList(self, fluidObjs):
        print("Legend: FluidName (w) | (i) -> (w)=written, (i)=ignored, unchanged coefficients")
        print("Writing fluids to JSON:", end="")
        for obj in fluidObjs:
            self.printStatusID(fluidObjs, obj)                 
            try: 
                self.toJSON(obj)
            except (TypeError, ValueError) as e:
                print("An error occurred for fluid: {0}".format(obj.name))
                print(obj)
                print(e)
                pass
        print(" ... done") 
        return 
        
        
#class FitGraphWriter(object):
#    """ 
#    A base class that defines all the variables needed 
#    in order to make a proper fit. You can copy this code
#    put in your data and add some documentation for where the
#    information came from. 
#    """
#    def __init__(self):
#        self.verbose = True
#        
#    def inCoolProp(self,name):
#        #print FluidsList() 
#        result = name in FluidsList()
#        if not result:
#            try:
#                CP.PropsU('Tmin','T',0,'P',0,name,"SI")
#                return True
#            except ValueError as e:
#                print e
#                return False
#            
#    def getFluidList(self):
#        containerList = []
##        containerList += [TherminolD12()]
##        containerList += [TherminolVP1(), Therminol66(), Therminol72()]
##        containerList += [DowthermJ(), DowthermQ()]
##        containerList += [Texatherm22(),  NitrateSalt(), SylthermXLT()]
##        containerList += [HC50(), HC40(), HC30(), HC20(), HC10()]
##        containerList += [AS10(), AS20(), AS30(), AS40(), AS55()]
##        containerList += [ZS10(), ZS25(), ZS40(), ZS45(), ZS55()]
#        return containerList
#    
#    def relError(A=[],B=[],PCT=False):
#        """
#        Returns the relative Error from (A-B)/B, if PCT is True, it returns percent.
#        """
#        result = (np.array(A)-np.array(B))/np.array(B);
#        if PCT:
#            return result * 100. 
#        else:
#            return result
#        
#        
#    def makePlots(self, fluid):
#        # row and column sharing for test plots
#        #matplotlib.pyplot.subplots_adjust(top=0.85)
#        f, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = matplotlib.pyplot.subplots(3, 2, sharex='col')
#        f.set_size_inches(matplotlib.pyplot.figaspect(1.2)*1.5)
#        #f.suptitle("Fit for "+str(data.Desc), fontsize=14)
#        
##        ### This is the actual fitting
#        tData = data.T
#        tDat1 = numpy.linspace(numpy.min(tData)+1, numpy.max(tData)-1, 10)
#        Pin = 1e20 # Dummy pressure
#        inCP =liqObj.inCoolProp(data.Name)
#        print "Fluid in CoolProp: "+str(inCP)
#        print
        
    
         
