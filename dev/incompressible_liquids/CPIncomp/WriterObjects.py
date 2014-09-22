from __future__ import division, print_function
import numpy as np

import hashlib, os, json, sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from CPIncomp.DataObjects import SolutionData
from CPIncomp.BaseObjects import IncompressibleData, IncompressibleFitter
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator
from matplotlib.backends.backend_pdf import PdfPages
import itertools
import matplotlib
import csv
from CoolProp.BibtexParser import BibTeXerClass
from warnings import warn

class SolutionDataWriter(object):
    """
    A base class that defines all the variables needed
    in order to make a proper fit. You can copy this code
    put in your data and add some documentation for where the
    information came from.
    """
    def __init__(self):
        self.bibtexer = BibTeXerClass('../../Web/fluid_properties/Incompressibles.bib')

    def fitAll(self, fluidObject=SolutionData()):

        if fluidObject.Tbase is None:
            fluidObject.Tbase = (fluidObject.Tmin + fluidObject.Tmax) / 2.0

        if fluidObject.xbase is None:
            fluidObject.xbase = (fluidObject.xmin + fluidObject.xmax) / 2.0

        tData = fluidObject.temperature.data
        xData = fluidObject.concentration.data
        tBase = fluidObject.Tbase
        xBase = fluidObject.xbase

        # Set the standard order for polynomials
        std_xorder = 3+1
        std_yorder = 5+1
        std_coeffs = np.zeros((std_xorder,std_yorder))

        errList = (ValueError, AttributeError, TypeError, RuntimeError)

        if fluidObject.density.coeffs is None:
            try:
                fluidObject.density.setxyData(tData,xData)
                fluidObject.density.coeffs = np.copy(std_coeffs)
                fluidObject.density.type   = IncompressibleData.INCOMPRESSIBLE_POLYNOMIAL
                fluidObject.density.fitCoeffs(tBase,xBase)
            except errList as ve:
                if fluidObject.density.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(fluidObject.name,'density',ve))
                pass

        if fluidObject.specific_heat.coeffs is None:
            try:
                fluidObject.specific_heat.setxyData(tData,xData)
                fluidObject.specific_heat.coeffs = np.copy(std_coeffs)
                fluidObject.specific_heat.type   = IncompressibleData.INCOMPRESSIBLE_POLYNOMIAL
                fluidObject.specific_heat.fitCoeffs(tBase,xBase)
            except errList as ve:
                if fluidObject.specific_heat.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(fluidObject.name,'specific heat',ve))
                pass

        if fluidObject.conductivity.coeffs is None:
            try:
                fluidObject.conductivity.setxyData(tData,xData)
                fluidObject.conductivity.coeffs = np.copy(std_coeffs)
                fluidObject.conductivity.type   = IncompressibleData.INCOMPRESSIBLE_POLYNOMIAL
                fluidObject.conductivity.fitCoeffs(tBase,xBase)
            except errList as ve:
                if fluidObject.conductivity.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(fluidObject.name,'conductivity',ve))
                pass

        if fluidObject.viscosity.coeffs is None:
            try:
                fluidObject.viscosity.setxyData(tData,xData)
                tried = False
                if len(fluidObject.viscosity.yData)==1:# and np.isfinite(fluidObject.viscosity.data).sum()<10:
                    fluidObject.viscosity.coeffs = np.array([+5e+2, -6e+1, +1e+1])
                    fluidObject.viscosity.type   = IncompressibleData.INCOMPRESSIBLE_EXPONENTIAL
                    fluidObject.viscosity.fitCoeffs(tBase,xBase)
                    if fluidObject.viscosity.coeffs is None or IncompressibleFitter.allClose(fluidObject.viscosity.coeffs, np.array([+5e+2, -6e+1, +1e+1])): # Fit failed
                        tried = True
                if len(fluidObject.viscosity.yData)>1 or tried:
                    #fluidObject.viscosity.coeffs = np.zeros(np.round(np.array(std_coeffs.shape) * 1.5))
                    fluidObject.viscosity.coeffs = np.copy(std_coeffs)
                    fluidObject.viscosity.type   = IncompressibleData.INCOMPRESSIBLE_EXPPOLYNOMIAL
                    fluidObject.viscosity.fitCoeffs(tBase,xBase)
            except errList as ve:
                if fluidObject.viscosity.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(fluidObject.name,'viscosity',ve))
                pass

        if fluidObject.saturation_pressure.coeffs is None:
            try:
                fluidObject.saturation_pressure.setxyData(tData,xData)
                tried = False
                if len(fluidObject.saturation_pressure.yData)==1:# and np.isfinite(fluidObject.saturation_pressure.data).sum()<10:
                    fluidObject.saturation_pressure.coeffs = np.array([-5e+3, +6e+1, -1e+1])
                    fluidObject.saturation_pressure.type   = IncompressibleData.INCOMPRESSIBLE_EXPONENTIAL
                    fluidObject.saturation_pressure.fitCoeffs(tBase,xBase)
                    if fluidObject.saturation_pressure.coeffs is None or IncompressibleFitter.allClose(fluidObject.saturation_pressure.coeffs, np.array([-5e+3, +6e+1, -1e+1])): # Fit failed
                        tried = True
                if len(fluidObject.saturation_pressure.yData)>1 or tried:
                    #fluidObject.saturation_pressure.coeffs = np.zeros(np.round(np.array(std_coeffs.shape) * 1.5))
                    fluidObject.saturation_pressure.coeffs = np.copy(std_coeffs)
                    fluidObject.saturation_pressure.type   = IncompressibleData.INCOMPRESSIBLE_EXPPOLYNOMIAL
                    fluidObject.saturation_pressure.fitCoeffs(tBase,xBase)
            except errList as ve:
                if fluidObject.saturation_pressure.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(fluidObject.name,'saturation pressure',ve))
                pass

        # reset data for getArray and read special files
        if fluidObject.xid!=fluidObject.ifrac_pure and fluidObject.xid!=fluidObject.ifrac_undefined:
            if fluidObject.T_freeze.coeffs is None:
                fluidObject.T_freeze.setxyData([0.0],xData)
                try:
                    if len(fluidObject.T_freeze.xData)==1:# and np.isfinite(fluidObject.T_freeze.data).sum()<10:
                        fluidObject.T_freeze.coeffs = np.array([+7e+2, -6e+1, +1e+1])
                        fluidObject.T_freeze.type   = IncompressibleData.INCOMPRESSIBLE_EXPONENTIAL
                    else:
                        fluidObject.specific_heat.coeffs = np.copy(std_coeffs)
                        fluidObject.T_freeze.type   = IncompressibleData.INCOMPRESSIBLE_EXPPOLYNOMIAL
                    fluidObject.T_freeze.fitCoeffs(tBase,xBase)
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


    def get_json_file(self,name):
        return os.path.join("json","{0}.json".format(name))

    def get_report_file(self,name):
        return os.path.join("report","{0}_fitreport.pdf".format(name))


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
        stdFmt = "1.{0}e".format(int(data.significantDigits-1))
        #pr = np.finfo(float64).eps * 10.0
        #pr = np.finfo(float64).precision - 2 # stay away from numerical precision
        #json.encoder.FLOAT_REPR = lambda o: format(np.around(o,decimals=pr), stdFmt)
        json.encoder.FLOAT_REPR = lambda o: format(o, stdFmt)
        dump = json.dumps(jobj, indent = 2, sort_keys = True)
        json.encoder.FLOAT_REPR = original_float_repr

        #print dump
        hashes = self.load_hashes()
        hash   = self.get_hash(dump)

        name = jobj['name']

        if name not in hashes or \
          hashes[name] != hash or \
          not os.path.isfile(self.get_json_file(name)): # update hashes and write file

            hashes[name] = hash
            self.write_hashes(hashes)
            path = self.get_json_file(name)
            if not os.path.exists(os.path.dirname(path)):
                os.makedirs(os.path.dirname(path))
            with open(path, 'w') as fp:
                fp.write(dump)

            if not quiet: print(" ({0})".format("w"), end="")
        else:
            if not quiet: print(" ({0})".format("i"), end="")


        # Update the object:
        self.fromJSON(data=data)


    def fromJSON(self,data=SolutionData()):

        path = os.path.join("json",data.name+'.json')

        if not os.path.isfile(path):
            return None

        with open(path) as json_file:
            jobj = json.load(json_file)

            data.name        = jobj['name']  # Name of the current fluid
            data.description = jobj['description']  # Description of the current fluid
            data.reference   = jobj['reference']  # Reference data for the current fluid

            data.Tmax = jobj['Tmax']          # Maximum temperature in K
            data.Tmin = jobj['Tmin']          # Minimum temperature in K
            data.xmax = jobj['xmax']          # Maximum concentration
            data.xmin = jobj['xmin']          # Minimum concentration
            data.xid = jobj['xid']            # Concentration is mole, mass or volume-based
            data.TminPsat = jobj['TminPsat']      # Minimum saturation temperature in K
            data.Tbase = jobj['Tbase']        # Base value for temperature fits
            data.xbase = jobj['xbase']        # Base value for concentration fits

            #data.temperature  # Temperature for data points in K
            #data.concentration  # Concentration data points in weight fraction
            data.density.fromJSON(jobj['density']) # Density in kg/m3
            data.specific_heat.fromJSON(jobj['specific_heat']) # Heat capacity in J/(kg.K)
            data.viscosity.fromJSON(jobj['viscosity']) # Dynamic viscosity in Pa.s
            data.conductivity.fromJSON(jobj['conductivity']) # Thermal conductivity in W/(m.K)
            data.saturation_pressure.fromJSON(jobj['saturation_pressure']) # Saturation pressure in Pa
            data.T_freeze.fromJSON(jobj['T_freeze']) # Freezing temperature in K
            data.mass2input.fromJSON(jobj['mass2input']) # dd
            data.volume2input.fromJSON(jobj['volume2input']) # dd
            data.mole2input.fromJSON(jobj['mole2input']) # dd

            return data


    def printStatusID(self, fluidObjs, obj):
        #obj = fluidObjs[num]
        if obj==fluidObjs[0]:
            print(" {0}".format(obj.name), end="")
        elif obj==fluidObjs[-1]:
            print(", {0}".format(obj.name), end="")
        else:
            print(", {0}".format(obj.name), end="")

        sys.stdout.flush()
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


    def readFluidList(self, fluidObjs):
        print("Reading fluids:", end="")
        for obj in fluidObjs:
            self.printStatusID(fluidObjs, obj)
            try:
                self.fromJSON(obj)
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

    def writeReportList(self, fluidObjs, pdfFile=None):
        print("Writing fitting reports:", end="")

        if not pdfFile is None:
            if not os.path.exists(os.path.dirname(pdfFile)):
                os.makedirs(os.path.dirname(pdfFile))
            with PdfPages(pdfFile) as pdfObj:
                for obj in fluidObjs:
                    matplotlib.pyplot.close("all")
                    self.printStatusID(fluidObjs, obj)
                    try:
                        self.makeFitReportPage(obj,pdfObj=pdfObj)
                    except (TypeError, ValueError) as e:
                        print("An error occurred for fluid: {0}".format(obj.name))
                        print(obj)
                        print(e)
                        pass
        else:
            for obj in fluidObjs:
                matplotlib.pyplot.close("all")
                self.printStatusID(fluidObjs, obj)
                try:
                    self.makeFitReportPage(obj)
                except (TypeError, ValueError) as e:
                    print("An error occurred for fluid: {0}".format(obj.name))
                    print(obj)
                    print(e)
                    pass

        print(" ... done")
        return


    #####################################
    # Plotting routines
    #####################################
    def relError(self, A=[], B=[], PCT=False):
        """
        Returns the absolute relative Error from either
        (B-A)/B or (A-B)/A for abs(B)>abs(A) and abs(A)>abs(B),
        respectively. If PCT is True, it returns it in percent.
        """
        A_a = np.array(A)
        B_a = np.array(B)

        abl = np.absolute(B_a)
        eps = np.ones_like(abl) * np.finfo(float).eps
        #div = np.amax(np.hstack((abl,eps)), axis=1)*np.sign(B_a)

        pos = np.isfinite(B_a)
        pos2 = (B_a>eps)
        result = np.ones_like(A_a)*np.NAN

        result[pos & pos2] = (A_a[pos & pos2]-B_a[pos & pos2])/B_a[pos & pos2]

        if PCT:
            return result * 100.
        else:
            return result

    ############################################################
    # Define the general purpose routines for plotting
    ############################################################
    def wireFrame2D(self,xz,yz,linesX=5,linesY=None,color='black',ax=None,plot=False):
        """
        xz is a 2D array that holds x-values for constant z1 and z2.
        yz is a 2D array that holds y-values for the same z1 and z2.
        xz and yz have to have the same size.
        The first dimension of xz should be greater than
        or equal to the lines input.
        """
        if xz.ndim!=2:
            raise ValueError("xz has to be a 2D array.")
        if yz.ndim!=2:
            raise ValueError("yz has to be a 2D array.")
        if xz.shape!=yz.shape:
            raise ValueError("xz and yz have to have the same shape: {0} != {1}".format(xz.shape,yz.shape))
        if linesY is None and not linesX is None:
            linesY = linesX
        if linesX is None and not linesY is None:
            linesX = linesY
        if linesY is None and linesX is None:
            raise ValueError("You have to provide linesX or linesY")

        xl,yl = xz.shape
        x_index = np.round(np.linspace(0, xl-1, linesX))
        y_index = np.round(np.linspace(0, yl-1, linesY))
        x_toPlot = []
        y_toPlot = []
        for i in x_index:
            x_toPlot += [xz.T[i]]
            y_toPlot += [yz.T[i]]
        for i in y_index:
            x_toPlot += [xz[i]]
            y_toPlot += [yz[i]]
        if plot==False:
            return x_toPlot,y_toPlot
        if ax is None:
            raise ValueError("You have to give an axis to plot.")
        for i in range(len(x_toPlot)):
            ax.plot(x_toPlot[i],y_toPlot[i],color=color)
        return x_toPlot,y_toPlot


    def plotValues(self,axVal,axErr,solObj=SolutionData(),dataObj=IncompressibleData(),func=None,old=None):
        """
        Plots two data series using the same axis. You can
        choose if you prefer points or a line for the reference
        data. This can be used to show that we have experimental
        data or a reference equation.
        You can use the old input to call CoolProp with this ID.
        You can use this feature to visualise changes for new
        data fits. Primarily intended to highlight changes from
        v4 to v5 of CoolProp.
        """

        if dataObj.type==dataObj.INCOMPRESSIBLE_NOT_SET \
          or dataObj.source==dataObj.SOURCE_NOT_SET:
            return

        # TODO: Improve this work-around
        xFunction = False
        try:
            if solObj.T_freeze.coeffs.shape==dataObj.coeffs.shape:
                if np.all(solObj.T_freeze.coeffs==dataObj.coeffs):
                    xFunction = True
        except AttributeError as ae:
            if False: print(ae)
            pass

        points = 30

        dataFormatter = {}
        tData = None
        xData = None
        pData = None
        zData = None
        zError= None

        if dataObj.source==dataObj.SOURCE_DATA or dataObj.source==dataObj.SOURCE_EQUATION:
            dataFormatter['color']  = 'blue'
            dataFormatter['marker'] = 'o'
            dataFormatter['ls']     = 'none'


            dataObj.setxyData(solObj.temperature.data,solObj.concentration.data)

            tData = dataObj.xData
            xData = dataObj.yData
            pData = 1e7 # 100 bar
            zData = dataObj.data

            if not func is None and not zData is None:
                r,c = zData.shape
                zError= np.zeros((r,c))
                for i in range(r):
                    for j in range(c):
                        zError[i,j]= func(tData[i],pData,xData[j])

                zError = self.relError(zData, zError) * 1e2


                ## Find the column with the largest single error
                #maxVal = np.amax(zError, axis=0) # largest error per column
                #col2plot = np.argmax(maxVal) # largest error in row
                ## Find the column with the largest total error
                #totVal = np.sum(zError, axis=0) # summed error per column
                #col2plot = np.argmax(totVal) # largest error in row
                # Find the column with the largest average error
                if xFunction:
                    #avgVal = np.average(zError, axis=1) # summed error per column
                    #set2plot = np.argmax(avgVal) # largest error in row
                    set2plot = int(np.round(r/2.0))
                    tData = np.array([tData[set2plot]])
                    zData = zData[set2plot]
                    zError= zError[set2plot]
                else:
                    #avgVal = np.average(zError, axis=0) # summed error per column
                    #set2plot = np.argmax(avgVal) # largest error in row
                    set2plot = int(np.round(c/2.0))
                    xData = np.array([xData[set2plot]])
                    zData = zData.T[set2plot]
                    zError= zError.T[set2plot]
            else:
                raise ValueError("You have to provide data and a fitted function.")

        elif dataObj.source==dataObj.SOURCE_COEFFS:
            dataFormatter['color']  = 'blue'
            #dataFormatter['marker'] = 'o'
            dataFormatter['ls']     = 'solid'

            if xFunction:
                xData = np.linspace(solObj.xmin, solObj.xmax, num=points)
                if solObj.xid==solObj.ifrac_pure: tData = np.array([0.0])
                else: tData = np.array([solObj.Tmin+solObj.Tmax])/2.0
            else:
                tData = np.linspace(solObj.Tmin, solObj.Tmax, num=points)
                if solObj.xid==solObj.ifrac_pure: xData = np.array([0.0])
                else: xData = np.array([solObj.xmin+solObj.xmax])/2.0

            pData = 1e7 # 100 bar

            #zData= np.zeros((len(tData),len(xData)))
            #for i in range(len(tData)):
            #    for j in range(len(xData)):
            #        zData[i,j] = func(tData[i],pData,xData[j])
            #r,c = zData.shape
            #if r==1 or c==1:
            #    zData = np.array(zData.flat)
            #else:
            #    raise ValueError("Cannot plot non-flat arrays!")


        # Copy the arrays
        tFunc = tData
        xFunc = xData
        pFunc = pData
        zFunc = None
        zMiMa = None
        xFree = xData
        tFree = tData
        zFree = None

        if not func is None:
            if len(tFunc)<points and len(tFunc)>1:
                tFunc = np.linspace(solObj.Tmin, solObj.Tmax, num=points)
            if len(xFunc)<points and len(xFunc)>1:
                xFunc = np.linspace(solObj.xmin, solObj.xmax, num=points)

            zFunc = np.zeros((len(tFunc),len(xFunc)))
            for i in range(len(tFunc)):
                for j in range(len(xFunc)):
                    zFunc[i,j] = func(tFunc[i],pFunc,xFunc[j])
            r,c = zFunc.shape
            if r==1 or c==1:
                zFunc = np.array(zFunc.flat)
            else:
                raise ValueError("Cannot plot non-flat arrays!")

            if xFunction:
                tMiMa = np.array([solObj.Tmin, solObj.Tmax])
                xMiMa = xFunc
            else:
                tMiMa = tFunc
                xMiMa = np.array([solObj.xmin, solObj.xmax])

            zMiMa = np.zeros((len(tMiMa),len(xMiMa)))
            for i in range(len(tMiMa)):
                for j in range(len(xMiMa)):
                    zMiMa[i,j] = func(tMiMa[i],pFunc,xMiMa[j])

            if not xFunction: # add the freezing front
                if solObj.T_freeze.type!=IncompressibleData.INCOMPRESSIBLE_NOT_SET:
                    cols = len(tMiMa)
                    conc = np.linspace(solObj.xmin, solObj.xmax, num=cols)
                    tFree = np.zeros_like(conc)
                    zFree = np.zeros_like(conc)
                    for i in range(cols):
                        tFree[i] = solObj.Tfreeze(10.0, p=pFunc, x=conc[i])
                        zFree[i] = func(tFree[i],pFunc,conc[i])
                    #zMiMa = np.hstack((zMiMa,temp.reshape((len(conc),1))))

        fitFormatter = {}
        fitFormatter['color'] = 'red'
        fitFormatter['ls'] = 'solid'

        errorFormatter = {}
        errorFormatter['color'] = dataFormatter['color']
        errorFormatter['marker'] = 'o'
        errorFormatter['ls'] = 'none'
        errorFormatter['alpha'] = 0.25

        pData = None
        pFree = None
        if xFunction:
            pData = xData
            pFunc = xFunc
            pMiMa = xMiMa
            zMiMa = zMiMa.T
            #pFree = xFree
            #zFree = zFree.T
        else:
            pData = tData - 273.15
            pFunc = tFunc - 273.15
            pMiMa = tMiMa - 273.15
            #zMiMa = zMiMa
            pFree = tFree - 273.15

        if not zData is None and not axVal is None:
            axVal.plot(pData, zData, label='data', **dataFormatter)

        if not zFunc is None and not axVal is None:
            axVal.plot(pFunc, zFunc, label='function' , **fitFormatter)
            if solObj.xid!=solObj.ifrac_pure and not xFunction:
                axVal.set_title("showing x={0:3.2f}".format(xFunc[0]))
            else:
                axVal.set_title(" ")

        if not zMiMa is None and not axVal is None:
            axVal.plot(pMiMa, zMiMa, alpha=0.25, ls=':', color=fitFormatter["color"])

        if not zFree is None and not axVal is None:
            axVal.plot(pFree, zFree, alpha=0.25, ls=':', color=fitFormatter["color"])

        if not zError is None and not axErr is None:
            #formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
            #axErr.yaxis.set_major_formatter(formatter)
            #axErr.yaxis.get_major_formatter().useOffset=False
            axErr.plot(pData, zError, label='error' , **errorFormatter)

        elif not axErr is None:
            errorFormatter['alpha'] = 0.00
            axErr.plot([pData[0],pData[-1]], [0,0], **errorFormatter)
            #axErr.xaxis.set_visible(False)
            #axErr.yaxis.set_visible(False)
            #axErr.plot(pData, zFunc, label='function' , **fitFormatter)





        #else:
        #    plt.setp(axErr.get_yticklabels(), visible=False)
        #    plt.setp(axErr.yaxis.get_label(), visible=False)


    def printFluidInfo(self,ax,solObj=SolutionData()):
        """
        Prints some fluid information on top of the fitting report.
        """
        #ax = subplot(111, frame_on=False)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)

        #cellText  = [["1","2"],["3","4"]]
        #rowLabels = ["Name","Source"]
        ## Add a table at the bottom of the axes
        #the_table = ax.table(cellText=cellText)

        annotateSettingsTitle = {}
        #annotateSettingsLabel['xycoords']=('figure fraction', 'figure fraction')
        annotateSettingsTitle['ha']         = 'center'
        annotateSettingsTitle['va']         = 'baseline'
        annotateSettingsTitle['fontsize']   = 'xx-large'
        annotateSettingsTitle['fontweight'] = 'bold'

        annotateSettingsLabel = {}
        #annotateSettingsLabel['xycoords']=('figure fraction', 'figure fraction')
        annotateSettingsLabel['ha']         = 'right'
        annotateSettingsLabel['va']         = 'baseline'
        #annotateSettingsLabel['fontsize']   = 'large'
        annotateSettingsLabel['fontweight'] = 'semibold'

        annotateSettingsText = {}
        #annotateSettingsText['xycoords']=('figure fraction', 'figure fraction')
        annotateSettingsText['ha']         = 'left'
        annotateSettingsText['va']         = 'baseline'
        #annotateSettingsText['fontsize']   = 'large'
        annotateSettingsText['fontweight'] = 'medium'


        #ax.set_title('Fitting Report for {0}'.format(solObj.name))
        ax.text(0.5, 0.8, 'Fitting Report for {0}'.format(solObj.name), **annotateSettingsTitle)

        def myAnnotate(label,text,x=0.0,y=0.0):
            dx = 0.005
            ax.text(x-dx,y, label, **annotateSettingsLabel)
            ax.text(x+dx,y, text,  **annotateSettingsText)

        #ax.annotate(r'Enthalpy [$\mathdefault{10^5\!J/kg}$]', xy=(0.25, 0.03), **annotateSettings)
        #ax.annotate(r'Enthalpy [$\mathdefault{10^5\!J/kg}$]', xy=(0.25, 0.03), **annotateSettings)

        dx = 0.50
        dy = 0.10

        xStart = 0.175 ; x = xStart
        yStart = 0.575; y = yStart

        #myAnnotate('Name: ',solObj.name,x=x,y=y); x += .0; y -= dy
        myAnnotate('Description: ',solObj.description,x=x,y=y); x += .0; y -= dy

        # TODO: Debug bibtexer
        refs = solObj.reference.split(",")
        maxLength = 75
        for i in range(len(refs)):
            refs[i] = refs[i].strip()
            try:
                refs[i] = self.bibtexer.getEntry(key=refs[i], fmt='plaintext').strip()
            except Exception as e:
                warn("Your string \"{0}\"was not a valid Bibtex key, I will use it directly: {1}".format(refs[i],e))
                pass
            if len(refs[i])>maxLength:
                refs[i] = refs[i][0:maxLength-3]+u'...'

            if i==0:
                myAnnotate('Source: ',refs[i],x=x,y=y); x += .0 #; y -= 2*dy
            elif i==1:
                myAnnotate(   '     ',refs[i],x=x,y=y-dy); x += .0 #; y -= 2*dy

        y -= 2*dy
        yRestart = y
        myAnnotate('Temperature: ',u'{0} \u00B0C to {1} \u00B0C'.format(solObj.Tmin-273.15, solObj.Tmax-273.15),x=x,y=y); x += .0; y -= dy
        conc = False
        if solObj.xid==solObj.ifrac_mass: conc=True
        if solObj.xid==solObj.ifrac_volume: conc=True
        if solObj.xid==solObj.ifrac_mole: conc=True
        if conc==True:
            myAnnotate('Composition: ',u'{0} % to {1} %, {2}'.format(solObj.xmin*100., solObj.xmax*100., solObj.xid),x=x,y=y)
        else:
            myAnnotate('Composition: ','pure fluid',x=x,y=y)
        x += .0; y -= dy

        if solObj.density.source!=solObj.density.SOURCE_NOT_SET:
            myAnnotate('Density: ',u'{0} to {1} {2}'.format(solObj.density.source, solObj.density.type, solObj.density.coeffs.shape),x=x,y=y)
        else:
            myAnnotate('Density: ','no information',x=x,y=y)
        x += .0; y -= dy
        if solObj.specific_heat.source!=solObj.specific_heat.SOURCE_NOT_SET:
            myAnnotate('Spec. Heat: ',u'{0} to {1} {2}'.format(solObj.specific_heat.source, solObj.specific_heat.type, solObj.specific_heat.coeffs.shape),x=x,y=y)
        else:
            myAnnotate('Spec. Heat: ','no information',x=x,y=y)
        x += .0; y -= dy

        x = xStart + dx; y = yRestart
        if solObj.conductivity.source!=solObj.conductivity.SOURCE_NOT_SET:
            myAnnotate('Th. Cond.: ',u'{0} to {1} {2}'.format(solObj.conductivity.source, solObj.conductivity.type, solObj.conductivity.coeffs.shape),x=x,y=y)
        else:
            myAnnotate('Th. Cond.: ','no information',x=x,y=y)
        x += .0; y -= dy
        if solObj.viscosity.source!=solObj.viscosity.SOURCE_NOT_SET:
            myAnnotate('Viscosity: ',u'{0} to {1} {2}'.format(solObj.viscosity.source, solObj.viscosity.type, solObj.viscosity.coeffs.shape),x=x,y=y)
        else:
            myAnnotate('Viscosity: ','no information',x=x,y=y)
        x += .0; y -= dy
        if solObj.saturation_pressure.source!=solObj.saturation_pressure.SOURCE_NOT_SET:
            myAnnotate('Psat: ',u'{0} to {1} {2}'.format(solObj.saturation_pressure.source, solObj.saturation_pressure.type, solObj.saturation_pressure.coeffs.shape),x=x,y=y)
        else:
            myAnnotate('Psat: ','no information',x=x,y=y)
        x += .0; y -= dy
        if solObj.T_freeze.source!=solObj.T_freeze.SOURCE_NOT_SET:
            myAnnotate('Tfreeze: ',u'{0} to {1} {2}'.format(solObj.T_freeze.source, solObj.T_freeze.type, solObj.T_freeze.coeffs.shape),x=x,y=y)
        else:
            myAnnotate('Tfreeze: ','no information',x=x,y=y)
        x += .0; y -= dy


        #ax5.set_xlabel(ur'$\mathregular{Temperature\/(\u00B0C)}$')

        #x += dx; y = yStart
        #myAnnotate('Name: ',solObj.name,x=x,y=y); x += .0; y -= dy

        ax.set_xlim((0,1))
        ax.set_ylim((0,1))


    def printFitDetails(self):
        pass



    def makeFitReportPage(self, solObj=SolutionData(), pdfObj=None, quiet=False):
        """
        Creates a whole page with some plots and basic information
        for both fit quality, reference data, data sources and
        more.
        """

        # First we determine some basic settings and check the JSON file

        json_path = self.get_json_file(solObj.name)
        report_path = self.get_report_file(solObj.name)

        if os.path.isfile(json_path):
            json_time = os.path.getctime(json_path)
        else:
            json_time = 0

        if os.path.isfile(report_path):
            report_time = os.path.getctime(report_path)
        else:
            report_time = 0

        if json_time<report_time and pdfObj is None:
            #print("The JSON file {0} has not been updated, skipping.".format(self.get_json_file(solObj.name)))
            if not quiet: print(" ({0})".format("i"), end="")
            return 1

#         from os import path
# from time import ctime
# from datetime import datetime, timedelta
#
# two_days_ago = datetime.now() - timedelta(days=2)
# filetime = datetime.fromtimestamp(path.getctime(file))
#
# if filetime < two_days_ago:
#   print "File is more than two days old"




        gs = gridspec.GridSpec(4, 2, wspace=None, hspace=None, height_ratios=[2.5,3,3,3])
        #gs.update(top=0.75, hspace=0.05)
        div = 22
        fig = plt.figure(figsize=(210/div,297/div))
        table_axis         = plt.subplot(gs[0,:], frame_on=False)




        # Info text settings
        infoText = {}
        infoText['ha']         = 'center'
        infoText['va']         = 'baseline'
        infoText['fontsize']   = 'smaller'
        infoText['xycoords']   =('axes fraction', 'axes fraction')

        # Setting the labels
        #errLabel  = ur'$\mathdefault{rel.\/Error\/[\u2030]}$'
        errLabel  = r'$\mathdefault{rel.\/Error\/[\%]}$'
        tempLabel = ur'$\mathdefault{Temperature\/(\u00B0C)}$'

        density_axis       = plt.subplot(gs[1,0])
        density_error      = density_axis.twinx()
        density_axis.set_ylabel(r'Density [$\mathdefault{kg/m^3\!}$]')
        density_axis.set_xlabel(tempLabel)
        density_error.set_ylabel(errLabel)
        if solObj.density.source!=solObj.density.SOURCE_NOT_SET:
            self.plotValues(density_axis,density_error,solObj=solObj,dataObj=solObj.density,func=solObj.rho)
        else:
            raise ValueError("Density data has to be provided!")

        capacity_axis      = plt.subplot(gs[1,1])
        capacity_error     = capacity_axis.twinx()
        capacity_axis.set_ylabel(r'Heat Capacity [$\mathdefault{J/kg/K}$]')
        capacity_axis.set_xlabel(tempLabel)
        capacity_error.set_ylabel(errLabel)
        if solObj.specific_heat.source!=solObj.specific_heat.SOURCE_NOT_SET:
            self.plotValues(capacity_axis,capacity_error,solObj=solObj,dataObj=solObj.specific_heat,func=solObj.c)
        else:
            raise ValueError("Specific heat data has to be provided!")

        # Optional plots, might not all be shown
        conductivity_axis  = plt.subplot(gs[2,0])#, sharex=density_axis)
        conductivity_error = conductivity_axis.twinx()
        conductivity_axis.set_ylabel(r'Thermal Conductivity [$\mathdefault{W/m/K}$]')
        conductivity_axis.set_xlabel(tempLabel)
        conductivity_error.set_ylabel(errLabel)
        if solObj.conductivity.source!=solObj.conductivity.SOURCE_NOT_SET:
            self.plotValues(conductivity_axis,conductivity_error,solObj=solObj,dataObj=solObj.conductivity,func=solObj.cond)
        else:
            #conductivity_axis.xaxis.set_visible(False)
            #conductivity_axis.yaxis.set_visible(False)
            #conductivity_error.xaxis.set_visible(False)
            #conductivity_error.yaxis.set_visible(False)
            conductivity_axis.annotate("No conductivity information",xy=(0.5,0.5),**infoText)

        viscosity_axis     = plt.subplot(gs[2,1])#, sharex=capacity_axis)
        viscosity_error    = viscosity_axis.twinx()
        viscosity_axis.set_yscale('log')
        viscosity_axis.set_ylabel(r'Dynamic Viscosity [$\mathdefault{Pa\/s}$]')
        viscosity_axis.set_xlabel(tempLabel)
        viscosity_error.set_ylabel(errLabel)
        if solObj.viscosity.source!=solObj.viscosity.SOURCE_NOT_SET:
            self.plotValues(viscosity_axis,viscosity_error,solObj=solObj,dataObj=solObj.viscosity,func=solObj.visc)
        else:
            #viscosity_axis.xaxis.set_visible(False)
            #viscosity_axis.yaxis.set_visible(False)
            #viscosity_error.xaxis.set_visible(False)
            #viscosity_error.yaxis.set_visible(False)
            viscosity_axis.annotate("No viscosity information",xy=(0.5,0.5),**infoText)

        saturation_axis  = plt.subplot(gs[3,0])#, sharex=density_axis)
        saturation_error = saturation_axis.twinx()
        saturation_axis.set_yscale('log')
        saturation_axis.set_ylabel(r'Saturation Pressure [$\mathdefault{Pa}$]')
        saturation_axis.set_xlabel(tempLabel)
        saturation_error.set_ylabel(errLabel)
        if solObj.saturation_pressure.source != solObj.saturation_pressure.SOURCE_NOT_SET: # exists
            self.plotValues(saturation_axis,saturation_error,solObj=solObj,dataObj=solObj.saturation_pressure,func=solObj.psat)
        else:
            #saturation_axis.xaxis.set_visible(False)
            #saturation_axis.yaxis.set_visible(False)
            #saturation_error.xaxis.set_visible(False)
            #saturation_error.yaxis.set_visible(False)
            saturation_axis.annotate("No saturation state information",xy=(0.5,0.5),**infoText)

        Tfreeze_axis     = plt.subplot(gs[3,1])#, sharex=capacity_axis)
        Tfreeze_error    = Tfreeze_axis.twinx()
        Tfreeze_axis.set_ylabel(r'Freezing Temperature [$\mathdefault{K}$]')
        Tfreeze_axis.set_xlabel("{0} fraction".format(solObj.xid.title()))
        Tfreeze_error.set_ylabel(errLabel)
        if solObj.T_freeze.source != solObj.T_freeze.SOURCE_NOT_SET: # exists
            self.plotValues(Tfreeze_axis,Tfreeze_error,solObj=solObj,dataObj=solObj.T_freeze,func=solObj.Tfreeze)
        else:
            #Tfreeze_axis.xaxis.set_visible(False)
            #Tfreeze_axis.yaxis.set_visible(False)
            #Tfreeze_error.xaxis.set_visible(False)
            #Tfreeze_error.yaxis.set_visible(False)
            Tfreeze_axis.annotate("No freezing point information",xy=(0.5,0.5),**infoText)
            Tfreeze_axis.set_xlabel("Fraction")

        #saturation_axis   = plt.subplot2grid((3,2), (2,0))
        #Tfreeze_axis       = plt.subplot2grid((3,2), (2,0))
        #mass2input_axis   = plt.subplot2grid((3,2), (2,0))
        #volume2input_axis = plt.subplot2grid((3,2), (2,0))

        # Set a minimum error level and do some more formatting
        minAbsErrorScale = 0.05 # in per cent
        for a in fig.axes:
            if a.get_ylabel()==errLabel:
                mi,ma = a.get_ylim()
                if mi>-minAbsErrorScale: a.set_ylim(bottom=-minAbsErrorScale)
                if ma< minAbsErrorScale: a.set_ylim(   top= minAbsErrorScale)
            a.xaxis.set_major_locator(MaxNLocator(5))
            #a.yaxis.set_major_locator(MaxNLocator(7))



        # print headlines etc.
        self.printFluidInfo(table_axis, solObj)
        # Prepare the legend
        legenddict = {}
        for a in fig.axes:
            handles, labels = a.get_legend_handles_labels()
            for i in range(len(labels)):
                legenddict[labels[i]] = handles[i]

        legKey  = ["Legend: "]
        legVal  = [Rectangle((0, 0), 1, 1, alpha=0.0)]
        legKey += legenddict.keys()
        legVal += legenddict.values()
        legKey += [" "]
        legVal += [Rectangle((0, 0), 1, 1, alpha=0.0)]

        table_axis.legend(
          legVal, legKey,
          bbox_to_anchor=(0.0, -0.03, 1., -0.03),
          ncol=len(legKey), mode="expand", borderaxespad=0.,
          numpoints=1)
        #table_axis.legend(handles, labels, bbox_to_anchor=(0.0, -0.1), loc=2, ncol=3)

        gs.tight_layout(fig)#, rect=[0, 0, 1, 0.75])
        # Fine-tune figure; make subplots close to each other
        # and hide x ticks for all but bottom plot.
        #fig.subplots_adjust(wspace=0)
        #plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)

        if json_time<report_time:
            if not quiet: print(" ({0})".format("i"), end="")
        else:
            if not os.path.exists(os.path.dirname(report_path)):
                os.makedirs(os.path.dirname(report_path))
            plt.savefig(report_path)
            if not quiet: print(" ({0})".format("w"), end="")

        if not pdfObj is None: pdfObj.savefig(fig)

        plt.close()
        pass


    def makeSolutionPlots(self, solObjs=[SolutionData()], pdfObj=None):
        """
        Creates a whole page with some plots and basic information
        for both fit quality, reference data, data sources and
        more.
        """
        # First we determine some basic settings
        water=None
        solutions=[]

        for i in range(len(solObjs)-1):
            if solObjs[i].xid==SolutionData.ifrac_mass or \
               solObjs[i].xid==SolutionData.ifrac_mole or \
               solObjs[i].xid==SolutionData.ifrac_volume:
                solutions += [solObjs[i]]
            #elif solObjs[i].xid==SolutionData.ifrac_pure:
            #    purefluids += [doneObjs[i]]
            elif solObjs[i].name=="NBS":
                water = solObjs[i]
                solutions += [solObjs[i]]

        if water is None: raise ValueError("No water found, reference values missing.")

        # Set temperature data for all fluids
        dataDict = {}
        dataList = []
        obj = SolutionData()
        for i in range(len(solutions)-1):
            obj = solutions[i]
            T = np.arange(np.max([275,np.round(obj.Tmin)]), np.min([300,np.round(obj.Tmax)]), 1)
            P = 100e5
            x = obj.xmin
            dataDict["name"] = obj.name
            dataDict["desc"] = obj.description
            dataDict["T"] = T
            dataDict["P"] = P
            dataDict["x"] = x
            if obj.density.type!=IncompressibleData.INCOMPRESSIBLE_NOT_SET:
                dataDict["D"] = np.array([obj.rho(Ti, P, x) for Ti in T])
            else: dataDict["D"] = None
            if obj.specific_heat.type!=IncompressibleData.INCOMPRESSIBLE_NOT_SET:
                dataDict["C"] = np.array([obj.c(Ti, P, x) for Ti in T])
            else: dataDict["C"] = None
            if obj.conductivity.type!=IncompressibleData.INCOMPRESSIBLE_NOT_SET:
                dataDict["L"] = np.array([obj.cond(Ti, P, x) for Ti in T])
            else: dataDict["L"] = None
            if obj.viscosity.type!=IncompressibleData.INCOMPRESSIBLE_NOT_SET:
                dataDict["V"] = np.array([obj.visc(Ti, P, x) for Ti in T])
            else: dataDict["V"] = None
            dataList.append(dataDict.copy())

        rat = np.array([5.5,3.05])
        mul =  2.25
        #fig = plt.figure(figsize=(297/div,210/div))
        fig = plt.figure(figsize=rat*mul)

        ax = fig.add_subplot(121)
        ax.set_xlabel(r'Temperature [$\mathdefault{K}$]')
        ax.set_ylabel(r'Density [$\mathdefault{kg/m^3\!}$]')
        fig.suptitle(r'Aqueous solutions with a concentration of 0.0',fontsize='x-large',fontweight='bold')

        #obj = water
        #T = np.linspace(obj.Tmin, obj.Tmax, num=int(obj.Tmax-obj.Tmin))
        #D =

        #import matplotlib.pyplot as plt
        #from itertools import cycle
        lines   = ["-","--","-.", ":"]
        colours = ['r', 'g', 'b', 'c', 'm']

        linecycler   = itertools.cycle(lines)
        colourcycler = itertools.cycle(colours)

        for i in range(len(dataList)-1):
            obj = dataList[i]
            if not obj["T"] is None and not obj["D"] is None:
                if obj["x"]==0.0 and np.any(np.isfinite(obj["D"])) and obj["name"]!="NBS":
                    lc = next(colourcycler)
                    if lc==colours[0]: ls = next(linecycler)
                    ax.plot(obj["T"],obj["D"],label="{0}: {1}".format(obj["name"],obj["desc"]),ls=ls,color=lc)
                #if not np.any(np.isfinite(obj["D"])):
                #    print("Name: {0}, Dmin: {1}, Dmax: {1}".format(obj["name"],np.min(obj["D"]),np.max(obj["D"])))

        obj = (item for item in dataList if item["name"] == "NBS").next()
        ax.plot(obj["T"],obj["D"],label="{0}: {1}".format(obj["name"],obj["desc"]),ls='-',color='black')

        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=1)#, prop={'size':'smaller'})
        plt.tight_layout(rect=(0, 0, 1, 0.95))
        plt.savefig("all_solutions_00.pdf")
        return 0



    #####################################
    # Table generation routines
    #####################################
    # See http://stackoverflow.com/questions/11347505/what-are-some-approaches-to-outputting-a-python-data-structure-to-restructuredte
    def make_table(self, grid):
        max_cols = [max(out) for out in map(list, zip(*[[len(item) for item in row] for row in grid]))]
        rst = self.table_div(max_cols, 1)
        for i, row in enumerate(grid):
            header_flag = False
            if i == 0 or i == len(grid)-1: header_flag = True
            rst += self.normalize_row(row,max_cols)
            rst += self.table_div(max_cols, header_flag )
        return rst

    def table_div(self, max_cols, header_flag=1, indent=2):
        out = ""
        for i in range(indent):
            out += " "
        if header_flag == 1:
            style = "="
        else:
            style = "-"
        for max_col in max_cols:
            out += max_col * style + " "
        out += "\n"
        return out

    def normalize_row(self, row, max_cols, indent=2):
        r = ""
        for i in range(indent):
            r += " "

        for i, max_col in enumerate(max_cols):
            r += row[i] + (max_col  - len(row[i]) + 1) * " "
        return r + "\n"


    def writeTextToFile(self, path,text):
        #print("Writing to file: {0}".format(path))
        if not os.path.exists(os.path.dirname(path)):
            os.makedirs(os.path.dirname(path))
        with open(path, 'w') as f:
            f.write(text)
        return True

    def writeTxtTableToFile(self, path,table,head=""):
        if not head == "":
            return self.writeTextToFile(path+".txt", head+"\n\n"+self.make_table(table))
        return self.writeTextToFile(path+".txt", self.make_table(table))

    def writeCsvTableToFile(self, path,table):
        if not os.path.exists(os.path.dirname(path+".csv")):
            os.makedirs(os.path.dirname(path+".csv"))
        with open(path+".csv", 'wb') as f:
            writer = csv.writer(f)
            writer.writerows(table)
        return True

    # Interface
    def writeTableToFile(self, path,table):
        self.writeCsvTableToFile(path,table)
        self.writeTxtTableToFile(path,table)
        return True



    def getReportLink(self, name):
        reportFile = os.path.join("..","_static","fluid_properties","incompressible","report","{0}_fitreport.pdf".format(name))
        return self.d(name,reportFile)

    def checkForNumber(self, number):
        try:
            n = float(number)
        except:
            n = np.NAN
            pass
        return n

    def d(self, text,target):
#         try:
#             if os.path.isfile(target):
#                 link = ":download:`{0}<{1}>`".format(text,target)
#             else:
#                 link = "{0}".format(text)
#         except:
#             link = "{0}".format(text)
#             pass
        # TODO: Fix this!
        link = ":download:`{0}<{1}>`".format(text,target)
        return link

    def m(self, math):
        text = ":math:`{0}`".format(math)
        return text

    def c(self, number):
        #text = "{0:5.2f} |degC|".format(self.checkForNumber(number)-273.15)
        text = "{0:5.2f}".format(self.checkForNumber(number)-273.15)
        return text

    def x(self, number):
        text = "{0:3.2f}".format(self.checkForNumber(number))
        return text


    def generateRstTable(self, solObjs=[SolutionData()], path="table"):

        solObjs = sorted(solObjs, key=lambda x: x.name)
        xmin = np.array([x.xmin for x in solObjs])
        xmax = np.array([x.xmax for x in solObjs])

        if np.any(xmin>0.0) and np.any(xmax<1.0): use_x = True
        else: use_x = False

        header = ['Name', 'Description', 'Reference', \
          self.m('T_{min}')+" (|degC|)", self.m('T_{max}')+" (|degC|)"]
        if use_x: header.extend([self.m('x_{min}'), self.m('x_{max}')])

        testTable = []
        testTable.append(header) # Headline
        for fluid in solObjs:
            testTable.append([self.getReportLink(fluid.name), fluid.description, fluid.reference, self.c(fluid.Tmin), self.c(fluid.Tmax)])
            if use_x: testTable[-1].extend([self.x(fluid.xmin), self.x(fluid.xmax)])

        self.writeTableToFile(path, testTable)
        return True





