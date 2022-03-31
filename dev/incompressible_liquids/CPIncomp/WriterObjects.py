#!/usr/bin/env python
# -*- coding: utf8 -*-
from __future__ import division, print_function
import numpy as np

import matplotlib
import copy
matplotlib.use("agg")

import hashlib, os, json, sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from CPIncomp.DataObjects import SolutionData
from CPIncomp.BaseObjects import IncompressibleData, IncompressibleFitter
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator
from matplotlib.backends.backend_pdf import PdfPages
import itertools
from CoolProp.BibtexParser import BibTeXerClass
from warnings import warn

# See: https://docs.python.org/2/library/csv.html#csv-examples
import csv, codecs, io


class UTF8Recoder:
    """
    Iterator that reads an encoded stream and reencodes the input to UTF-8
    """

    def __init__(self, f, encoding):
        self.reader = codecs.getreader(encoding)(f)

    def __iter__(self):
        return self

    def next(self):
        return next(self.reader).encode("utf-8")


class UnicodeReader:
    """
    A CSV reader which will iterate over lines in the CSV file "f",
    which is encoded in the given encoding.
    """

    def __init__(self, f, dialect=csv.excel, encoding="utf-8", **kwds):
        f = UTF8Recoder(f, encoding)
        self.reader = csv.reader(f, dialect=dialect, **kwds)

    def next(self):
        row = next(self.reader)
        return [unicode(s, "utf-8") for s in row]

    def __iter__(self):
        return self


class UnicodeWriter:
    """
    A CSV writer which will write rows to CSV file "f",
    which is encoded in the given encoding.
    """

    def __init__(self, f, dialect=csv.excel, encoding="utf-8", **kwds):
        # Redirect output to a queue
        self.queue = io.StringIO()
        self.writer = csv.writer(self.queue, dialect=dialect, **kwds)
        self.stream = f
        self.encoder = codecs.getincrementalencoder(encoding)()

    def writerow(self, row):
        self.writer.writerow(row)
        # Fetch UTF-8 output from the queue ...
        data = self.queue.getvalue()
        try: data = data.decode("utf-8")
        except: pass
        try: data = str(data, encoding ="utf-8")
        except: pass
        # ... and re-encode it into the target encoding
        data = self.encoder.encode(data)
        # write to the target stream
        self.stream.write(data)
        # empty queue
        self.queue.truncate(0)

    def writerows(self, rows):
        for row in rows:
            self.writerow(row)


class SolutionDataWriter(object):
    """
    A base class that defines all the variables needed
    in order to make a proper fit. You can copy this code
    put in your data and add some documentation for where the
    information came from.
    """

    def __init__(self):
        bibFile = os.path.join(os.path.dirname(__file__), '../../../CoolPropBibTeXLibrary.bib')
        self.bibtexer = BibTeXerClass(bibFile)

        matplotlib.rcParams['axes.formatter.useoffset'] = False

        # For standard report generation
        self.ext = "pdf"
        self.usebp = False
        self.ispage = True  # Do you want a page or a figure?
        self.isportrait = True
        self.resolveRef = True  # Resolve references and print text

        # Latex document mode
        #self.ext        = "pgf"
        #self.usebp      = True
        # self.ispage     = False # Do you want a page or a figure?
        #self.isportrait = True
        # self.resolveRef = False # Resolve references and print text

        if self.ext == "pgf" or matplotlib.rcParams['text.usetex']:
            self.usetex = True
            matplotlib.rcParams['text.usetex'] = True
            # preamble = [r'\usepackage{hyperref}', r'\usepackage{siunitx}']#, r'\usepackage[style=alphabetic,natbib=true,backend=biber]{biblatex}']
            preamble = [r'\usepackage{hyperref}', r'\usepackage{siunitx}']
            #matplotlib.rcParams['text.latex.preamble'] = list(matplotlib.rcParams['text.latex.preamble']) + preamble
            #matplotlib.rcParams["pgf.preamble"] = list(matplotlib.rcParams['pgf.preamble']) + preamble
            matplotlib.rcParams['text.latex.preamble'] = preamble
            matplotlib.rcParams["pgf.preamble"] = preamble
            self.percent = r'\si{\percent}'
            self.celsius = r'\si{\celsius}'
            self.errLabel = r'rel. Error (' + self.percent + r')'
            self.tempLabel = r'Temperature (' + self.celsius + r')'
            self.densLabel = r'Density (\si{\kilo\gram\per\cubic\metre})'
            self.heatLabel = r'Heat Capacity (\si{\joule\per\kilo\gram\per\kelvin})'
            self.condLabel = r'Thermal Conductivity (\si{\watt\per\metre\per\kelvin})'
            self.viscLabel = r'Dynamic Viscosity (\si{\pascal\second})'
            self.satPLabel = r'Saturation Pressure (\si{\pascal})'
            self.TfreLabel = r'Freezing Temperature (\si{\kelvin})'
            f = 0.85
        else:
            self.usetex = False
            matplotlib.rcParams['text.usetex'] = False
            self.percent = r'%'
            self.celsius = u'\u00B0C'
            self.errLabel = r'rel. Error (' + self.percent + r')'
            self.tempLabel = u'Temperature (' + self.celsius + u')'
            self.densLabel = r'Density ($\mathdefault{kg/m^3\!}$)'
            self.heatLabel = r'Heat Capacity ($\mathdefault{J/kg/K}$)'
            self.condLabel = r'Thermal Conductivity ($\mathdefault{W/m/K}$)'
            self.viscLabel = r'Dynamic Viscosity ($\mathdefault{Pa\/s}$)'
            self.satPLabel = r'Saturation Pressure ($\mathdefault{Pa}$)'
            self.TfreLabel = r'Freezing Temperature ($\mathdefault{K}$)'
            f = 1.00

#         self.percent   = ur'pc'
#         self.celsius   = ur'degC'
#         self.errLabel  = ur'rel. Error ('+self.percent+ur')'
#         self.tempLabel = ur'Temperature ('+self.celsius+ur')'
#         self.densLabel = ur'Density (kg/m3)'
#         self.heatLabel = ur'Heat Capacity (J/kg/K)'
#         self.condLabel = ur'Thermal Conductivity (W/m/K)'
#         self.viscLabel = ur'Dynamic Viscosity (Pa s)'
#         self.satPLabel = ur'Saturation Pressure (Pa)'
#         self.TfreLabel = ur'Freezing Temperature (K)'

        if self.usebp:
            from jopy.dataPlotters import BasePlotter
            self.bp = BasePlotter()
            ccycle = self.bp.getColourCycle(length=2)
            # ccycle.next() # skip the first one
            # ccycle.next() # skip the first one
            self.secondaryColour = next(ccycle)
            self.primaryColour = next(ccycle)
        else:
            self.primaryColour = 'blue'
            self.secondaryColour = 'red'

        baseSize = 297.0 * f  # A4 in mm
        ratio = 210.0 / 297.0  # A4 in mm
        mm_to_inch = 3.93700787401575 / 100.0  # factor mm to inch

        self.deta = 0.75  # factor for table
        if self.ispage:
            longSide = baseSize  # Make A4
            shortSide = baseSize * ratio
        else:
            # lofa = 1.0 - self.deta * 2.5/11.5 # Half of the original table
            longSide = baseSize  # * lofa #* 0.85 # Make smaller than A4
            shortSide = baseSize * ratio  # TODO: connect to gridspec and ylim

        if self.isportrait: self.figsize = (shortSide * mm_to_inch, longSide * mm_to_inch)
        else: self.figsize = (longSide * mm_to_inch, shortSide * mm_to_inch)

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
        std_xorder = 3 + 1
        std_yorder = 5 + 1
        std_coeffs = np.zeros((std_xorder, std_yorder))

        errList = (ValueError, AttributeError, TypeError, RuntimeError)

        if fluidObject.density.coeffs is None:
            try:
                fluidObject.density.setxyData(tData, xData)
                fluidObject.density.coeffs = np.copy(std_coeffs)
                fluidObject.density.type = IncompressibleData.INCOMPRESSIBLE_POLYNOMIAL
                fluidObject.density.fitCoeffs(tBase, xBase)
            except errList as ve:
                if fluidObject.density.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(fluidObject.name, 'density', ve))
                pass

        if fluidObject.specific_heat.coeffs is None:
            try:
                fluidObject.specific_heat.setxyData(tData, xData)
                fluidObject.specific_heat.coeffs = np.copy(std_coeffs)
                fluidObject.specific_heat.type = IncompressibleData.INCOMPRESSIBLE_POLYNOMIAL
                fluidObject.specific_heat.fitCoeffs(tBase, xBase)
            except errList as ve:
                if fluidObject.specific_heat.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(fluidObject.name, 'specific heat', ve))
                pass

        if fluidObject.conductivity.coeffs is None:
            try:
                fluidObject.conductivity.setxyData(tData, xData)
                fluidObject.conductivity.coeffs = np.copy(std_coeffs)
                fluidObject.conductivity.type = IncompressibleData.INCOMPRESSIBLE_POLYNOMIAL
                fluidObject.conductivity.fitCoeffs(tBase, xBase)
            except errList as ve:
                if fluidObject.conductivity.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(fluidObject.name, 'conductivity', ve))
                pass

        if fluidObject.viscosity.coeffs is None:
            try:
                fluidObject.viscosity.setxyData(tData, xData)
                tried = False
                if len(fluidObject.viscosity.yData) == 1:  # and np.isfinite(fluidObject.viscosity.data).sum()<10:
                    fluidObject.viscosity.coeffs = np.array([+5e+2, -6e+1, +1e+1])
                    fluidObject.viscosity.type = IncompressibleData.INCOMPRESSIBLE_EXPONENTIAL
                    fluidObject.viscosity.fitCoeffs(tBase, xBase)
                    if fluidObject.viscosity.coeffs is None or IncompressibleFitter.allClose(fluidObject.viscosity.coeffs, np.array([+5e+2, -6e+1, +1e+1])):  # Fit failed
                        tried = True
                if len(fluidObject.viscosity.yData) > 1 or tried:
                    #fluidObject.viscosity.coeffs = np.zeros(np.round(np.array(std_coeffs.shape) * 1.5))
                    fluidObject.viscosity.coeffs = np.copy(std_coeffs)
                    fluidObject.viscosity.type = IncompressibleData.INCOMPRESSIBLE_EXPPOLYNOMIAL
                    fluidObject.viscosity.fitCoeffs(tBase, xBase)
            except errList as ve:
                if fluidObject.viscosity.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(fluidObject.name, 'viscosity', ve))
                pass

        if fluidObject.saturation_pressure.coeffs is None:
            try:
                fluidObject.saturation_pressure.setxyData(tData, xData)
                tried = False
                if len(fluidObject.saturation_pressure.yData) == 1:  # and np.isfinite(fluidObject.saturation_pressure.data).sum()<10:
                    fluidObject.saturation_pressure.coeffs = np.array([-5e+3, +6e+1, -1e+1])
                    fluidObject.saturation_pressure.type = IncompressibleData.INCOMPRESSIBLE_EXPONENTIAL
                    fluidObject.saturation_pressure.fitCoeffs(tBase, xBase)
                    if fluidObject.saturation_pressure.coeffs is None or IncompressibleFitter.allClose(fluidObject.saturation_pressure.coeffs, np.array([-5e+3, +6e+1, -1e+1])):  # Fit failed
                        tried = True
                if len(fluidObject.saturation_pressure.yData) > 1 or tried:
                    #fluidObject.saturation_pressure.coeffs = np.zeros(np.round(np.array(std_coeffs.shape) * 1.5))
                    fluidObject.saturation_pressure.coeffs = np.copy(std_coeffs)
                    fluidObject.saturation_pressure.type = IncompressibleData.INCOMPRESSIBLE_EXPPOLYNOMIAL
                    fluidObject.saturation_pressure.fitCoeffs(tBase, xBase)
            except errList as ve:
                if fluidObject.saturation_pressure.DEBUG: print("{0}: Could not fit polynomial {1} coefficients: {2}".format(fluidObject.name, 'saturation pressure', ve))
                pass

        # reset data for getArray and read special files
        if fluidObject.xid != fluidObject.ifrac_pure and fluidObject.xid != fluidObject.ifrac_undefined:
            if fluidObject.T_freeze.coeffs is None:
                fluidObject.T_freeze.setxyData([0.0], xData)
                try:
                    if len(fluidObject.T_freeze.xData) == 1:  # and np.isfinite(fluidObject.T_freeze.data).sum()<10:
                        fluidObject.T_freeze.coeffs = np.array([+7e+2, -6e+1, +1e+1])
                        fluidObject.T_freeze.type = IncompressibleData.INCOMPRESSIBLE_EXPONENTIAL
                    else:
                        fluidObject.T_freeze.coeffs = np.copy(std_coeffs)
                        fluidObject.T_freeze.type = IncompressibleData.INCOMPRESSIBLE_EXPPOLYNOMIAL
                    fluidObject.T_freeze.fitCoeffs(tBase, xBase)
                except errList as ve:
                    if fluidObject.T_freeze.DEBUG: print("{0}: Could not fit {1} coefficients: {2}".format(fluidObject.name, "T_freeze", ve))
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

    def get_hash(self, data):
        return hashlib.sha224(data.encode()).hexdigest()

    def get_hash_file(self):
        return os.path.join(os.path.dirname(__file__), 'data', "hashes.json")

    def load_hashes(self):
        hashes_fname = self.get_hash_file()
        if os.path.exists(hashes_fname):
            hashes = json.load(open(hashes_fname, 'r'))
        else:
            hashes = dict()
        return hashes

    def write_hashes(self, hashes):
        hashes_fname = self.get_hash_file()
        fp = open(hashes_fname, 'w')
        fp.write(json.dumps(hashes))
        fp.close()
        return True

    def get_json_file(self, name):
        return os.path.join("json", "{0}.json".format(name))

    def get_report_file(self, name):
        return os.path.join("report", "{0}_fitreport.{1}".format(name, self.ext))

    def toJSON(self, data, quiet=False):
        jobj = {}

        jobj['name'] = data.name  # Name of the current fluid
        jobj['description'] = data.description  # Description of the current fluid
        jobj['reference'] = data.reference  # Reference data for the current fluid

        jobj['Tmax'] = data.Tmax         # Maximum temperature in K
        jobj['Tmin'] = data.Tmin         # Minimum temperature in K
        jobj['xmax'] = data.xmax         # Maximum concentration
        jobj['xmin'] = data.xmin         # Minimum concentration
        jobj['xid'] = data.xid           # Concentration is mole, mass or volume-based
        jobj['TminPsat'] = data.TminPsat     # Minimum saturation temperature in K
        jobj['Tbase'] = data.Tbase       # Base value for temperature fits
        jobj['xbase'] = data.xbase       # Base value for concentration fits

        # data.temperature  # Temperature for data points in K
        # data.concentration  # Concentration data points in weight fraction
        jobj['density'] = data.density.toJSON()  # Density in kg/m3
        jobj['specific_heat'] = data.specific_heat.toJSON()  # Heat capacity in J/(kg.K)
        jobj['viscosity'] = data.viscosity.toJSON()  # Dynamic viscosity in Pa.s
        jobj['conductivity'] = data.conductivity.toJSON()  # Thermal conductivity in W/(m.K)
        jobj['saturation_pressure'] = data.saturation_pressure.toJSON()  # Saturation pressure in Pa
        jobj['T_freeze'] = data.T_freeze.toJSON()  # Freezing temperature in K
        jobj['mass2input'] = data.mass2input.toJSON()  # dd
        jobj['volume2input'] = data.volume2input.toJSON()  # dd
        jobj['mole2input'] = data.mole2input.toJSON()  # dd

        # Quickfix to allow building the docs with Python 3, see #1786
        try:
            original_float_repr = json.encoder.FLOAT_REPR
            # print json.dumps(1.0001)
            stdFmt = "1.{0}e".format(int(data.significantDigits - 1))
            #pr = np.finfo(float64).eps * 10.0
            # pr = np.finfo(float64).precision - 2 # stay away from numerical precision
            #json.encoder.FLOAT_REPR = lambda o: format(np.around(o,decimals=pr), stdFmt)
            json.encoder.FLOAT_REPR = lambda o: format(o, stdFmt)
            dump = json.dumps(jobj, indent=2, sort_keys=True)
            json.encoder.FLOAT_REPR = original_float_repr
        except:
            # Assume Python3
            stdFmt = "1.{0}e".format(int(data.significantDigits - 1))
            class RoundingEncoder(json.JSONEncoder):
                def default(self, obj):
                    if isinstance(obj, float):
                        return format(o, stdFmt)
                    # Let the base class default method raise the TypeError
                    return json.JSONEncoder.default(self, obj)
            dump = json.dumps(jobj, indent=2, sort_keys=True, cls=RoundingEncoder)


        # print dump
        hashes = self.load_hashes()
        hash = self.get_hash(dump)

        name = jobj['name']

        if name not in hashes or \
          hashes[name] != hash or \
          not os.path.isfile(self.get_json_file(name)):  # update hashes and write file

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

    def fromJSON(self, data=SolutionData()):

        path = os.path.join("json", data.name + '.json')

        if not os.path.isfile(path):
            return None

        with open(path) as json_file:
            jobj = json.load(json_file)

            data.name = jobj['name']  # Name of the current fluid
            data.description = jobj['description']  # Description of the current fluid
            data.reference = jobj['reference']  # Reference data for the current fluid

            data.Tmax = jobj['Tmax']          # Maximum temperature in K
            data.Tmin = jobj['Tmin']          # Minimum temperature in K
            data.xmax = jobj['xmax']          # Maximum concentration
            data.xmin = jobj['xmin']          # Minimum concentration
            data.xid = jobj['xid']            # Concentration is mole, mass or volume-based
            data.TminPsat = jobj['TminPsat']      # Minimum saturation temperature in K
            data.Tbase = jobj['Tbase']        # Base value for temperature fits
            data.xbase = jobj['xbase']        # Base value for concentration fits

            # data.temperature  # Temperature for data points in K
            # data.concentration  # Concentration data points in weight fraction
            data.density.fromJSON(jobj['density'])  # Density in kg/m3
            data.specific_heat.fromJSON(jobj['specific_heat'])  # Heat capacity in J/(kg.K)
            data.viscosity.fromJSON(jobj['viscosity'])  # Dynamic viscosity in Pa.s
            data.conductivity.fromJSON(jobj['conductivity'])  # Thermal conductivity in W/(m.K)
            data.saturation_pressure.fromJSON(jobj['saturation_pressure'])  # Saturation pressure in Pa
            data.T_freeze.fromJSON(jobj['T_freeze'])  # Freezing temperature in K
            data.mass2input.fromJSON(jobj['mass2input'])  # dd
            data.volume2input.fromJSON(jobj['volume2input'])  # dd
            data.mole2input.fromJSON(jobj['mole2input'])  # dd

            return data

    def printStatusID(self, fluidObjs, obj):
        #obj = fluidObjs[num]
        if obj == fluidObjs[0]:
            print(" {0}".format(obj.name), end="")
        elif obj == fluidObjs[-1]:
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
                print(str(e))
                pass
        print(" ... done")
        return

    def writeReportList(self, fluidObjs, pdfFile=None):
        print("Writing fitting reports:", end="")

        if self.usetex and pdfFile is not None:
            warn("Unsetting PDF object, LaTeX output requested.")
            pdfFile = None

        if not pdfFile is None:
            if not os.path.exists(os.path.dirname(pdfFile)):
                os.makedirs(os.path.dirname(pdfFile))
            with PdfPages(pdfFile) as pdfObj:
                for obj in fluidObjs:
                    matplotlib.pyplot.close("all")
                    self.printStatusID(fluidObjs, obj)
                    try:
                        self.makeFitReportPage(obj, pdfObj=pdfObj)
                    except (TypeError, ValueError) as e:
                        print("An error occurred for fluid: {0}".format(obj.name))
                        print(obj)
                        print(e)
                        pass
        else:
            for obj in fluidObjs:
                matplotlib.pyplot.close("all")
                self.printStatusID(fluidObjs, obj)
                self.makeFitReportPage(obj)
                # TODO: Fix Python errors
                # try:
                #    self.makeFitReportPage(obj)
                # except (TypeError, ValueError) as e:
                #    print("An error occurred for fluid: {0}".format(obj.name))
                #    print(obj)
                #    print(e)
                #    pass

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
        pos2 = (B_a > eps)
        result = np.ones_like(A_a) * np.NAN

        result[pos & pos2] = (A_a[pos & pos2] - B_a[pos & pos2]) / B_a[pos & pos2]

        if PCT:
            return result * 100.
        else:
            return result

    ############################################################
    # Define the general purpose routines for plotting
    ############################################################
    def wireFrame2D(self, xz, yz, linesX=5, linesY=None, color='black', ax=None, plot=False):
        """
        xz is a 2D array that holds x-values for constant z1 and z2.
        yz is a 2D array that holds y-values for the same z1 and z2.
        xz and yz have to have the same size.
        The first dimension of xz should be greater than
        or equal to the lines input.
        """
        if xz.ndim != 2:
            raise ValueError("xz has to be a 2D array.")
        if yz.ndim != 2:
            raise ValueError("yz has to be a 2D array.")
        if xz.shape != yz.shape:
            raise ValueError("xz and yz have to have the same shape: {0} != {1}".format(xz.shape, yz.shape))
        if linesY is None and not linesX is None:
            linesY = linesX
        if linesX is None and not linesY is None:
            linesX = linesY
        if linesY is None and linesX is None:
            raise ValueError("You have to provide linesX or linesY")

        xl, yl = xz.shape
        x_index = np.round(np.linspace(0, xl - 1, linesX))
        y_index = np.round(np.linspace(0, yl - 1, linesY))
        x_toPlot = []
        y_toPlot = []
        for i in x_index:
            x_toPlot += [xz.T[i]]
            y_toPlot += [yz.T[i]]
        for i in y_index:
            x_toPlot += [xz[i]]
            y_toPlot += [yz[i]]
        if plot == False:
            return x_toPlot, y_toPlot
        if ax is None:
            raise ValueError("You have to give an axis to plot.")
        for i in range(len(x_toPlot)):
            ax.plot(x_toPlot[i], y_toPlot[i], color=color)
        return x_toPlot, y_toPlot

    def plotValues(self, axVal, axErr, solObj=SolutionData(), dataObj=IncompressibleData(), func=None, old=None):
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

        if dataObj.type == dataObj.INCOMPRESSIBLE_NOT_SET \
          or dataObj.source == dataObj.SOURCE_NOT_SET:
            return

        # TODO: Improve this work-around
        xFunction = False
        try:
            if solObj.T_freeze.coeffs.shape == dataObj.coeffs.shape:
                if np.all(solObj.T_freeze.coeffs == dataObj.coeffs):
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
        zError = None

        if dataObj.source == dataObj.SOURCE_DATA or dataObj.source == dataObj.SOURCE_EQUATION:

            dataFormatter['color'] = self.primaryColour
            dataFormatter['marker'] = 'o'
            dataFormatter['ls'] = 'none'

            dataObj.setxyData(solObj.temperature.data, solObj.concentration.data)

            tData = dataObj.xData
            xData = dataObj.yData
            pData = 1e7  # 100 bar
            zData = dataObj.data

            if not func is None and not zData is None:
                r, c = zData.shape
                zError = np.zeros((r, c))
                for i in range(r):
                    for j in range(c):
                        zError[i, j] = func(tData[i], pData, xData[j])

                zError = self.relError(zData, zError) * 1e2

                # Find the column with the largest single error
                # maxVal = np.amax(zError, axis=0) # largest error per column
                # col2plot = np.argmax(maxVal) # largest error in row
                # Find the column with the largest total error
                # totVal = np.sum(zError, axis=0) # summed error per column
                # col2plot = np.argmax(totVal) # largest error in row
                # Find the column with the largest average error
                if xFunction:
                    # avgVal = np.average(zError, axis=1) # summed error per column
                    # set2plot = np.argmax(avgVal) # largest error in row
                    set2plot = int(np.round(r / 2.0))
                    tData = np.array([tData[set2plot]])
                    zData = zData[set2plot]
                    zError = zError[set2plot]
                else:
                    # avgVal = np.average(zError, axis=0) # summed error per column
                    # set2plot = np.argmax(avgVal) # largest error in row
                    set2plot = int(np.round(c / 2.0))
                    xData = np.array([xData[set2plot]])
                    zData = zData.T[set2plot]
                    zError = zError.T[set2plot]
            else:
                raise ValueError("You have to provide data and a fitted function.")

        elif dataObj.source == dataObj.SOURCE_COEFFS:
            dataFormatter['color'] = self.primaryColour
            #dataFormatter['marker'] = 'o'
            dataFormatter['ls'] = 'solid'

            if xFunction:
                xData = np.linspace(solObj.xmin, solObj.xmax, num=points)
                if solObj.xid == solObj.ifrac_pure: tData = np.array([0.0])
                else: tData = np.array([solObj.Tmin + solObj.Tmax]) / 2.0
            else:
                tData = np.linspace(solObj.Tmin, solObj.Tmax, num=points)
                if solObj.xid == solObj.ifrac_pure: xData = np.array([0.0])
                else: xData = np.array([solObj.xmin + solObj.xmax]) / 2.0

            pData = 1e7  # 100 bar

            #zData= np.zeros((len(tData),len(xData)))
            # for i in range(len(tData)):
            #    for j in range(len(xData)):
            #        zData[i,j] = func(tData[i],pData,xData[j])
            #r,c = zData.shape
            # if r==1 or c==1:
            #    zData = np.array(zData.flat)
            # else:
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
            if len(tFunc) < points and len(tFunc) > 1:
                tFunc = np.linspace(solObj.Tmin, solObj.Tmax, num=points)
            if len(xFunc) < points and len(xFunc) > 1:
                xFunc = np.linspace(solObj.xmin, solObj.xmax, num=points)

            zFunc = np.zeros((len(tFunc), len(xFunc)))
            for i in range(len(tFunc)):
                for j in range(len(xFunc)):
                    zFunc[i, j] = func(tFunc[i], pFunc, xFunc[j])
            r, c = zFunc.shape
            if r == 1 or c == 1:
                zFunc = np.array(zFunc.flat)
            else:
                raise ValueError("Cannot plot non-flat arrays!")

            if xFunction:
                tMiMa = np.array([solObj.Tmin, solObj.Tmax])
                xMiMa = xFunc
            else:
                tMiMa = tFunc
                xMiMa = np.array([solObj.xmin, solObj.xmax])

            zMiMa = np.zeros((len(tMiMa), len(xMiMa)))
            for i in range(len(tMiMa)):
                for j in range(len(xMiMa)):
                    zMiMa[i, j] = func(tMiMa[i], pFunc, xMiMa[j])

            if not xFunction:  # add the freezing front
                if solObj.T_freeze.type != IncompressibleData.INCOMPRESSIBLE_NOT_SET:
                    cols = len(tMiMa)
                    conc = np.linspace(solObj.xmin, solObj.xmax, num=cols)
                    tFree = np.zeros_like(conc)
                    zFree = np.zeros_like(conc)
                    for i in range(cols):
                        tFree[i] = solObj.Tfreeze(10.0, p=pFunc, x=conc[i])
                        zFree[i] = func(tFree[i], pFunc, conc[i])
                    #zMiMa = np.hstack((zMiMa,temp.reshape((len(conc),1))))

        fitFormatter = {}
        fitFormatter['color'] = self.secondaryColour
        fitFormatter['ls'] = 'solid'

        errorFormatter = {}
        errorFormatter['color'] = dataFormatter['color']
        errorFormatter['marker'] = 'o'
        errorFormatter['ls'] = 'none'
        errorFormatter['alpha'] = 0.5

        boundsFormatter = fitFormatter.copy()
        boundsFormatter['ls'] = ':'
        boundsFormatter['alpha'] = 0.75

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
            axVal.plot(pFunc, zFunc, label='function', **fitFormatter)
            if solObj.xid != solObj.ifrac_pure and not xFunction:
                axVal.set_title("showing x={0:3.2f}".format(xFunc[0]), fontsize='small')
            else:
                axVal.set_title(" ", fontsize='small')

        if not zMiMa is None and not axVal is None:
            axVal.plot(pMiMa, zMiMa, label='bounds', **boundsFormatter)

        if not zFree is None and not axVal is None:
            axVal.plot(pFree, zFree, label='bounds', **boundsFormatter)

        if not zError is None and not axErr is None:
            #formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
            # axErr.yaxis.set_major_formatter(formatter)
            # axErr.yaxis.get_major_formatter().useOffset=False
            axErr.plot(pData, zError, label='error', **errorFormatter)

        elif not axErr is None:
            errorFormatter['alpha'] = 0.00
            axErr.plot([pData[0], pData[-1]], [0, 0], **errorFormatter)
            # axErr.xaxis.set_visible(False)
            # axErr.yaxis.set_visible(False)
            #axErr.plot(pData, zFunc, label='function' , **fitFormatter)

        # else:
        #    plt.setp(axErr.get_yticklabels(), visible=False)
        #    plt.setp(axErr.yaxis.get_label(), visible=False)

    def printFluidInfo(self, ax, solObj=SolutionData()):
        """
        Prints some fluid information on top of the fitting report.
        """
        #ax = subplot(111, frame_on=False)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)

        #cellText  = [["1","2"],["3","4"]]
        #rowLabels = ["Name","Source"]
        # Add a table at the bottom of the axes
        #the_table = ax.table(cellText=cellText)

        annotateSettingsTitle = {}
        #annotateSettingsLabel['xycoords']=('figure fraction', 'figure fraction')
        annotateSettingsTitle['ha'] = 'center'
        annotateSettingsTitle['va'] = 'baseline'
        annotateSettingsTitle['fontsize'] = 'large'
        annotateSettingsTitle['fontweight'] = 'bold'

        annotateSettingsLabel = {}
        #annotateSettingsLabel['xycoords']=('figure fraction', 'figure fraction')
        annotateSettingsLabel['ha'] = 'right'
        annotateSettingsLabel['va'] = 'baseline'
        annotateSettingsLabel['fontsize'] = 'small'
        annotateSettingsLabel['fontweight'] = 'semibold'

        annotateSettingsText = {}
        #annotateSettingsText['xycoords']=('figure fraction', 'figure fraction')
        annotateSettingsText['ha'] = 'left'
        annotateSettingsText['va'] = 'baseline'
        annotateSettingsText['fontsize'] = 'small'
        annotateSettingsText['fontweight'] = 'medium'

        #ax.set_title('Fitting Report for {0}'.format(solObj.name))
        yHead = 0.0
        if self.ispage:
            yHead = 0.8
            ax.text(0.5, yHead, 'Fitting Report for {0}'.format(solObj.name), **annotateSettingsTitle)
            yHead += 0.2

        def myAnnotate(label, text, x=0.0, y=0.0):
            dx = 0.005
            ax.text(x - dx, y, label, **annotateSettingsLabel)
            ax.text(x + dx, y, text, **annotateSettingsText)

        #ax.annotate(r'Enthalpy [$\mathdefault{10^5\!J/kg}$]', xy=(0.25, 0.03), **annotateSettings)
        #ax.annotate(r'Enthalpy [$\mathdefault{10^5\!J/kg}$]', xy=(0.25, 0.03), **annotateSettings)

        dx = 0.50
        dy = 0.10

        xStart = 0.175; x = xStart
        yStart = 0.575; y = yStart

        #myAnnotate('Name: ',solObj.name,x=x,y=y); x += .0; y -= dy
        if self.usetex:
            myAnnotate('Description: ', solObj.description.encode("latex"), x=x, y=y); x += .0; y -= dy
        else:
            myAnnotate('Description: ', solObj.description, x=x, y=y); x += .0; y -= dy

        # TODO: Debug bibtexer
        refs = solObj.reference.split(",")
        maxLength = 75
        for i in range(len(refs)):
            refs[i] = refs[i].strip()

            if self.resolveRef:
                try:
                    if self.usetex:
                        refs[i] = self.bibtexer.getEntry(key=refs[i], fmt='latex').strip()
                    else:
                        refs[i] = self.bibtexer.getEntry(key=refs[i], fmt='plaintext').strip()
                except Exception as e:
                    warn("Your string \"{0}\" was not a valid Bibtex key, I will use it directly: {1}".format(refs[i], e))
                    pass
            else:
                if self.usetex:
                    refs[i] = r'\cite{' + str(refs[i]) + r'}'
                else:
                    refs[i] = str(refs[i])

            if len(refs[i]) > maxLength and not self.usetex:
                refs[i] = refs[i][0:maxLength - 3] + u'...'
            # elif self.usetex:
            #    refs[i] = r'\truncate{'+str(maxLength)+r'pt}{'+str(refs[i])+r'}'
            # elif i>0 and (len(refs[i])+len(refs[i-1]))<maxLength:
            #    refs[i-1] = refs[i-1]+r", "+refs[i]
            #    refs[i] = ""

        for i in range(len(refs)):
            if i == 0:
                myAnnotate('Source: ', refs[i], x=x, y=y); x += .0  # ; y -= 2*dy
            elif i == 1:
                myAnnotate('     ', refs[i], x=x, y=y - dy); x += .0  # ; y -= 2*dy
            else:
                warn("Discarding all reference after the second line")

        y -= 2 * dy
        yRestart = y
        myAnnotate('Temperature: ', u'{0} {2} to {1} {2}'.format(solObj.Tmin - 273.15, solObj.Tmax - 273.15, self.celsius), x=x, y=y); x += .0; y -= dy
        conc = False
        if solObj.xid == solObj.ifrac_mass: conc = True
        if solObj.xid == solObj.ifrac_volume: conc = True
        if solObj.xid == solObj.ifrac_mole: conc = True
        if conc == True:
            myAnnotate('Composition: ', u'{0} {3} to {1} {3}, {2}'.format(solObj.xmin * 100., solObj.xmax * 100., solObj.xid, self.percent), x=x, y=y)
        else:
            myAnnotate('Composition: ', 'pure fluid', x=x, y=y)
        x += .0; y -= dy

        if solObj.density.source != solObj.density.SOURCE_NOT_SET:
            myAnnotate('Density: ', u'{0} to {1} {2}'.format(solObj.density.source, solObj.density.type, solObj.density.coeffs.shape), x=x, y=y)
        else:
            myAnnotate('Density: ', 'no information', x=x, y=y)
        x += .0; y -= dy
        if solObj.specific_heat.source != solObj.specific_heat.SOURCE_NOT_SET:
            myAnnotate('Spec. Heat: ', u'{0} to {1} {2}'.format(solObj.specific_heat.source, solObj.specific_heat.type, solObj.specific_heat.coeffs.shape), x=x, y=y)
        else:
            myAnnotate('Spec. Heat: ', 'no information', x=x, y=y)
        x += .0; y -= dy

        x = xStart + dx; y = yRestart
        if solObj.conductivity.source != solObj.conductivity.SOURCE_NOT_SET:
            myAnnotate('Th. Cond.: ', u'{0} to {1} {2}'.format(solObj.conductivity.source, solObj.conductivity.type, solObj.conductivity.coeffs.shape), x=x, y=y)
        else:
            myAnnotate('Th. Cond.: ', 'no information', x=x, y=y)
        x += .0; y -= dy
        if solObj.viscosity.source != solObj.viscosity.SOURCE_NOT_SET:
            myAnnotate('Viscosity: ', u'{0} to {1} {2}'.format(solObj.viscosity.source, solObj.viscosity.type, solObj.viscosity.coeffs.shape), x=x, y=y)
        else:
            myAnnotate('Viscosity: ', 'no information', x=x, y=y)
        x += .0; y -= dy
        if solObj.saturation_pressure.source != solObj.saturation_pressure.SOURCE_NOT_SET:
            myAnnotate('Psat: ', u'{0} to {1} {2}'.format(solObj.saturation_pressure.source, solObj.saturation_pressure.type, solObj.saturation_pressure.coeffs.shape), x=x, y=y)
        else:
            myAnnotate('Psat: ', 'no information', x=x, y=y)
        x += .0; y -= dy
        if solObj.T_freeze.source != solObj.T_freeze.SOURCE_NOT_SET:
            myAnnotate('Tfreeze: ', u'{0} to {1} {2}'.format(solObj.T_freeze.source, solObj.T_freeze.type, solObj.T_freeze.coeffs.shape), x=x, y=y)
        else:
            myAnnotate('Tfreeze: ', 'no information', x=x, y=y)
        x += .0; y -= dy

        # ax5.set_xlabel(ur'$\mathregular{Temperature\/(\u00B0C)}$')

        #x += dx; y = yStart
        #myAnnotate('Name: ',solObj.name,x=x,y=y); x += .0; y -= dy

        xmin = 0
        xmax = 1

        ymin = y + 1 * dy
        ymax = max(yStart, yHead)

        ax.set_xlim((xmin, xmax))
        ax.set_ylim((ymin, ymax))

        #xpos,ypos = ax.transAxes.inverted().transform(ax.transData.transform((0,y)))
        return [xmin, y + 0.5 * dy, xmax, y + 0.5 * dy]

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

        if json_time < report_time and pdfObj is None:
            #print("The JSON file {0} has not been updated, skipping.".format(self.get_json_file(solObj.name)))
            if not quiet: print(" ({0})".format("i"), end="")
            return 1

        # from os import path
        # from time import ctime
        # from datetime import datetime, timedelta
        #
        # two_days_ago = datetime.now() - timedelta(days=2)
        # filetime = datetime.fromtimestamp(path.getctime(file))
        #
        # if filetime < two_days_ago:
        #   print "File is more than two days old"

        if self.ispage:
            gs = gridspec.GridSpec(4, 2, wspace=None, hspace=None, height_ratios=[2.50, 3, 3, 3])
        else:  # Reduce height of upper graph by two fifth, TODO: connect to ylim
            gs = gridspec.GridSpec(4, 2, wspace=None, hspace=None, height_ratios=[2.50 * self.deta, 3, 3, 3])

        #gs.update(top=0.75, hspace=0.05)
        fig = plt.figure(figsize=self.figsize)

        table_axis = plt.subplot(gs[0, :], frame_on=False)

        # Info text settings
        infoText = {}
        infoText['ha'] = 'center'
        infoText['va'] = 'baseline'
        infoText['fontsize'] = 'smaller'
        infoText['xycoords'] = ('axes fraction', 'axes fraction')

        density_axis = plt.subplot(gs[1, 0])
        density_error = density_axis.twinx()
        density_axis.set_ylabel(self.densLabel)
        density_axis.set_xlabel(self.tempLabel)
        density_error.set_ylabel(self.errLabel)
        if solObj.density.source != solObj.density.SOURCE_NOT_SET:
            self.plotValues(density_axis, density_error, solObj=solObj, dataObj=solObj.density, func=solObj.rho)
        else:
            raise ValueError("Density data has to be provided!")

        capacity_axis = plt.subplot(gs[1, 1])
        capacity_error = capacity_axis.twinx()
        capacity_axis.set_ylabel(self.heatLabel)
        capacity_axis.set_xlabel(self.tempLabel)
        capacity_error.set_ylabel(self.errLabel)
        if solObj.specific_heat.source != solObj.specific_heat.SOURCE_NOT_SET:
            self.plotValues(capacity_axis, capacity_error, solObj=solObj, dataObj=solObj.specific_heat, func=solObj.c)
        else:
            raise ValueError("Specific heat data has to be provided!")

        # Optional plots, might not all be shown
        conductivity_axis = plt.subplot(gs[2, 0])  # , sharex=density_axis)
        conductivity_error = conductivity_axis.twinx()
        conductivity_axis.set_ylabel(self.condLabel)
        conductivity_axis.set_xlabel(self.tempLabel)
        conductivity_error.set_ylabel(self.errLabel)
        if solObj.conductivity.source != solObj.conductivity.SOURCE_NOT_SET:
            self.plotValues(conductivity_axis, conductivity_error, solObj=solObj, dataObj=solObj.conductivity, func=solObj.cond)
        else:
            # conductivity_axis.xaxis.set_visible(False)
            # conductivity_axis.yaxis.set_visible(False)
            # conductivity_error.xaxis.set_visible(False)
            # conductivity_error.yaxis.set_visible(False)
            conductivity_axis.annotate("No conductivity information", xy=(0.5, 0.5), **infoText)

        viscosity_axis = plt.subplot(gs[2, 1])  # , sharex=capacity_axis)
        viscosity_error = viscosity_axis.twinx()
        viscosity_axis.set_yscale('log')
        viscosity_axis.set_ylabel(self.viscLabel)
        viscosity_axis.set_xlabel(self.tempLabel)
        viscosity_error.set_ylabel(self.errLabel)
        if solObj.viscosity.source != solObj.viscosity.SOURCE_NOT_SET:
            self.plotValues(viscosity_axis, viscosity_error, solObj=solObj, dataObj=solObj.viscosity, func=solObj.visc)
        else:
            # viscosity_axis.xaxis.set_visible(False)
            # viscosity_axis.yaxis.set_visible(False)
            # viscosity_error.xaxis.set_visible(False)
            # viscosity_error.yaxis.set_visible(False)
            viscosity_axis.annotate("No viscosity information", xy=(0.5, 0.5), **infoText)

        saturation_axis = plt.subplot(gs[3, 0])  # , sharex=density_axis)
        saturation_error = saturation_axis.twinx()
        saturation_axis.set_yscale('log')
        saturation_axis.set_ylabel(self.satPLabel)
        saturation_axis.set_xlabel(self.tempLabel)
        saturation_error.set_ylabel(self.errLabel)
        if solObj.saturation_pressure.source != solObj.saturation_pressure.SOURCE_NOT_SET:  # exists
            self.plotValues(saturation_axis, saturation_error, solObj=solObj, dataObj=solObj.saturation_pressure, func=solObj.psat)
        else:
            # saturation_axis.xaxis.set_visible(False)
            # saturation_axis.yaxis.set_visible(False)
            # saturation_error.xaxis.set_visible(False)
            # saturation_error.yaxis.set_visible(False)
            saturation_axis.annotate("No saturation state information", xy=(0.5, 0.5), **infoText)

        Tfreeze_axis = plt.subplot(gs[3, 1])  # , sharex=capacity_axis)
        Tfreeze_error = Tfreeze_axis.twinx()
        Tfreeze_axis.set_ylabel(self.TfreLabel)
        Tfreeze_axis.set_xlabel("{0} fraction".format(solObj.xid.title()))
        Tfreeze_error.set_ylabel(self.errLabel)
        if solObj.T_freeze.source != solObj.T_freeze.SOURCE_NOT_SET:  # exists
            self.plotValues(Tfreeze_axis, Tfreeze_error, solObj=solObj, dataObj=solObj.T_freeze, func=solObj.Tfreeze)
        else:
            # Tfreeze_axis.xaxis.set_visible(False)
            # Tfreeze_axis.yaxis.set_visible(False)
            # Tfreeze_error.xaxis.set_visible(False)
            # Tfreeze_error.yaxis.set_visible(False)
            Tfreeze_axis.annotate("No freezing point information", xy=(0.5, 0.5), **infoText)
            Tfreeze_axis.set_xlabel("Fraction")

        #saturation_axis   = plt.subplot2grid((3,2), (2,0))
        #Tfreeze_axis       = plt.subplot2grid((3,2), (2,0))
        #mass2input_axis   = plt.subplot2grid((3,2), (2,0))
        #volume2input_axis = plt.subplot2grid((3,2), (2,0))

        # Set a minimum error level and do some more formatting
        minAbsErrorScale = 0.05  # in per cent
        for a in fig.axes:
            if a.get_ylabel() == self.errLabel:
                mi, ma = a.get_ylim()
                if mi > -minAbsErrorScale: a.set_ylim(bottom=-minAbsErrorScale)
                if ma < minAbsErrorScale: a.set_ylim(top=minAbsErrorScale)
            a.xaxis.set_major_locator(MaxNLocator(5))
            # a.yaxis.set_major_locator(MaxNLocator(7))

        # print headlines etc.
        lims = self.printFluidInfo(table_axis, solObj)
        # Prepare the legend
        legenddict = {}
        for a in fig.axes:
            handles, labels = a.get_legend_handles_labels()
            for i in range(len(labels)):
                legenddict[labels[i]] = handles[i]

        legKey = ["Legend: "]
        legVal = [Rectangle((0, 0), 1, 1, alpha=0.0)]
        legKey += legenddict.keys()
        legVal += legenddict.values()
        if self.usetex: legKey += ["~"]
        else: legKey += [" "]
        legVal += [Rectangle((0, 0), 1, 1, alpha=0.0)]

        # TODO: Fix this problem: ValueError: Can only output finite numbers in PDF
        if self.ext == "pdf" and self.usetex:
            warn("This is a dangerous combination, be prepared to experience problems with the PDF backend. It might help to manually change the number of columns.")
        else:
            table_axis.legend(
              legVal, legKey,
              bbox_to_anchor=tuple(lims),
              bbox_transform=table_axis.transData,
              ncol=len(legVal), loc=9,
              mode="expand", borderaxespad=0.,
              numpoints=1, fontsize='small')
            #table_axis.legend(handles, labels, bbox_to_anchor=(0.0, -0.1), loc=2, ncol=3)

        gs.tight_layout(fig)  # , rect=[0, 0, 1, 0.75])
        # Fine-tune figure; make subplots close to each other
        # and hide x ticks for all but bottom plot.
        # fig.subplots_adjust(wspace=0)
        #plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)

        if json_time < report_time:
            if not quiet: print(" ({0})".format("i"), end="")
        else:
            if not os.path.exists(os.path.dirname(report_path)):
                os.makedirs(os.path.dirname(report_path))

            if self.usebp and self.ext == "pgf":
                self.bp.savepgf(report_path, fig=fig, customReplace=["\\cite{", "\\citet{"])
            else:
                fig.savefig(report_path)

            if not quiet: print(" ({0})".format("w"), end="")

        if not pdfObj is None: pdfObj.savefig(fig)

        plt.close(fig)
        pass

    def makeSolutionPlots(self, solObjs=[SolutionData()], pdfObj=None):
        """
        Creates a whole page with some plots and basic information
        for both fit quality, reference data, data sources and
        more.
        """
        # First we determine some basic settings
        water = None
        solutions = []

        for i in range(len(solObjs) - 1):
            if solObjs[i].xid == SolutionData.ifrac_mass or \
               solObjs[i].xid == SolutionData.ifrac_mole or \
               solObjs[i].xid == SolutionData.ifrac_volume:
                solutions += [solObjs[i]]
            # elif solObjs[i].xid==SolutionData.ifrac_pure:
            #    purefluids += [doneObjs[i]]
            elif solObjs[i].name == "NBS":
                water = solObjs[i]
                solutions += [solObjs[i]]

        if water is None: raise ValueError("No water found, reference values missing.")

        # Set temperature data for all fluids
        dataDict = {}
        dataList = []
        obj = SolutionData()
        for i in range(len(solutions) - 1):
            obj = solutions[i]
            T = np.arange(np.max([275, np.round(obj.Tmin)]), np.min([300, np.round(obj.Tmax)]), 1)
            P = 100e5
            x = obj.xmin
            dataDict["name"] = obj.name
            dataDict["desc"] = obj.description
            dataDict["T"] = T
            dataDict["P"] = P
            dataDict["x"] = x
            if obj.density.type != IncompressibleData.INCOMPRESSIBLE_NOT_SET:
                dataDict["D"] = np.array([obj.rho(Ti, P, x) for Ti in T])
            else: dataDict["D"] = None
            if obj.specific_heat.type != IncompressibleData.INCOMPRESSIBLE_NOT_SET:
                dataDict["C"] = np.array([obj.c(Ti, P, x) for Ti in T])
            else: dataDict["C"] = None
            if obj.conductivity.type != IncompressibleData.INCOMPRESSIBLE_NOT_SET:
                dataDict["L"] = np.array([obj.cond(Ti, P, x) for Ti in T])
            else: dataDict["L"] = None
            if obj.viscosity.type != IncompressibleData.INCOMPRESSIBLE_NOT_SET:
                dataDict["V"] = np.array([obj.visc(Ti, P, x) for Ti in T])
            else: dataDict["V"] = None
            dataList.append(dataDict.copy())

        rat = np.array([5.5, 3.05])
        mul = 2.25
        #fig = plt.figure(figsize=(297/div,210/div))
        fig = plt.figure(figsize=rat * mul)

        ax = fig.add_subplot(121)
        ax.set_xlabel(self.tempLabel)
        ax.set_ylabel(self.densLabel)
        fig.suptitle(r'Aqueous solutions with a concentration of 0.0', fontsize='x-large', fontweight='bold')

        #obj = water
        #T = np.linspace(obj.Tmin, obj.Tmax, num=int(obj.Tmax-obj.Tmin))
        # D =

        #import matplotlib.pyplot as plt
        #from itertools import cycle
        lines = ["-", "--", "-.", ":"]
        colours = ['r', 'g', 'b', 'c', 'm']

        linecycler = itertools.cycle(lines)
        colourcycler = itertools.cycle(colours)

        for i in range(len(dataList) - 1):
            obj = dataList[i]
            if not obj["T"] is None and not obj["D"] is None:
                if obj["x"] == 0.0 and np.any(np.isfinite(obj["D"])) and obj["name"] != "NBS":
                    lc = next(colourcycler)
                    if lc == colours[0]: ls = next(linecycler)
                    ax.plot(obj["T"], obj["D"], label="{0}: {1}".format(obj["name"], obj["desc"]), ls=ls, color=lc)
                # if not np.any(np.isfinite(obj["D"])):
                #    print("Name: {0}, Dmin: {1}, Dmax: {1}".format(obj["name"],np.min(obj["D"]),np.max(obj["D"])))

        # Skip pure water
        #obj = (item for item in dataList if item["name"] == "NBS").next()
        #ax.plot(obj["T"], obj["D"], label="{0}: {1}".format(obj["name"], obj["desc"]), ls='-', color='black')

        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=1)  # , prop={'size':'smaller'})
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
            if i == 0 or i == len(grid) - 1: header_flag = True
            rst += self.normalize_row(row, max_cols)
            rst += self.table_div(max_cols, header_flag)
        return rst

    def table_div(self, max_cols, header_flag=1, indent=2):
        out = u""
        for i in range(indent):
            out += u" "
        if header_flag == 1:
            style = u"="
        else:
            style = u"-"
        for max_col in max_cols:
            out += max_col * style + u" "
        out += u"\n"
        return out

    def normalize_row(self, row, max_cols, indent=2):
        r = u""
        for i in range(indent):
            r += u" "

        for i, max_col in enumerate(max_cols):
            r += row[i] + (max_col - len(row[i]) + 1) * u" "
        return r + u"\n"

    def writeTextToFile(self, path, text):
        #print("Writing to file: {0}".format(path))
        if not os.path.exists(os.path.dirname(path)):
            os.makedirs(os.path.dirname(path))
        with open(path, 'w') as f:
            f.write(text)

        return True

    def writeTxtTableToFile(self, path, table, head=u""):
        if not head == u"":
            return self.writeTextToFile(path + ".txt", head + u"\n\n" + self.make_table(table))
        return self.writeTextToFile(path + ".txt", self.make_table(table))

    def writeCsvTableToFile(self, path, table):
        if not os.path.exists(os.path.dirname(path + ".csv")):
            os.makedirs(os.path.dirname(path + ".csv"))
        with codecs.open(path + ".csv", 'w', encoding='utf-8') as f:
            writer = csv.writer(f)
            # writer = UnicodeWriter(f)
            writer.writerows(table)
        return True

    # Interface
    def writeTableToFile(self, path, table):
        self.writeCsvTableToFile(path, table)
        # self.writeTxtTableToFile(path,table)
        return True

    def getReportLink(self, name):
        reportFile = os.path.join("..", "_static", "fluid_properties", "Incompressibles_reports", "{0}_fitreport.{1}".format(name, self.ext))
        return self.d(name, reportFile)

    def getCitation(self, keys):
        return u":cite:`{0}`".format(keys)

    def checkForNumber(self, number):
        try:
            n = float(number)
        except:
            n = np.NAN
            pass
        return n

    def d(self, text, target):
#         try:
#             if os.path.isfile(target):
#                 link = ":download:`{0}<{1}>`".format(text,target)
#             else:
#                 link = "{0}".format(text)
#         except:
#             link = "{0}".format(text)
#             pass
        # TODO: Fix this!
        link = u":download:`{0}<{1}>`".format(text, target)
        return link

    def m(self, math):
        text = u":math:`{0}`".format(math)
        return text

    def c(self, number):
        #text = "{0:5.2f} |degC|".format(self.checkForNumber(number)-273.15)
        text = u"{0:5.2f}".format(self.checkForNumber(number) - 273.15)
        return text

    def x(self, number):
        text = u"{0:3.2f}".format(self.checkForNumber(number))
        return text

    def generateTexTable(self, solObjs=[SolutionData()], path="table"):

        solObjs = sorted(solObjs, key=lambda x: x.name)
        xmin = np.array([x.xmin for x in solObjs])
        xmax = np.array([x.xmax for x in solObjs])

        # TODO: This is a dirty hack!
        if np.any(xmin > 0.0) and np.any(xmax < 1.0): use_x = True
        else: use_x = False

        header = [r'Name', r'Description', r'Reference', r'{$T_\text{min}$ (\si{\celsius})}', r'{$T_\text{max}$ (\si{\celsius})}', r'{$T_\text{base}$ (\si{\kelvin})}']
        if use_x: header.extend([r'{$x_\text{min}$}', r'{$x_\text{max}$}'])

        testTable = []
        testTable.append(header)  # Headline

        testSection = []

        figPath = r"appendices/IncompressibleFluidsReports/"

#         def getFigureText(fld):
#             text  = "\\begin{pgffigure} \n"
#             text += "\\fbox{\\resizebox{\\textwidth}{!}{ \n"
#             text += "\\subimport{"+str(figPath)+"}{"+str(fld)+"_fitreport.pgf} \n"
#             text += "}} \n"
#             text += "\\end{pgffigure} \n"
#             return text
        def getFigureText(fld):
            text = "{\\centering"
            text += "\\resizebox{\\textwidth}{!}{ \n"
            text += "\\subimport{" + str(figPath) + "}{" + str(fld) + "_fitreport.pgf} \n"
            text += "}} \\vfill \n"
            return text

        def getFigureLabel(fld):
            return "subsec:fit" + str(fld)

        for fluid in solObjs:
            # \hyperref[label_name]{''link text''}
            testTable.append([
                fluid.name + ", p.~\\pageref{" + getFigureLabel(fluid.name) + "}",
                fluid.description.encode("latex"),
                r'\cite{' + str(fluid.reference) + '}',
                self.c(fluid.Tmin),
                self.c(fluid.Tmax),
                self.c(fluid.Tbase + 273.15)
            ])
            if use_x: testTable[-1].extend([self.x(fluid.xmin), self.x(fluid.xmax)])

            testSection.append([
                "\\subsection{Fitted functions for " + str(fluid.name) + "}",
                "\\label{" + getFigureLabel(fluid.name) + "}",
                getFigureText(fluid.name)
            ])

        text = "\\cr \\topcrule \n"
        i = -1
        for row in testTable:
            for i__ in range(len(row)):
                try: row[i__] = str(row[i__], encoding ="utf-8")
                except: pass
            i += 1
            if i < 1:
                tmp = r"\headcol " + " & ".join(row)
            elif i % 2 == 0:
                tmp = r"\rowcol " + " & ".join(row)
            else:
                tmp = " & ".join(row)
            text += tmp + " \\\\\n"
            if i == 0: text += "\\midcrule \n"

        if i % 2 == 0:
            text += "\\bottomcrule \\cr \n"
        else:
            text += "\\bottomrule \\cr \n"

        with open(path + ".tex", 'w') as f:
            f.write(text)


        text = "\n"
        for i, row in enumerate(testSection):
            for i__ in range(len(row)):
                try: row[i__] = str(row[i__], encoding ="utf-8")
                except: pass
            tmp = "\n".join(row)
            text += tmp + " \n\n"

        with open(path + "-section.tex", 'w') as f:
            f.write(text)

        return True

    def generateStatsTable(self, objLists=[[SolutionData()]], labelList=[]):

        if len(objLists) != len(labelList):
            raise ValueError("Wrong length")

        header = [u'Property']
        header.extend(labelList)

        testTable = []
        testTable.append(header)  # Headline

        def getEntries(solObj):
            data = {}
            data["name"] = dict(type=solObj.name, nrms=solObj.name)

            kt = [("density", solObj.density),
              ("specific_heat", solObj.specific_heat),
              ("conductivity", solObj.conductivity),
              ("viscosity", solObj.viscosity),
              ("saturation_pressure", solObj.saturation_pressure),
              ("T_freeze", solObj.T_freeze)]

            for k, _ in kt: data[k] = {}

            for k, t in kt:
                try:
                    data[k]["type"] = t.type
                except:
                    data[k]["type"] = None
                try:
                    data[k]["nrms"] = t.NRMS
                except:
                    data[k]["nrms"] = None

            return data

        errcolumns = []
        typcolumns = []
        for group in objLists:
            errcolumn = {}
            typcolumn = {}
            for obj in group:
                dat = getEntries(obj)
                for k in dat:
                    ent = errcolumn.get(k, [])
                    ent.append(dat[k]["nrms"])
                    errcolumn[k] = ent
                    ent = typcolumn.get(k, [])
                    ent.append(dat[k]["type"])
                    typcolumn[k] = ent
            for k in errcolumn:
                if k != "name":
                    for i, _ in enumerate(errcolumn[k]):
                        try:
                            errcolumn[k][i] = float(errcolumn[k][i])
                        except:
                            errcolumn[k][i] = np.NAN
                        # try:
                        #    typcolumn[k][i] = float(typcolumn[k][i])
                        # except:
                        #    typcolumn[k][i] = np.NAN

            errcolumns.append(errcolumn)
            typcolumns.append(typcolumn)

        keys = errcolumns[0].keys()

        errTable = copy.copy(testTable)
        typTable = copy.copy(testTable)

        for k in keys:
            if k != "name":
                # Error lines
                maxLine = []
                avgLine = []
                minLine = []

                minLine.append("{0:25s}".format(k))
                avgLine.append("{0:25s}".format(""))
                maxLine.append("{0:25s}".format(""))

                # Function type lines
                polyLine = []
                expoLine = []
                logeLine = []
                expPLine = []

                polyLine.append("{0:25s}".format(k))
                expoLine.append("{0:25s}".format("exp"))
                logeLine.append("{0:25s}".format("logexp"))
                expPLine.append("{0:25s}".format("exppoly"))

                #print(k+": ", end="")
                for e, t in zip(errcolumns, typcolumns):
                    try: mi = np.nanargmin(e[k])
                    except: mi = 0
                    try: ma = np.nanargmax(e[k])
                    except: ma = 0
                    try: av = np.nanmean(e[k])
                    except: av = 0.0
                    #print("min: {0}({1}), avg: {2}, max: {3}({4})".format(c[k][mi],c["name"][mi],av,c[k][ma],c["name"][ma]),end="")
                    minLine.append("{0:5.3f} ({1:5s})".format(e[k][mi] * 100.0, e["name"][mi]))
                    avgLine.append("{0:5.3f}  {1:5s} ".format(av * 100.0, ""))
                    maxLine.append("{0:5.3f} ({1:5s})".format(e[k][ma] * 100.0, e["name"][ma]))

                    polyLine.append("{0:3d}".format(t[k].count(IncompressibleData.INCOMPRESSIBLE_POLYNOMIAL)))
                    expoLine.append("{0:3d}".format(t[k].count(IncompressibleData.INCOMPRESSIBLE_EXPONENTIAL)))
                    logeLine.append("{0:3d}".format(t[k].count(IncompressibleData.INCOMPRESSIBLE_LOGEXPONENTIAL)))
                    expPLine.append("{0:3d}".format(t[k].count(IncompressibleData.INCOMPRESSIBLE_EXPPOLYNOMIAL)))
                    print(t[k])

                #print(" ")
                errTable.append(minLine)
                errTable.append(avgLine)
                errTable.append(maxLine)
                errTable.append([r"\midrule"])

                typTable.append(polyLine)
                typTable.append(expoLine)
                typTable.append(logeLine)
                typTable.append(expPLine)
                typTable.append([r"\midrule"])

        errString = ""
        for lin in errTable:
            errString += " & ".join(lin) + "\\\\ \n"
        print(errString)

        typString = ""
        for lin in typTable:
            typString += " & ".join(lin) + "\\\\ \n"
        print(typString)

        return True
#
#
#
#         figPath = r"appendices/IncompressibleFluidsReports/"
#
# #         def getFigureText(fld):
# #             text  = "\\begin{pgffigure} \n"
# #             text += "\\fbox{\\resizebox{\\textwidth}{!}{ \n"
# #             text += "\\subimport{"+str(figPath)+"}{"+str(fld)+"_fitreport.pgf} \n"
# #             text += "}} \n"
# #             text += "\\end{pgffigure} \n"
# #             return text
#         def getFigureText(fld):
#             text  = "{\\centering"
#             text += "\\resizebox{\\textwidth}{!}{ \n"
#             text += "\\subimport{"+str(figPath)+"}{"+str(fld)+"_fitreport.pgf} \n"
#             text += "}} \\vfill \n"
#             return text
#
#         def getFigureLabel(fld):
#             return "subsec:fit"+str(fld)
#
#         for fluid in solObjs:
#             #\hyperref[label_name]{''link text''}
#             testTable.append([
#                 fluid.name+", p.~\\pageref{"+getFigureLabel(fluid.name)+"}",
#                 fluid.description.encode("latex"),
#                 r'\cite{'+str(fluid.reference)+'}',
#                 self.c(fluid.Tmin),
#                 self.c(fluid.Tmax),
#                 self.c(fluid.Tbase+273.15)
#             ])
#             if use_x: testTable[-1].extend([self.x(fluid.xmin), self.x(fluid.xmax)])
#
#             testSection.append([
#                 "\\subsection{Fitted functions for "+str(fluid.name)+"}",
#                 "\\label{"+getFigureLabel(fluid.name)+"}",
#                 getFigureText(fluid.name)
#             ])
#
#         text = "\\cr \\topcrule \n"
#         i = -1
#         for row in testTable:
#             i += 1
#             if i < 1:
#                 tmp   = r"\headcol "+" & ".join(row)
#             elif i % 2 == 0:
#                 tmp   = r"\rowcol " +" & ".join(row)
#             else:
#                 tmp   =              " & ".join(row)
#             text += tmp + " \\\\\n"
#             if i == 0: text += "\\midcrule \n"
#
#         if i % 2 == 0:
#             text += "\\bottomcrule \\cr \n"
#         else:
#             text += "\\bottomrule \\cr \n"
#
#         with open(path+".tex", 'w') as f:
#             f.write(text.encode('utf-8'))
#
#         text = "\n"
#         for i,row in enumerate(testSection):
#             tmp   = "\n".join(row)
#             text += tmp + " \n\n"
#
#         with open(path+"-section.tex", 'w') as f:
#             f.write(text.encode('utf-8'))
#
#         return True

    def generateRstTable(self, solObjs=[SolutionData()], path="table"):

        solObjs = sorted(solObjs, key=lambda x: x.name)
        xmin = np.array([x.xmin for x in solObjs])
        xmax = np.array([x.xmax for x in solObjs])

        # TODO: This is a dirty hack!
        if np.any(xmin > 0.0) and np.any(xmax < 1.0): use_x = True
        else: use_x = False

        header = [u'Name', u'Description', u'Reference', \
          self.m(u'T_\\text{min}') + u" (°C)", self.m(u'T_\\text{max}') + u" (°C)", self.m(u'T_\\text{base}') + u" (K)"]
        if use_x: header.extend([self.m(u'x_\\text{min}'), self.m(u'x_\\text{max}')])

        testTable = []
        testTable.append(header)  # Headline
        for fluid in solObjs:
            testTable.append([
                self.getReportLink(fluid.name),
                fluid.description,
                self.getCitation(fluid.reference),
                self.c(fluid.Tmin),
                self.c(fluid.Tmax),
                self.c(fluid.Tbase + 273.15)
            ])
            if use_x: testTable[-1].extend([self.x(fluid.xmin), self.x(fluid.xmax)])

        self.writeTableToFile(path, testTable)
        return True
