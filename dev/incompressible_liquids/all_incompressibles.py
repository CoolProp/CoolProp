from __future__ import division, absolute_import, print_function
from CPIncomp.WriterObjects import SolutionDataWriter
from CPIncomp import getExampleNames, getPureFluids, getSolutionFluids, getCoefficientFluids, getDigitalFluids, getSecCoolFluids, getMelinderFluids
import sys
from CPIncomp.DataObjects import SolutionData
import argparse
import os
import numpy as np

def getTime(path):
    if os.path.isfile(path):
        return os.path.getctime(path)
    else:
        return 0

def mergePdfIfNewer(singlePdfs, combined):
    from PyPDF2 import PdfFileMerger, PdfFileReader
    singles_time = np.array([])
    for fi in singlePdfs:
        singles_time = np.append(singles_time, [getTime(fi)])
    combined_time = getTime(combined)
    if np.any(singles_time>combined_time):
        allcombined = PdfFileMerger()
        for fl in singlePdfs:
            allcombined.append(PdfFileReader(fl,"rb"))
        allcombined.write(combined)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-t","--test", action='store_true', help="only run a subset of fluid")
    parser.add_argument("-nf","--nofit", action='store_true', help="Do not fit the data, but read the JSON files")
    parser.add_argument("-nr","--noreports", action='store_true', help="Do not write the fitting reports")
    parser.add_argument("-ns","--nosummary", action='store_true', help="Do not generate the summary figures")
    #parser.add_argument("-f","--fluid", help="Only process the fluid FLUID")
    parser.add_argument("-fr","--fullreport", action='store_true', help="Generate one file with all fitting reports")

    args = parser.parse_args()
#    if args.verbosity:
#     print "verbosity turned on"

    if args.test:      runTest    = True
    else:              runTest    = False
    if args.nofit:     runFitting = False
    else:              runFitting = True
    if args.noreports: runReports = False
    else:              runReports = True
    if args.nosummary: runSummary = False
    else:              runSummary = True
    if args.fullreport:runFullreport = True
    else:              runFullreport = False
    #if args.fluid:     onlyFluid  = args.fluid
    #else:              onlyFluid  = None

    print("")
    print("Processing the incompressible fluids for CoolProp")
    print("Legend: FluidName (w) | (i) -> (w)=written, (i)=ignored, unchanged coefficient or reports")
    print("")

    writer = SolutionDataWriter()
    doneObjs = []

    # To debug single fluids
    if runTest:
        solObjs = []
        from CPIncomp.SecCoolFluids import SecCoolSolutionData,SecCoolIceData,ThermogenVP1869
        from CPIncomp.PureFluids import Texatherm22
        solObjs += [SecCoolSolutionData(sFile='Melinder, Ammonia'            ,sFolder='xMass',name='MAM2',desc='Melinder, Ammonia'            ,ref='Melinder-BOOK-2010, SecCool software')]
        solObjs += [SecCoolIceData(sFile='IceNA'   ,sFolder='xMass',name='IceNA',desc='Ice slurry with NaCl' ,ref='Danish Technological Institute, SecCool software')]
        #solObjs = [Freezium()]
        #solObjs[0].density.DEBUG = True
        #solObjs[0].specific_heat.DEBUG = True
        #solObjs[0].conductivity.DEBUG = True
        #solObjs[0].viscosity.DEBUG = True
        #solObjs[0].T_freeze.DEBUG = True
        #writer.fitSecCoolList(solObjs)
        solObjs = [ThermogenVP1869()]#,Therminol72()]
        solObjs[0].viscosity.DEBUG=True
        #solObjs[0].saturation_pressure.DEBUG=True
        #
        ##from CPIncomp.ExampleObjects import SecCoolExample
        ##solObjs = [SecCoolExample()]
        #writer.fitFluidList(solObjs)
        writer.fitSecCoolList(solObjs)
        writer.writeFluidList(solObjs)
        writer.writeReportList(solObjs)
        sys.exit(0)

    # treat the examples first
    fluidObjs = getExampleNames(obj=True)
    examplesToFit = ["ExamplePure","ExampleSolution","ExampleDigital","ExampleDigitalPure"]

    print("\nProcessing example fluids")
    for obj in fluidObjs:
        if obj.name in examplesToFit:
            if runFitting: writer.fitAll(obj)
            else: writer.fromJSON(obj)
    doneObjs += fluidObjs[:]
    if runFitting: writer.writeFluidList(doneObjs)
    if runReports:
        # TODO: The new method for multipage PDFs produces larger files, why?
        writer.writeReportList(doneObjs, pdfFile="all_examples.pdf")
        #singleNames = [writer.get_report_file(fl.name) for fl in doneObjs]
        #mergePdfIfNewer(singleNames, "all_examples.pdf")


    # If the examples did not cause any errors,
    # we can proceed to the real data.
    doneObjs = []

    print("\nProcessing fluids with given coefficients")
    fluidObjs = getCoefficientFluids()
    doneObjs += fluidObjs[:]

    print("\nProcessing digital fluids")
    fluidObjs = getDigitalFluids()
    if runFitting: writer.fitFluidList(fluidObjs)
    else: writer.readFluidList(fluidObjs)
    doneObjs += fluidObjs[:]

    print("\nProcessing Melinder fluids")
    fluidObjs = getMelinderFluids()
    doneObjs += fluidObjs[:]

    print("\nProcessing pure fluids")
    fluidObjs = getPureFluids()
    if runFitting: writer.fitFluidList(fluidObjs)
    else: writer.readFluidList(fluidObjs)
    doneObjs += fluidObjs[:]

    print("\nProcessing solutions")
    fluidObjs = getSolutionFluids()
    if runFitting: writer.fitFluidList(fluidObjs)
    else: writer.readFluidList(fluidObjs)
    doneObjs += fluidObjs[:]

    print("\nProcessing SecCool fluids")
    fluidObjs = getSecCoolFluids()
    if runFitting: writer.fitSecCoolList(fluidObjs)
    else: writer.readFluidList(fluidObjs)
    doneObjs += fluidObjs[:]

    print("\nAll {0} fluids processed, all coefficients should be set.".format(len(doneObjs)))
    print("Checking the list of fluid objects.")
    #doneObjs += getCoefficientObjects()[:]
    doneObjs = sorted(doneObjs, key=lambda x: x.name)

    purefluids = []
    solMass    = []
    solMole    = []
    solVolu    = []
    errors     = []



    for i in range(0,len(doneObjs)):

        curObj = doneObjs[i]

        if i < len(doneObjs)-2:
            nexObj = doneObjs[i+1]
        else:
            nexObj = doneObjs[0]

        if curObj.name==nexObj.name:
            print("Conflict between {0} and {1}, aborting".format(curObj,nexObj))
            raise ValueError("Two elements have the same name, that does not work: {0}".format(curObj.name))
        else:
            if curObj.xid==SolutionData.ifrac_mass:
                solMass    += [curObj]
            elif curObj.xid==SolutionData.ifrac_mole:
                solMole    += [curObj]
            elif curObj.xid==SolutionData.ifrac_volume:
                solVolu    += [curObj]
            elif curObj.xid==SolutionData.ifrac_pure:
                purefluids += [curObj]
            else:
                errors += [curObj]

    if len(errors)>0:
        raise ValueError("There was a problem processing the fluid(s): {0}".format([error.name for error in errors]))

    solutions  = solMass
    solutions += solMole
    solutions += solVolu

    if runFitting: print("All checks passed, going to write parameters to disk.")
    if runFitting: writer.writeFluidList(doneObjs)

    if runReports:
        print("Creating the fitting reports for the different groups.")
        #writer.writeReportList(doneObjs)
        #doneObjs.sort(key=lambda x: (x.xid ,x.name))
        if len(purefluids)>0:
            print("Processing {0:2d} pure fluids   - ".format(len(purefluids)), end="")
            writer.writeReportList(purefluids)#, pdfFile="all_pure.pdf")
            #singleNames = [writer.get_report_file(fl.name) for fl in purefluids]
            #mergePdfIfNewer(singleNames, "all_pure.pdf")

        if len(solutions)>0:
            print("Processing {0:2d} solutions     - ".format(len(solutions)), end="")
            writer.writeReportList(solutions)#, pdfFile="all_solutions.pdf")
            #singleNames = [writer.get_report_file(fl.name) for fl in solutions]
            #mergePdfIfNewer(singleNames, "all_solutions.pdf")

    if runFullreport:
        print("Creating the full fitting report.")

        combined_name = "all_incompressibles.pdf"

        if len(doneObjs)>0:
            singles_time = np.array([])
            for fl in doneObjs:
                singles_time = np.append(singles_time, [getTime(writer.get_report_file(fl.name))])
            combined_time = getTime(combined_name)
            if np.any(singles_time>combined_time):
                print("Processing {0:2d} fluids        - ".format(len(doneObjs)), end="")
                writer.writeReportList(doneObjs, pdfFile=combined_name)

    if runSummary:
        writer.makeSolutionPlots(solObjs=doneObjs, pdfObj=None)

    print("All done, bye")
    sys.exit(0)

