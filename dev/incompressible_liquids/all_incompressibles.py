from __future__ import division, absolute_import, print_function
from CPIncomp.WriterObjects import SolutionDataWriter
from CPIncomp import getExampleNames, getPureFluids, getCoefficientFluids,\
    getDigitalFluids, getSecCoolFluids, getMelinderFluids
import sys
from CPIncomp.DataObjects import SolutionData


if __name__ == '__main__':   
    
    runTest     = False    
    runFitting  = False
    runReports  = False
    runSummary  = True
        
    writer = SolutionDataWriter()
    doneObjs = []
    
    # To debug single fluids
    if runTest:
        solObjs = []
        from CPIncomp.SecCoolFluids import SecCoolSolutionData,SecCoolIceData
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
        solObjs += [Texatherm22()]#,Therminol72()]
        solObjs[0].viscosity.DEBUG=True
        #solObjs[0].saturation_pressure.DEBUG=True
        #
        ##from CPIncomp.ExampleObjects import SecCoolExample
        ##solObjs = [SecCoolExample()]
        writer.fitFluidList(solObjs)
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
    if runReports: writer.writeReportList(doneObjs, pdfFile="all_examples.pdf")
        
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
    for i in range(len(doneObjs)-1):
        if doneObjs[i].name==doneObjs[i+1].name:
            print("Conflict between {0} and {1}, aborting".format(doneObjs[i],doneObjs[i+1]))
            raise ValueError("Two elements have the same name, that does not work: {0}".format(doneObjs[i].name))
        else:
            if doneObjs[i].xid==SolutionData.ifrac_mass:
                solMass    += [doneObjs[i]]
            elif doneObjs[i].xid==SolutionData.ifrac_mole:
                solMole    += [doneObjs[i]]
            elif doneObjs[i].xid==SolutionData.ifrac_volume:
                solVolu    += [doneObjs[i]]
            elif doneObjs[i].xid==SolutionData.ifrac_pure:
                purefluids += [doneObjs[i]]
            else: 
                errors += [doneObjs[i]]
    solutions  = solMass
    solutions += solMole
    solutions += solVolu
        
    if runFitting: print("All checks passed, going to write parameters to disk.")
    if runFitting: writer.writeFluidList(doneObjs)
    
    if runReports: 
        print("Creating the fitting reports for the different groups.")
        #writer.writeReportList(doneObjs)
        #doneObjs.sort(key=lambda x: (x.xid ,x.name))
        if len(purefluids)>0 and runReports:
            print("Processing {0:2d} pure fluids   - ".format(len(purefluids)), end="")
            writer.writeReportList(purefluids, pdfFile="all_pure.pdf")
        if len(solutions)>0 and runReports:
            print("Processing {0:2d} solutions     - ".format(len(solutions)), end="")
            writer.writeReportList(solutions, pdfFile="all_solutions.pdf")  
        if len(errors)>0 and runReports:
            print("Processing {0:2d} faulty fluids - ".format(len(errors)), end="")
            writer.writeReportList(errors, pdfFile="all_errors.pdf")
    
    if runSummary:
        writer.makeSolutionPlots(solObjs=doneObjs, pdfObj=None)

    print("All done, bye")
    sys.exit(0)

    