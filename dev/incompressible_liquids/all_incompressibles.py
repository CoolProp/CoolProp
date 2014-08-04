from __future__ import division, absolute_import, print_function
from CPIncomp.WriterObjects import SolutionDataWriter
from CPIncomp import getExampleNames, getPureFluids, getCoefficientFluids,\
    getDigitalFluids, getSecCoolFluids, getMelinderFluids
import sys
from CPIncomp.DataObjects import SolutionData


if __name__ == '__main__':   
    
    runTest     = False
    runExamples = False
    
    runCoeffs   = True #Processing fluids with given coefficients
    runDigital  = True #Processing digital fluids
    runMelinder = True #Processing Melinder fluids
    runPure     = True #Processing pure fluids
    runSecCool  = True #Processing SecCool fluids
    
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
    if runExamples:
        fluidObjs = getExampleNames(obj=True)
        examplesToFit = ["ExamplePure","ExampleSolution","ExampleDigital","ExampleDigitalPure"]
        
        print("\nProcessing example fluids")
        for obj in fluidObjs:
            if obj.name in examplesToFit:
                writer.fitAll(obj)
        doneObjs += fluidObjs[:]       
        writer.writeFluidList(doneObjs)
        writer.writeReportList(doneObjs, pdfFile="all_examples.pdf")
        
    # If the examples did not cause any errors, 
    # we can proceed to the real data.
    doneObjs = []
    
    if runCoeffs:
        print("\nProcessing fluids with given coefficients")
        fluidObjs = getCoefficientFluids()
        doneObjs += fluidObjs[:]
    
    if runDigital:
        print("\nProcessing digital fluids")
        fluidObjs = getDigitalFluids()
        writer.fitFluidList(fluidObjs)
        doneObjs += fluidObjs[:]
    
    if runMelinder:
        print("\nProcessing Melinder fluids")
        fluidObjs = getMelinderFluids()
        doneObjs += fluidObjs[:]
    
    if runPure:
        print("\nProcessing pure fluids")
        fluidObjs = getPureFluids()
        writer.fitFluidList(fluidObjs)
        doneObjs += fluidObjs[:]
    
    if runSecCool:
        print("\nProcessing SecCool fluids")
        fluidObjs = getSecCoolFluids()
        writer.fitSecCoolList(fluidObjs)
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
        
    print("All checks passed, going to write parameters to disk.")
    writer.writeFluidList(doneObjs)
    
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
       
    writer.makeSolutionPlots(solObjs=doneObjs, pdfObj=None)

    print("All done, bye")
    sys.exit(0)

    