"""Driver for the incompressible-fluid fitting pipeline.

Fits all fluid definitions from the CPIncomp package and writes the
per-fluid coefficient files to json/. Optionally also generates fitting
reports (PDF), summary figures and the RST/TeX property tables used by the
online documentation.

Fitting and JSON output need only numpy and scipy. The reports, figures and
tables additionally need matplotlib and a built CoolProp Python package --
skip them with:

    python all_incompressibles.py -nr -ns -nt -nst

See README.md in this directory for how to add a new fluid.
"""
import argparse
import os

import numpy as np

from CPIncomp import getExampleNames, getPureFluids, getSolutionFluids, getCoefficientFluids, getDigitalFluids, getSecCoolFluids, getMelinderFluids
from CPIncomp.DataObjects import SolutionData
from CPIncomp.WriterObjects import SolutionDataWriter


def getTime(path):
    if os.path.isfile(path):
        return os.path.getctime(path)
    return 0


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-nf", "--nofit", action='store_true', help="Do not fit the data, but read the JSON files")
    parser.add_argument("-nr", "--noreports", action='store_true', help="Do not write the fitting reports")
    parser.add_argument("-ns", "--nosummary", action='store_true', help="Do not generate the summary figures")
    parser.add_argument("-nt", "--notables", action='store_true', help="Do not write the fluid tables")
    parser.add_argument("-nst", "--nostats", action='store_true', help="Do not process statistical parameters")
    args = parser.parse_args()

    runFitting = not args.nofit
    runReports = not args.noreports
    runSummary = not args.nosummary
    runTables = not args.notables
    runStats = not args.nostats

    print("")
    print("Processing the incompressible fluids for CoolProp")
    print("Legend: FluidName (w) | (i) -> (w)=written, (i)=ignored, unchanged coefficient or reports")
    print("")

    writer = SolutionDataWriter()

    # Process the examples first: they are small and fail fast if the
    # pipeline itself is broken.
    print("\nProcessing example fluids")
    doneObjs = getExampleNames(obj=True)
    examplesToFit = ["ExamplePure", "ExampleSolution", "ExampleDigital", "ExampleDigitalPure"]
    for obj in doneObjs:
        if obj.name in examplesToFit:
            if runFitting: writer.fitAll(obj)
            else: writer.fromJSON(obj)
    if runFitting: writer.writeFluidList(doneObjs)
    if runReports:
        if writer.usetex: combined_name = None
        else: combined_name = os.path.join(os.path.abspath("report"), "all_examples.pdf")
        writer.writeReportList(doneObjs, pdfFile=combined_name)

    # If the examples did not cause any errors,
    # we can proceed to the real data.
    doneObjs = []

    print("\nProcessing fluids with given coefficients")
    doneObjs += getCoefficientFluids()

    print("\nProcessing digital fluids")
    fluidObjs = getDigitalFluids()
    if runFitting: writer.fitFluidList(fluidObjs)
    else: writer.readFluidList(fluidObjs)
    doneObjs += fluidObjs

    print("\nProcessing Melinder fluids")
    doneObjs += getMelinderFluids()

    print("\nProcessing pure fluids")
    fluidObjs = getPureFluids()
    if runFitting: writer.fitFluidList(fluidObjs)
    else: writer.readFluidList(fluidObjs)
    doneObjs += fluidObjs

    print("\nProcessing solutions")
    fluidObjs = getSolutionFluids()
    if runFitting: writer.fitFluidList(fluidObjs)
    else: writer.readFluidList(fluidObjs)
    doneObjs += fluidObjs

    print("\nProcessing SecCool fluids")
    fluidObjs = getSecCoolFluids()
    if runFitting: writer.fitSecCoolList(fluidObjs)
    else: writer.readFluidList(fluidObjs)
    doneObjs += fluidObjs

    print("\nAll {0} fluids processed, all coefficients should be set.".format(len(doneObjs)))
    print("Checking the list of fluid objects.")
    doneObjs = sorted(doneObjs, key=lambda x: x.name)

    names = [obj.name for obj in doneObjs]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise ValueError("Two fluids have the same name, that does not work: {0}".format(duplicates))

    purefluids = [obj for obj in doneObjs if obj.xid == SolutionData.ifrac_pure]
    solMass = [obj for obj in doneObjs if obj.xid == SolutionData.ifrac_mass]
    solMole = [obj for obj in doneObjs if obj.xid == SolutionData.ifrac_mole]
    solVolu = [obj for obj in doneObjs if obj.xid == SolutionData.ifrac_volume]
    known = purefluids + solMass + solMole + solVolu
    errors = [obj for obj in doneObjs if obj not in known]
    if errors:
        raise ValueError("There was a problem processing the fluid(s): {0}".format([error.name for error in errors]))

    if runFitting:
        print("All checks passed, going to write parameters to disk.")
        writer.writeFluidList(doneObjs)

    if runReports:
        if writer.usetex:
            combined_name = None
            combined_time = 0
        else:
            combined_name = os.path.join(os.path.abspath("report"), "all_incompressibles.pdf")
            combined_time = getTime(combined_name)

        singles_time = np.array([getTime(writer.get_json_file(fl.name)) for fl in doneObjs])
        if np.any(singles_time > combined_time):
            print("Processing {0:2d} fluids - ".format(len(doneObjs)), end="")
            writer.writeReportList(doneObjs, pdfFile=combined_name)
        else:
            print("No newer file found, aborting.")

    if runSummary:
        writer.makeSolutionPlots(solObjs=doneObjs, pdfObj=None)

    if runTables:
        FLUID_INFO_FOLDER = os.path.abspath("tables")
        objLists = [purefluids, solMass, solMole, solVolu]
        filLists = [os.path.join(FLUID_INFO_FOLDER, "Incompressibles_" + name)
                    for name in ["pure-fluids", "mass-based-fluids", "mole-based-fluids", "volume-based-fluids"]]
        for objList, filList in zip(objLists, filLists):
            writer.generateRstTable(objList, filList)
            writer.generateTexTable(objList, filList)

    if runStats:
        lists = [purefluids, solMass, solVolu]
        labels = ["Pure", "Mass", "Volume"]
        writer.generateStatsTable(lists, labels)

    print("All done, bye")
