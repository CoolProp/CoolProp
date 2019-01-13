"""
CPIncomp - Jorrit's collection of routines for fitting incompressible liquids
=====

readme.md       - General instructions and copyright information / credits.

"""
from __future__ import division, absolute_import, print_function
import inspect
from . import DataObjects, ExampleObjects, PureFluids, CoefficientFluids, DigitalFluids, MelinderFluids, SolutionFluids
from .SecCoolFluids import SecCoolSolutionData


def getBaseClassNames():
    """
    Returns a list of names of the abstract base
    classes that should not be instantiated. Can
    be used to build an ignore list.
    """
    ignList = []
    for i in inspect.getmembers(DataObjects):
        if inspect.isclass(i[1]):
            ignList.append(i[0])
    return ignList


def getExampleNames(obj=False):
    """
    Returns a list of names of the example fluid
    classes that should not be treated like the
    normal fluids. Can be used to build an ignore
    list.
    Use the obj switch to return the objects instead
    of the names.
    """
    ignList = getBaseClassNames()
    outList = []
    for i in inspect.getmembers(ExampleObjects):
        if inspect.isclass(i[1]):
            if i[0] not in ignList:
                if obj: outList.append(i[1]())
                else: outList.append(i[0])
    return outList


def getIgnoreNames():
    """
    Returns a list of names of classes that should
    not be treated like the normal fluids.
    """
    ignList = []
    ignList += getBaseClassNames()
    ignList += getExampleNames()
    return ignList


def getCoefficientFluids():
    """
    Returns a list of CoefficientData objects, which
    already contain all coefficients. These objects
    can be written to JSON without further processing.
    """
    classes = []
    ignList = getIgnoreNames()
    for name, obj in inspect.getmembers(CoefficientFluids):
        if inspect.isclass(obj):
            # print(name)
            if not name in ignList:  # Ignore the base classes
                classes.append(obj())
    return classes


def getDigitalFluids():
    """
    Returns a list of DigitalData objects, which
    contain data for the fitting. These objects
    only hold the data and you still have to call
    the fitting routines.
    a) Data in these classes is based on equations
    that cannot be converted to the forms available
    in CoolProp.
    b) There are no equations available, but the
    fluid data can be accessed from python in the
    form of another library.
    c) There are data files that contain the experimental
    data. The fit has to be done after loading the fluids.
    """
    classes = []
    ignList = getIgnoreNames()
    for name, obj in inspect.getmembers(DigitalFluids):
        if inspect.isclass(obj):
            # print(name)
            if not name in ignList:  # Ignore the base classes
                classes.append(obj())
    return classes


def getMelinderFluids():
    """
    Returns a list of CoefficientData objects, which
    already contain all coefficients. These objects
    can be written to JSON without further processing.
    All coefficients are taken from the same reference:
    "Properties of Secondary Working Fluids for Indirect Systems"
    written by Aake Melinder and published in 2010 by IIR
    """
    classes = []
    ignList = getIgnoreNames()
    for name, obj in inspect.getmembers(MelinderFluids):
        if inspect.isclass(obj):
            # print(name)
            if not name in ignList:  # Ignore the base classes
                classes.append(obj())
    return classes


def getPureFluids():
    """
    Returns a list of SolutionData objects, which
    contain data for fitting pure fluids. These
    objects only hold the data and you still have
    to call the fitting routines.
    """
    classes = []
    ignList = getIgnoreNames()
    for name, obj in inspect.getmembers(PureFluids):
        if inspect.isclass(obj):
            # print(name)
            if not name in ignList:  # Ignore the base classes
                classes.append(obj())
    return classes


def getSolutionFluids():
    """
    Returns a list of SolutionData objects, which
    contain data for fitting solutions. These
    objects only hold the data and you still have
    to call the fitting routines.
    """
    classes = []
    ignList = getIgnoreNames()
    for name, obj in inspect.getmembers(SolutionFluids):
        if inspect.isclass(obj):
            # print(name)
            if not name in ignList:  # Ignore the base classes
                classes.append(obj())
    return classes


def getSecCoolFluids():
    """
    Returns a list of DigitalData objects, which
    contain data for the fits. All objects here
    implement the fitFluid() function, which can
    be called to set the coefficients before writing
    the JSON files.
    """
    return SecCoolSolutionData.factory()


def get_version():
    return 0.5


if __name__ == "__main__":
    print('You are using version %s of the Python package for incompressible liquids in CoolProp.' % (get_version()))
    print()
