"""
CPIncomp - fitting routines for the incompressible fluids in CoolProp
=====

Fluids are defined as classes in the fluid modules (PureFluids.py,
SolutionFluids.py, ...) and discovered automatically: every class in those
modules that is not a base class or an example is instantiated and fitted.
Adding a fluid therefore needs no registration step -- see README.md in
dev/incompressible_liquids for a walk-through.
"""
import inspect

from . import DataObjects, ExampleObjects, PureFluids, CoefficientFluids, DigitalFluids, MelinderFluids, SolutionFluids
from .SecCoolFluids import SecCoolSolutionData


def getBaseClassNames():
    """
    Returns a list of names of the abstract base
    classes that should not be instantiated. Can
    be used to build an ignore list.
    """
    return [name for name, obj in inspect.getmembers(DataObjects) if inspect.isclass(obj)]


def getExampleNames(obj=False):
    """
    Returns a list of names of the example fluid
    classes that should not be treated like the
    normal fluids. Can be used to build an ignore
    list.
    Use the obj switch to return the objects instead
    of the names.
    """
    baseNames = getBaseClassNames()
    if obj:
        return [cls() for name, cls in inspect.getmembers(ExampleObjects) if inspect.isclass(cls) and name not in baseNames]
    return [name for name, cls in inspect.getmembers(ExampleObjects) if inspect.isclass(cls) and name not in baseNames]


def getIgnoreNames():
    """
    Returns a list of names of classes that should
    not be treated like the normal fluids.
    """
    return getBaseClassNames() + getExampleNames()


def _instantiateFluids(module):
    """One instance of every fluid class in module, skipping base classes and examples."""
    ignoreNames = getIgnoreNames()
    return [cls() for name, cls in inspect.getmembers(module) if inspect.isclass(cls) and name not in ignoreNames]


def getCoefficientFluids():
    """
    Returns a list of CoefficientData objects, which
    already contain all coefficients. These objects
    can be written to JSON without further processing.
    """
    return _instantiateFluids(CoefficientFluids)


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
    return _instantiateFluids(DigitalFluids)


def getMelinderFluids():
    """
    Returns a list of CoefficientData objects, which
    already contain all coefficients. These objects
    can be written to JSON without further processing.
    All coefficients are taken from the same reference:
    "Properties of Secondary Working Fluids for Indirect Systems"
    written by Aake Melinder and published in 2010 by IIR
    """
    return _instantiateFluids(MelinderFluids)


def getPureFluids():
    """
    Returns a list of SolutionData objects, which
    contain data for fitting pure fluids. These
    objects only hold the data and you still have
    to call the fitting routines.
    """
    return _instantiateFluids(PureFluids)


def getSolutionFluids():
    """
    Returns a list of SolutionData objects, which
    contain data for fitting solutions. These
    objects only hold the data and you still have
    to call the fitting routines.
    """
    return _instantiateFluids(SolutionFluids)


def getSecCoolFluids():
    """
    Returns a list of DigitalData objects, which
    contain data for the fits. All objects here
    implement the fitFluid() function, which can
    be called to set the coefficients before writing
    the JSON files.
    """
    return SecCoolSolutionData.factory()
