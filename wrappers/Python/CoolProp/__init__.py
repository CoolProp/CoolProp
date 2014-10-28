from __future__ import absolute_import

try:
    import matplotlib
except ImportError:
    raise ImportError('It was not possible to import the required package matplotlib.  Please install it using your package manager')

try:
    import numpy
except ImportError:
    raise ImportError('It was not possible to import the required package numpy.  Please install it using your package manager')

from .CoolProp import AbstractState
from . import CoolProp
from . import HumidAirProp
from . import State
from .constants import *

__fluids__ = CoolProp.get_global_param_string('fluids_list').split(',')
__incompressibles_pure__ = CoolProp.get_global_param_string('incompressible_list_pure').split(',')
__incompressibles_solution__ = CoolProp.get_global_param_string('incompressible_list_solution').split(',')
__version__ = CoolProp.get_global_param_string('version')
__gitrevision__ = CoolProp.get_global_param_string('gitrevision')

def get(s):
    """
    This is just a shorthand function for getting a parameter from
    ``CoolProp.get_global_param_string``
    """
    return CoolProp.get_global_param_string(s)
def test():
    """
    Run the tests in the test folder
    """
    from .tests import runner
    runner.run()

def get_include_directory():
    """
    Get the include directory for CoolProp header files that are needed if you want
    to compile anything else that uses the CoolProp Cython extension type
    
    Returns
    -------
    include_directory: The path to the include folder for CoolProp
    """
    import os
    head, file = os.path.split(__file__)
    return os.path.join(head, 'include')
    
def copy_BibTeX_library(file = None, folder = None):
    """
    Copy the CoolProp BibTeX library file to the file given by ``file``, or the folder given by ``folder``
    
    If no inputs are provided, the file will be copied to the current working 
    directory
    
    Parameters
    ----------
    file : string
        Provide if you want to put the file into a given file
    folder : string
        Provide if you want to put the CoolPropBibTeXLibrary.bib file into the given folder
    
    """
    import os, shutil
    path_to_bib = os.path.join(os.path.split(__file__)[0],'CoolPropBibTeXLibrary.bib')
    if file is None and folder is None:
        shutil.copy2(path_to_bib,os.path.abspath(os.curdir))
    elif file and folder is None:
        shutil.copy2(path_to_bib,file)
    elif folder and file is None:
        shutil.copy2(path_to_bib,os.path.join(folder,file))
    else:
        raise ValueError('can only provide one of file or folder')
        
    