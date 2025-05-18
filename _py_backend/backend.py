from setuptools import build_meta as _orig
from setuptools.build_meta import *
import os 
import sys
from pathlib import Path

#

def build_sdist(wheel_directory, config_settings=None):
    return  _orig.build_sdist(wheel_directory, config_settings)
    
def build_wheel(wheel_directory, config_settings=None, metadata_directory=None):
    sys.path.append(str(Path(__file__).parent))
    os.chdir('wrappers/Python')
    return  _orig.build_wheel(wheel_directory, config_settings, metadata_directory)
