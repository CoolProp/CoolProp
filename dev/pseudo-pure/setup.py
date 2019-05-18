from distutils.core import setup
from Cython.Build import cythonize

import sys, numpy
import Cython
Cython.Compiler.Options.annotate = True
sys.argv += ['build_ext', '--inplace']

setup(
    name="",
    ext_modules=cythonize('summer.pyx'),  # accepts a glob pattern
    include_dirs=[numpy.get_include()]
)
