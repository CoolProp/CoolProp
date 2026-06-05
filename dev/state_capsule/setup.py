"""Build the frozen `State` compat shim.

Needs only the C-ABI struct header (CoolProp/detail/state_capi.h) at compile
time -- it is link-free (the function table is fetched from CoolProp._capi at
runtime).  Run: python setup.py build_ext --inplace
"""
import os

from setuptools import Extension, setup
from Cython.Build import cythonize

HERE = os.path.dirname(os.path.abspath(__file__))
INCLUDE = os.path.abspath(os.path.join(HERE, "..", "..", "include"))  # worktree include/

exts = [
    Extension("State", ["State.pyx"], language="c++", include_dirs=[INCLUDE]),
    # cimport_check cimports State.pxd -> proves cdef-level cimport compatibility
    Extension("cimport_check", ["cimport_check.pyx"], language="c++", include_dirs=[INCLUDE]),
]
setup(ext_modules=cythonize(exts, language_level=3))
