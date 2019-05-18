
from distutils.core import setup, Extension
import subprocess, shutil, os, sys

sys.argv += ['build_ext', '--inplace', '--reswig']

if '--reswig' in sys.argv:
    import subprocess
    subprocess.check_output(['swig', '-python', '-outcurrentdir', '-c++', '-I../../../CoolProp', 'Helmholtz.i'])
    sys.argv.remove('--reswig')

commons = dict()

helm_module = Extension('_helmholtz',
                        sources=['Helmholtz_wrap.cxx', '../../../CoolProp/Helmholtz.cpp', '../../../CoolProp/CoolPropTools.cpp'],
                        **commons
                        )

setup(name='EOSTerms',
       version='0.0.0',
       author="Ian Bell",
       author_email='ian.h.bell@gmail.com',
       description=""" helmholtz energy terms for EOS fitting """,
       ext_modules=[helm_module],
       )
