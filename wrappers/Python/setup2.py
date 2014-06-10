from __future__ import print_function
#Check for cython >= 0.17 due to the use of STL containers
try:
    import Cython
except ImportError:
    raise ImportError("Cython not found")
major,minor = Cython.__version__.split('.')[0:2]
#Convert major to integer
major = int(major)
iEnd = 0
while minor[iEnd].isdigit():
    iEnd+=1
    if iEnd == len(minor):
        break

minor = int(minor[0:iEnd])
if not(major > 0 or minor >= 17):
    raise ImportError('Cython version >= 0.17 required due to the use of STL wrappers.  Please update your version of cython')

if minor >= 20:
    _profiling_enabled = True
else:
    _profiling_enabled = False
        
from distutils.core import setup, Extension
import subprocess, shutil, os, sys, glob
from Cython.Build import cythonize
from Cython.Distutils import build_ext
from Cython.Distutils.extension import Extension as CyExtension
from distutils.sysconfig import get_python_inc
from distutils.ccompiler import new_compiler 
from distutils.dep_util import newer_group

def find_cpp_sources(root = os.path.join('..','..','src'), extensions = ['.cpp'], skip_files = None):
    file_listing = []
    for path, dirs, files in os.walk(root):
        for file in files:
            n,ext = os.path.splitext(file)
            fname = os.path.relpath(os.path.join(path, file))
            if skip_files is not None and fname in skip_files: continue
            if ext in extensions:
                file_listing.append(fname)
    return file_listing
    
#This will generate HTML to show where there are still pythonic bits hiding out
Cython.Compiler.Options.annotate = True

# Two options for the location of the C++ files.  Either in the normal ../../CoolProp folder
# or they are in the CoolPropSource folder in this folder

# For PYPI purposes
if 'sdist' in sys.argv:
    CProot = '.'
    shutil.copy2(os.path.join('..','..','version.txt'),'version.txt')
    shutil.copytree(os.path.join('..','..','CoolProp'),'CoolPropSource')
else:
    CProot = os.path.join('..','..')
    
if os.path.exists('CoolPropSource'):
    CPSourceDir = 'CoolPropSource'
    CProot = '.'
else:
    CPSourceDir = os.path.join(CProot,'src')
    CPIncludeDir = os.path.join(CProot,'include')
    
version = open(os.path.join(CProot,'version.txt'),'r').read().strip()

if __name__=='__main__':
    
    #This will generate HTML to show where there are still pythonic bits hiding out
    Cython.Compiler.Options.annotate = True

    ## If the file is run directly without any parameters, build and install
    if len(sys.argv)==1:
        #sys.argv += ['build_ext','--inplace']
       sys.argv += ['build','--compiler=mingw32','install']
#         sys.argv += ['build','install']
        #sys.argv += ['install']

    Sources = find_cpp_sources(CPSourceDir,'*.cpp')
    
    ### Include folders for build
    include_dirs = [os.path.join(CPIncludeDir)]

    if _profiling_enabled:
        cython_directives = dict(profile = True,
                                 embed_signature = True)
    else:
        cython_directives = dict(embed_signature = True)
        
    common_args = dict(include_dirs = include_dirs,
                       language='c++',
                       cython_c_in_temp = True,
                       cython_directives = cython_directives
                       )
                        
    AbstractState_module = CyExtension('CoolProp5.AbstractState',
                        [os.path.join('CoolProp5','AbstractState.pyx')]+Sources,
                        **common_args)
                        
    CoolProp_module = CyExtension('CoolProp5.CoolProp',
                        [os.path.join('CoolProp5','CoolProp.pyx')]+Sources,
                        **common_args)
                        
    constants_module = CyExtension('CoolProp5.constants',
                        [os.path.join('CoolProp5','constants.pyx')],
                        **common_args)
                            
#     # Collect all the header files in the main folder into an include folder
#     try:
#         os.mkdir(os.path.join('CoolProp','include'))
#     except:
#         pass

#     for folder in [os.path.join(CPSourceDir)]:
#         for header in glob.glob(os.path.join(folder,'*.h')):
#             pth,fName = os.path.split(header)
#             shutil.copy2(header,os.path.join('CoolProp','include',fName))
    
    #shutil.copy2(os.path.join(CPSourceDir,'CoolPropBibTeXLibrary.bib'),os.path.join('CoolProp','CoolPropBibTeXLibrary.bib'))
    
    setup (name = 'CoolProp5',
           version = version, #look above for the definition of version variable - don't modify it here
           author = "Ian Bell",
           author_email='ian.h.bell@gmail.com',
           url='http://coolprop.sourceforge.net',
           description = """Open-source thermodynamic and transport properties database""",
           packages = ['CoolProp5','CoolProp5.Plots','CoolProp5.tests','CoolProp5.GUI'],
           ext_modules = [CoolProp_module, AbstractState_module, constants_module],
           package_dir = {'CoolProp5':'CoolProp5',},
           #package_data = {'CoolProp5':['State.pxd','CoolProp.pxd','constants_header.pxd','include/*.h','include/rapidjson/*.h','include/rapidjson/internal/*.h','CoolPropBibTeXLibrary.bib']},
           cmdclass={'build_ext': build_ext},
           
           classifiers = [
            "Programming Language :: Python",
            "Development Status :: 4 - Beta",
            "Environment :: Other Environment",
            "Intended Audience :: Developers",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
            "Topic :: Software Development :: Libraries :: Python Modules"
            ],
           )
    
    sys.path.pop(0)
    import CoolProp5.CoolProp as CP5
    print(CP5.PropsSI('L','T',300,'D',1e-10,'Water'))
    
#     #Clean up the include folder
#     shutil.rmtree(os.path.join('CoolProp','include'), ignore_errors = True)
#     os.remove(os.path.join('CoolProp','CoolPropBibTeXLibrary.bib'))
#     
#     for file in glob.glob(os.path.join('CoolProp','__init__.*')):
#         try:
#             os.remove(file)
#         except:
#             pass
#     
#     if 'sdist' in sys.argv:
#         shutil.rmtree('CoolPropSource')
#         os.remove('version.txt')
#     touch('setup.py')
    
#    try:
#        import nose
#        import CoolProp
#        CoolProp.test()
#    except ImportError:
#        print("Could not run tests, nose not installed")
        
        
    
