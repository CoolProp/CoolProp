from __future__ import print_function

def copy_files():
    import shutil
    shutil.copytree(os.path.join(CProot, 'include'), os.path.join('CoolProp5','include'))
    shutil.copy2(os.path.join(CProot, 'CoolPropBibTeXLibrary.bib'), os.path.join('CoolProp5', 'CoolPropBibTeXLibrary.bib'))
    
def remove_files():
    import shutil
    shutil.rmtree(os.path.join('CoolProp5','include'), ignore_errors = True)
    os.remove(os.path.join('CoolProp5', 'CoolPropBibTeXLibrary.bib'))
    
if __name__=='__main__':

    import subprocess, shutil, os, sys, glob
    
    # Check if a sdist build for pypi
    pypi = os.path.exists('.use_this_directory_as_root')
    
    """
    Modes of operation:
    1) Building the source distro (generate_headers.py must have been run before making the repo)
    2) Installing from source (generate_headers.py must have been run before making the repo)
    3) Installing from git repo (need to make sure to run generate_headers.py)
    4) 
    """ 
    
    # Determine whether or not to use Cython - default is to use cython unless the file .build_without_cython is found in the current working directory
    USE_CYTHON = not os.path.exists('.build_without_cython')
    cy_ext = 'pyx' if USE_CYTHON else 'cpp'

    if USE_CYTHON:
        # Check for cython >= 0.21 due to the use of cpdef enum
        try:
            import Cython
        except ImportError:
            raise ImportError("Cython not found, please install it.  You can do a pip install Cython")
            
        from pkg_resources import parse_version
        if parse_version(Cython.__version__) < parse_version('0.20'):
            raise ImportError('Your version of Cython (%s) must be >= 0.20 .  Please update your version of cython' % (Cython.__version__,))

        if parse_version(Cython.__version__) >= parse_version('0.20'):
            _profiling_enabled = True
        else:
            _profiling_enabled = False
            
        if _profiling_enabled:
            cython_directives = dict(profile = True,
                                     embed_signature = True)
        else:
            cython_directives = dict(embed_signature = True)
    else:
        cython_directives = {}

    # Determine the path to the root of the repository, the folder that contains the CMakeLists.txt file 
    # for normal builds, or the main directory for sdist builds
    if pypi:
        CProot = '.'
    else:
        if os.path.exists(os.path.join('..','..','CMakeLists.txt')):
            # Good working directory
            CProot = os.path.join('..','..')
        else:
            raise ValueError('Could not run script from this folder(' + os.path.abspath(os.path.curdir) + '). Run from wrappers/Python folder')
    
        # Generate the headers - does nothing if up to date - but only if not pypi
        subprocess.check_call('python generate_headers.py', 
                              shell = True, 
                              stdout = sys.stdout, 
                              cwd = os.path.join(CProot, 'dev')
                              )
    
                    
    # Read the version from a bare string stored in file in root directory
    version = open(os.path.join(CProot,'.version'),'r').read().strip()

    setup_kwargs = {}
    from setuptools import setup, Extension, find_packages
    if USE_CYTHON:
        import Cython.Compiler
        from Cython.Distutils.extension import Extension
        from Cython.Build import cythonize
        from Cython.Distutils import build_ext
        
        # This will always generate HTML to show where there are still pythonic bits hiding out
        Cython.Compiler.Options.annotate = True
        
        setup_kwargs['cmdclass'] = dict(build_ext = build_ext)
        
        print('Cython will be used; cy_ext is ' + cy_ext)
    else:
        print('Cython will not be used; cy_ext is ' + cy_ext)

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
        
    # Set variables for C++ sources and include directories
    sources = find_cpp_sources(os.path.join(CProot,'src'), '*.cpp') 
    include_dirs = [os.path.join(CProot, 'include'), os.path.join(CProot, 'externals', 'Eigen')]

    ## If the file is run directly without any parameters, clean, build and install
    if len(sys.argv)==1:
       sys.argv += ['clean', 'build', 'install']
        
    common_args = dict(include_dirs = include_dirs,
                       language='c++')
   
    if USE_CYTHON:
        common_args.update(dict(cython_c_in_temp = True,
                                cython_directives = cython_directives
                                )
                           )
                        
    AbstractState_module = Extension('CoolProp5.AbstractState',
                        [os.path.join('CoolProp5','AbstractState.' + cy_ext)] + sources,
                        **common_args)
    CoolProp_module = Extension('CoolProp5.CoolProp',
                        [os.path.join('CoolProp5','CoolProp.' + cy_ext)] + sources,
                        **common_args)
    constants_module = Extension('CoolProp5.constants',
                        [os.path.join('CoolProp5','constants.' + cy_ext)],
                        **common_args)
     
     
    if not pypi:
        copy_files()

    ext_modules = [CoolProp_module, AbstractState_module, constants_module]
    
    if USE_CYTHON:
        ext_modules = cythonize(ext_modules)
        
    try:
        setup (name = 'CoolProp5',
               version = version, # look above for the definition of version variable - don't modify it here
               author = "Ian Bell",
               author_email='ian.h.bell@gmail.com',
               url='http://www.coolprop.org',
               description = """Open-source thermodynamic and transport properties database""",
               packages = find_packages(),
               ext_modules = ext_modules,
               package_data = {'CoolProp5':['State.pxd',
                                            'CoolProp.pxd',
                                            'CoolPropBibTeXLibrary.bib'] + find_cpp_sources(os.path.join('include'), '*.h')},
               classifiers = [
                "Programming Language :: Python",
                "Development Status :: 4 - Beta",
                "Environment :: Other Environment",
                "Intended Audience :: Developers",
                "License :: OSI Approved :: MIT License",
                "Operating System :: OS Independent",
                "Topic :: Software Development :: Libraries :: Python Modules"
                ],
               **setup_kwargs
               )
    except BaseException as E:
        if not pypi:
            remove_files()
        raise
    else:
        if not pypi:
            remove_files()
