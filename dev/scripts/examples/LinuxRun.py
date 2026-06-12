from __future__ import print_function
import subprocess, os, sys
from example_generator import *
import shutil
import codecs


def tee_call(call, file, **kwargs):
    callee = subprocess.Popen(call,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           **kwargs)
    stdout, stderr = callee.communicate()
    print(stdout, stderr)
    file.write(stdout.decode('utf-8'))
    file.write(stderr.decode('utf-8'))
    # if callee.poll() != 0:
    #    raise ValueError('Return code is non-zero')


def copyfiles(lang, ext):
    shutil.copy2(lang + '/Example.' + ext, '../../../Web/coolprop/wrappers/' + lang + '/Example.' + ext)
    shutil.copy2(lang + '/Example.out', '../../../Web/coolprop/wrappers/' + lang + '/Example.out')


if __name__ == '__main__':

    # Each COOLPROP_*_MODULE build below recompiles the entire CoolProp core
    # from scratch, so the core is compiled once per wrapper language (4x per
    # run).  A compiler cache avoids redoing that work when the core source is
    # unchanged.  This helps most ACROSS runs (docs-only changes, the daily/
    # weekly cron on an unchanged master, re-runs): a later build of the same
    # module is then a full cache hit.  It does NOT share core objects between
    # the four modules within a single run -- the Octave/Java/R cmake blocks add
    # their header dirs via global include_directories(), so those -I flags land
    # on the core .cpp compile lines and make each module's core objects
    # distinct.  Scoping those includes to the swig target would unlock within-
    # run sharing (tracked separately).  Empty string => no-op when ccache isn't
    # installed, so local runs without it behave exactly as before.
    CCACHE = (' -DCMAKE_CXX_COMPILER_LAUNCHER=ccache -DCMAKE_C_COMPILER_LAUNCHER=ccache'
              if shutil.which('ccache') else '')
    if CCACHE:
        # Zero the stats so the ccache -s at the end reads as a per-run hit rate.
        subprocess.call('ccache -z', shell=True)

    # C++
    #kwargs = dict(stdout = sys.stdout, stderr = sys.stderr, shell = True)
    #subprocess.check_call('cmake ../../../.. -DCOOLPROP_MY_MAIN=Example.cpp -DCMAKE_VERBOSE_MAKEFILE=ON', **kwargs)
    #subprocess.check_call('cmake --build .', **kwargs)

    if not os.path.exists('Python'): os.mkdir('Python')
    P = Python()
    code = P.parse()
    P.write('Python/Example.py', code)
    with codecs.open('Python/Example.out', 'w', encoding='utf-8') as fp:
        tee_call('python Example.py', fp, shell=True, cwd='Python')
    copyfiles('Python', 'py')

    if not os.path.exists('Octave'): os.mkdir('Octave')
    O = Octave()
    O.write('Octave/Example.m', O.parse())
    kwargs = dict(stdout=sys.stdout, stderr=sys.stderr, shell=True, cwd='Octave')
    subprocess.check_call('cmake ../../../.. -G Ninja -DCOOLPROP_OCTAVE_MODULE=ON -DCMAKE_VERBOSE_MAKEFILE=ON -DCOOLPROP_NO_EXAMPLES=ON' + CCACHE, **kwargs)
    subprocess.check_call('cmake --build .', **kwargs)
    with codecs.open('Octave/Example.out', 'w', encoding='utf-8') as fp:
       tee_call(r'octave Example.m', fp, shell=True, cwd='Octave')
    copyfiles('Octave', 'm')
    #
    if not os.path.exists('Java'): os.mkdir('Java')
    J = Java()
    J.write('Java/Example.java', J.parse())
    kwargs = dict(stdout=sys.stdout, stderr=sys.stderr, shell=True, cwd='Java')
    subprocess.check_call('cmake ../../../.. -G Ninja -DCOOLPROP_JAVA_MODULE=ON -DCMAKE_VERBOSE_MAKEFILE=ON -DCOOLPROP_NO_EXAMPLES=ON' + CCACHE, **kwargs)
    subprocess.check_call('cmake --build .', **kwargs)
    subprocess.check_call(r'javac *.java', **kwargs)
    with codecs.open('Java/Example.out', 'w', encoding='utf-8') as fp:
        tee_call(r'java -Djava.library.path="'+os.path.abspath('Java')+'" Example', fp, shell=True, cwd='Java')
    copyfiles('Java', 'java')
    
    if not os.path.exists('Csharp'): os.mkdir('Csharp')
    C = Csharp()
    C.write('Csharp/Example.cs', C.parse())
    kwargs = dict(stdout=sys.stdout, stderr=sys.stderr, shell=True, cwd='Csharp')
    subprocess.check_call('cmake ../../../.. -G Ninja -DCOOLPROP_CSHARP_MODULE=ON -DCMAKE_VERBOSE_MAKEFILE=ON -DCOOLPROP_NO_EXAMPLES=ON' + CCACHE, **kwargs)
    subprocess.check_call('cmake --build .', **kwargs)
    subprocess.check_call(r'mcs -out:Example *.cs', **kwargs)
    with codecs.open('Csharp/Example.out', 'w', encoding='utf-8') as fp:
        tee_call(r'mono Example', fp, shell=True, cwd='Csharp')
    copyfiles('Csharp', 'cs')
    
    if not os.path.exists('R'): os.mkdir('R')
    RR = R()
    RR.write('R/Example.R', RR.parse())
    kwargs = dict(stdout=sys.stdout, stderr=sys.stderr, shell=True, cwd='R')
    subprocess.check_call('cmake ../../../.. -G Ninja -DCOOLPROP_R_MODULE=ON -DCMAKE_VERBOSE_MAKEFILE=ON -DCOOLPROP_NO_EXAMPLES=ON -DR_BIN=/usr/bin' + CCACHE, **kwargs)
    subprocess.check_call('cmake --build .', **kwargs)
    with codecs.open('R/Example.out', 'w', encoding='utf-8') as fp:
        tee_call(r'DYLD_LIBRARY_PATH=/opt/refprop Rscript Example.R', fp, shell=True, cwd='R')
    copyfiles('R', 'R')

    if CCACHE:
        # Surface the hit rate in the CI log so a regression (e.g. a flag that
        # silently busts the cache across modules) is visible.
        subprocess.call('ccache -s', shell=True)
    #
    #if not os.path.exists('MATLAB'): os.mkdir('MATLAB')
    #M = MATLAB()
    #M.write('MATLAB/Example.m', M.parse())
    #kwargs = dict(stdout = sys.stdout, stderr = sys.stderr, shell = True, cwd = 'MATLAB')
    #subprocess.check_call('PATH=${HOME}/swig-matlab/bin:$PATH cmake ../../../.. -DCOOLPROP_MATLAB_SWIG_MODULE=ON -DSWIG_DIR=${HOME}/swig-matlab/bin -DCMAKE_VERBOSE_MAKEFILE=ON', **kwargs)
    #subprocess.check_call('PATH=${HOME}/swig-matlab/bin:$PATH cmake --build .', **kwargs)
    #retcode = subprocess.call('matlab -nosplash -nojvm -nodesktop -nodisplay -r "result = runtests(\'Example\'); exit(result.Failed)" -logfile Example.out', shell = True, cwd = 'MATLAB')
    #copyfiles('MATLAB','m')
