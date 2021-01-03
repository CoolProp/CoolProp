from __future__ import print_function
import subprocess, os
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
    subprocess.check_call('cmake ../../../.. -DCOOLPROP_OCTAVE_MODULE=ON -DCMAKE_VERBOSE_MAKEFILE=ON -DCOOLPROP_NO_EXAMPLES=ON', **kwargs)
    subprocess.check_call('cmake --build .', **kwargs)
    with codecs.open('Octave/Example.out', 'w', encoding='utf-8') as fp:
        tee_call(r'octave Example.m', fp, shell=True, cwd='Octave')
    copyfiles('Octave', 'm')

    if not os.path.exists('Java'): os.mkdir('Java')
    J = Java()
    J.write('Java/Example.java', J.parse())
    kwargs = dict(stdout=sys.stdout, stderr=sys.stderr, shell=True, cwd='Java')
    subprocess.check_call('cmake ../../../.. -DCOOLPROP_JAVA_MODULE=ON -DCMAKE_VERBOSE_MAKEFILE=ON -DCOOLPROP_NO_EXAMPLES=ON', **kwargs)
    subprocess.check_call('cmake --build .', **kwargs)
    subprocess.check_call(r'javac *.java', **kwargs)
    with codecs.open('Java/Example.out', 'w', encoding='utf-8') as fp:
        tee_call(r'java -Djava.library.path="'+os.path.abspath('Java')+'" Example', fp, shell=True, cwd='Java')
    copyfiles('Java', 'java')

    if not os.path.exists('Csharp'): os.mkdir('Csharp')
    C = Csharp()
    C.write('Csharp/Example.cs', C.parse())
    kwargs = dict(stdout=sys.stdout, stderr=sys.stderr, shell=True, cwd='Csharp')
    subprocess.check_call('cmake ../../../.. -DCOOLPROP_CSHARP_MODULE=ON -DCMAKE_VERBOSE_MAKEFILE=ON -DCOOLPROP_NO_EXAMPLES=ON', **kwargs)
    subprocess.check_call('cmake --build .', **kwargs)
    subprocess.check_call(r'mcs -out:Example *.cs', **kwargs)
    with codecs.open('Csharp/Example.out', 'w', encoding='utf-8') as fp:
        tee_call(r'mono Example', fp, shell=True, cwd='Csharp')
    copyfiles('Csharp', 'cs')

    if not os.path.exists('R'): os.mkdir('R')
    RR = R()
    RR.write('R/Example.R', RR.parse())
    kwargs = dict(stdout=sys.stdout, stderr=sys.stderr, shell=True, cwd='R')
    subprocess.check_call('cmake ../../../.. -DCOOLPROP_R_MODULE=ON -DCMAKE_VERBOSE_MAKEFILE=ON -DCOOLPROP_NO_EXAMPLES=ON -DR_BIN=/usr/bin', **kwargs)
    subprocess.check_call('cmake --build .', **kwargs)
    with codecs.open('R/Example.out', 'w', encoding='utf-8') as fp:
        tee_call(r'DYLD_LIBRARY_PATH=/opt/refprop Rscript Example.R', fp, shell=True, cwd='R')
    copyfiles('R', 'R')

    #if not os.path.exists('MATLAB'): os.mkdir('MATLAB')
    #M = MATLAB()
    #M.write('MATLAB/Example.m', M.parse())
    #kwargs = dict(stdout = sys.stdout, stderr = sys.stderr, shell = True, cwd = 'MATLAB')
    #subprocess.check_call('PATH=${HOME}/swig-matlab/bin:$PATH cmake ../../../.. -DCOOLPROP_MATLAB_SWIG_MODULE=ON -DSWIG_DIR=${HOME}/swig-matlab/bin -DCMAKE_VERBOSE_MAKEFILE=ON', **kwargs)
    #subprocess.check_call('PATH=${HOME}/swig-matlab/bin:$PATH cmake --build .', **kwargs)
    #retcode = subprocess.call('matlab -nosplash -nojvm -nodesktop -nodisplay -r "result = runtests(\'Example\'); exit(result.Failed)" -logfile Example.out', shell = True, cwd = 'MATLAB')
    # copyfiles('MATLAB','m')
