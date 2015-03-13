import subprocess, os
from example_generator import *
            
def tee_call(call, file, **kwargs):
    callee = subprocess.Popen(call,
                           stdout = subprocess.PIPE,
                           stderr = subprocess.PIPE,
                           **kwargs)
    stdout, stderr = callee.communicate()
    print stdout, stderr
    file.write(stdout)
    file.write(stderr)
    if callee.poll() != 0:
        raise ValueError('Return code is non-zero')
            
if __name__=='__main__':

    #C++
    #kwargs = dict(stdout = sys.stdout, stderr = sys.stderr, shell = True)
    #subprocess.check_call('cmake ../../../.. -DCOOLPROP_MY_MAIN=Example.cpp -DCMAKE_VERBOSE_MAKEFILE=ON', **kwargs)
    #subprocess.check_call('cmake --build .', **kwargs)
    """
    if not os.path.exists('Octave'): os.mkdir('Octave')
    P = Python()
    code = P.parse()
    P.write('Python/Example.py', code)
    with open('Python/Example.out','w') as fp:
        tee_call(r'python Example.py', fp, shell = True, cwd = 'Python')

    if not os.path.exists('Octave'): os.mkdir('Octave')
    O = Octave()
    code = O.parse()
    O.write('Octave/Example.m', code)
    kwargs = dict(stdout = sys.stdout, stderr = sys.stderr, shell = True, cwd = 'Octave')
    subprocess.check_call('cmake ../../../.. -DCOOLPROP_OCTAVE_MODULE=ON -DCMAKE_VERBOSE_MAKEFILE=ON', **kwargs)
    subprocess.check_call('cmake --build .', **kwargs)
    with open('Octave/Example.out','w') as fp:
        tee_call(r'octave Example', fp, shell = True, cwd = 'Octave')
    
    if not os.path.exists('Java'): os.mkdir('Java')
    J = Java()
    code = J.parse()
    J.write('Java/Example.java', code)
    kwargs = dict(stdout = sys.stdout, stderr = sys.stderr, shell = True, cwd = 'Java')
    subprocess.check_call('cmake ../../../.. -DCOOLPROP_JAVA_MODULE=ON -DCMAKE_VERBOSE_MAKEFILE=ON', **kwargs)
    subprocess.check_call('cmake --build .', **kwargs)
    subprocess.check_call(r'javac *.java', **kwargs)
    with open('Java/Example.out','w') as fp:
        tee_call(r'java Example', fp, shell = True, cwd = 'Java')
    
    if not os.path.exists('Csharp'): os.mkdir('Csharp')
    C = Csharp()
    C.write('Csharp/Example.cs', C.parse())
    kwargs = dict(stdout = sys.stdout, stderr = sys.stderr, shell = True, cwd = 'Csharp')
    subprocess.check_call('cmake ../../../.. -DCOOLPROP_CSHARP_MODULE=ON -DCMAKE_VERBOSE_MAKEFILE=ON', **kwargs)
    subprocess.check_call('cmake --build .', **kwargs)
    subprocess.check_call(r'mcs -out:Example *.cs', **kwargs)
    with open('Csharp/Example.out','w') as fp:
        tee_call(r'mono Example', fp, shell = True, cwd = 'Csharp')
    """
    if not os.path.exists('MATLAB'): os.mkdir('MATLAB')
    #M = MATLAB()
    #M.write('MATLAB/Example.m', M.parse())
    kwargs = dict(stdout = sys.stdout, stderr = sys.stderr, shell = True, cwd = 'MATLAB')
    subprocess.check_call('PATH=${HOME}/swig-matlab/bin:$PATH cmake ../../../.. -DCOOLPROP_MATLAB_SWIG_MODULE=ON -DSWIG_DIR=${HOME}/swig-matlab/bin -DCMAKE_VERBOSE_MAKEFILE=ON', **kwargs)
    subprocess.check_call('PATH=${HOME}/swig-matlab/bin:$PATH cmake --build .', **kwargs)
    retcode = subprocess.call('matlab -nosplash -nojvm -nodesktop -nodisplay -r "result = runtests(\'Example\'); exit(result.Failed)" -logfile Example.out', shell = True, cwd = 'MATLAB')
    