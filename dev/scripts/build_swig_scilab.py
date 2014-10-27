import subprocess, wget, os, shutil, sys, glob

if not os.path.exists('swig-scilab'):
	subprocess.call('git clone https://github.com/swig/swig.git -b gsoc2012-scilab swig-scilab', shell = True, stdout = sys.stdout, stderr = sys.stderr)
else:
	subprocess.call('git pull', shell = True, cwd = 'swig-scilab', stdout = sys.stdout, stderr = sys.stderr)
os.chdir('swig-scilab')
if not glob.glob('pcre-*.tar.gz'):
    for rev in ['8.34','8.35','8.36']:
        try:
            wget.download('ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-'+rev+'.tar.gz'); break
        except:
            pass

prefix = '--prefix=${PWD}/swig-scilab-bin'
if '--windows' in sys.argv:
    compilers = "CXX=i686-w64-mingw32-g++ CC=i686-w64-mingw32-gcc "
    prefix += ' LDFLAGS="-static-libgcc -static-libstdc++ -static"'
else:
    compilers = ''
commons = dict(shell = True, stdout = sys.stdout, stderr = sys.stderr)

subprocess.check_call(compilers+'Tools/pcre-build.sh', **commons)
subprocess.check_call(compilers+'./autogen.sh', **commons)
subprocess.check_call(compilers+'./configure --disable-ccache --with-scilab-inc=${SCILAB_HOME}/include --with-scilab=${SCILAB_HOME}/bin/scilab-cli'+' '+prefix, **commons)
subprocess.check_call(compilers+'make', **commons)
subprocess.check_call(compilers+'make install', **commons)

if '--windows' in sys.argv:
    subprocess.check_call('cp swig.exe swig3.0.exe', cwd='swig-scilab-bin/bin', **commons)
    subprocess.check_call('cp swig.exe swig2.0.exe', cwd='swig-scilab-bin/bin', **commons)
else:
    subprocess.check_call('cp swig swig3.0', cwd='swig-scilab-bin/bin', **commons)
    subprocess.check_call('cp swig swig2.0', cwd='swig-scilab-bin/bin', **commons)
