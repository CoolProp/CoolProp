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

if '--windows' is in sys.argv:
    env = {}
else:
    env = dict(CC = 'i686-w64-mingw32-g++', CC = 'i686-w64-mingw32-gcc')
commons = dict(shell = True, stdout = sys.stdout, stderr = sys.stderr, env = env)
subprocess.check_call('Tools/pcre-build.sh', **commons)
subprocess.check_call('./autogen.sh', **commons)
subprocess.check_call('./configure --disable-ccache --with-scilab-inc=${SCILAB_HOME}/include --with-scilab=${SCILAB_HOME}/bin/scilab-cli --prefix=${PWD}/swig-scilab-bin', **commons)
subprocess.check_call('make', **commons)
subprocess.check_call('make install', **commons)
subprocess.check_call('cp swig swig3.0', cwd='swig-scilab-bin/bin', **commons)
subprocess.check_call('cp swig swig2.0', cwd='swig-scilab-bin/bin', **commons)