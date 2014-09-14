import subprocess, wget, os, shutil, sys, glob

if not os.path.exists('swig-scilab'):
	subprocess.call('git clone https://github.com/swig/swig.git -b gsoc2012-scilab swig-scilab', shell = True, stdout = sys.stdout, stderr = sys.stderr)
else:
	subprocess.call('git pull', shell = True, cwd = 'swig-scilab', stdout = sys.stdout, stderr = sys.stderr)
os.chdir('swig-scilab')
if not glob.glob('pcre-*.tar.gz'):
	wget.download('ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-8.33.tar.gz')
subprocess.check_call('Tools/pcre-build.sh', shell = True, stdout = sys.stdout, stderr = sys.stderr)
subprocess.check_call('./autogen.sh', shell = True, stdout = sys.stdout, stderr = sys.stderr)
subprocess.check_call('./configure --disable-ccache --with-scilab-inc=${SCILAB_HOME}/include --with-scilab=${SCILAB_HOME}/bin/scilab-cli --prefix=${PWD}/swig-scilab-bin', shell = True, stdout = sys.stdout, stderr = sys.stderr)
subprocess.check_call('make', shell = True, stdout = sys.stdout, stderr = sys.stderr)
subprocess.check_call('make install', shell = True, stdout = sys.stdout, stderr = sys.stderr)
