import subprocess, wget, os, shutil, sys, glob

if not os.path.exists('swig-matlab'):
	subprocess.call('git clone https://github.com/KrisThielemans/swig swig-matlab', shell = True, stdout = sys.stdout, stderr = sys.stderr)
else:
	subprocess.call('git pull', shell = True, cwd = 'swig-matlab', stdout = sys.stdout, stderr = sys.stderr)

os.chdir('swig-matlab')
if not glob.glob('pcre-*.tar.gz'):
    for rev in ['8.34','8.35','8.36']:
        try:
            wget.download('ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-'+rev+'.tar.gz'); break
        except:
            pass
        
subprocess.check_call('Tools/pcre-build.sh', shell = True, stdout = sys.stdout, stderr = sys.stderr)
subprocess.check_call('./autogen.sh', shell = True, stdout = sys.stdout, stderr = sys.stderr)
subprocess.check_call('./configure --disable-ccache --with-matlab=/usr/local/MATLAB/R2014a --prefix=${PWD}/swig-matlab-bin', shell = True, stdout = sys.stdout, stderr = sys.stderr)
subprocess.check_call('make', shell = True, stdout = sys.stdout, stderr = sys.stderr)
subprocess.check_call('make install', shell = True, stdout = sys.stdout, stderr = sys.stderr)
subprocess.check_call('cp swig swig3.0', shell = True, stdout = sys.stdout, stderr = sys.stderr, cwd='swig-matlab-bin/bin')
