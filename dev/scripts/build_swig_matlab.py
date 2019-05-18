import subprocess, wget, os, shutil, sys, glob

if '--windows' in sys.argv:
    compilers = " CXX=i686-w64-mingw32-g++ CC=i686-w64-mingw32-gcc "
    host = " --host=x86_64-unknown-linux-gnu "
    extra = ' LDFLAGS="-static-libgcc -static-libstdc++ -static"'
else:
    compilers = ''
    extra = ''
    host = ''
commons = dict(shell=True, stdout=sys.stdout, stderr=sys.stderr)

if not os.path.exists('swig-matlab'):
    subprocess.call('git clone https://github.com/jaeandersson/swig swig-matlab', **commons)
else:
    subprocess.call('git pull', shell=True, cwd='swig-matlab', stdout=sys.stdout, stderr=sys.stderr)

os.chdir('swig-matlab')
if '--no-pcre' not in sys.argv:
    if not glob.glob('pcre-*.tar.gz'):
        for rev in ['8.34', '8.35', '8.36']:
            try:
                wget.download('ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-' + rev + '.tar.gz'); break
            except:
                pass
    subprocess.check_call('Tools/pcre-build.sh' + compilers + host, **commons)

prefix = '${PWD}/swig-matlab-bin'
if '--prefix' in sys.argv:
    prefix = sys.argv[sys.argv.index('--prefix') + 1]

subprocess.check_call(compilers + './autogen.sh', **commons)
subprocess.check_call('./configure --disable-ccache --with-matlab=/usr/local/MATLAB/R2014a --prefix=' + prefix + extra + compilers + host, **commons)
subprocess.check_call(compilers + 'make', **commons)
subprocess.check_call(compilers + 'make install', **commons)

if '--windows' in sys.argv:
    subprocess.check_call('cp swig.exe swig3.0.exe', cwd='swig-matlab-bin/bin', **commons)
    subprocess.check_call('cp swig.exe swig2.0.exe', cwd='swig-matlab-bin/bin', **commons)
else:
    subprocess.check_call('cp swig swig3.0', cwd='swig-matlab-bin/bin', **commons)
    subprocess.check_call('cp swig swig2.0', cwd='swig-matlab-bin/bin', **commons)
