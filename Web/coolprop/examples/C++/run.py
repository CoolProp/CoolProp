import subprocess, sys
kwargs = dict(stdout = sys.stdout, stderr = sys.stderr, shell = True)
subprocess.check_call('cmake ../../../.. -DCOOLPROP_MY_MAIN=Example.cpp -DCMAKE_VERBOSE_MAKEFILE=ON', **kwargs)
subprocess.check_call('cmake --build .', **kwargs)