# Hack due to hang in plot_directive for some unknown reason
import subprocess, sys
subprocess.call('python methane-ethane.py', shell=True, stdout=sys.stdout, stderr=sys.stderr, cwd='../fluid_properties')
