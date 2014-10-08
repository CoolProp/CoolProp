#!/usr/bin/env python
# -*- coding: utf8 -*-
import os.path, glob, subprocess, sys, time
#
if len(sys.argv) < 2: 
    full_rebuild = False
if len(sys.argv)== 2: 
    if   sys.argv[1]=="True": full_rebuild = True
    elif sys.argv[1]=="1"   : full_rebuild = True
    else: full_rebuild = False
    print str(sys.argv[1])
    print "full_rebuild was set to: {0}".format(full_rebuild)
if len(sys.argv) > 2: 
    full_rebuild = False
    print "Cannot process more than one parameter: {0}".format(str(sys.argv))
#
def touch(fname):
    if os.path.exists(fname): os.utime(fname, None)
    else: open(fname, 'a').close()
#
def get_ftime(fname):
    if os.path.isfile(fname): return os.path.getctime(fname)
    else: return 0
#
web_dir      = os.path.abspath(os.path.join(os.path.dirname(__file__),'..'))
script_dir   = os.path.abspath(os.path.join(web_dir,'scripts'))
touch_file   = os.path.abspath(os.path.join(script_dir,'last_run'))
#
cur_time = time.time()
lim_days = 7
lim_time = cur_time - 60*60*24*lim_days # seconds
fil_time = get_ftime(touch_file)
#
if fil_time < lim_time and not full_rebuild:
    print "The static files have not been updated in {0} days, forcing an update now.".format(lim_days)
    full_rebuild = True
#
print "Executing the normal scripts for generating the static files."
script_files = glob.glob(os.path.join(script_dir,'*.py')) # Avoid recursion
script_files = [os.path.abspath(f) for f in script_files if not os.path.abspath(f)==os.path.abspath(__file__)]
for script in script_files:
    print "Executing {0}".format(script)
    subprocess.call('python {0}'.format(os.path.basename(script)), cwd=script_dir, shell=True)
#
if full_rebuild:
    print "Executing the computationally expensive scripts for generating the static files."
    touch(touch_file)
    script_files = glob.glob(os.path.join(script_dir,'*.sh'))
    for script in script_files:
	print "Executing {0}".format(script)
	subprocess.call('./{0}'.format(os.path.basename(script)), cwd=script_dir, shell=True)
else:
    print "Skipping the computationally expensive scripts for generating the static files."