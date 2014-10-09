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
fil_time = get_ftime(touch_file)
#
reg_hour   = time.strftime("%H")
reg_minute = time.strftime("%M")
sch_hour   = 13 #scheduled hour = 4am Boulder = 1pm CPH
sch_minute =  7 #scheduled minute = 7 past
#
lim_days = 1
lim_time = cur_time - 60*60*24*lim_days # seconds
if int(reg_hour)==sch_hour and sch_minute+2>reg_minute and sch_minute-2<reg_minute and not full_rebuild:
    print "This is a scheduled rebuild at {0}:{1}.".format(sch_hour,sch_minute)
    if fil_time < lim_time: full_rebuild = True
    else: print "It looks like the files have been rebuilt during the last day, reduced rebuild."
#
lim_days = 1.5
lim_time = cur_time - 60*60*24*lim_days # seconds
if fil_time < lim_time and not full_rebuild:
    print "The static files have not been updated in {0} days, forcing an update now.".format(lim_days)
    full_rebuild = True

#req_dir = [os.path.abspath(os.path.join(web_dir,'_static','fluid_properties','Incompressibles_reports'))]
#req_fil = [os.path.abspath(os.path.join(web_dir,'fluid_properties','Mixtures.csv')),
#  os.path.abspath(os.path.join(web_dir,'fluid_properties','PurePseudoPure.csv')),
#  os.path.abspath(os.path.join(web_dir,'fluid_properties','Incompressibles_pure-fluids.csv'))]
#
#for d in req_dir:
#    if not os.path.exists(d) and not full_rebuild:
#        print "The required directory {0} is missing, trying to rebuild it.".format(d)
#        full_rebuild = True
#for f in req_fil:
#    if not os.path.exists(f): 
#        print "The required file {0} is missing, trying to rebuild it.".format(f)
#        full_rebuild = True

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
	subprocess.call('chmod +x {0}'.format(os.path.basename(script)), cwd=script_dir, shell=True)
	subprocess.call('./{0}'.format(os.path.basename(script)), cwd=script_dir, shell=True)
else:
    print "Skipping the computationally expensive scripts for generating the static files."