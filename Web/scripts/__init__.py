#!/usr/bin/env python
# -*- coding: utf8 -*-
import os.path, glob, subprocess, sys, time, datetime, pytz
#
if len(sys.argv) < 2:
    full_rebuild = False
if len(sys.argv) == 2:
    if sys.argv[1] == "True": full_rebuild = True
    elif sys.argv[1] == "1": full_rebuild = True
    else: full_rebuild = False
if len(sys.argv) > 2:
    full_rebuild = False
    print("Cannot process more than one parameter: {0}".format(str(sys.argv)))
#


def touch(fname):
    if os.path.exists(fname): os.utime(fname, None)
    else: open(fname, 'a').close()
#


def get_ftime(fname):
    if os.path.isfile(fname): return os.path.getctime(fname)
    else: return 0


#
web_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
script_dir = os.path.abspath(os.path.join(web_dir, 'scripts'))
touch_file = os.path.abspath(os.path.join(script_dir, 'last_run'))
root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
#
cur_time = time.time()
fil_time = get_ftime(touch_file)
#
# Static execution time
#reg_hour   = time.strftime("%H")
#reg_minute = time.strftime("%M")
# sch_hour   = 12 #scheduled hour = 3am Boulder = 12pm CPH
# sch_minute =  7 #scheduled minute = 7 past
#
# Dynamically calculated execution (includes daylight saving time etc
masterTime = pytz.timezone('US/Pacific')
#slaveTime  = pytz.timezone('Europe/Copenhagen')
now_master = datetime.datetime.now(masterTime)
run_master = datetime.datetime.strptime("03:00:00", '%H:%M:%S')
#
now_master = datetime.time(now_master.hour, now_master.minute, now_master.second)
run_master = datetime.time(run_master.hour, run_master.minute, run_master.second)
run_master_end = datetime.time(run_master.hour, run_master.minute + 5, run_master.second)
#
lim_days = 0.90
lim_time = cur_time - 60 * 60 * 24 * lim_days  # seconds
#
if now_master >= run_master and \
   now_master <= run_master_end and \
   not full_rebuild:
    print("This is a scheduled rebuild at {0}.".format(run_master))
    if fil_time < lim_time: full_rebuild = True
    else: print("It looks like the files have been rebuilt during the last day.")
#
lim_days = 3
lim_time = cur_time - 60 * 60 * 24 * lim_days  # seconds
if fil_time < lim_time and not full_rebuild:
    print("The static files have not been updated in {0} days, forcing an update now.".format(lim_days))
    full_rebuild = True

#req_dir = [os.path.abspath(os.path.join(web_dir,'_static','fluid_properties','Incompressibles_reports'))]
# req_fil = [os.path.abspath(os.path.join(web_dir,'fluid_properties','Mixtures.csv')),
#  os.path.abspath(os.path.join(web_dir,'fluid_properties','PurePseudoPure.csv')),
#  os.path.abspath(os.path.join(web_dir,'fluid_properties','Incompressibles_pure-fluids.csv'))]
#
# for d in req_dir:
#    if not os.path.exists(d) and not full_rebuild:
#        print "The required directory {0} is missing, trying to rebuild it.".format(d)
#        full_rebuild = True
# for f in req_fil:
#    if not os.path.exists(f):
#        print "The required file {0} is missing, trying to rebuild it.".format(f)
#        full_rebuild = True
# print "Executing the normal scripts for generating the static files."
# script_files = glob.glob(os.path.join(script_dir,'*.py')) # Avoid recursion
# script_files = [os.path.abspath(f) for f in script_files if not os.path.abspath(f)==os.path.abspath(__file__)]
# for script in script_files:
#     print "Executing {0}".format(script)
#     subprocess.call('python {0}'.format(os.path.basename(script)), cwd=script_dir, shell=True)
#


def run_script(path):
    if os.path.exists(path):
        file_path = os.path.dirname(path)
        file_name = os.path.basename(path)
        file_extension = path.split(".")[-1]
        #file_name, file_extension = os.path.splitext(path)

        if file_extension.lower() == "py":
            subprocess.check_call('python -u {0}'.format(file_name), cwd=file_path, shell=True)
        elif file_extension.lower() == "sh" or file_extension.lower() == "bsh":
            subprocess.check_call('chmod +x {0}'.format(file_name), cwd=file_path, shell=True)
            subprocess.check_call('./{0}'.format(file_name), cwd=file_path, shell=True)
        else:
            print("Unknown file extension in {0}".format(path))
    else:
        print("Could not find the file {0}".format(path))


# Inject the version of CoolProp into the doxygen configuration files
# Put it at the end, overwrites prior value
import CoolProp
with open(os.path.join(root_dir, 'Doxyfile'), 'a+') as fp:
    fp.write('\n\n PROJECT_NUMBER         = ' + CoolProp.__version__ + '\n')

# The normal tasks that are carried out each time the script runs
normal_tasks = ["../../dev/scripts/examples/LinuxRun.py", "coolprop.tabular.speed.py", "fluid_properties.phase_envelope.py", "fluid_properties.PurePseudoPure.py", "fluid_properties.Mixtures.py", "coolprop.parametric_table.py", "coolprop.configuration.py", "logo_2014.py", "fluid_properties.REFPROPcomparison.py"]
# The expensive tasks that are fired when full_rebuild is True
expensive_tasks = ["fluid_properties.Consistency.py", "fluid_properties.Incompressibles.sh"]
print("Executing the normal scripts for generating static files.")
for script in normal_tasks:
    print("Executing {0}".format(script))
    run_script(os.path.normpath(os.path.join(script_dir, script)))
#
if full_rebuild:
    print("Executing the computationally expensive scripts for generating the static files.")
    for script in expensive_tasks:
        print("Executing {0}".format(script))
        run_script(os.path.join(script_dir, script))
    touch(touch_file)
else:
    print("Skipping the computationally expensive scripts for generating the static files.")
