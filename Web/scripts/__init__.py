#!/usr/bin/env python
# -*- coding: utf8 -*-
import os.path, glob, subprocess, sys, time, datetime, pytz

# Start with detecting the full build
def detect_full_rebuild(args):
    if len(sys.argv) >= 2:
        arg = str(sys.argv[1]).lower()
        if arg == "true": return True
        if arg == "1": return True
    return False


full_rebuild = detect_full_rebuild(sys.argv)
print("Detected rebuild argument: full_rebuild = {}".format(full_rebuild))

# File system functions
def touch(fname):
    if os.path.exists(fname): os.utime(fname, None)
    else: open(fname, 'a').close()
#
def get_ftime(fname):
    if os.path.isfile(fname): return os.path.getctime(fname)
    else: return 0
#
# Directory settings
script_root_dir = os.path.abspath(os.path.dirname(__file__))
repo_root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
task_list = []

def add_if_exists(fname):
    if os.path.isfile(fname):
        task_list.append(fname)
        print("Added '{}' to the task list.".format(fname))
        return True
    return False

def add_to_task_list(fname):
    if add_if_exists(os.path.abspath(fname))
        return True
    fname = os.path.join(script_root_dir, fname)
    if add_if_exists(os.path.abspath(fname))
        return True
    fname = os.path.join(repo_root_dir, fname)
    if add_if_exists(os.path.abspath(fname))
        return True
    return False

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
normal_tasks = ["dev/scripts/examples/LinuxRun.py", "coolprop.tabular.speed.py", "fluid_properties.phase_envelope.py", "fluid_properties.PurePseudoPure.py", "fluid_properties.Mixtures.py", "coolprop.parametric_table.py", "coolprop.configuration.py", "logo_2014.py", "fluid_properties.REFPROPcomparison.py"]
# The expensive tasks that are fired when full_rebuild is True
expensive_tasks = ["fluid_properties.Consistency.py", "fluid_properties.Incompressibles.sh"]

print("Adding the normal scripts to the task list.")
selected_tasks = normal_tasks[:]
if full_rebuild:
    print("Adding the computationally expensive scripts to the task list.")
    selected_tasks += expensive_tasks[:]

for fname in selected_tasks:
    add_to_task_list(fname)

print("Processing the selected tasks to generate the static files.")
for fname in task_list:
    print("Executing {0}".format(fname))
    run_script(os.path.normpath(fname))
