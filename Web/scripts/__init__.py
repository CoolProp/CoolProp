#!/usr/bin/env python
# -*- coding: utf8 -*-
import os.path, glob, subprocess, sys, time, datetime, pytz

# Start with detecting the full build
def detect_full_rebuild():
    if len(sys.argv) >= 2:
        arg = str(sys.argv[1]).lower()
        if arg == "true": return True
        if arg == "1": return True
    return False


full_rebuild = detect_full_rebuild()
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

def add_to_task_list(fname_in):
    fname = fname_in
    if add_if_exists(os.path.abspath(fname)):
        return True
    fname = os.path.join(script_root_dir, fname_in)
    if add_if_exists(os.path.abspath(fname)):
        return True
    fname = os.path.join(repo_root_dir, fname_in)
    if add_if_exists(os.path.abspath(fname)):
        return True
    print("Error: Could not find '{}'.".format(fname_in))
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
            subprocess.check_call('sed -i "s/\r//g" {0}'.format(file_name), cwd=file_path, shell=True)
            subprocess.check_call('./{0}'.format(file_name), cwd=file_path, shell=True)
        else:
            print("Unknown file extension in {0}".format(path))
    else:
        print("Could not find the file {0}".format(path))


# Inject the version of CoolProp into the doxygen configuration files
# Put it at the end, overwrites prior value
import CoolProp
with open(os.path.join(repo_root_dir, 'Doxyfile'), 'a+') as fp:
    fp.write('\n\n PROJECT_NUMBER         = {}\n'.format(CoolProp.__version__))

# The normal tasks that are carried out each time the script runs
print("Adding the normal scripts to the task list.")
add_to_task_list("dev/scripts/examples/LinuxRun.py")
add_to_task_list("coolprop.tabular.speed.py")
add_to_task_list("fluid_properties.phase_envelope.py")
add_to_task_list("fluid_properties.PurePseudoPure.py")
add_to_task_list("fluid_properties.Mixtures.py")
add_to_task_list("coolprop.parametric_table.py")
add_to_task_list("coolprop.configuration.py")
add_to_task_list("logo_2014.py")
add_to_task_list("fluid_properties.REFPROPcomparison.py")

# The expensive tasks that are fired when full_rebuild is True
if full_rebuild:
    print("Adding the computationally expensive scripts to the task list.")
    add_to_task_list("fluid_properties.Consistency.py")
    add_to_task_list("fluid_properties.Incompressibles.sh")

# Run all the files in the task list
print("Processing the selected tasks to generate the static files.")
for fname in task_list:
    print("Executing {0}".format(fname))
    run_script(os.path.normpath(fname))
