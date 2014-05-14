from __future__ import print_function
import os, subprocess, sys
import subprocess
import itertools
import platform

def validate_pair(ob):
    try:
        if not (len(ob) == 2):
            print("Unexpected result:", ob, file=sys.stderr)
            raise ValueError
    except:
        return False
    return True

def consume(iter):
    try:
        while True: next(iter)
    except StopIteration:
        pass    
    
# See http://stackoverflow.com/questions/1214496/how-to-get-environment-from-a-subprocess-in-python
def get_environment_from_batch_command(env_cmd, opts = None, initial=None):
    """
    Take a command (either a single command or list of arguments)
    and return the environment created after running that command.
    Note that if the command must be a batch file or .cmd file, or the
    changes to the environment will not be captured.

    If initial is supplied, it is used as the initial environment passed
    to the child process.
    """
    if not isinstance(env_cmd, (list, tuple)):
        env_cmd = [env_cmd]
    # construct the command that will alter the environment
    env_cmd = subprocess.list2cmdline(env_cmd)
    # create a tag so we can tell in the output when the proc is done
    tag = 'Done running command'
    # construct a cmd.exe command to do accomplish this
    cmd = 'cmd.exe /s /c "{env_cmd} {opts} && echo "{tag}" && set"'.format(**vars())
    # launch the process
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell = True, env=initial)
    if proc.returncode == 0:
        return None
    # parse the output sent to stdout
    lines = proc.stdout
    # consume whatever output occurs until the tag is reached
    consume(itertools.takewhile(lambda l: tag not in l, lines))
    # define a way to handle each KEY=VALUE line
    handle_line = lambda l: l.rstrip().split('=',1)
    # parse key/values into pairs
    pairs = map(handle_line, lines)
    # make sure the pairs are valid
    valid_pairs = filter(validate_pair, pairs)
    # Upper case for the pair
    for i in range(len(valid_pairs)):
        valid_pairs[i][0] = valid_pairs[i][0].upper()
    # construct a dictionary of the pairs
    result = dict(valid_pairs)
    # let the process finish
    proc.communicate()
    
    return result
    
def find_VS_installations():
    # Only on windows
    if platform.system() != 'Windows': return None
        
    installations = []
    
    # Use environmental variables to file candidate program files folders on this system
    program_file_paths = set([os.environ[k] for k in ['PROGRAMFILES','PROGRAMFILES(X86)','ProgramW6432']])
    for ver in ['8.0','9.0','10.0','11.0','12.0','13.0','14.0','15.0']:
        for PF_path in program_file_paths:
            for bitness in ['x86','amd64']:
                    
                bat_path = os.path.join(PF_path,'Microsoft Visual Studio '+ver,'VC','vcvarsall.bat')
                bin_path = os.path.join(PF_path,'Microsoft Visual Studio '+ver,'VC','bin')

                # Skip if this compiler bat file doesn't exist
                if not os.path.exists(bat_path): 
                    continue
                
                env = get_environment_from_batch_command(bat_path, opts = bitness)
                
                if 'LIB' not in env: continue
                
                if bitness =='x86':
                    bits = 32
                elif bitness == 'amd64':
                    bits = 64
                
                installations.append(dict(ver = ver, bits = bits, bin_path = bin_path, env = env))
    return installations
                    
def find_cpp_sources(root = os.path.join('..','..','src'), extensions = ['.cpp'], skip_files = None):
    file_listing = []
    for path, dirs, files in os.walk(root):
        for file in files:
            n,ext = os.path.splitext(file)
            fname = os.path.relpath(os.path.join(path, file))
            if skip_files is not None and fname in skip_files:
                continue
                
            if ext in extensions:
                file_listing.append(fname)
    return file_listing
    
if __name__=='__main__':
    print("don't run me")