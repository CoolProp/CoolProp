from __future__ import print_function
import os
import json
import sys
import subprocess
import util 
import platform

def stdout_check_call(*args, **kwargs):
    print('calling: '+args[0])
    subprocess.check_call(*args, shell = True, stdout = sys.stdout, stderr = sys.stdout, **kwargs)

def find_VS_compiler(compilers, bits, ver):
    for c in compilers:
        if int(c['bits']) == bits and ver == c['ver']:
            return c
    raise ValueError('no VS compiler found')
    
def newest_VS_compiler(compilers, bits):
    joined = [(float(c['ver']),c) for c in compilers if int(c['bits']) == bits]
    joined = sorted(joined)
    joined.reverse()
    return joined[0][1]
        
def build_target(target, project_root, cpp_sources, project_name = 'Java'):
    """
    
    """
    if 'name_suffix' in target:
        suffix = '-' + target['name_suffix']
    else:
        suffix = ''
    build_name = platform.system() + '-' + target['compiler'] + '-' + str(target['bitness'])+'bit' + suffix
    build_path = os.path.join(project_name, build_name)
    try: 
        os.makedirs(build_path) 
    except: 
        pass
    dist_path = os.path.join(project_name, build_name)
    
    # Add any other sources needed
    if 'sources' in target:
        cpp_sources += target['sources']
        
    # Replace macros in sources
    for i,source in enumerate(cpp_sources):
        cpp_sources[i] = source.replace("%PACKAGE_HOME%", project_root)
    
    # Check that the platform agrees with this platform
    if platform.system() not in target['platform']:
        print('skipped this run since this platform is invalid:'+str(target))
    
    # Do any sort of preliminary processing needed based on the pre key
    if 'pre' in target:
        if not isinstance(target['pre'], basestring): 
            raise ValueError("target['pre'] is not a string: "+str(target['pre']))
        if target['pre'] != '':
            stdout_check_call(target['pre'])
    
    # Get compiler name
    compiler = target['compiler']
    if not isinstance(compiler, basestring): 
        raise ValueError("compiler is not a string: "+str(compiler))
        
    # Do the compilation step
    VS_keys = ['VS+','VS8.0','VS9.0','VS10.0','VS11.0','VS12.0','VS13.0','VS14.0','VS15.0']
    
    if compiler in VS_keys:
        vs_installs = util.find_VS_installations()
        if compiler == 'VS+':
            selected_compiler = newest_VS_compiler(vs_installs, int(target['bitness']))
        else:
            selected_compiler = find_VS_compiler(vs_installs, int(target['bitness']), compiler.lstrip('VS'))
            
        ### *************** COMPILE ********************
        call_dict = dict(c_flags = target['c_flags'], 
                         cpp_sources = ' '.join(cpp_sources),
                         bin_path = '"' + selected_compiler['bin_path'] + '"'
                         )
        
        for file in cpp_sources:
            call_dict['file'] = file
            
            path,ofile = os.path.split(file)
            
            # Get the output file name in the build folder
            call_dict['ofile'] = os.path.join(build_path, ofile.rsplit('.',1)[0]+'.obj')
            # Build the string to be run
            compile_string = 'cl /nologo /c {c_flags} /Fo{ofile} {file}'.format(**call_dict)
            # Actually run the build command using the environment for the selected compiler
            stdout_check_call(compile_string, env = selected_compiler['env'])
        
        if 'link_type' in target:
            if target['link_type'] == 'DLL':
                ### *************************** LINK ***************************
                call_dict = dict(link_flags = target['link_flags'], 
                                 bin_path = '"' + selected_compiler['bin_path'] + '"',
                                 build_path = build_path,
                                 link_fname = target['link_fname']
                                 )
                link_string = 'link /DLL /nologo {build_path}/*.obj {link_flags} /OUT:{build_path}/{link_fname}'.format(**call_dict)
                # Actually run the link command using the environment for the selected compiler
                stdout_check_call(link_string, env = selected_compiler['env'])
        
    elif compiler == 'GCC':
        pass
    elif compiler == 'clang':
        pass
    else:
        raise ValueError()
    
    # Do the post step
    if 'post' in target:
        if not isinstance(target['post'], basestring): 
            raise ValueError("target['post'] is not a string: "+str(target['post']))
        if target['post'] != '':
            stdout_check_call(target['post'])

def build_all_targets(project_root, project_name):
    json_file_name = os.path.relpath(os.path.join(project_root,'build.json'))
    
    if not os.path.exists(json_file_name): 
        raise OSError('File [{json_file_name}] does not exist'.format(**vars()))
        
    # Load the JSON file
    project = json.load(open(json_file_name,'r'))
    
    # Find all the C++ sources that will need to be compiled for sure
    cpp_sources = util.find_cpp_sources(os.path.join('..','..','src'))
    
    for target in project:
        build_target(target, project_root = project_root, cpp_sources = cpp_sources, project_name = project_name)

if __name__=='__main__':
#     build_all_targets(project_root = '../../wrappers/EES', project_name = 'EES')
#     build_all_targets(project_root = '../../wrappers/C#', project_name = 'Csharp')
#     build_all_targets(project_root = '../../wrappers/Java', project_name = 'Java')
    build_all_targets(project_root = '../../wrappers/SharedLibrary', project_name = 'SharedLibrary')