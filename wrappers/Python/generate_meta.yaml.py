import jinja2
from jinja2 import Environment
import os,sys
import requests
import json
from distutils.version import LooseVersion #, StrictVersion
import codecs

""" A simple script to create a conda recipe from the PyPI release and build the package"""

template = """
package:
  name: coolprop
  version: {{ version }}

{% if pypi %}
source:
  fn: {{ fil }}
  url: {{ url }}
  md5: {{ md5 }}
{% endif %}

{% if local %}
source:
  path: .
{% endif %}
 

  # If this is a new build for the same version, increment the build
  # number. If you do not include this key, it defaults to 0.
  # number: 1

requirements:
  build:
    - python
    - setuptools{% for pkg in pkgs %}
    - {{ pkg -}}
{% endfor %}

  run:
    - python{% for pkg in pkgs %}
    - {{ pkg -}}
{% endfor %}

test:
  # Python imports
  imports:
    - CoolProp
    #- CoolProp.GUI
    #- CoolProp.Plots
    - CoolProp.tests

  # commands:
    # You can put test commands to be run here.  Use this to test that the
    # entry points work.


  # You can also put a file called run_test.py in the recipe that will be run
  # at test time.

  # requires:
    # Put any additional test requirements here.  For example
    # - nose

about:
  home: {{ home }}
  license: {{ license }}
  summary: {{ summary }}

"""

target_dir = os.path.join(os.path.dirname(__file__),'..','..')

#loader = jinja2.FileSystemLoader(template_dir)
#environment = jinja2.Environment(loader=loader)
#template = environment.get_template(os.path.join(template_dir,'conda_'+target+'.tpl'))
#
template =Environment().from_string(template)


pypi = False 
local = not pypi

pkgs = ["numpy", "scipy", "matplotlib", "pandas", "cython"]
target = 'meta.yaml'

if pypi:
    # Get the additional information from PyPI 
    r = requests.get('https://pypi.python.org/pypi/CoolProp/json')
    if(r.ok):
        item = json.loads(r.text or r.content)
        version = item['info']['version']
        #version = sorted(item['releases'].keys())[-1]
        home = item['info']['home_page']
        license = 'MIT'
        summary = item['info']['summary']
        
        for u in item['urls']:
            if u['python_version'] != 'source': continue
            fil = u['filename']
            url = u['url']
            md5 = u['md5_digest']
            continue
        
        

if local:
    coolprop_dir = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..'))
    version = open(os.path.join(coolprop_dir,'.version'),'r').read().strip()
    home = "http://www.coolprop.org"
    license = "MIT"
    summary = "Open-source thermodynamic and transport properties database"
    fil = None
    url = None
    md5 = None




f = codecs.open(os.path.join(target_dir,target),mode='wb',encoding='utf-8')
f.write(template.render(
  pypi=pypi,
  local=local,
  version=version,
  fil=fil,
  url=url,
  md5=md5,
  pkgs=pkgs,
  home = home,
  license = license,
  summary = summary
  ))
f.close()

bat_template = """
pushd wrappers\Python
"%PYTHON%" setup.py install
if errorlevel 1 exit 1
popd 

:: Add more build steps here, if they are necessary.

:: See
:: http://docs.continuum.io/conda/build.html
:: for a list of environment variables that are set during the build process.
"""
target = "bld.bat"
f = codecs.open(os.path.join(target_dir,target),mode='wb',encoding='utf-8')
f.write(bat_template)
f.close()

bsh_template = """
#!/bin/bash
pushd wrappers/Python
$PYTHON setup.py install
popd 

# Add more build steps here, if they are necessary.

# See
# http://docs.continuum.io/conda/build.html
# for a list of environment variables that are set during the build process.
"""
target = "build.sh"
f = codecs.open(os.path.join(target_dir,target),mode='wb',encoding='utf-8')
f.write(bat_template)
f.close()

runner_template = """
from __future__ import print_function
import sys, shutil, subprocess, os, errno
def run_command(cmd):
    '''given shell command, returns communication tuple of stdout and stderr'''
    return subprocess.Popen(cmd, 
                            stdout=subprocess.PIPE, 
                            stderr=subprocess.PIPE, 
                            stdin=subprocess.PIPE).communicate()
#subprocess.check_call('conda build .', shell = True, stdout = sys.stdout, stderr = sys.stderr)
filename = str(os.path.abspath(run_command('conda build --output .')[0])).strip()
tar = str(os.path.abspath(os.path.join(os.path.dirname(__file__),'conda','Python_conda'))).strip()
try:
    os.makedirs(tar)
except Exception as e:
    if os.path.isdir(tar): pass
    else: raise
print("Copying: "+str(filename)+" to "+str(tar)) 
shutil.copy(filename,tar)
sys.exit(0)
"""
target = "runner.py"
f = codecs.open(os.path.join(target_dir,target),mode='wb',encoding='utf-8')
f.write(runner_template)
f.close()
sys.exit(0)
