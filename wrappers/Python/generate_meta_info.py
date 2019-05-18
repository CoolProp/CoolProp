from jinja2 import Environment
import os, sys
import requests
import json
from distutils.version import LooseVersion  # , StrictVersion
import codecs

""" A simple script to create a conda recipe and the infrastructure files for PyPI"""

first_line = "# CAUTION: This file is generated automatically, any customisation will be lost.\n"

python_dir = os.path.abspath(os.path.dirname(__file__))
target_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))

pypi = False
local = not pypi

run_pkgs = ["numpy", "scipy", "matplotlib", "pandas"]
dev_pkgs = run_pkgs + ["cython"]
tst_pkgs = dev_pkgs + ["nose"]

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
    coolprop_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    sys.path.append(os.path.join(coolprop_dir, 'dev'))
    import generate_headers
    # Generate the headers - does nothing if up to date - but only if not pypi
    generate_headers.generate()
    del generate_headers
    version = open(os.path.join(coolprop_dir, '.version'), 'r').read().strip()
    home = "http://www.coolprop.org"
    license = "MIT"
    summary = "Open-source thermodynamic and transport properties database"
    fil = None
    url = None
    md5 = None

local_info = dict(
  pypi=pypi,
  local=local,
  version=version,
  fil=fil,
  url=url,
  md5=md5,
  run_pkgs=run_pkgs,
  dev_pkgs=dev_pkgs,
  tst_pkgs=tst_pkgs,
  home=home,
  license=license,
  summary=summary
)

#######################


template = """{% for pkg in run_pkgs %}{{ pkg }}
{% endfor %}
"""
target = "requirements.txt"
template = Environment().from_string(template)
f = codecs.open(os.path.join(python_dir, target), mode='wb', encoding='utf-8')
f.write(first_line)
f.write(template.render(**local_info))
f.close()


template = """
package:
  name: coolprop
  version: {{ version }}

{% if pypi %}source:
  fn: {{ fil }}
  url: {{ url }}
  md5: {{ md5 }}
{% endif %}
{% if local %}source:
  path: .
{% endif %}

#build:
#  script: python setup.py install [not win]
#  script: "%PYTHON%" setup.py install & if errorlevel 1 exit 1 [win]

  # If this is a new build for the same version, increment the build
  # number. If you do not include this key, it defaults to 0.
  # number: 1

requirements:
  build:
    - python
    - setuptools{% for pkg in dev_pkgs %}
    - {{ pkg -}}
{% endfor %}

  run:
    - python{% for pkg in run_pkgs %}
    - {{ pkg -}}
{% endfor %}

test:
  # Python imports
  imports:
    - CoolProp
#    #- CoolProp.GUI
#    #- CoolProp.Plots
#    - CoolProp.tests

  # commands:
    # You can put test commands to be run here.  Use this to test that the
    # entry points work.


  # You can also put a file called run_test.py in the recipe that will be run
  # at test time.

  requires:
    # Put any additional test requirements here.  For example
    # - nose{% for pkg in tst_pkgs %}
    - {{ pkg -}}
{% endfor %}

about:
  home: {{ home }}
  license: {{ license }}
  summary: {{ summary }}

"""
target = 'meta.yaml'
template = Environment().from_string(template)
f = codecs.open(os.path.join(target_dir, target), mode='wb', encoding='utf-8')
f.write(first_line)
f.write(template.render(**local_info))
f.close()

template = """
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
f = codecs.open(os.path.join(target_dir, target), mode='wb', encoding='utf-8')
f.write(":: " + first_line)
f.write(template)
f.close()

template = """
pushd wrappers/Python
$PYTHON setup.py install
popd

# Add more build steps here, if they are necessary.

# See
# http://docs.continuum.io/conda/build.html
# for a list of environment variables that are set during the build process.
"""
target = "build.sh"
f = codecs.open(os.path.join(target_dir, target), mode='wb', encoding='utf-8')
f.write("#!/bin/bash\n" + first_line)
f.write(template)
f.close()


template = """
from __future__ import print_function
import sys, shutil, subprocess, os, stat
#
def run_command(cmd):
    '''given shell command, returns communication tuple of stdout and stderr'''
    print(str(__file__)+": "+' '.join(cmd))
    return subprocess.Popen(cmd,
      stdout=subprocess.PIPE,
      stderr=subprocess.PIPE,
      stdin=subprocess.PIPE).communicate()
#
tar = os.path.abspath(os.path.join(os.path.dirname(__file__),'install_root')).strip()
ver =  sys.version_info
cmd = ['conda','build','--python',str(ver[0])+'.'+str(ver[1])]
print(run_command(['conda', 'clean', '-y', '-lts'])[0].decode("utf-8").strip())
filename = os.path.abspath(run_command(cmd+['--output','.'])[0].decode("utf-8").strip())
tar = os.path.join(tar,'Python_conda',os.path.basename(os.path.dirname(filename))).strip()
try:
    subprocess.check_call(cmd+['.'], stdout=sys.stdout, stderr=sys.stderr)
except Exception as e:
    print("conda build failed: "+str(e))
    pass
try:
    os.makedirs(tar)
except Exception as e:
    if os.path.isdir(tar): pass
    else: raise
try:
    print("Copying: "+str(filename)+" to "+str(tar))
    shutil.copy(filename,tar)
except Exception as e:
    print("Copy operation failed: "+str(e))
    pass
sys.exit(0)
"""
target = "runner.py"
f = codecs.open(os.path.join(target_dir, target), mode='wb', encoding='utf-8')
f.write(template)
f.close()
sys.exit(0)
