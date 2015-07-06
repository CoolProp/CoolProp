import jinja2
import os
import requests
import json
from distutils.version import LooseVersion #, StrictVersion
import codecs

""" A simple script to create a conda recipe from the PyPI release and build the package"""

template_dir = os.path.dirname(__file__)
target_dir = os.path.abspath(os.path.join(os.path.dirname(__file__),'coolprop'))

pkgs = ["numpy", "scipy", "matplotlib", "pandas", "cython"]

loader = jinja2.FileSystemLoader(template_dir)
environment = jinja2.Environment(loader=loader)

target = 'meta.yaml'
template = environment.get_template(os.path.join(template_dir,'conda_'+target+'.tpl'))

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
    
    f = codecs.open(os.path.join(target_dir,target),mode='wb',encoding='utf-8')
    f.write(template.render(
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
        
    
    # release = item['releases'][version]
    
    # for d in release:
        # if d['python_version'] == 'source':
            # fil = 
            # url =
            # md5 = d['md5_digest']

    

    

##tag = sorted(tags.keys())[-1]

##for tag in sorted(tags.keys()):
##    print tag
##     r = requests.get('https://api.github.com/repos/coolprop/coolprop/git/tags/'+tags[tag])
##     if(r.ok):
##         items = json.loads(r.text or r.content)
##         print str(items)

##def cmp(x,y): return LooseVersion(x).__cmp__(y)
##tag = sorted(tags.keys(),cmp=cmp)[-1]
#tag = sorted(tags.keys())[-1]
##from pkg_resources import parse_version
##>>> parse_version('1.4') > parse_version('1.4-rc2')
#if tag[0]=='v': version = tag[1:]
#else: version = tag

#f = codecs.open(os.path.join(target_dir,target),mode='wb',encoding='utf-8')
##f = open(name,mode='w')
#f.write(template.render(version=version,tag=tag,pkgs=pkgs))
#f.close()