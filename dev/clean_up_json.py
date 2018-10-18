import json, sys, glob
sys.path.append('..')
from package_json import json_options

for fluid in glob.glob('fluids/*.json'):

    print(fluid)
    j = json.load(open(fluid,'r'))

    fp = open(fluid,'w')
    fp.write(json.dumps(j,**json_options))
    fp.close()
