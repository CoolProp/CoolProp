import json, sys, glob, os
here = os.path.dirname(os.path.abspath(__file__))
sys.path.append(here+'/..')
from package_json import json_options

for fluid in glob.glob(here+'/fluids/*.json'):

    print(fluid)
    j = json.load(open(fluid, 'r'))

    fp = open(fluid, 'w')
    fp.write(json.dumps(j, **json_options))
    fp.close()
