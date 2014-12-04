from __future__ import print_function
import urllib, json, sys

def generate(milestone):
    
    # Find the milestone number for the given name
    milestones_json = json.loads(urllib.urlopen('https://api.github.com/repos/CoolProp/CoolProp/milestones').read())

    # Map between name and number
    title_to_number_map = {stone['title']: stone['number'] for stone in milestones_json}

    # Find the desired number
    number = title_to_number_map[milestone]

    # Get the issues associated with the milestone
    issues = json.loads(urllib.urlopen('https://api.github.com/repos/CoolProp/CoolProp/issues?state=all&milestone='+str(number)).read())

    # Make sure all issues are closed in this milestone
    for issue in issues:
        if issue['state'] != 'closed': raise ValueError('This issue is still open: ' + issue['title'])
        
    rst = 'Issues Closed:\n\n'+'\n'.join(['* `#{n:d} <http://github.com/CoolProp/CoolProp/issues/{n:d}>`_ : {t:s}'.format(n = issue['number'], t = issue['title']) for issue in issues])
    
    return rst

if __name__=='__main__':
    sys.argv += ['v5.0.2']
    if len(sys.argv) != 2:
        raise ValueError('This script should be called like this: python milestone2rst.py v5')
        
    print(generate(sys.argv[1]))