from __future__ import print_function
import urllib, json, sys

def generate_issues(milestone):
    
    # Find the milestone number for the given name
    milestones_json = json.loads(urllib.urlopen('https://api.github.com/repos/CoolProp/CoolProp/milestones').read())

    print(milestones_json)

    # Map between name and number
    title_to_number_map = {stone['title']: stone['number'] for stone in milestones_json}

    # Find the desired number
    number = title_to_number_map[milestone]

    # Get the issues associated with the milestone
    issues = json.loads(urllib.urlopen('https://api.github.com/repos/CoolProp/CoolProp/issues?state=all&per_page=1000&milestone='+str(number)).read())

    # Make sure all issues are closed in this milestone
    for issue in issues:
        if issue['state'] != 'closed': raise ValueError('This issue is still open: ' + issue['title'])
        
    rst = 'Issues Closed:\n\n'+'\n'.join(['* `#{n:d} <https://github.com/CoolProp/CoolProp/issues/{n:d}>`_ : {t:s}'.format(n = issue['number'], t = issue['title']) for issue in issues])
    
    return rst

def generate_PR(milestone):
    
    # Find the milestone number for the given name
    milestones_json = json.loads(urllib.urlopen('https://api.github.com/repos/CoolProp/CoolProp/milestones').read())

    # Map between name and number
    title_to_number_map = {stone['title']: stone['number'] for stone in milestones_json}

    # Find the desired number
    number = title_to_number_map[milestone]

    # Get the merged pull requests associated with the milestone
    PR = json.loads(urllib.urlopen('https://api.github.com/repos/CoolProp/CoolProp/pulls?state=closed&per_page=1000&milestone='+str(number)).read())

    for issue in PR:
        print(issue)
        
    rst = 'Pull Requests merged:\n\n'+'\n'.join(['* `#{n:d} <https://github.com/CoolProp/CoolProp/pull/{n:d}>`_ : {t:s}'.format(n = issue['number'], t = issue['title'].encode('utf-8')) for issue in PR if issue['milestone']['title'] == milestone])
    
    return rst

if __name__=='__main__':
    if len(sys.argv) != 2:
        raise ValueError('This script should be called like this: python milestone2rst.py v5')
        
    print(generate_issues(sys.argv[1]))
    print('')
    print(generate_PR(sys.argv[1]))