from __future__ import print_function
import urllib, json, sys, os

try:
    from urllib.request import urlopen
except ImportError:
    from urllib import urlopen

def get_milestones(milestone):
    fname = milestone+'-milestones.json'
    if not os.path.exists(fname):
        # Find the milestone number for the given name
        milestones_json = json.loads(urlopen('https://api.github.com/repos/CoolProp/CoolProp/milestones').read())
        with open(fname,'w') as fp:
            fp.write(json.dumps(milestones_json, indent = 2))
    with open(fname,'r') as fp:
        return json.load(fp)


def get_PR_JSON(milestone, number):
    # Get the merged pull requests associated with the milestone
    fname = milestone+'-PR.json'
    if not os.path.exists(fname):
        # Find the milestone number for the given name
        PR = json.loads(urlopen('https://api.github.com/repos/CoolProp/CoolProp/pulls?state=closed&per_page=1000&milestone='+str(number)).read())
        with open(fname,'w') as fp:
            fp.write(json.dumps(PR, indent = 2))
    with open(fname,'r') as fp:
        return json.load(fp)


def get_issues_JSON(milestone, number):
    # Get the merged pull requests associated with the milestone
    fname = milestone+'-issues.json'
    if not os.path.exists(fname):
        # Find the milestone number for the given name
        issues = json.loads(urlopen('https://api.github.com/repos/CoolProp/CoolProp/issues?state=all&per_page=1000&milestone='+str(number)).read())
        with open(fname,'w') as fp:
            fp.write(json.dumps(issues, indent = 2))
    with open(fname,'r') as fp:
        return json.load(fp)


def generate_issues(milestone):

    milestones_json = get_milestones(milestone)

    # Map between name and number
    title_to_number_map = {stone['title']: stone['number'] for stone in milestones_json}

    # Find the desired number
    number = title_to_number_map[milestone]

    PR = get_PR_JSON(milestone, number)
    pr_numbers = [issue['number'] for issue in PR]

    # Get the issues associated with the milestone
    issues = get_issues_JSON(milestone,     number)

    # Make sure all issues are closed in this milestone
    l = 0
    for issue in issues:
        if issue['state'] != 'closed': raise ValueError('This issue is still open: ' + issue['title'])
        l = l + 1

    for i in reversed(range(l)):
        if issues[i]['number'] in pr_numbers:
            issues.pop(i)

    rst = 'Issues Closed:\n\n'+'\n'.join(['* `#{n:d} <https://github.com/CoolProp/CoolProp/issues/{n:d}>`_ : {t:s}'.format(n = issue['number'], t = issue['title']) for issue in issues])

    return rst


def generate_PR(milestone):

    # Find the milestone number for the given name
    milestones_json = get_milestones(milestone)

    # Map between name and number
    title_to_number_map = {stone['title']: stone['number'] for stone in milestones_json}

    # Find the desired number
    number = title_to_number_map[milestone]

    PR = get_PR_JSON(milestone, number)

    rst = 'Pull Requests merged:\n\n'
    for issue in PR:
        if issue['milestone'] is not None and issue['milestone']['title'] == milestone:
            rst += '* `#{n:d} <https://github.com/CoolProp/CoolProp/pull/{n:d}>`_ : {t:s}\n'.format(n = issue['number'], t = issue['title'])

    return rst


if __name__=='__main__':
    if len(sys.argv) != 2:
        raise ValueError('This script should be called like this: python milestone2rst.py v5')

    print(generate_issues(sys.argv[1]))
    print('')
    print(generate_PR(sys.argv[1]))
