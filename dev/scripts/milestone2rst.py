from __future__ import print_function
import urllib, json, sys, os

try:
    from urllib.request import urlopen, Request
except ImportError:
    from urllib import urlopen, Request


def _log_msg(msg: str):
    """A pseudo logger that prints to stdout without line breaks"""
    try:
        msg = msg.encode(encoding="ascii", errors="ignore")
        msg = msg.decode()
    except AttributeError:
        pass
    try:
        print(msg, end='', flush=True)
    except Exception as e:
        print(str(e))


def _make_request(url: str):
    """Takes a URL, makes a request and returns a JSON dict"""
    request = Request(url)
    # use a global variable for the authentication token
    try:
        request.add_header("Authorization", "token {0}".format(oauth2_token))
    except:
        pass

    response = urlopen(request).read()

    # Decode bytes to string (using UTF-8)
    try:
        response = response.decode()
    except AttributeError:
        pass

    return json.loads(response, encoding="utf-8")


REPO_NAME = "CoolProp/CoolProp"
REPO_URL = "https://api.github.com/repos"
SEARCH_URL = "https://api.github.com/search"


def get_latest_tag_and_date():
    # Get latest release
    _REQUEST_URL = "/".join([REPO_URL, REPO_NAME, 'releases/latest'])
    _release = _make_request(_REQUEST_URL)
    _tag_name = _release["tag_name"]
    # Get all tags for further processing
    _REQUEST_URL = "/".join([REPO_URL, REPO_NAME, 'tags'])
    _tags = _make_request(_REQUEST_URL)
    # Add dates to commits:
    _tag_date = None
    for _t in _tags:
        if _t["name"] == _tag_name:
            _REQUEST_URL = "/".join([REPO_URL, REPO_NAME, 'commits', _t["commit"]["sha"]])
            _commit = _make_request(_REQUEST_URL)
            _date = _commit["commit"]["author"]["date"]
            _tag_date = _date
            _release["tag_date"] = _tag_date
            break
    _log_msg("Determining last tag/release date using the GitHub API: {0} ({1}) \n".format(_tag_name, _tag_date))
    return _release


def get_issues_closed_since(tag_date: str):
    """Finds all issues that have been closed after the tag_data"""
    BASE_URL = "/".join([SEARCH_URL, "issues"])
    POST_VARS = dict(page=1, per_page=1000, sort="created", order="asc")
    QUER_VARS = dict(repo=REPO_NAME, closed=">="+tag_date)
    QUER_VARS["is"] = "issue"

    QUER_VARS_LIST = []
    for k,v in QUER_VARS.items():
        QUER_VARS_LIST.append("{0}:{1}".format(k, v))
    QUER_VARS_STRING = "+".join(QUER_VARS_LIST)

    POST_VARS_LIST = []
    for k,v in POST_VARS.items():
        POST_VARS_LIST.append("{0}={1}".format(k, v))
    POST_VARS_STRING = "&".join(POST_VARS_LIST)

    _REQUEST_URL = BASE_URL + "?" + POST_VARS_STRING + "&q=" + QUER_VARS_STRING

    _issues_dict = _make_request(_REQUEST_URL)
    return _issues_dict


def check_issues_for_labels_and_milestone(ms: str, _issues_dict: dict):
    """Check whether the issues have the correct milestone information or are labeled for exclusion"""
    _no_label_or_ms = []
    _wrong_milestone = []
    for _i in _issues_dict["items"]:
        _num = _i["number"]
        _labels = [_l["name"] for _l in _i["labels"]]
        _milestone = _i["milestone"]

        if "duplicate" in _labels:
            continue
        if "invalid" in _labels:
            continue
        if "wontfix" in _labels:
            continue

        if _milestone is not None:
            if _milestone["title"] == ms:
                continue
            _wrong_milestone.append(str(_num))
        else:
            _no_label_or_ms.append(str(_num))

    if (len(_no_label_or_ms) + len(_wrong_milestone)) > 0:
        _log_msg("The following issues seem to have missing information:\n")
        _log_msg(' - No excluded label and no milestone information ({0}): \n'.format(len(_no_label_or_ms)))
        for _num in _no_label_or_ms:
            _log_msg('     #{0} - https://github.com/CoolProp/CoolProp/issues/{0} \n'.format(_num))
        _log_msg(' - Closed issues that belong to a different milestone ({0}): \n'.format(len(_wrong_milestone)))
        for _num in _wrong_milestone:
            _log_msg('     #{0} - https://github.com/CoolProp/CoolProp/issues/{0} \n'.format(_num))
        return False
    return True


def get_milestones(milestone):
    fname = milestone + '-milestones.json'
    if not os.path.exists(fname):
        # Find the milestone number for the given name
        milestones_json = json.loads(urlopen(REPO_URL + '/milestones').read())
        with open(fname, 'w') as fp:
            fp.write(json.dumps(milestones_json, indent=2))
    with open(fname, 'r') as fp:
        return json.load(fp)


def get_PR_JSON(milestone, number):
    # Get the merged pull requests associated with the milestone
    fname = milestone + '-PR.json'
    if not os.path.exists(fname):
        # Find the milestone number for the given name
        PR = json.loads(urlopen(REPO_URL + '/pulls?state=closed&per_page=1000&milestone=' + str(number)).read())
        with open(fname, 'w') as fp:
            fp.write(json.dumps(PR, indent=2))
    with open(fname, 'r') as fp:
        return json.load(fp)


def get_issues_JSON(milestone, number):
    # Get the issues associated with the milestone
    fname = milestone + '-issues.json'
    if not os.path.exists(fname):
        # Find the milestone number for the given name
        issues = json.loads(urlopen(REPO_URL + '/issues?state=all&per_page=1000&milestone=' + str(number)).read())
        with open(fname, 'w') as fp:
            fp.write(json.dumps(issues, indent=2))
    with open(fname, 'r') as fp:
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
    issues = get_issues_JSON(milestone, number)

    # Make sure all issues are closed in this milestone
    l = 0
    for issue in issues:
        if issue['state'] != 'closed': raise ValueError('This issue is still open: ' + issue['title'])
        l = l + 1

    for i in reversed(range(l)):
        if issues[i]['number'] in pr_numbers:
            issues.pop(i)

    rst = 'Issues Closed:\n\n' + '\n'.join(['* `#{n:d} <https://github.com/CoolProp/CoolProp/issues/{n:d}>`_ : {t:s}'.format(n=issue['number'], t=issue['title']) for issue in issues])

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
            rst += '* `#{n:d} <https://github.com/CoolProp/CoolProp/pull/{n:d}>`_ : {t:s}\n'.format(n=issue['number'], t=issue['title'])

    return rst


if __name__ == '__main__':

    if len(sys.argv) != 3:
        raise ValueError(
            'This script should be called with a command (check or changelog) and a version number ("v5.3.4") as parameter')


    if sys.argv[1] == "check":
        release_json = get_latest_tag_and_date()
        issues_json = get_issues_closed_since(release_json["tag_date"])
        succ = check_issues_for_labels_and_milestone(sys.argv[2], issues_json)
        if succ:
             _log_msg("All issues seem to have the correct labels and milestones, congrats!\n")
        #else:
        #    _log_msg("There are some issues that seem to have incorrect label and milestone information.\n")

    elif sys.argv[1] == "changelog":
        print(generate_issues(sys.argv[2]))
        print('')
        print(generate_PR(sys.argv[2]))







    #sys.exit(0)


    #if len(sys.argv) != 2:
    #    raise ValueError('This script should be called like this: python milestone2rst.py v5')
    #
    #print(generate_issues(sys.argv[1]))
    #print('')
    #print(generate_PR(sys.argv[1]))
