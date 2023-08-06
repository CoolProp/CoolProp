from __future__ import print_function
import requests, json, sys, os
from typing import Dict, List


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
    headers = {'user-agent': 'CoolProp/0.0.1'}
    try:
        headers["Authorization"] = "token {0}".format(oauth2_token)
    except:
        pass
    #_log_msg("Sending request to the GitHub API: {0}\n".format(url))
    r = requests.get(url, headers=headers)
    #print(url)
    #print(r.json())
    return r.json()

def _make_request_urllib(url: str):
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
SITE_URL = "https://github.com"
BASE_PATH = os.path.dirname(__file__)

def query_search_api(QUER_VARS: Dict, POST_VARS: Dict = dict(sort="created", order="asc")):
    """Finds all issues that have been closed after the tag_data"""
    # Prepare URL
    BASE_URL = "/".join([SEARCH_URL, "issues"])
    QUER_VARS["repo"] = REPO_NAME
    # Prepare query
    QUER_VARS_LIST = []
    for k,v in QUER_VARS.items():
        QUER_VARS_LIST.append("{0}:{1}".format(k, v))
    QUER_VARS_STRING = "+".join(QUER_VARS_LIST)
    POST_VARS["q"] = QUER_VARS_STRING
    # Prepare loop variables
    _result_dict = None
    if "page" in POST_VARS.keys():
        pageCounter = POST_VARS["page"]
        isRunning = False
    else:
        pageCounter = 1
        isRunning = True
    # Loop over the results
    result_dict = None
    while True:
        POST_VARS["page"] = pageCounter
        POST_VARS_LIST = ["{0}={1}".format(k, v) for k,v in POST_VARS.items()]
        POST_VARS_STRING = "&".join(POST_VARS_LIST)
        REQUEST_URL = BASE_URL + "?" + POST_VARS_STRING
        tmp_dict = _make_request(REQUEST_URL)
        isRunning &= len(tmp_dict["items"]) > 0
        if result_dict is None:
            result_dict = tmp_dict
        elif isRunning:
            result_dict["items"] += tmp_dict["items"]
        pageCounter += 1
        if not isRunning:
            break
    
    return result_dict

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


def get_issues_closed_since(tag_date: str, what: str):
    """Finds all issues that have been closed after the tag_data"""
    POST_VARS = dict(sort="created", order="asc")
    QUER_VARS = {}
    QUER_VARS["closed"] = ">="+tag_date
    QUER_VARS["is"] = what
    return query_search_api(QUER_VARS, POST_VARS)


def check_issues_for_labels_and_milestone(ms: str, _issues_dict: dict):
    """Check whether the issues have the correct milestone information or are labeled for exclusion"""
    _no_label_or_ms = []
    _wrong_milestone = []
    print("Processing {} items".format(len(_issues_dict["items"])))
    for _i in _issues_dict["items"]:
    
        #_thetype = None
        #if "issues" in _i["url"]:
        #    _thetype = "issue"
        #if "pulls" in _i["url"]:
        #    _thetype = "pull request"
        #print("Checking {} #{}".format(_thetype, _i["number"]))
    
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
            _wrong_milestone.append(_num)
        else:
            _no_label_or_ms.append(_num)

    _LINK_URL = "/".join([SITE_URL, REPO_NAME, 'issues'])

    if (len(_no_label_or_ms) + len(_wrong_milestone)) > 0:
        _log_msg("The following issues/pull requests seem to have missing information:\n")
        _log_msg(' - No excluded label and no milestone information ({0}): \n'.format(len(_no_label_or_ms)))
        for _num in _no_label_or_ms:
            _log_msg('     #{0} - {1}/{0} \n'.format(_num, _LINK_URL))
        _log_msg(' - Closed items that belong to a different milestone ({0}): \n'.format(len(_wrong_milestone)))
        for _num in _wrong_milestone:
            _log_msg('     #{0} - {1}/{0} \n'.format(_num, _LINK_URL))
        return False
    return True


def get_milestones(milestone):
    fname = milestone + '-milestones.json'
    if True or not os.path.exists(fname):
        _REQUEST_URL = "/".join([REPO_URL, REPO_NAME, 'milestones'])
        # Find the milestone number for the given name
        milestones_json = _make_request(_REQUEST_URL)
        with open(fname, 'w') as fp:
            fp.write(json.dumps(milestones_json, indent=2))
    with open(fname, 'r') as fp:
        return json.load(fp)


def get_milestone_JSON(milestone: str, what: str):
    # Get the merged pull requests associated with the milestone
    fname = "{}-{}.json".format(milestone, what)
    if True or not os.path.exists(fname):
        POST_VARS = dict(sort="created", order="asc")
        QUER_VARS = {}
        QUER_VARS["milestone"] = milestone
        QUER_VARS["is"] = what
        #QUER_VARS["is"] = "closed"
        #QUER_VARS["state"] = "all"
        result = query_search_api(QUER_VARS, POST_VARS)
        with open(fname, 'w') as fp:
            fp.write(json.dumps(result, indent=2))
    with open(fname, 'r') as fp:
        return json.load(fp)

def get_milestone_items(milestone: str, what: str):
    # Get the items associated with the milestone
    items = get_milestone_JSON(milestone, what)["items"]
    print("Found {} {}s associated with {}".format(len(items), what, milestone))
    for item in items:
        if item["state"] != "closed":
            raise ValueError("This {} is still open: {} - {}".format(what, item["number"], item["title"]))
    return items

def generate_issues(milestone):
    # Get the items associated with the milestone
    issues = get_milestone_items(milestone, "issue")
    rst = 'Issues closed:\n\n' + '\n'.join(['* `#{n:d} <https://github.com/CoolProp/CoolProp/issues/{n:d}>`_ : {t:s}'.format(n=issue['number'], t=issue['title']) for issue in issues])
    return rst


def generate_prs(milestone):
    # Get the items associated with the milestone
    issues = get_milestone_items(milestone, "pr")
    rst = 'Pull requests merged:\n\n' + '\n'.join(['* `#{n:d} <https://github.com/CoolProp/CoolProp/pull/{n:d}>`_ : {t:s}'.format(n=issue['number'], t=issue['title']) for issue in issues])
    return rst


if __name__ == '__main__':

    if len(sys.argv) != 3:
        raise ValueError(
            'This script should be called with a command (check or changelog) and a version number ("v5.3.4") as parameter')


    if sys.argv[1] == "check":
        release_json = get_latest_tag_and_date()

        issues_json = get_issues_closed_since(release_json["tag_date"], what="issue")
        succ = check_issues_for_labels_and_milestone(sys.argv[2], issues_json)
        if succ:
             _log_msg("All issues seem to have the correct labels and milestones, congrats!\n")

        issues_json = get_issues_closed_since(release_json["tag_date"], what="pr")
        succ = check_issues_for_labels_and_milestone(sys.argv[2], issues_json)
        if succ:
             _log_msg("All pull requests seem to have the correct labels and milestones, congrats!\n")

    elif sys.argv[1] == "changelog":
        issues_rst = generate_issues(sys.argv[2])
        print(issues_rst)
        with open(os.path.join(BASE_PATH,"snippet_issues.rst.txt"), 'w') as fp:
            fp.write(issues_rst)

        print('')

        issues_rst = generate_prs(sys.argv[2])
        print(issues_rst)
        with open(os.path.join(BASE_PATH,"snippet_pulls.rst.txt"), 'w') as fp:
            fp.write(issues_rst)




    #sys.exit(0)


    #if len(sys.argv) != 2:
    #    raise ValueError('This script should be called like this: python milestone2rst.py v5')
    #
    #print(generate_issues(sys.argv[1]))
    #print('')
    #print(generate_PR(sys.argv[1]))
