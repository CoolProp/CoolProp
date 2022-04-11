import argparse
import re
import requests
from packaging import version
from pathlib import Path

ROOT_DIR = Path(__file__).parent.parent


def parse_pypi_version(pypi=False):
    if pypi:
        response = requests.get('https://pypi.org/pypi/CoolProp/json')
    else:
        response = requests.get('https://test.pypi.org/pypi/CoolProp/json')
    response.raise_for_status()
    data = response.json()
    releases = [version.parse(v) for v in data['releases'].keys()]
    return releases


def parse_cmake_version_info():
    with open(ROOT_DIR / 'CMakeLists.txt', 'r') as f:
        content = f.read()

    no_comments_lines = []
    for line in content.splitlines():
        l = line.strip().split('#')[0]
        if l:
            no_comments_lines.append(l)
    content = "\n".join(no_comments_lines)

    m_major = re.search(r'set\s*\(COOLPROP_VERSION_MAJOR (\d+)\)', content)
    m_minor = re.search(r'set\s*\(COOLPROP_VERSION_MINOR (\d+)\)', content)
    m_patch = re.search(r'set\s*\(COOLPROP_VERSION_PATCH (\d+)\)', content)
    m_rev = re.search(r'set\s*\(COOLPROP_VERSION_REVISION "*(.*?)"*\)', content)

    coolprop_version = ''
    if m_major:
        COOLPROP_VERSION_MAJOR = m_major.groups()[0]
        coolprop_version += COOLPROP_VERSION_MAJOR

    if m_minor:
        COOLPROP_VERSION_MINOR = m_minor.groups()[0]
        coolprop_version += "." + COOLPROP_VERSION_MINOR

    if m_patch:
        COOLPROP_VERSION_PATCH = m_patch.groups()[0]
        coolprop_version += "." + COOLPROP_VERSION_PATCH

    if m_rev:
        COOLPROP_VERSION_REV = m_rev.groups()[0]
        coolprop_version += "-" + COOLPROP_VERSION_REV

    return coolprop_version


def replace_setup_py(new_v: version.Version):
    fp = ROOT_DIR / 'wrappers/Python/setup.py'

    with open(fp, 'r') as f:
        content = f.read()

    with open(fp, 'w') as f:
        f.write(content.replace('version=version,', f"version='{new_v}',"))

    print(f"Replaced version '{new_v}' in {fp}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Find the right version from pypi/testpypi")

    parser.add_argument("--cmake-only", default=False,
                        action='store_true',
                        help="Print the version from cmake only")

    parser.add_argument("--pypi", default=False,
                        action='store_true',
                        help="Check pypi instead of testpypi")

    parser.add_argument("--current", default=False,
                        action='store_true',
                        help="Check current version instead of incrementing by one")

    parser.add_argument("--replace-setup-py", default=False,
                        action='store_true',
                        help="Do replacement in setup.py")

    args = parser.parse_args()
    current_v = parse_cmake_version_info()
    if args.cmake_only:
        print(current_v, end="")
        exit(0)

    current_v = version.Version(current_v)

    releases = parse_pypi_version(pypi=args.pypi)

    matched_releases = [v for v in releases
                        if v.base_version == current_v.base_version]

    new_v = current_v.base_version
    if matched_releases:
        max_v = max(matched_releases)
        if max_v.pre:
            pre_iden, pre_v = max_v.pre
            if args.current:
                new_v += f"{pre_iden}{pre_v}"
            else:
                new_v += f"{pre_iden}{pre_v + 1}"
        else:
            new_v += ".post1"
    else:
        new_v = str(current_v)

    new_v = version.Version(new_v)
    if args.replace_setup_py:
        remote = "PyPi" if args.pypi else "TestPyPi"
        print(f"Found next available version on {remote}: {new_v}")

        replace_setup_py(new_v=new_v)
    else:
        print(new_v, end="")
