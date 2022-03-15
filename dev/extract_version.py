import re
from pathlib import Path

ROOT_DIR = Path(__file__).parent.parent

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

print(f"{coolprop_version}")
