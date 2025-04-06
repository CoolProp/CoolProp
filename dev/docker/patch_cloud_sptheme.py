from pathlib import Path
import cloud_sptheme
path = Path(cloud_sptheme.__file__).parent / 'ext' / 'index_styling.py'
assert path.exists()

contents = path.open('r').read().replace("from jinja2 import Markup as literal, escape","import jinja2.utils; literal = jinja2.utils.markupsafe.Markup; escape = jinja2.utils.markupsafe.escape")
with path.open('w') as fp:
    fp.write(contents)
print('finished patching cloud_sptheme')
