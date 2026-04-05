from __future__ import print_function

import json
import os.path
import re
import sys

import CoolProp.CoolProp as CP

web_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
json_path = os.path.abspath(os.path.join(web_dir, '..', 'dev', 'mixtures', 'predefined_mixtures.json'))
csvfile = os.path.join(web_dir, 'fluid_properties', 'PredefinedMixtures.csv')

with open(json_path, 'r') as f:
    mixtures = json.load(f)

mixtures = sorted(mixtures, key=lambda m: m['name'])


def check_mixture(name):
    """Return None if the mixture loads successfully, or a short error string."""
    try:
        CP.PropsSI('H', 'T', 300, 'P', 101325, name + '.mix')
        return None
    except Exception as e:
        msg = str(e)
        if 'was not found' in msg:
            return 'Not found in predefined mixture library'
        m = re.search(r'binary pair \[([^\]]+)\]', msg)
        if m:
            return 'Missing binary interaction parameters for pair [{}]'.format(m.group(1))
        # Fallback: first sentence of the error
        return msg.split(';')[0].strip()


errors = 0
with open(csvfile, 'w') as fp:
    fp.write('Name,Components,Mole fractions,Notes\n')
    for mix in mixtures:
        name = mix['name']
        fluids = ' + '.join(mix['fluids'])
        fracs = ' + '.join('{:.6g}'.format(x) for x in mix['mole_fractions'])
        note = check_mixture(name)
        if note:
            errors += 1
            sys.stdout.write('WARN {}: {}\n'.format(name, note))
        fp.write('{},{},{},{}\n'.format(name, fluids, fracs, note or ''))

print('Wrote {} predefined mixtures ({} with errors) to {}'.format(
    len(mixtures), errors, csvfile))
