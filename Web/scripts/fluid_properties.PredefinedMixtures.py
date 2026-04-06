from __future__ import print_function

import json
import os.path
import re
import sys

import CoolProp.CoolProp as CP

web_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
json_path = os.path.abspath(os.path.join(web_dir, '..', 'dev', 'mixtures', 'predefined_mixtures.json'))
csvfile = os.path.join(web_dir, 'fluid_properties', 'PredefinedMixtures.csv')
countfile = os.path.join(web_dir, 'fluid_properties', 'PredefinedMixturesCount.rst')

with open(json_path, 'r') as f:
    mixtures = json.load(f)

mixtures = sorted(mixtures, key=lambda m: m['name'])

# Build a CAS number -> fluid name map for translating binary-pair error messages
_cas_to_name = {}
for _f in CP.FluidsList():
    try:
        _cas_to_name[CP.get_fluid_param_string(_f, 'CAS')] = _f
    except Exception:
        pass

_binary_pair_re = re.compile(r'binary pair \[([^\]]+)\]')
_not_found_key_re = re.compile(r'key \[([^\]]+)\] was not found')


def _pair_names(match):
    """Convert a binary-pair regex match to a human-readable 'A + B' string."""
    return ' + '.join(_cas_to_name.get(c.strip(), c.strip())
                      for c in match.group(1).split(','))


def check_mixture(name, components, fracs):
    """Return None if the mixture loads successfully, or a short error string."""
    try:
        CP.PropsSI('H', 'T', 300, 'P', 101325, name + '.mix')
        return None
    except Exception as e:
        msg = str(e)

        # Missing binary interaction parameters reported directly by the .mix loader
        m = _binary_pair_re.search(msg)
        if m:
            return 'Missing binary interaction parameters: {}'.format(_pair_names(m))

        # Predefined mixture definition not found – dig into components for root cause
        if 'was not found' in msg:
            # Check for missing pure fluids first
            missing = [f for f in components
                       if not _fluid_exists(f)]
            if missing:
                return 'Missing pure fluid(s): {}'.format(' + '.join(missing))

            # All pure fluids exist; try loading by components to detect missing BIP
            try:
                AS = CP.AbstractState('HEOS', '&'.join(components))
                AS.set_mole_fractions(fracs)
                AS.update(CP.PT_INPUTS, 101325, 300)
                # Loaded fine by components – just not registered as a named mixture
                return 'Not registered as predefined mixture'
            except Exception as e2:
                msg2 = str(e2)
                m2 = _binary_pair_re.search(msg2)
                if m2:
                    return 'Missing binary interaction parameters: {}'.format(_pair_names(m2))
                m3 = _not_found_key_re.search(msg2)
                if m3:
                    return 'Missing pure fluid: {}'.format(m3.group(1))
                return msg2.split(';')[0].replace(',', ' ').strip()

        # Fallback: first clause of the error
        return msg.split(';')[0].replace(',', ' ').strip()


def _fluid_exists(name):
    try:
        CP.get_fluid_param_string(name, 'name')
        return True
    except Exception:
        return False


errors = 0
with open(csvfile, 'w') as fp:
    fp.write('Name,Components,Mole fractions,Notes\n')
    for mix in mixtures:
        name = mix['name']
        fluids = ' + '.join(mix['fluids'])
        fracs = ' + '.join('{:.6g}'.format(x) for x in mix['mole_fractions'])
        note = check_mixture(name, mix['fluids'], mix['mole_fractions'])
        if note:
            errors += 1
            sys.stdout.write('WARN {}: {}\n'.format(name, note))
        fp.write('{},{},{},{}\n'.format(name, fluids, fracs, note or ''))

with open(countfile, 'w') as fp:
    fp.write('.. |predefined_mixture_count| replace:: {}\n'.format(len(mixtures)))

print('Wrote {} predefined mixtures ({} with errors) to {}'.format(
    len(mixtures), errors, csvfile))
