from __future__ import print_function

import CoolProp.CoolProp as CP, json
jj = json.loads(CP.get_config_as_json_string())
# UTF-8 encoding gives us headroom for non-ASCII characters in
# config_key_description (em-dashes, units, etc.) without depending
# on the docs container's locale being UTF-8 by default.
with open('../coolprop/configuration_keys.rst.in', 'w', encoding='utf-8') as fp:
    for key in sorted(jj.keys()):
        fp.write('``' + key + '``: ' + CP.config_key_description(key) + '\n\n')
