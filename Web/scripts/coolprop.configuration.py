from __future__ import print_function

import CoolProp.CoolProp as CP, json
jj = json.loads(CP.get_config_as_json_string())
with open('../coolprop/configuration_keys.rst.in', 'w') as fp:
    for key in sorted(jj.keys()):
        fp.write('``' + key + '``: ' + CP.config_key_description(key) + '\n\n')
