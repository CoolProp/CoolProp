from __future__ import print_function
import os.path
import CoolProp
import subprocess
import sys

web_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
root_dir = os.path.abspath(os.path.join(web_dir, '..'))
fluids_path = os.path.join(web_dir, 'fluid_properties', 'fluids')
plots_path = os.path.join(web_dir, 'fluid_properties', 'fluids', 'Superancillaryplots')

template = r"""
import matplotlib
matplotlib.use('Agg')  # Force mpl to use a non-GUI backend
import matplotlib.pyplot as plt
import json
import numpy as np
import pandas as pd
import CoolProp.CoolProp as CP

CP.set_config_bool(CP.ENABLE_SUPERANCILLARIES, False)

from pathlib import Path
CP.set_config_string(CP.ALTERNATIVE_REFPROP_PATH, str(Path('~/REFPROP10').expanduser()))

AS = CP.AbstractState('HEOS', "{fluid}")
jSuper = json.loads(CP.get_fluid_param_string("{fluid}", "JSON"))[0]['EOS'][0]['SUPERANCILLARY']
superanc = CP.SuperAncillary(json.dumps(jSuper))
RPname = AS.fluid_param_string("REFPROP_name")

# Load extended precision calcs from the release on github
chk = json.load(open(f'/Users/ianbell/Documents/Code/fastchebpure/outputcheck/{fluid}_check.json'))
df = pd.DataFrame(chk['data'])
# df.info() # uncomment to see what fields are available

Tcrit_num = AS.get_fluid_parameter_double(0, "SUPERANC::Tcrit_num")
T = df['T / K'].to_numpy()
Theta = (Tcrit_num-T)/Tcrit_num

fig, axes = plt.subplots(3, 1, sharex=True, figsize=(6,12))


plt.sca(axes[0])
rhoL_anc = np.zeros_like(T)
superanc.eval_sat_many(T, 'D', 0, rhoL_anc)
err = np.abs(df["rho'(mp) / mol/m^3"]/rhoL_anc-1)
plt.plot(Theta, err, 'o', label=r'$\Upsilon$:SA')

errCP = np.abs(df["rho'(mp) / mol/m^3"]/CP.PropsSI('Dmolar', 'T', T, 'Q', 0, f'HEOS::{fluid}')-1)
plt.plot(Theta, errCP, '^', label=r'$\Upsilon$:HEOS')

try:
    errRP = np.abs(df["rho'(mp) / mol/m^3"]/CP.PropsSI('Dmolar', 'T', T, 'Q', 0, f'REFPROP::{{RPname}}')-1)
    plt.plot(Theta, errRP, 'x', label=r'$\Upsilon$:REFPROP')
except BaseException as BE:
    print(BE)

plt.legend(loc='best')
plt.ylabel(r"$|\rho_{{\rm \Upsilon}}'/\rho_{{\rm ep}}'-1|$")
plt.yscale('log')



plt.sca(axes[1])
rhoV_anc = np.zeros_like(T)
superanc.eval_sat_many(T, 'D', 1, rhoV_anc)
err = np.abs(df["rho''(mp) / mol/m^3"]/rhoV_anc-1)
plt.plot(Theta, err, 'o', label=r'$\Upsilon$:SA')

errCP = np.abs(df["rho''(mp) / mol/m^3"]/CP.PropsSI('Dmolar', 'T', T, 'Q', 1, f'HEOS::{fluid}')-1)
plt.plot(Theta, errCP, '^', label=r'$\Upsilon$:HEOS')

try:
    errRP = np.abs(df["rho''(mp) / mol/m^3"]/CP.PropsSI('Dmolar', 'T', T, 'Q', 1, f'REFPROP::{{RPname}}')-1)
    plt.plot(Theta, errRP, 'x', label=r'$\Upsilon$:REFPROP')
except BaseException as BE:
    print(BE)

plt.legend(loc='best')
plt.ylabel(r"$|\rho_{{\rm \Upsilon}}''/\rho_{{\rm ep}}''-1|$")
plt.yscale('log')



plt.sca(axes[2])
p_anc = np.zeros_like(T)
superanc.eval_sat_many(T, 'P', 1, p_anc)
err = np.abs(df["p(mp) / Pa"]/p_anc-1)
plt.plot(Theta, err, 'o', label=r'$\Upsilon$:SA')

errCP = np.abs(df["p(mp) / Pa"]/CP.PropsSI('P', 'T', T, 'Q', 1, f'HEOS::{fluid}')-1)
plt.plot(Theta, errCP, '^', label=r'$\Upsilon$:HEOS')

try:
    errRP = np.abs(df["p(mp) / Pa"]/CP.PropsSI('P', 'T', T, 'Q', 1, f'REFPROP::{{RPname}}')-1)
    plt.plot(Theta, errRP, 'x', label=r'$\Upsilon$:REFPROP')
except BaseException as BE:
    print(BE)
plt.legend(loc='best')

# print(CP.PropsSI('gas_constant', 'T', T[0], 'Q', 1, f'HEOS::{fluid}'))
# print(CP.PropsSI('gas_constant', 'T', T[0], 'Q', 1, f'REFPROP::{fluid}'))

plt.ylabel(r"$|p_{{\rm \Upsilon}}/p_{{\rm ep}}-1|$")
plt.yscale('log')

plt.sca(axes[2])
plt.xlabel(r'$\Theta=(T_{{\rm crit,num}}-T)/T_{{\rm crit,num}}$')
plt.xscale('log')

plt.suptitle('Superancillary v. Extended Precision')
plt.savefig('{fluid:s}.png', dpi = 30)
plt.savefig('{fluid:s}.pdf')
plt.close()
"""
if not os.path.exists(plots_path):
    os.makedirs(plots_path)

for fluid in CoolProp.__fluids__:
    print('fluid:', fluid)
    file_string = template.format(fluid=fluid)
    file_path = os.path.join(plots_path, fluid + '.py')
    print('Writing to', file_path)
    with open(file_path, 'w') as fp:
        fp.write(file_string)
    print('calling:', 'python "' + fluid + '.py"', 'in', plots_path)
    subprocess.check_call('python -u "' + fluid + '.py"', cwd=plots_path, stdout=sys.stdout, stderr=sys.stderr, shell=True)
