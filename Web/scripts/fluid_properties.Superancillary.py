from __future__ import print_function
import os.path
import CoolProp
import urllib.request
import subprocess
import sys
from pathlib import Path
from zipfile import ZipFile

web_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
root_dir = os.path.abspath(os.path.join(web_dir, '..'))
fluids_path = os.path.join(web_dir, 'fluid_properties', 'fluids')
plots_path = os.path.join(web_dir, 'fluid_properties', 'fluids', 'Superancillaryplots')

outputversion = '2025.04.26'
if not Path(f'{outputversion}.zip').exists():
    print('Downloading the chebyshev output file to ', Path('.').absolute())
    urllib.request.urlretrieve(f'https://github.com/CoolProp/fastchebpure/archive/refs/tags/{outputversion}.zip', f'{outputversion}.zip')

    with ZipFile(f'{outputversion}.zip') as z:
        z.extractall('.')

    with (Path('.') / f'fastchebpure-{outputversion}' / '.gitignore').open('w') as fp:
        fp.write("*")

outputcheck = (Path('.') / f'fastchebpure-{outputversion}' / 'outputcheck').absolute()


template = r"""
import matplotlib
matplotlib.use('Agg')  # Force mpl to use a non-GUI backend
import matplotlib.pyplot as plt
import json
import numpy as np
import pandas as pd
import CoolProp.CoolProp as CP

CP.set_config_bool(CP.ENABLE_SUPERANCILLARIES, False)

AS = CP.AbstractState('HEOS', "{fluid}")

# Skip pseudo-pure fluids, pure fluids only; pseudo-pure do not have superancillaries
if AS.fluid_param_string("pure") != "true":
    quit()

jEOS = json.loads(CP.get_fluid_param_string("{fluid}", "JSON"))[0]['EOS'][0]
if 'SUPERANCILLARY' not in jEOS:
    fig = plt.figure()
    fig.text(0.5, 0.5, 'Superancillary not available')
    plt.savefig('{fluid:s}.png', dpi = 300)
    plt.savefig('{fluid:s}.pdf')
    plt.close()
    quit()
else:
    jSuper = jEOS['SUPERANCILLARY']

superanc = CP.SuperAncillary(json.dumps(jSuper))
RPname = AS.fluid_param_string("REFPROP_name")

# Load extended precision calcs from the release on github
chk = json.load(open('{outputcheck}/{fluid}_check.json'))
df = pd.DataFrame(chk['data'])
# df.info() # uncomment to see what fields are available

Tcrit_num = AS.get_fluid_parameter_double(0, "SUPERANC::Tcrit_num")
T = df['T / K'].to_numpy()
Theta = (Tcrit_num-T)/Tcrit_num

fig, axes = plt.subplots(3, 1, sharex=True, figsize=(3.5,7))


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
plt.tight_layout(pad=0.2, rect=(0,0,1,0.95))
plt.savefig('{fluid:s}.png', dpi = 300)
plt.savefig('{fluid:s}.pdf')
plt.close()
"""
if not os.path.exists(plots_path):
    os.makedirs(plots_path)

from pathlib import Path

for fluid in CoolProp.__fluids__:
    print('fluid:', fluid)
    file_string = template.format(fluid=fluid, outputcheck=outputcheck)
    file_path = os.path.join(plots_path, fluid + '.py')
    print('Writing to', file_path)
    with open(file_path, 'w') as fp:
        fp.write(file_string)
    print('calling:', 'python "' + fluid + '.py"', 'in', plots_path)
    subprocess.check_call('python -u "' + fluid + '.py"', cwd=plots_path, stdout=sys.stdout, stderr=sys.stderr, shell=True)
