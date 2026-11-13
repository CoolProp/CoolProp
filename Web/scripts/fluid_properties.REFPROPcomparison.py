from __future__ import print_function
import os.path
import CoolProp, CoolProp.CoolProp as CP
import subprocess
import sys

web_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
root_dir = os.path.abspath(os.path.join(web_dir, '..'))
fluids_path = os.path.join(web_dir, 'fluid_properties', 'fluids')
plots_path = os.path.join(web_dir, 'fluid_properties', 'fluids', 'REFPROPplots')

template = """
from __future__ import division, print_function
import matplotlib
matplotlib.use('Agg') # Use a non-GUI backend
import numpy as np, matplotlib.pyplot as plt

import CoolProp
CP = CoolProp.CoolProp

fluid = '{fluid:s}'

fig, ax = plt.subplots()
plt.ylim(10**-18, 10**2)

not_in_REFPROP = False
try:
    if CP.get_fluid_param_string(fluid, "REFPROP_name") == 'N/A':
        not_in_REFPROP = True
    else:
        RPfluid = 'REFPROP::'  + CP.get_fluid_param_string(fluid, "REFPROP_name")
        CAS = CP.get_fluid_param_string(RPfluid, "CAS")
except (RuntimeError,ValueError) as E:
    not_in_REFPROP = True

if not_in_REFPROP:
    ax.set_xlim(0,1)
    xlims = ax.get_xlim()
    ylims = ax.get_ylim()
    ax.plot([xlims[0],xlims[1]],[ylims[0],ylims[1]],lw = 3,c = 'r')
    ax.plot([xlims[0],xlims[1]],[ylims[1],ylims[0]],lw = 3,c = 'r')

    x = 0.5
    y = (ylims[0]*ylims[1])**0.5

    ax.text(x,y,'Not\\nin\\nREFPROP',ha='center',va ='center',bbox = dict(fc = 'white'))
else:
    RPfluid = 'REFPROP::'  + CP.get_fluid_param_string(fluid, "REFPROP_name")
    symbols = ["o", "v", "^", "<", ">","8", "s","p","*","h","H","+","x"]

    T = np.min([1.01*CP.PropsSI(fluid, 'Tcrit'),CP.PropsSI(fluid, 'Tmax')])
    rhoc = CP.PropsSI(fluid, 'rhomolar_critical')

    # Transport properties are only compared when both CoolProp and REFPROP have a
    # model for the fluid.  Probe with a scalar call at a single-phase state; an
    # array call hides the error message.  The match ignores case; any other
    # error is not swallowed.
    no_model_msgs = ('Viscosity model is not available',
                     'Thermal conductivity model is not available',
                     'Transport equations are not available')
    def has_model(key, name):
        try:
            CP.PropsSI(key, 'T', T, 'Dmolar', 0.5*rhoc, name)
        except ValueError as E:
            if any(msg.lower() in str(E).lower() for msg in no_model_msgs):
                return False
            raise
        return True

    allkeys = ['P','V','L','Cpmolar','Cvmolar']
    normalkeys = list(allkeys)
    skipped = dict()
    for key, label in (('V', 'viscosity'), ('L', 'thermal conductivity')):
        missing_in = [n for n, f in (('CoolProp', fluid), ('REFPROP', RPfluid)) if not has_model(key, f)]
        if missing_in:
            normalkeys.remove(key)
            skipped.setdefault(' or '.join(missing_in), []).append(label)

    # REFPROP rejects states above twice its own upper pressure limit, so end the
    # density range where CoolProp reaches that pressure (or at 2*rhoc, if sooner)
    p_limit = 2*CP.PropsSI('pmax', RPfluid)*0.995
    if not 0 < p_limit < np.inf:
        raise ValueError('REFPROP upper pressure limit is not usable for ' + RPfluid + ': ' + str(p_limit))
    rho_max = 2*rhoc
    if CP.PropsSI('P', 'T', T, 'Dmolar', rho_max, fluid) > p_limit:
        lo, hi = 1e-10, rho_max
        for _ in range(60):
            mid = 0.5*(lo + hi)
            if CP.PropsSI('P', 'T', T, 'Dmolar', mid, fluid) > p_limit:
                hi = mid
            else:
                lo = mid
        rho_max = lo

    # Normal properties
    rho = np.linspace(1e-10, rho_max)
    RPdata = CP.PropsSI(normalkeys, 'T', T, 'Dmolar', rho, RPfluid)
    CPdata = CP.PropsSI(normalkeys, 'T', T, 'Dmolar', rho, fluid)
    for i, key in enumerate(normalkeys):
        plt.plot(rho/rhoc, np.abs(RPdata[:,i]/CPdata[:,i]-1)*100, lw = 0, label = key, marker = symbols[allkeys.index(key)%len(symbols)])

    # Special properties
    keys = ['Hmolar','Smolar']
    for i, key in enumerate(keys):
        RPdata = CP.PropsSI(key, 'T', T, 'Dmolar', rho, RPfluid) - CP.PropsSI(key, 'T', T, 'Dmolar', 1, RPfluid)
        CPdata = CP.PropsSI(key, 'T', T, 'Dmolar', rho, fluid) - CP.PropsSI(key, 'T', T, 'Dmolar', 1, fluid)
        plt.plot(rho/rhoc, np.abs(RPdata/CPdata-1)*100, lw = 0, label = key, marker = symbols[(i+len(allkeys))%len(symbols)])

    ax.legend(loc='best', ncol = 2)
    if skipped:
        fig.subplots_adjust(bottom = 0.14 + 0.04*len(skipped))
        note = ['Not plotted: ' + ', '.join(labels) + ' (no model in ' + src + ')' for src, labels in skipped.items()]
        fig.text(0.5, 0.02, '\\n'.join(note), ha = 'center', va = 'bottom', fontsize = 8)

plt.xlabel(r'Reduced density [$\\rho/\\rho_c$]')
plt.ylabel(r'Relative deviation $(y_{{CP}}/y_{{RP}}-1)\\times 100$ [%]')

ax.set_yscale('log')
plt.title('Comparison between CoolProp and REFPROP({rpv:s}) along T = 1.01*Tc')
plt.savefig(fluid+'.png', dpi = 100)
plt.savefig(fluid+'.pdf')
plt.close('all')

"""

if not os.path.exists(plots_path):
    os.makedirs(plots_path)

with open(os.path.join(plots_path, 'matplotlibrc'), 'w') as fp:
    fp.write("backend : agg\n")

for fluid in CoolProp.__fluids__:
    print('fluid:', fluid)
    file_string = template.format(fluid=fluid, rpv=CP.get_global_param_string("REFPROP_version"))
    file_path = os.path.join(plots_path, fluid + '.py')
    print('Writing to', file_path)
    with open(file_path, 'w') as fp:
        fp.write(file_string)
    subprocess.check_call('python "' + fluid + '.py"', cwd=plots_path, stdout=sys.stdout, stderr=sys.stderr, shell=True)
