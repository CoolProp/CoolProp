from __future__ import print_function

import os.path
import CoolProp
import subprocess
import sys

web_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
root_dir = os.path.abspath(os.path.join(web_dir, '..')) 
fluids_path = os.path.join(web_dir,'fluid_properties','fluids')
plots_path = os.path.join(web_dir,'fluid_properties','fluids','REFPROPplots')

template = """
from __future__ import division, print_function

import numpy as np, matplotlib.pyplot as plt

import CoolProp
CP = CoolProp.CoolProp

fluid = '{fluid:s}'

fig, ax = plt.subplots()
plt.ylim(10**-18, 10**2)

if CP.get_fluid_param_string(fluid, "REFPROP_name") == 'N/A': 
    ax.set_xlim(0,1)
    # Not in REFPROP
    xlims = ax.get_xlim()
    ylims = ax.get_ylim()
    ax.plot([xlims[0],xlims[1]],[ylims[0],ylims[1]],lw = 3,c = 'r')
    ax.plot([xlims[0],xlims[1]],[ylims[1],ylims[0]],lw = 3,c = 'r')
    
    x = 0.5
    y = (ylims[0]*ylims[1])**0.5
        
    ax.text(x,y,'Not\\nin\\nREFPROP',ha='center',va ='center',bbox = dict(fc = 'white'))
else:
    RPfluid = 'REFPROP::'  + CP.get_fluid_param_string(fluid, "REFPROP_name")
    
    T = 1.1*CP.PropsSI(fluid, 'Tcrit')
    rhoc = CP.PropsSI(fluid, 'rhomolar_critical')

    # Normal properties
    rho = np.linspace(1e-10, 2*rhoc)
    keys = ['P','V','L','Cpmolar','Cvmolar']
    data = dict()
    for i, key in enumerate(keys):
        RPdata = CP.PropsSI(key, 'T', T, 'Dmolar', rho, RPfluid)
        CPdata = CP.PropsSI(key, 'T', T, 'Dmolar', rho, fluid)
        plt.plot(rho/rhoc, np.abs(RPdata/CPdata-1)*100, label = key, dashes = [1, 1+1.5*i], lw = 1.5)
        
    # Special properties
    rho = np.linspace(1e-10, 2*rhoc)
    keys = ['Hmolar','Smolar']
    data = dict()
    for i,key in enumerate(keys):
        RPdata = CP.PropsSI(key, 'T', T, 'Dmolar', rho, RPfluid) - CP.PropsSI(key, 'T', T, 'Dmolar', 1, RPfluid)
        CPdata = CP.PropsSI(key, 'T', T, 'Dmolar', rho, fluid) - CP.PropsSI(key, 'T', T, 'Dmolar', 1, fluid)
        plt.plot(rho/rhoc, np.abs(RPdata/CPdata-1)*100, label = key, dashes = [1, 1+1.5*i], lw = 1.5)

    ax.legend(loc='best', ncol = 2)
    
plt.xlabel(r'Reduced density [$\\rho/\\rho_c$]')
plt.ylabel(r'Relative deviation $(y_{{CP}}/y_{{RP}}-1)\\times 100$ [%]')

ax.set_yscale('log')
plt.savefig(fluid+'.png', dpi = 100)
plt.savefig(fluid+'.pdf')
plt.close('all')

"""
if not os.path.exists(plots_path):
    os.makedirs(plots_path)
    
for fluid in CoolProp.__fluids__:
    print('fluid:', fluid)
    file_string = template.format(fluid = fluid)
    file_path = os.path.join(plots_path, fluid + '.py')
    print('Writing to', file_path)
    with open(file_path, 'w') as fp:
        fp.write(file_string)
    subprocess.check_call('python "' + fluid + '.py"', cwd = plots_path, stdout = sys.stdout, stderr = sys.stderr, shell = True)