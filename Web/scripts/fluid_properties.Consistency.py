from __future__ import print_function
import os.path
import CoolProp
import subprocess
import sys

web_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
root_dir = os.path.abspath(os.path.join(web_dir, '..'))
fluids_path = os.path.join(web_dir, 'fluid_properties', 'fluids')
plots_path = os.path.join(web_dir, 'fluid_properties', 'fluids', 'Consistencyplots')

template = """from __future__ import division, print_function
import matplotlib
matplotlib.use('Agg') #Force mpl to use a non-GUI backend

import matplotlib.pyplot as plt
from CoolProp.Plots.ConsistencyPlots import ConsistencyFigure

ff = ConsistencyFigure('{fluid:s}')
ff.savefig('{fluid:s}.png', dpi = 30)
ff.savefig('{fluid:s}.pdf')
plt.close()
del ff
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
