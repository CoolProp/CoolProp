
import sys, shutil
sys.path.append('../../dev/TTSE/')
import AbstractStateTTSE

shutil.copy2('../../dev/TTSE/AbstractStateTTSE.py', '../coolprop/speed_script.py')

with open('../coolprop/tabular_data.rst.in', 'w') as fp:
    fp.write(AbstractStateTTSE.generate_rst())
