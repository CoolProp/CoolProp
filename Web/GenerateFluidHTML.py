"""
This script will generate the fluid HTML and RST files in the Fluids folder
"""
from __future__ import print_function

import os,sys,glob
import CoolProp
import subprocess

def generate_notebook(fluid):
    """
    Generate a new copy of the template notebook with the fluid name changed
    """
    lines = open('FluidTemplate.ipynb','r').readlines()
    
    line_indices_with_fluid = []
    for i in range(len(lines)):
        line = lines[i].strip()
        if line.startswith("\"Fluid = \'") and line.endswith('\'\\n\",'):
            line_indices_with_fluid.append(i)
    
    # Make sure only one line matches the above criteria
    assert(len(line_indices_with_fluid)==1)
    i = line_indices_with_fluid[0]
    # Get the current fluid in the line
    current_fluid = lines[i].split("\"Fluid = \'")[1].split('\'\\n\",')[0]
    
    # Actually replace the fluid in the line
    lines[i] = lines[i].replace(current_fluid, fluid)
    
    # Write back to temporary file
    fp = open('temp.ipynb','w')
    fp.write(''.join(lines))
    fp.close()
    
stub_template = """{fluid:s}
{line:s}

View this page as an `IPython notebook <http://nbviewer>`_

.. raw:: html
    :file: {fluid}.html
"""

def make_rst_stub(fluid):
    """
    Make a stub file that will raw import the generated HTML
    """
    s = stub_template.format(fluid = fluid, line = '='*len(fluid))
    fp = open(os.path.join('Fluids',fluid+'.rst'),'w')
    fp.write(s)
    fp.close()
    
def make_html_file(fluid):
    """
    Make the 
    """

    call = 'runipy temp.ipynb Fluids\{fluid:s}.ipynb --quiet --html Fluids\{fluid:s}.html --template output_toggle --skip-exceptions'.format(fluid=fluid)
    print('About to make HTML for '+fluid+'; call:', call)
    subprocess.check_output(call, shell = True)
    
index_template = """{header:s}
{line:s}

.. toctree::
    :maxdepth: 1

    {fluids:s}
"""

def index_file(Fluids, header):
    FluidList = '\n    '.join([Fluid+'.rst' for Fluid in Fluids])
    return index_template.format(fluids = FluidList,
                                 header = header,
                                 line = '='*len(header))

# Make this function do everything so that we can use a multiprocessing Pool
def do_fluid(fluid):
    pass
    
if __name__=='__main__':
    
    if not os.path.exists('Fluids'):
        print('making Fluids folder')
        os.mkdir('Fluids')
        
    pure_fluids,pseudo_pure_fluids = [],[]
    for fluid in sorted(CoolProp.__fluids__):
        make_rst_stub(fluid)
        generate_notebook(fluid)
        make_html_file(fluid)
        if CoolProp.CoolProp.IsFluidType(fluid,'PureFluid'):
            pure_fluids.append(fluid)
        else:
            pseudo_pure_fluids.append(fluid)
    
    fp = open(os.path.join('Fluids','pure_fluids.rst'),'w')
    fp.write(index_file(pure_fluids,"Pure Fluids"))
    fp.close()
    
    fp = open(os.path.join('Fluids','pseudo_pure_fluids.rst'),'w')
    fp.write(index_file(pseudo_pure_fluids,"Pseudo-pure Fluids"))
    fp.close()