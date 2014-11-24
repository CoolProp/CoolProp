import CoolProp
import os

fluid_template = """.. _fluid_{fluid:s}:

{fluid_stars:s}
{fluid:s}
{fluid_stars:s}

References
==========

{references:s}
{aliases:s}

Fluid Information
=================

.. csv-table::
   :header-rows: 1
   :widths: 40, 60
   :file: {fluid:s}-info.csv
   
REFPROP Validation Data
=======================

.. note::

    You can download the script that generated the following figure here: :download:`(link to script)<REFPROPplots/{fluid:s}.py>`, right-click the link and then save as... or the equivalent in your browser.  You can also download this figure :download:`as a PDF<REFPROPplots/{fluid:s}.pdf>`. 

.. image:: REFPROPplots/{fluid:s}.png

Consistency Plots
=================

"""

table_template = """ Parameter, Value
Molar mass [kg/mol],{mm:s}
CAS, {CAS:s}
ASHRAE, {ASHRAE:s}
Triple point temperature [K],{Tt:s}
Triple point pressure [Pa], {pt:s}
Critical point temperature [K], {Tc:s}
Critical point pressure [Pa], {pc:s}
Critical point density [kg/m3], {rhoc_mass:s}
Critical point density [mol/m3], {rhoc_molar:s}
"""
class FluidInfoTableGenerator(object):
    
    def __init__(self, name):
    
        self.name = name
        
    def write(self, path):
        def tos(n):
            ''' convert number to nicely formatted string '''
            n = str(n)
            if 'e' in n:
                n = n.replace('e',':math:`\times 10^{')
                n += '}`'
            else:
                return n
        molar_mass = CoolProp.CoolProp.PropsSI(self.name,'molemass')
        Tt = CoolProp.CoolProp.PropsSI(self.name,'Ttriple')
        Tc = CoolProp.CoolProp.PropsSI(self.name,'Tcrit')
        pc = CoolProp.CoolProp.PropsSI(self.name,'pcrit')
        pt = CoolProp.CoolProp.PropsSI(self.name,'ptriple')
        rhoc_mass = CoolProp.CoolProp.PropsSI(self.name,'rhomass_critical')
        rhoc_molar = CoolProp.CoolProp.PropsSI(self.name,'rhomolar_critical')
        CAS = CoolProp.CoolProp.get_fluid_param_string(self.name, "CAS")
        ASHRAE = CoolProp.CoolProp.get_fluid_param_string(self.name, "ASHRAE")
            
        args = dict(mm = tos(molar_mass),
                    Tt = tos(Tt),
                    pt = tos(pt),
                    Tc = tos(Tc),
                    rhoc_mass = tos(rhoc_mass),
                    rhoc_molar = tos(rhoc_molar),
                    pc = tos(pc),
                    CAS = tos(CAS),
                    ASHRAE = tos(ASHRAE))
        out = table_template.format(**args)
        
        with open(os.path.join(path, self.name+'-info.csv'),'w') as fp:
            print 'writing', os.path.join(path, self.name+'-info.csv')
            fp.write(out)
            
class FluidGenerator(object):
    def __init__(self, fluid):
        self.fluid = fluid
        
    def write(self, path):
        
        # Write CSV table data for fluid information
        ITG = FluidInfoTableGenerator(self.fluid)
        ITG.write(path)
        
        aliases = ', '.join(['``' + a.strip() + '``' for a in CoolProp.CoolProp.get_fluid_param_string(self.fluid, 'aliases').strip().split(',') if a])
        if aliases:
            aliases = 'Aliases\n=======\n\n'+aliases + '\n'
        
        # Write RST file for fluid
        out = fluid_template.format(aliases = aliases,
                                    fluid = self.fluid,
                                    fluid_stars = '*'*len(self.fluid),
                                    references = ''
                                    )
        
        with open(os.path.join(path, self.fluid+'.rst'),'w') as fp:
            print 'writing', os.path.join(path, self.fluid+'.rst')
            fp.write(out)