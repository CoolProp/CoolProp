import CoolProp
CP = CoolProp.CoolProp
import os
import codecs

web_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
root_dir = os.path.abspath(os.path.join(web_dir, '..'))

fluid_template = u""".. _fluid_{fluid:s}:

{fluid_stars:s}
{fluid:s}
{fluid_stars:s}

{references:s}
{aliases:s}

Fluid Information
=================

.. csv-table::
   :header-rows: 1
   :widths: 40, 60
   :delim: ;
   :file: {fluid:s}-info.csv

REFPROP Validation Data
=======================

.. note::

    This figure compares the results generated from CoolProp and those generated from REFPROP.  They are all results obtained in the form :math:`Y(T,\\rho)`, where :math:`Y` is the parameter of interest and which for all EOS is a direct evaluation of the EOS

    You can download the script that generated the following figure here: :download:`(link to script)<REFPROPplots/{fluid:s}.py>`, right-click the link and then save as... or the equivalent in your browser.  You can also download this figure :download:`as a PDF<REFPROPplots/{fluid:s}.pdf>`.

.. image:: REFPROPplots/{fluid:s}.png

Consistency Plots
=================

The following figure shows all the flash routines that are available for this fluid.  A red + is a failure of the flash routine, a black dot is a success.  Hopefully you will only see black dots.  The red curve is the maximum temperature curve, and the blue curve is the melting line if one is available for the fluid.

In this figure, we start off with a state point given by T,P and then we calculate each of the other possible output pairs in turn, and then try to re-calculate T,P from the new input pair.  If we don't arrive back at the original T,P values, there is a problem in the flash routine in CoolProp.  For more information on how these figures were generated, see :py:mod:`CoolProp.Plots.ConsistencyPlots`

.. note::

    You can download the script that generated the following figure here: :download:`(link to script)<Consistencyplots/{fluid:s}.py>`, right-click the link and then save as... or the equivalent in your browser.  You can also download this figure :download:`as a PDF<Consistencyplots/{fluid:s}.pdf>`.

.. image:: Consistencyplots/{fluid:s}.png

"""

table_template = """ Parameter, Value
**General**;
Molar mass [kg/mol];{mm:s}
CAS number; {CAS:s}
ASHRAE class; {ASHRAE:s}
Formula; {formula:s}
Acentric factor; {acentric:s}
InChI; {inchi:s}
InChIKey; {inchikey:s}
SMILES; {smiles:s}
ChemSpider ID; {ChemSpider_id:s}
2D image; .. image:: {twoDurl:s}
**Limits**;
Maximum temperature [K];{Tmax:s}
Maximum pressure [Pa];{pmax:s}
**Triple point**;
Triple point temperature [K];{Tt:s}
Triple point pressure [Pa]; {pt:s}
**Critical point**;
Critical point temperature [K]; {Tc:s}
Critical point density [kg/m3]; {rhoc_mass:s}
Critical point density [mol/m3]; {rhoc_molar:s}
Critical point pressure [Pa]; {pc:s}
{reducing_string:s}
"""

reducing_template = """**Reducing point**;
Reducing point temperature [K]; {Tr:s}
Reducing point density [mol/m3]; {rhor_molar:s}
"""

bibtex_keys = ['EOS', 'CP0', 'CONDUCTIVITY', 'VISCOSITY', 'MELTING_LINE', 'SURFACE_TENSION']
bibtex_map = {'EOS': 'Equation of State',
              'CP0': 'Ideal gas specific heat',
              'CONDUCTIVITY': 'Thermal Conductivity',
              'VISCOSITY': 'Viscosity',
              'MELTING_LINE': 'Melting Line',
              'SURFACE_TENSION': 'Surface Tension'}

from pybtex.database.input import bibtex
parser = bibtex.Parser()
bibdata = parser.parse_file(os.path.join(root_dir, "CoolPropBibTeXLibrary.bib"))

from CoolProp.BibtexParser import BibTeXerClass
BTC = BibTeXerClass(os.path.join(root_dir, "CoolPropBibTeXLibrary.bib"))

# See http://stackoverflow.com/questions/19751402/does-pybtex-support-accent-special-characters-in-bib-file/19754245#19754245
import pybtex
style = pybtex.plugin.find_plugin('pybtex.style.formatting', 'plain')()
backend = pybtex.plugin.find_plugin('pybtex.backends', 'html')()
parser = pybtex.database.input.bibtex.Parser()


def entry2html(entry):
    for e in entry:
        return e.text.render(backend).replace('{', '').replace('}', '').replace('\n', ' ')


def generate_bibtex_string(fluid):
    string = ''
    for key in bibtex_keys:
        header_string = ''
        sect_strings = []
        try:
            # get the item
            bibtex_key = CoolProp.CoolProp.get_BibTeXKey(fluid, key).strip()
            for thekey in bibtex_key.split(','):
                if thekey.strip() in bibdata.entries.keys():
                    html = BTC.getEntry(key=thekey.strip(), fmt='html')
                    if len(sect_strings) == 0:
                        sect = bibtex_map[key]
                        header_string = sect + '\n' + '-' * len(sect) + '\n\n'
                    sect_strings.append('.. raw:: html\n\n   ' + html + '\n\n')
        except ValueError as E:
            print("error: %s" % E)
        string += header_string + '\n\n.. raw:: html\n\n    <br><br> \n\n'.join(sect_strings)
    return string


class FluidInfoTableGenerator(object):

    def __init__(self, name):

        self.name = name

    def write(self, path):
        def tos(n):
            ''' convert number to nicely formatted string '''
            n = str(n)
            if 'e' in n:
                n = n.replace('e', ':math:`\times 10^{')
                n += '}`'
            else:
                return n
            return n
        molar_mass = CoolProp.CoolProp.PropsSI(self.name, 'molemass')
        Tt = CoolProp.CoolProp.PropsSI(self.name, 'Ttriple')
        Tc = CoolProp.CoolProp.PropsSI(self.name, 'Tcrit')
        Tr = CoolProp.CoolProp.PropsSI(self.name, 'T_reducing')
        pc = CoolProp.CoolProp.PropsSI(self.name, 'pcrit')
        pt = CoolProp.CoolProp.PropsSI(self.name, 'ptriple')
        if pt is None:
            pt = "Unknown"
        Tmax = CoolProp.CoolProp.PropsSI(self.name, 'Tmax')
        pmax = CoolProp.CoolProp.PropsSI(self.name, 'pmax')
        acentric = CoolProp.CoolProp.PropsSI(self.name, 'acentric')
        rhoc_mass = CoolProp.CoolProp.PropsSI(self.name, 'rhomass_critical')
        rhoc_molar = CoolProp.CoolProp.PropsSI(self.name, 'rhomolar_critical')
        rhor_molar = CoolProp.CoolProp.PropsSI(self.name, 'rhomolar_reducing')

        CAS = CoolProp.CoolProp.get_fluid_param_string(self.name, "CAS")
        ASHRAE = CoolProp.CoolProp.get_fluid_param_string(self.name, "ASHRAE34")
        formula = CoolProp.CoolProp.get_fluid_param_string(self.name, "formula")
        if formula:
            formula = ':math:`' + formula + '`'
        else:
            formula = 'Not applicable'
        formula = formula.replace('_{1}', '')
        InChI = CoolProp.CoolProp.get_fluid_param_string(self.name, "INCHI")
        InChiKey = CoolProp.CoolProp.get_fluid_param_string(self.name, "INCHIKEY")
        smiles = CoolProp.CoolProp.get_fluid_param_string(self.name, "SMILES")
        ChemSpider_id = CoolProp.CoolProp.get_fluid_param_string(self.name, "CHEMSPIDER_ID")
        twoDurl = CoolProp.CoolProp.get_fluid_param_string(self.name, "2DPNG_URL")

        # Generate (or not) the reducing data
        reducing_data = ''
        if abs(Tr - Tc) > 1e-3:
            reducing_data = reducing_template.format(Tr=tos(Tr),
                                                     rhor_molar=tos(rhor_molar))

        args = dict(mm=tos(molar_mass),
                    Tt=tos(Tt),
                    pt=tos(pt),
                    Tc=tos(Tc),
                    rhoc_mass=tos(rhoc_mass),
                    rhoc_molar=tos(rhoc_molar),
                    pc=tos(pc),
                    acentric=tos(acentric),
                    CAS=tos(CAS),
                    ASHRAE=tos(ASHRAE),
                    Tmax=tos(Tmax),
                    pmax=tos(pmax),
                    reducing_string=reducing_data,
                    formula=formula,
                    inchi=InChI,
                    inchikey=InChiKey,
                    smiles=smiles,
                    ChemSpider_id=ChemSpider_id,
                    twoDurl=twoDurl
                    )
        out = table_template.format(**args)

        with open(os.path.join(path, self.name + '-info.csv'), 'w') as fp:
            print("writing %s" % os.path.join(path, self.name + '-info.csv'))
            fp.write(out)


class FluidGenerator(object):
    def __init__(self, fluid):
        self.fluid = fluid

    def write(self, path):

        # Write CSV table data for fluid information
        ITG = FluidInfoTableGenerator(self.fluid)
        ITG.write(path)

        del_old = CP.get_config_string(CP.LIST_STRING_DELIMITER)

        CP.set_config_string(CP.LIST_STRING_DELIMITER, '|')
        try:
            aliases = ', '.join(['``' + a.strip() + '``' for a in CoolProp.CoolProp.get_fluid_param_string(self.fluid, 'aliases').strip().split('|') if a])
        finally:
            CP.set_config_string(CP.LIST_STRING_DELIMITER, del_old)

        if aliases:
            aliases = 'Aliases\n=======\n\n' + aliases + '\n'

        references = generate_bibtex_string(self.fluid)
        if references:
            references = 'References\n==========\n' + references + '\n'

        # Write RST file for fluid
        out = fluid_template.format(aliases=aliases,
                                    fluid=self.fluid,
                                    fluid_stars='*' * len(self.fluid),
                                    references=references
                                    )

        with codecs.open(os.path.join(path, self.fluid + '.rst'), 'w', encoding='utf-8') as fp:
            print("writing %s" % os.path.join(path, self.fluid + '.rst'))
            fp.write(out)
