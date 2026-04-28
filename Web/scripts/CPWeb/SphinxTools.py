import CoolProp
CP = CoolProp.CoolProp
import os
import codecs
import json
import re
from datetime import datetime, timedelta, timezone
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

web_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
root_dir = os.path.abspath(os.path.join(web_dir, '..'))

fluid_template = u""".. _fluid_{fluid:s}:

{fluid_stars:s}
{fluid:s}
{fluid_stars:s}

{references:s}
{aliases:s}

{molecule_viewer_rst:s}Fluid Information
=================

.. csv-table::
   :header-rows: 1
   :escape: @
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

Superancillary Plots
====================

The following figure shows the accuracy of the superancillary functions relative to extended precision calculations carried out in C++ with the teqp library. The results of the iterative calculations with REFPROP and CoolProp are also shown.

.. note::

    You can download the script that generated the following figure here: :download:`(link to script)<Superancillaryplots/{fluid:s}.py>`, right-click the link and then save as... or the equivalent in your browser.  You can also download this figure :download:`as a PDF<Superancillaryplots/{fluid:s}.pdf>`.

.. image:: Superancillaryplots/{fluid:s}.png

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

_INCHIKEY_RE = re.compile(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$')

# Populated during FluidGenerator.write() runs; consumed by
# write_geometry_status_report() to emit a single human-friendly summary
# of which fluids got 3D conformers, 2D fallback, or no structure at all.
GEOMETRY_STATUS = {
    'has_3d': [],            # list of fluid names with 3D conformer
    'has_2d_only': [],       # list of fluid names with 2D fallback only
    'pubchem_missing': [],   # list of (fluid, inchikey) PubChem returned nothing for
    'no_inchikey': [],       # list of fluid names without a valid InChIKey
}

# PubChem cache freshness — the manifest sits next to the SDF files and records
# (timestamp, status) per (inchikey, record_type), so freshness survives git
# checkouts (mtime does not). Stale entries discovered during a build are
# accumulated here for the end-of-build report.
_MANIFEST_FILENAME = 'manifest.json'
_DEFAULT_MAX_AGE_DAYS = 180

CACHE_FRESHNESS = {
    'stale': [],              # list of manifest keys older than max_age_days
    'max_age_days': _DEFAULT_MAX_AGE_DAYS,
    'refresh_mode': '',       # 'all' | 'stale' | '' (resolved at first lookup)
}


def _now_iso():
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _parse_iso(s):
    # tolerate both "...+00:00" and trailing "Z"
    return datetime.fromisoformat(s.replace('Z', '+00:00'))


def _is_stale(iso_str, max_age_days):
    try:
        ts = _parse_iso(iso_str)
    except (ValueError, AttributeError):
        return True
    return datetime.now(timezone.utc) - ts > timedelta(days=max_age_days)


def _load_manifest(cache_dir):
    path = os.path.join(cache_dir, _MANIFEST_FILENAME)
    if not os.path.exists(path):
        return {}
    try:
        with open(path, 'r', encoding='utf-8') as fh:
            return json.load(fh)
    except (json.JSONDecodeError, OSError):
        return {}


def _save_manifest(cache_dir, manifest):
    path = os.path.join(cache_dir, _MANIFEST_FILENAME)
    with open(path, 'w', encoding='utf-8') as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
        fh.write('\n')


def _resolve_refresh_settings(max_age_days, refresh_mode):
    if max_age_days is None:
        max_age_days = int(os.environ.get('COOLPROP_PUBCHEM_MAX_AGE_DAYS',
                                          _DEFAULT_MAX_AGE_DAYS))
    if refresh_mode is None:
        refresh_mode = os.environ.get('COOLPROP_PUBCHEM_REFRESH', '').lower()
    if refresh_mode not in ('', 'all', 'stale'):
        print(f'  warning: ignoring COOLPROP_PUBCHEM_REFRESH={refresh_mode!r} '
              f'(expected "all", "stale", or unset)')
        refresh_mode = ''
    CACHE_FRESHNESS['max_age_days'] = max_age_days
    CACHE_FRESHNESS['refresh_mode'] = refresh_mode
    return max_age_days, refresh_mode


def _migrate_legacy_entry(cache_path, fail_flag):
    """Synthesize a manifest entry from a pre-manifest cache file or .failed flag."""
    if os.path.exists(cache_path):
        mtime = datetime.fromtimestamp(os.path.getmtime(cache_path), timezone.utc)
        return {'fetched': mtime.replace(microsecond=0).isoformat(), 'status': 'ok'}
    if os.path.exists(fail_flag):
        mtime = datetime.fromtimestamp(os.path.getmtime(fail_flag), timezone.utc)
        return {'fetched': mtime.replace(microsecond=0).isoformat(), 'status': 'missing'}
    return None


def fetch_pubchem_sdf(inchikey, cache_dir, max_age_days=None, refresh_mode=None):
    """
    Fetch SDF from PubChem for *inchikey*, trying 3-D conformer first then 2-D.
    Results are cached in *cache_dir* and tracked in ``manifest.json`` with
    fetch timestamps so freshness survives git checkouts.

    Cache freshness (see ``CACHE_FRESHNESS``) is enforced via:
      - ``max_age_days`` / env ``COOLPROP_PUBCHEM_MAX_AGE_DAYS`` (default 180)
      - ``refresh_mode`` / env ``COOLPROP_PUBCHEM_REFRESH``: ``""`` (use cache,
        warn on stale), ``"stale"`` (refetch entries older than max_age_days),
        ``"all"`` (refetch every entry).

    Returns ``(sdf_string, is_3d)`` or ``(None, None)`` when unavailable.
    """
    max_age_days, refresh_mode = _resolve_refresh_settings(max_age_days, refresh_mode)
    os.makedirs(cache_dir, exist_ok=True)
    manifest = _load_manifest(cache_dir)

    try:
        for record_type, is_3d in (('3d', True), ('2d', False)):
            key = f'{inchikey}_{record_type}'
            cache_path = os.path.join(cache_dir, f'{key}.sdf')
            fail_flag = cache_path + '.failed'
            entry = manifest.get(key)
            if entry is None:
                entry = _migrate_legacy_entry(cache_path, fail_flag)
                if entry is not None:
                    manifest[key] = entry

            stale = entry is not None and _is_stale(entry['fetched'], max_age_days)
            force_refetch = (refresh_mode == 'all') or (refresh_mode == 'stale' and stale)

            if entry is not None and not force_refetch:
                if stale and key not in CACHE_FRESHNESS['stale']:
                    CACHE_FRESHNESS['stale'].append(key)
                if entry['status'] == 'ok' and os.path.exists(cache_path):
                    with open(cache_path, 'r', encoding='utf-8') as fh:
                        return fh.read(), is_3d
                if entry['status'] == 'missing':
                    continue
                # status=ok but file missing → fall through and refetch

            url = (f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/'
                   f'{inchikey}/SDF?record_type={record_type}')
            session = requests.Session()
            retry = Retry(total=3, backoff_factor=1,
                          status_forcelist=[429, 500, 502, 503, 504],
                          allowed_methods=['GET'])
            session.mount('https://', HTTPAdapter(max_retries=retry))
            try:
                resp = session.get(url, timeout=15)
                if resp.status_code == 404:
                    open(fail_flag, 'w').close()
                    manifest[key] = {'fetched': _now_iso(), 'status': 'missing'}
                    continue
                resp.raise_for_status()
                content = resp.text
                with open(cache_path, 'w', encoding='utf-8') as fh:
                    fh.write(content)
                if os.path.exists(fail_flag):
                    os.remove(fail_flag)
                manifest[key] = {'fetched': _now_iso(), 'status': 'ok'}
                action = 'refreshed' if force_refetch else 'fetched'
                print(f'  {action} {record_type.upper()} SDF for {inchikey}')
                return content, is_3d
            except Exception as exc:
                print(f'  PubChem {record_type.upper()} SDF unavailable for {inchikey}: {exc}')
                open(fail_flag, 'w').close()
                manifest[key] = {'fetched': _now_iso(), 'status': 'missing'}
        return None, None
    finally:
        _save_manifest(cache_dir, manifest)


def generate_3dmol_rst(fluid, sdf_data, is_3d):
    """
    Return an RST string containing a ``.. raw:: html`` block with an embedded
    interactive py3Dmol viewer.  The SDF data is inlined as a JS template
    literal so no external file serving is required.
    """
    # Escape SDF content for safe embedding inside a JS template literal
    sdf_js = sdf_data.replace('\\', '\\\\').replace('`', '\\`').replace('${', '\\${')
    dim_label = '3D conformer' if is_3d else '2D structure'

    html = (
        f'<div class="molecule-viewer" style="text-align:center;margin:1em 0;">\n'
        f'  <div id="mol3d_{fluid}" style="width:400px;height:320px;'
        f'position:relative;display:inline-block;border:1px solid #ccc;border-radius:6px;"></div>\n'
        f'  <p style="margin-top:0.4em;color:#666;font-size:0.85em;">\n'
        f'    {fluid} \u2014 {dim_label} (interactive: click and drag to rotate)\n'
        f'  </p>\n'
        f'</div>\n'
        f'<script>\n'
        f'(function() {{\n'
        f'  var sdf = `{sdf_js}`;\n'
        f'  function init() {{\n'
        f'    if (typeof $3Dmol === "undefined") {{ setTimeout(init, 150); return; }}\n'
        f'    var v = $3Dmol.createViewer(\n'
        f'      document.getElementById("mol3d_{fluid}"),\n'
        f'      {{backgroundColor: "white"}});\n'
        f'    v.addModel(sdf, "sdf");\n'
        f'    v.setStyle({{}}, {{stick: {{radius: 0.15}}, sphere: {{scale: 0.3}}}});\n'
        f'    v.zoomTo();\n'
        f'    v.resize();\n'
        f'    v.render();\n'
        f'  }}\n'
        f'  if (document.readyState === "loading") {{\n'
        f'    document.addEventListener("DOMContentLoaded", init);\n'
        f'  }} else {{\n'
        f'    init();\n'
        f'  }}\n'
        f'}})();\n'
        f'</script>'
    )

    indented = '\n'.join('   ' + line for line in html.split('\n'))
    return (
        'Molecular Structure\n'
        '===================\n'
        '\n'
        '.. raw:: html\n'
        '\n'
        f'{indented}\n'
        '\n'
    )


def _freshness_section(cache_dir):
    """Return Markdown lines summarising PubChem cache freshness, plus a
    boolean indicating whether any stale entries were found."""
    manifest = _load_manifest(cache_dir)
    max_age = CACHE_FRESHNESS['max_age_days']
    stale = []
    for key, entry in manifest.items():
        if _is_stale(entry.get('fetched', ''), max_age):
            stale.append((key, entry.get('fetched', '?'), entry.get('status', '?')))
    stale.sort(key=lambda t: t[1])  # oldest first

    lines = ['## PubChem cache freshness', '']
    if not manifest:
        lines += ['_No `manifest.json` found — cache freshness cannot be checked._', '']
        return lines, False
    lines += [
        f'- **Cache entries tracked:** {len(manifest)}',
        f'- **Max age before refresh recommended:** {max_age} days',
        f'- **Stale entries:** {len(stale)}',
        '',
    ]
    if stale:
        lines += [
            'To refresh stale entries, rebuild with '
            '`COOLPROP_PUBCHEM_REFRESH=stale`. To refresh everything, use '
            '`COOLPROP_PUBCHEM_REFRESH=all`. Override the threshold with '
            '`COOLPROP_PUBCHEM_MAX_AGE_DAYS=<n>`.',
            '',
            '| Key | Fetched (UTC) | Status |',
            '| --- | --- | --- |',
        ]
        lines += [f'| `{k}` | {ts} | {st} |' for k, ts, st in stale]
        lines.append('')
    return lines, bool(stale)


def write_geometry_status_report(out_path, cache_dir=None):
    """
    Write a human-friendly Markdown summary of per-fluid 3D molecule-viewer
    status to *out_path*.  Buckets are populated by FluidGenerator.write() as
    it runs, so call this *after* the per-fluid generation loop.

    If *cache_dir* is given (the directory containing the SDF cache and
    ``manifest.json``), the report also includes a PubChem cache freshness
    section.
    """
    s = GEOMETRY_STATUS
    total = sum(len(s[k]) for k in ('has_3d', 'has_2d_only', 'pubchem_missing', 'no_inchikey'))
    lines = [
        '# Molecule geometry status',
        '',
        'Generated by `Web/scripts/fluid_properties.PurePseudoPure.py` from '
        '`fetch_pubchem_sdf()` outcomes. Buckets are mutually exclusive and '
        'cover every fluid in `CoolProp.__fluids__`.',
        '',
        f'- **Total fluids:** {total}',
        f'- **3D conformer:** {len(s["has_3d"])}',
        f'- **2D fallback only:** {len(s["has_2d_only"])}',
        f'- **No PubChem record:** {len(s["pubchem_missing"])}',
        f'- **No valid InChIKey** (pseudo-pures, ortho/para variants): {len(s["no_inchikey"])}',
        '',
        '## 2D fallback only',
        '',
        'PubChem returned a 2D structure but no 3D conformer.',
        '',
    ]
    if s['has_2d_only']:
        lines += [f'- {name}' for name in sorted(s['has_2d_only'])]
    else:
        lines.append('_(none)_')
    lines += ['', '## No PubChem record', '',
              'PubChem returned 404 for both 3D and 2D requests.', '']
    if s['pubchem_missing']:
        lines += [f'- {name} — `{key}`' for name, key in sorted(s['pubchem_missing'])]
    else:
        lines.append('_(none)_')
    lines += ['', '## No valid InChIKey', '',
              'These fluids have `INCHIKEY = "N/A"` (or empty) and are skipped '
              'before any PubChem lookup. Typical cases: pseudo-pure mixtures '
              '(R407C, R410A, SES36, …) and ortho/para hydrogen variants.', '']
    if s['no_inchikey']:
        lines += [f'- {name}' for name in sorted(s['no_inchikey'])]
    else:
        lines.append('_(none)_')
    lines += ['', '## 3D conformer (success)', '']
    if s['has_3d']:
        lines += [f'- {name}' for name in sorted(s['has_3d'])]
    else:
        lines.append('_(none)_')
    lines.append('')

    has_stale = False
    if cache_dir is not None:
        freshness_lines, has_stale = _freshness_section(cache_dir)
        lines += freshness_lines

    with codecs.open(out_path, 'w', encoding='utf-8') as fp:
        fp.write('\n'.join(lines))
    print(f'wrote molecule geometry status report: {out_path}')
    if has_stale:
        n = len([1 for entry in _load_manifest(cache_dir).values()
                 if _is_stale(entry.get('fetched', ''), CACHE_FRESHNESS['max_age_days'])])
        print(f'  WARNING: {n} PubChem cache entries are older than '
              f'{CACHE_FRESHNESS["max_age_days"]} days — see report for details. '
              f'Refresh with COOLPROP_PUBCHEM_REFRESH=stale.')
    return has_stale


def formula2RST(formula):
    """
    See: https://docutils.sourceforge.io/docs/ref/rst/roles.html#subscript
    """
    return formula.replace('_{', r'\ :sub:`').replace('}',r'`\ ').replace(r'\ :sub:`1`\ ', '')


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
                l, r = n.split('e')
                n = rf' :math:`{l}@\times 10^{{{r}}}`'
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
            formula = formula2RST(formula)
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

        # Generate interactive 3-D molecule viewer (py3Dmol via PubChem SDF)
        molecule_viewer_rst = ''
        inchikey = CoolProp.CoolProp.get_fluid_param_string(self.fluid, "INCHIKEY")
        if inchikey and _INCHIKEY_RE.match(inchikey):
            sdf_cache = os.path.join(path, 'molecule_sdf')
            sdf_data, is_3d = fetch_pubchem_sdf(inchikey, sdf_cache)
            if sdf_data:
                molecule_viewer_rst = generate_3dmol_rst(self.fluid, sdf_data, is_3d)
                bucket = 'has_3d' if is_3d else 'has_2d_only'
                GEOMETRY_STATUS[bucket].append(self.fluid)
            else:
                GEOMETRY_STATUS['pubchem_missing'].append((self.fluid, inchikey))
        else:
            GEOMETRY_STATUS['no_inchikey'].append(self.fluid)

        # Write RST file for fluid
        out = fluid_template.format(aliases=aliases,
                                    fluid=self.fluid,
                                    fluid_stars='*' * len(self.fluid),
                                    references=references,
                                    molecule_viewer_rst=molecule_viewer_rst
                                    )

        with codecs.open(os.path.join(path, self.fluid + '.rst'), 'w', encoding='utf-8') as fp:
            print("writing %s" % os.path.join(path, self.fluid + '.rst'))
            fp.write(out)
