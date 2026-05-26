from __future__ import print_function
import os.path
import CoolProp
import pandas
import subprocess
import sys
from CoolProp.Plots import _consistency_report as rpt

web_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
root_dir = os.path.abspath(os.path.join(web_dir, '..'))
fluids_path = os.path.join(web_dir, 'fluid_properties', 'fluids')

backend = os.environ.get('COOLPROP_CONSISTENCY_BACKEND', 'HEOS')
subdir = 'Consistencyplots' if backend == 'HEOS' else 'Consistencyplots_' + backend
plots_path = os.path.join(fluids_path, subdir)
# REFPROP stub dir is always present so the fluid-page include never breaks.
refprop_subdir = 'Consistencyplots_REFPROP'
refprop_path = os.path.join(fluids_path, refprop_subdir)

template = """from __future__ import division, print_function
import matplotlib
matplotlib.use('Agg')  # Force mpl to use a non-GUI backend

import matplotlib.pyplot as plt
from CoolProp.Plots.ConsistencyPlots import ConsistencyFigure
from CoolProp.Plots import _consistency_report as rpt

backend = {backend!r}
fluid = {fluid!r}
csv_relpath = {csv_relpath!r}
intro_rst = {intro_rst!r}

ff = ConsistencyFigure(fluid, backend=backend)
ff.savefig(fluid + '.png', dpi=30)
ff.savefig(fluid + '.pdf')
plt.close()

rpt.write_csv(ff.errors, fluid + '-consistency.csv')
summary = rpt.summarize(ff.errors)
summary.insert(0, 'fluid', fluid)
summary.to_csv(fluid + '-summary.csv', index=False)
rpt.write_rst_fragment(ff.errors, fluid, backend, csv_relpath,
                       fluid + '-report.rst', intro_rst=intro_rst)
del ff
"""

force = os.environ.get('COOLPROP_FORCE_CONSISTENCY', '').lower() in ('1', 'true', 'yes')

if not os.path.exists(plots_path):
    os.makedirs(plots_path)
if not os.path.exists(refprop_path):
    os.makedirs(refprop_path)


def refprop_intro(fluid):
    """Image + downloads RST that the fluid page does NOT carry for REFPROP."""
    return (
        '.. image:: {sub}/{fluid}.png\n\n'
        ':download:`REFPROP consistency plot (PDF) <{sub}/{fluid}.pdf>`\n'
        .format(sub=refprop_subdir, fluid=fluid))


def regenerate_fragment_from_csv(fluid):
    """Cheap rebuild of the RST fragment + summary from a cached failures CSV."""
    csv_path = os.path.join(plots_path, fluid + '-consistency.csv')
    errors = pandas.read_csv(csv_path)
    csv_relpath = subdir + '/' + fluid + '-consistency.csv'
    intro_rst = refprop_intro(fluid) if backend == 'REFPROP' else ''
    rpt.write_rst_fragment(errors, fluid, backend, csv_relpath,
                           os.path.join(plots_path, fluid + '-report.rst'), intro_rst=intro_rst)
    summary = rpt.summarize(errors)
    summary.insert(0, 'fluid', fluid)
    summary.to_csv(os.path.join(plots_path, fluid + '-summary.csv'), index=False)


build_failures = []

for fluid in CoolProp.__fluids__:
    png_path = os.path.join(plots_path, fluid + '.png')
    csv_path = os.path.join(plots_path, fluid + '-consistency.csv')

    if os.path.exists(png_path) and not force and os.path.exists(csv_path):
        print('fluid:', fluid, '- cached; regenerating fragment from CSV')
        regenerate_fragment_from_csv(fluid)
    else:
        print('fluid:', fluid, '- generating (backend=%s)' % backend)
        csv_relpath = subdir + '/' + fluid + '-consistency.csv'
        intro_rst = refprop_intro(fluid) if backend == 'REFPROP' else ''
        file_string = template.format(backend=backend, fluid=fluid,
                                      csv_relpath=csv_relpath, intro_rst=intro_rst)
        file_path = os.path.join(plots_path, fluid + '.py')
        with open(file_path, 'w') as fp:
            fp.write(file_string)
        try:
            subprocess.check_call('python -u "' + fluid + '.py"', cwd=plots_path,
                                  stdout=sys.stdout, stderr=sys.stderr, shell=True)
        except subprocess.CalledProcessError as exc:
            print('BUILD FAILED for', fluid, ':', exc)
            build_failures.append((fluid, str(exc)))
            continue

    # During the default (HEOS) build, ensure a REFPROP stub fragment exists so the
    # fluid page's REFPROP include never points at a missing file.
    if backend == 'HEOS':
        refprop_frag = os.path.join(refprop_path, fluid + '-report.rst')
        if not os.path.exists(refprop_frag):
            rpt.write_stub_fragment(refprop_frag,
                                    'REFPROP consistency plots were not generated in this build.')

# Aggregate per-fluid summaries into the consolidated report page.
summaries = []
for fluid in CoolProp.__fluids__:
    spath = os.path.join(plots_path, fluid + '-summary.csv')
    if os.path.exists(spath):
        s = pandas.read_csv(spath)
        if len(s):
            summaries.append(s)
combined = pandas.concat(summaries, ignore_index=True) if summaries else pandas.DataFrame()

report_name = 'ConsistencyReport.rst' if backend == 'HEOS' else 'ConsistencyReport_%s.rst' % backend
report_path = os.path.join(web_dir, 'fluid_properties', report_name)
rpt.write_consolidated_rst(combined, report_path, backend,
                           build_failures=build_failures, orphan=(backend != 'HEOS'))
print('Wrote consolidated report:', report_path)
