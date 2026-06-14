from __future__ import print_function
import os.path
import concurrent.futures
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
dpi = {dpi!r}

ff = ConsistencyFigure(fluid, backend=backend)
ff.savefig(fluid + '.png', dpi=dpi)
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
# The figure is a 15x23" 5x3 grid; dpi=30 (the historical default) was unreadable.
dpi = int(os.environ.get('COOLPROP_CONSISTENCY_DPI', '100'))

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


def build_one_fluid(fluid):
    """Generate (or cheaply regenerate from cache) the consistency artifacts for
    one fluid and guarantee its include fragments exist.  Returns
    ``(fluid, failure_msg_or_None)``.

    Runs in a worker thread; the expensive rendering is an isolated subprocess,
    so concurrent calls scale across cores without GIL contention.  Any
    exception is captured into the failure message rather than raised — one
    bad fluid must never abort the whole loop (the stub fragment written in the
    ``finally`` keeps the fluid page's ``.. include::`` from dangling, exactly
    as the old serial fall-through did)."""
    failure = None
    try:
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
                                          csv_relpath=csv_relpath, intro_rst=intro_rst, dpi=dpi)
            file_path = os.path.join(plots_path, fluid + '.py')
            with open(file_path, 'w') as fp:
                fp.write(file_string)
            # Capture output and emit it as one contiguous block so concurrent
            # fluids' logs don't interleave line-by-line.
            proc = subprocess.run([sys.executable, '-u', fluid + '.py'], cwd=plots_path,
                                  stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                  universal_newlines=True)
            if proc.stdout:
                sys.stdout.write('----- %s -----\n%s\n' % (fluid, proc.stdout.rstrip('\n')))
                sys.stdout.flush()
            if proc.returncode != 0:
                failure = 'exit code %d' % proc.returncode
                print('BUILD FAILED for', fluid, ':', failure)
    except Exception as exc:  # noqa: BLE001 — a single fluid must not kill the build
        failure = str(exc)
        print('BUILD FAILED for', fluid, ':', failure)
    finally:
        if failure is not None:
            # Drop stale per-fluid data artifacts on failure. These dirs are
            # restored from cache, so a prior SUCCESSFUL run may have left a
            # <fluid>-summary.csv (which the consolidated report aggregates
            # unconditionally below) and a <fluid>-consistency.csv (whose
            # presence, with the png, marks the fluid "cached" so the next run
            # skips it). Leaving them would make this build report stale success
            # for a fluid that actually failed, and mask the failure on reruns.
            # Removing them excludes the fluid from this run's aggregation and
            # forces a clean regeneration next run; build_failures still records
            # it for the report's failure list.
            for stale in (fluid + '-summary.csv', fluid + '-consistency.csv'):
                stale_path = os.path.join(plots_path, stale)
                try:
                    if os.path.exists(stale_path):
                        os.remove(stale_path)
                except OSError as exc:
                    # Best effort: failing to delete a stale artifact must not
                    # escape this finally (it would surface via ex.map and abort
                    # every other fluid). A leftover stale file is the lesser
                    # evil and is logged.
                    print('WARN: could not remove stale %s: %s' % (stale_path, exc))

        # Guarantee the per-fluid HEOS report fragment exists. The fluid page
        # includes it unconditionally, so a crash before the subprocess wrote it
        # (or a partial cache) must not leave a dangling `.. include::`.
        report_frag = os.path.join(plots_path, fluid + '-report.rst')
        if not os.path.exists(report_frag):
            rpt.write_stub_fragment(report_frag,
                                    'Consistency data could not be generated for this fluid in this build.')

        # During the default (HEOS) build, ensure a REFPROP stub fragment exists so
        # the fluid page's REFPROP include never points at a missing file.
        if backend == 'HEOS':
            refprop_frag = os.path.join(refprop_path, fluid + '-report.rst')
            if not os.path.exists(refprop_frag):
                rpt.write_stub_fragment(refprop_frag,
                                        'REFPROP consistency plots were not generated in this build.')

    return fluid, failure


# Render fluids concurrently: each is an isolated, single-threaded subprocess, so
# a worker per core saturates the (4-core) CI builder instead of running 136
# fluids strictly serially.  Worker count is overridable via
# COOLPROP_CONSISTENCY_JOBS; default to the core count.  Each worker only ever
# touches files named for its own fluid, so there are no shared-file races.
fluids = list(CoolProp.__fluids__)
jobs_env = os.environ.get('COOLPROP_CONSISTENCY_JOBS')
try:
    max_workers = int(jobs_env) if jobs_env else (os.cpu_count() or 1)
except ValueError:
    max_workers = os.cpu_count() or 1
max_workers = max(1, min(max_workers, len(fluids)))
print('Building consistency plots for %d fluids with %d worker(s)' % (len(fluids), max_workers))

build_failures = []
with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as ex:
    # map() preserves input order, so build_failures stays deterministic (fluid
    # order) for the consolidated report.
    for fluid, failure in ex.map(build_one_fluid, fluids):
        if failure is not None:
            build_failures.append((fluid, failure))

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
