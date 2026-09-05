"""Build the per-fluid flash-consistency plots and reports for the documentation.

Environment:

  COOLPROP_CONSISTENCY_BACKEND  backend to evaluate (default HEOS).  Anything other
                                than HEOS writes to fluids/Consistencyplots_<backend>/
                                and to an orphan ConsistencyReport_<backend>.rst, and
                                has its fluid list filtered to what the backend
                                actually carries.
  COOLPROP_CONSISTENCY_FLUIDS   comma-separated subset (default: every fluid).
  COOLPROP_CONSISTENCY_JOBS     worker count (default: core count).
  COOLPROP_CONSISTENCY_DPI      PNG resolution (default 100).
  COOLPROP_FORCE_CONSISTENCY    rebuild even when a cached PNG + CSV exist.
  COOLPROP_REFPROP_ROOT         REFPROP install directory; inherited by the per-fluid
                                worker subprocesses, so it is all that a REFPROP run
                                needs.  The documentation container has no REFPROP,
                                which is why the REFPROP pass is opt-in rather than
                                part of the default build.

  COOLPROP_REFPROP_ROOT=~/REFPROP10 COOLPROP_CONSISTENCY_BACKEND=REFPROP \
      python fluid_properties.Consistency.py
"""
from __future__ import print_function
import os.path
import concurrent.futures
import signal
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
                if proc.returncode < 0:
                    # Killed by a signal: the library died in native code rather than
                    # raising.  That is a different class from a flash that threw, and
                    # it is the one a reader must not skim past, so it is labelled and
                    # banner-printed.  The build still continues -- one fluid must not
                    # cost the documentation -- and the fluid is listed in the report.
                    sig = -proc.returncode
                    try:
                        signame = signal.Signals(sig).name
                    except ValueError:
                        signame = 'signal %d' % sig
                    failure = 'CRASHED: killed by signal %d (%s)' % (sig, signame)
                    print('!' * 78)
                    print('!! %s CRASHED the %s backend in native code (%s)' % (fluid, backend, signame))
                    print('!! Nothing was measured for it.  Reproduce with:')
                    print('!!   COOLPROP_CONSISTENCY_BACKEND=%s COOLPROP_CONSISTENCY_FLUIDS=%s '
                          'COOLPROP_FORCE_CONSISTENCY=1 python fluid_properties.Consistency.py'
                          % (backend, fluid))
                    print('!' * 78)
                else:
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


def split_by_availability(backend, fluids):
    """``(usable, unavailable)`` for ``backend``, where ``unavailable`` is a list of
    ``(fluid, reason)``.

    ``CoolProp.__fluids__`` is the HEOS fluid list, so it is authoritative only for
    the HEOS backend.  REFPROP ships a different (largely overlapping) set: about
    eight CoolProp fluids have no REFPROP .FLD, and rendering those would spend a
    subprocess each to produce an identical "could not load" failure and then bury
    the real findings under them in the consolidated report.  Probe once, up front,
    and report the gap as a gap rather than as failures."""
    if backend == 'HEOS':
        return list(fluids), []
    usable, unavailable = [], []
    for fluid in fluids:
        try:
            AS = CoolProp.AbstractState(backend, fluid)
            # Touch a fluid constant: constructing the state is not always enough to
            # force the backend to actually load the fluid files.
            AS.keyed_output(CoolProp.iT_critical)
            usable.append(fluid)
        except (KeyboardInterrupt, SystemExit):
            raise  # an interrupt must abort the build, not mark this fluid unavailable
        except BaseException as exc:  # noqa: BLE001 -- a backend may raise anything here
            # `or ['']`: str(exc) can be empty, and ''.splitlines() is [], so indexing
            # it would raise IndexError from inside this handler at module scope and
            # take the whole build down over one blank message.
            unavailable.append((fluid, (str(exc).splitlines() or [''])[0]))
    return usable, unavailable


# Render fluids concurrently: each is an isolated, single-threaded subprocess, so
# a worker per core saturates the (4-core) CI builder instead of running 136
# fluids strictly serially.  Worker count is overridable via
# COOLPROP_CONSISTENCY_JOBS; default to the core count.  Each worker only ever
# touches files named for its own fluid, so there are no shared-file races.
# COOLPROP_CONSISTENCY_FLUIDS restricts the run to a comma-separated subset.  The
# full sweep is 130-odd subprocesses, so without this there is no cheap way to
# exercise a backend end to end (the REFPROP path especially, which CI cannot run).
_subset = [f.strip() for f in os.environ.get('COOLPROP_CONSISTENCY_FLUIDS', '').split(',') if f.strip()]
candidates = _subset or list(CoolProp.__fluids__)
fluids, unavailable = split_by_availability(backend, candidates)
if unavailable:
    print('%d of %d fluids are not available in the %s backend and will be skipped:'
          % (len(unavailable), len(candidates), backend))
    for fluid, reason in unavailable:
        print('  %-24s %s' % (fluid, reason))
if not fluids:
    # Every fluid failed to load: the backend itself is missing (no REFPROP on this
    # machine, say).  Stop with one clear message instead of writing a report that
    # claims every fluid is broken.  Only the (empty) output directories exist at
    # this point -- no fragment, summary or consolidated report has been written, so
    # any previously cached artifacts are left exactly as they were.
    sys.exit('No fluid could be loaded with the %s backend; is the backend installed '
             'and on the search path? (REFPROP: set COOLPROP_REFPROP_ROOT)' % backend)
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
for fluid in fluids:
    spath = os.path.join(plots_path, fluid + '-summary.csv')
    if os.path.exists(spath):
        s = pandas.read_csv(spath)
        if len(s):
            summaries.append(s)
combined = pandas.concat(summaries, ignore_index=True) if summaries else pandas.DataFrame()

report_name = 'ConsistencyReport.rst' if backend == 'HEOS' else 'ConsistencyReport_%s.rst' % backend
report_path = os.path.join(web_dir, 'fluid_properties', report_name)
if _subset:
    # The consolidated page aggregates the fluids this run built.  Writing it from a
    # deliberately partial run would replace the full report with a handful of rows
    # and read as "everything else is clean".
    print('Subset run (COOLPROP_CONSISTENCY_FLUIDS); leaving %s untouched' % report_name)
    raise SystemExit(0)
rpt.write_consolidated_rst(combined, report_path, backend,
                           build_failures=build_failures, unavailable=unavailable,
                           orphan=(backend != 'HEOS'))
print('Wrote consolidated report:', report_path)

# Repeat any crash at the very end: one line among 130 fluids scrolls out of a CI
# log, and a native-code death is the finding most likely to be missed and least
# safe to miss.
crashed = [(fluid, msg) for fluid, msg in build_failures if str(msg).startswith('CRASHED')]
if crashed:
    print('!' * 78)
    print('!! %d fluid(s) crashed the %s backend in native code:' % (len(crashed), backend))
    for fluid, msg in crashed:
        print('!!   %-24s %s' % (fluid, msg))
    print('!' * 78)
