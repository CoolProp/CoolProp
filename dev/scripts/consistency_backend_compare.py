#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Run the flash-consistency grid on two backends and compare them fluid by fluid.

The docs build renders :class:`CoolProp.Plots.ConsistencyPlots.ConsistencyFigure`
once per fluid on a single backend.  This script runs the *same* grid on two
backends over the fluids both of them carry, and emits

* ``points/<fluid>__<backend>.csv``  -- every failing state point (the per-fluid
  CSV the docs build publishes, one per backend);
* ``backend_compare_points.csv``     -- all of the above concatenated;
* ``backend_compare_pairs.csv``      -- one row per (fluid, input pair): failure
  counts and mean flash time for each backend, plus the time ratio;
* ``backend_compare.json``           -- the same data plus per-fluid totals and
  the deduplicated exception messages, for a report generator.

Usage::

    COOLPROP_REFPROP_ROOT=~/REFPROP10 \\
        python dev/scripts/consistency_backend_compare.py --out /tmp/cmp

    # or, if REFPROP lives somewhere CoolProp cannot guess:
    python dev/scripts/consistency_backend_compare.py --out /tmp/cmp \\
        --refprop-path ~/REFPROP10/

Each fluid is rendered in its own subprocess (``--one-fluid``), so a backend that
segfaults or leaks costs one fluid rather than the run.  Both backends are timed
back to back *inside the same worker*, so the per-pair time ratio is meaningful
even when several workers compete for cores; the absolute times are not.
"""
from __future__ import print_function, division, absolute_import

import argparse
import concurrent.futures
import json
import os
import signal
import subprocess
import sys
import time

import pandas

# The grid is unchanged from the docs build; these are the panels it evaluates.
# (ConsistencyFigure crosses out `not_implemented_solvers` itself.)
CLASSES = ['INCONSISTENT', 'EXCEPTION', 'BAD_PHASE']


def _apply_refprop_path(path):
    """Point CoolProp at a REFPROP install for this process and its children."""
    if not path:
        return
    import CoolProp.CoolProp as CP
    path = os.path.expanduser(path)
    # ALTERNATIVE_REFPROP_PATH wants a trailing separator; COOLPROP_REFPROP_ROOT is
    # what the worker subprocesses inherit, so set both from the one flag.
    CP.set_config_string(CP.ALTERNATIVE_REFPROP_PATH, path.rstrip(os.sep) + os.sep)
    os.environ['COOLPROP_REFPROP_ROOT'] = path.rstrip(os.sep)


def loadable(backend, fluid):
    """True when ``backend`` can actually instantiate ``fluid``.

    Constructing the state is not always enough to force the fluid files to load,
    so a fluid constant is read as well.
    """
    import CoolProp.CoolProp as CP
    try:
        AS = CP.AbstractState(backend, fluid)
        AS.keyed_output(CP.iT_critical)
        return True
    except (KeyboardInterrupt, SystemExit):
        raise  # an interrupt must abort the run, not mark this fluid unavailable
    except BaseException:  # noqa: BLE001 -- a backend may raise anything here
        return False


def run_one_backend(fluid, backend, out_dir):
    """Render one (fluid, backend) grid.  Returns the per-backend result dict."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from CoolProp.Plots.ConsistencyPlots import ConsistencyFigure, not_implemented_solvers
    from CoolProp.Plots import _consistency_report as rpt

    tic = time.time()
    ff = ConsistencyFigure(fluid, backend=backend)
    wall = time.time() - tic

    errors = ff.errors
    failures = errors[errors['cls'].isin(CLASSES)] if len(errors) and 'cls' in errors else errors.head(0)

    points_dir = os.path.join(out_dir, 'points')
    os.makedirs(points_dir, exist_ok=True)
    csv_path = os.path.join(points_dir, '{0}__{1}.csv'.format(fluid, backend))
    rpt.write_csv(errors, csv_path)

    # Per-pair counts and timings.  The timings come off the axes rather than the
    # errors frame: ConsistencyFigure drops the GOOD rows (and with them the timing
    # of every flash that actually solved) before handing back `errors`.
    pairs = {}
    for ax, pair in zip(ff.axes_list, ff.pairs):
        if pair in not_implemented_solvers:
            continue  # panel is crossed out for every backend
        # NB: do NOT test the timings for this instead.  A panel that early-returned
        # from setup_failure_frame has no timing but does have a failure row, and
        # skipping it here would drop that row from the per-pair table while leaving
        # it in the totals.
        rows = failures[failures['pair'] == pair] if len(failures) else failures
        counts = {cls.lower(): int((rows['cls'] == cls).sum()) if len(rows) else 0 for cls in CLASSES}
        # Split the exceptions: 'reference' rows are grid points the backend would not
        # define at all (the p,T / Q,T setup flash threw), so they never tested the pair
        # flash.  Counting them as flash failures makes a strict backend look broken.
        reference = 0
        if len(rows) and 'type' in rows:
            reference = int(((rows['cls'] == 'EXCEPTION') & (rows['type'] == 'reference')).sum())
        # 'setup' rows (a panel that could not be built) count with the pair-flash
        # failures, not with the domain refusals: something went wrong on our side.
        pairs[pair] = dict(
            inconsistent=counts['inconsistent'],
            exceptions=counts['exception'],
            reference=reference,
            bad_phase=counts['bad_phase'],
            # GOOD-only means: a backend that throws early on many points would look
            # fast under the all-points mean that the plot annotation shows.
            t1=ax.mean_elapsed_1phase_good,
            t2=ax.mean_elapsed_2phase_good,
            t1_all=ax.mean_elapsed_1phase,
            t2_all=ax.mean_elapsed_2phase,
        )

    # Deduplicated exception texts, most frequent first, for the report.
    # Deduplicated failure texts for the report, most frequent first.  BAD_PHASE is
    # included alongside EXCEPTION: its message names the two phase labels that
    # disagreed, and one systematic mislabelling can account for tens of thousands of
    # points, which is worth seeing as one line rather than one count.
    messages = []
    if len(failures) and 'err' in failures:
        named = failures[failures['cls'].isin(['EXCEPTION', 'BAD_PHASE'])]
        if len(named):
            # A frame with no 'type' column at all (a fluid whose only failures are
            # BAD_PHASE rows, which set none) must still group: assigning the bare
            # string here is what makes that arm work -- .fillna on it would raise,
            # and the raise would surface as "this backend failed" for the fluid.
            kind = named['type'].fillna('update') if 'type' in named else 'update'
            grouped = named.assign(_kind=kind).groupby(['cls', 'err', '_kind']).size()
            messages = [dict(cls=str(k[0]), err=str(k[1]), kind=str(k[2]), count=int(v))
                        for k, v in grouped.items()]
            messages.sort(key=lambda m: -m['count'])

    plt.close(ff.fig)
    del ff
    return dict(ok=True, wall=wall, pairs=pairs, messages=messages,
                totals=dict(
                    {cls.lower(): int((failures['cls'] == cls).sum()) if len(failures) else 0
                     for cls in CLASSES},
                    reference=int(((failures['cls'] == 'EXCEPTION') & (failures['type'] == 'reference')).sum())
                    if (len(failures) and 'type' in failures) else 0),
                csv=os.path.relpath(csv_path, out_dir))


def one_fluid(fluid, backends, out_dir):
    """Render every backend for one fluid; a backend that dies is recorded, not raised."""
    result = dict(fluid=fluid, backends={})
    for backend in backends:
        try:
            result['backends'][backend] = run_one_backend(fluid, backend, out_dir)
        except (KeyboardInterrupt, SystemExit):
            raise  # an interrupt must abort the run, not be recorded as a backend failure
        except BaseException as exc:  # noqa: BLE001 -- one backend must not lose the fluid
            result['backends'][backend] = dict(ok=False, error='{0}: {1}'.format(type(exc).__name__, exc),
                                               wall=None, pairs={}, messages=[], totals={})
    return result


def describe_exit(returncode):
    """``(crashed, human text)`` for a worker's exit status.

    A negative return code means the process was killed by a signal -- a segfault
    or abort inside the native library, not a Python-level failure.  That is a
    different animal from a flash that raised: nothing was measured, and the cause
    is a bug in the library rather than a property of a state point.  It gets its
    own label so it can never be read as ordinary noise.
    """
    if returncode < 0:
        sig = -returncode
        try:
            name = signal.Signals(sig).name
        except ValueError:
            name = 'signal %d' % sig
        return True, 'CRASHED: killed by signal %d (%s)' % (sig, name)
    if returncode > 0:
        return False, 'worker exited %d' % returncode
    return False, 'worker exited 0 without writing a result'


def _worker(fluid, backends, out_dir, refprop_path):
    """Parent-side half of a worker: run one fluid in a fresh interpreter."""
    # exist_ok: several worker threads reach this concurrently.
    summary_dir = os.path.join(out_dir, 'summary')
    os.makedirs(summary_dir, exist_ok=True)
    summary_path = os.path.join(summary_dir, fluid + '.json')
    # Drop any summary left by an earlier run: the exit check below treats a present
    # summary as proof the worker got that far, and a stale one would let a crashed
    # worker report the previous run's numbers as this run's.
    if os.path.exists(summary_path):
        os.remove(summary_path)
    cmd = [sys.executable, '-u', os.path.abspath(__file__),
           '--one-fluid', fluid, '--out', out_dir, '--backends', ','.join(backends)]
    if refprop_path:
        cmd += ['--refprop-path', refprop_path]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          universal_newlines=True)
    if proc.returncode != 0 or not os.path.exists(summary_path):
        # The subprocess died before it could write its summary.  Record that as the
        # fluid's result rather than dropping it: a missing fluid in the report would
        # read as "nothing to see here".
        crashed, what = describe_exit(proc.returncode)
        tail = (proc.stdout or '').strip().splitlines()[-6:]
        return dict(fluid=fluid, crashed=crashed,
                    backends={b: dict(ok=False, crashed=crashed, wall=None, pairs={}, messages=[],
                                      totals={}, error='%s%s' % (what, (': ' + ' | '.join(tail)) if tail else ''))
                              for b in backends})
    with open(summary_path) as fp:
        return json.load(fp)


def aggregate(results, backends, out_dir, unavailable):
    """Write the concatenated CSVs and the JSON aggregate."""
    pair_rows = []
    for res in results:
        fluid = res['fluid']
        pairs = sorted({p for b in backends for p in res['backends'].get(b, {}).get('pairs', {})})
        for pair in pairs:
            row = dict(fluid=fluid, pair=pair)
            for b in backends:
                d = res['backends'].get(b, {}).get('pairs', {}).get(pair)
                prefix = b + '_'
                if d is None:
                    for k in ('inconsistent', 'exceptions', 'reference', 'bad_phase', 't1', 't2'):
                        row[prefix + k] = None
                    continue
                for k in ('inconsistent', 'exceptions', 'reference', 'bad_phase', 't1', 't2'):
                    row[prefix + k] = d[k]
            # Ratio of the *second* backend to the first, per phase region.
            a, b2 = backends[0], backends[1] if len(backends) > 1 else backends[0]
            for region in ('t1', 't2'):
                num, den = row.get(b2 + '_' + region), row.get(a + '_' + region)
                row['ratio_' + region] = (num / den) if (num and den) else None
            pair_rows.append(row)

    pairs_df = pandas.DataFrame(pair_rows)
    pairs_df.to_csv(os.path.join(out_dir, 'backend_compare_pairs.csv'), index=False)

    frames = []
    for res in results:
        for b in backends:
            rel = res['backends'].get(b, {}).get('csv')
            if not rel:
                continue
            path = os.path.join(out_dir, rel)
            if os.path.exists(path):
                df = pandas.read_csv(path)
                if len(df):
                    frames.append(df)
    # Always write, even with nothing to write: leaving a previous run's file in
    # place would present stale points as this run's result.
    points_path = os.path.join(out_dir, 'backend_compare_points.csv')
    if frames:
        pandas.concat(frames, ignore_index=True).to_csv(points_path, index=False)
    else:
        from CoolProp.Plots._consistency_report import _CSV_COLS
        pandas.DataFrame(columns=_CSV_COLS).to_csv(points_path, index=False)

    payload = dict(backends=backends, fluids=results, unavailable=unavailable,
                   generated=time.strftime('%Y-%m-%d %H:%M:%S'))
    with open(os.path.join(out_dir, 'backend_compare.json'), 'w') as fp:
        json.dump(payload, fp)
    return pairs_df


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--out', required=True, help='output directory')
    parser.add_argument('--backends', default='HEOS,REFPROP',
                        help='comma-separated backends, reference first (default: HEOS,REFPROP)')
    parser.add_argument('--fluids', default='', help='comma-separated subset (default: all common fluids)')
    parser.add_argument('--jobs', type=int, default=0, help='worker processes (default: core count)')
    parser.add_argument('--refprop-path', default='', help='REFPROP install directory')
    parser.add_argument('--one-fluid', default='', help=argparse.SUPPRESS)
    args = parser.parse_args(argv)

    backends = [b for b in args.backends.split(',') if b]
    out_dir = os.path.abspath(os.path.expanduser(args.out))
    os.makedirs(out_dir, exist_ok=True)
    _apply_refprop_path(args.refprop_path)

    if args.one_fluid:
        result = one_fluid(args.one_fluid, backends, out_dir)
        summary_dir = os.path.join(out_dir, 'summary')
        os.makedirs(summary_dir, exist_ok=True)
        with open(os.path.join(summary_dir, args.one_fluid + '.json'), 'w') as fp:
            json.dump(result, fp)
        return 0

    import CoolProp
    candidates = [f for f in args.fluids.split(',') if f] or list(CoolProp.__fluids__)
    fluids, unavailable = [], []
    for fluid in candidates:
        missing = [b for b in backends if not loadable(b, fluid)]
        if missing:
            unavailable.append(dict(fluid=fluid, missing_from=missing))
        else:
            fluids.append(fluid)
    print('%d of %d fluids are available in all of %s' % (len(fluids), len(candidates), ', '.join(backends)))
    for entry in unavailable:
        print('  skipping %-24s (not in %s)' % (entry['fluid'], ', '.join(entry['missing_from'])))
    if not fluids:
        sys.exit('No fluid is available in every requested backend.')

    jobs = args.jobs or (os.cpu_count() or 1)
    jobs = max(1, min(jobs, len(fluids)))
    print('Rendering %d fluid(s) x %d backend(s) with %d worker(s)' % (len(fluids), len(backends), jobs))

    results = []
    tic = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as ex:
        futures = {ex.submit(_worker, f, backends, out_dir, args.refprop_path): f for f in fluids}
        for n, fut in enumerate(concurrent.futures.as_completed(futures), 1):
            res = fut.result()
            results.append(res)
            dead = [b for b in backends if not res['backends'].get(b, {}).get('ok')]
            mark = 'ok'
            if res.get('crashed'):
                mark = '*** CRASHED: ' + ','.join(dead)
            elif dead:
                mark = 'FAILED: ' + ','.join(dead)
            print('[%3d/%3d] %-24s %s' % (n, len(fluids), res['fluid'], mark))
    results.sort(key=lambda r: r['fluid'])
    pairs_df = aggregate(results, backends, out_dir, unavailable)
    print('Done in %.1f s; %d pair rows -> %s' % (time.time() - tic, len(pairs_df), out_dir))

    # A worker killed by a signal is a native-code bug, and one line among a hundred
    # and thirty scrolls away.  Repeat them all at the end, and exit non-zero, so a
    # crash cannot be mistaken for a clean run by a human or by a shell.
    crashed = [(r['fluid'], b, r['backends'][b].get('error', ''))
               for r in results for b in backends
               if r['backends'].get(b, {}).get('crashed')]
    failed = [(r['fluid'], b, r['backends'][b].get('error', ''))
              for r in results for b in backends
              if not r['backends'].get(b, {}).get('ok') and not r['backends'][b].get('crashed')]
    if failed:
        print('\n%d (fluid, backend) run(s) produced no result:' % len(failed))
        for fluid, backend, err in failed:
            print('  %-24s %-10s %s' % (fluid, backend, str(err)[:120]))
    if crashed:
        bar = '=' * 78
        print('\n' + bar)
        print(' %d WORKER CRASH(ES) -- a backend died in native code; nothing was measured' % len(crashed))
        print(bar)
        for fluid, backend, err in crashed:
            print('  %-24s %-10s %s' % (fluid, backend, str(err)[:120]))
        print('\n Reproduce one on its own to get the stack:')
        # Absolute paths: relpath against whatever cwd the run happened to have
        # produces a line of ../../.. that nobody can paste.
        print('   %s %s --out %s --fluids %s%s\n'
              % (sys.executable, os.path.abspath(__file__), out_dir, crashed[0][0],
                 (' --refprop-path ' + args.refprop_path) if args.refprop_path else ''))
        print(bar)
        return 2
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
