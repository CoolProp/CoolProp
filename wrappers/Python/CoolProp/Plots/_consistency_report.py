# -*- coding: utf-8 -*-
"""Reporting helpers for the consistency plots.

Pure functions over the ``errors`` DataFrame produced by
``CoolProp.Plots.ConsistencyPlots.ConsistencyFigure``. No CoolProp import so the
module is unit-testable in isolation.
"""
from __future__ import print_function, division, absolute_import
import math
import codecs
import datetime as _dt
import pandas

# Map raw class labels (as stored in the errors DataFrame) to summary columns.
_CLASS_COL = {'INCONSISTENT': 'inconsistent', 'EXCEPTION': 'exceptions', 'BAD_PHASE': 'bad_phase'}
_SUMMARY_COLS = ['pair', 'inconsistent', 'exceptions', 'bad_phase', 'total_failures']
_CSV_COLS = ['fluid', 'backend', 'pair', 'phase_region', 'cls',
             'in1', 'val1', 'in2', 'val2', 'P', 'T', 'dev', 'err']


def format_time(seconds):
    """Human-readable auto-scaled duration: s / ms / microseconds."""
    if seconds is None:
        return 'n/a'
    try:
        if not math.isfinite(seconds):
            return 'n/a'
    except TypeError:
        return 'n/a'
    if seconds >= 1.0:
        return '{0:.3g} s'.format(seconds)
    if seconds >= 1e-3:
        return '{0:.3g} ms'.format(seconds * 1e3)
    return '{0:.3g} µs'.format(seconds * 1e6)


def summarize(errors):
    """Per-input-pair failure counts. Returns a DataFrame with columns
    ``_SUMMARY_COLS``, sorted by ``total_failures`` descending."""
    if errors is None or len(errors) == 0 or 'cls' not in errors:
        return pandas.DataFrame(columns=_SUMMARY_COLS)
    df = errors[errors['cls'].isin(_CLASS_COL.keys())]
    if len(df) == 0:
        return pandas.DataFrame(columns=_SUMMARY_COLS)
    counts = df.groupby(['pair', 'cls']).size().unstack(fill_value=0)
    counts = counts.rename(columns=_CLASS_COL)
    for col in ['inconsistent', 'exceptions', 'bad_phase']:
        if col not in counts:
            counts[col] = 0
    counts['total_failures'] = counts[['inconsistent', 'exceptions', 'bad_phase']].sum(axis=1)
    counts = counts.reset_index()[_SUMMARY_COLS]
    return counts.sort_values('total_failures', ascending=False).reset_index(drop=True)


def write_csv(errors, path):
    """Write the full list of failing points to CSV (header-only if none)."""
    if errors is None or len(errors) == 0 or 'cls' not in errors:
        pandas.DataFrame(columns=_CSV_COLS).to_csv(path, index=False)
        return
    df = errors[errors['cls'].isin(_CLASS_COL.keys())]
    df.reindex(columns=_CSV_COLS).to_csv(path, index=False)


def _cell(v):
    if v is None:
        return ''
    if isinstance(v, float):
        if math.isnan(v):
            return ''
        return '{0:.6g}'.format(v)
    return str(v).replace('\n', ' ').replace('|', '\\|').strip()


def _select_sample(errors, sample_cap):
    """Up to ``sample_cap`` representative failing rows per (pair, cls).
    EXCEPTION rows are de-duplicated by message; INCONSISTENT rows are ordered
    largest-deviation first when ``dev`` is available."""
    if errors is None or len(errors) == 0 or 'cls' not in errors:
        return pandas.DataFrame()
    df = errors[errors['cls'].isin(_CLASS_COL.keys())]
    chunks = []
    for (_pair, cls), grp in df.groupby(['pair', 'cls']):
        if cls == 'EXCEPTION' and 'err' in grp:
            grp = grp.drop_duplicates(subset='err')
        if 'dev' in grp and grp['dev'].notna().any():
            grp = grp.sort_values('dev', ascending=False, na_position='last')
        chunks.append(grp.head(sample_cap))
    if not chunks:
        return df.head(0)
    return pandas.concat(chunks)


def _rst_list_table(sample, indent='   '):
    cols = ['pair', 'cls', 'phase_region', 'P', 'T', 'in1', 'val1', 'in2', 'val2', 'err']
    headers = ['Pair', 'Class', 'Region', 'P [Pa]', 'T [K]', 'In1', 'Val1', 'In2', 'Val2', 'Error']
    lines = [indent + '.. list-table::', indent + '   :header-rows: 1', '']

    def row(cells):
        block = [indent + '   * - ' + _cell(cells[0])]
        for c in cells[1:]:
            block.append(indent + '     - ' + _cell(c))
        return block

    lines += row(headers)
    for _, r in sample.iterrows():
        lines += row([r.get(c, '') for c in cols])
    return '\n'.join(lines) + '\n'


def write_rst_fragment(errors, fluid, backend, csv_relpath, out_path,
                       sample_cap=20, intro_rst=''):
    """Per-fluid RST fragment: counts line + CSV download + collapsible sample table.

    ``intro_rst`` is prepended verbatim (used by the REFPROP path to carry the
    image/download directives that the fluid page does not provide)."""
    summary = summarize(errors)
    with codecs.open(out_path, 'w', encoding='utf-8') as fp:
        if intro_rst:
            fp.write(intro_rst)
            if not intro_rst.endswith('\n'):
                fp.write('\n')
            fp.write('\n')
        if len(summary) == 0:
            fp.write('**Flash consistency ({0}):** no failures.\n'.format(backend))
            return
        tot_inc = int(summary['inconsistent'].sum())
        tot_exc = int(summary['exceptions'].sum())
        tot_bad = int(summary['bad_phase'].sum())
        npairs = int((summary['total_failures'] > 0).sum())
        fp.write('**Flash consistency ({0}):** {1} inconsistent, {2} exceptions, '
                 '{3} bad-phase across {4} input pair(s).\n\n'
                 .format(backend, tot_inc, tot_exc, tot_bad, npairs))
        fp.write(':download:`Download full failure list (CSV) <{0}>`\n\n'.format(csv_relpath))
        sample = _select_sample(errors, sample_cap)
        fp.write('.. dropdown:: Failing state points (sample, up to {0} per pair/class)\n\n'
                 .format(sample_cap))
        fp.write(_rst_list_table(sample))


def write_stub_fragment(out_path, message):
    """Write a one-line placeholder fragment so an `.. include::` never breaks."""
    with codecs.open(out_path, 'w', encoding='utf-8') as fp:
        fp.write(message.rstrip() + '\n')


def write_consolidated_rst(combined_summary, out_path, backend,
                           build_failures=None, date=None, orphan=False):
    """Cross-fluid failure summary page (one row per fluid/pair with failures)."""
    build_failures = build_failures or []
    date = date or _dt.date.today().isoformat()
    with codecs.open(out_path, 'w', encoding='utf-8') as fp:
        if orphan:
            fp.write(':orphan:\n\n')
        title = 'Consistency Failure Report ({0})'.format(backend)
        fp.write(title + '\n' + '=' * len(title) + '\n\n')
        fp.write('Generated {0}. One row per fluid / input-pair with at least one '
                 'flash-consistency failure, sorted by total failures.\n\n'.format(date))
        if combined_summary is None or len(combined_summary) == 0:
            fp.write('No failures recorded across all fluids. \U0001F389\n\n')
        else:
            cs = combined_summary[combined_summary['total_failures'] > 0].copy()
            cs = cs.sort_values('total_failures', ascending=False)
            if len(cs) == 0:
                fp.write('No failures recorded across all fluids. \U0001F389\n\n')
            else:
                fp.write('.. list-table::\n   :header-rows: 1\n\n')
                fp.write('   * - Fluid\n     - Pair\n     - Inconsistent\n'
                         '     - Exceptions\n     - Bad-phase\n     - Total\n')
                for _, r in cs.iterrows():
                    fp.write('   * - :ref:`{0} <fluid_{0}>`\n'.format(r['fluid']))
                    fp.write('     - {0}\n'.format(r['pair']))
                    fp.write('     - {0}\n'.format(int(r['inconsistent'])))
                    fp.write('     - {0}\n'.format(int(r['exceptions'])))
                    fp.write('     - {0}\n'.format(int(r['bad_phase'])))
                    fp.write('     - {0}\n'.format(int(r['total_failures'])))
                fp.write('\n')
        if build_failures:
            fp.write('Build failures\n--------------\n\n')
            for fluid, msg in build_failures:
                fp.write('* {0}: {1}\n'.format(fluid, str(msg).replace('\n', ' ')))
            fp.write('\n')
