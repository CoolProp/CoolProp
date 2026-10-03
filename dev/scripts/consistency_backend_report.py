#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Render the JSON emitted by ``consistency_backend_compare.py`` as a standalone
HTML report.

Kept separate from the comparison driver so a re-render never costs another
render of the grid::

    python dev/scripts/consistency_backend_report.py \\
        --json /tmp/cmp/backend_compare.json --out /tmp/cmp/report.html

The page is self-contained (one file, no external assets beyond a Google Fonts
link) so it can be opened locally, attached to an issue, or published as-is.
"""
from __future__ import print_function, division, absolute_import

import argparse
import html
import json
import os
import re
import sys

# Numeric literals in a REFPROP/CoolProp error message vary point to point; the
# message class does not.  Collapse them so "Density above upper limit: D =
# 15.6157 mol/L" and its 400 neighbours count as one finding.
_NUM = re.compile(r'[-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?')


# CoolProp::phases, in declaration order (include/CoolProp/DataStructures.h).
_PHASES = ['liquid', 'supercritical', 'supercritical_gas', 'supercritical_liquid',
           'critical_point', 'gas', 'twophase', 'unknown', 'not_imposed']
_BAD_PHASE = re.compile(r'^phase (\d+) instead of (\d+)$')

# The input-pair names CoolProp's REFPROP backend prepends to its error text
# ("DmolarSmolar: [DSFLSH error 207] ...").  Matched by name so that an exception
# type prefixed by err_text() is never mistaken for one.
_PAIR_TAGS = ['Dmolar', 'Hmolar', 'Smolar', 'Umolar', 'P', 'T', 'Q']
_PAIR_TAGS = sorted({a + b for a in _PAIR_TAGS for b in _PAIR_TAGS} | set(_PAIR_TAGS),
                    key=len, reverse=True)


def _phase_name(index):
    try:
        return _PHASES[int(index)]
    except (ValueError, IndexError):
        return str(index)


def message_class(err):
    """Normalise one failure string to its class: strip REFPROP's fixed-width
    padding, drop the leading input-pair tag, and mask the numbers.

    A bad-phase message is the exception: its two numbers are phase enum indices,
    the whole content of the message, so they are named instead of masked --
    masking them would collapse every distinct mislabelling into one line.
    """
    text = ' '.join(str(err).split())
    # Strip the input-pair tag CoolProp's REFPROP backend prepends ("DmolarSmolar: ..."),
    # but NOT the exception type that err_text() deliberately prepends for a
    # non-ValueError: that prefix is the whole point of err_text, and erasing it here
    # would file a library defect under the same class as an ordinary out-of-range
    # report.
    # The negative lookahead covers ANY capitalised identifier, not just names ending
    # Error/Exception: err_text() prepends type(exc).__name__ for every non-ValueError,
    # so StopIteration/KeyboardInterrupt/SystemExit are all producible and were being
    # stripped.  CoolProp's input-pair tags (DmolarSmolar, HmolarP, P2T) are also
    # capitalised, so they are matched by name instead.
    text = re.sub(r'^(?:' + '|'.join(_PAIR_TAGS) + r'): ', '', text)
    bad = _BAD_PHASE.match(text)
    if bad:
        return 'phase {0} instead of {1}'.format(_phase_name(bad.group(1)), _phase_name(bad.group(2)))
    return _NUM.sub('#', text)


def collect(payload):
    """Fold the raw per-fluid records into everything the page renders."""
    backends = payload['backends']
    ref, other = backends[0], backends[1] if len(backends) > 1 else backends[0]

    fluids, pair_names = [], []
    totals = {b: dict(inconsistent=0, exception=0, bad_phase=0, reference=0) for b in backends}
    classes = {b: {} for b in backends}

    for rec in sorted(payload['fluids'], key=lambda r: r['fluid'].lower()):
        row = dict(fluid=rec['fluid'], backends={})
        for b in backends:
            d = rec['backends'].get(b, {})
            t = d.get('totals') or {}
            for k in totals[b]:
                totals[b][k] += int(t.get(k, 0) or 0)
            ref_count = int(t.get('reference', 0) or 0)
            row['backends'][b] = dict(
                ok=bool(d.get('ok')),
                error=d.get('error'),
                crashed=bool(d.get('crashed')),
                wall=d.get('wall'),
                inconsistent=int(t.get('inconsistent', 0) or 0),
                # "exceptions" everywhere below means the pair flash under test threw.
                # Grid points the backend would not define at all are counted apart.
                exceptions=int(t.get('exception', 0) or 0) - ref_count,
                reference=ref_count,
                bad_phase=int(t.get('bad_phase', 0) or 0),
                pairs=d.get('pairs') or {},
            )
            for msg in (d.get('messages') or []):
                key = (message_class(msg['err']), msg.get('kind', 'update'), msg.get('cls', 'EXCEPTION'))
                slot = classes[b].setdefault(key, dict(count=0, fluids=set()))
                slot['count'] += int(msg['count'])
                slot['fluids'].add(rec['fluid'])
            for pair in (d.get('pairs') or {}):
                if pair not in pair_names:
                    pair_names.append(pair)
        fluids.append(row)

    # Per-pair speed: the median over fluids of (other / ref) mean solve time,
    # taken per phase region.  Median, not mean, because a handful of fluids sit
    # orders of magnitude out and would otherwise set the number by themselves.
    def median(values):
        vals = sorted(v for v in values if v is not None)
        if not vals:
            return None
        mid = len(vals) // 2
        return vals[mid] if len(vals) % 2 else 0.5 * (vals[mid - 1] + vals[mid])

    pair_stats = []
    for pair in pair_names:
        entry = dict(pair=pair)
        for region in ('t1', 't2'):
            ratios, ref_times, other_times = [], [], []
            for row in fluids:
                a = row['backends'][ref]['pairs'].get(pair, {}).get(region)
                b = row['backends'][other]['pairs'].get(pair, {}).get(region)
                if a:
                    ref_times.append(a)
                if b:
                    other_times.append(b)
                if a and b:
                    ratios.append(b / a)
            entry[region] = dict(ratio=median(ratios), n=len(ratios),
                                 ref=median(ref_times), other=median(other_times))
        for b in backends:
            entry[b] = dict(
                inconsistent=sum(r['backends'][b]['pairs'].get(pair, {}).get('inconsistent', 0) for r in fluids),
                exceptions=sum((r['backends'][b]['pairs'].get(pair, {}).get('exceptions', 0)
                                - r['backends'][b]['pairs'].get(pair, {}).get('reference', 0))
                               for r in fluids),
                reference=sum(r['backends'][b]['pairs'].get(pair, {}).get('reference', 0) for r in fluids),
                bad_phase=sum(r['backends'][b]['pairs'].get(pair, {}).get('bad_phase', 0) for r in fluids),
            )
        pair_stats.append(entry)

    # Per-fluid median speed ratio, for the fluid table's speed column.
    for row in fluids:
        ratios = []
        for pair in pair_names:
            a = row['backends'][ref]['pairs'].get(pair, {}).get('t1')
            b = row['backends'][other]['pairs'].get(pair, {}).get('t1')
            if a and b:
                ratios.append(b / a)
        row['ratio'] = median(ratios)
        for b in backends:
            # "total" drives the clean / disagreement counts, so it is the defect
            # total: domain refusals are excluded deliberately.
            row['backends'][b]['total'] = (row['backends'][b]['inconsistent']
                                           + row['backends'][b]['exceptions']
                                           + row['backends'][b]['bad_phase'])

    class_lists = {b: sorted(({'cls': k[0], 'kind': k[1], 'failure': k[2], 'count': v['count'],
                               'nfluids': len(v['fluids']), 'fluids': sorted(v['fluids'])}
                              for k, v in classes[b].items()),
                             key=lambda d: -d['count'])
                   for b in backends}

    # `ok` False means the worker died before it could measure anything, so its
    # zero counts are absence of data, not absence of failures.  Excluded from the
    # clean lists and surfaced in their own section -- a segfaulting backend
    # producing a perfect score is the worst way for this page to be wrong.
    clean = {b: [r['fluid'] for r in fluids
                 if r['backends'][b]['ok'] and r['backends'][b]['total'] == 0]
             for b in backends}
    incomplete = [(r['fluid'], b, r['backends'][b].get('error') or 'no result recorded',
                   bool(r['backends'][b].get('crashed')))
                  for r in fluids for b in backends if not r['backends'][b]['ok']]
    return dict(backends=backends, ref=ref, other=other, fluids=fluids, pairs=pair_names,
                pair_stats=pair_stats, totals=totals, classes=class_lists, clean=clean,
                incomplete=incomplete,
                unavailable=payload.get('unavailable', []), generated=payload.get('generated', ''))


CSS = """
:root {
  --ground:#eef1f4; --surface:#ffffff; --surface-2:#f6f8f9; --line:#d3dadf;
  --ink:#131a20; --ink-2:#4d5b66; --ink-3:#78868f;
  --accent:#0d666f; --accent-soft:#dcebec;
  --inc:#9a6512; --exc:#a8342a; --bad:#5f479b;
  --inc-soft:#f4e9d6; --exc-soft:#f6e0dd; --bad-soft:#e6e1f2;
  --fast:#1c6b52; --slow:#a8342a;
}
:root:not([data-theme="light"]) { color-scheme: light dark; }
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --ground:#0d1317; --surface:#141d23; --surface-2:#18232a; --line:#2a373f;
    --ink:#e3eaee; --ink-2:#a3b1ba; --ink-3:#7b8b95;
    --accent:#5cc3ce; --accent-soft:#12333a;
    --inc:#e0aa5c; --exc:#f0857a; --bad:#b09ce8;
    --inc-soft:#33290f; --exc-soft:#3a1f1c; --bad-soft:#241f38;
    --fast:#5fc79f; --slow:#f0857a;
  }
}
:root[data-theme="dark"] {
  --ground:#0d1317; --surface:#141d23; --surface-2:#18232a; --line:#2a373f;
  --ink:#e3eaee; --ink-2:#a3b1ba; --ink-3:#7b8b95;
  --accent:#5cc3ce; --accent-soft:#12333a;
  --inc:#e0aa5c; --exc:#f0857a; --bad:#b09ce8;
  --inc-soft:#33290f; --exc-soft:#3a1f1c; --bad-soft:#241f38;
  --fast:#5fc79f; --slow:#f0857a;
}

* { box-sizing:border-box; }
body {
  margin:0; background:var(--ground); color:var(--ink);
  font-family:"IBM Plex Sans","Helvetica Neue",Arial,sans-serif;
  font-size:15px; line-height:1.6;
}
.wrap { max-width:1180px; margin:0 auto; padding:0 24px 96px; }
h1,h2,h3 { font-family:"IBM Plex Sans Condensed","IBM Plex Sans",Arial,sans-serif;
  text-wrap:balance; margin:0; letter-spacing:-0.01em; }
h1 { font-size:clamp(30px,4.4vw,46px); font-weight:600; line-height:1.08; }
h2 { font-size:24px; font-weight:600; margin-bottom:6px; }
h3 { font-size:17px; font-weight:600; margin-bottom:4px; }
p { margin:0 0 14px; max-width:68ch; color:var(--ink-2); }
p.lede { color:var(--ink); font-size:17px; max-width:66ch; }
a { color:var(--accent); }
code, .mono { font-family:"IBM Plex Mono",ui-monospace,SFMono-Regular,Menlo,monospace; }

header.masthead { border-bottom:1px solid var(--line); background:var(--surface);
  padding:44px 0 26px; margin-bottom:34px; }
.eyebrow { font-family:"IBM Plex Mono",monospace; font-size:11px; letter-spacing:0.16em;
  text-transform:uppercase; color:var(--accent); margin:0 0 12px; }
.provenance { display:grid; grid-template-columns:repeat(auto-fit,minmax(178px,1fr));
  gap:14px 26px; margin-top:26px; padding-top:20px; border-top:1px solid var(--line); }
.provenance div { display:flex; flex-direction:column; gap:2px; }
.provenance dt, .provenance .k { font-size:10.5px; letter-spacing:0.13em; text-transform:uppercase;
  color:var(--ink-3); font-family:"IBM Plex Mono",monospace; }
.provenance .v { font-size:14px; color:var(--ink); font-family:"IBM Plex Mono",monospace; }

section { margin:0 0 46px; }
section > p:first-of-type { margin-top:8px; }
.tiles { display:grid; grid-template-columns:repeat(auto-fit,minmax(210px,1fr)); gap:14px; margin:22px 0 8px; }
.tile { background:var(--surface); border:1px solid var(--line); border-radius:3px; padding:16px 18px;
  display:flex; flex-direction:column; gap:2px; border-top:3px solid var(--accent); }
.tile.inc { border-top-color:var(--inc); } .tile.exc { border-top-color:var(--exc); }
.tile.bad { border-top-color:var(--bad); }
.tile .n { font-size:30px; font-weight:600; font-variant-numeric:tabular-nums;
  font-family:"IBM Plex Sans Condensed",sans-serif; line-height:1.1; }
.tile .lbl { font-size:11px; letter-spacing:0.11em; text-transform:uppercase; color:var(--ink-3);
  font-family:"IBM Plex Mono",monospace; }
.tile .sub { font-size:13px; color:var(--ink-2); }

.tablewrap { overflow-x:auto; border:1px solid var(--line); border-radius:3px; background:var(--surface); }
table { border-collapse:collapse; width:100%; font-size:13.5px; }
th, td { padding:7px 12px; text-align:right; white-space:nowrap; border-bottom:1px solid var(--line); }
th:first-child, td:first-child { text-align:left; }
thead th { position:sticky; top:0; background:var(--surface-2); z-index:1;
  font-family:"IBM Plex Mono",monospace; font-size:10.5px; letter-spacing:0.08em;
  text-transform:uppercase; color:var(--ink-3); font-weight:500; }
thead th.sortable { cursor:pointer; user-select:none; }
thead th.sortable:hover { color:var(--accent); }
thead th .arrow { opacity:0.45; }
tbody tr:hover { background:var(--surface-2); }
td.num { font-variant-numeric:tabular-nums; font-family:"IBM Plex Mono",monospace; }
td.name { font-family:"IBM Plex Mono",monospace; }
.zero { color:var(--ink-3); }
.chip { display:inline-block; min-width:2.4em; padding:1px 7px; border-radius:2px;
  font-family:"IBM Plex Mono",monospace; font-size:12.5px; font-variant-numeric:tabular-nums; }
.chip.inc { background:var(--inc-soft); color:var(--inc); }
.chip.exc { background:var(--exc-soft); color:var(--exc); }
.chip.bad { background:var(--bad-soft); color:var(--bad); }
.grp { border-left:1px solid var(--line); }
.faster { color:var(--fast); } .slower { color:var(--slow); }
tr.detail td { background:var(--surface-2); padding:0; }
tr.detail table { font-size:12.5px; }
tr.detail th { background:transparent; }
tr.expandable td:first-child { cursor:pointer; }
tr.expandable td:first-child::before { content:"\\25B8"; color:var(--ink-3); display:inline-block;
  width:1.1em; transition:transform .12s; }
tr.expandable.open td:first-child::before { transform:rotate(90deg); }

.msg { font-family:"IBM Plex Mono",monospace; font-size:12.5px; white-space:normal;
  text-align:left; color:var(--ink); max-width:none; }
td.msg { white-space:normal; }
.fluidlist { font-family:"IBM Plex Mono",monospace; font-size:12px; color:var(--ink-2);
  white-space:normal; }
.banner { border:1px solid var(--exc); border-left:5px solid var(--exc); background:var(--exc-soft);
  color:var(--exc); padding:13px 17px; border-radius:3px; margin:0 0 24px; font-size:14.5px; }
.banner a { color:var(--exc); font-weight:600; }
.note { border-left:3px solid var(--accent); background:var(--surface); padding:12px 16px;
  margin:18px 0; border-radius:0 3px 3px 0; }
.note p:last-child { margin-bottom:0; }
ul { color:var(--ink-2); max-width:68ch; padding-left:20px; }
li { margin-bottom:6px; }
.controls { display:flex; gap:10px; align-items:center; margin:14px 0 10px; flex-wrap:wrap; }
.controls input { font-family:"IBM Plex Mono",monospace; font-size:13px; padding:6px 10px;
  border:1px solid var(--line); border-radius:3px; background:var(--surface); color:var(--ink); }
.controls input:focus-visible, thead th.sortable:focus-visible { outline:2px solid var(--accent); outline-offset:2px; }
.count { font-size:12.5px; color:var(--ink-3); font-family:"IBM Plex Mono",monospace; }
footer { border-top:1px solid var(--line); padding-top:18px; color:var(--ink-3); font-size:13px; }
@media (prefers-reduced-motion: reduce) { * { transition:none !important; } }
"""


def fmt_time(seconds):
    if not seconds:
        return '—'
    if seconds >= 1e-3:
        return '{0:.3g} ms'.format(seconds * 1e3)
    return '{0:.3g} µs'.format(seconds * 1e6)


def build_html(data, meta):
    ref, other = data['ref'], data['other']
    esc = html.escape
    # '<' is escaped so a '</script>' inside any string (worker error text is the
    # last six lines of arbitrary subprocess output) cannot close the data block
    # early -- which would leave JSON.parse throwing and the fluid table silently
    # empty while the rest of the page still looked right.
    payload = json.dumps(dict(fluids=data['fluids'], pairs=data['pairs'],
                              ref=ref, other=other)).replace('<', '\\u003c')

    def tiles(backend, kind):
        t = data['totals'][backend]
        exc = t['exception'] - t['reference']
        total = t['inconsistent'] + exc + t['bad_phase']
        return (
            '<div class="tile {k}"><span class="lbl">{b} failing points</span>'
            '<span class="n">{n}</span>'
            '<span class="sub">{i:,} inconsistent &middot; {e:,} exception &middot; {p:,} bad-phase'
            '<br>+ {r:,} grid points it declines to define</span></div>'
        ).format(k=kind, b=esc(backend), n='{0:,}'.format(total),
                 i=t['inconsistent'], e=exc, p=t['bad_phase'], r=t['reference'])

    rows = []
    for st in data['pair_stats']:
        def cell(b, key, kind):
            v = st[b][key]
            return ('<span class="chip {k}">{v}</span>'.format(k=kind, v=v) if v
                    else '<span class="zero">0</span>')

        def ratio_cell(r, grp=False):
            klass = 'num grp' if grp else 'num'
            if not r:
                return '<td class="{c} zero">—</td>'.format(c=klass)
            tone = 'slower' if r > 1.15 else ('faster' if r < 0.87 else '')
            return '<td class="{c} {t}">{v}&times;</td>'.format(c=klass, t=tone, v='{0:.2f}'.format(r))
        rows.append(
            '<tr><td class="name">{pair}</td>'
            '<td class="num">{a}</td><td class="num">{b}</td><td class="num">{c}</td>'
            '<td class="num grp">{d}</td><td class="num">{e}</td><td class="num">{f}</td>'
            '<td class="num">{r}</td>'
            '<td class="num grp">{g}</td><td class="num">{h}</td>{i}{j}</tr>'.format(
                pair=esc(st['pair']),
                a=cell(ref, 'inconsistent', 'inc'), b=cell(ref, 'exceptions', 'exc'),
                c=cell(ref, 'bad_phase', 'bad'),
                d=cell(other, 'inconsistent', 'inc'), e=cell(other, 'exceptions', 'exc'),
                f=cell(other, 'bad_phase', 'bad'),
                r=('<span class="zero">' + str(st[other]['reference']) + '</span>'
                   if not st[other]['reference'] else str(st[other]['reference'])),
                g=fmt_time(st['t1']['ref']), h=fmt_time(st['t1']['other']),
                i=ratio_cell(st['t1']['ratio'], grp=True), j=ratio_cell(st['t2']['ratio'])))
    pair_rows = '\n'.join(rows)

    def class_table(backend, kind, limit=18, failure='EXCEPTION'):
        out = []
        kinds = {'update', 'setup'} if kind == 'update' else {kind}
        for c in [c for c in data['classes'][backend]
                  if c['kind'] in kinds and c['failure'] == failure][:limit]:
            fl = ', '.join(c['fluids'][:8]) + ('\u2026' if len(c['fluids']) > 8 else '')
            out.append('<tr><td class="num"><span class="chip {k}">{n}</span></td>'
                       '<td class="num">{nf}</td>'
                       '<td class="msg">{m}<div class="fluidlist">{f}</div></td></tr>'
                       .format(k=('bad' if failure == 'BAD_PHASE' else ('exc' if kind == 'update' else 'inc')),
                               n='{0:,}'.format(c['count']), nf=c['nfluids'],
                               m=esc(c['cls']), f=esc(fl)))
        return '\n'.join(out) or '<tr><td colspan="3">None recorded.</td></tr>'

    unavailable = ', '.join(esc(u['fluid']) for u in data['unavailable']) or 'none'

    crashes = [row for row in data['incomplete'] if row[3]]
    if data['incomplete']:
        rows_i = '\n'.join(
            '<tr><td class="name">{0}</td><td class="name">{1}</td>'
            '<td class="num">{2}</td><td class="msg">{3}</td></tr>'
            .format(esc(f), esc(b),
                    '<span class="chip exc">crash</span>' if crashed else '<span class="zero">error</span>',
                    esc(err))
            for f, b, err, crashed in data['incomplete'])
        incomplete_section = (
            '<section id="incomplete">\n  <h2>Runs that did not complete</h2>\n'
            '  <p>These (fluid, backend) grids produced no measurement, so they count as '
            'nothing anywhere else on this page &mdash; not as clean, and not as failing. '
            'A <strong>crash</strong> means the worker process was killed by a signal: the '
            'backend died in native code, which is a library bug rather than a property of '
            'a state point.</p>\n'
            '  <div class="tablewrap"><table><thead><tr><th>Fluid</th><th>Backend</th>'
            '<th>Kind</th><th style="text-align:left">What happened</th></tr></thead><tbody>\n'
            + rows_i + '\n</tbody></table></div>\n</section>\n')
    else:
        incomplete_section = ''

    # A crash is the one thing on this page that must not wait for the reader to
    # scroll: every other number is silently conditioned on it.
    if crashes:
        crash_banner = (
            '<div class="banner"><strong>{0} worker crash(es).</strong> A backend died in '
            'native code on {1}, so those grids measured nothing and every total below '
            'excludes them. <a href="#incomplete">See what crashed &rarr;</a></div>'
        ).format(len(crashes),
                 esc(', '.join(sorted({'%s (%s)' % (f, b) for f, b, _e, _c in crashes}))[:300]))
    else:
        crash_banner = ''

    return """<title>{title}</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans+Condensed:wght@600&family=IBM+Plex+Sans:wght@400;500;600&display=swap">
<style>{css}</style>

<header class="masthead"><div class="wrap">
  {crash_banner}
  <p class="eyebrow">CoolProp &middot; flash consistency</p>
  <h1>{ref} vs {other} flash audit</h1>
  <p class="lede">{summary}</p>
  <div class="provenance">
    <div><span class="k">Fluids compared</span><span class="v">{nfluids}</span></div>
    <div><span class="k">Input pairs</span><span class="v">{npairs}</span></div>
    <div><span class="k">CoolProp</span><span class="v">{cpver}</span></div>
    <div><span class="k">REFPROP</span><span class="v">{rpver}</span></div>
    <div><span class="k">Grid per panel</span><span class="v">{grid}</span></div>
    <div><span class="k">Generated</span><span class="v">{generated}</span></div>
  </div>
</div></header>

<div class="wrap">

<section>
  <h2>Where the failures are</h2>
  <p>Every panel of the consistency grid flashes from p,T to the panel's input pair and back.
     A point is <em>inconsistent</em> when the round trip lands somewhere else, an
     <em>exception</em> when the flash throws, and <em>bad&#8209;phase</em> when it lands in the
     right state but labels the phase differently than the p,T flash did. Grid points a
     backend declines to define at all &mdash; where the setting-up p,T flash itself throws
     &mdash; never test a flash routine, so they are counted separately throughout.</p>
  <div class="tiles">
    {tile_ref}
    {tile_other}
    <div class="tile"><span class="lbl">Fluids clean in {ref}</span><span class="n">{nclean_ref}</span>
      <span class="sub">of {nfluids}, with no failing point on any panel</span></div>
    <div class="tile"><span class="lbl">Fluids clean in {other}</span><span class="n">{nclean_other}</span>
      <span class="sub">of {nfluids}, on the same grid</span></div>
  </div>
  <div class="note">
    <p>{headline}</p>
  </div>
</section>

<section>
  <h2>By input pair</h2>
  <p>Failing points summed over every fluid, and the median per-point solve time of the
     flashes that succeeded. The ratio is {other} relative to {ref}: above 1 means {other}
     is slower.</p>
  <div class="tablewrap"><table>
    <thead>
      <tr>
        <th rowspan="2">Input pair</th>
        <th colspan="3">{ref} failing points</th>
        <th colspan="4" class="grp">{other} failing points</th>
        <th colspan="2" class="grp">Median solve time (1&#8209;phase)</th>
        <th colspan="2">Time ratio {other}/{ref}</th>
      </tr>
      <tr>
        <th>inc</th><th>exc</th><th>bad</th>
        <th class="grp">inc</th><th>exc</th><th>bad</th><th>declined</th>
        <th class="grp">{ref}</th><th>{other}</th>
        <th class="grp">1&#8209;phase</th><th>2&#8209;phase</th>
      </tr>
    </thead>
    <tbody>{pair_rows}</tbody>
  </table></div>
</section>

<section>
  <h2>By fluid</h2>
  <p>Click a fluid to open its per-pair breakdown. Sort by any column.</p>
  <div class="controls">
    <input id="filter" type="search" placeholder="filter fluids&hellip;" aria-label="Filter fluids">
    <span class="count" id="count"></span>
  </div>
  <div class="tablewrap"><table id="fluidtable">
    <thead><tr>
      <th class="sortable" data-key="fluid" tabindex="0">Fluid <span class="arrow"></span></th>
      <th class="sortable" data-key="ref_inc" tabindex="0">{ref} inc <span class="arrow"></span></th>
      <th class="sortable" data-key="ref_exc" tabindex="0">exc <span class="arrow"></span></th>
      <th class="sortable" data-key="ref_bad" tabindex="0">bad <span class="arrow"></span></th>
      <th class="sortable grp" data-key="oth_inc" tabindex="0">{other} inc <span class="arrow"></span></th>
      <th class="sortable" data-key="oth_exc" tabindex="0">exc <span class="arrow"></span></th>
      <th class="sortable" data-key="oth_bad" tabindex="0">bad <span class="arrow"></span></th>
      <th class="sortable" data-key="oth_ref" tabindex="0">declined <span class="arrow"></span></th>
      <th class="sortable grp" data-key="ratio" tabindex="0">median time ratio <span class="arrow"></span></th>
    </tr></thead>
    <tbody id="fluidbody"></tbody>
  </table></div>
</section>

<section>
  <h2>What {ref} throws</h2>
  <p>Message text with the numeric literals masked, so one class counts once however many
     state points hit it.</p>
  <div class="tablewrap"><table>
    <thead><tr><th style="text-align:right">points</th><th>fluids</th><th style="text-align:left">message class</th></tr></thead>
    <tbody>{classes_ref}</tbody>
  </table></div>
</section>

<section>
  <h2>What {other} throws</h2>
  <p>The pair flash under test, failing on a state point that {other} itself defined.</p>
  <div class="tablewrap"><table>
    <thead><tr><th style="text-align:right">points</th><th>fluids</th><th style="text-align:left">message class</th></tr></thead>
    <tbody>{classes_other}</tbody>
  </table></div>

  <h3 style="margin-top:30px">&hellip;and what it declines to define</h3>
  <p>These are the p,T (or Q,T) flashes that set up the comparison, refused before any pair
     flash ran. They are a statement about {other}'s validated domain, not about its flash
     routines &mdash; {ref} answers the same points by extrapolating past the equation's
     stated limits. Counted apart everywhere above.</p>
  <div class="tablewrap"><table>
    <thead><tr><th style="text-align:right">points</th><th>fluids</th><th style="text-align:left">message class</th></tr></thead>
    <tbody>{declined_other}</tbody>
  </table></div>
</section>

{incomplete_section}
<section>
  <h2>Phase-label disagreements</h2>
  <p>Compared inside one backend: the p,T flash and the pair flash reach the same state,
     agreeing on density, pressure and temperature to within 0.1%, and then label its
     phase differently. {ref} records none. The two-phase check is still suppressed for
     {other} &mdash; the phase-index convention genuinely differs inside the dome
     (GitHub #1057) &mdash; so everything below is single-phase.</p>
  <div class="tablewrap"><table>
    <thead><tr><th style="text-align:right">points</th><th>fluids</th><th style="text-align:left">pair flash said &hellip; where p,T said &hellip;</th></tr></thead>
    <tbody>{badphase_other}</tbody>
  </table></div>
</section>

<section>
  <h2>How to read this</h2>
  <ul>
    <li><strong>Grid.</strong> {grid} per panel, {npairs} input pairs per fluid &mdash; the same grid
        the documentation build renders. Four panels are crossed out for both backends
        (<span class="mono">SmolarUmolar</span>, <span class="mono">HmolarUmolar</span>,
        <span class="mono">HmolarT</span>, <span class="mono">TUmolar</span>).</li>
    <li><strong>Timing.</strong> Mean wall time per point over the flashes that <em>solved</em>.
        A backend that throws early on many points would look fast on an all-points mean.
        Both backends run back to back inside one worker process, so the ratio survives the
        CPU contention of a parallel run; the absolute times do not &mdash; treat them as
        indicative, not as a benchmark.</li>
    <li><strong>Is the REFPROP side charged for CoolProp's wrapper?</strong> Barely. Timing the
        raw <span class="mono">TDFLSHdll</span> / <span class="mono">TPFLSHdll</span> calls
        against <span class="mono">AbstractState::update</span> on identical inputs in C++,
        with no Python in the path, puts CoolProp's own per-call work at
        <strong>0.04&ndash;0.17 &micro;s, or 2&ndash;5% of the measured time</strong>
        (Water, CO<sub>2</sub>, n&#8209;Propane, Nitrogen; 20,000 points each). The
        <span class="mono">WMOLdll</span> call CoolProp makes on every update costs 2&ndash;3 ns.
        What these columns compare is the flash routines, not the binding.</li>
    <li><strong>Bad-phase.</strong> Compared within one backend: the p,T flash and the pair
        flash must agree on the phase label for the same state. The two-phase check stays
        suppressed for REFPROP (the phase-index convention genuinely differs inside the dome,
        GitHub #1057); the single-phase check now runs for both.</li>
    <li><strong>Not in REFPROP.</strong> <span class="mono">{unavailable}</span> &mdash; present in
        CoolProp's fluid list with no REFPROP equivalent, so they are out of scope here.</li>
  </ul>
</section>

<footer>Generated by <span class="mono">dev/scripts/consistency_backend_compare.py</span> and
<span class="mono">dev/scripts/consistency_backend_report.py</span>. Full failing-point CSVs
accompany this page.</footer>
</div>

<script id="data" type="application/json">{payload}</script>
<script>
(function () {{
  var D = JSON.parse(document.getElementById('data').textContent);
  var body = document.getElementById('fluidbody');
  var countEl = document.getElementById('count');
  var sortKey = 'oth_total', sortDir = -1, filter = '';

  function flat(row) {{
    var a = row.backends[D.ref], b = row.backends[D.other];
    return {{
      fluid: row.fluid, row: row,
      ref_inc: a.inconsistent, ref_exc: a.exceptions, ref_bad: a.bad_phase,
      oth_inc: b.inconsistent, oth_exc: b.exceptions, oth_bad: b.bad_phase, oth_ref: b.reference,
      ok: (a.ok !== false) && (b.ok !== false),
      oth_total: b.inconsistent + b.exceptions + b.bad_phase + a.inconsistent + a.exceptions + a.bad_phase,
      ratio: row.ratio
    }};
  }}
  var rows = D.fluids.map(flat);

  function num(v, cls) {{
    if (!v) return '<td class="num zero">0</td>';
    return '<td class="num"><span class="chip ' + cls + '">' + v + '</span></td>';
  }}
  function ratioCell(r) {{
    if (!r) return '<td class="num grp zero">&mdash;</td>';
    var tone = r > 1.15 ? 'slower' : (r < 0.87 ? 'faster' : '');
    return '<td class="num grp ' + tone + '">' + r.toFixed(2) + '&times;</td>';
  }}
  function fmt(t) {{
    if (!t) return '&mdash;';
    return t >= 1e-3 ? (t * 1e3).toPrecision(3) + ' ms' : (t * 1e6).toPrecision(3) + ' \\u00b5s';
  }}

  function detail(row) {{
    var out = '<tr class="detail"><td colspan="9"><div class="tablewrap"><table><thead><tr>' +
      '<th>pair</th><th>' + D.ref + ' inc</th><th>exc</th><th>bad</th>' +
      '<th class="grp">' + D.other + ' inc</th><th>exc</th><th>bad</th><th>declined</th>' +
      '<th class="grp">' + D.ref + ' 1ph</th><th>' + D.other + ' 1ph</th><th class="grp">ratio</th>' +
      '</tr></thead><tbody>';
    D.pairs.forEach(function (p) {{
      var a = row.backends[D.ref].pairs[p], b = row.backends[D.other].pairs[p];
      if (!a && !b) return;
      a = a || {{}}; b = b || {{}};
      var r = (a.t1 && b.t1) ? b.t1 / a.t1 : null;
      out += '<tr><td class="name">' + p + '</td>' +
        num(a.inconsistent, 'inc') + num((a.exceptions || 0) - (a.reference || 0), 'exc') +
        num(a.bad_phase, 'bad') +
        '<td class="num grp">' + (b.inconsistent ? '<span class="chip inc">' + b.inconsistent + '</span>' : '<span class="zero">0</span>') + '</td>' +
        num((b.exceptions || 0) - (b.reference || 0), 'exc') + num(b.bad_phase, 'bad') +
        '<td class="num zero">' + (b.reference || 0) + '</td>' +
        '<td class="num grp">' + fmt(a.t1) + '</td><td class="num">' + fmt(b.t1) + '</td>' +
        ratioCell(r) + '</tr>';
    }});
    return out + '</tbody></table></div></td></tr>';
  }}

  function render() {{
    var shown = rows.filter(function (r) {{
      return !filter || r.fluid.toLowerCase().indexOf(filter) >= 0;
    }});
    shown.sort(function (x, y) {{
      var a = x[sortKey], b = y[sortKey];
      if (typeof a === 'string') return sortDir * a.localeCompare(b);
      a = (a === null || a === undefined) ? -1 : a;
      b = (b === null || b === undefined) ? -1 : b;
      return sortDir * (a - b) || x.fluid.localeCompare(y.fluid);
    }});
    body.innerHTML = shown.map(function (r) {{
      // A grid that produced no measurement must not render as a row of zeroes next to
      // the genuinely clean fluids -- that is the one reading the section text forbids.
      if (!r.ok) {{
        return '<tr class="expandable" data-fluid="' + r.fluid + '"><td class="name">' + r.fluid + '</td>' +
          '<td class="num" colspan="7"><span class="chip exc">no result &mdash; see "Runs that did not complete"</span></td></tr>';
      }}
      return '<tr class="expandable" data-fluid="' + r.fluid + '"><td class="name">' + r.fluid + '</td>' +
        num(r.ref_inc, 'inc') + num(r.ref_exc, 'exc') + num(r.ref_bad, 'bad') +
        '<td class="num grp">' + (r.oth_inc ? '<span class="chip inc">' + r.oth_inc + '</span>' : '<span class="zero">0</span>') + '</td>' +
        num(r.oth_exc, 'exc') + num(r.oth_bad, 'bad') +
        '<td class="num zero">' + (r.oth_ref || 0) + '</td>' + ratioCell(r.ratio) + '</tr>';
    }}).join('');
    countEl.textContent = shown.length + ' of ' + rows.length + ' fluids';
  }}

  body.addEventListener('click', function (ev) {{
    var tr = ev.target.closest('tr.expandable');
    if (!tr) return;
    if (tr.classList.contains('open')) {{
      tr.classList.remove('open');
      if (tr.nextElementSibling && tr.nextElementSibling.classList.contains('detail')) {{
        tr.nextElementSibling.remove();
      }}
      return;
    }}
    var rec = rows.filter(function (r) {{ return r.fluid === tr.dataset.fluid; }})[0];
    tr.classList.add('open');
    tr.insertAdjacentHTML('afterend', detail(rec.row));
  }});

  document.querySelectorAll('th.sortable').forEach(function (th) {{
    function go() {{
      var k = th.dataset.key;
      if (sortKey === k) {{ sortDir = -sortDir; }} else {{ sortKey = k; sortDir = k === 'fluid' ? 1 : -1; }}
      document.querySelectorAll('th.sortable .arrow').forEach(function (a) {{ a.textContent = ''; }});
      th.querySelector('.arrow').textContent = sortDir > 0 ? '\\u2191' : '\\u2193';
      render();
    }}
    th.addEventListener('click', go);
    th.addEventListener('keydown', function (e) {{ if (e.key === 'Enter' || e.key === ' ') {{ e.preventDefault(); go(); }} }});
  }});
  document.getElementById('filter').addEventListener('input', function (e) {{
    filter = e.target.value.trim().toLowerCase();
    render();
  }});
  render();
}})();
</script>
""".format(
        title=esc(meta['title']), css=CSS, ref=esc(ref), other=esc(other),
        summary=esc(meta['summary']), nfluids=len(data['fluids']), npairs=len(data['pairs']),
        cpver=esc(meta['coolprop']), rpver=esc(meta['refprop']), grid=esc(meta['grid']),
        generated=esc(data['generated']),
        tile_ref=tiles(ref, 'exc'), tile_other=tiles(other, 'inc'),
        nclean_ref=len(data['clean'][ref]), nclean_other=len(data['clean'][other]),
        headline=meta['headline'],
        pair_rows=pair_rows,
        classes_ref=class_table(ref, 'update'), classes_other=class_table(other, 'update'),
        declined_other=class_table(other, 'reference'),
        badphase_other=class_table(other, 'update', failure='BAD_PHASE'),
        unavailable=unavailable, incomplete_section=incomplete_section,
        crash_banner=crash_banner, payload=payload)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--json', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--title', default='HEOS vs REFPROP Flash Audit')
    parser.add_argument('--summary', default='')
    parser.add_argument('--coolprop', default='')
    parser.add_argument('--refprop', default='')
    parser.add_argument('--grid', default='40x40 (1-phase), 20x20 (2-phase)')
    parser.add_argument('--headline', default='', help='HTML for the callout under the tiles')
    args = parser.parse_args(argv)

    with open(args.json) as fp:
        payload = json.load(fp)
    data = collect(payload)
    # The single biggest failure class is worth stating up front: one systematic
    # mislabelling can account for a five-figure point count, and a reader who only
    # sees the total will read it as a five-figure number of independent defects.
    # Only real failures: the domain refusals ('reference') are the largest class of
    # all, and leading with them would headline the one thing this report argues is
    # not a defect.
    biggest = None
    for b in (data['other'], data['ref']):
        for c in data['classes'][b]:
            if c['kind'] == 'reference':
                continue
            if biggest is None or c['count'] > biggest[1]['count']:
                biggest = (b, c)
    headline = args.headline
    if not headline and biggest:
        headline = ('Largest single failure class: <strong>{0:,} points</strong> across {1} fluids '
                    'on the {2} backend, all of them <span class="mono">{3}</span>. One systematic '
                    'cause, not {0:,} independent defects &mdash; read the totals with that in mind.').format(
            biggest[1]['count'], biggest[1]['nfluids'], html.escape(biggest[0]),
            html.escape(biggest[1]['cls'][:160]))
    meta = dict(title=args.title, summary=args.summary or (
        'Every CoolProp fluid REFPROP also carries, put through the same flash-consistency '
        'grid on both backends.'),
        coolprop=args.coolprop, refprop=args.refprop, grid=args.grid, headline=headline)
    with open(os.path.abspath(os.path.expanduser(args.out)), 'w') as fp:
        fp.write(build_html(data, meta))
    print('Wrote', args.out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
