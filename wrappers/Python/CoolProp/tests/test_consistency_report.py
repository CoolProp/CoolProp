import pandas
import pytest
from CoolProp.Plots import _consistency_report as rpt


def _sample_errors():
    return pandas.DataFrame([
        dict(fluid='Water', backend='HEOS', pair='HmolarP', phase_region='1phase',
             cls='EXCEPTION', in1='Hmolar', val1=1.0, in2='P', val2=2.0, P=2.0, T=300.0,
             dev=float('nan'), err='boom'),
        dict(fluid='Water', backend='HEOS', pair='HmolarP', phase_region='1phase',
             cls='EXCEPTION', in1='Hmolar', val1=1.1, in2='P', val2=2.1, P=2.1, T=301.0,
             dev=float('nan'), err='boom'),
        dict(fluid='Water', backend='HEOS', pair='DmolarP', phase_region='2phase',
             cls='INCONSISTENT', in1='Dmolar', val1=5.0, in2='P', val2=6.0, P=6.0, T=350.0,
             dev=0.4, err=''),
        dict(fluid='Water', backend='HEOS', pair='DmolarP', phase_region='2phase',
             cls='BAD_PHASE', in1='Dmolar', val1=7.0, in2='P', val2=8.0, P=8.0, T=360.0,
             dev=float('nan'), err='phase 6 instead of 5'),
    ])


def test_summarize_counts():
    s = rpt.summarize(_sample_errors())
    by_pair = {r['pair']: r for _, r in s.iterrows()}
    assert by_pair['HmolarP']['exceptions'] == 2
    assert by_pair['HmolarP']['total_failures'] == 2
    assert by_pair['DmolarP']['inconsistent'] == 1
    assert by_pair['DmolarP']['bad_phase'] == 1
    assert by_pair['DmolarP']['total_failures'] == 2
    # sorted by total_failures descending
    assert list(s['total_failures']) == sorted(s['total_failures'], reverse=True)


def test_summarize_empty():
    s = rpt.summarize(pandas.DataFrame())
    assert list(s.columns) == ['pair', 'inconsistent', 'exceptions', 'bad_phase', 'total_failures']
    assert len(s) == 0


def test_format_time_units():
    assert rpt.format_time(2.0).endswith(' s')
    assert rpt.format_time(2e-3).endswith(' ms')
    assert rpt.format_time(2e-6).endswith(' µs')
    assert rpt.format_time(None) == 'n/a'
    assert rpt.format_time(float('nan')) == 'n/a'


def test_write_csv_columns_and_rows(tmp_path):
    path = tmp_path / 'errors.csv'
    rpt.write_csv(_sample_errors(), str(path))
    out = pandas.read_csv(path)
    assert list(out.columns) == rpt._CSV_COLS
    assert len(out) == 4  # all four failing rows


def test_write_csv_empty(tmp_path):
    path = tmp_path / 'empty.csv'
    rpt.write_csv(pandas.DataFrame(), str(path))
    out = pandas.read_csv(path)
    assert list(out.columns) == rpt._CSV_COLS
    assert len(out) == 0


def test_fragment_no_failures(tmp_path):
    out = tmp_path / 'frag.rst'
    rpt.write_rst_fragment(pandas.DataFrame(), 'Water', 'HEOS',
                           'Consistencyplots/Water-consistency.csv', str(out))
    text = out.read_text(encoding='utf-8')
    assert 'no failures' in text
    assert '.. dropdown::' not in text


def test_fragment_with_failures(tmp_path):
    out = tmp_path / 'frag.rst'
    rpt.write_rst_fragment(_sample_errors(), 'Water', 'HEOS',
                           'Consistencyplots/Water-consistency.csv', str(out))
    text = out.read_text(encoding='utf-8')
    assert 'Flash consistency (HEOS)' in text
    assert '2 exceptions' in text
    assert ':download:' in text
    assert 'Consistencyplots/Water-consistency.csv' in text
    assert '.. dropdown::' in text
    assert '.. list-table::' in text


def test_fragment_intro_prepended(tmp_path):
    out = tmp_path / 'frag.rst'
    rpt.write_rst_fragment(_sample_errors(), 'Water', 'REFPROP',
                           'Consistencyplots_REFPROP/Water-consistency.csv', str(out),
                           intro_rst='.. image:: Consistencyplots_REFPROP/Water.png\n')
    text = out.read_text(encoding='utf-8')
    assert text.startswith('.. image:: Consistencyplots_REFPROP/Water.png')


def test_sample_cap_and_dedup():
    # 5 EXCEPTION rows, same err -> dedup to 1; INCONSISTENT capped at 2
    rows = []
    for i in range(5):
        rows.append(dict(fluid='X', backend='HEOS', pair='HmolarP', phase_region='1phase',
                         cls='EXCEPTION', in1='Hmolar', val1=i, in2='P', val2=i, P=i, T=i,
                         dev=float('nan'), err='same'))
    for i in range(5):
        rows.append(dict(fluid='X', backend='HEOS', pair='DmolarP', phase_region='1phase',
                         cls='INCONSISTENT', in1='Dmolar', val1=i, in2='P', val2=i, P=i, T=i,
                         dev=float(i), err=''))
    sample = rpt._select_sample(pandas.DataFrame(rows), sample_cap=2)
    exc = sample[sample['cls'] == 'EXCEPTION']
    inc = sample[sample['cls'] == 'INCONSISTENT']
    assert len(exc) == 1                 # deduped on err message
    assert len(inc) == 2                 # capped
    assert list(inc['dev']) == [4.0, 3.0]  # largest-deviation first


def test_consolidated_page(tmp_path):
    combined = pandas.DataFrame([
        dict(fluid='Water', pair='HmolarP', inconsistent=0, exceptions=2, bad_phase=0, total_failures=2),
        dict(fluid='R134a', pair='DmolarP', inconsistent=1, exceptions=0, bad_phase=1, total_failures=2),
        dict(fluid='Air', pair='PT', inconsistent=0, exceptions=0, bad_phase=0, total_failures=0),
    ])
    out = tmp_path / 'ConsistencyReport.rst'
    rpt.write_consolidated_rst(combined, str(out), 'HEOS',
                               build_failures=[('Foo', 'kaboom')], date='2026-05-26')
    text = out.read_text(encoding='utf-8')
    assert 'Consistency Failure Report (HEOS)' in text
    assert '2026-05-26' in text
    assert ':ref:`Water <fluid_Water>`' in text
    assert 'Air' not in text.split('Build failures')[0]  # zero-failure fluid omitted from table
    assert 'Build failures' in text
    assert 'Foo: kaboom' in text


def test_consolidated_empty(tmp_path):
    out = tmp_path / 'ConsistencyReport.rst'
    rpt.write_consolidated_rst(pandas.DataFrame(), str(out), 'HEOS', date='2026-05-26')
    text = out.read_text(encoding='utf-8')
    assert 'No failures recorded' in text


def test_write_stub(tmp_path):
    out = tmp_path / 'stub.rst'
    rpt.write_stub_fragment(str(out), 'REFPROP consistency not generated in this build.')
    assert out.read_text(encoding='utf-8').strip() == 'REFPROP consistency not generated in this build.'


def test_instrumentation_columns_and_timing():
    # Tiny grid so this runs fast; asserts the errors DataFrame is self-describing.
    from CoolProp.Plots.ConsistencyPlots import ConsistencyFigure
    import matplotlib
    matplotlib.use('Agg')
    ff = ConsistencyFigure('Water', backend='HEOS',
                           NT_1phase=4, Np_1phase=4, NT_2phase=3, NQ_2phase=3)
    err = ff.errors
    # fluid/backend stamped
    assert (err['fluid'] == 'Water').all() if len(err) else True
    assert (err['backend'] == 'HEOS').all() if len(err) else True
    # required columns exist
    for col in ['cls', 'pair', 'phase_region']:
        assert col in err.columns
    # at least one panel recorded a single-phase mean time
    means = [a.mean_elapsed_1phase for a in ff.axes_list
             if getattr(a, 'mean_elapsed_1phase', None) is not None]
    assert means and all(m > 0 for m in means)
    import matplotlib.pyplot as plt
    plt.close(ff.fig)


def test_panel_timing_annotation():
    from CoolProp.Plots.ConsistencyPlots import ConsistencyFigure
    import matplotlib
    matplotlib.use('Agg')
    ff = ConsistencyFigure('Water', backend='HEOS',
                           NT_1phase=4, Np_1phase=4, NT_2phase=3, NQ_2phase=3)
    # Every implemented panel that ran should carry a mean-time annotation.
    annotated = 0
    for ax in ff.axes_list:
        texts = [t.get_text() for t in ax.ax.texts]
        if any('mean t/pt' in s for s in texts):
            annotated += 1
    assert annotated > 0
    import matplotlib.pyplot as plt
    plt.close(ff.fig)


if __name__ == '__main__':
    import pytest
    pytest.main([__file__, '-v'])
