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


def test_consolidated_unavailable_section(tmp_path):
    """Fluids the backend does not carry are a coverage gap, listed apart from
    the build failures so the failure list stays about real defects."""
    out = tmp_path / 'ConsistencyReport_REFPROP.rst'
    rpt.write_consolidated_rst(pandas.DataFrame(), str(out), 'REFPROP',
                               build_failures=[('Foo', 'kaboom')],
                               unavailable=[('SES36', 'Could not load these fluids: SES36')],
                               date='2026-05-26', orphan=True)
    text = out.read_text(encoding='utf-8')
    assert text.startswith(':orphan:')
    assert 'Not available in this backend' in text
    assert '1 fluid(s) in the CoolProp fluid list have no REFPROP equivalent' in text
    assert 'SES36: Could not load these fluids: SES36' in text
    # The two lists stay distinct.
    assert 'SES36' not in text.split('Not available in this backend')[0]
    assert 'Foo: kaboom' in text.split('Not available in this backend')[0]


def test_consolidated_no_unavailable_section_when_empty(tmp_path):
    out = tmp_path / 'ConsistencyReport.rst'
    rpt.write_consolidated_rst(pandas.DataFrame(), str(out), 'HEOS', date='2026-05-26')
    assert 'Not available in this backend' not in out.read_text(encoding='utf-8')


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


def test_lowest_valid_T_prefers_the_equation_minimum():
    """The grid floor is the colder of the two limits the backend reports, not the
    triple point alone: REFPROP publishes a true triple point below the range its
    equation is fitted over."""
    from CoolProp.Plots.ConsistencyPlots import lowest_valid_T
    import CoolProp.CoolProp as CP

    class Fake(object):
        def __init__(self, T_triple, T_min):
            self._t, self._m = T_triple, T_min

        def keyed_output(self, key):
            if key == CP.iT_triple:
                return self._t
            if key == CP.iT_min:
                return self._m
            raise KeyError(key)

    assert lowest_valid_T(Fake(89.54, 120.0)) == 120.0   # REFPROP R14
    assert lowest_valid_T(Fake(273.16, 251.165)) == 273.16  # REFPROP Water (extended range)
    assert lowest_valid_T(Fake(200.0, 200.0)) == 200.0   # HEOS: the two coincide

    class NoTmin(Fake):
        def keyed_output(self, key):
            if key == CP.iT_triple:
                return self._t
            raise ValueError('T_min not available for this backend')

    # A backend with no T_min at all must degrade to the triple point, not throw.
    assert lowest_valid_T(NoTmin(150.0, None)) == 150.0


def test_setup_failure_frame_is_one_visible_exception():
    """A panel that cannot be set up must be reported, not silently empty --
    an empty frame is indistinguishable from a clean panel."""
    from CoolProp.Plots.ConsistencyPlots import ConsistencyFigure
    import matplotlib
    matplotlib.use('Agg')
    ff = ConsistencyFigure('Water', backend='HEOS',
                           NT_1phase=3, Np_1phase=3, NT_2phase=3, NQ_2phase=3)
    df = ff.axes_list[0].setup_failure_frame('no usable P_min')
    assert len(df) == 1
    assert df['cls'].iloc[0] == 'EXCEPTION'
    assert df['pair'].iloc[0] == ff.axes_list[0].pair
    assert 'no usable P_min' in df['err'].iloc[0]
    import matplotlib.pyplot as plt
    plt.close(ff.fig)


def test_good_only_timing_is_recorded_separately():
    """The annotation keeps the all-points mean; the GOOD-only twin exists for a
    backend-to-backend comparison, where counting throw-fast points as fast lies."""
    from CoolProp.Plots.ConsistencyPlots import ConsistencyFigure
    import matplotlib
    matplotlib.use('Agg')
    ff = ConsistencyFigure('Water', backend='HEOS',
                           NT_1phase=4, Np_1phase=4, NT_2phase=3, NQ_2phase=3)
    good = [a.mean_elapsed_1phase_good for a in ff.axes_list
            if getattr(a, 'mean_elapsed_1phase_good', None) is not None]
    assert good and all(m > 0 for m in good)
    import matplotlib.pyplot as plt
    plt.close(ff.fig)


if __name__ == '__main__':
    import pytest
    pytest.main([__file__, '-v'])
