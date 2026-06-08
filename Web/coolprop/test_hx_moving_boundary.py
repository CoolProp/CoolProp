"""Regression guard for the provider-based moving-boundary HX solver.

Runs from the notebook's directory so `import hx_moving_boundary` resolves.
Exercisable with `python3 -m pytest test_hx_moving_boundary.py -v` or directly
with `python3 test_hx_moving_boundary.py`.
"""
from hx_moving_boundary import make_evaporator, ORACLE_Q


def test_heos_matches_oracle():
    # The Python port must reproduce dev/reference/HX.py's Q to tight tolerance.
    hx = make_evaporator('HEOS', A=4.0)
    Q = hx.run()
    assert abs(Q - ORACLE_Q) / ORACLE_Q < 1e-6


def test_backends_agree():
    Qh = make_evaporator('HEOS', A=4.0).run()
    Qs = make_evaporator('SVDSBTL&HEOS', A=4.0).run()
    assert abs(Qh - Qs) / Qh < 1e-5


def test_effectiveness_physical():
    hx = make_evaporator('HEOS', A=4.0)
    Q = hx.run()
    eps = Q / hx.Qmax
    assert 0.0 < eps < 1.0


if __name__ == '__main__':
    test_heos_matches_oracle()
    test_backends_agree()
    test_effectiveness_physical()
    print('all tests passed')
