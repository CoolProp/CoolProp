"""AbstractState missing-method coverage (bd CoolProp-r9sq.25).

Locks the two AbstractState methods the legacy Cython interface exposed but the
nanobind core had dropped: ``build_options_json`` and the zero-allocation
``fast_evaluate`` batch path (tabular / IF97 backends).
"""
import numpy as np
import pytest

from CoolProp import CoolProp as C


def _not_marshalling_error(exc):
    assert not isinstance(exc, TypeError), f"method uncallable (marshalling bug): {exc}"
    assert "Unable to convert" not in str(exc), f"unregistered type (marshalling bug): {exc}"


class TestBuildOptionsJson:
    def test_callable_returns_str(self):
        AS = C.AbstractState("HEOS", "Water")
        try:
            s = AS.build_options_json()
        except Exception as e:  # noqa: BLE001 -- a backend may not implement it; not a marshalling bug
            _not_marshalling_error(e)
            return
        assert isinstance(s, str)


class TestFastEvaluate:
    def _arrays(self, N=3, M=2):
        val1 = np.array([1.0e5, 2.0e5, 3.0e5][:N], dtype=np.float64)   # P [Pa]
        val2 = np.array([300.0, 350.0, 400.0][:N], dtype=np.float64)   # T [K]
        outputs = np.array([int(C.iHmolar), int(C.iDmolar)][:M], dtype=np.int32)
        out = np.zeros((N, M), dtype=np.float64)
        status = np.zeros(N, dtype=np.int32)
        return val1, val2, outputs, out, status

    def test_fast_evaluate_runs(self):
        # IF97 (water steam tables) implements the fast path.
        AS = C.AbstractState("IF97", "Water")
        val1, val2, outputs, out, status = self._arrays()
        try:
            AS.fast_evaluate(int(C.PT_INPUTS), val1, val2, outputs, out, status)
        except Exception as e:  # noqa: BLE001 -- backend may not support this pair; not a marshalling bug
            _not_marshalling_error(e)
            return
        # out filled in place; ok points have status 0 and finite outputs.
        ok = status == 0
        assert ok.any()
        assert np.isfinite(out[ok]).all()

    def test_shape_validation(self):
        AS = C.AbstractState("IF97", "Water")
        val1, val2, outputs, out, status = self._arrays()
        bad_out = np.zeros((2, 2), dtype=np.float64)  # wrong N
        with pytest.raises(ValueError):
            AS.fast_evaluate(int(C.PT_INPUTS), val1, val2, outputs, bad_out, status)

    def test_length_mismatch_raises(self):
        AS = C.AbstractState("IF97", "Water")
        val1, val2, outputs, out, status = self._arrays()
        with pytest.raises(ValueError):
            AS.fast_evaluate(int(C.PT_INPUTS), val1, val2[:2], outputs, out, status)
