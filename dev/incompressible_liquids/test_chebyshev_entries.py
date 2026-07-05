"""Checks on the committed ``*_cheb`` caloric entries in json/*.json.

Three layers:
- schema/shape sanity for every committed entry (stdlib-level checks);
- physical sanity: density and cp evaluated from the Chebyshev fit are
  strictly positive across each fluid's whole (T, x) domain -- this is the
  net that catches a sparse-grid fit oscillating wildly between data;
- agreement with the committed polynomial entry: exact (basis conversions)
  or fit-level (tabular refits);
- a golden-master refit (DowJ, LiBr, AKF), mirroring
  test_fitting_regression.py, so numpy/scipy drift in the Chebyshev fitter
  is caught before release;
- entropy machinery: h/s built from the committed cheb cp entries agree
  with adaptive quadrature (the production integral path validated in
  prototype_chebyshev_caloric.py, applied to the committed data).
"""

import glob
import json
import os

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("scipy")

from numpy.polynomial import chebyshev as ncheb
from scipy.integrate import quad

from CPIncomp import ChebyshevFits

JSON_DIR = os.path.join(os.path.dirname(__file__), "json")
CALORIC = ChebyshevFits.CALORIC_PROPERTIES


def _fluids():
    for path in sorted(glob.glob(os.path.join(JSON_DIR, "*.json"))):
        with open(path) as fh:
            yield os.path.basename(path)[:-5], json.load(fh)


def _domain_grid(fluid, entry):
    T0, T1 = entry["Trange"]
    Ts = np.linspace(T0, T1, 25)
    xmin, xmax = float(fluid.get("xmin", 0.0)), float(fluid.get("xmax", 0.0))
    xs = np.linspace(xmin, xmax, 5) if xmax > xmin else np.array([0.0])
    return Ts, xs


def _poly_eval(committed, Tbase, xbase, T, x):
    P = np.asarray(committed["coeffs"], dtype=float)
    if P.ndim == 1:
        P = P.reshape(-1, 1)
    T = np.asarray(T, dtype=float)
    return np.polynomial.polynomial.polyval2d(T - Tbase, np.full_like(T, x) - xbase, P)


def test_every_caloric_property_has_a_cheb_entry():
    missing = []
    for name, fluid in _fluids():
        for prop in CALORIC:
            has_poly = fluid.get(prop, {}).get("type") == "polynomial"
            if has_poly and prop + "_cheb" not in fluid:
                missing.append((name, prop))
    assert not missing, missing


def test_cheb_entries_schema():
    for name, fluid in _fluids():
        for prop in CALORIC:
            entry = fluid.get(prop + "_cheb")
            if entry is None:
                continue
            assert entry["type"] == "chebyshev", name
            T0, T1 = entry["Trange"]
            assert 0.0 < T0 < T1, (name, prop, entry["Trange"])
            assert entry["fit_source"] in ("tabular_data", "basis_conversion"), name
            coeffs = np.asarray(entry["coeffs"], dtype=float)
            assert coeffs.ndim == 2 and coeffs.size, (name, prop)
            assert np.all(np.isfinite(coeffs)), (name, prop)


def test_density_and_cp_positive_across_domain():
    bad = []
    for name, fluid in _fluids():
        for prop in CALORIC:
            entry = fluid.get(prop + "_cheb")
            if entry is None:
                continue
            Ts, xs = _domain_grid(fluid, entry)
            for x in xs:
                vals = ChebyshevFits.evaluate(entry["coeffs"], Ts, x, entry["Trange"], entry["xbase"])
                if np.min(vals) <= 0.0:
                    bad.append((name, prop, float(x), float(np.min(vals))))
    assert not bad, bad


def test_conversions_reproduce_committed_polynomial_exactly():
    for name, fluid in _fluids():
        Tbase = float(fluid.get("Tbase", 0.0) or 0.0)
        for prop in CALORIC:
            entry = fluid.get(prop + "_cheb")
            committed = fluid.get(prop, {})
            if entry is None or entry["fit_source"] != "basis_conversion" or committed.get("type") != "polynomial":
                continue
            Ts, xs = _domain_grid(fluid, entry)
            for x in xs:
                cheb = ChebyshevFits.evaluate(entry["coeffs"], Ts, x, entry["Trange"], entry["xbase"])
                poly = _poly_eval(committed, Tbase, entry["xbase"], Ts, x)
                rel = np.max(np.abs(cheb - poly) / np.maximum(np.abs(poly), 1e-30))
                assert rel < 1e-9, (name, prop, float(x), rel)


def test_tabular_fits_describe_their_data():
    # For refitted entries, compare against the actual data (and against the
    # committed polynomial AT the data points). Comparing the two fits away
    # from the data would only measure how differently they extrapolate into
    # data-free corners of sparse grids, where neither is authoritative.
    from add_chebyshev_entries import collect_fluid_objects, raw_grids

    objects = collect_fluid_objects()
    for name, fluid in _fluids():
        Tbase = float(fluid.get("Tbase", 0.0) or 0.0)
        for prop in CALORIC:
            entry = fluid.get(prop + "_cheb")
            committed = fluid.get(prop, {})
            if entry is None or entry["fit_source"] != "tabular_data":
                continue
            assert entry["NRMS"] is None or entry["NRMS"] < 0.06, (name, prop, entry["NRMS"])
            rawT, rawX, rawGrid = raw_grids(objects.get(name), prop)
            if rawGrid is None:
                continue
            rawT = np.asarray(rawT, dtype=float).ravel()
            rawX = np.asarray(rawX if rawX is not None else [0.0], dtype=float).ravel()
            grid = np.asarray(rawGrid, dtype=float).reshape(len(rawT), len(rawX))
            spread = float(np.nanmax(grid) - np.nanmin(grid)) or 1.0
            for jx, x in enumerate(rawX):
                col = grid[:, jx]
                mask = np.isfinite(col)
                if not mask.any():
                    continue
                cheb = ChebyshevFits.evaluate(entry["coeffs"], rawT[mask], x, entry["Trange"], entry["xbase"])
                # fit describes the data
                assert np.max(np.abs(cheb - col[mask])) / spread < 0.15, (name, prop, float(x))
                # and stays close to the committed fit where the data lives
                # (except the synthetic Example* templates, whose committed
                # polynomial is a known-poor fit -- 13.8% NRMS, see DATA_AUDIT)
                if committed.get("type") == "polynomial" and not name.startswith("Example"):
                    poly = _poly_eval(committed, Tbase, entry["xbase"], rawT[mask], x)
                    assert np.max(np.abs(cheb - poly)) / spread < 0.20, (name, prop, float(x))


def test_refit_golden_master():
    # Mirrors test_fitting_regression.py: rebuild the cheb entries for three
    # representative fluids from raw data and compare with what's committed.
    from add_chebyshev_entries import collect_fluid_objects, raw_grids

    objects = collect_fluid_objects()
    for name in ("DowJ", "LiBr", "AKF"):
        with open(os.path.join(JSON_DIR, name + ".json")) as fh:
            fluid = json.load(fh)
        for prop in CALORIC:
            committed = fluid.get(prop + "_cheb")
            if committed is None:
                continue
            rawT, rawX, rawGrid = raw_grids(objects[name], prop)
            rebuilt = ChebyshevFits.build_entry(fluid, prop, rawT, rawX, rawGrid)
            assert rebuilt is not None, (name, prop)
            a, b = np.asarray(rebuilt["coeffs"]), np.asarray(committed["coeffs"])
            assert a.shape == b.shape, (name, prop)
            scale = np.max(np.abs(b))
            assert np.allclose(a, b, rtol=1e-2, atol=1e-6 * scale), (name, prop)


def test_entropy_from_committed_cheb_matches_quadrature():
    # The production integral path (Lobatto re-expansion of cp/T, then the
    # Chebyshev antiderivative) applied to the committed entries.
    for name in ("Water", "LiBr", "MEG", "LiqNa"):
        with open(os.path.join(JSON_DIR, name + ".json")) as fh:
            fluid = json.load(fh)
        entry = fluid.get("specific_heat_cheb")
        assert entry is not None, name
        T0f, T1f = entry["Trange"]
        x = 0.5 * (float(fluid.get("xmin", 0.0)) + float(fluid.get("xmax", 0.0)))
        coeffs = np.asarray(entry["coeffs"], dtype=float)
        dx = np.power(x - entry["xbase"], np.arange(coeffs.shape[1]))
        cp1d = ncheb.Chebyshev(coeffs @ dx, domain=[T0f, T1f])

        degree = coeffs.shape[0] + 8
        while True:
            aux = ncheb.Chebyshev.interpolate(lambda T: cp1d(T) / T, degree, domain=[T0f, T1f])
            tail = np.max(np.abs(aux.coef[-3:])) / np.max(np.abs(aux.coef))
            if tail < 1e-13 or degree >= 96:
                break
            degree += 16
        s_fun = aux.integ()

        Ta, Tb = T0f + 0.1 * (T1f - T0f), T1f - 0.1 * (T1f - T0f)
        ds = s_fun(Tb) - s_fun(Ta)
        truth, _ = quad(lambda T: cp1d(T) / T, Ta, Tb, epsabs=1e-13, epsrel=1e-13, limit=200)
        assert abs(ds - truth) / abs(truth) < 1e-12, (name, ds, truth)
