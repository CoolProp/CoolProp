"""Chebyshev fits for the caloric properties (density, specific heat).

Emits the optional ``<property>_cheb`` JSON entries that sit alongside the
classic centered-polynomial entries::

    "specific_heat_cheb": {
        "type": "chebyshev",
        "Trange": [Tmin, Tmax],   # the fit domain in T, explicit
        "xbase": 0.21,            # centering of the composition direction
        "coeffs": [[...], ...],   # rows k: Chebyshev-in-T coefficients,
                                  # cols j: multiplied by (x - xbase)^j
        "NRMS": 0.001
    }

so that ``value(T, x) = sum_k sum_j C[k][j] * T_k(u(T)) * (x - xbase)^j``
with ``u(T) = (2*T - (Tmax + Tmin)) / (Tmax - Tmin)``.

Two provenance paths, mirroring how the polynomial fits are produced:

- fluids with tabular data are least-squares fitted from the raw grid
  (NaN-masked, T-degree capped by the number of temperature points);
- coefficient-only fluids (Melinder book, food correlations, LiBr, ...)
  have no raw data -- their committed polynomial IS the ground truth, so
  the Chebyshev coefficients are an EXACT basis conversion of it
  (interpolation at Chebyshev-Lobatto nodes at the same degree).

The math was validated end-to-end in prototype_chebyshev_caloric.py
(entropy from these fits agrees with adaptive quadrature to ~1e-16; see
NOTES_thermodynamic_consistency.md).
"""

import numpy as np
from numpy.polynomial import chebyshev as ncheb

# T-degree for data fits: enough for every smooth fluid (the prototype's
# fit-quality tables plateau at or before 8), capped by the data.
DEFAULT_MAX_DEGREE_T = 8
# A Chebyshev fit in T needs at least a handful of temperatures; freeze
# curves and conversion tables with fewer points are skipped (they are not
# caloric properties anyway).
MIN_TEMPERATURE_POINTS = 3

CALORIC_PROPERTIES = ("density", "specific_heat")


def _scaled_T(T, Trange):
    return (2.0 * np.asarray(T, dtype=float) - (Trange[1] + Trange[0])) / (Trange[1] - Trange[0])


def evaluate(coeffs, T, x, Trange, xbase):
    """Evaluate a 2D chebyshev-in-T x monomial-in-(x - xbase) fit."""
    coeffs = np.asarray(coeffs, dtype=float)
    dx = np.power(x - xbase, np.arange(coeffs.shape[1]))
    return ncheb.chebval(_scaled_T(T, Trange), coeffs @ dx)


def fit_from_data(Tvec, xvec, grid, xbase, deg_x, deg_T=DEFAULT_MAX_DEGREE_T):
    """Least-squares 2D fit from a raw data grid.

    Returns (coeffs, Trange, NRMS) or None if the grid cannot support a fit.
    """
    Tvec = np.asarray(Tvec, dtype=float).ravel()
    xvec = np.asarray(xvec, dtype=float).ravel()
    grid = np.asarray(grid, dtype=float).reshape(len(Tvec), len(xvec))
    finite_T_rows = int(np.isfinite(grid).any(axis=1).sum())
    if finite_T_rows < MIN_TEMPERATURE_POINTS:
        return None
    Trange = (float(np.min(Tvec)), float(np.max(Tvec)))
    if Trange[1] <= Trange[0]:
        return None

    TT, XX = np.meshgrid(Tvec, xvec, indexing="ij")
    T, x, z = TT.ravel(), XX.ravel(), grid.ravel()
    mask = np.isfinite(z)
    points = int(mask.sum())

    # Fit only as many coefficients as the FINITE data supports: sparse
    # SecCool grids pad out-of-range cells with NaN, and some hardcoded
    # arrays (e.g. LiqNa cp) carry NaN gaps.
    deg_T_cap = min(deg_T, finite_T_rows - 1)
    deg_x = min(deg_x, max(len(xvec) - 1, 0))
    while deg_x > 0 and 2 * (deg_x + 1) > points:
        deg_x -= 1
    if 2 * (deg_x + 1) > points:
        return None

    u = _scaled_T(T[mask], Trange)
    xpow = np.vander(x[mask] - xbase, deg_x + 1, increasing=True)
    spread = float(np.nanmax(grid) - np.nanmin(grid))

    # Select the T-degree by generalized cross-validation, GCV = n*RSS/(n-p)^2:
    # tabulated data is typically rounded to 3-4 digits, and a high-order fit
    # of a near-linear property just rings between the data points (0.5 kg/m3
    # wiggles on ExamplePure's density). GCV keeps the order low for noisy
    # low-structure data while still choosing high order where the data
    # genuinely bends (the ice slurries' cp).
    best = None
    for candidate in range(1, deg_T_cap + 1):
        n_coeffs = (candidate + 1) * (deg_x + 1)
        if n_coeffs >= points:  # interpolation has no leftover DOF to validate
            break
        cheb_cols = ncheb.chebvander(u, candidate)
        design = (cheb_cols[:, :, None] * xpow[:, None, :]).reshape(points, -1)
        flat, *_ = np.linalg.lstsq(design, z[mask], rcond=None)
        rss = float(np.sum(np.square(design @ flat - z[mask])))
        gcv = points * rss / (points - n_coeffs) ** 2
        if best is None or gcv < best[0]:
            best = (gcv, candidate, flat, rss)
    if best is None:
        # too few points for any validated fit: fall back to the exactly-
        # determined lowest order that the points allow
        candidate = max(min(deg_T_cap, points // (deg_x + 1) - 1), 1)
        cheb_cols = ncheb.chebvander(u, candidate)
        design = (cheb_cols[:, :, None] * xpow[:, None, :]).reshape(points, -1)
        flat, *_ = np.linalg.lstsq(design, z[mask], rcond=None)
        best = (0.0, candidate, flat, float(np.sum(np.square(design @ flat - z[mask]))))

    _, chosen_deg_T, flat, rss = best
    coeffs = flat.reshape(chosen_deg_T + 1, deg_x + 1)
    residual_rms = float(np.sqrt(rss / points))
    nrms = residual_rms / spread if spread > 0 else residual_rms
    return coeffs, Trange, nrms


def convert_polynomial(poly_coeffs, Tbase, xbase, Trange):
    """EXACT basis conversion of a centered 2D polynomial fit.

    The committed polynomial ``sum_k sum_j P[k][j] (T-Tbase)^k (x-xbase)^j``
    is, per x-column, a degree-K polynomial in T; interpolating it at K+1
    Chebyshev-Lobatto nodes reproduces it exactly (same function space), so
    this adds no fitting error. The x-direction (monomial in x - xbase) is
    carried over unchanged.
    """
    P = np.asarray(poly_coeffs, dtype=float)
    if P.ndim == 1:
        P = P.reshape(-1, 1)
    deg_T = P.shape[0] - 1
    out = np.zeros_like(P)
    for j in range(P.shape[1]):
        col = ncheb.Chebyshev.interpolate(
            lambda T: np.polynomial.polynomial.polyval(np.asarray(T) - Tbase, P[:, j]),
            deg_T, domain=list(Trange))
        out[:deg_T + 1, j] = col.coef
    return out


def _positive_on_domain(coeffs, Trange, xbase, xmin, xmax):
    """Density and heat capacity are strictly positive; a fit that swings
    non-positive anywhere on the fluid's (T, x) domain is oscillating
    through data-free regions and must not ship."""
    Ts = np.linspace(Trange[0], Trange[1], 41)
    xs = np.linspace(xmin, xmax, 9) if xmax > xmin else [xbase]
    return all(np.min(evaluate(coeffs, Ts, x, Trange, xbase)) > 0.0 for x in xs)


def _fit_covers_fluid_range(rawT, rawGrid, Tmin, Tmax):
    """A tabular refit is only trustworthy when the finite data spans
    (essentially all of) the fluid's advertised temperature range --
    otherwise the entry would extrapolate freely where the committed
    low-order polynomial extrapolates gently (LiqNa's cp data stops at
    1000 K of an advertised 2500 K, see DATA_AUDIT.md)."""
    grid = np.asarray(rawGrid, dtype=float).reshape(np.size(rawT), -1)
    T = np.asarray(rawT, dtype=float).ravel()[np.isfinite(grid).any(axis=1)]
    if T.size < MIN_TEMPERATURE_POINTS:
        return False
    return (T.max() - T.min()) >= 0.9 * (Tmax - Tmin)


def build_entry(fluid_json, prop, rawT=None, rawX=None, rawGrid=None):
    """Build the ``<prop>_cheb`` entry dict for one caloric property.

    fluid_json is the committed per-fluid JSON dict (source of the
    polynomial coefficients, Tbase/xbase, Tmin/Tmax). If a raw data grid
    covering the fluid's range is supplied, the entry is fitted from it,
    lowering the fit orders until the result is positive over the whole
    domain; otherwise (no data, poor coverage, or no positive fit) it is
    an exact basis conversion of the committed polynomial. Returns None
    when nothing usable exists.
    """
    committed = fluid_json.get(prop, {})
    xbase = float(fluid_json.get("xbase", 0.0) or 0.0)
    xmin = float(fluid_json.get("xmin", 0.0) or 0.0)
    xmax = float(fluid_json.get("xmax", 0.0) or 0.0)
    has_poly = committed.get("type") == "polynomial" and committed.get("coeffs") not in (None, "null")

    if (rawGrid is not None and rawT is not None
            and _fit_covers_fluid_range(rawT, rawGrid, float(fluid_json["Tmin"]), float(fluid_json["Tmax"]))):
        deg_x0 = max(len(committed["coeffs"][0]) - 1, 0) if has_poly \
            else (min(5, np.size(rawX) - 1) if rawX is not None and np.size(rawX) > 1 else 0)
        # Highest orders first; the first fit that is positive over the whole
        # domain wins. Sparse grids (SecCool -1 padding) often need lower
        # orders than dense ones to avoid oscillating between data bands.
        candidates = sorted(
            ((dT, dx) for dT in range(DEFAULT_MAX_DEGREE_T, 1, -1) for dx in range(deg_x0, -1, -1)),
            key=lambda p: (p[0] + p[1], p[0]), reverse=True)
        for deg_T, deg_x in candidates:
            fitted = fit_from_data(rawT, rawX if rawX is not None else [0.0], rawGrid, xbase, deg_x, deg_T)
            if fitted is None:
                continue
            coeffs, Trange, nrms = fitted
            if np.all(np.isfinite(coeffs)) and _positive_on_domain(coeffs, Trange, xbase, xmin, xmax):
                return {
                    "type": "chebyshev",
                    "Trange": [Trange[0], Trange[1]],
                    "xbase": xbase,
                    "coeffs": coeffs.tolist(),
                    "NRMS": nrms,
                    "fit_source": "tabular_data",
                }

    # Fall back to the exact re-basis of the committed polynomial: it is the
    # shipped ground truth, including its (gentle) extrapolation behavior.
    if not has_poly:
        return None
    Trange = (float(fluid_json["Tmin"]), float(fluid_json["Tmax"]))
    if Trange[1] <= Trange[0]:
        return None
    Tbase = float(fluid_json.get("Tbase", 0.0) or 0.0)
    coeffs = convert_polynomial(committed["coeffs"], Tbase, xbase, Trange)
    if not np.all(np.isfinite(coeffs)):
        return None
    return {
        "type": "chebyshev",
        "Trange": [Trange[0], Trange[1]],
        "xbase": xbase,
        "coeffs": coeffs.tolist(),
        "NRMS": committed.get("NRMS"),
        "fit_source": "basis_conversion",
    }
