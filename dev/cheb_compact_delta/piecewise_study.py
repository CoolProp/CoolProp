"""Piecewise subdivision study: for each IAPWS-95 term family, at several tau,
how many adaptive subintervals of u = (delta - 1) / (delta + 1) are needed
to represent the per-term pressure contribution g(u) at max per-piece degree
n_max in {4, 8, 12, 16} with tail tolerance 1e-12.

This tells us whether the compactified-delta piecewise architecture can hit
the single-digit-microsecond budget (which requires small n per piece).
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from iapws95 import IAPWS95Terms
from compact_cheb import values_to_coeffs, chebyshev_lobatto_nodes
from convergence_study import (
    _single_power_term_fn,
    _single_gauss_term_fn,
    _single_nonan_term_fn,
)

OUTDIR = os.path.dirname(os.path.abspath(__file__))


def adaptive_subdivide(f, a, b, n_max, tol_rel=1e-12, tol_abs=1e-300,
                       min_width=1e-10, max_pieces=2000):
    """Return list of (aa, bb, coeffs) where each piece has a Chebyshev series
    of degree <= n_max satisfying tail-decay tol. Piece-level convergence
    criterion: ||last tail_len coeffs|| < tol_rel * max(|c|) + tol_abs."""
    stack = [(a, b)]
    pieces = []
    while stack:
        if len(pieces) + len(stack) > max_pieces:
            raise RuntimeError(f"Exceeded {max_pieces} pieces")
        aa, bb = stack.pop()
        nodes = chebyshev_lobatto_nodes(n_max, aa, bb)
        try:
            vals = np.array([f(x) for x in nodes])
        except Exception:
            vals = np.full(n_max + 1, np.nan)
        if not np.all(np.isfinite(vals)):
            if bb - aa < min_width:
                pieces.append((aa, bb, np.full(n_max + 1, np.nan)))
                continue
            mid = 0.5 * (aa + bb)
            stack.append((mid, bb))
            stack.append((aa, mid))
            continue
        coeffs = values_to_coeffs(vals)
        abs_c = np.abs(coeffs)
        cmax = abs_c.max() if abs_c.max() > 0 else 1.0
        tail_len = max(2, (n_max + 1) // 4)
        tail = abs_c[-tail_len:]
        if np.max(tail) < tol_rel * cmax + tol_abs:
            pieces.append((aa, bb, coeffs))
            continue
        if bb - aa < min_width:
            pieces.append((aa, bb, coeffs))
            continue
        mid = 0.5 * (aa + bb)
        stack.append((mid, bb))
        stack.append((aa, mid))
    pieces.sort(key=lambda p: p[0])
    return pieces


def study_family(eos, family_name, term_fn_list, labels, tau_values,
                 n_max_values=(4, 8, 12, 16), a=-0.99999, b=0.99999):
    print(f"\n### {family_name} ###")
    header = f"{'term':<20} {'tau':>8} " + "  ".join([f"n={nm:<2d} pieces"
                                                        for nm in n_max_values])
    print(header)
    results = []
    for idx, (term_fn, lbl) in enumerate(zip(term_fn_list, labels)):
        for tau in tau_values:
            def f(u, _t=tau, _fn=term_fn):
                return _fn(u, _t)
            row = f"{lbl:<20} {tau:>8.3f} "
            piece_counts_for_label = []
            for nm in n_max_values:
                try:
                    pieces = adaptive_subdivide(f, a, b, n_max=nm, tol_rel=1e-12)
                    row += f"  {len(pieces):>6d}    "
                    piece_counts_for_label.append(len(pieces))
                except Exception as e:
                    row += f"  ERR       "
                    piece_counts_for_label.append(-1)
            results.append((lbl, tau, piece_counts_for_label))
            print(row)
    return results


def main():
    eos = IAPWS95Terms()
    tau_values = [1.5, 1.0, 0.8]  # subcritical, critical, supercritical

    # Pick representative terms per family
    # (idx into the power block: one l=0, one each l=1,2,3,4,6)
    pwr_reps = [
        (0, "power_l0_i0"),
        (7, "power_l1_i7"),
        (22, "power_l2_i22"),
        (42, "power_l3_i42"),
        (46, "power_l4_i46"),
        (50, "power_l6_i50"),
    ]
    pwr_fns = [_single_power_term_fn(eos, i) for i, _ in pwr_reps]
    pwr_lbls = [lbl for _, lbl in pwr_reps]

    study_family(eos, "Power terms", pwr_fns, pwr_lbls, tau_values)

    gauss_fns = [_single_gauss_term_fn(eos, i) for i in range(3)]
    gauss_lbls = [f"gauss_i{i}" for i in range(3)]
    study_family(eos, "Gaussian bell terms", gauss_fns, gauss_lbls, tau_values)

    nonan_fns = [_single_nonan_term_fn(eos, i) for i in range(2)]
    nonan_lbls = [f"nonan_i{i}" for i in range(2)]
    study_family(eos, "Non-analytic terms", nonan_fns, nonan_lbls, tau_values)

    # Also study the WORST-CASE combined per-term function (a typical fluid
    # term may be the pressure contribution for just the polyexp+gaussian sum)
    def combined_polyexp_gauss(u, tau, eos=eos):
        from compact_cheb import delta_of_u
        delta = delta_of_u(u)
        pe = eos.alphar_01_power_terms(tau, delta).sum() - \
            sum(eos.alphar_01_power_terms(tau, delta)[eos.l_p == 0])
        gg = eos.alphar_01_gauss_terms(tau, delta).sum()
        return delta * (pe + gg)  # delta * (polyexp+gauss α^r_δ) contribution to P/(ρ_c RT)

    print(f"\n### Combined (44 polyexp + 3 gauss, delta * sum) ###")
    header = f"{'function':<20} {'tau':>8} " + "  ".join([f"n={nm:<2d} pieces"
                                                            for nm in (4, 8, 12, 16, 24)])
    print(header)
    for tau in tau_values:
        def f(u, _t=tau):
            return combined_polyexp_gauss(u, _t)
        row = f"{'polyexp+gauss':<20} {tau:>8.3f} "
        for nm in (4, 8, 12, 16, 24):
            pieces = adaptive_subdivide(f, -0.99999, 0.99999, n_max=nm, tol_rel=1e-12)
            row += f"  {len(pieces):>6d}    "
        print(row)


if __name__ == "__main__":
    main()
