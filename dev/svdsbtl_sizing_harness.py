#!/usr/bin/env python3
"""SVDSBTL per-region sizing harness (CoolProp-6oe).

For a given (fluid, source, input_pair, region), this script measures
the data the sizing decision actually depends on:

  1. σ_k / σ_1 singular-value decay
  2. on-grid rank-vs-err curve (reconstruction error of the rank-r SVD
     at the training grid points)
  3. off-grid rank-vs-err curve (the metric that matches runtime
     behavior: rank-r SVD + natural-cubic Hermite interpolation,
     evaluated at random (xi, eta) probes against HEOS truth)
  4. NT / NR sensitivity at fixed rank — finds the grid-saturation
     point past which more samples stop buying accuracy

Output: one CSV per (fluid, region, property) plus a matplotlib figure
per region.  Reports a recommended (NT, NR, rank) triple for each
(region, property) given a target off-grid max-error budget.

This harness deliberately reimplements the region geometry in Python
rather than calling the C++ SurfacePresets.  Python is fast enough for
the offline analysis (n_probes = 2000-5000 per cell takes seconds),
and the standalone implementation lets us experiment with axis
transforms / boundary curves that haven't been productionized yet.

Usage:
    /Users/ianbell/Code/CoolProp/.venv/bin/python dev/svdsbtl_sizing_harness.py Argon V
    /Users/ianbell/Code/CoolProp/.venv/bin/python dev/svdsbtl_sizing_harness.py R245fa L --region NC --primary pow1/3
    /Users/ianbell/Code/CoolProp/.venv/bin/python dev/svdsbtl_sizing_harness.py Water L --interactive

Region tags:
    L      LIQUID  (parent subcritical liquid leg)
    V      VAPOR   (parent subcritical vapor leg)
    S      SUPER   (parent supercritical)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

import CoolProp.CoolProp as CP
from CoolProp import AbstractState


# ---------------------- region geometry ----------------------

def h_sat(state, p, Q):
    try:
        state.update(CP.PQ_INPUTS, p, Q)
        return state.hmass()
    except Exception:
        return np.nan


def h_pt(state, p, T):
    try:
        state.update(CP.PT_INPUTS, p, T)
        return state.hmass()
    except Exception:
        return np.nan


def h_floor_robust(state, p, T_min):
    """Production-matching h_floor that walks T up over the melting line."""
    try:
        T_melt = state.melting_line(CP.iT, CP.iP, p)
        T_eff = max(T_min, T_melt + 0.01)
    except Exception:
        T_eff = T_min
    return h_pt(state, p, T_eff)


def make_region(state, region):
    """Returns (p_lo, p_hi, h_lo_fn(p), h_hi_fn(p)) for the named region.

    Geometry mirrors src/SBTL/SurfacePresets.cpp.  LIQUID/VAPOR/SUPER
    are the parent regions; the NC sub-regions are addressed by
    --p-lo / --p-hi CLI overrides (no separate region keys here so the
    harness stays general-purpose).
    """
    p_trip = state.keyed_output(CP.iP_triple) * 1.01
    pc = state.p_critical()
    p_max = state.pmax() * 0.99
    T_max = state.Tmax()
    T_min = max(state.Ttriple(), state.Tmin())
    if region == "L":  # LIQUID
        return p_trip, 0.999 * pc, lambda p: h_floor_robust(state, p, T_min), lambda p: h_sat(state, p, 0)
    if region == "V":  # VAPOR
        return p_trip, 0.999 * pc, lambda p: h_sat(state, p, 1), lambda p: h_pt(state, p, T_max - 0.5)
    if region == "S":  # SUPER
        return 1.001 * pc, p_max, lambda p: h_floor_robust(state, p, T_min), lambda p: h_pt(state, p, T_max - 0.5)
    raise ValueError(f"unknown region {region!r}")


# ---------------------- grid construction ----------------------

def cheb_eta(NT, pad=0.001):
    """Chebyshev-cosine on η with 0.001 pad (matches make_grid_axes in SVDSurfaceFactory.cpp)."""
    t = np.arange(NT) / (NT - 1)
    return pad + (1.0 - 2.0 * pad) * 0.5 * (1.0 - np.cos(np.pi * t))


def primary_axis(primary, p_lo, p_hi):
    """Returns (fwd: xi -> p, inv: p -> xi) for the named primary axis."""
    if primary == "log":
        log_lo, log_hi = np.log(p_lo), np.log(p_hi)
        return (lambda xi: np.exp(log_lo + xi * (log_hi - log_lo)),
                lambda p: (np.log(p) - log_lo) / (log_hi - log_lo))
    if primary == "pow1/3":
        # POWER: forward xi = 1 - cbrt((a_hi - p)/(a_hi - a_lo)) — crowds at a_hi
        return (lambda xi: p_hi - (p_hi - p_lo) * (1.0 - xi) ** 3,
                lambda p: 1.0 - np.cbrt((p_hi - p) / (p_hi - p_lo)))
    if primary == "pow1/2":
        return (lambda xi: p_hi - (p_hi - p_lo) * (1.0 - xi) ** 2,
                lambda p: 1.0 - np.sqrt((p_hi - p) / (p_hi - p_lo)))
    raise ValueError(f"unknown primary axis {primary!r}")


# ---------------------- sampling + SVD ----------------------

PROPS = {
    "rho": (lambda s: s.rhomass(), "log"),
    "T":   (lambda s: s.T(),       "identity"),
    "s":   (lambda s: s.smass(),   "identity"),
    "u":   (lambda s: s.umass(),   "identity"),
}


def sample_grid(state, geom, NT, NR, prop, primary):
    """Returns M (NT × NR) sample matrix indexed by (i=η, j=ξ)."""
    p_lo, p_hi, h_lo_fn, h_hi_fn = geom
    eta_grid = cheb_eta(NT)
    xi_grid = np.linspace(0, 1, NR)
    fwd, _ = primary_axis(primary, p_lo, p_hi)
    p_grid = fwd(xi_grid)
    get, _ = PROPS[prop]
    M = np.full((NT, NR), np.nan)
    for j, p in enumerate(p_grid):
        b_lo, b_hi = h_lo_fn(p), h_hi_fn(p)
        if not (np.isfinite(b_lo) and np.isfinite(b_hi) and b_lo < b_hi):
            continue
        for i, eta in enumerate(eta_grid):
            h = b_lo + eta * (b_hi - b_lo)
            try:
                state.update(CP.HmassP_INPUTS, h, p)
                v = get(state)
                if np.isfinite(v) and v > 0:
                    M[i, j] = v
            except Exception:
                pass
    return M, eta_grid, xi_grid, p_grid


def svd_analyze(M, prop_transform="log"):
    """Compute σ decay and on-grid rank-vs-err.  Returns (s, sigma_ratio,
    [(rank, max, p99, p50), ...])."""
    M_t = np.log(M) if prop_transform == "log" else M.copy()
    mask = np.isfinite(M_t)
    M_filled = M_t.copy()
    for i in range(M_t.shape[0]):
        row = M_t[i, :]
        good = np.isfinite(row)
        if good.any() and (~good).any():
            M_filled[i, ~good] = np.median(row[good])
        elif not good.any():
            M_filled[i, :] = 0.0
    U, s, Vt = np.linalg.svd(M_filled, full_matrices=False)
    ranks_full = sorted(set([1, 2, 3, 5, 8, 10, 12, 15, 20, 25, 30, 40, 50, len(s)]))
    ranks = [r for r in ranks_full if r <= len(s)]
    out = []
    for r in ranks:
        M_r = U[:, :r] @ np.diag(s[:r]) @ Vt[:r, :]
        if prop_transform == "log":
            M_orig = np.exp(M_t); M_r_orig = np.exp(M_r)
        else:
            M_orig = M_t; M_r_orig = M_r
        e = np.abs(M_r_orig - M_orig) / np.abs(M_orig)
        e_finite = e[mask]
        out.append((r, float(e_finite.max()), float(np.quantile(e_finite, 0.99)), float(np.quantile(e_finite, 0.50))))
    return s, s / s[0], out


def build_offgrid_eval(M, eta_grid, xi_grid, rank, prop_transform):
    """Mirror of SVDEvaluator: per-mode natural-cubic Hermite interpolation."""
    M_t = np.log(M) if prop_transform == "log" else M.copy()
    for i in range(M_t.shape[0]):
        good = np.isfinite(M_t[i, :])
        if good.any() and (~good).any(): M_t[i, ~good] = np.median(M_t[i, good])
        elif not good.any(): M_t[i, :] = 0.0
    U, s, Vt = np.linalg.svd(M_t, full_matrices=False)
    U_r, S_r, V_r = U[:, :rank], s[:rank], Vt[:rank, :]
    u_sp = [CubicSpline(eta_grid, U_r[:, k], bc_type="natural") for k in range(rank)]
    v_sp = [CubicSpline(xi_grid, V_r[k, :], bc_type="natural") for k in range(rank)]

    def ev(xi, eta):
        uv = np.array([u_sp[k](eta) for k in range(rank)])
        vv = np.array([v_sp[k](xi) for k in range(rank)])
        val = float(np.sum(S_r * uv * vv))
        return np.exp(val) if prop_transform == "log" else val
    return ev


def offgrid_probe(state, geom, primary, eta_grid, xi_grid, M_dict, ranks, prop, n_probes=3000, seed=0):
    """Probes (xi, eta) ∈ (eps, 1-eps)^2 at the chosen ranks.

    Returns (errs, n_valid, probes) where:
      errs    : dict[rank] -> rel-err array (one entry per valid probe)
      n_valid : number of valid probes (HEOS-evaluable + property > 0)
      probes  : dict with keys 'xi', 'eta', 'p', 'h' — all length n_probes,
                NaN where invalid.  Lets the caller plot worst-cell location.
    """
    p_lo, p_hi, h_lo_fn, h_hi_fn = geom
    fwd, _ = primary_axis(primary, p_lo, p_hi)
    rng = np.random.default_rng(seed)
    eps = 1e-3
    xis = rng.uniform(eps, 1 - eps, n_probes)
    etas = rng.uniform(eps, 1 - eps, n_probes)
    p_arr = fwd(xis)
    get, transform = PROPS[prop]
    truth = np.full(n_probes, np.nan)
    h_arr = np.full(n_probes, np.nan)
    for i in range(n_probes):
        p = p_arr[i]
        b_lo, b_hi = h_lo_fn(p), h_hi_fn(p)
        if not (np.isfinite(b_lo) and np.isfinite(b_hi) and b_lo < b_hi):
            continue
        h = b_lo + etas[i] * (b_hi - b_lo)
        h_arr[i] = h
        try:
            state.update(CP.HmassP_INPUTS, h, p)
            v = get(state)
            if np.isfinite(v) and v > 0:
                truth[i] = v
        except Exception:
            pass
    valid = np.isfinite(truth)
    out = {}
    M = M_dict[prop]
    for r in ranks:
        evaluator = build_offgrid_eval(M, eta_grid, xi_grid, r, transform)
        # Returned arrays are length-n_probes with NaN at invalid probes
        # so xi/eta/p/h indexes line up with errs[r].
        e_full = np.full(n_probes, np.nan)
        for i in range(n_probes):
            if not valid[i]: continue
            try:
                pred = evaluator(xis[i], etas[i])
                e_full[i] = abs(pred / truth[i] - 1)
            except Exception:
                pass
        out[r] = e_full
    probes = {"xi": xis, "eta": etas, "p": p_arr, "h": h_arr}
    return out, int(valid.sum()), probes


# ---------------------- recommender ----------------------

def recommend(offgrid_curve, budget_max=1e-6):
    """Given dict[rank] -> err array, return (rank, achieved_max).  rank is
    the smallest tested rank whose max-err meets the budget; if none
    meet it, returns the rank with the smallest max-err."""
    items = sorted(offgrid_curve.items())
    for r, e in items:
        if e.max() <= budget_max:
            return r, float(e.max())
    r_best, e_best = min(items, key=lambda kv: kv[1].max())
    return r_best, float(e_best.max())


# ---------------------- driver ----------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fluid")
    ap.add_argument("region", choices=["L", "V", "S"], help="L=LIQUID, V=VAPOR, S=SUPER (parent)")
    ap.add_argument("--source", default="HEOS")
    ap.add_argument("--primary", default="log", choices=["log", "pow1/2", "pow1/3"])
    ap.add_argument("--p-lo-mult", type=float, default=None, help="Override p_lo = mult · pc (for NC sub-bands)")
    ap.add_argument("--p-hi-mult", type=float, default=None, help="Override p_hi = mult · pc")
    ap.add_argument("--NT-list", type=int, nargs="+", default=[100, 200, 400])
    ap.add_argument("--NR-list", type=int, nargs="+", default=[400, 800])
    ap.add_argument("--ranks", type=int, nargs="+", default=[5, 10, 15, 20, 25, 30, 40])
    ap.add_argument("--probes", type=int, default=3000)
    ap.add_argument("--budget", type=float, default=1e-6, help="Target max off-grid relative error")
    ap.add_argument("--props", nargs="+", default=["rho", "T", "s", "u"])
    ap.add_argument("--interactive", action="store_true")
    ap.add_argument("--out", default="/tmp/svdsbtl_sizing")
    args = ap.parse_args()

    if args.interactive:
        import matplotlib
        matplotlib.use("MacOSX")
    import matplotlib.pyplot as plt

    state = AbstractState(args.source, args.fluid)
    pc = state.p_critical()
    geom = make_region(state, args.region)
    if args.p_lo_mult is not None or args.p_hi_mult is not None:
        p_lo = (args.p_lo_mult * pc) if args.p_lo_mult is not None else geom[0]
        p_hi = (args.p_hi_mult * pc) if args.p_hi_mult is not None else geom[1]
        geom = (p_lo, p_hi, geom[2], geom[3])
    p_lo, p_hi, *_ = geom

    label = f"{args.fluid} {args.region} {args.source} [{p_lo/pc:.3f}, {p_hi/pc:.3f}]·pc  primary={args.primary}"
    print(f"\n=== {label} ===")
    print(f"probes={args.probes}  budget={args.budget:.0e}  props={args.props}")

    out_dir = Path(args.out); out_dir.mkdir(parents=True, exist_ok=True)
    # Three rows: (0) off-grid rank curve, (1) σ-decay, (2) worst-cell (ξ, η) localization
    fig, axes = plt.subplots(3, len(args.props), figsize=(4*len(args.props), 11), squeeze=False)

    rows = []  # CSV rows
    for col, prop in enumerate(args.props):
        get, transform = PROPS[prop]
        print(f"\n--- {prop} ({transform}) ---")
        per_grid = []
        for NT in args.NT_list:
            for NR in args.NR_list:
                M, eta_grid, xi_grid, p_grid = sample_grid(state, geom, NT, NR, prop, args.primary)
                ok = np.sum(np.isfinite(M))
                s, sr, on_grid = svd_analyze(M, transform)
                off, n_valid, probes = offgrid_probe(state, geom, args.primary, eta_grid, xi_grid, {prop: M}, args.ranks, prop, n_probes=args.probes)
                # Strip NaN per rank for the recommender
                off_clean = {r: e[np.isfinite(e)] for r, e in off.items()}
                r_rec, max_at_rec = recommend(off_clean, args.budget)
                per_grid.append((NT, NR, s, sr, on_grid, off, off_clean, probes, r_rec, max_at_rec, n_valid))
                print(f"  NT={NT:>3} NR={NR:>3}  cells_ok={ok:>5}/{M.size}  σ_20/σ_1={sr[min(19, len(sr)-1)]:.2e}  rec rank={r_rec} (max={max_at_rec:.2e})")
                for r in args.ranks:
                    e = off_clean[r]
                    rows.append((args.fluid, args.region, args.source, args.primary, prop, NT, NR, r, e.max(), np.quantile(e, 0.99), np.median(e), n_valid))

        # Row 0: rank-vs-err for each (NT, NR) grid
        ax_off = axes[0, col]
        for tup in per_grid:
            NT, NR = tup[0], tup[1]
            off_clean = tup[6]
            maxes = [off_clean[r].max() for r in args.ranks]
            ax_off.semilogy(args.ranks, maxes, marker="o", lw=1.0, label=f"{NT}×{NR}")
        ax_off.axhline(args.budget, color="gray", ls=":", lw=0.8, label=f"budget={args.budget:.0e}")
        ax_off.set_xlabel("rank r")
        ax_off.set_ylabel("off-grid max rel err")
        ax_off.set_title(f"{prop} — off-grid (Hermite eval)")
        ax_off.legend(fontsize=8)
        ax_off.grid(True, alpha=0.3, which="both")

        # Row 1: σ decay for the largest grid
        ax_sig = axes[1, col]
        s_big = per_grid[-1][2]; sr_big = per_grid[-1][3]
        ax_sig.semilogy(np.arange(1, len(sr_big)+1), sr_big, lw=1.2, color="C0")
        ax_sig.axhline(1e-12, color="gray", ls=":", lw=0.5)
        ax_sig.set_xlabel("mode k")
        ax_sig.set_ylabel(r"$\sigma_k / \sigma_1$")
        ax_sig.set_title(f"{prop} — σ decay  ({per_grid[-1][0]}×{per_grid[-1][1]})")
        ax_sig.set_xlim(0, min(50, len(sr_big)))
        ax_sig.set_ylim(1e-16, 2.0)
        ax_sig.grid(True, alpha=0.3, which="both")

        # Row 2: worst-cell (ξ, η) localization at the LARGEST grid + HIGHEST rank
        ax_loc = axes[2, col]
        tup = per_grid[-1]
        NT, NR = tup[0], tup[1]
        off_full = tup[5]
        off_clean = tup[6]
        probes = tup[7]
        r_top = args.ranks[-1]
        e = off_full[r_top]
        ok = np.isfinite(e)
        xis_ok = probes["xi"][ok]; etas_ok = probes["eta"][ok]; e_ok = e[ok]
        log_e = np.log10(np.clip(e_ok, 1e-16, 1.0))
        sc = ax_loc.scatter(xis_ok, etas_ok, c=log_e, s=4, alpha=0.6, cmap="viridis", vmin=-12, vmax=-3)
        thr = np.quantile(e_ok, 0.99)
        worst = ok & (e >= thr)
        ax_loc.scatter(probes["xi"][worst], probes["eta"][worst], facecolors="none", edgecolors="red", s=30, lw=0.8, label=f"top 1% (≥{thr:.1e})")
        thr99 = np.quantile(e_ok, 0.999)
        worst9 = ok & (e >= thr99)
        ax_loc.scatter(probes["xi"][worst9], probes["eta"][worst9], facecolors="magenta", edgecolors="black", s=70, marker="*", lw=0.4, label=f"top 0.1% (≥{thr99:.1e})")
        ax_loc.set_xlabel("ξ (primary)"); ax_loc.set_ylabel("η (secondary)")
        ax_loc.set_xlim(0, 1); ax_loc.set_ylim(0, 1)
        ax_loc.set_title(f"{prop} worst @ rank {r_top}  ({NT}×{NR})")
        ax_loc.legend(fontsize=7, loc="upper right")
        ax_loc.grid(True, alpha=0.3)
        plt.colorbar(sc, ax=ax_loc, fraction=0.04, pad=0.02, label=r"log10|Δ/truth|")

    fig.suptitle(label, fontsize=10)
    plt.tight_layout()
    png = out_dir / f"sizing_{args.fluid}_{args.region}_{args.source}_{args.primary.replace('/','_')}.png"
    fig.savefig(png, dpi=110, bbox_inches="tight")
    print(f"\nPNG: {png}")

    csv = out_dir / f"sizing_{args.fluid}_{args.region}_{args.source}_{args.primary.replace('/','_')}.csv"
    import csv as _csv
    with open(csv, "w", newline="") as f:
        w = _csv.writer(f)
        w.writerow(["fluid", "region", "source", "primary", "prop", "NT", "NR", "rank", "max_err", "p99_err", "p50_err", "n_valid_probes"])
        for r in rows:
            w.writerow(r)
    print(f"CSV: {csv}")

    # Final recommendation table per prop, over best (NT, NR)
    print(f"\nRECOMMENDATIONS  (smallest (NT,NR,rank) meeting budget={args.budget:.0e}):")
    print(f"{'prop':>5} {'best (NT,NR)':>14} {'rank':>5} {'achieved':>12}")
    for prop in args.props:
        # Find smallest grid + rank combo that meets budget
        candidates = [r for r in rows if r[4] == prop and r[8] <= args.budget]
        if not candidates:
            best = min((r for r in rows if r[4] == prop), key=lambda r: r[8])
            print(f"{prop:>5} {f'{best[5]}x{best[6]}':>14} {best[7]:>5} {best[8]:>12.2e} (budget NOT met)")
        else:
            best = min(candidates, key=lambda r: r[5] * r[6] * r[7])  # minimize NT·NR·rank
            print(f"{prop:>5} {f'{best[5]}x{best[6]}':>14} {best[7]:>5} {best[8]:>12.2e}")

    if args.interactive:
        plt.show()


if __name__ == "__main__":
    main()
