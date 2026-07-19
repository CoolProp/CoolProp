#!/usr/bin/env python3
"""
pxy_teqp.py -- robust binary VLE envelopes (P-x-y at fixed T, or T-x-y at fixed P)
using teqp's continuation tracer, driven by CoolProp-format JSON binary-interaction
parameters (mixture_binary_pairs.json + mixture_departure_functions.json).

Why teqp continuation rather than point-wise QT flashing:
  A point-wise Q=0 flash swept over composition falls onto spurious roots for
  wide-boiling / near-critical pairs -- trivial roots (x==y at the wrong pressure)
  and "beyond-critical" roots (imposing a bubble point on a feed past the mixture
  critical composition, which has no bubble point).  teqp follows the phase boundary
  as a single connected locus and terminates cleanly at the mixture critical, so it
  never produces those artifacts.  See the H2-CO2 260 K case, where QT flashing
  reported a bogus ~340 MPa critical while the true model critical is ~67 MPa.

The multifluid model reads the SAME JSON files CoolProp uses, so you can point it at
a CoolProp source tree (its dev/fluids + dev/mixtures) and, optionally, at modified
BIP/departure files to test new binary models before they are baked into a build.

Requires: teqp (`pip install teqp`).  CoolProp Python is used when available for exact
pure-fluid saturation endpoints; a teqp-only fallback is provided.

Examples
--------
# H2-N2 at 70.4 K using a CoolProp tree that already has the new binary parameters:
python pxy_teqp.py --pair "Nitrogen&Hydrogen" --T 70.4 \
    --root /path/to/CoolProp --logp --out n2h2.png --csv n2h2.csv

# Override the binary parameters explicitly (e.g. an un-built new model):
python pxy_teqp.py --pair "CarbonDioxide&Hydrogen" --T 260 \
    --root /path/to/CoolProp \
    --bip  /path/to/CoolProp/dev/mixtures/mixture_binary_pairs.json \
    --dep  /path/to/CoolProp/dev/mixtures/mixture_departure_functions.json \
    --logp --out h2co2.png

# Isobaric T-x-y at 5 MPa:
python pxy_teqp.py --pair "Nitrogen&Oxygen" --P 5e6 --root /path/to/CoolProp --out n2o2_Txy.png
"""
import argparse
import os
import numpy as np


# ----------------------------------------------------------------------------- model
def build_binary_model(pair, root, bip="", dep=""):
    """Build a teqp multifluid model for a binary pair from CoolProp-format JSON.

    pair : "A&B" or "A,B" (CoolProp fluid names).
    root : CoolProp source tree containing dev/fluids/*.json.
    bip  : optional path to a mixture_binary_pairs.json to use instead of the default.
    dep  : optional path to a mixture_departure_functions.json to use instead.
    """
    import teqp

    comps = [c.strip() for c in pair.replace("&", ",").split(",") if c.strip()]
    if len(comps) != 2:
        raise ValueError(f"expected exactly two components, got {comps!r}")
    model = teqp.build_multifluid_model(comps, root, bip or "", {}, dep or "")
    return model, comps


def _pure_saturation(fluid, T, root):
    """Return (rhoL, rhoV) [mol/m^3] for a pure fluid at temperature T.

    Prefers CoolProp (exact); falls back to teqp.pure_VLE_T with critical-density-
    scaled guesses so the helper still works in a teqp-only environment.
    """
    try:
        import CoolProp.CoolProp as CP

        rhoL = CP.PropsSI("Dmolar", "T", T, "Q", 0, fluid)
        rhoV = CP.PropsSI("Dmolar", "T", T, "Q", 1, fluid)
        if rhoL > rhoV > 0:
            return float(rhoL), float(rhoV)
    except Exception:
        pass

    import teqp

    mp = teqp.build_multifluid_model([fluid], root)
    rhoc = 1.0 / mp.get_vcvec()[0]
    for fL in (2.8, 2.5, 3.1, 2.2, 2.0):
        for fV in (1e-2, 1e-3, 5e-2, 1e-1, 3e-1):
            try:
                r = mp.pure_VLE_T(T, fL * rhoc, fV * rhoc, 30)
                rhoL, rhoV = float(r[0]), float(r[1])
                if rhoL > rhoV > 0:
                    return rhoL, rhoV
            except Exception:
                continue
    raise RuntimeError(f"could not converge pure saturation for {fluid} at {T} K")


def _start_component(model, comps, T):
    """Pick a subcritical pure component (Tc > T) to start the isotherm trace from.

    Prefer the most condensable (highest Tc); its saturation endpoint is the most
    robust place to begin.  Raise if both components are supercritical at T (no VLE).
    """
    Tc = list(model.get_Tcvec())
    candidates = [i for i in (0, 1) if Tc[i] > T]
    if not candidates:
        raise RuntimeError(
            f"both components supercritical at {T} K (Tc = {Tc}); no VLE to trace"
        )
    return max(candidates, key=lambda i: Tc[i])


# ----------------------------------------------------------------------------- traces
def _extract(J, comps, key_prefix):
    """Turn a teqp trace (list of dicts) into arrays keyed by component name."""
    x0 = np.array([p[f"xL_0 / mole frac."] for p in J])
    y0 = np.array([p[f"xV_0 / mole frac."] for p in J])
    P = np.array([p["pL / Pa"] for p in J])
    T = np.array([p["T / K"] for p in J])
    return {
        "T": T,
        "P": P,
        "x": {comps[0]: x0, comps[1]: 1.0 - x0},  # liquid mole fractions
        "y": {comps[0]: y0, comps[1]: 1.0 - y0},  # vapour mole fractions
        "comps": comps,
    }


def _critical(res):
    """Locate the mixture critical: an *interior* point where liquid and vapour
    compositions merge (x==y with 0<x<1).  The pure trace endpoints trivially have
    x==y and must be excluded.  If no interior merge exists (e.g. an envelope that
    spans pure-to-pure with both components subcritical), report closed=False."""
    c0 = res["comps"][0]
    x, y = res["x"][c0], res["y"][c0]
    d = np.abs(x - y)
    interior = (x > 5e-3) & (x < 1 - 5e-3) & (y > 5e-3) & (y < 1 - 5e-3)
    idxs = np.where(interior)[0]
    if len(idxs) == 0:
        return {"index": None, "closed": False, "T": None, "P": None, "x": None}
    i = int(idxs[np.argmin(d[idxs])])
    return {
        "index": i,
        "T": float(res["T"][i]),
        "P": float(res["P"][i]),
        "x": {k: float(v[i]) for k, v in res["x"].items()},
        "closed": bool(d[i] < 5e-3),
    }


def trace_pxy(pair, T, root, bip="", dep=""):
    """Trace an isothermal P-x-y envelope.  Returns a result dict (see _extract)."""
    model, comps = build_binary_model(pair, root, bip, dep)
    idx = _start_component(model, comps, T)
    rhoL, rhoV = _pure_saturation(comps[idx], T, root)
    rhovecL = np.zeros(2)
    rhovecV = np.zeros(2)
    rhovecL[idx] = rhoL
    rhovecV[idx] = rhoV
    J = model.trace_VLE_isotherm_binary(float(T), rhovecL, rhovecV)
    res = _extract(J, comps, "T")
    res["mode"] = "pxy"
    res["fixed"] = ("T", float(T))
    res["critical"] = _critical(res)
    res["start"] = comps[idx]
    return res


def trace_txy(pair, P, root, bip="", dep=""):
    """Trace an isobaric T-x-y envelope.  Returns a result dict (see _extract)."""
    model, comps = build_binary_model(pair, root, bip, dep)
    # start from a component that has a saturation temperature at this pressure
    Tc = list(model.get_Tcvec())
    start = None
    for idx in sorted((0, 1), key=lambda i: -Tc[i]):
        try:
            import CoolProp.CoolProp as CP

            Tsat = CP.PropsSI("T", "P", P, "Q", 0, comps[idx])
        except Exception:
            # crude fallback: just below Tc
            Tsat = 0.98 * Tc[idx]
        try:
            rhoL, rhoV = _pure_saturation(comps[idx], Tsat, root)
        except Exception:
            continue
        rhovecL = np.zeros(2)
        rhovecV = np.zeros(2)
        rhovecL[idx] = rhoL
        rhovecV[idx] = rhoV
        try:
            J = model.trace_VLE_isobar_binary(float(P), float(Tsat), rhovecL, rhovecV)
            start = comps[idx]
            break
        except Exception:
            continue
    if start is None:
        raise RuntimeError(f"could not start an isobar trace for {pair} at {P} Pa")
    res = _extract(J, comps, "P")
    res["mode"] = "txy"
    res["fixed"] = ("P", float(P))
    res["critical"] = _critical(res)
    res["start"] = start
    return res


# ------------------------------------------------------------------------------ plot
def plot_envelope(res, axis_component=None, logy=False, title=None, ax=None):
    """Plot a bubble/dew envelope.  axis_component chooses the x-axis composition
    (defaults to the second component, typically the more volatile one)."""
    import matplotlib.pyplot as plt

    comps = res["comps"]
    comp = axis_component or comps[1]
    xb = res["x"][comp]
    yb = res["y"][comp]
    if res["mode"] == "pxy":
        vb = res["P"] / 1e6
        ylabel = "pressure  $P$ / MPa"
        fixed = f"{res['fixed'][1]:g} K"
    else:
        vb = res["T"]
        ylabel = "temperature  $T$ / K"
        fixed = f"{res['fixed'][1] / 1e6:g} MPa"

    if ax is None:
        _, ax = plt.subplots(figsize=(7.5, 6.2))
    ax.plot(xb, vb, "-", color="#1f77b4", lw=2.4, label="bubble (liquid)")
    ax.plot(yb, vb, "--", color="#d62728", lw=2.2, label="dew (vapour)")
    cr = res["critical"]
    if cr["closed"]:
        cval = cr["P"] / 1e6 if res["mode"] == "pxy" else cr["T"]
        ax.plot(cr["x"][comp], cval, "k*", ms=13,
                label=f"critical ({cr['x'][comp]:.3f}, "
                      f"{cval:.3g} {'MPa' if res['mode']=='pxy' else 'K'})")
    if logy:
        ax.set_yscale("log")
    ax.set_xlim(0, 1)
    ax.grid(alpha=0.3, which="both")
    ax.set_xlabel(f"{comp} mole fraction  $x, y$", fontsize=12)
    ax.set_ylabel(ylabel + ("  (log)" if logy else ""), fontsize=12)
    ax.legend(fontsize=9, loc="best")
    ax.set_title(title or f"{comps[0]}–{comps[1]} "
                 f"{'P' if res['mode']=='pxy' else 'T'}–x–y at {fixed}", fontsize=12)
    return ax


def write_csv(res, path):
    comps = res["comps"]
    with open(path, "w") as f:
        f.write(f"T_K,P_Pa,xL_{comps[0]},xL_{comps[1]},yV_{comps[0]},yV_{comps[1]}\n")
        for i in range(len(res["P"])):
            f.write(f"{res['T'][i]:.6f},{res['P'][i]:.6f},"
                    f"{res['x'][comps[0]][i]:.6f},{res['x'][comps[1]][i]:.6f},"
                    f"{res['y'][comps[0]][i]:.6f},{res['y'][comps[1]][i]:.6f}\n")


# ------------------------------------------------------------------------------- CLI
def main():
    ap = argparse.ArgumentParser(description="Binary VLE envelope via teqp + CoolProp JSON")
    ap.add_argument("--pair", required=True, help='e.g. "CarbonDioxide&Hydrogen"')
    ap.add_argument("--root", required=True, help="CoolProp source tree (has dev/fluids)")
    ap.add_argument("--T", type=float, help="isotherm temperature [K] -> P-x-y")
    ap.add_argument("--P", type=float, help="isobar pressure [Pa] -> T-x-y")
    ap.add_argument("--bip", default="", help="override mixture_binary_pairs.json path")
    ap.add_argument("--dep", default="", help="override mixture_departure_functions.json path")
    ap.add_argument("--axis", default=None, help="component for the x-axis (default: 2nd)")
    ap.add_argument("--logp", action="store_true", help="log pressure/temperature axis")
    ap.add_argument("--out", default=None, help="output PNG path")
    ap.add_argument("--csv", default=None, help="output CSV path")
    args = ap.parse_args()

    if (args.T is None) == (args.P is None):
        ap.error("specify exactly one of --T (P-x-y) or --P (T-x-y)")

    if args.T is not None:
        res = trace_pxy(args.pair, args.T, args.root, args.bip, args.dep)
    else:
        res = trace_txy(args.pair, args.P, args.root, args.bip, args.dep)

    cr = res["critical"]
    npts = len(res["P"])
    print(f"traced {npts} points; start pure = {res['start']}")
    if cr["closed"]:
        if res["mode"] == "pxy":
            print(f"mixture critical: P = {cr['P']/1e6:.2f} MPa, "
                  f"x = {cr['x']}")
        else:
            print(f"mixture critical: T = {cr['T']:.2f} K, x = {cr['x']}")
    else:
        print("no critical closure in traced range (both pures subcritical, "
              "or trace terminated at a composition bound)")

    if args.csv:
        write_csv(res, args.csv)
        print(f"wrote {args.csv}")
    if args.out:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        plot_envelope(res, axis_component=args.axis, logy=args.logp)
        plt.tight_layout()
        plt.savefig(args.out, dpi=145)
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
