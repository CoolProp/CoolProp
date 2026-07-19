"""
Add the binary hydrogen mixture models of Beckmueller, Thol, Bell, Lemmon & Span (2021),
"New Equations of State for Binary Hydrogen-Mixtures Containing Methane, Nitrogen, Carbon
Monoxide, and Carbon Dioxide", J. Phys. Chem. Ref. Data (NIST pub 930967), for the pairs
H2 + {CH4, N2, CO, CO2}.

The models use the GERG-2008 mathematical form (composition-dependent reducing functions +
binary-specific departure functions with polynomial + exponential + Gaussian-bell-shaped terms),
which CoolProp already implements as the "Gaussian+Exponential" departure type.  So this is a
pure data addition: new departure functions + updated binary reducing parameters (Fij = 1).

Data source: Supplementary Material, Table S2 (reducing) and Table S3 (departure).  The pure-fluid
critical parameters used by the reducing functions (Table S1) already match CoolProp exactly,
including CoolProp's updated normal-hydrogen reference (Tc = 33.145 K, rhoc = 15.508 mol/dm^3).

Term ordering for CoolProp's Gaussian+Exponential: the "power" terms (polynomial l=0 AND
exponential l>0) come first (count = Npower, evaluated by add_Power which applies exp(-delta^l)),
then the Gaussian-bell-shaped terms (eta/epsilon/beta/gamma).  Table S3 already lists the terms in
that order.  The table columns are (eta, beta, gamma, epsilon); they are mapped BY NAME here.
"""
import json, os, textwrap

HERE = os.path.dirname(__file__)


def _entry_spans(text):
    """Yield (start, end, obj) char spans of each top-level object in a JSON array text."""
    dec = json.JSONDecoder()
    i = text.index("[") + 1
    n = len(text)
    while i < n:
        while i < n and text[i] in " \t\r\n,":
            i += 1
        if i >= n or text[i] == "]":
            break
        obj, j = dec.raw_decode(text, i)
        yield i, j, obj
        i = j

BECK = "Beckmueller-JPCRD-2021"

# --- Table S3: departure functions (power terms first, then Gaussian) -------------------------
# Each: n, t, d, l (0 for polynomial/Gaussian, >0 for exponential), eta/epsilon/beta/gamma
# (nonzero only for the Gaussian terms), and Npower = number of leading power terms.
departures = {
    "Methane-Hydrogen-Beckmueller2021": dict(
        Npower=4,
        n=[1.280, -0.774, -0.914, 0.360, -2.450, 4.462, -0.972, -2.070],
        t=[0.340, 0.420, 0.800, 3.290, 0.050, 0.254, 0.410, 2.630],
        d=[1, 2, 1, 2, 1, 2, 3, 1],
        l=[0, 0, 1, 1, 0, 0, 0, 0],
        eta=[0, 0, 0, 0, 0.83, 0.34, 0.57, 0.44],
        beta=[0, 0, 0, 0, 0.97, 0.20, 0.26, 0.18],
        gamma=[0, 0, 0, 0, 0.77, 1.74, 0.79, 0.40],
        epsilon=[0, 0, 0, 0, 1.36, 1.44, 1.69, 1.53],
    ),
    "Nitrogen-Hydrogen-Beckmueller2021": dict(
        Npower=4,
        n=[-1.659, -0.268, -0.502, 0.323, 2.976, 5.084, -3.866, -5.474],
        t=[0.724, 0.097, 2.953, 3.500, 2.939, 0.694, 2.329, 1.066],
        d=[1, 2, 1, 2, 1, 1, 1, 1],
        l=[1, 1, 2, 2, 0, 0, 0, 0],
        eta=[0, 0, 0, 0, 1.86, 0.03, 1.85, 0.13],
        beta=[0, 0, 0, 0, 0.93, 0.36, 0.99, 0.24],
        gamma=[0, 0, 0, 0, 1.30, 0.80, 1.50, 0.23],
        epsilon=[0, 0, 0, 0, 2.51, 0.39, 2.51, 1.22],
    ),
    "CarbonMonoxide-Hydrogen-Beckmueller2021": dict(
        Npower=2,
        n=[-0.521, -0.387, -2.590, 4.350],
        t=[2.250, 0.473, 0.585, 0.091],
        d=[1, 2, 1, 2],
        l=[0, 0, 0, 0],
        eta=[0, 0, 0.647, 0.344],
        beta=[0, 0, 0.751, 0.660],
        gamma=[0, 0, 1.86, 2.23],
        epsilon=[0, 0, 1.380, 0.773],
    ),
    "CarbonDioxide-Hydrogen-Beckmueller2021": dict(
        Npower=2,
        n=[3.560, -1.036, -4.835, 12.380, -2.650, -3.300],
        t=[1.470, 1.170, 1.950, 0.310, 1.412, 2.280],
        d=[1, 2, 1, 2, 3, 1],
        l=[0, 0, 0, 0, 0, 0],
        eta=[0, 0, 0.580, 0.200, 0.292, 0.120],
        beta=[0, 0, 0.465, 0.820, 0.520, 1.0],
        gamma=[0, 0, 0.17, 2.11, 1.49, 1.73],
        epsilon=[0, 0, 0.52, 0.15, 0.24, 0.15],
    ),
}

# --- Table S2: reducing parameters, given in the paper's (X + H2) order ------------------------
# betaT, gammaT, betaV, gammaV for order (i=X, j=H2).
reducing_XH2 = {
    "Methane": dict(betaT=1.010, gammaT=1.440, betaV=1.086, gammaV=0.804),
    "Nitrogen": dict(betaT=1.027, gammaT=1.240, betaV=0.993, gammaV=0.773),
    "CarbonMonoxide": dict(betaT=1.078, gammaT=1.105, betaV=1.037, gammaV=1.040),
    "CarbonDioxide": dict(betaT=0.979, gammaT=1.961, betaV=1.198, gammaV=0.842),
}

H2 = "1333-74-0"
CAS = {"Methane": "74-82-8", "Nitrogen": "7727-37-9", "CarbonMonoxide": "630-08-0", "CarbonDioxide": "124-38-9"}
FUNC = {k: f"{k}-Hydrogen-Beckmueller2021" for k in reducing_XH2}


def build():
    # sanity: term lengths + Npower consistency
    for name, d in departures.items():
        L = len(d["n"])
        for k in ("t", "d", "l", "eta", "beta", "gamma", "epsilon"):
            assert len(d[k]) == L, f"{name}: {k} length {len(d[k])} != {L}"
        assert 0 <= d["Npower"] <= L
    print("departure term counts OK")

    # ---- departure functions JSON: append 4 new entries in the file's existing style ----
    dep_path = os.path.join(HERE, "mixture_departure_functions.json")
    dep_txt = open(dep_path).read()
    existing = {o["Name"] for _, _, o in _entry_spans(dep_txt)}

    def fmt_arr(a):
        return "[" + ", ".join(("%g" % v if isinstance(v, float) else str(v)) for v in a) + "]"

    new_blocks = []
    for name, d in departures.items():
        if name in existing:
            print("  (already present)", name); continue
        L = len(d["n"])
        block = (
            "{\n"
            f'"Name" : "{name}",\n'
            f'"BibTeX" : "{BECK}",\n'
            '"aliases" : [],\n'
            '"type" : "Gaussian+Exponential",\n'
            f'"Npower" : {d["Npower"]},\n'
            f'"n" : {fmt_arr(d["n"])},\n'
            f'"t" : {fmt_arr(d["t"])},\n'
            f'"d" : {fmt_arr(d["d"])},\n'
            f'"l" : {fmt_arr(d["l"])},\n'
            f'"eta" : {fmt_arr([float(x) for x in d["eta"]])},\n'
            f'"beta" : {fmt_arr([float(x) for x in d["beta"]])},\n'
            f'"gamma" : {fmt_arr([float(x) for x in d["gamma"]])},\n'
            f'"epsilon" : {fmt_arr([float(x) for x in d["epsilon"]])}\n'
            "}"
        )
        new_blocks.append(block)
        print("  + departure", name, f"(Npower={d['Npower']}, {L} terms)")
    if new_blocks:
        assert dep_txt.rstrip().endswith("]")
        head = dep_txt.rstrip()[:-1].rstrip()  # drop trailing ']'
        dep_txt = head + ",\n" + ",\n".join(new_blocks) + "\n]\n"
        open(dep_path, "w").write(dep_txt)

    # ---- binary pairs JSON: replace ONLY the 4 target entries, leave the rest byte-identical ----
    bip_path = os.path.join(HERE, "mixture_binary_pairs.json")
    bip_txt = open(bip_path).read()
    targets = {frozenset({CAS[X], H2}): X for X in reducing_XH2}
    edits = []  # (start, end, new_text)
    for s, e_, obj in _entry_spans(bip_txt):
        key = frozenset({obj.get("CAS1"), obj.get("CAS2")})
        if key not in targets:
            continue
        X = targets[key]
        red = reducing_XH2[X]
        reversed_order = obj["CAS1"] == H2  # CoolProp stores this pair as (H2, X)
        bt, gt, bv, gv = red["betaT"], red["gammaT"], red["betaV"], red["gammaV"]
        if reversed_order:
            bt, bv = 1.0 / bt, 1.0 / bv  # reducing symmetry: beta_ij = 1/beta_ji, gamma unchanged
        obj["F"] = 1.0
        obj["betaT"], obj["gammaT"], obj["betaV"], obj["gammaV"] = bt, gt, bv, gv
        obj["function"] = FUNC[X]
        obj["BibTeX"] = BECK
        new_text = textwrap.indent(json.dumps(obj, indent=2, sort_keys=True), "  ").lstrip()
        edits.append((s, e_, new_text))
        print(f"  ~ binary {obj['Name1']}-{obj['Name2']} (reversed={reversed_order}) -> "
              f"betaT={bt:.9f} gammaT={gt} betaV={bv:.9f} gammaV={gv} function={FUNC[X]}")
    for s, e_, new_text in sorted(edits, reverse=True):
        bip_txt = bip_txt[:s] + new_text + bip_txt[e_:]
    open(bip_path, "w").write(bip_txt)


if __name__ == "__main__":
    build()
    print("done")
