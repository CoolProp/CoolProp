#!/usr/bin/env python3
"""
Generate the SRK and Peng-Robinson critical-point constants in
src/Backends/Cubics/CubicConstants.h from their closed forms (Bell and Deiters,
Ind. Eng. Chem. Res., 2021), in 60-digit decimal arithmetic, rounded once to the
nearest double and written as hex-float literals.

#2745 typed these constants as decimal digits and got three of four wrong; never
transcribe them by hand.  Usage:

    python3 dev/scripts/cubic_constants.py          # print the literals
    python3 dev/scripts/cubic_constants.py --check  # verify the header matches

Standard library only.
"""
import argparse
import re
import sys
from decimal import Decimal, getcontext
from pathlib import Path

getcontext().prec = 60
HEADER = Path(__file__).resolve().parents[2] / "src" / "Backends" / "Cubics" / "CubicConstants.h"


def cbrt(x: Decimal) -> Decimal:
    y = Decimal(float(x) ** (1.0 / 3.0))
    for _ in range(10):  # Newton; quadratic convergence from a 16-digit start
        y -= (y * y * y - x) / (3 * y * y)
    return y


def constants() -> dict:
    c2 = cbrt(Decimal(2))
    out = {
        "SRK_Omega_a": 1 / (9 * (c2 - 1)),
        "SRK_Omega_b": (c2 - 1) / 3,
    }
    # Peng-Robinson: eta = b/v_c = 1/(1 + (4 - sqrt 8)^(1/3) + (4 + sqrt 8)^(1/3)).  In units with
    # b = 1 and RT = 1 (so v_c = 1/eta), dp/dv = 0 gives a = g^2/((v - 1)^2 (2v + 2)) with
    # g = v^2 + 2v - 1, then p_c = 1/(v - 1) - a/g, Omega_b = p_c, Omega_a = a p_c, Z_c = p_c v_c.
    s8 = Decimal(8).sqrt()
    v = 1 + cbrt(4 - s8) + cbrt(4 + s8)
    g = v * v + 2 * v - 1
    a = g * g / ((v - 1) ** 2 * (2 * v + 2))
    p = 1 / (v - 1) - a / g
    out.update({"PR_Omega_a": a * p, "PR_Omega_b": p, "PR_Zc": p * v})
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--check", action="store_true", help=f"verify the literals in {HEADER.name}")
    args = parser.parse_args()
    # float(Decimal) rounds correctly to the nearest double
    literals = {name: float(value).hex() for name, value in constants().items()}
    if not args.check:
        for name, value in constants().items():
            print(f"{name:12s} {literals[name]:24s} {str(value)[:40]}")
        return 0
    text = HEADER.read_text()
    bad = 0
    for name, lit in literals.items():
        m = re.search(rf"double {name}\(\)\s*\{{\s*return\s+([0-9a-fA-FxXpP.+-]+);", text)
        if m is None or float.fromhex(m.group(1)) != float.fromhex(lit):
            print(f"{name}: header has {m.group(1) if m else 'nothing'}, expected {lit}")
            bad += 1
    print(f"{HEADER.name}: {'OK' if bad == 0 else f'{bad} mismatch(es)'}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
