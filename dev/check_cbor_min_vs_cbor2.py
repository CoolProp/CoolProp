#!/usr/bin/env python3
"""Pin the vendored CBOR encoder (dev/cbor_min.py) against the reference cbor2.

`generate_headers.py` uses the stdlib-only `cbor_min` to embed the fluid data
as CBOR (no third-party build dependency). This check guarantees `cbor_min`
stays byte-for-byte identical to `cbor2` for the JSON value model, so the
hand-rolled encoder can never silently diverge from the reference. It needs
cbor2 installed, so it runs ONLY in a dedicated CI job on a standard runner —
never on the C++ build path.

Exits non-zero on any divergence.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cbor_min  # noqa: E402

try:
    import cbor2  # noqa: E402
except ImportError:
    sys.exit("cbor2 is required for this check: pip install cbor2")

# Edge cases exercising every cbor_min encoding branch (int length widths,
# negatives, doubles, string/array/map length widths, unicode, nesting, empties).
CASES = [
    None, True, False,
    0, 1, 23, 24, 255, 256, 65535, 65536, 2**32 - 1, 2**32, 2**64 - 1,
    -1, -24, -25, -256, -65536, -(2**32), -(2**64),
    0.0, -0.0, 1.5, -246119.46072729214, 3.141592653589793, 1e308, 5e-324, 2.2250738585072014e-308,
    "", "a", "x" * 23, "x" * 24, "x" * 256, "x" * 70000,
    "café ☃ 𝄞 — unicode",
    [], [1, 2, 3], [None, True, "x", 1.5, [1, [2, [3]]]],
    {}, {"k": "v"}, {"a": 1, "b": [1, 2], "c": {"d": -0.5}},
    {str(i): i for i in range(300)},   # map length > 255 (2-byte length)
    list(range(1000)),                 # array length > 255
]


def check(obj, label):
    a = cbor_min.dumps(obj)
    b = cbor2.dumps(obj)
    if a != b:
        print(f"MISMATCH on {label}: cbor_min={len(a)}B cbor2={len(b)}B")
        for i, (x, y) in enumerate(zip(a, b)):
            if x != y:
                print(f"  first diff at byte {i}: cbor_min=0x{x:02x} cbor2=0x{y:02x}")
                break
        return 1
    return 0


def main():
    fails = 0
    for i, obj in enumerate(CASES):
        fails += check(obj, f"case[{i}] {type(obj).__name__}")
    # Whole edge-case set as one document.
    fails += check(CASES, "all-cases-array")

    # The real embedded fluid data, if it has been generated.
    # The real embedded fluid data is the point of this check — require it.
    # A missing file is a hard failure (not a silent edge-case-only pass): in CI
    # this step runs after the build, which generates all_fluids.json; locally,
    # run dev/generate_headers.py first.
    src = os.path.join(os.path.dirname(__file__), "all_fluids.json")
    if not os.path.exists(src):
        print(f"FAIL: {src} not found — run dev/generate_headers.py first so the "
              "real embedded fluid data is compared, not just synthetic cases.")
        return 1
    with open(src) as f:
        real = json.load(f)
    n = check(real, "all_fluids.json")
    fails += n
    if n == 0:
        count = len(real) if isinstance(real, list) else 1
        print(f"all_fluids.json: {count} fluids compared — byte-identical")

    if fails:
        print(f"\n{fails} divergence(s) between cbor_min and cbor2 — FAIL")
        return 1
    print("cbor_min is byte-identical to cbor2 across all cases — OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
