#!/usr/bin/env python3
"""Fail if a CI smoke-test screenshot looks blank (uniform pixels).

Usage:
    check-screenshot-rendered.py PATH [THRESHOLD]

Default THRESHOLD is 20.0 on the 0-255 grayscale std-dev scale. Calibration
points (gh-2825 / gh-2826):

    blank Windows WebView2 (subtle gradient)        std-dev ≈ 13.9
    rendered macOS GUI (real Tauri+React content)   std-dev ≈ 37.7

A threshold of 20 sits well between them. Adjust per-OS by passing a
different second argument from the workflow.

Exits 0 if rendered, 1 if blank, 2 on usage error.
"""
import sys
from pathlib import Path

try:
    from PIL import Image  # workflow installs Pillow before invoking us
except ImportError as exc:
    print(f"check-screenshot-rendered: Pillow not available ({exc})", file=sys.stderr)
    sys.exit(2)


def std_dev(path: Path) -> float:
    img = Image.open(path).convert("L")  # grayscale
    pixels = list(img.getdata())
    n = len(pixels)
    if n == 0:
        return 0.0
    mean = sum(pixels) / n
    var = sum((p - mean) ** 2 for p in pixels) / n
    return var ** 0.5


def main() -> int:
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print(__doc__, file=sys.stderr)
        return 2
    path = Path(sys.argv[1])
    threshold = float(sys.argv[2]) if len(sys.argv) == 3 else 20.0
    if not path.is_file():
        print(f"check-screenshot-rendered: {path} not found", file=sys.stderr)
        return 2
    s = std_dev(path)
    print(f"{path}: std-dev = {s:.3f}  (threshold {threshold:.3f})")
    if s < threshold:
        print(f"FAIL: screenshot appears blank/uniform — render likely didn't paint")
        return 1
    print("OK: screenshot has plausible variance")
    return 0


if __name__ == "__main__":
    sys.exit(main())
