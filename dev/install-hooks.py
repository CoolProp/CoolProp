#!/usr/bin/env python3
"""Point this repo's git hooks at dev/hooks/ so the pre-commit clang-format
check ships in the repo and stays version-controlled.

Run once per clone (works on macOS, Linux, and Windows):
    python dev/install-hooks.py
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def main() -> int:
    repo_root = Path(
        subprocess.run(
            ["git", "rev-parse", "--show-toplevel"],
            check=True, capture_output=True, text=True,
        ).stdout.strip()
    )
    subprocess.run(["git", "config", "core.hooksPath", "dev/hooks"], check=True, cwd=repo_root)
    print("core.hooksPath -> dev/hooks")
    print("Active hooks:")
    for p in sorted((repo_root / "dev" / "hooks").iterdir()):
        if p.is_file():
            print(f"  {p.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
