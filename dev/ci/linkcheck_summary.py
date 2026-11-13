#!/usr/bin/env python3
"""Summarize the docs link check and group its output.txt by type.

Usage: linkcheck_summary.py OUTPUT_TXT [LOG ...]

Prints one count per line for the CI log, then rewrites OUTPUT_TXT in place
with its lines grouped under a header per type (and sorted within each group),
so the uploaded artifact is easy to scan.  The LOG files (the homepage and
Sphinx build output) are only searched for Python warning lines.  Reporting
only: the caller ignores the exit status, and a missing file is not an error.
"""
import os
import sys

# (title, test) in output order; the first test that matches a line wins.
GROUPS = [
    ("Broken", lambda t: "[broken]" in t),
    ("Timed out", lambda t: "[timeout]" in t),
    ("Permanent redirects (301)", lambda t: "[redirected permanently]" in t),
    ("Temporary redirects (302, incl. DOIs)", lambda t: "[redirected" in t),
]
OTHER = "Other"


def read_lines(path):
    if not os.path.isfile(path):
        return None
    with open(path, encoding="utf-8", errors="replace") as f:
        return [ln.rstrip("\n") for ln in f if ln.strip()]


def main(argv):
    if not argv:
        print(__doc__)
        return 2
    output = read_lines(argv[0])
    if output is not None:
        # Skip headers from an earlier run so a repeated run gives the same result.
        output = [ln for ln in output if not (ln.startswith("== ") and ln.endswith(" =="))]
    logs = [read_lines(p) for p in argv[1:]]

    grouped = {title: [] for title, _ in GROUPS}
    grouped[OTHER] = []
    for line in output or []:
        title = next((t for t, test in GROUPS if test(line)), OTHER)
        grouped[title].append(line)

    def warn(name):
        return sum(name in ln for lines in logs for ln in lines or [])

    broken = len(grouped["Broken"]) + len(grouped["Timed out"])
    rows = [
        ("Broken or timed-out links (incl. transients)", broken),
        ("Permanent redirects (301)", len(grouped["Permanent redirects (301)"])),
        ("Temporary redirects incl. DOIs (302)", len(grouped["Temporary redirects (302, incl. DOIs)"])),
        ("UserWarning lines", warn("UserWarning")),
        ("RuntimeWarning lines", warn("RuntimeWarning")),
        ("Total output.txt lines", len(output) if output is not None else 0),
    ]
    width = max(len(name) for name, _ in rows)
    print("::group::Link check summary")
    for name, n in rows:
        print(f"{name:<{width}} : {n}")
    if output is None:
        print("(no output.txt, so the link counts are zero)")
    print("::endgroup::")

    if output is not None:
        with open(argv[0], "w", encoding="utf-8", newline="\n") as f:
            for title in [t for t, _ in GROUPS] + [OTHER]:
                if grouped[title]:
                    f.write(f"== {title}: {len(grouped[title])} ==\n")
                    f.writelines(ln + "\n" for ln in sorted(grouped[title]))
                    f.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
