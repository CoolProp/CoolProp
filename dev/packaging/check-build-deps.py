#!/usr/bin/env python3
"""Check that the packaging recipes agree on their build dependencies.

The same build dependencies are spelled out in four places, because each
build system insists on reading its own file:

    dev/packaging/obs/debian.control    Build-Depends, for the Debian package
    dev/packaging/obs/coolprop.dsc      Build-Depends, what OBS and sbuild read
    dev/packaging/obs/coolprop.spec     BuildRequires, for the RPM package
    .github/workflows/packaging_offline.yml   the apt-get line in CI

Nothing made those four agree, and that is not hypothetical: the commit that
added Eigen and fmt to three of them missed coolprop.dsc, which is precisely
the file the Debian OBS build reads, so the bug it was fixing stayed live on
Debian and CI still went green.  This script is the missing coupling.  It is
run by the packaging workflow, so dropping a dependency from one recipe now
fails CI instead of failing somebody's distribution build months later.

Run it with no arguments to check.  Exit status 0 means the recipes agree.
"""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent

CONTROL = ROOT / "dev/packaging/obs/debian.control"
DSC = ROOT / "dev/packaging/obs/coolprop.dsc"
SPEC = ROOT / "dev/packaging/obs/coolprop.spec"
WORKFLOW = ROOT / ".github/workflows/packaging_offline.yml"

# COOLPROP_VENDOR_THIRD_PARTY=OFF, which every distribution build passes, makes
# cmake/dependencies.cmake resolve these two with find_package instead of CPM,
# so both must be installed at build time.  Names differ per ecosystem.
DEBIAN_THIRD_PARTY = {"libeigen3-dev", "libfmt-dev"}
RPM_THIRD_PARTY = {"eigen3-devel", "fmt-devel"}


def parse_build_depends(text):
    """Return the package names in an RFC-822 Build-Depends field.

    Handles both shapes we use: the one-line form in coolprop.dsc and the
    continuation-line form in debian.control, where following lines start
    with whitespace.  Version constraints such as "(>= 3.4)" are stripped;
    this compares which packages are declared, not their bounds.
    """
    match = re.search(
        r"^Build-Depends:(.*?)(?=^\S)", text, re.MULTILINE | re.DOTALL
    )
    if not match:
        return None
    names = set()
    for entry in match.group(1).split(","):
        entry = entry.strip()
        if not entry:
            continue
        # Drop a version constraint and any architecture qualifier.
        name = re.split(r"[\s(\[]", entry, maxsplit=1)[0].strip()
        if name:
            names.add(name)
    return names


def parse_build_requires(text):
    """Return the package names on BuildRequires lines in an RPM spec.

    Lines inside %if blocks are collected too: this asks what the spec can
    declare, not what one distribution ends up with.
    """
    names = set()
    for line in text.splitlines():
        if not line.startswith("BuildRequires:"):
            continue
        entry = line.split(":", 1)[1].strip()
        # "eigen3-devel >= 3.4" and "cmake >= 3.14" both reduce to the name.
        names.add(entry.split()[0])
    return names


def parse_workflow_apt_packages(text):
    """Return the packages the offline workflow installs with apt-get."""
    match = re.search(r"^\s*sudo apt-get install -y (.+)$", text, re.MULTILINE)
    if not match:
        return None
    return set(match.group(1).split())


def main():
    problems = []

    control_deps = parse_build_depends(CONTROL.read_text(encoding="utf-8"))
    dsc_deps = parse_build_depends(DSC.read_text(encoding="utf-8"))
    spec_reqs = parse_build_requires(SPEC.read_text(encoding="utf-8"))
    workflow_pkgs = parse_workflow_apt_packages(
        WORKFLOW.read_text(encoding="utf-8")
    )

    # A parse that finds nothing must be an error, never a silent pass.  This
    # is the fail-open that would make every check below vacuously true.
    for label, value in (
        ("debian.control Build-Depends", control_deps),
        ("coolprop.dsc Build-Depends", dsc_deps),
        ("coolprop.spec BuildRequires", spec_reqs or None),
        ("packaging_offline.yml apt-get install", workflow_pkgs),
    ):
        if not value:
            problems.append(
                "could not find {0}; this checker cannot verify anything "
                "until its parser is fixed".format(label)
            )
    if problems:
        report(problems)
        return 1

    # 1. The two Debian declarations must be identical.  They describe one
    #    package built one way; OBS reads the .dsc and dpkg reads the control.
    if control_deps != dsc_deps:
        only_control = sorted(control_deps - dsc_deps)
        only_dsc = sorted(dsc_deps - control_deps)
        if only_control:
            problems.append(
                "in debian.control but not coolprop.dsc: "
                + ", ".join(only_control)
            )
        if only_dsc:
            problems.append(
                "in coolprop.dsc but not debian.control: " + ", ".join(only_dsc)
            )

    # 2. Every recipe must declare the third-party packages that
    #    COOLPROP_VENDOR_THIRD_PARTY=OFF needs at build time.
    for label, declared, required in (
        ("debian.control", control_deps, DEBIAN_THIRD_PARTY),
        ("coolprop.dsc", dsc_deps, DEBIAN_THIRD_PARTY),
        ("coolprop.spec", spec_reqs, RPM_THIRD_PARTY),
    ):
        missing = sorted(required - declared)
        if missing:
            problems.append(
                "{0} builds with COOLPROP_VENDOR_THIRD_PARTY=OFF but does not "
                "declare: {1}".format(label, ", ".join(missing))
            )

    # 3. What CI installs must be declared by the Debian recipes, so that the
    #    workflow cannot pass by installing something no recipe asks for.
    undeclared = sorted(workflow_pkgs - control_deps)
    if undeclared:
        problems.append(
            "packaging_offline.yml installs packages no recipe declares: "
            + ", ".join(undeclared)
        )
    uninstalled = sorted(DEBIAN_THIRD_PARTY - workflow_pkgs)
    if uninstalled:
        problems.append(
            "packaging_offline.yml does not install: " + ", ".join(uninstalled)
        )

    if problems:
        report(problems)
        return 1

    print("ok: the packaging recipes agree on their build dependencies")
    print("  debian.control / coolprop.dsc: " + ", ".join(sorted(control_deps)))
    print("  coolprop.spec:                 " + ", ".join(sorted(spec_reqs)))
    return 0


def report(problems):
    print("FAIL: the packaging recipes disagree", file=sys.stderr)
    for problem in problems:
        print("  - " + problem, file=sys.stderr)
    print(
        "\nEvery build dependency has to be declared in debian.control, "
        "coolprop.dsc\nand coolprop.spec, and installed by "
        "packaging_offline.yml.  See dev/packaging/README.md.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    sys.exit(main())
