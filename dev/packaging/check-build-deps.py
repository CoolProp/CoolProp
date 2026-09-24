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

It checks two things that must not drift: the build dependencies above, and
the -DCOOLPROP_* options the three recipes configure CoolProp with, since
adopting an option in two recipes and not the third is the same mistake in a
different place.

Run it with no arguments to check.  Exit status 0 means the recipes agree.
"""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent

CONTROL = ROOT / "dev/packaging/obs/debian.control"
DSC = ROOT / "dev/packaging/obs/coolprop.dsc"
SPEC = ROOT / "dev/packaging/obs/coolprop.spec"
RULES = ROOT / "dev/packaging/obs/debian.rules"
WORKFLOW = ROOT / ".github/workflows/packaging_offline.yml"

# COOLPROP_VENDOR_THIRD_PARTY=OFF, which every distribution build passes, makes
# cmake/dependencies.cmake resolve these two with find_package instead of CPM,
# so both must be installed at build time.  Names differ per ecosystem, and so
# does the syntax of a version bound, hence one row per dependency.
#
# The bound matters as much as the name.  cmake/dependencies.cmake hard-fails
# below Eigen 3.4 and the exported CoolPropConfig.cmake repeats that floor, so
# a recipe naming eigen3 without the bound describes a build that resolves
# happily and then dies in the middle of compiling.
THIRD_PARTY = (
    # debian name, rpm name, the bound every recipe must state (None = any)
    ("libeigen3-dev", "eigen3-devel", ">= 3.4"),
    ("libfmt-dev", "fmt-devel", None),
)
DEBIAN_THIRD_PARTY = {row[0] for row in THIRD_PARTY}


def parse_build_depends(text):
    """Return {package name: version bound} for an RFC-822 Build-Depends field.

    Handles both shapes we use: the one-line form in coolprop.dsc and the
    continuation-line form in debian.control, where following lines start
    with whitespace.  The bound is normalised to "<operator> <version>", or
    None where the entry states none, because dropping it would let a recipe
    lose "(>= 3.4)" without this checker noticing.
    """
    match = re.search(
        r"^Build-Depends:(.*?)(?=^\S)", text, re.MULTILINE | re.DOTALL
    )
    if not match:
        return None
    declared = {}
    for entry in match.group(1).split(","):
        entry = entry.strip()
        if not entry:
            continue
        # The name ends at the first space, bracket or architecture qualifier.
        name = re.split(r"[\s(\[]", entry, maxsplit=1)[0].strip()
        if not name:
            continue
        bound = re.search(r"\(\s*([<>=]+)\s*([^)]+?)\s*\)", entry)
        declared[name] = (
            "{0} {1}".format(bound.group(1), bound.group(2)) if bound else None
        )
    return declared


def parse_build_requires(text):
    """Return {package name: version bound} for a spec's BuildRequires lines.

    Lines inside %if blocks are collected too: this asks what the spec can
    declare, not what one distribution ends up with.  As above, the bound is
    kept rather than discarded.
    """
    declared = {}
    for line in text.splitlines():
        if not line.startswith("BuildRequires:"):
            continue
        entry = line.split(":", 1)[1].strip()
        parts = entry.split(None, 1)
        if not parts:
            continue
        declared[parts[0]] = parts[1].strip() if len(parts) > 1 else None
    return declared


def parse_workflow_apt_packages(text):
    """Return the packages the offline workflow installs with apt-get."""
    match = re.search(r"^\s*sudo apt-get install -y (.+)$", text, re.MULTILINE)
    if not match:
        return None
    return set(match.group(1).split())


def parse_distro_cmake_options(text):
    """Return the -DCOOLPROP_* options of the distribution configure call.

    The three recipes are a YAML workflow, an RPM spec and a makefile, so the
    only structure they share is a configure command continued over several
    lines with trailing backslashes.  This therefore drops comments, rejoins
    those continuations into whole commands, and picks the command that passes
    COOLPROP_REQUIRE_VENDORED_DEPS, which only an offline packaging build does.
    Dropping comments first matters: all three files describe these options in
    prose above the command, and a comment that quoted one would otherwise be
    mistaken for the command itself.

    Raises ValueError rather than returning a partial answer, because an option
    this cannot read is an option it would silently stop comparing.
    """
    commands = []
    current = ""
    for raw in text.splitlines():
        if raw.lstrip().startswith("#"):
            continue
        line = raw.rstrip()
        if line.endswith("\\"):
            current += line[:-1] + " "
            continue
        current += line
        commands.append(current)
        current = ""
    if current:
        commands.append(current)

    matching = [c for c in commands if "COOLPROP_REQUIRE_VENDORED_DEPS" in c]
    if not matching:
        return None
    if len(matching) > 1:
        raise ValueError(
            "found {0} configure commands passing "
            "COOLPROP_REQUIRE_VENDORED_DEPS; expected exactly one".format(
                len(matching)
            )
        )

    command = matching[0]
    # A value may be bare, double quoted or single quoted; a quoted one can
    # contain spaces, as -DCOOLPROP_PC_REQUIRES="eigen3 fmt" will when step 2
    # of GH #3388 lands.
    pairs = re.findall(
        r"""-D(COOLPROP_\w+)=("[^"]*"|'[^']*'|\S+)""", command
    )
    # Every -DCOOLPROP_ in the command must have been read.  Without this, a
    # value shape the pattern above cannot match would simply disappear from
    # the comparison and the checker would pass while the recipes disagree.
    if len(pairs) != command.count("-DCOOLPROP_"):
        raise ValueError(
            "read {0} of {1} -DCOOLPROP_ options; one of them has a value "
            "this parser cannot read".format(
                len(pairs), command.count("-DCOOLPROP_")
            )
        )
    return {name: value.strip("\"'") for name, value in pairs}


def format_declared(declared):
    """Render a {name: bound} mapping as one readable line."""
    return ", ".join(
        name if bound is None else "{0} {1}".format(name, bound)
        for name, bound in sorted(declared.items())
    )


def describe_bound(bound):
    """Render a version bound for a diagnostic, including its absence."""
    return "without a version bound" if bound is None else "with " + bound


def main():
    problems = []

    control_deps = parse_build_depends(CONTROL.read_text(encoding="utf-8"))
    dsc_deps = parse_build_depends(DSC.read_text(encoding="utf-8"))
    spec_reqs = parse_build_requires(SPEC.read_text(encoding="utf-8"))
    workflow_pkgs = parse_workflow_apt_packages(
        WORKFLOW.read_text(encoding="utf-8")
    )

    options = {}
    for label, path in (
        ("packaging_offline.yml", WORKFLOW),
        ("coolprop.spec", SPEC),
        ("debian.rules", RULES),
    ):
        try:
            options[label] = parse_distro_cmake_options(
                path.read_text(encoding="utf-8")
            )
        except ValueError as error:
            problems.append("{0}: {1}".format(label, error))
            options[label] = None
    workflow_opts = options["packaging_offline.yml"]
    spec_opts = options["coolprop.spec"]
    rules_opts = options["debian.rules"]

    # A parse that finds nothing must be an error, never a silent pass.  This
    # is the fail-open that would make every check below vacuously true.
    for label, value in (
        ("debian.control Build-Depends", control_deps),
        ("coolprop.dsc Build-Depends", dsc_deps),
        ("coolprop.spec BuildRequires", spec_reqs or None),
        ("packaging_offline.yml apt-get install", workflow_pkgs),
        ("packaging_offline.yml cmake options", workflow_opts),
        ("coolprop.spec cmake options", spec_opts),
        ("debian.rules cmake options", rules_opts),
    ):
        if not value:
            problems.append(
                "could not find {0}; this checker cannot verify anything "
                "until its parser is fixed".format(label)
            )
    if problems:
        report(problems)
        return 1

    # 1. The two Debian declarations must be identical, bounds included.  They
    #    describe one package built one way; OBS reads the .dsc and dpkg reads
    #    the control, so a bound present in only one of them is a real skew.
    for name in sorted(set(control_deps) | set(dsc_deps)):
        if name not in dsc_deps:
            problems.append(
                "in debian.control but not coolprop.dsc: " + name
            )
        elif name not in control_deps:
            problems.append(
                "in coolprop.dsc but not debian.control: " + name
            )
        elif control_deps[name] != dsc_deps[name]:
            problems.append(
                "{0} is {1} in debian.control but {2} in coolprop.dsc".format(
                    name,
                    describe_bound(control_deps[name]),
                    describe_bound(dsc_deps[name]),
                )
            )

    # 2. Every recipe must declare the third-party packages that
    #    COOLPROP_VENDOR_THIRD_PARTY=OFF needs at build time, with the version
    #    bound CMake enforces.  A name without its bound is not enough.
    for label, declared, column in (
        ("debian.control", control_deps, 0),
        ("coolprop.dsc", dsc_deps, 0),
        ("coolprop.spec", spec_reqs, 1),
    ):
        for row in THIRD_PARTY:
            name = row[column]
            required_bound = row[2]
            if name not in declared:
                problems.append(
                    "{0} builds with COOLPROP_VENDOR_THIRD_PARTY=OFF but does "
                    "not declare {1}".format(label, name)
                )
            elif required_bound is not None and declared[name] != required_bound:
                problems.append(
                    "{0} declares {1} {2}, but cmake/dependencies.cmake "
                    "requires {3}".format(
                        label,
                        name,
                        describe_bound(declared[name]),
                        required_bound,
                    )
                )

    # 3. What CI installs must be declared by the Debian recipes, so that the
    #    workflow cannot pass by installing something no recipe asks for.
    undeclared = sorted(workflow_pkgs - set(control_deps))
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

    # 4. The three recipes must configure CoolProp identically.  Adopting an
    #    option in two of them and not the third is how this checker's own
    #    subject matter went wrong in the first place.
    for label, opts in (("coolprop.spec", spec_opts), ("debian.rules", rules_opts)):
        for name in sorted(set(workflow_opts) | set(opts)):
            in_ci = workflow_opts.get(name)
            in_recipe = opts.get(name)
            if in_ci != in_recipe:
                describe = (
                    "does not set {0}".format(name)
                    if in_recipe is None
                    else "sets {0}={1}".format(name, in_recipe)
                )
                expected = (
                    "packaging_offline.yml does not set it"
                    if in_ci is None
                    else "packaging_offline.yml uses {0}".format(in_ci)
                )
                problems.append(
                    "{0} {1}, but {2}".format(label, describe, expected)
                )

    if problems:
        report(problems)
        return 1

    print("ok: the recipes agree on build dependencies and configure options")
    print("  debian.control / coolprop.dsc: " + format_declared(control_deps))
    print("  coolprop.spec:                 " + format_declared(spec_reqs))
    print("  shared cmake options:          "
          + ", ".join("{0}={1}".format(k, v)
                      for k, v in sorted(workflow_opts.items())))
    return 0


def report(problems):
    print("FAIL: the packaging recipes disagree", file=sys.stderr)
    for problem in problems:
        print("  - " + problem, file=sys.stderr)
    print(
        "\nEvery build dependency has to be declared in debian.control, "
        "coolprop.dsc and\ncoolprop.spec, and installed by "
        "packaging_offline.yml.  Every -DCOOLPROP_ option\nhas to be passed "
        "the same way by packaging_offline.yml, coolprop.spec and\n"
        "debian.rules.  See dev/packaging/README.md.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    sys.exit(main())
