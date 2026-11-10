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
DEPENDENCIES = ROOT / "cmake/dependencies.cmake"
CMAKELISTS = ROOT / "CMakeLists.txt"
GITATTRIBUTES = ROOT / ".gitattributes"
CHANGELOG = ROOT / "dev/packaging/obs/debian.changelog"
PYPROJECT = ROOT / "pyproject.toml"
SERVICE = ROOT / "dev/packaging/obs/_service"
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
    # debian name, rpm name, whether the Eigen floor below applies to it
    ("libeigen3-dev", "eigen3-devel", True),
    ("libfmt-dev", "fmt-devel", False),
)
DEBIAN_THIRD_PARTY = {row[0] for row in THIRD_PARTY}


def parse_upstream_version(text):
    """Return the version make-release-tarball.sh will name the archive with.

    CMakeLists.txt is the single source of truth: the four components there
    decide the tarball name, so every recipe has to follow it rather than
    restate it.  Returns None if any component cannot be read, because a
    half-parsed version would be compared against and silently agree.
    """
    # Comments go first and every pattern is anchored to the start of a line,
    # because make-release-tarball.sh parses the same file that way.  An
    # unanchored search reads "#set(COOLPROP_VERSION_PATCH 1)" sitting above a
    # live "set(COOLPROP_VERSION_PATCH 2)", so the two tools would answer
    # differently about the same CMakeLists and this whole check would be void:
    # the recipes would stay at the old version, CI would pass, and OBS would
    # get a tarball whose directory %autosetup cannot find.
    live = "\n".join(
        line for line in text.splitlines() if not line.lstrip().startswith("#")
    )

    parts = {}
    for name in ("MAJOR", "MINOR", "PATCH"):
        match = re.search(
            r"^[ \t]*set[ \t]*\([ \t]*COOLPROP_VERSION_{0}[ \t]+"
            r"([0-9]+)[ \t]*\)".format(name),
            live,
            re.M,
        )
        if not match:
            return None
        parts[name] = match.group(1)
    # REVISION is "dev" between releases and empty on a release tag, so an
    # empty match is a legitimate answer here and only a missing LINE is not.
    revision = re.search(
        r"^[ \t]*set[ \t]*\([ \t]*COOLPROP_VERSION_REVISION[ \t]*"
        r"([A-Za-z0-9]*)[ \t]*\)",
        live,
        re.M,
    )
    if revision is None:
        return None
    return "{0}.{1}.{2}".format(
        parts["MAJOR"], parts["MINOR"], parts["PATCH"]
    ), revision.group(1)


def distro_version(numeric, revision):
    """Render the version a distribution package should carry.

    A pre-release has to sort BELOW the release it precedes.  Both RPM and dpkg
    spell that with a tilde, and neither treats a bare suffix that way: dpkg
    puts 8.0.1dev ABOVE 8.0.1, so a snapshot named that way would shadow the
    release and block the upgrade to it.
    """
    return numeric if not revision else "{0}~{1}".format(numeric, revision)


def parse_spec_versions(text):
    """Return (upstream_version, Version) from the RPM spec."""
    upstream = re.search(r"^%global\s+upstream_version\s+(\S+)", text, re.M)
    version = re.search(r"^Version:\s*(\S+)", text, re.M)
    return (
        upstream.group(1) if upstream else None,
        version.group(1) if version else None,
    )


def parse_dsc_version(text):
    """Return (Version, tarball name) from the Debian source control file.

    The tarball is read from the Files: stanza specifically, not from the first
    indented line that happens to look like "token number token".  dpkg writes
    Checksums-Sha1: and Checksums-Sha256: BEFORE Files:, and their entries have
    that same shape, so a looser search reads a checksum line and never looks
    at the name this check exists to verify.  Exactly one entry is required:
    a second one would otherwise go unexamined.
    """
    version = re.search(r"^Version:\s*(\S+)", text, re.M)

    tarball = None
    lines = text.splitlines()
    for index, line in enumerate(lines):
        if line.rstrip() != "Files:":
            continue
        entries = []
        for entry in lines[index + 1:]:
            if not entry[:1].isspace():
                break
            fields = entry.split()
            if len(fields) >= 3:
                entries.append(fields[2])
        # One source tarball, or this parser cannot say which one matters.
        tarball = entries[0] if len(entries) == 1 else None
        break

    return (
        version.group(1) if version else None,
        tarball,
    )


def parse_soname_major(spec_text, dsc_text, control_text):
    """Return the shared-library major each Debian and RPM file hardcodes.

    coolprop.spec defines it as %define sover, and the Debian package name
    carries it in the binary package name (libcoolprop8).  All of them have to
    track COOLPROP_VERSION_MAJOR, which cmake/CoolPropLibrary.cmake uses as the
    target SOVERSION, or %files claims a libCoolProp.so.N that was never built.
    """
    sover = re.search(r"^%define\s+sover\s+([0-9]+)", spec_text, re.M)
    binaries = re.findall(r"libcoolprop([0-9]+)", dsc_text + control_text)
    return (
        sover.group(1) if sover else None,
        sorted(set(binaries)),
    )


def parse_changelog_version(text):
    """Return the version of the newest debian.changelog entry."""
    match = re.match(r"^\S+\s+\(([^)]+)\)", text)
    return match.group(1) if match else None


def packaging_is_lf_only(text):
    """Is dev/packaging declared eol=lf in .gitattributes?

    The repository sets "* text=auto", which converts text files to the
    checkout platform's native line ending.  On Windows that gives
    coolprop.spec CRLF, and rpmbuild then writes the %prep body into a shell
    script where the stray CR is executed as a command:

        /var/tmp/rpm-tmp.XXXX: line 46: $'\r': command not found

    That failed a real OBS build.  CI cannot detect it, because a Linux
    checkout normalises to LF either way, so what is checkable is that the
    rule exists at all.
    """
    return bool(
        re.search(
            r"^dev/packaging/\*\*\s+text\s+eol=lf\s*$", text, re.M
        )
    )


def lto_is_disabled(spec_text, rules_text):
    """Do both recipes switch link-time optimisation off?

    CoolProp embeds dev/all_fluids.cbor with incbin, which emits a top-level
    __asm__ holding a .incbin directive.  The assembler finds that file
    through the -I paths of the compile step.  Under LTO the asm is streamed
    into the LTO objects and re-assembled at link time by lto-wrapper, which
    runs from /tmp without those -I paths, so the link dies with

        /tmp/ccXXXXXX.s:46: Error: file not found: all_fluids.cbor

    That failed real OBS builds on Tumbleweed x86_64 and i586, because
    openSUSE and Fedora both put -flto=auto in %optflags.  Debian does not
    enable LTO by default, so its guard is precautionary, but both are checked
    so that neither can be dropped without the other being reconsidered.

    Returns a list of the recipes that are missing their guard.
    """
    missing = []
    if not re.search(r"^%define\s+_lto_cflags\s+%\{nil\}\s*$", spec_text, re.M):
        missing.append("coolprop.spec (%define _lto_cflags %{nil})")
    if not re.search(
        r"^export\s+DEB_BUILD_MAINT_OPTIONS\s*=.*\boptimize=-lto\b",
        rules_text,
        re.M,
    ):
        missing.append("debian.rules (optimize=-lto)")
    return missing


# The four places dev/packaging/obs/_service names a version.  All of them
# have to say the same thing, and that thing is the RELEASE this development
# series becomes, not the dev snapshot: a snapshot is never published as a
# GitHub release, so there is nothing for download_url to fetch.
# dev/generate_headers.py runs during the build and annotates with builtin
# generics (list[Path] at dev/generate_headers.py:364), which is a syntax error
# before Python 3.9.  This constant is the requirement, stated once.
#
# It is deliberately NOT read from pyproject.toml.  Every floor below used to
# be compared against that file, which meant lowering that one line satisfied
# the whole check while generate_headers.py still needed 3.9.
REQUIRED_PYTHON = (3, 9)

SERVICE_VERSION_SITES = (
    "the release tag in the download path",
    "the tarball name in the download path",
    "the download_url filename",
    "the verify_file file",
)


def parse_service_versions(text):
    """Return {where: version} for every version _service names.

    _service is the bootstrap path a packager runs by hand with
    "osc service manualrun": it downloads a published release tarball from
    GitHub and verifies its sha256.  Nothing else in this script looked at
    it, so its version could sit at an old release while CMakeLists.txt moved
    on, and the first symptom would be a packager fetching the wrong tarball
    or a 404.

    A site that cannot be read is simply absent from the returned dict, and
    the caller treats that as a failure.  Quietly returning fewer sites than
    there are would make this check pass by finding nothing, which is the
    fail-open this whole script exists to prevent.
    """
    found = {}

    # Drop XML comments first.  parse_upstream_version strips "#" lines from
    # CMakeLists for the same reason: re.search takes the FIRST hit, so a
    # commented-out old block above the live one is read instead of it, and
    # the check then reports on a version that is not in effect.
    text = re.sub(r"<!--.*?-->", "", text, flags=re.S)

    path = re.search(r'<param\s+name="path">([^<]*)</param>', text)
    if path:
        target = path.group(1).strip()
        tag = re.search(r"/releases/download/v([0-9][^/]*)/", target)
        if tag:
            found[SERVICE_VERSION_SITES[0]] = tag.group(1)
        name = re.search(r"/coolprop-([0-9][^/]*?)\.tar\.gz$", target)
        if name:
            found[SERVICE_VERSION_SITES[1]] = name.group(1)

    for site, param in (
        (SERVICE_VERSION_SITES[2], "filename"),
        (SERVICE_VERSION_SITES[3], "file"),
    ):
        match = re.search(
            r'<param\s+name="{0}">\s*coolprop-([0-9].*?)\.tar\.gz\s*</param>'.format(
                param
            ),
            text,
        )
        if match:
            found[site] = match.group(1)

    return found


def parse_python_floors(
    pyproject_text, cmakelists_text, spec_text, control_text
):
    """Return the minimum Python each file demands, as (major, minor) tuples.

    The floor is written in four places: requires-python in pyproject.toml, the
    find_package(Python ...) call that makes CMake skip an older interpreter,
    the spec's BuildRequires and Debian's Build-Depends so the build root
    actually has one.  Leap 15.x ships 3.6 and did fail on exactly this.

    A python requirement that states NO floor counts as (0, 0) rather than
    being skipped.  Skipping it was a fail-open: dropping ">= 3.9" from a line,
    or reverting "python311" to "python3-base", made that line invisible and
    min() over the remaining branches still answered 3.9, so the gate passed
    while the build root could be too old.  That is the regression this is
    here to catch, not a wrong version number.

    A file whose requirement cannot be found at all is absent from the result,
    and the caller treats that as a failure.
    """
    found = {}

    match = re.search(
        r'^requires-python\s*=\s*"[><=~^ ]*([0-9]+)\.([0-9]+)',
        pyproject_text,
        re.M,
    )
    if match:
        found["pyproject.toml requires-python"] = (
            int(match.group(1)),
            int(match.group(2)),
        )

    match = re.search(
        r"^[ \t]*find_package\([ \t]*Python[ \t]+([0-9]+)\.([0-9]+)",
        cmakelists_text,
        re.M,
    )
    if match:
        found["the find_package(Python ...) floor in CMakeLists.txt"] = (
            int(match.group(1)),
            int(match.group(2)),
        )

    # The spec asks for Python once per distribution family.  The weakest
    # branch decides whether some build root ends up too old, so every branch
    # is read and the lowest wins.
    spec_floors = []
    for line in spec_text.splitlines():
        match = re.match(
            r"^BuildRequires:[ \t]+(python[0-9][-A-Za-z0-9_]*)[ \t]*(.*)$", line
        )
        if not match:
            continue
        name, rest = match.group(1), match.group(2).strip()
        bound = re.match(r"^>=[ \t]*([0-9]+)\.([0-9]+)", rest)
        if bound:
            spec_floors.append((int(bound.group(1)), int(bound.group(2))))
            continue
        # A versioned package name such as python311 carries its own floor.
        in_name = re.match(r"^python([0-9])([0-9]+)(?:-[A-Za-z0-9_]+)?$", name)
        if in_name:
            spec_floors.append((int(in_name.group(1)), int(in_name.group(2))))
            continue
        # Neither a bound nor a version in the name: this branch promises
        # nothing, which is weaker than any number.
        spec_floors.append((0, 0))
    if spec_floors:
        found["the weakest python BuildRequires in coolprop.spec"] = min(
            spec_floors
        )

    # Debian's floor matters for the same reason: Ubuntu 20.04 still ships
    # Python 3.8.  Read it from the parsed Build-Depends rather than by
    # searching the file, so a bound mentioned in a comment or in some binary
    # package's Depends: cannot be mistaken for the real one.
    declared = parse_build_depends(control_text)
    if declared:
        control_floors = []
        for name, bound in declared.items():
            if not re.match(r"^python[0-9]", name):
                continue
            version = re.match(r"^>=\s*([0-9]+)\.([0-9]+)", bound or "")
            control_floors.append(
                (int(version.group(1)), int(version.group(2)))
                if version
                else (0, 0)
            )
        if control_floors:
            found["the python3 Build-Depends in debian.control"] = min(
                control_floors
            )

    return found


def parse_eigen_floor(text):
    """Return the Eigen version cmake/dependencies.cmake refuses to go below.

    Read rather than restated, because a floor written down here as well would
    be one more copy to drift: raise it in the CMake and a hardcoded constant
    would keep passing three recipes that all understate it, which is the exact
    failure this script exists to prevent.
    """
    match = re.search(
        r"Eigen3_VERSION\s+VERSION_LESS\s+([0-9][0-9.]*)", text
    )
    return match.group(1) if match else None


def split_bound(bound):
    """Return (operator, version tuple) for a bound, or None if unparsable."""
    if bound is None:
        return None
    match = re.match(r"^([<>=]+)\s*([0-9][0-9.]*)$", bound.strip())
    if not match:
        return None
    return match.group(1), tuple(int(n) for n in match.group(2).split("."))


def satisfies_floor(bound, floor):
    """Is this bound at least as strict as a ">= floor" requirement?

    Compares versions numerically, so ">= 3.4.0" satisfies a floor of "3.4"
    and ">= 3.5" does too.  A missing, non-">=" or lower bound does not.
    """
    parsed = split_bound(bound)
    if parsed is None:
        return False
    operator, version = parsed
    if operator != ">=":
        return False
    wanted = tuple(int(n) for n in floor.split("."))
    length = max(len(version), len(wanted))
    version = version + (0,) * (length - len(version))
    wanted = wanted + (0,) * (length - len(wanted))
    return version >= wanted


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
        # An entry may offer alternatives ("a | b"); dpkg satisfies the build
        # from the first installable one, and only the first is recorded here.
        # The bound must be searched inside that alternative alone: searching
        # the whole entry would credit b's bound to a, so "libeigen3-dev |
        # other (>= 3.4)" would read as a floor this recipe does not have.
        alternative = entry.split("|")[0].strip()
        # The name ends at the first space, bracket or architecture qualifier.
        name = re.split(r"[\s(\[]", alternative, maxsplit=1)[0].strip()
        if not name:
            continue
        bound = re.search(r"\(\s*([<>=]+)\s*([^)]+?)\s*\)", alternative)
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
        # "eigen3-devel >= 3.4" is one package with a bound, but RPM also
        # allows "gcc-c++ make" as two packages on one line.  Tell them apart
        # by whether the tail starts with a comparison operator, so the second
        # package is recorded rather than mistaken for a version.
        if len(parts) > 1 and re.match(r"^[<>=]", parts[1].strip()):
            declared[parts[0]] = parts[1].strip()
        else:
            for name in entry.split():
                declared.setdefault(name, None)
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
    eigen_floor = parse_eigen_floor(DEPENDENCIES.read_text(encoding="utf-8"))
    lf_declared = packaging_is_lf_only(
        GITATTRIBUTES.read_text(encoding="utf-8")
    )
    lto_missing = lto_is_disabled(
        SPEC.read_text(encoding="utf-8"), RULES.read_text(encoding="utf-8")
    )
    service_versions = parse_service_versions(
        SERVICE.read_text(encoding="utf-8")
    )
    python_floors = parse_python_floors(
        PYPROJECT.read_text(encoding="utf-8"),
        CMAKELISTS.read_text(encoding="utf-8"),
        SPEC.read_text(encoding="utf-8"),
        CONTROL.read_text(encoding="utf-8"),
    )
    upstream = parse_upstream_version(CMAKELISTS.read_text(encoding="utf-8"))
    spec_upstream, spec_version = parse_spec_versions(
        SPEC.read_text(encoding="utf-8")
    )
    dsc_version, dsc_tarball = parse_dsc_version(DSC.read_text(encoding="utf-8"))
    changelog_version = parse_changelog_version(
        CHANGELOG.read_text(encoding="utf-8")
    )
    spec_sover, deb_sonames = parse_soname_major(
        SPEC.read_text(encoding="utf-8"),
        DSC.read_text(encoding="utf-8"),
        CONTROL.read_text(encoding="utf-8"),
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
        ("the Eigen floor in cmake/dependencies.cmake", eigen_floor),
        ("the version in CMakeLists.txt", upstream),
        ("%global upstream_version in coolprop.spec", spec_upstream),
        ("Version: in coolprop.spec", spec_version),
        ("Version: in coolprop.dsc", dsc_version),
        ("the tarball name in coolprop.dsc", dsc_tarball),
        ("the version in debian.changelog", changelog_version),
        ("%define sover in coolprop.spec", spec_sover),
        ("the libcoolprop<N> package name", deb_sonames or None),
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
            floor_applies = row[2]
            if name not in declared:
                problems.append(
                    "{0} builds with COOLPROP_VENDOR_THIRD_PARTY=OFF but does "
                    "not declare {1}".format(label, name)
                )
            elif floor_applies and not satisfies_floor(
                declared[name], eigen_floor
            ):
                problems.append(
                    "{0} declares {1} {2}, but cmake/dependencies.cmake "
                    "refuses to build below {3}".format(
                        label,
                        name,
                        describe_bound(declared[name]),
                        eigen_floor,
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

    # 5. Every recipe's version must follow CMakeLists.txt, which is what
    #    make-release-tarball.sh names the archive from.  Five files restate
    #    it, and the failures are late ones, inside the build service: rpmbuild
    #    cannot find the directory %autosetup names, and a stale soname makes
    #    %files claim a libCoolProp.so.N that was never built.
    #
    #    Note on the changelog: OBS's debtransform rewrites the top changelog
    #    entry to match the .dsc Version rather than failing, so that one is
    #    checked to keep the files honest for anyone reading or building them
    #    outside OBS, not because OBS would reject it.
    numeric, revision = upstream
    tarball_version = numeric + revision
    wanted = distro_version(numeric, revision)

    if spec_upstream != tarball_version:
        problems.append(
            "coolprop.spec builds from coolprop-{0}.tar.gz, but CMakeLists.txt "
            "produces coolprop-{1}.tar.gz".format(spec_upstream, tarball_version)
        )
    if spec_version != wanted:
        problems.append(
            "coolprop.spec has Version: {0}, but CMakeLists.txt means "
            "{1}".format(spec_version, wanted)
        )
    # A full match, not a prefix: "8.0.1~dev-" would otherwise pass, and dpkg
    # rejects an empty revision.  An epoch ("1:8.0.1~dev-1") is also refused;
    # it is legitimate dpkg syntax and the escape hatch for walking a version
    # back, so if one is ever needed this check is what to relax, deliberately.
    if not re.fullmatch(re.escape(wanted) + r"-[A-Za-z0-9.+~]+", dsc_version):
        problems.append(
            "coolprop.dsc has Version: {0}, but CMakeLists.txt means {1} with a "
            "non-empty Debian revision suffix".format(dsc_version, wanted)
        )
    if dsc_tarball != "coolprop-{0}.tar.gz".format(tarball_version):
        problems.append(
            "coolprop.dsc names {0}, but make-release-tarball.sh produces "
            "coolprop-{1}.tar.gz".format(dsc_tarball, tarball_version)
        )
    if changelog_version != dsc_version:
        problems.append(
            "debian.changelog is {0} but coolprop.dsc is {1}".format(
                changelog_version, dsc_version
            )
        )

    expected_python = (
        "pyproject.toml requires-python",
        "the find_package(Python ...) floor in CMakeLists.txt",
        "the weakest python BuildRequires in coolprop.spec",
        "the python3 Build-Depends in debian.control",
    )
    unread = [w for w in expected_python if w not in python_floors]
    if unread:
        problems.append(
            "could not read the Python floor from {0}; nothing then checks "
            "that the build root has an interpreter new enough to run "
            "dev/generate_headers.py".format(", ".join(unread))
        )
    else:
        too_low = sorted(
            "{0} allows {1}.{2}".format(where, *python_floors[where])
            for where in expected_python
            if python_floors[where] < REQUIRED_PYTHON
        )
        if too_low:
            problems.append(
                "dev/generate_headers.py needs Python {0}.{1}, but {2}".format(
                    REQUIRED_PYTHON[0], REQUIRED_PYTHON[1], "; ".join(too_low)
                )
            )

    missing_sites = [
        site for site in SERVICE_VERSION_SITES if site not in service_versions
    ]
    if missing_sites:
        problems.append(
            "_service: could not read {0}; the version there is then checked "
            "by nothing".format(", ".join(missing_sites))
        )
    else:
        wrong = sorted(
            "{0} says {1}".format(site, service_versions[site])
            for site in SERVICE_VERSION_SITES
            if service_versions[site] != numeric
        )
        if wrong:
            problems.append(
                "_service must name the release this series becomes, {0}, but "
                "{1}".format(numeric, "; ".join(wrong))
            )

    if lto_missing:
        problems.append(
            "link-time optimisation is not disabled in {0}; incbin's .incbin "
            "directive cannot be resolved when lto-wrapper re-assembles it at "
            "link time, and the build fails with 'file not found: "
            "all_fluids.cbor'".format(" and ".join(lto_missing))
        )

    if not lf_declared:
        problems.append(
            ".gitattributes does not declare 'dev/packaging/** text eol=lf', "
            "so a Windows checkout gives these recipes CRLF and rpmbuild fails "
            "in %prep on the stray carriage return"
        )

    major = numeric.split(".")[0]
    if spec_sover != major:
        problems.append(
            "coolprop.spec defines sover {0}, but COOLPROP_VERSION_MAJOR is "
            "{1}, which is the SOVERSION the library is built with".format(
                spec_sover, major
            )
        )
    wrong_sonames = [n for n in deb_sonames if n != major]
    if wrong_sonames:
        problems.append(
            "the Debian files name libcoolprop{0}, but COOLPROP_VERSION_MAJOR "
            "is {1}".format(", libcoolprop".join(wrong_sonames), major)
        )

    if problems:
        report(problems)
        return 1

    print("ok: the recipes agree on build dependencies and configure options")
    print("  debian.control / coolprop.dsc: " + format_declared(control_deps))
    print("  coolprop.spec:                 " + format_declared(spec_reqs))
    print("  version:                       {0} (tarball coolprop-{1}.tar.gz)".format(
        wanted, tarball_version))
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
