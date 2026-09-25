# Linux distribution packaging

This directory holds everything needed to build CoolProp packages for the main
Linux distributions on the [openSUSE Build Service][obs] (OBS).  The background,
the decision record and the staging plan live in
[GH #3388](https://github.com/CoolProp/CoolProp/issues/3388); this file is the
operational half.

[obs]: https://build.opensuse.org

```
dev/packaging/
  make-release-tarball.sh     build the offline source tarball
  vendor-deps.sh              fill externals/cpm/ so a build needs no network
  check-build-deps.py         assert the recipes below declare the same build
                              dependencies and configure options (run by CI)
  obs/
    _service                  bootstrap-only: fetch + checksum a release tarball
    coolprop.spec             RPM recipe   (openSUSE, SLE, Fedora, RHEL/EPEL)
    coolprop.dsc              Debian source control
    debian.control            Debian binary packages
    debian.rules              Debian build rules
    debian.changelog          Debian changelog
    debian.copyright          Debian copyright file (Policy 12.5 requires one)
    debian.libcoolprop8.install
    debian.libcoolprop-dev.install
```

## Why an offline tarball

Every distribution build system (OBS, `sbuild`, `mock`, `makepkg` in a clean
chroot) builds inside a sandbox with **no network at all**.  CoolProp fetches
ten dependencies with CPM at configure time, so a plain `git archive` of the
tree cannot be built by any of them.

`make-release-tarball.sh` produces a tarball that can: it exports the tree at a
given revision, runs `vendor-deps.sh` to copy every CPM dependency into
`externals/cpm/`, and writes `dev/gitrevision.txt` so that
`dev/generate_headers.py` does not need git either.

`cmake/dependencies.cmake` picks `externals/cpm/<name>/` up on its own.  Nothing
changes for a normal developer checkout, which has no `externals/cpm/`.

```bash
# from a clean checkout, with network access
dev/packaging/make-release-tarball.sh --ref v8.0.1
# -> dist/coolprop-8.0.1.tar.gz
#    dist/coolprop-8.0.1.tar.gz.sha256

# without --ref it follows the working tree, which between releases is a
# snapshot and is named accordingly:
dev/packaging/make-release-tarball.sh
# -> dist/coolprop-8.0.1dev.tar.gz
```

The tarball is reproducible for a given revision: entries sorted, ownership
zeroed, timestamps taken from the commit date.  Rebuilding it does not churn
the checksum that `_service` and the release notes pin.

### Verifying it really is offline

`vendor-deps.sh` finishes by reconfiguring with
`-DCOOLPROP_REQUIRE_VENDORED_DEPS=ON`, which makes
`cmake/dependencies.cmake` abort if **any** package resolved from outside
`externals/cpm/`.  That check walks `CPM_PACKAGES`, the list CPM itself keeps,
rather than a hand-written list of dependency names, so a dependency somebody
adds to `cmake/dependencies.cmake` later cannot slip past it.

Both packaging recipes pass `-DCOOLPROP_REQUIRE_VENDORED_DEPS=ON` as well, so
the failure mode in the OBS sandbox is a clear configure error rather than a
download timeout.

## The FHS install layout

The FHS layout itself comes from GH #3311, which landed in master while this
branch was open.  A packaging build asks for it with four options:

- `-DCOOLPROP_INSTALL_CMAKE_PACKAGE=ON` installs the system layout below
- `-DCOOLPROP_INSTALL_LEGACY_LAYOUT=OFF` drops CoolProp's release-artifact
  folders (`shared_library/Linux/64bit_GNU_11/`, which is what the SourceForge
  uploads want, and which has no business in `/usr`)
- `-DCOOLPROP_VENDOR_THIRD_PARTY=OFF` keeps the bundled Eigen and fmt out of
  `/usr/include`; leaving it ON is right for a relocatable SDK and is the one
  thing a distribution will reject
- `-DCOOLPROP_INSTALL_FLAT_HEADERS=OFF` keeps the deprecated flat shims
  (`AbstractState.h`, `Solvers.h`, `Exceptions.h`, `Ice.h` and friends, GH
  #1280) out of the include root.  They forward to the canonical
  `<CoolProp/*.h>` and their names are far too generic for `/usr/include`.
  Leaving it ON is right for an SDK, which is why that is the default; for an
  RPM it is not optional, because the recipes ship only
  `%{_includedir}/CoolProp/` and rpmbuild aborts on "Installed (but
  unpackaged) files found"

`COOLPROP_VENDOR_THIRD_PARTY=OFF` is not free: it makes
`cmake/dependencies.cmake` resolve Eigen and fmt with
`find_package(... CONFIG REQUIRED)` instead of CPM, so both have to be
installed at **build** time, not just depended on at install time.  A recipe
that turns the option off without declaring them fails in the configure step
with "Could not find a package configuration file provided by Eigen3".  That is
why `coolprop.spec` carries `BuildRequires: eigen3-devel` and
`BuildRequires: fmt-devel`, `debian.control` carries `libeigen3-dev` and
`libfmt-dev` in `Build-Depends`, and `packaging_offline.yml` installs the same
two packages before it takes the network away.

That gives:

| What | Where |
|---|---|
| shared library | `${CMAKE_INSTALL_LIBDIR}`, so `/usr/lib/x86_64-linux-gnu` on Debian and `/usr/lib64` on Fedora |
| headers | `${CMAKE_INSTALL_INCLUDEDIR}/CoolProp/` |
| pkg-config | `${CMAKE_INSTALL_LIBDIR}/pkgconfig/coolprop.pc` |
| CMake package | `${CMAKE_INSTALL_LIBDIR}/cmake/CoolProp/` |

Of these, only the pkg-config file comes from this branch; GH #3311 supplies
the rest.  This branch originally carried its own `cmake/CoolPropInstall.cmake`
doing the same job, and that was deleted when #3311 merged rather than left to
fight with it: two `install(EXPORT)` rules over the same targets is not a merge
conflict CMake reports, it is one that ships.

So a downstream project can finally do either of:

```cmake
find_package(CoolProp REQUIRED)
target_link_libraries(my_app PRIVATE CoolProp::CoolProp)
```

```bash
c++ my_app.cpp $(pkg-config --cflags --libs coolprop)
```

Two deliberate differences from the release archives:

- **The flat legacy headers are not installed.**  `AbstractState.h`,
  `Solvers.h`, `MatrixMath.h`, `Exceptions.h`, `Ice.h` and about two dozen more
  are deprecation shims (GH #1280, to be removed at v9) that forward to
  `CoolProp/<name>.h`.  Putting names that generic directly into
  `/usr/include` would collide with other packages.  Distribution consumers
  include `<CoolProp/CoolProp.h>`.
- **`CMAKE_INSTALL_PREFIX` is honoured.**  CoolProp used to force it to
  `<source>/install_root` with `CACHE ... FORCE`, which silently swallowed
  `-DCMAKE_INSTALL_PREFIX=/usr`.  It now only does that when no prefix was
  asked for, so `%cmake` and `dh_auto_configure` work while a plain
  `cmake -B build -S .` still collects artefacts in `install_root/` the way the
  wrapper build jobs expect.  The check runs before `project()`, where
  `CMAKE_INSTALL_PREFIX` is still empty unless the user put it in the cache from
  the command line; `packaging_offline.yml` asserts all three outcomes.

### Known limitation: Eigen and fmt version skew

The installed C++ headers are not self-contained: they include `<Eigen/Dense>`,
and `<fmt/format.h>` unless the consumer defines `NO_FMTLIB`.  What follows
from that depends on which way `COOLPROP_VENDOR_THIRD_PARTY` is set, and the
two cases are genuinely different:

- **ON** (the default, and what an SDK-style install wants): CoolProp builds
  against its CPM-pinned Eigen 5.0.1 and fmt 12.0.0 and installs them under
  `<prefix>/include/CoolProp/third_party/`, so a consumer compiles against
  those same copies.  Library and consumer agree by construction.
- **OFF** (what every distribution build passes): CoolProp builds against the
  distribution's Eigen and fmt and installs neither, so a consumer compiles
  against those same system copies.  Library and consumer agree here too.

So a packaged build has no library-versus-consumer skew.  The hazard it does
have is a different one, and it is worth stating precisely rather than as
"version skew": CoolProp is developed and tested against its pinned Eigen
5.0.1, and a distribution build compiles it against whatever that distribution
ships instead.

That is not one number, and it is not safe to assume it.  A rolling
distribution may already carry the same 5.x CoolProp pins, in which case there
is nothing to worry about; a stable release may sit near the 3.4 floor, which
is a different major.  One data point from a real OBS build root, 2026-09-24:
openSUSE Tumbleweed installed `eigen3-devel 5.0.1` and `fmt-devel 12.1.0`.
Eigen there is exactly the version CoolProp pins; fmt is one minor ahead of the
pinned 12.0.0, so even on a rolling distribution the skew is not zero, it is
just small.  So the exposure has to be read per target, from the
Eigen version that target actually packages, rather than assumed from this
file.  `cmake/dependencies.cmake` enforces the 3.4 floor and the recipes
declare it, but "configures and compiles" is not "behaves identically", and
the test suite has only ever been run against the pinned version.

`coolprop.pc` used to ship an empty `Requires:`, which meant
`pkg-config --cflags coolprop` emitted no include path for Eigen, so a **C++**
consumer could not compile at all: Debian keeps Eigen in `/usr/include/eigen3`
rather than `/usr/include`.  That is fixed.  The OFF build now writes
`Requires.private: eigen3 fmt`, and the ON build instead names the bundled
copies under `<includedir>/CoolProp/third_party/` in `Cflags`.  Private rather
than public because both are header-only here, so pkg-config always emits
their include flags and keeps their link flags behind `--static`.

The CI consumer test now compiles `dev/ci/cmake-consumer/cpp_api.cpp` through
pkg-config to keep it that way.  It was C-only before, deliberately, which is
exactly why a `.pc` no C++ consumer could use went unnoticed.

Declaring the dependency does not settle the version question above.  It only
stops the file being unusable.

For the **C API** (`<CoolProp/CoolPropLib.h>`) this does not arise; it is a
plain C interface over a compiled library.

For the **C++ API** it is a real hazard, and it is why step 2 of GH #3388,
`COOLPROP_USE_SYSTEM_DEPS=ON` with `find_package()` for the six dependencies
that distributions already package, has to land before the `-dev` package can
be called production-ready.  Until then:

- `coolprop.pc` names `eigen3` and `fmt` without a version bound, so
  pkg-config accepts whatever the distribution ships.  The floor is enforced by
  CMake and by the recipes' own build dependencies instead.
- The `-dev` packages still depend on `libeigen3-dev` / `eigen3-devel` and the
  fmt equivalents, so the headers a consumer needs are at least present.

### Known limitation: the packages are built without LTO

CoolProp embeds its fluid database (`dev/all_fluids.cbor`) into the library
with incbin.  That works by emitting a top-level `__asm__` holding a `.incbin`
directive, and the assembler locates the file through the `-I` paths the
compiler passes it, one of which is `dev/`.

Link-time optimisation breaks that.  Under `-flto`, top-level assembly is
streamed into the LTO objects and re-assembled at link time by `lto-wrapper`,
which runs from `/tmp` with its own option set.  The original `-I` is gone, so
the link fails:

```
/tmp/ccXXXXXX.s:46: Error: file not found: all_fluids.cbor
lto-wrapper: fatal error: make returned 2 exit status
ld: error: lto-wrapper failed
```

openSUSE and Fedora both put `-flto=auto` in `%optflags`, so this hits every
RPM target; it failed real Tumbleweed builds on x86_64 and i586 alike, at the
final link after a full compile.  Debian does not enable LTO by default.

So `coolprop.spec` sets `%define _lto_cflags %{nil}` and `debian.rules` adds
`optimize=-lto`, and `check-build-deps.py` fails if either guard is dropped.
The cost is the cross-object inlining LTO would have bought.

The fix that would let LTO back in is to hand incbin an absolute path, so the
re-assembly can still find the file no matter which directory it runs from.
That has been confirmed to work, but it changes the main build system rather
than a packaging recipe, so it belongs in its own change where CoolProp's test
suite can be run against it with LTO enabled.

### Known limitation: `debian.copyright` is not per-dependency yet

`debian.copyright` is **incomplete, and says so in its own header**.  It is
good enough for our own OBS repositories and is not good enough for a
distribution archive.  Two gaps:

**The in-git third-party code is enumerated only as far as a grep reached.**
The tarball is a `git archive` of the whole tree, and that tree carries code
under licences other than MIT.  Stanzas exist for the ones found - ExternalMedia
under `wrappers/Modelica/src/` (**Modelica License 2**, not MIT), the BSD-licensed
CMake find-modules under `dev/cmake/Modules/` taken from Ceres and GDCM,
`wrappers/Lua/lualib.mk` (ISC), `wrappers/MATLAB/`, `externals/incbin`,
`externals/miniz-3.1.1`, `cmake/CPM.cmake` and `wrappers/Rust` - each verified
against the file it names.  What is *not* established is that the set is
complete, so `Files: *` may still cover third-party code.  Closing this needs a
`licensecheck(1)` pass over an unpacked tarball, which is a review task rather
than a scripted one.

**The vendored sources are not enumerated per dependency.**  For the trees
under `externals/cpm/`, the file points at the licence inside each rather than
restating terms that could drift.  Those licences are **not uniform and not all
MIT** - Eigen is MPL-2.0, msgpack-c and Catch2 are BSL-1.0, valijson is
BSD-2-Clause, and the terms for IF97, REFPROP-headers, boost-headers and
multicomplex have not been established at all.  A blanket "MIT or compatible"
claim over that set would simply be false.

Both are step 6 concerns in the GH #3388 staging (archive inclusion), not step
4 ones (our own OBS repos), which is why they are recorded rather than done.
Landing `COOLPROP_USE_SYSTEM_DEPS=ON` shrinks the second one first, since a
dependency taken from the distribution is not vendored and needs no stanza.
Narrowing what the tarball ships (a `.gitattributes` `export-ignore` over the
wrapper trees a distribution build does not compile) would shrink the first,
and is worth considering on its own merits.

## Setting up OBS

### 1. Account and project

1. Sign up at <https://build.opensuse.org> (free for open source).  An OBS
   account is also an openSUSE account; it is not tied to SUSE employment.
2. Everyone gets a home project, `home:<username>`.  Use that for the first
   run-through.
3. Once the packages build, request a proper namespace such as
   `science:CoolProp` by opening a request against the `science` project, or
   keep `home:<username>:CoolProp` and point users at it.  A devel project
   under `science:` is the more discoverable home for a numerics library.

Install the client locally:

```bash
# openSUSE
sudo zypper install osc
# Debian/Ubuntu
sudo apt install osc
# Fedora
sudo dnf install osc
```

`osc` reads credentials from `~/.config/osc/oscrc`; `osc` will prompt on the
first command.  Prefer an API token (Profile -> Tokens on the web UI) over
storing the account password.

### 2. Create the package

```bash
osc checkout home:<username>
cd home:<username>
osc mkpac CoolProp
cd CoolProp

cp /path/to/CoolProp/dev/packaging/obs/* .
cp /path/to/CoolProp/dist/coolprop-8.0.1dev.tar.gz .

osc add *
osc commit -m "Initial CoolProp packaging"
```

### Versions

`CMakeLists.txt` is the single source of truth.  `make-release-tarball.sh` names
the archive from it, so between releases it emits `coolprop-8.0.1dev.tar.gz`,
and the recipes follow.

The distribution version is **not** the same string as the tarball version, and
that is deliberate.  A pre-release has to sort *below* the release it precedes,
and a bare suffix does the opposite:

```
$ dpkg --compare-versions 8.0.1~dev lt 8.0.1   # true
$ dpkg --compare-versions 8.0.1dev  lt 8.0.1   # FALSE
```

so a snapshot called `8.0.1dev` would outrank the release and block the upgrade
to it.  Both RPM and dpkg spell a pre-release with `~`.  Hence:

| | between releases | on a release tag |
|---|---|---|
| tarball | `coolprop-8.0.1dev.tar.gz` | `coolprop-8.0.1.tar.gz` |
| `coolprop.spec` | `%global upstream_version 8.0.1dev`, `Version: 8.0.1~dev` | both `8.0.1` |
| `coolprop.dsc` | `Version: 8.0.1~dev-1` | `8.0.1-1` |
| `debian.changelog` | `(8.0.1~dev-1)` | `(8.0.1-1)` |

`dev/packaging/check-build-deps.py` checks all four against `CMakeLists.txt` on
every CI run, so they cannot drift apart silently.  To cut a release, set
`COOLPROP_VERSION_REVISION` to empty and update the four; the checker tells you
if you miss one.

**`_service` is still a placeholder**: 64 zeros for the checksum and a release
asset URL that does not exist yet, so `osc service manualrun` fails closed
rather than fetching something unchecked.  It is the tagged-release path and is
not used for a snapshot trial, where the tarball is added with `osc add` --
by hand for a one-off, or by the CI job in section 6, which does the same
thing on every push.

Do not do both in the same package directory.  `dev/packaging/obs-upload.sh`
handles this by removing every `*.tar.gz` from the checkout before copying the
new one in, and by never uploading `_service`.  `coolprop.dsc` carries no
`Debtransform-Tar:` line, so OBS's `debtransform` finds the source archive by
scanning the directory, and it aborts with "Too many files looking like a
usable source tarball" when a hand-added `coolprop-8.0.1dev.tar.gz` sits beside
a `coolprop-8.0.1.tar.gz` that `_service` fetched.

**If you check out on Windows, mind the line endings.**  `.gitattributes` sets
`* text=auto`, so a text file is converted to the platform's native ending on
checkout.  For these recipes that is fatal rather than cosmetic: rpmbuild
writes the `%prep` body into a shell script, and a CR there is executed as a
command, which ends the build with

```
/var/tmp/rpm-tmp.XXXXXX: line 46: $'\r': command not found
error: Bad exit status from /var/tmp/rpm-tmp.XXXXXX (%prep)
```

`dev/packaging/**` and `dev/ci/**` are therefore pinned to `eol=lf`, since they
are only ever consumed by Linux build systems.

If you already have a checkout with CRLF in these files, note that
`git add --renormalize .` alone does not fix it: that updates the index, not
the files on disk, so the copy you would upload to OBS still has the CRLF and
still fails in `%prep`.  Delete the files and check them out again, which
applies the new `eol=lf` attribute:

```bash
# Delete the tracked files, then check them out again.  git ls-files is used
# rather than "rm -rf dev/packaging" so that untracked files in those
# directories, such as a tarball you already built, are left alone.
git ls-files -z dev/packaging dev/ci | xargs -0 rm -f
git checkout -- dev/packaging dev/ci
```

Re-cloning works too.  Verify with `file dev/packaging/obs/coolprop.spec`,
which must not say "CRLF line terminators".

### rpmlint fails the build after a successful compile

openSUSE runs rpmlint after packaging and aborts when the accumulated "badness"
exceeds 1000. A single error can be worth 10000, so a build that compiled,
linked, installed and produced all three RPMs still ends in
`failed "build coolprop.spec"`. The log looks like a toolchain failure and is
not one; read past the compiler warnings to the `RPMLINT report:` block near
the end.

One rule has already caught us, and will again if anyone tidies the spec:

- **A package holding one shared library must be named after its SONAME.** Ours
  is `libCoolProp.so.8`, so the package is `libCoolProp8`, capitals included.
  Naming it `libcoolprop8` gives
  `E: shlib-policy-name-error (Badness: 10000) libCoolProp8`, where the name in
  brackets is what rpmlint wanted. The `-devel` package is not covered by that
  rule and stays lower case, which is also what the install instructions above
  tell users to type.

  **The RPM and the Debian package therefore have different names on purpose**,
  and that is not something to tidy up:

  | | Package | Rule |
  |---|---|---|
  | `debian.control` | `libcoolprop8` | Debian Policy allows lower case only |
  | `coolprop.spec` | `libCoolProp8` | openSUSE derives the name from the SONAME |

  The two policies genuinely contradict each other, so each file follows its
  own and neither is wrong. Making them agree would mean renaming the library
  itself to `libcoolprop.so.8`, which is not a packaging detail: `coolprop.pc`
  advertises `-lCoolProp`, the CMake package exports `CoolProp::CoolProp`, and
  the same library is `CoolProp.dll` on Windows and `libCoolProp.dylib` on
  macOS. That would break every existing consumer to make two package names
  match, so the mixed case stays.

Two warnings are known and accepted, because warnings score nothing against the
threshold:

- `no-binary` on the `coolprop` package, which ships only the licence and the
  README. Making it `noarch` is not possible from the main preamble without
  making the shared library noarch too.
- `macro-in-comment`, if a `%` macro is written unescaped in a spec comment.
  Macros expand inside comments as well, so double the percent sign: `%%files`.

`osc build` runs rpmlint locally, so section 4 below catches all of this before
a push.

### 3. Choose build targets

In the web UI, **Repositories -> Add from a distribution**.

**Take the exact repository names from that dialog rather than from this
file.** Distribution releases go end of life on their own schedule, so any list
written down here is wrong within a year or two, and a packager who copies a
stale name selects a target that no longer receives updates. What follows is
therefore the shape of a reasonable set, not a list to copy:

| Pick | Why |
|---|---|
| `openSUSE_Tumbleweed` | rolling openSUSE, and the one to develop against |
| the current openSUSE Leap | the regular-release openSUSE users actually run |
| the two newest Fedora releases | Fedora supports roughly the last two |
| the current CentOS Stream | stands in for RHEL and its rebuilds; needs EPEL, see below |
| Debian stable, and oldstable while it is still supported | what most Debian users have |
| the Ubuntu LTS releases still in standard support | Ubuntu users mostly track LTS |
| `Arch` | possible, though the AUR below is the better home for it |

Enable `x86_64` everywhere and `aarch64` where the distribution offers it.

`eigen3-devel` and `fmt-devel` are not in the RHEL or CentOS Stream base
repositories, they come from EPEL, so a CentOS Stream target is unresolvable
until EPEL is added to it; see "Adding EPEL to a CentOS Stream target" below.
The openSUSE and Fedora targets carry both packages themselves.

Debian and Ubuntu targets build from `coolprop.dsc` plus the `debian.*` files;
the RPM targets build from `coolprop.spec`.  OBS decides per repository, so the
two recipes live side by side in one package directory.

#### Editing the target list as XML

The web UI dialog writes into the project metadata, and that metadata can be
edited directly instead.  This is the faster route once more than one or two
targets are involved, and it is the only route that shows the two settings the
dialog hides (a second `<path>`, and a per-target architecture switch).

```bash
osc meta prj -e home:<username>            # opens $EDITOR on the live metadata
osc meta prj home:<username> > prj.xml     # or: dump, edit, upload
osc meta prj home:<username> -F prj.xml
```

Each target is one `<repository>` block.  The `name` is what appears in the
download URL and in `osc build`, the `<path>` says which distribution supplies
the build root, and each `<arch>` is one architecture OBS will schedule:

```xml
<repository name="openSUSE_Tumbleweed">
  <path project="openSUSE:Factory" repository="snapshot"/>
  <arch>x86_64</arch>
  <arch>aarch64</arch>
</repository>
```

The project and repository names inside `<path>` are OBS names, not
distribution names, and they change as releases come and go.  Add one target
through **Repositories -> Add from a distribution** first, then dump the
metadata and copy the shape: that way the names come from the server rather
than from anyone's memory.

Removing a target is deleting its block.  OBS does not retire a repository when
the distribution behind it goes end of life, it keeps scheduling builds against
a base that no longer receives updates, so end-of-life targets have to be taken
out by hand.

#### Adding EPEL to a CentOS Stream target

`eigen3-devel` and `fmt-devel` are not in the CentOS Stream base repositories,
which is why a CentOS Stream target reports `unresolvable` rather than failing
to compile: OBS cannot assemble a build root that satisfies the `BuildRequires`
in the spec, so nothing is ever built.

The fix is a second `<path>` in that repository block, pointing at whichever
project on the server mirrors EPEL for that release.  Paths are searched in the
order they are listed, so the base distribution goes first and EPEL after it:

```xml
<repository name="CentOS_Stream">
  <path project="CENTOS_STREAM_PROJECT" repository="standard"/>
  <path project="EPEL_PROJECT_FOR_THAT_RELEASE" repository="standard"/>
  <arch>x86_64</arch>
</repository>
```

Take both names from the **Add from a distribution** dialog, which lists what
the server actually carries.  Until the second path is there, the target stays
unresolvable no matter what the spec says.

#### Turning off one architecture without dropping the target

Deleting an `<arch>` line stops that architecture being scheduled, but it also
discards the build history behind it.  To keep the target and silence one
architecture, for instance while an `aarch64` failure is being investigated,
disable it instead:

```xml
<build>
  <disable repository="openSUSE_Tumbleweed" arch="aarch64"/>
</build>
```

The `<build>` element sits at the end of the project metadata, after the
`<repository>` blocks.  The same element is accepted in the package metadata
(`osc meta pkg -e home:<username> coolprop`), which is the better place for a
switch that is about this package rather than about the project.

A disabled architecture reports `disabled` in `osc results`.  That is not
counted as a failure by the CI job in section 6, deliberately: a target that is
switched off on purpose should not turn the GitHub check red.  Do note that
switching off every target would leave nothing to check, which the job treats
as a failure in its own right.

#### Project-wide build settings

Macros, package preferences and similar build-root settings live in a separate
document, the project config:

```bash
osc meta prjconf -e home:<username>
```

Nothing in this packaging currently needs an entry there.  It is worth knowing
about because build failures that look like a broken spec, for example two
packages both providing the same dependency, are usually resolved with a
`Prefer:` line in the project config rather than by changing the spec.

### 4. Build locally before pushing

`osc build` runs the same build in a local chroot, which is far faster than
waiting on the server and gives you the full log:

```bash
osc build openSUSE_Tumbleweed x86_64 coolprop.spec
osc build Debian_13 x86_64 coolprop.dsc
```

Substitute whichever target names you actually enabled above; the Debian one
here is only an example.

The first run downloads a base chroot and needs root (`osc build` uses `sudo`
for that, or set `su-wrapper` in `oscrc`).

This stays a local step, and the CI job in section 6 does not use `osc build`.
The reasons are worth recording, because it looks at first like the obvious way
to test packaging on a GitHub runner:

- it is not offline.  `osc build` fetches the build configuration and the
  package list from the OBS API before it starts, so it needs the same
  credentials and the same reachable server as pushing does.  It gives no
  independence from OBS, only a different place to run the compiler.
- it needs a privileged chroot on the runner.
- Ubuntu's `obs-build`, which supplies the distribution definitions, is a
  snapshot several years old, so the newest targets are not described in it.

A container running `rpmbuild` plus `rpmlint` directly, with no OBS involved,
is the cheaper way to reproduce an RPM target on a runner if that is ever
wanted.  It needs no credentials, so it also works on pull requests from forks.

### 5. What users then do

OBS publishes a signed repository per target.  For openSUSE:

```bash
sudo zypper addrepo https://download.opensuse.org/repositories/home:/<username>/openSUSE_Tumbleweed/home:<username>.repo
sudo zypper install libcoolprop-devel
```

For Debian/Ubuntu the project page's "Go to download repository" link gives the
exact `deb` line and the signing key.

### 6. Building on OBS from CI (the testing loop)

The `obs-build` job in `.github/workflows/packaging_offline.yml` is how the
packaging is exercised day to day. It runs after `offline-build`, downloads the
tarball that job produced, pushes it **and every recipe file in
`dev/packaging/obs/`** to the OBS package with `osc`, waits for the builds, and
fails the run if any target failed. The OBS outcome therefore shows up as an
ordinary check on the commit, next to the other CI jobs.

That replaces the hand upload. The tarball and the spec go up together, so a
spec change is always tested against the code it belongs with. Uploading a new
spec on top of an old tarball is the mistake that cost a full round of
debugging on the `-m32`/`-m64` failures, and this makes it impossible.

It shares a workflow with the offline build on purpose. The obvious gain is
that the tarball is built once rather than twice, but the one that matters more
is that OBS builds the **byte-identical** artifact the offline job just proved,
rather than a second tarball that ought to be the same. `_service` states that
as the goal; downloading the artifact is what delivers it. A useful side effect
is ordering: a tarball that already failed to build offline never reaches OBS,
so an obvious break costs seconds rather than a slow OBS round trip.

The two jobs are not redundant. `offline-build` unpacks the tarball, installs
the dependencies the distribution recipes declare, builds with the network
switched off, checks the installed layout, and compiles a consumer against the
staged install through both `pkg-config` and `find_package`. OBS builds real
RPMs and DEBs on real distributions but does not test a downstream consumer,
and it needs credentials and an external service. Keeping the fast, credential
free job as the gate and OBS as the distribution reality check is the point of
the split.

Setup, once. The job is **inert until `OBS_PROJECT` is set**, so nothing
happens in a fork or in a clone with no OBS package to push to. A manual run
can supply the project instead of the variable, but only somebody with write
access can start one.

| Kind | Name | Value |
|---|---|---|
| Variable | `OBS_PROJECT` | e.g. `home:jowr`. Setting this is what switches the job on for automatic runs. |
| Variable | `OBS_PACKAGE` | Optional, defaults to `coolprop`. |
| Variable | `OBS_APIURL` | Optional, defaults to `https://api.opensuse.org`. |
| Secret | `OBS_USERNAME` | Your OBS account name. |
| Secret | `OBS_PASSWORD` | Your OBS password. |

Set these under **Settings -> Secrets and variables -> Actions**, on the
Variables and Secrets tabs respectively. `osc` reads the credentials from the
environment, so neither value appears in a command line.

No `oscrc` is written on the runner, but only because the script sets
`OSC_CONFIG=/dev/null`. That is not a tidiness measure. osc resolves its config
file before it consults the environment, so on a machine with no
`~/.config/osc/oscrc` it announces that it is going to create one and blocks
waiting for a username on stdin; the environment credentials are never reached.
Pointing `OSC_CONFIG` at `/dev/null` makes osc start from an empty config and
use them.

If `OBS_PROJECT` is set but the credentials are missing, the run fails with a
clear message rather than skipping. That is deliberate: a packaging gate that
quietly does nothing is worse than one that is switched off, because it still
shows a green tick.

Two behaviours worth knowing before the first run:

- **An unresolvable target fails the run.** `osc results --fail-on-error`
  treats `failed`, `broken` and `unresolvable` as failures, so a
  `BuildRequires:` that does not exist on a target (the reason CentOS 8 Stream
  and possibly Leap's `gcc13-c++` fall over) turns the check red rather than
  passing quietly. `excluded` and `disabled` are not failures, so a repository
  that does not build a given architecture is fine.
- **A rebuild that OBS never reacts to fails the run.** After the commit there
  is a window where OBS still reports the previous build as finished, and
  reading results in that window would pass on stale data. So when the sources
  changed, the script snapshots the build state before committing and then
  waits until OBS either marks the package dirty or reports something different
  from that snapshot, and treats running out of time as a failure. Comparing
  against a snapshot rather than looking for failure codes is deliberate: this
  package already has failing targets, and any check that simply looked for a
  failed state would match the leftovers from the previous build on its first
  poll and report those instead. When the sources did NOT change, there is
  nothing for OBS to react to, and the run reports the existing results without
  waiting.
- **A result list with no rows fails the run.** `--fail-on-error` decides from
  the rows it iterated, so zero rows exits 0. A package with no build targets
  enabled, or a mistyped `OBS_PACKAGE`, would otherwise look like a clean
  build, so the row count is asserted separately.

The job installs `osc` from PyPI rather than from apt, and that is not a
preference. Ubuntu 24.04 ships osc 0.169.1, while reading credentials from
environment variables arrived in osc 1.6.0 and `results --fail-on-error` in
1.8.0. With the apt version the job would prompt for a password and would never
report a failing target, so the pin must stay at or above those versions.

The `obs-build` jobs are serialised on a concurrency group keyed by the
project and package they write to, because a single OBS package cannot hold two
commits at once. Runs aimed at different packages, which is what the
`obs_project` input makes possible, do not hold each other up. They queue
rather than cancelling each other; see the second limit below for why that
matters.

When it runs, and when it does not:

- **On packaging changes, not on code changes.** The workflow's `paths:` filter
  covers `CMakeLists.txt`, `cmake/**`, `dev/packaging/**`,
  `dev/generate_headers.py` and the workflow file itself. A change under `src/`
  or `include/` does not trigger an OBS build. That keeps the loop fast while
  the packaging is what is being worked on; widen the filter if OBS should also
  be a check on the library itself.
- **On pushes to master, main and develop, and on pull requests to them.** The
  pull request case is the day to day loop: every push to a branch with an open
  pull request gets an OBS result. A branch with no pull request open does not,
  and needs `workflow_dispatch`.
- **Not on release tags.** The workflow itself does fire on `v*`, but
  `obs-build` skips them. On a release tag the OBS package belongs to the
  release path in section 7, whose services fetch a published tarball and check
  its sha256 against the same project and package. Letting CI commit its own
  tarball at that moment would leave whoever revives that path debugging
  sources that were quietly replaced.
- **Not for a pull request from a fork.** A fork gets no secrets, so the job is
  skipped rather than failing a contributor's pull request on a credential they
  were never going to have. Note a skipped job reads as neutral, so do not make
  `obs-build` a required check without thinking that through.

#### Triggering a build by hand

**Actions -> Packaging (offline build) -> Run workflow**, then pick the branch.
Three inputs:

| Input | Default | What it does |
|---|---|---|
| `push_to_obs` | ticked | Untick to run only the offline build and leave OBS alone. |
| `obs_project` | blank | A project for this run only. Blank uses the `OBS_PROJECT` variable. |
| `obs_package` | blank | Likewise for the package name. Blank uses `OBS_PACKAGE`, or `coolprop`. |

The whole workflow runs either way, because `obs-build` consumes the tarball
artifact that `offline-build` produces, so there is no path to OBS that does not
build the tarball first. A manual run therefore costs one offline build even
when all you wanted was the OBS half.

Two inputs rather than a variable is what makes a one-off build to a staging
project possible without editing the variable that every automatic run uses.
The credentials are deliberately not overridable: they are secrets, and a
dispatch input would record the value in the run's parameters in plain text.

A manual run ignores the `paths:` filter, so it is also the way to get an OBS
build for a change that the filter does not cover, or for a branch with no pull
request open.

Two behaviours specific to manual runs:

- **Asking for an OBS build with no project configured fails the run**, rather
  than skipping the job the way an automatic run does. Skipping is right when it
  means "this repository has no OBS package", but somebody who ticked
  `push_to_obs` asked for an OBS build in so many words, and answering that with
  a skipped job and a green run is the quiet no-op this setup exists to avoid.
  The first step in the job says what is missing, before it installs anything.
- **Dispatching against a tag fails the run**, rather than skipping quietly,
  for the same reason. On a release tag the package belongs to the release path
  in section 7, so the job refuses and says so. An automatic run on a tag is
  still skipped, because nobody asked for anything there.
- **The project and package names are checked for shape, not just emptiness.**
  Any non-empty input wins over the variable, so a single space left behind by
  a copy-paste would otherwise pass an empty-looking name and silently override
  a perfectly good `OBS_PROJECT`. A name with a space in it is rejected.
- **The credentials are checked in the same step.** A repository that was never
  set up finds out before `osc` is installed rather than several steps later.

Two limits to be aware of before leaning on this:

- **Everything that triggers it automatically pushes to the same OBS
  package**, including every push to an open pull request, which is the highest
  churn case of all. `OBS_PROJECT` names one project, and `home:jowr` is a
  published repository that users install from, so any of these replaces what
  is published there until the next run. While that is only ever your own home
  project this is a nuisance rather than a hazard, but point `OBS_PROJECT` at a
  throwaway staging project if iteration and a repository other people consume
  ever become the same place. A manual run escapes this with the `obs_project`
  input; an automatic one cannot.
- **Overlapping runs settle in completion order, not push order.** The
  `obs-build` jobs queue on one concurrency group, and a job joins that group
  only once `offline-build` has finished, so a run whose offline build was slow
  can reach OBS after a newer one. Each run still commits its own tarball and
  reports the result for the sources it committed, so no check ever lies; but
  the package can be left holding the older sources until the next run. Queuing
  rather than cancelling is deliberate, so that every commit gets a verdict
  instead of the overtaken one being cancelled.

### 7. Rebuilding automatically from GitHub (the release path)

`.obs/workflows.yml` in the repository root tells OBS what to do when GitHub
calls it: on a release tag it re-runs the package's source services, which
fetch the tarball from the GitHub release and verify its sha256.

**The file alone does nothing.** OBS only reads it if a token and a webhook
exist, and both are created by hand, once:

1. On OBS, create a workflow token and give it a GitHub personal access token
   so OBS can read the repository and report back. Do this in the browser,
   under **Profile -> Manage Your Tokens**, choosing the `workflow` operation
   and pasting the GitHub token into the form.

   It returns a token id and a secret. Keep the secret; it is shown once.

   `osc` can do the same thing:

   ```bash
   osc token --create --operation workflow --scm-token <github-pat>
   ```

   Prefer the browser anyway. The GitHub token is a command-line argument
   here, so it lands in your shell history and is readable in
   `/proc/<pid>/cmdline` by anyone else on the machine for as long as the
   command runs. That matters little on a personal laptop and rather more on
   a shared build host.

2. On GitHub, under **Settings -> Webhooks -> Add webhook** for the repository:

   - Payload URL: `https://build.opensuse.org/trigger/workflow?id=<token-id>`
   - Content type: `application/json`
   - Secret: the secret from step 1
   - Events: "Let me select individual events", then **Pushes**. OBS derives
     its `tag_push` event from a push whose ref is under `refs/tags/`, so the
     separate "Branch or tag creation" event is not needed.

3. `.obs/workflows.yml` has to be on the repository's default branch. OBS reads
   it from there, not from the branch that triggered the event.

Setting the token and the webhook up does not yet make the trigger do
anything. Three things are in the way, and all three have to be dealt with
before a tag rebuilds the package:

1. **The services are skipped.** Both services in `dev/packaging/obs/_service`
   carry `mode="manual"`. A server-side services run, which is what
   `trigger_services` asks for, skips exactly the modes `localonly`,
   `disabled`, `manual` and `buildtime`; the list is in the OBS backend, in
   `bs_service` under "collect services to run". So the webhook logs
   `Skip download_url` and `Skip verify_file` and does nothing else. Dropping
   the mode makes them run, at the price of re-downloading the tarball on
   every source change, which is the thing the mode was added to avoid.
2. **There is no tarball at the URL.** `_service` fetches a GitHub release
   asset under `releases/download/<tag>/`, and nothing publishes one.
   `release_all_files.yml` does build the tarball with
   `make-release-tarball.sh`, but it uploads it as a workflow artifact and
   rsyncs it to SourceForge, and its `tags: ['v*']` trigger is commented out,
   so no tag runs it at all.
3. **`_service` names the wrong version.** It still pins v8.0.1 with an
   all-zero placeholder checksum.

None of this can publish a wrong tarball. With the services skipped nothing
happens at all, and if they were un-skipped today the download would 404 or
`verify_file` would reject the checksum. Wiring the release job to attach the
tarball and rewrite `_service` is item 3 under "Keeping it alive" below.

There is also an ordering point that survives all three fixes: `tag_push`
fires when the tag is pushed, which is before the release workflow has
uploaded anything. The dependable trigger is therefore that workflow calling
OBS after its own upload, with this webhook as the manual re-run path.

Things to check on the first tag rather than assume:

- **How often it fires.** `.obs/workflows.yml` has no tag filter, so every
  tag triggers a services run, including the `gui-v*` ones. Once the three
  items above are fixed, that means an unrelated tag re-fetches whatever
  version `_service` names at the time. The file explains why a filter was
  removed rather than left in unverified, but this is the cost of that
  choice and is worth watching.
- **Which project the workflow targets.** It currently names `home:jowr`, a
  personal project. That is the right place while the packaging is being
  proven, and the wrong place afterwards: a project owned by the CoolProp
  organisation should own the published packages, and the workflow file should
  be updated with it.

There is deliberately no pull-request workflow here. OBS can branch a package
and build it per pull request, but this package's source is a release tarball
rather than the git checkout, so such a build would compile whatever tarball
the package already holds and report green without having tested the pull
request at all. CoolProp's own CI covers pull requests; see
`.github/workflows/packaging_offline.yml`, which builds the tarball with the
network switched off.

## Keeping it alive

The single clearest lesson from the 2014 attempt recorded in GH #3388 is that
**a packaging path nobody exercises is already dead**.  `wrappers/DEB/` broke
the day a `cmake-format` pass removed a space, and nobody found out for twelve
years, because no CI job ran it.

So:

1. `.github/workflows/packaging_offline.yml` builds CoolProp from a
   freshly-made offline tarball, with the network blocked, on every push that
   touches the build system or this directory.  If dependency handling regresses,
   that job goes red on the commit that caused it, not months later.
2. The `prepare_sources` job in `.github/workflows/release_all_files.yml`
   builds the tarball as a second output alongside the existing
   `CoolProp_sources.zip`, so both land in `binaries/source/` and get rsynced
   to SourceForge with everything else.  That zip is not a substitute: it has
   no vendored dependencies (so it cannot build in a chroot), unpacks as
   `source/` rather than `coolprop-<version>/`, carries no version in its name
   and is not reproducible.
3. Still to wire up: pushing the tarball to OBS with `osc` from that same job,
   and attaching it to the GitHub release.  Version numbers appear in four
   places that must agree, so have CI rewrite them rather than a human:
   `coolprop.spec` (`Version:`), `coolprop.dsc` (`Version:`, `Files:`),
   `debian.changelog` (top entry) and `_service`.  OBS's `set_version` service
   can do the spec and changelog; the `.dsc` needs the tarball name too.
   `dpkg-source` rejects a changelog version that does not match the tarball,
   which is one of the things that would have stopped the 2014 scripts even if
   their version parser had worked.
4. Never let a packaging check fail open.  `vendor-deps.sh` aborts if CPM
   reports no packages at all rather than declaring an empty tree vendored, and
   `make-release-tarball.sh` asserts that each version component it parsed is a
   number instead of building `coolprop_..orig.tar.gz` the way the 2014 script
   did.

## Cheaper wins worth doing alongside

- **AUR**: a `PKGBUILD` is about thirty lines and the Arch community maintains
  it.  Covers Arch for nearly nothing.
- **conda-forge**: a feedstock plus the autotick bot.  For CoolProp's actual
  audience, engineers and scientists on managed clusters and mixed operating
  systems, this probably reaches more users than Debian and Fedora combined,
  and maintenance is approving bot pull requests.

Both are listed as step 5 in GH #3388.
