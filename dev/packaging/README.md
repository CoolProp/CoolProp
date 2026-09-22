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

`-DCOOLPROP_SYSTEM_INSTALL=ON` switches the install rules from CoolProp's
release-artifact folders (`shared_library/Linux/64bit_GNU_11/`, which is what
the SourceForge uploads want) to a normal system layout:

| What | Where |
|---|---|
| shared library | `${CMAKE_INSTALL_LIBDIR}`, so `/usr/lib/x86_64-linux-gnu` on Debian and `/usr/lib64` on Fedora |
| headers | `${CMAKE_INSTALL_INCLUDEDIR}/CoolProp/` |
| pkg-config | `${CMAKE_INSTALL_LIBDIR}/pkgconfig/coolprop.pc` |
| CMake package | `${CMAKE_INSTALL_LIBDIR}/cmake/CoolProp/` |

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
and `<fmt/format.h>` unless the consumer defines `NO_FMTLIB`.  CoolProp builds
against its own CPM-pinned Eigen 5.0.1 and fmt 12.0.0 and does not install
them, so a consumer compiles those same headers against whatever the
distribution ships, which is a different Eigen major on most of them.

For the **C API** (`<CoolProp/CoolPropLib.h>`) this does not arise; it is a
plain C interface over a compiled library.

For the **C++ API** it is a real hazard, and it is why step 2 of GH #3388,
`COOLPROP_USE_SYSTEM_DEPS=ON` with `find_package()` for the six dependencies
that distributions already package, has to land before the `-dev` package can
be called production-ready.  Until then:

- `coolprop.pc` ships with an empty `Requires:` and `CoolPropConfig.cmake`
  with no `find_dependency()` calls, because naming `eigen3` there would claim
  a compatibility that has not been established.
- Both are one `cmake -D` away when it has been:
  `-DCOOLPROP_PC_REQUIRES="eigen3 fmt"` and
  `-DCOOLPROP_EXPORTED_DEPENDENCIES="Eigen3;fmt"`.
- The `-dev` packages still depend on `libeigen3-dev` / `eigen3-devel` and the
  fmt equivalents, so the headers a consumer needs are at least present.

### Known limitation: `debian.copyright` is not per-dependency yet

`debian.copyright` names the copyright holders and licences for the upstream
tree and for the third-party code carried in git (`externals/incbin`,
`externals/miniz-3.1.1`, `cmake/CPM.cmake`, `wrappers/Rust`).  For the sources
the release tarball vendors under `externals/cpm/`, it points at the licence
file inside each tree instead of restating the terms.

That is honest but not sufficient for a distribution archive.  Those licences
are **not uniform and not all MIT** - Eigen is MPL-2.0, msgpack-c and Catch2
are BSL-1.0, valijson is BSD-2-Clause - so an archive submission needs a
`Files:` stanza per vendored dependency with its real short name and text.
Writing one blanket "MIT or compatible" claim over the set would be wrong, and
a copyright file is the last place to guess.

This is a step 6 concern in the GH #3388 staging (archive inclusion), not a
step 4 one (our own OBS repos), which is why it is recorded rather than done.
Landing `COOLPROP_USE_SYSTEM_DEPS=ON` shrinks the problem first, since a
dependency taken from the distribution is not vendored and needs no stanza.

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
cp /path/to/CoolProp/dist/coolprop-8.0.1.tar.gz .

osc add *
osc commit -m "Initial CoolProp packaging"
```

**The recipes carry placeholder values, on purpose, and nothing in CI checks
them yet.** Fix these by hand for the first run, then automate it as described
under "Keeping it alive" below:

| File | Placeholder | Why it is there |
|---|---|---|
| `coolprop.spec` | `Version: 8.0.1` | the tree is `8.0.1dev`, so `make-release-tarball.sh` currently emits `coolprop-8.0.1dev.tar.gz` and `%autosetup -n coolprop-%{version}` will not find it |
| `coolprop.dsc` | `Version: 8.0.1-1`, and a `Files:` line of zeros | OBS's `debtransform` normally rewrites the `Files:` stanza, so the zeros are expected to work, but that is unverified here |
| `debian.changelog` | `(8.0.1-1)` | `dpkg-source` rejects a changelog version that does not match the tarball |
| `_service` | 64 zeros for the checksum, and a release asset URL that does not exist yet | so `osc service manualrun` fails closed rather than fetching something unchecked |

Tag a release, or set `COOLPROP_VERSION_REVISION` to empty in `CMakeLists.txt`,
and all four line up at `8.0.1`.

### 3. Choose build targets

In the web UI, **Repositories -> Add from a distribution**.  A reasonable
starting set:

| Repository | Gives you |
|---|---|
| `openSUSE_Tumbleweed` | rolling openSUSE |
| `openSUSE_Leap_15.6` | current Leap |
| `Fedora_41`, `Fedora_42` | Fedora |
| `CentOS_9_Stream` | RHEL-compatible |
| `Debian_12`, `Debian_13` | Debian stable and next |
| `xUbuntu_24.04`, `xUbuntu_22.04` | Ubuntu LTS |
| `Arch` | Arch (community-maintained is usually better; see AUR below) |

Enable `x86_64` everywhere and `aarch64` where the distribution offers it.
Each repository can be added from the command line too, by editing the project
metadata with `osc meta prj -e home:<username>`.

Debian and Ubuntu targets build from `coolprop.dsc` plus the `debian.*` files;
the RPM targets build from `coolprop.spec`.  OBS decides per repository, so the
two recipes live side by side in one package directory.

### 4. Build locally before pushing

`osc build` runs the same build in a local chroot, which is far faster than
waiting on the server and gives you the full log:

```bash
osc build openSUSE_Tumbleweed x86_64 coolprop.spec
osc build Debian_12 x86_64 coolprop.dsc
```

The first run downloads a base chroot and needs root (`osc build` uses `sudo`
for that, or set `su-wrapper` in `oscrc`).

### 5. What users then do

OBS publishes a signed repository per target.  For openSUSE:

```bash
sudo zypper addrepo https://download.opensuse.org/repositories/home:/<username>/openSUSE_Tumbleweed/home:<username>.repo
sudo zypper install libcoolprop-devel
```

For Debian/Ubuntu the project page's "Go to download repository" link gives the
exact `deb` line and the signing key.

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
