#
# spec file for package coolprop
#
# Built on the openSUSE Build Service for openSUSE, SLE, Fedora and RHEL/EPEL.
# See dev/packaging/README.md in the CoolProp source tree for the whole setup.
#
# The tarball this consumes is produced by dev/packaging/make-release-tarball.sh
# and carries every third-party dependency under externals/cpm/, because OBS
# builds inside a sandbox with no network access.
#

# CoolProp embeds its fluid database (dev/all_fluids.cbor) into the library with
# incbin, which emits a top-level __asm__ carrying a .incbin directive.  The
# assembler resolves that filename through the -I paths the compiler hands it.
#
# Link-time optimisation breaks that.  With -flto, top-level asm is streamed
# into the LTO objects and re-assembled at link time by lto-wrapper, which runs
# from /tmp with its own option set, so the -I .../dev of the compile step is
# gone and the assembler cannot find the file:
#
#     /tmp/ccXXXXXX.s:46: Error: file not found: all_fluids.cbor
#     lto-wrapper: fatal error: make returned 2 exit status
#
# openSUSE and Fedora both put -flto=auto in %%optflags, so this hits every RPM
# target (reproduced on Tumbleweed x86_64 and i586 alike).  LTO is therefore
# switched off for this package.  The fix that would keep LTO is to give incbin
# an absolute path so the re-assembly can still find the file; that belongs
# upstream in the build system, not in a packaging recipe.
%define _lto_cflags %{nil}

# Shared library major version.  It follows COOLPROP_VERSION_MAJOR, which
# CMakeLists.txt sets as the target SOVERSION, so that two CoolProp majors can
# be installed side by side during a transition.
%define sover 8

Name:           coolprop
# The tarball name and the RPM version are deliberately not the same string.
#
# make-release-tarball.sh names the archive from CMakeLists.txt, which carries
# COOLPROP_VERSION_REVISION=dev between releases, so a snapshot unpacks as
# coolprop-8.0.1dev/.  That suffix cannot be the RPM Version: "8.0.1dev" sorts
# ABOVE "8.0.1", so a snapshot would shadow the release it precedes and block
# the upgrade.  RPM spells a pre-release with "~", which sorts below.
#
# Hence two values: upstream_version follows the tarball, Version follows RPM's
# ordering rules.  On a release tag COOLPROP_VERSION_REVISION is empty and both
# become plain 8.0.1.  dev/packaging/check-build-deps.py checks they agree with
# CMakeLists.txt and with the Debian files.
%global upstream_version 8.0.1dev
Version:        8.0.1~dev
Release:        0%{?dist}
Summary:        Thermophysical property library for pure fluids, mixtures and humid air
License:        MIT
URL:            https://www.coolprop.org
Source0:        coolprop-%{upstream_version}.tar.gz

BuildRequires:  cmake >= 3.14
# CoolProp uses std::filesystem, so it needs GCC 9 or newer; CMakeLists.txt
# states that floor and refuses anything older.  Leap 15.x still installs GCC
# 7.5 as its default compiler, where the build died on
# "fatal error: filesystem: No such file or directory", so that distribution
# gets a newer compiler by name, exactly as it gets a newer Python.
%if 0%{?suse_version} && 0%{?suse_version} < 1600
BuildRequires:  gcc13-c++
%global cp_cc   %{_bindir}/gcc-13
%global cp_cxx  %{_bindir}/g++-13
%else
BuildRequires:  gcc-c++
%global cp_cc   %{_bindir}/gcc
%global cp_cxx  %{_bindir}/g++
%endif
BuildRequires:  pkgconfig
# COOLPROP_VENDOR_THIRD_PARTY=OFF below resolves Eigen and fmt with
# find_package, so they must be present at build time as well as at
# install time.  Without these two the configure step fails with
# "Could not find a package configuration file provided by Eigen3".
BuildRequires:  eigen3-devel >= 3.4
BuildRequires:  fmt-devel
# dev/generate_headers.py runs at build time to turn the fluid JSON into the
# generated headers that get compiled into the library.  It needs Python 3.9 or
# newer, which is what pyproject.toml declares and what CMakeLists.txt asks
# find_package for: the script uses builtin generics such as list[Path].
#
# Leap 15.x still ships Python 3.6 as its default python3, so there the build
# died with "TypeError: 'type' object is not subscriptable".  Ask that
# distribution for a separate, newer interpreter; CMake picks it up because the
# 3.9 floor makes it skip the 3.6 one and keep looking.
# cp_python is the interpreter the build is told to use, by absolute path.
# Leaving CMake to discover one is not good enough on Leap 15.x: python311
# installs /usr/bin/python3.11 and leaves /usr/bin/python3 at 3.6, and
# FindPython only probes the interpreter names its own version list knows, so
# an older CMake can fail to consider 3.11 at all and reject the 3.6 it does
# find.  Naming the path removes the guesswork.  CMake still checks that the
# interpreter satisfies the 3.9 floor, so a wrong path fails loudly.
%if 0%{?suse_version}
Group:          Development/Libraries/C and C++
%if 0%{?suse_version} < 1600
BuildRequires:  python311
%global cp_python %{_bindir}/python3.11
%else
BuildRequires:  python3-base >= 3.9
%global cp_python %{_bindir}/python3
%endif
%else
BuildRequires:  python3 >= 3.9
%global cp_python %{_bindir}/python3
%endif

%description
CoolProp is a C++ library that implements pure and pseudo-pure fluid equations
of state and transport properties for 122 components, mixture properties using
high-accuracy Helmholtz energy formulations, correlations of incompressible
fluids and brines, and psychrometric routines for humid air.

This package contains the documentation and licence only; see
libCoolProp%{sover} for the shared library and libcoolprop-devel to build
against it.

# The capitals are not a style choice.  openSUSE requires a package holding a
# single shared library to be named after that library's SONAME, which here is
# libCoolProp.so.8, giving libCoolProp8.  rpmlint scores a mismatch at badness
# 10000 against a threshold of 1000, so it fails the build AFTER everything has
# compiled and linked, which reads like a toolchain problem and is not one.
# Do not "normalise" this to lower case.  The -devel package below is free of
# that rule and stays lower case, matching what the README tells users to
# install.
%package -n libCoolProp%{sover}
Summary:        Thermophysical property library for pure fluids, mixtures and humid air
%if 0%{?suse_version}
Group:          System/Libraries
%endif

%description -n libCoolProp%{sover}
CoolProp is a C++ library that implements pure and pseudo-pure fluid equations
of state and transport properties for 122 components, mixture properties using
high-accuracy Helmholtz energy formulations, correlations of incompressible
fluids and brines, and psychrometric routines for humid air.

This package contains the shared library.  The fluid data is compiled into the
binary, so there is no separate data package and no runtime data path to
configure.

%package -n libcoolprop-devel
Summary:        Development files for CoolProp
Requires:       libCoolProp%{sover} = %{version}-%{release}
Requires:       eigen3-devel
Requires:       fmt-devel
%if 0%{?suse_version}
Group:          Development/Libraries/C and C++
%else
# Owns %%{_libdir}/cmake on Fedora and RHEL; see the %%files section below.
Requires:       cmake-filesystem
%endif

%description -n libcoolprop-devel
Headers, the pkg-config file and the CMake package configuration needed to
build against CoolProp.

Consumers include <CoolProp/CoolProp.h> for the C++ API or
<CoolProp/CoolPropLib.h> for the C API.  The flat legacy headers that the
upstream release archives carry (AbstractState.h and friends, directly in the
include root) are not shipped here: their names are too generic to put into
%{_includedir}.

%prep
%autosetup -n coolprop-%{upstream_version}

%build
# COOLPROP_REQUIRE_VENDORED_DEPS makes the configure step fail if anything
# would have to be downloaded, rather than letting it fail later and less
# clearly inside the OBS sandbox.
#
# COOLPROP_VENDOR_THIRD_PARTY=OFF keeps CoolProp's pinned Eigen and fmt out of
# /usr/include: a private copy of Eigen in a system include directory is the one
# thing a distribution will not accept.  The consumer therefore compiles the
# installed C++ headers against the distribution's own Eigen.
#
# COOLPROP_PC_REQUIRES is not set here because it does not need to be: with
# COOLPROP_VENDOR_THIRD_PARTY=OFF the CMake defaults it to "eigen3 fmt", which
# is what the build actually compiled against.  Set it only to override that.
%cmake \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DCOOLPROP_SHARED_LIBRARY=ON \
    -DCOOLPROP_INSTALL_CMAKE_PACKAGE=ON \
    -DCOOLPROP_INSTALL_LEGACY_LAYOUT=OFF \
    -DCOOLPROP_INSTALL_FLAT_HEADERS=OFF \
    -DCOOLPROP_VENDOR_THIRD_PARTY=OFF \
    -DCOOLPROP_REQUIRE_VENDORED_DEPS=ON \
    -DCOOLPROP_NO_EXAMPLES=ON \
    -DPython_EXECUTABLE=%{cp_python} \
    -DCMAKE_C_COMPILER=%{cp_cc} \
    -DCMAKE_CXX_COMPILER=%{cp_cxx}

%cmake_build

%install
%cmake_install

%post   -n libCoolProp%{sover} -p /sbin/ldconfig
%postun -n libCoolProp%{sover} -p /sbin/ldconfig

%files
%license LICENSE
%doc README.md

%files -n libCoolProp%{sover}
%{_libdir}/libCoolProp.so.%{sover}
# libCoolProp.so.8.0.1 on a release, libCoolProp.so.8.0.1dev on a snapshot:
# CMakeLists.txt appends COOLPROP_VERSION_REVISION to the library VERSION.
%{_libdir}/libCoolProp.so.%{sover}.*

%files -n libcoolprop-devel
%{_includedir}/CoolProp/
%{_libdir}/libCoolProp.so
%{_libdir}/pkgconfig/coolprop.pc
# Fedora ships cmake-filesystem, which owns %%{_libdir}/cmake, so depend on it
# there (claiming the directory as well would be dual ownership and a review
# flag).  openSUSE has no such package, so claim it there instead; an unowned
# directory is an rpmlint error and OBS will bounce it.
%if 0%{?suse_version}
%dir %{_libdir}/cmake
%endif
%{_libdir}/cmake/CoolProp/

%changelog
* Mon Sep 21 2026 CoolProp developers <coolprop@coolprop.org> - 8.0.1~dev-0
- Initial packaging for the openSUSE Build Service (GH #3388).
