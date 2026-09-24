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
BuildRequires:  gcc-c++
BuildRequires:  pkgconfig
# COOLPROP_VENDOR_THIRD_PARTY=OFF below resolves Eigen and fmt with
# find_package, so they must be present at build time as well as at
# install time.  Without these two the configure step fails with
# "Could not find a package configuration file provided by Eigen3".
BuildRequires:  eigen3-devel >= 3.4
BuildRequires:  fmt-devel
# dev/generate_headers.py runs at build time to turn the fluid JSON into the
# generated headers that get compiled into the library.
%if 0%{?suse_version}
Group:          Development/Libraries/C and C++
BuildRequires:  python3-base
%else
BuildRequires:  python3
%endif

%description
CoolProp is a C++ library that implements pure and pseudo-pure fluid equations
of state and transport properties for 122 components, mixture properties using
high-accuracy Helmholtz energy formulations, correlations of incompressible
fluids and brines, and psychrometric routines for humid air.

This package contains the documentation and licence only; see
libcoolprop%{sover} for the shared library and libcoolprop-devel to build
against it.

%package -n libcoolprop%{sover}
Summary:        Thermophysical property library for pure fluids, mixtures and humid air
%if 0%{?suse_version}
Group:          System/Libraries
%endif

%description -n libcoolprop%{sover}
CoolProp is a C++ library that implements pure and pseudo-pure fluid equations
of state and transport properties for 122 components, mixture properties using
high-accuracy Helmholtz energy formulations, correlations of incompressible
fluids and brines, and psychrometric routines for humid air.

This package contains the shared library.  The fluid data is compiled into the
binary, so there is no separate data package and no runtime data path to
configure.

%package -n libcoolprop-devel
Summary:        Development files for CoolProp
Requires:       libcoolprop%{sover} = %{version}-%{release}
Requires:       eigen3-devel
Requires:       fmt-devel
%if 0%{?suse_version}
Group:          Development/Libraries/C and C++
%else
# Owns %{_libdir}/cmake on Fedora and RHEL; see the %files section below.
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
    -DCOOLPROP_NO_EXAMPLES=ON

%cmake_build

%install
%cmake_install

%post   -n libcoolprop%{sover} -p /sbin/ldconfig
%postun -n libcoolprop%{sover} -p /sbin/ldconfig

%files
%license LICENSE
%doc README.md

%files -n libcoolprop%{sover}
%{_libdir}/libCoolProp.so.%{sover}
# libCoolProp.so.8.0.1 on a release, libCoolProp.so.8.0.1dev on a snapshot:
# CMakeLists.txt appends COOLPROP_VERSION_REVISION to the library VERSION.
%{_libdir}/libCoolProp.so.%{sover}.*

%files -n libcoolprop-devel
%{_includedir}/CoolProp/
%{_libdir}/libCoolProp.so
%{_libdir}/pkgconfig/coolprop.pc
# Fedora ships cmake-filesystem, which owns %{_libdir}/cmake, so depend on it
# there (claiming the directory as well would be dual ownership and a review
# flag).  openSUSE has no such package, so claim it there instead; an unowned
# directory is an rpmlint error and OBS will bounce it.
%if 0%{?suse_version}
%dir %{_libdir}/cmake
%endif
%{_libdir}/cmake/CoolProp/

%changelog
* Mon Sep 21 2026 CoolProp developers <coolprop@coolprop.org> - 8.0.1-0
- Initial packaging for the openSUSE Build Service (GH #3388).
