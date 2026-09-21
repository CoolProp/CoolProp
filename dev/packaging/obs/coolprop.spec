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
Version:        8.0.1
Release:        0
Summary:        Thermophysical property library for pure fluids, mixtures and humid air
License:        MIT
URL:            https://www.coolprop.org
Source0:        coolprop-%{version}.tar.gz

BuildRequires:  cmake >= 3.14
BuildRequires:  gcc-c++
BuildRequires:  pkgconfig
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

%package devel
Summary:        Development files for CoolProp
Requires:       libcoolprop%{sover} = %{version}-%{release}
Requires:       eigen3-devel
%if 0%{?suse_version}
Group:          Development/Libraries/C and C++
Requires:       fmt-devel
%else
Requires:       fmt-devel
%endif

%description devel
Headers, the pkg-config file and the CMake package configuration needed to
build against CoolProp.

Consumers include <CoolProp/CoolProp.h> for the C++ API or
<CoolProp/CoolPropLib.h> for the C API.  The flat legacy headers that the
upstream release archives carry (AbstractState.h and friends, directly in the
include root) are not shipped here: their names are too generic to put into
%{_includedir}.

%prep
%autosetup -n coolprop-%{version}

%build
# COOLPROP_REQUIRE_VENDORED_DEPS makes the configure step fail if anything
# would have to be downloaded, rather than letting it fail later and less
# clearly inside the OBS sandbox.
#
# COOLPROP_PC_REQUIRES and COOLPROP_EXPORTED_DEPENDENCIES are deliberately left
# unset.  CoolProp is compiled against its own pinned Eigen and fmt, so telling
# a consumer to find the distribution's copies would claim a compatibility
# nobody has established yet.  Set them to "eigen3 fmt" and "Eigen3;fmt" once
# COOLPROP_USE_SYSTEM_DEPS exists (step 2 of GH #3388); see
# dev/packaging/README.md.
%cmake \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DCOOLPROP_SHARED_LIBRARY=ON \
    -DCOOLPROP_SYSTEM_INSTALL=ON \
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

%files devel
%{_includedir}/CoolProp/
%{_libdir}/libCoolProp.so
%{_libdir}/pkgconfig/coolprop.pc
%{_libdir}/cmake/CoolProp/

%changelog
