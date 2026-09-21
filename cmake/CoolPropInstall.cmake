# FHS install layout for CoolProp, enabled with COOLPROP_SYSTEM_INSTALL=ON.
#
# The historical CoolProp install rules drop the library into release-artifact
# folders such as shared_library/Linux/64bit_GNU_11/, which is exactly what the
# SourceForge uploads want.  A Linux distribution package needs the opposite: a
# library in /usr/lib<multiarch-or-64>, headers in /usr/include, a pkg-config
# file and a CMake package config, so that a downstream project can write
#
#     find_package(CoolProp REQUIRED)          # -> CoolProp::CoolProp
#     pkg-config --cflags --libs coolprop
#
# and get a working build.  Neither of those worked before; see GH #3388.
#
# This file is included from CMakeLists.txt with the library target
# (${LIB_NAME}) already created, so it only adds install rules.

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# --------------------------------------------------------------------------
# Public headers
# --------------------------------------------------------------------------
#
# Ship the include/CoolProp/ tree (GH #1280), minus detail/json.h and
# detail/msgpack.h, which are not self-contained: they pull in nlohmann,
# valijson and msgpack.hpp, none of which this package installs.
# dev/ci/check-installed-headers.sh compile-checks what ends up here.
#
# The flat back-compat shims (include/*.h forwarding to include/CoolProp/*.h,
# to be removed at v9) are deliberately NOT installed in this layout.  Their
# names are generic -- Solvers.h, MatrixMath.h, Exceptions.h, Configuration.h,
# Ice.h -- and dropping about thirty of those straight into /usr/include would
# collide with other packages.  A distribution consumer includes
# <CoolProp/CoolProp.h>; the release-artifact layout still carries the shims
# for everyone building against an unpacked SourceForge archive.
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/CoolProp
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" FILES_MATCHING PATTERN "*.h"
        REGEX "detail/(json|msgpack)\\.h$" EXCLUDE)

# --------------------------------------------------------------------------
# The library itself
# --------------------------------------------------------------------------
#
# ${LIB_NAME} was given the whole private include list (Eigen, fmt, msgpack,
# ...) as PUBLIC include directories, which is right for an in-tree build but
# cannot be exported: those paths are build-machine scratch directories that do
# not exist for whoever installs the package.  Rewrite the interface so the
# build tree keeps them ($<BUILD_INTERFACE:...>) and the installed package
# points at the installed header tree instead ($<INSTALL_INTERFACE:...>).
get_target_property(_coolprop_iface_includes ${LIB_NAME}
                    INTERFACE_INCLUDE_DIRECTORIES)
set(_coolprop_build_includes "")
if(_coolprop_iface_includes)
  foreach(_dir IN LISTS _coolprop_iface_includes)
    list(APPEND _coolprop_build_includes "$<BUILD_INTERFACE:${_dir}>")
  endforeach()
endif()
set_property(TARGET ${LIB_NAME} PROPERTY INTERFACE_INCLUDE_DIRECTORIES
             ${_coolprop_build_includes}
             "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")

# An alias so that a project pulling CoolProp in with add_subdirectory() links
# against the same name it would get from find_package(CoolProp).
if(NOT TARGET CoolProp::CoolProp)
  add_library(CoolProp::CoolProp ALIAS ${LIB_NAME})
endif()

# No INCLUDES DESTINATION here: the $<INSTALL_INTERFACE:...> set above already
# puts ${CMAKE_INSTALL_INCLUDEDIR} into the exported interface, and naming it
# twice just lists the same path twice in CoolPropTargets.cmake.
install(
  TARGETS ${LIB_NAME}
  EXPORT CoolPropTargets
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")

# --------------------------------------------------------------------------
# CMake package config: lib/cmake/CoolProp/CoolPropConfig.cmake
# --------------------------------------------------------------------------

set(COOLPROP_CMAKE_CONFIG_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/CoolProp")

install(
  EXPORT CoolPropTargets
  FILE CoolPropTargets.cmake
  NAMESPACE CoolProp::
  DESTINATION "${COOLPROP_CMAKE_CONFIG_DESTINATION}")

# The installed headers include <Eigen/Dense>, and <fmt/format.h> unless the
# consumer defines NO_FMTLIB, so a consumer has to find those packages too.
# CoolProp itself builds against its own CPM-fetched copies and does not
# install them, which means only whoever is packaging knows what the consumer
# should find.  Hence a variable rather than a hard-coded list, empty by
# default so that a plain developer install does not demand packages that are
# not on the machine.
set(COOLPROP_EXPORTED_DEPENDENCIES
    ""
    CACHE STRING
    "CMake packages that CoolPropConfig.cmake should find_dependency(), e.g. \"Eigen3;fmt\"")

set(COOLPROP_CONFIG_DEPENDENCIES "")
foreach(_dep IN LISTS COOLPROP_EXPORTED_DEPENDENCIES)
  string(APPEND COOLPROP_CONFIG_DEPENDENCIES "find_dependency(${_dep})\n")
endforeach()

configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/CoolPropConfig.cmake.in"
  "${PROJECT_BINARY_DIR}/CoolPropConfig.cmake"
  INSTALL_DESTINATION "${COOLPROP_CMAKE_CONFIG_DESTINATION}")

# SameMajorVersion: CoolProp's SOVERSION is the major version, so any 8.x
# satisfies a request for 8.y as far as ABI is concerned.
write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/CoolPropConfigVersion.cmake"
  VERSION "${COOLPROP_VERSION_MAJOR}.${COOLPROP_VERSION_MINOR}.${COOLPROP_VERSION_PATCH}"
  COMPATIBILITY SameMajorVersion)

install(FILES "${PROJECT_BINARY_DIR}/CoolPropConfig.cmake"
              "${PROJECT_BINARY_DIR}/CoolPropConfigVersion.cmake"
        DESTINATION "${COOLPROP_CMAKE_CONFIG_DESTINATION}")

# --------------------------------------------------------------------------
# pkg-config: lib/pkgconfig/coolprop.pc
# --------------------------------------------------------------------------
#
# Express libdir and includedir relative to ${prefix} where GNUInstallDirs gave
# us a relative directory; an absolute CMAKE_INSTALL_LIBDIR/INCLUDEDIR is passed
# through unchanged.
#
# Note that prefix= itself is the absolute CMAKE_INSTALL_PREFIX baked in at
# install time, so a relocated tree needs the consumer's help.  `pkg-config
# --define-prefix` is NOT that help here: it recomputes the prefix by stripping
# two components from the .pc file's own directory, which is wrong as soon as
# libdir is multiarch (lib/x86_64-linux-gnu/pkgconfig is three components) and
# silently yields -I<root>/lib/include.  `--define-variable=prefix=<root>`
# works.  Distribution packages install to the prefix they were configured
# with, so this only bites someone relocating a tarball by hand.
if(IS_ABSOLUTE "${CMAKE_INSTALL_LIBDIR}")
  set(COOLPROP_PC_LIBDIR "${CMAKE_INSTALL_LIBDIR}")
else()
  set(COOLPROP_PC_LIBDIR "\${exec_prefix}/${CMAKE_INSTALL_LIBDIR}")
endif()
if(IS_ABSOLUTE "${CMAKE_INSTALL_INCLUDEDIR}")
  set(COOLPROP_PC_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}")
else()
  set(COOLPROP_PC_INCLUDEDIR "\${prefix}/${CMAKE_INSTALL_INCLUDEDIR}")
endif()

# Anything including <CoolProp/AbstractState.h> needs Eigen, and fmt too unless
# it defines NO_FMTLIB, because those headers are not self-contained.  CoolProp
# does not install its own (CPM-fetched) copies, so a package that ships this
# layout has to name the distribution's own packages here.  Left empty by
# default: on a plain developer install there may well be no eigen3.pc on the
# machine, and an unsatisfiable Requires makes pkg-config fail outright.
set(COOLPROP_PC_REQUIRES
    ""
    CACHE STRING
    "pkg-config Requires: line, e.g. \"eigen3 fmt\" when packaging for a distribution")

# ${CMAKE_DL_LIBS} is "dl" on glibc and empty elsewhere.  It only matters when
# someone links the static library, hence Libs.private.
if(CMAKE_DL_LIBS)
  set(COOLPROP_PC_LIBS_PRIVATE "-l${CMAKE_DL_LIBS}")
else()
  set(COOLPROP_PC_LIBS_PRIVATE "")
endif()

configure_file("${PROJECT_SOURCE_DIR}/cmake/coolprop.pc.in"
               "${PROJECT_BINARY_DIR}/coolprop.pc" @ONLY)

install(FILES "${PROJECT_BINARY_DIR}/coolprop.pc"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")

message(STATUS "COOLPROP_SYSTEM_INSTALL: installing to an FHS layout")
message(STATUS "  libraries -> ${CMAKE_INSTALL_FULL_LIBDIR}")
message(STATUS "  headers   -> ${CMAKE_INSTALL_FULL_INCLUDEDIR}")
