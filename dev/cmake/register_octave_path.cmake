# Add "addpath('<PREFIX>/<SUBDIR>');" to ~/.octaverc (once), so the installed
# CoolProp.oct is on Octave's path at every startup.
#
# Run as a standalone `cmake -P` process from the Octave install step in the
# top-level CMakeLists.txt, with -DPREFIX=<install prefix> and
# -DSUBDIR=<install dir below the prefix>. Any failure here (no home
# directory, unreadable or unwritable ~/.octaverc, and so on) ends only this
# process with an error, which the caller turns into a warning, so a failed
# registration never aborts `cmake --install`.
#
# HOME and DESTDIR are read from the environment of the installing user, not
# of whoever ran the original `cmake ..`.
if(NOT DEFINED PREFIX OR NOT DEFINED SUBDIR)
  message(FATAL_ERROR "register_octave_path.cmake requires -DPREFIX=... -DSUBDIR=...")
endif()

# Single quotes are doubled so a prefix containing one stays a valid Octave string.
string(REPLACE "'" "''" prefix "${PREFIX}")
set(line "addpath('${prefix}/${SUBDIR}');")
message(STATUS "Octave path line for CoolProp: ${line}")

if(DEFINED ENV{DESTDIR} AND NOT "$ENV{DESTDIR}" STREQUAL "")
  message(STATUS "DESTDIR staging detected ($ENV{DESTDIR}); skipping ~/.octaverc registration")
  return()
endif()

if(DEFINED ENV{HOME})
  file(TO_CMAKE_PATH "$ENV{HOME}" home)
elseif(DEFINED ENV{USERPROFILE})
  file(TO_CMAKE_PATH "$ENV{USERPROFILE}" home)
else()
  message(FATAL_ERROR "Neither HOME nor USERPROFILE is set")
endif()
if(NOT IS_DIRECTORY "${home}")
  message(FATAL_ERROR "${home} (from HOME/USERPROFILE) is not a directory")
endif()

set(rc "${home}/.octaverc")
set(existing "")
if(EXISTS "${rc}")
  file(READ "${rc}" existing)
endif()
string(FIND "${existing}" "${line}" pos)
if(NOT pos EQUAL -1)
  message(STATUS "${rc} already registers CoolProp")
  return()
endif()

set(sep "")
if(NOT existing STREQUAL "" AND NOT existing MATCHES "\n$")
  set(sep "\n")
endif()
file(APPEND "${rc}" "${sep}${line}\n")
message(STATUS "Added CoolProp to ${rc}")
