#######################################
#       COMPILER CONFIGURATION        #
#-------------------------------------#
# Compiler flags, C++ standard,       #
# platform-specific settings          #
#######################################

# Force C++17 standard (lambdas are used in CPStrings.h, std::make_unique in DataStructures.cpp)
set(CMAKE_CXX_STANDARD 17)

# Address Sanitizer configuration for Clang
if(COOLPROP_ASAN)
  # https://stackoverflow.com/a/64294837 (CC BY-SA 4.0)
  if(isMultiConfig)
    if(NOT "Asan" IN_LIST CMAKE_CONFIGURATION_TYPES)
      list(APPEND CMAKE_CONFIGURATION_TYPES Asan)
    endif()
  else()
    set(allowedBuildTypes Asan Debug Release RelWithDebInfo MinSizeRel)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "${allowedBuildTypes}")

    if(CMAKE_BUILD_TYPE AND NOT CMAKE_BUILD_TYPE IN_LIST allowedBuildTypes)
      message(FATAL_ERROR "Invalid build type: ${CMAKE_BUILD_TYPE}")
    endif()
  endif()

  set(CMAKE_C_FLAGS_ASAN
      "${CMAKE_C_FLAGS_RelWithDebInfo} -fsanitize=address -fno-omit-frame-pointer" CACHE STRING
      "Flags used by the C compiler for Asan build type or configuration." FORCE)

  set(CMAKE_CXX_FLAGS_ASAN
      "${CMAKE_CXX_FLAGS_RelWithDebInfo} -fsanitize=address -fno-omit-frame-pointer" CACHE STRING
      "Flags used by the C++ compiler for Asan build type or configuration." FORCE)

  set(CMAKE_EXE_LINKER_FLAGS_ASAN
      "${CMAKE_EXE_LINKER_FLAGS_RelWithDebInfo} -fsanitize=address" CACHE STRING
      "Linker flags to be used to create executables for Asan build type." FORCE)

  set(CMAKE_SHARED_LINKER_FLAGS_ASAN
      "${CMAKE_SHARED_LINKER_FLAGS_RelWithDebInfo} -fsanitize=address" CACHE STRING
      "Linker lags to be used to create shared libraries for Asan build type." FORCE)
endif()

# macOS libc++ vs libstdc++ configuration
# See:
# https://stackoverflow.com/questions/52509602/cant-compile-c-program-on-a-mac-after-upgrade-to-mojave
# https://support.enthought.com/hc/en-us/articles/204469410-OS-X-GCC-Clang-and-Cython-in-10-9-Mavericks
# https://github.com/pandas-dev/pandas/pull/24274/files
# https://github.com/explosion/thinc/pull/84/files
# https://github.com/jlfaucher/builder/commit/d144d3a695949f90c5e2acff4dfd94fdcf8dcdfa
# https://github.com/CoolProp/CoolProp/issues/1778
# https://gitlab.kitware.com/cmake/cmake/issues/18396
if(DEFINED DARWIN_USE_LIBCPP)
  if(DARWIN_USE_LIBCPP)
    set(CMAKE_OSX_DEPLOYMENT_TARGET "10.9" CACHE STRING "Minimum OS X deployment version")
    set(OSX_COMPILE_FLAGS "${OSX_COMPILE_FLAGS} -stdlib=libc++")
    set(OSX_COMPILE_FLAGS "${OSX_COMPILE_FLAGS} -mmacosx-version-min=10.9")
    set(OSX_LINK_FLAGS "${OSX_LINK_FLAGS} -lc++")
    set(OSX_LINK_FLAGS "${OSX_LINK_FLAGS} -nodefaultlibs")
    set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
  else()
    set(CMAKE_OSX_DEPLOYMENT_TARGET "10.5" CACHE STRING "Minimum OS X deployment version")
    set(OSX_COMPILE_FLAGS "${OSX_COMPILE_FLAGS} -stdlib=libstdc++")
    set(OSX_COMPILE_FLAGS "${OSX_COMPILE_FLAGS} -mmacosx-version-min=10.5")
    set(OSX_LINK_FLAGS "${OSX_LINK_FLAGS} -lstdc++")
    set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libstdc++")
  endif()
  message(STATUS "DARWIN_USE_LIBCPP was set added some flags:")
  message(STATUS "  OSX_COMPILE_FLAGS: ${OSX_COMPILE_FLAGS}")
  message(STATUS "     OSX_LINK_FLAGS: ${OSX_LINK_FLAGS}")
else()
  if("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
    message(STATUS "OSX build detected:")
    message(STATUS "  You might want to pass the -DDARWIN_USE_LIBCPP=ON/OFF parameter")
    message(STATUS "  to enable or disable different C++ standard libraries.")
    message(STATUS "  You can also specify the environment variable MACOSX_DEPLOYMENT_TARGET=10.9 to force clang builds.")
  endif()
endif()

# Add definitions to silence warnings in MSVC2017 related to shared ptr code
if(MSVC AND MSVC_VERSION GREATER_EQUAL 1910)
  add_definitions(-D_SILENCE_CXX17_OLD_ALLOCATOR_MEMBERS_DEPRECATION_WARNING)
  add_definitions(-D_SILENCE_CXX17_RESULT_OF_DEPRECATION_WARNING)
endif()

# Handle MSVC release/debug configurations
if(COOLPROP_RELEASE AND COOLPROP_DEBUG)
  message(FATAL_ERROR "You cannot set both COOLPROP_RELEASE and COOLPROP_DEBUG")
endif()

if(COOLPROP_RELEASE)
  if(NOT MSVC)
    set(CMAKE_BUILD_TYPE "Release")
  else()
    # Multi-config generator (Visual Studio)
    # Can't set CMAKE_BUILD_TYPE, it's controlled by the IDE
    message(STATUS "Building in release mode with MSVC")
  endif()
endif()

# Handle MSVC static/dynamic runtime library linking
if(COOLPROP_MSVC_STATIC)
  set(COOLPROP_MSVC_ALL "/MTd" "/MT" "/MDd" "/MD")
  foreach(flag_var CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
                   CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
    foreach(TO_REPLACE ${COOLPROP_MSVC_ALL})
      string(REPLACE "${TO_REPLACE}" "" ${flag_var} "${${flag_var}}")
    endforeach()
    if(COOLPROP_MSVC_DEBUG)
      set(${flag_var} "${${flag_var}} /MTd")
    else()
      set(${flag_var} "${${flag_var}} /MT")
    endif()
  endforeach()
endif()
