#######################################
#           BUILD OPTIONS             #
#-------------------------------------#
# These options are available to be   #
# modified in the build process.      #
# packages may want to modify these   #
# to suit, or just leave as defaults. #
#######################################

# Library build types
option(COOLPROP_STATIC_LIBRARY "Build CoolProp as a static library (.lib, .a)" OFF)
option(COOLPROP_SHARED_LIBRARY "Build CoolProp as a shared library (.dll, .so)" OFF)
option(COOLPROP_OBJECT_LIBRARY "Build CoolProp objects, but do not link them (.obj, .o)" OFF)

# Platform-specific packages
option(COOLPROP_WINDOWS_PACKAGE "Build the Windows installer" OFF)
option(COOLPROP_DEBIAN_PACKAGE "Build Debian package" OFF)
option(COOLPROP_IOS_TARGET "Build for iOS" OFF)

# Build configuration
option(BUILD_TESTING "Enable testing for this given builder" OFF)
option(COOLPROP_RELEASE "Optimize the builds with the release specs" OFF)
option(COOLPROP_DEBUG "Make a debug build" OFF)
option(COOLPROP_NO_EXAMPLES "Do not generate example code, does only apply to some wrappers." OFF)

# Bitness options
option(FORCE_BITNESS_32 "Force a 32bit build regardless of the host" OFF)
option(FORCE_BITNESS_64 "Force a 64bit build regardless of the host" OFF)
option(FORCE_BITNESS_NATIVE "Force a native bitness build regardless of the host" OFF)

# MSVC-specific options
option(COOLPROP_MSVC_STATIC "Statically link Microsoft Standard library removes dependency on MSVCRXXX.dll." OFF)
option(COOLPROP_MSVC_DYNAMIC "Dynamically link Microsoft Standard library to integrate with other builds." OFF)
option(COOLPROP_MSVC_DEBUG "Link the debug version of Microsoft Standard library to the debug builds." ON)

# Library calling conventions (32-bit only)
option(COOLPROP_STDCALL_LIBRARY "Build CoolProp as a 32bit shared library with stdcall" OFF)
option(COOLPROP_CDECL_LIBRARY "Build CoolProp as a 32bit shared library with cdecl" OFF)
option(COOLPROP_EXTERNC_LIBRARY "Overwrite the export settings to force extern C" OFF)

# Language wrapper modules
option(COOLPROP_EES_MODULE "Build the EES module" OFF)
option(COOLPROP_PRIME_MODULE "Build Mathcad Prime module" OFF)
option(COOLPROP_MATHCAD15_MODULE "Build Mathcad 15 module" OFF)
option(COOLPROP_OCTAVE_MODULE "Build Octave module" OFF)
option(COOLPROP_PYTHON_MODULE "Build Python module" OFF)
option(COOLPROP_CSHARP_MODULE "Build C# module" OFF)
option(COOLPROP_VBDOTNET_MODULE "Build VB.NET module" OFF)
option(COOLPROP_R_MODULE "Build R module" OFF)
option(COOLPROP_JAVA_MODULE "Build Java module" OFF)
option(COOLPROP_ANDROID_MODULE "Build Android module" OFF)
option(COOLPROP_PHP_MODULE "Build PHP module" OFF)
option(COOLPROP_LIBREOFFICE_MODULE "Build LibreOffice extension" OFF)
option(COOLPROP_JAVASCRIPT_MODULE "Build Javascript module" OFF)
option(COOLPROP_MATHEMATICA_MODULE "Build Mathematica module" OFF)
option(COOLPROP_SMATH_MODULE "Build SMath module" OFF)
option(COOLPROP_SMATH_WORK_INPLACE "Build SMath wrapper in source directory" OFF)

# Example and test modules
option(COOLPROP_MAIN_MODULE "Build main example" OFF)
option(COOLPROP_CATCH_MODULE "Build Catch test runner" OFF)
option(COOLPROP_CPP_EXAMPLE_TEST "Build C++ example test" OFF)
option(COOLPROP_SNIPPETS "Build code snippets" OFF)

# Development options
option(COOLPROP_VXWORKS_MAKEFILE "Generate VxWorks makefile" OFF)
option(COOLPROP_VXWORKS_LIBRARY_MODULE "Build VxWorks library module" OFF)
option(COOLPROP_VXWORKS_LIBRARY "Build VxWorks library" OFF)
option(COOLPROP_FPIC "Add -fPIC flag" OFF)
option(COOLPROP_IWYU "Enable include-what-you-use" OFF)
option(COOLPROP_ASAN "Enable address sanitizer" OFF)
option(COOLPROP_CLANG_ADDRESS_SANITIZER "Build with Clang address sanitizer" OFF)
option(COOLPROP_LAZY_LOAD_SUPERANCILLARIES "Lazy load superancillaries" OFF)
