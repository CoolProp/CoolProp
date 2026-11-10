// Compiles Boost.CharConv's from_chars into CoolProp.  The library is not
// header-only for floating point, and CoolProp ships no separately built Boost
// libraries, so its source comes from the boost-headers subset (which carries
// libs/charconv/src for this purpose) and is built here.  parse_double_C() in
// CPstrings.cpp wraps it to read numbers without regard to the C locale.
//
// Living in src/ means every build that globs CoolProp's sources
// (the main CMakeLists.txt and wrappers/Python/CMakeLists.txt) picks it up.
//
// BOOST_CHARCONV_SOURCE is what Boost's own build defines when compiling this
// file.  Without it, boost/charconv/config.hpp makes MSVC auto-link a
// boost_charconv .lib that CoolProp does not ship, and the link fails.
#define BOOST_CHARCONV_SOURCE
#include "libs/charconv/src/from_chars.cpp"  // NOLINT(bugprone-suspicious-include)
