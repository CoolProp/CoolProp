// Compiles Boost.CharConv's from_chars into CoolProp.  The library is not
// header-only for floating point, and CoolProp ships no separately built Boost
// libraries, so its source comes from the boost-headers subset (which carries
// libs/charconv/src for this purpose) and is built here.  The lexer in
// Expression.cpp uses it to read number literals without regard to the C locale.
//
// Living in src/expression/ means every build that globs CoolProp's sources
// (the main CMakeLists.txt and wrappers/Python/CMakeLists.txt) picks it up.
#include "libs/charconv/src/from_chars.cpp"  // NOLINT(bugprone-suspicious-include)
