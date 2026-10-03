
#ifndef COOLPROP_STRINGS_H
#define COOLPROP_STRINGS_H

#include <algorithm>
#include <cctype>
#include <system_error>
#include <vector>

#if !defined(NO_FMTLIB)
#    ifndef FMT_HEADER_ONLY
#        define FMT_HEADER_ONLY
#    endif
#    include "fmt/format.h"  // For addition of the string formatting functions and macros from fmtlib
#    include "fmt/printf.h"  // For sprintf
#    undef FMT_HEADER_ONLY
#else
#    include <vector>
#    include <string>
#endif

#include "CoolProp/Exceptions.h"

#if !defined(__powerpc__)
/// Copy string to wstring
/// Dangerous if the string has non-ASCII characters; from http://stackoverflow.com/a/8969776/1360263
inline void StringToWString(const std::string& s, std::wstring& ws) {
    ws = std::wstring(s.begin(), s.end());
}
#endif

/// The following code for the trim functions was taken from http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring
// trim from start
#ifdef HAS_MOVE_SEMANTICS  //More robust c++11 detection https://stackoverflow.com/questions/10717502/is-there-a-preprocessor-directive-for-detecting-c11x-support
// #ifdef __cplusplus >= 201103L || (defined(_MSC_VER) && _MSC_VER >= 1900)
inline std::string& strlstrip(std::string& s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}
#else
inline std::string& strlstrip(std::string& s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) { return !std::isspace(ch); }));
    return s;
}
#endif
// trim from end
#ifdef HAS_MOVE_SEMANTICS  //More robust c++11 detection https://stackoverflow.com/questions/10717502/is-there-a-preprocessor-directive-for-detecting-c11x-support
// #ifdef __cplusplus >= 201103L || (defined(_MSC_VER) && _MSC_VER >= 1900)
inline std::string& strrstrip(std::string& s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}
#else
inline std::string& strrstrip(std::string& s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) { return !std::isspace(ch); }).base(), s.end());
    return s;
}
#endif
// trim from both ends
inline std::string& strstrip(std::string& s) {
    return strlstrip(strrstrip(s));
}
/// Simple string function to check for end of string being equal to given string
inline bool endswith(const std::string& s1, const std::string& s2) {
    // Impossible to match a string longer than the given string
    if (s2.size() > s1.size()) {
        return false;
    }
    long lhs = static_cast<long>(s1.rfind(s2));
    long rhs = static_cast<long>(s1.size()) - static_cast<long>(s2.size());
    return lhs == rhs;
}

#if defined(NO_FMTLIB)
// Missing string formatting function, this old guy is needed for ancient gcc compilers on PowerPC for VxWorks
inline std::string format(const char* fmt, ...);
#elif FMT_VERSION >= 50000
template <typename... Args>
inline std::string format(const char* format_str, const Args&... args) {
    return fmt::sprintf(format_str, args...);
}
#else
inline std::string format(const char* format, fmt::ArgList args) {
    return fmt::sprintf(format, args);
}
FMT_VARIADIC(std::string, format, const char*)
#endif

// Missing string split - like in Python
std::vector<std::string> strsplit(const std::string& s, char del);

inline std::string upper(std::string str) {
    std::transform(str.begin(), str.end(), str.begin(), ::toupper);
    return str;
}

inline std::string lower(std::string str) {
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    return str;
}

std::string strjoin(const std::vector<std::string>& strings, const std::string& delim);

/// A convenience function that return true if a string begins with the given other string
inline bool strstartswith(const std::string& s, const std::string& other) {
    return s.find(other) == 0;
};

/// Parse a decimal number at the start of [first, last) with '.' as the decimal
/// separator, whatever the C locale.  strtod, atof and std::stod follow
/// LC_NUMERIC, so a host program that calls setlocale(LC_NUMERIC, "de_DE")
/// would make them stop at the '.'.
///
/// Accepts an optional leading '+' or '-', then the same decimal syntax as
/// strtod (digits [. digits] [(e|E)[+|-]digits], plus inf and nan).  Unlike
/// strtod it does not skip leading whitespace and does not accept hex.
///
/// Returns std::errc() and sets value and end on success.  Returns
/// std::errc::invalid_argument, with end == first, if no number starts at first.
/// Returns std::errc::result_out_of_range, with end past the number, if it
/// overflows or underflows a double; value is then unchanged.
std::errc parse_double_C(const char* first, const char* last, double& value, const char*& end);

/// Convert a whole string to a double, accepting a FORTRAN-style 'D' or 'd'
/// exponent and surrounding whitespace.  Independent of the C locale.  Throws
/// CoolProp::ValueError if any of the string is not part of the number.
double string2double(const std::string& s);

#endif
