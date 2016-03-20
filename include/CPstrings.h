
#ifndef COOLPROP_STRINGS_H
#define COOLPROP_STRINGS_H

    #include <iterator>
    #include <algorithm>
    #include <functional>

    #define FMT_HEADER_ONLY
    #include "externals/cppformat/cppformat/format.h" // For addition of the string formatting functions and macros from cppformat
    #undef FMT_HEADER_ONLY

    #if !defined(__powerpc__)
    /// Copy string to wstring
    /// Dangerous if the string has non-ASCII characters; from http://stackoverflow.com/a/8969776/1360263 
    inline void StringToWString(const std::string &s, std::wstring &ws)
    {
        ws = std::wstring(s.begin(), s.end());
    }
    #endif

    /// The following code for the trim functions was taken from http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring
    // trim from start
    inline std::string &strlstrip(std::string &s) {
            s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
            return s;
    }
    // trim from end
    inline std::string &strrstrip(std::string &s) {
            s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
            return s;
    }
    // trim from both ends
    inline std::string &strstrip(std::string &s) {
            return strlstrip(strrstrip(s));
    }
    /// Simple string function to check for end of string being equal to given string
    inline bool endswith(const std::string &s1, const std::string &s2){
        long lhs = static_cast<long>(s1.rfind(s2));
        long rhs = static_cast<long>(s1.size()) - static_cast<long>(s2.size());
        return lhs == rhs;
    }

    // Missing std::string formatting function - provided by the cppformat library
    inline std::string format(const char *format, fmt::ArgList args) {
      return fmt::sprintf(format, args);
    }
    FMT_VARIADIC(std::string, format, const char *)

    // Missing string split - like in Python
    std::vector<std::string> strsplit(const std::string &s, char del);

    inline std::string upper(std::string str)
    {
        std::transform(str.begin(), str.end(), str.begin(), ::toupper);
        return str;
    }
	
	inline std::string lower(std::string str)
    {
        std::transform(str.begin(), str.end(), str.begin(), ::tolower);
        return str;
    }

    std::string strjoin(const std::vector<std::string> &strings, const std::string &delim);

#endif