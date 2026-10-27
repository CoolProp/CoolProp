#include "CoolProp/detail/strings.h"
#include <cctype>
#include <memory>
using std::shared_ptr;
#include <vector>
#include <string>

// Boost.CharConv's from_chars is compiled into CoolProp by src/boost_charconv.cpp,
// so suppress MSVC's auto-link to a separately built boost_charconv library.
#ifndef BOOST_CHARCONV_NO_LIB
#    define BOOST_CHARCONV_NO_LIB
#endif
#include <boost/charconv/from_chars.hpp>

std::string strjoin(const std::vector<std::string>& strings, const std::string& delim) {
    // Empty input vector
    if (strings.empty()) {
        return "";
    }

    std::string output = strings[0];
    for (unsigned int i = 1; i < strings.size(); i++) {
        output += format("%s%s", delim.c_str(), strings[i].c_str());
    }
    return output;
}

std::vector<std::string> strsplit(const std::string& s, char del) {
    std::vector<std::string> v;
    std::string::const_iterator i1 = s.begin(), i2;
    while (true) {
        i2 = std::find(i1, s.end(), del);
        v.emplace_back(i1, i2);
        if (i2 == s.end()) break;
        i1 = i2 + 1;
    }
    return v;
}

#if defined(NO_FMTLIB)
std::string format(const char* fmt, ...) {
    const int size = 512;
    struct deleter
    {
        static void delarray(char* p) {
            delete[] p;
        }
    };  // to use delete[]
    shared_ptr<char> buffer(new char[size], deleter::delarray);  // I'd prefer unique_ptr, but it's only available since c++11
    va_list vl;
    va_start(vl, fmt);
    int nsize = vsnprintf(buffer.get(), size, fmt, vl);
    if (size <= nsize) {                                     //fail delete buffer and try again
        buffer.reset(new char[++nsize], deleter::delarray);  //+1 for /0
        nsize = vsnprintf(buffer.get(), nsize, fmt, vl);
    }
    va_end(vl);
    return buffer.get();
}
#endif

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>
#    include "CoolProp/detail/tools.h"
#    include "CoolProp/CoolProp.h"

TEST_CASE("Test endswith function", "[endswith]") {
    REQUIRE(endswith("aaa", "-PengRobinson") == false);
    REQUIRE(endswith("Ethylbenzene", "-PengRobinson") == false);
    REQUIRE(endswith("Ethylbenzene-PengRobinson", "-PengRobinson") == true);
    REQUIRE(endswith("Ethylbenzene", "Ethylbenzene") == true);
}

#endif

std::errc parse_double_C(const char* first, const char* last, double& value, const char*& end) {
    end = first;
    const char* p = first;
    // from_chars takes '-' but not '+'; strtod took both.
    if (p < last && *p == '+') {
        ++p;
        if (p < last && *p == '-') return std::errc::invalid_argument;
    }
    const boost::charconv::from_chars_result r = boost::charconv::from_chars(p, last, value);
    if (r.ec == std::errc() || r.ec == std::errc::result_out_of_range) end = r.ptr;
    return r.ec;
}

double string2double(const std::string& s) {
    std::string mys = s;  // copy
    // Replace the first D or d with e (FORTRAN-style exponent)
    std::size_t pos = mys.find('D');
    if (pos != std::string::npos) mys.replace(pos, 1, "e");
    pos = mys.find('d');
    if (pos != std::string::npos) mys.replace(pos, 1, "e");

    const char* first = mys.c_str();
    const char* last = first + mys.size();
    // Skip surrounding whitespace, e.g. the trailing space left by Windows'
    // "set COOLPROP_X=0.25 && ..." idiom.
    while (first < last && std::isspace(static_cast<unsigned char>(*first)) != 0)
        ++first;
    while (last > first && std::isspace(static_cast<unsigned char>(last[-1])) != 0)
        --last;
    double val = 0.0;
    const char* end = nullptr;
    const std::errc ec = parse_double_C(first, last, val, end);
    if (ec == std::errc::result_out_of_range) {
        throw CoolProp::ValueError(format("Number is out of range for a double:%s", s.c_str()));
    }
    if (ec != std::errc() || end != last) {
        // Found a character that is not able to be converted to number
        throw CoolProp::ValueError(format("Unable to convert this string to a number:%s", s.c_str()));
    }
    return val;
}
