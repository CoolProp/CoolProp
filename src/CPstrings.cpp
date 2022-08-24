#include "CPstrings.h"
#include "crossplatform_shared_ptr.h"
#include <cstdio>
#include <vector>
#include <string>

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
        v.push_back(std::string(i1, i2));
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
    };                                                           // to use delete[]
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

#    include "crossplatform_shared_ptr.h"
#    include <catch2/catch_all.hpp>
#    include "CoolPropTools.h"
#    include "CoolProp.h"

TEST_CASE("Test endswith function", "[endswith]") {
    REQUIRE(endswith("aaa", "-PengRobinson") == false);
    REQUIRE(endswith("Ethylbenzene", "-PengRobinson") == false);
    REQUIRE(endswith("Ethylbenzene-PengRobinson", "-PengRobinson") == true);
    REQUIRE(endswith("Ethylbenzene", "Ethylbenzene") == true);
}

#endif
