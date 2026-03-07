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

std::vector<std::string> strsplit_brace_aware(const std::string& s, char del) {
    std::vector<std::string> v;
    int depth = 0;
    std::string current;
    for (std::size_t i = 0; i < s.size(); ++i) {
        char c = s[i];
        if (c == '{') {
            ++depth;
            current += c;
        } else if (c == '}') {
            if (depth == 0) {
                throw CoolProp::ValueError("strsplit_brace_aware: unmatched '}' in string: " + s);
            }
            --depth;
            current += c;
        } else if (c == del && depth == 0) {
            v.push_back(current);
            current.clear();
        } else {
            current += c;
        }
    }
    v.push_back(current);
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

TEST_CASE("Test strsplit_brace_aware", "[brace_aware]") {
    SECTION("No braces — behaves like strsplit") {
        auto v = strsplit_brace_aware("Water&Ethanol", '&');
        REQUIRE(v.size() == 2);
        REQUIRE(v[0] == "Water");
        REQUIRE(v[1] == "Ethanol");
    }
    SECTION("& inside JSON value is not a split point") {
        auto v = strsplit_brace_aware("Water{\"k\":\"a&b\"}&Ethanol", '&');
        REQUIRE(v.size() == 2);
        REQUIRE(v[0] == "Water{\"k\":\"a&b\"}");
        REQUIRE(v[1] == "Ethanol");
    }
    SECTION("Bare JSON token (mixture-level config)") {
        auto v = strsplit_brace_aware("Water&Ethanol&{\"mixing_rule\":\"HV\"}", '&');
        REQUIRE(v.size() == 3);
        REQUIRE(v[2] == "{\"mixing_rule\":\"HV\"}");
    }
    SECTION("Single fluid with JSON suffix") {
        auto v = strsplit_brace_aware("Water{\"EOS\":\"Wagner-JPCRD-2002\"}", '&');
        REQUIRE(v.size() == 1);
        REQUIRE(v[0] == "Water{\"EOS\":\"Wagner-JPCRD-2002\"}");
    }
    SECTION("Unmatched } throws") {
        REQUIRE_THROWS_AS(strsplit_brace_aware("Water}", '&'), CoolProp::ValueError);
    }
}

#endif
