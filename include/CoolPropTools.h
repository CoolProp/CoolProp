#ifndef COOLPROPTOOLS_H
#define COOLPROPTOOLS_H

#ifndef _CRT_SECURE_NO_WARNINGS
#    define _CRT_SECURE_NO_WARNINGS
#endif

#include "PlatformDetermination.h"
#include "Exceptions.h"
#include <string>
#include <vector>
#include <cctype>
#include <map>

#include "CPstrings.h"
#include "CPnumerics.h"
#include "CPfilepaths.h"

#ifndef __has_feature           // Optional of course.
#    define __has_feature(x) 0  // Compatibility with non-clang compilers.
#endif

#ifdef __EMSCRIPTEN__
#    define thread_local
#endif

// see http://stackoverflow.com/questions/18298280/how-to-declare-a-variable-as-thread-local-portably
#ifndef thread_local
#    if __STDC_VERSION__ >= 201112 && !defined __STDC_NO_THREADS__
#        define thread_local _Thread_local
#    elif defined _WIN32 && (defined _MSC_VER || defined __ICL || defined __DMC__ || defined __BORLANDC__)
#        define thread_local __declspec(thread)
#    elif defined(__ISAPPLE__) && (defined(__llvm__) || defined(__clang__)) && !__has_feature(cxx_thread_local)
#        define thread_local
/* note that ICC (linux) and Clang are covered by __GNUC__ */
#    elif defined __GNUC__ || defined __SUNPRO_C || defined __xlC__
#        define thread_local __thread
#    else
#        error "Cannot define thread_local"
//          #define thread_local
#    endif
#endif

#define COOLPROPDBL_MAPS_TO_DOUBLE
#ifdef COOLPROPDBL_MAPS_TO_DOUBLE
typedef double CoolPropDbl;
#else
typedef long double CoolPropDbl;
#endif

/// Define the deprecated macro to give compile-time warnings
#ifdef __GNUC__
#    define DEPRECATED(func) func __attribute__((deprecated))
#elif defined(_MSC_VER)
#    define DEPRECATED(func) __declspec(deprecated) func
#else
#    pragma message("WARNING: You need to implement DEPRECATED for this compiler")
#    define DEPRECATED(func) func
#endif

class Dictionary
{
   private:
    typedef std::map<std::string, double> numbers_map;
    numbers_map numbers;
    typedef std::map<std::string, std::string> strings_map;
    strings_map strings;
    typedef std::map<std::string, std::vector<double>> double_vectors_map;
    double_vectors_map double_vectors;
    typedef std::map<std::string, std::vector<std::string>> string_vectors_map;
    string_vectors_map string_vectors;

   public:
    Dictionary(){};
    bool is_empty(void) const {
        return numbers.empty() && strings.empty() && double_vectors.empty() && string_vectors.empty();
    }
    void add_string(const std::string& s1, const std::string& s2) {
        strings.insert(std::pair<std::string, std::string>(s1, s2));
    }
    void add_number(const std::string& s1, double d) {
        numbers.erase(s1);
        numbers.insert(std::pair<std::string, double>(s1, d));
    }
    bool has_number(const std::string& s1) {
        return numbers.find(s1) != numbers.end();
    }
    void add_double_vector(const std::string& s1, const std::vector<double>& d) {
        double_vectors.insert(std::pair<std::string, std::vector<double>>(s1, d));
    }
    void add_string_vector(const std::string& s1, const std::vector<std::string>& d) {
        string_vectors.insert(std::pair<std::string, std::vector<std::string>>(s1, d));
    }
    std::string get_string(const std::string& s) const {
        strings_map::const_iterator i = strings.find(s);
        if (i != strings.end()) {
            return i->second;
        } else {
            throw CoolProp::ValueError(format("%s could not be matched in get_string", s.c_str()));
        }
    };
    double get_double(const std::string& s) const {
        numbers_map::const_iterator i = numbers.find(s);
        if (i != numbers.end()) {
            return i->second;
        } else {
            throw CoolProp::ValueError(format("%s could not be matched in get_number", s.c_str()));
        }
    };
    /// Get a double, or return the default value if not found
    double get_double(const std::string& s, const double default_value) const {
        numbers_map::const_iterator i = numbers.find(s);
        if (i != numbers.end()) {
            return i->second;
        } else {
            return default_value;
        }
    };
    double get_number(const std::string& s) const {
        return get_double(s);
    };
    const std::vector<double>& get_double_vector(const std::string& s) const {
        double_vectors_map::const_iterator i = double_vectors.find(s);
        if (i != double_vectors.end()) {
            return i->second;
        } else {
            throw CoolProp::ValueError(format("%s could not be matched in get_double_vector", s.c_str()));
        }
    };
    const std::vector<std::string>& get_string_vector(const std::string& s) const {
        string_vectors_map::const_iterator i = string_vectors.find(s);
        if (i != string_vectors.end()) {
            return i->second;
        } else {
            throw CoolProp::ValueError(format("%s could not be matched in get_string_vector", s.c_str()));
        }
    };
};
/// Utility function to clear a std::map of pointers
//http://stackoverflow.com/questions/569110/why-is-memory-still-accessible-after-stdmapclear-is-called
template <typename M>
void freeClear(M& amap) {
    for (typename M::iterator it = amap.begin(); it != amap.end(); ++it) {
        delete it->second;
    }
    amap.clear();
}

#define CATCH_ALL_ERRORS_RETURN_HUGE(x) \
    try {                               \
        x                               \
    } catch (const std::exception& e) { \
        return _HUGE;                   \
    } catch (...) {                     \
        return _HUGE;                   \
    }

enum miniz_mode
{
    MINIZ_COMPRESS,
    MINIZ_DECOMPRESS
};
void miniz(const std::string& inFile, const std::string& outFile, miniz_mode mode);
#endif
