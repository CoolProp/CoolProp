#include <chrono>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <vector>
#include <type_traits>
#include <CoolPropLib.h>

using Clock = std::chrono::high_resolution_clock;
using time_point = std::chrono::time_point<Clock>;

//template <typename T>
//std::enable_if_t<is_cont<T>::value, std::string>
//print_vector(const T& container) {
//    std::string str = "{";
//    for (auto it = std::begin(container); it != std::end(container); ++it) {
//        str += print_vector(*it);
//        if (isLast(it, container)) {
//            str += '}';
//        } else {
//            str += ',';
//        }
//    }
//
//    return str;
//}
//
//template <typename T>
//std::enable_if_t<!is_cont<T>::value, std::string>
//print_vector(const T& value) {
//    return std::to_string(value);
//}

template <typename T>
std::string _print(const T& item) {
    std::stringstream stm;
    stm << item;
    return stm.str();
}

template <typename T>
std::string _print(const std::pair<T,T> &items) {
    std::stringstream stm;
    stm << _print(items.first) << " - " << _print(items.second);
    return stm.str();
}

template <typename T>
std::string _print(const std::vector<T> &items) {
    std::stringstream stm;
    stm << '{';
    for (auto it = std::begin(items); it != std::end(items); ++it) {
        stm << _print(*it);
        if (it == std::end(items)-1) {
            stm << '}';
        } else {
            stm << ',';
        }
    }
    return stm.str();
}




#define RP_LN 1000

bool get_refprop_fluids(std::vector<std::vector<std::string>>& fluids) {
    long ret = 0;
    char rp_name[RP_LN];
    int rp_length;
    for (auto& f_v : fluids) { 
        ret = get_fluid_param_string(f_v[0].c_str(), "REFPROP_name", &rp_name[0], RP_LN);
        if (ret == 1) {
            f_v.push_back(rp_name);
        } else {
            return false;
        }
    }
    return true;
}

std::vector<std::pair<std::string,std::string>> get_state_strings() {
    std::vector<std::pair<std::string, std::string>> string_vector;
    //std::vector<std::string> backends = { "HEOS", "REFPROP" };
    std::vector<std::string> fluids = {"Water", "Ammonia", "R134a", "R290", "R407c" };

    long ret = 0;
    char rp_name[RP_LN];
    int rp_length;
    for (auto& f_v : fluids) {
        string_vector.push_back({ { "HEOS" }, { f_v } });
        ret = get_fluid_param_string(f_v.c_str(), "REFPROP_name", &rp_name[0], RP_LN);
        if (ret == 1) {
            string_vector.push_back({{ "REFPROP" }, { std::string(rp_name) }});
        }
    }
    return string_vector;
}

std::vector<long> get_pointers() {
    //EXPORT_CODE long CONVENTION AbstractState_factory(const char* backend, const char* fluids, long *errcode, char *message_buffer, const long buffer_length);
    return { 1 };
}

std::vector<std::string> strings = {"one", "two", "three", "four", "five"};

std::chrono::duration<double> vector_push_back(const size_t n) {
    time_point start, end;
    start = Clock::now();
    
    std::vector<std::string> v;
    for (size_t i = 0; i < n; ++i) {
        v.push_back(strings[i % strings.size()]);
    }
    
    end = Clock::now();
    return end - start;
}

std::chrono::duration<double> vector_push_back_with_reserve(const size_t n) {
    time_point start, end;
    start = Clock::now();
    
    std::vector<std::string> v;
    v.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        v.push_back(strings[i % strings.size()]);
    }
    
    end = Clock::now();
    return end - start;
}

std::chrono::duration<double> vector_element_assignment(const size_t n) {
    time_point start, end;
    start = Clock::now();
    
    std::vector<std::string> v(n);
    for (size_t i = 0; i < n; ++i) {
        v[i] = strings[i % strings.size()];
    }
    
    end = Clock::now();
    return end - start;
}

std::chrono::duration<double> vector_emplace_back(const size_t n) {
    time_point start, end;
    start = Clock::now();
    
    std::vector<std::string> v;
    for (size_t i = 0; i < n; ++i) {
        v.emplace_back(strings[i % strings.size()]);
    }
    
    end = Clock::now();
    return end - start;
}

std::chrono::duration<double> vector_emplace_back_with_reserve(const size_t n) {
    time_point start, end;
    start = Clock::now();
    
    std::vector<std::string> v;
    v.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        v.emplace_back(strings[i % strings.size()]);
    }
    
    end = Clock::now();
    return end - start;
}

int main() {
    //const size_t n = 10000;
    //std::cout << "vector push_back: " << vector_push_back(n).count() << "\n";
    //std::cout << "vector push_back with reserve: " << vector_push_back(n).count() << "\n";
    //std::cout << "vector element assignment: " << vector_element_assignment(n).count() << "\n";
    //std::cout << "vector emplace_back: " << vector_emplace_back(n).count() << "\n";
    //std::cout << "vector emplace_back with reserve: " << vector_emplace_back_with_reserve(n).count() << "\n";
    //std::cout << "Fetching REFPROP names: " << get_refprop_fluids(fluids);
    std::cout << "Fluid names: " << _print(get_state_strings());
    return 1;
}
