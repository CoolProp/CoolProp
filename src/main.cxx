#include <chrono>
#include <iostream>
#include <vector>

using Clock = std::chrono::high_resolution_clock;
using time_point = std::chrono::time_point<Clock>;

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
    const size_t n = 10000;
    std::cout << "vector push_back: " << vector_push_back(n).count() << "\n";
    std::cout << "vector push_back with reserve: " << vector_push_back(n).count() << "\n";
    std::cout << "vector element assignment: " << vector_element_assignment(n).count() << "\n";
    std::cout << "vector emplace_back: " << vector_emplace_back(n).count() << "\n";
    std::cout << "vector emplace_back with reserve: " << vector_emplace_back_with_reserve(n).count() << "\n";
}
