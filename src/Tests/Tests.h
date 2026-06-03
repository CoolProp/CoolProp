#ifndef TESTS_H
#define TESTS_H

#include <vector>
#include <string>

void run_tests();
int run_fast_tests();
int run_not_slow_tests();
int run_user_defined_tests(const std::vector<std::string>& tests_or_tags);

#endif