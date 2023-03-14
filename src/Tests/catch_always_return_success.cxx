// This file is used in tests to always return success
// so that address sanitizer's return code will be the
// return code that is fed back to buildbot, not the
// number of failing catch tests

#define CATCH_CONFIG_RUNNER
#include <catch2/catch_all.hpp>
#include <iostream>

int main(int argc, char* argv[]) {
    int result = Catch::Session().run(argc, argv);
    std::cout << "Result is:" << result << std::endl;
    return EXIT_SUCCESS;
}
