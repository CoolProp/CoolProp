// This file is used in tests to always return success
// so that address saniter's return code will be the 
// return code that is fed back to buildbot, not the
// number of failing catch tests
#include "catch.hpp"

int main (int argc, char * argv[]) {
    int result = Catch::Session().run( argc, argv );
    std::cout << "Result is:" << result << std::endl;
    return EXIT_SUCCESS;
}