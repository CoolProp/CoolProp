#include "catch.hpp"

int main (int argc, char * argv[]) {
    int result = Catch::Session().run( argc, argv );
    std::cout << "Result is:" << result << std::endl;
    return EXIT_SUCCESS;
}