#include "CoolProp.h"
#include <iostream>
#include <stdlib.h>
using namespace CoolProp;
int main() {
    // First type (slowest, due to most string processing, exposed in DLL)
    std::cout << PropsSI("Dmolar", "T", 298, "P", 1e5, "Propane[0.5]&Ethane[0.5]") << std::endl;  // Default backend is HEOS
    std::cout << PropsSI("Dmolar", "T", 298, "P", 1e5, "HEOS::Propane[0.5]&Ethane[0.5]") << std::endl;
    return EXIT_SUCCESS;
}