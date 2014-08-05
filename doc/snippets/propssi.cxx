#include "CoolProp.h"
#include <iostream>
using namespace CoolProp;
int main()
{
    // First type (slowest, due to most string processing, exposed in DLL)
    std::cout << PropsSI("Dmolar","T",298,"P",1e5,"Propane[0.5]&Ethane[0.5]") << std::endl; // Default backend is HEOS
    std::cout << PropsSI("Dmolar","T",298,"P",1e5,"HEOS::Propane[0.5]&Ethane[0.5]") << std::endl;
    std::cout << PropsSI("Dmolar","T",298,"P",1e5,"REFPROP::Propane[0.5]&Ethane[0.5]") << std::endl;

    std::vector<double> z(2,0.5);
    // Second type (C++ only, a bit faster)
    std::cout << PropsSI("Dmolar","T",298,"P",1e5,"Propane&Ethane", z) << std::endl;
    std::cout << PropsSI("Dmolar","T",298,"P",1e5,"HEOS::Propane&Ethane", z) << std::endl;
    std::cout << PropsSI("Dmolar","T",298,"P",1e5,"REFPROP::Propane&Ethane", z) << std::endl;
    
    return EXIT_SUCCESS;
}