#include "CoolProp.h"
using namespace CoolProp;
int main()
{
    // First type (slowest, most string processing, exposed in DLL)
    double r0A = PropsSI("Dmolar","T",298,"P",1e5,"Propane[0.5]&Ethane[0.5]"); // Default backend is HEOS
    double r0B = PropsSI("Dmolar","T",298,"P",1e5,"HEOS::Propane[0.5]&Ethane[0.5]");
    double r0C = PropsSI("Dmolar","T",298,"P",1e5,"REFPROP::Propane[0.5]&Ethane[0.5]");

    std::vector<double> z(2,0.5);
    // Second type (C++ only, a bit faster)
    double r1A = PropsSI("Dmolar","T",298,"P",1e5,"Propane&Ethane", z);
    double r1B = PropsSI("Dmolar","T",298,"P",1e5,"HEOS::Propane&Ethane", z);
    double r1C = PropsSI("Dmolar","T",298,"P",1e5,"REFPROP::Propane&Ethane", z);
    
    return EXIT_SUCCESS;
}