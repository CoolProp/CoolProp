#include "CoolProp.h"
#include <iostream>
#include <stdlib.h>
using namespace CoolProp;
int main()
{
    std::cout << PropsSI("P", "D", 3.096728622, "U", 2043578.583, "n-Dodecane") << std::endl; // Breaks without Brent modification

    return EXIT_SUCCESS;
}
