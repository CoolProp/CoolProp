#include "CoolProp.h"
#include "AbstractState.h"
#include <iostream>
#include "crossplatform_shared_ptr.h"
using namespace CoolProp;
int main() {
    shared_ptr<AbstractState> Water(AbstractState::factory("HEOS", "Water"));
    Water->update(PQ_INPUTS, 101325, 0);  // SI units
    std::cout << "T: " << Water->T() << " K" << std::endl;
    std::cout << "rho': " << Water->rhomass() << " kg/m^3" << std::endl;
    std::cout << "rho': " << Water->rhomolar() << " mol/m^3" << std::endl;
    std::cout << "h': " << Water->hmass() << " J/kg" << std::endl;
    std::cout << "h': " << Water->hmolar() << " J/mol" << std::endl;
    std::cout << "s': " << Water->smass() << " J/kg/K" << std::endl;
    std::cout << "s': " << Water->smolar() << " J/mol/K" << std::endl;
    return EXIT_SUCCESS;
}