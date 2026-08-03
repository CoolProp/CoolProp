#include "AbstractState.h"
#include "CoolProp.h"
#include <iostream>
#include <iomanip>

using namespace CoolProp;

int main() {
    std::cout << "Testing mass quality calculations with R410A.MIX at 278.15 K and 0.5 quality" << std::endl;
    std::cout << "=================================================" << std::endl;

    try {
        shared_ptr<AbstractState> AS(AbstractState::factory("REFPROP", "R410A.MIX"));

        AS->update(QT_INPUTS, 0.5, 278.15);
        std::cout << "  Viscosity: " << std::setprecision(12) << AS->viscosity() << " Pa.s" << std::endl;
        std::cout << "  Pressure: " << std::setprecision(12) << AS->p() << " Pa" << std::endl;

        AS->update(QT_INPUTS, 0.0, 278.15);
        std::cout << "  Viscosity from QT (0.0, 278.15): " << std::setprecision(12) << AS->viscosity() << " Pa.s" << std::endl;

        AS->update(QT_INPUTS, 1.0, 278.15);
        std::cout << "  Viscosity from QT (1.0, 278.15): " << std::setprecision(12) << AS->viscosity() << " Pa.s" << std::endl;

        AS->update(QmassT_INPUTS, 0.0, 278.15);
        std::cout << "  Viscosity from QmassT (0.0, 278.15): " << std::setprecision(12) << AS->viscosity() << " Pa.s" << std::endl;

        AS->update(QmassT_INPUTS, 1.0, 278.15);
        std::cout << "  Viscosity from QmassT (1.0, 278.15): " << std::setprecision(12) << AS->viscosity() << " Pa.s" << std::endl;

        AS->update(QmassT_INPUTS, 0.5, 278.15);
        std::cout << "  Viscosity from QmassT (0.5, 278.15): " << std::setprecision(12) << AS->viscosity() << " Pa.s" << std::endl;

        AS->update(QT_INPUTS, 0.5, 278.15);
        std::cout << "  Viscosity from QT (0.5, 278.15): " << std::setprecision(12) << AS->viscosity() << " Pa.s" << std::endl;

        AS->update(QmassT_INPUTS, 0.5, 278.15);
        std::cout << "  Viscosity from QmassT (0.5, 278.15): " << std::setprecision(12) << AS->viscosity() << " Pa.s" << std::endl;

        // -----------------------------------------

        AS->update(PQ_INPUTS, 9e5, 0.5);
        std::cout << "  Viscosity: " << std::setprecision(12) << AS->viscosity() << " Pa.s" << std::endl;

        AS->update(PQ_INPUTS, 9e5, 0.0);
        std::cout << "  Viscosity from PQ (9e5, 0.0): " << std::setprecision(12) << AS->viscosity() << " Pa.s" << std::endl;

        AS->update(PQ_INPUTS, 9e5, 1.0);
        std::cout << "  Viscosity from PQ (9e5, 1.0): " << std::setprecision(12) << AS->viscosity() << " Pa.s" << std::endl;

        AS->update(PQmass_INPUTS, 9e5, 0.0);
        std::cout << "  Viscosity from PQmass (9e5, 0.0): " << std::setprecision(12) << AS->viscosity() << " Pa.s" << std::endl;

        AS->update(PQmass_INPUTS, 9e5, 1.0);
        std::cout << "  Viscosity from PQmass (9e5, 1.0): " << std::setprecision(12) << AS->viscosity() << " Pa.s" << std::endl;

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    }
}
