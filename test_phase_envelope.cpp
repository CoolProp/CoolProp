#include "AbstractState.h"
#include "CoolProp.h"
#include <iostream>
#include <iomanip>

using namespace CoolProp;

int main() {
    try {
        std::cout << "Testing mass quality calculations with R410A.MIX at 278.15 K and 0.5 quality" << std::endl;
        std::cout << "=================================================" << std::endl;
        
        // Create abstract state for R410A.MIX
        shared_ptr<AbstractState> AS(AbstractState::factory("REFPROP", "R410A.MIX"));
        
        AS->update(QT_INPUTS, 0.5, 278.15);  // Saturated vapor (q=1 in molar basis)
        double p_sat = AS->p();
        std::cout << "  Saturation pressure (molar basis): " << std::setprecision(12) << AS->p() << " Pa" << std::endl;        
        std::cout << "  Mole density: " << std::setprecision(12) << AS->rhomolar()/1e3 << " kmol/m3" << std::endl;
        std::cout << "  Mass density: " << std::setprecision(12) << AS->rhomass() << " kg/m3" << std::endl;
        std::cout << "  Molar quality: " << std::setprecision(12) << AS->Q() << " -" << std::endl;
        std::cout << "  Mass quality: " << std::setprecision(12) << AS->Qmass() << " -" << std::endl;
        auto mf = AS->get_mole_fractions(); std::cout << "  Mole fractions: "; for (const auto& x : mf) std::cout << std::setprecision(12) << x << " "; std::cout << std::endl;
        auto mf_mass = AS->get_mass_fractions(); std::cout << "  Mass fractions: "; for (const auto& x : mf_mass) std::cout << std::setprecision(12) << x << " "; std::cout << std::endl;

        AS->update(QmassT_INPUTS, 0.5, 278.15);  // Saturated vapor (q=1 in mass basis)
        std::cout << "  Saturation pressure (mass basis): " << std::setprecision(12) << AS->p() << " Pa" << std::endl;        
        std::cout << "  Mole density: " << std::setprecision(12) << AS->rhomolar()/1e3 << " kmol/m3" << std::endl;
        std::cout << "  Mass density: " << std::setprecision(12) << AS->rhomass() << " kg/m3" << std::endl;
        std::cout << "  Molar quality: " << std::setprecision(12) << AS->Q() << " -" << std::endl;
        std::cout << "  Mass quality: " << std::setprecision(12) << AS->Qmass() << " -" << std::endl;

        AS->update(PQ_INPUTS, 1e5, 0.5);  // Saturated vapor (q=1 in mass basis)
        std::cout << "  Saturation temperature (molar basis): " << std::setprecision(12) << AS->T() << " K" << std::endl;        
        std::cout << "  Mole density: " << std::setprecision(12) << AS->rhomolar()/1e3 << " kmol/m3" << std::endl;
        std::cout << "  Mass density: " << std::setprecision(12) << AS->rhomass() << " kg/m3" << std::endl;
        std::cout << "  Molar quality: " << std::setprecision(12) << AS->Q() << " -" << std::endl;
        std::cout << "  Mass quality: " << std::setprecision(12) << AS->Qmass() << " -" << std::endl;

        AS->update(PQmass_INPUTS, 1e5, 0.5);  // Saturated vapor (q=1 in mass basis)
        std::cout << "  Saturation temperature (mass basis): " << std::setprecision(12) << AS->T() << " K" << std::endl;        
        std::cout << "  Mole density: " << std::setprecision(12) << AS->rhomolar()/1e3 << " kmol/m3" << std::endl;
        std::cout << "  Mass density: " << std::setprecision(12) << AS->rhomass() << " kg/m3" << std::endl;
        std::cout << "  Molar quality: " << std::setprecision(12) << AS->Q() << " -" << std::endl;
        std::cout << "  Mass quality: " << std::setprecision(12) << AS->Qmass() << " -" << std::endl;

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    }
}


