#include "AbstractState.h"
#include "CoolProp.h"
#include <iostream>
#include <iomanip>

using namespace CoolProp;

int main() {
    try {
        std::cout << "Testing phase envelope for R410A.MIX at 278.15 K" << std::endl;
        std::cout << "=================================================" << std::endl;
        
        // Create abstract state for R410A.MIX
        shared_ptr<AbstractState> AS(AbstractState::factory("REFPROP", "R410A.MIX"));
        
        // Set quality flag: 1 = molar basis (default), 2 = mass basis
        //AS->set_quality_flag(1);  // Using molar basis
        //std::cout << "Quality flag: " << AS->get_quality_flag() << " (1=molar, 2=mass)" << std::endl;

        // Find data at 278.15 K
        double target_T = 278.15;
        AS->update(QT_INPUTS, 0.5, target_T);  // Saturated vapor (q=1 in molar basis)
        double p_sat = AS->p();
        std::cout << "  Saturation pressure (molar basis): " << std::setprecision(12) << AS->p() << " Pa" << std::endl;        
        std::cout << "  Vapor density: " << std::setprecision(12) << AS->rhomolar() << " kmol/m3" << std::endl;
        std::cout << "  Vapor density: " << std::setprecision(12) << AS->rhomass() << " kg/m3" << std::endl;
        std::cout << "  Quality (molar): " << std::setprecision(12) << AS->Q() << " -" << std::endl;
        std::cout << "  Quality (mass): " << std::setprecision(12) << AS->Qmass() << " -" << std::endl;
        auto mf = AS->get_mole_fractions(); std::cout << "  mole fractions: "; for (const auto& x : mf) std::cout << std::setprecision(12) << x << " "; std::cout << std::endl;
        auto mf_mass = AS->get_mass_fractions(); std::cout << "  mass fractions: "; for (const auto& x : mf_mass) std::cout << std::setprecision(12) << x << " "; std::cout << std::endl;

        //AS->set_quality_flag(2);  // Switch to mass basis
        AS->update(QmassT_INPUTS, 0.5, target_T);  // Saturated vapor (q=1 in mass basis)
        std::cout << "  Saturation pressure (mass basis): " << std::setprecision(12) << AS->p() << " Pa" << std::endl;        
        std::cout << "  Vapor density: " << std::setprecision(12) << AS->rhomolar() << " kmol/m3" << std::endl;
        std::cout << "  Vapor density: " << std::setprecision(12) << AS->rhomass() << " kg/m3" << std::endl;
        std::cout << "  Quality (molar): " << std::setprecision(12) << AS->Q() << " -" << std::endl;
        std::cout << "  Quality (mass): " << std::setprecision(12) << AS->Qmass() << " -" << std::endl;

        //AS->set_quality_flag(1);  // Switch to mole basis
        AS->update(PQ_INPUTS, 1e5, 0.5);  // Saturated vapor (q=1 in mass basis)
        std::cout << "  Saturation temperature (molar basis): " << std::setprecision(12) << AS->T() << " K" << std::endl;        
        std::cout << "  Vapor density: " << std::setprecision(12) << AS->rhomolar() << " kmol/m3" << std::endl;
        std::cout << "  Vapor density: " << std::setprecision(12) << AS->rhomass() << " kg/m3" << std::endl;
        std::cout << "  Quality (molar): " << std::setprecision(12) << AS->Q() << " -" << std::endl;
        std::cout << "  Quality (mass): " << std::setprecision(12) << AS->Qmass() << " -" << std::endl;

        //AS->set_quality_flag(2);  // Switch to mass basis
        AS->update(PQmass_INPUTS, 1e5, 0.5);  // Saturated vapor (q=1 in mass basis)
        std::cout << "  Saturation temperature (mass basis): " << std::setprecision(12) << AS->T() << " K" << std::endl;        
        std::cout << "  Vapor density: " << std::setprecision(12) << AS->rhomolar() << " kmol/m3" << std::endl;
        std::cout << "  Vapor density: " << std::setprecision(12) << AS->rhomass() << " kg/m3" << std::endl;
        std::cout << "  Quality (molar): " << std::setprecision(12) << AS->Q() << " -" << std::endl;
        std::cout << "  Quality (mass): " << std::setprecision(12) << AS->Qmass() << " -" << std::endl;



        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    }
}


