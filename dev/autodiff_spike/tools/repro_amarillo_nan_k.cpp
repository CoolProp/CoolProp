#include "CoolProp/AbstractState.h"
#include <cstdio>
#include <memory>
int main() {
    const std::string fl = "Methane&Nitrogen&CarbonDioxide&Ethane&Propane&IsoButane&n-Butane&Isopentane&n-Pentane&n-Hexane";
    const std::vector<double> z = {0.906724, 0.031284, 0.004676, 0.045279, 0.00828, 0.001037, 0.001563, 0.000321, 0.000443, 0.000393};
    for (auto s : {std::pair<double, double>{176.9762835, 15145363.4}, {180.1382565, 8593226.381}}) {
        std::unique_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("GERG2008", fl));
        AS->set_mole_fractions(z);
        try {
            AS->update(CoolProp::PT_INPUTS, s.second, s.first);
            std::printf("T=%.4f p=%.6g  rho=%.6g Q=%g phase=%d\n", s.first, s.second, AS->rhomolar(), AS->Q(), (int)AS->phase());
        } catch (std::exception& e) {
            std::printf("T=%.4f p=%.6g  THREW: %s\n", s.first, s.second, e.what());
        }
    }
}
