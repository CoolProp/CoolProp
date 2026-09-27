#include "CoolProp/AbstractState.h"
#include <cstdio>
#include <memory>
int main() {
    const std::string fluids = "Methane&Ethane&n-Propane&n-Butane&IsoButane&Nitrogen&CarbonDioxide";
    const std::vector<double> z = {0.9188, 0.0532, 0.0193, 0.0010, 0.0012, 0.0064, 0.0001};
    for (double T : {223.10, 223.15, 223.20, 223.28, 223.35, 223.40, 207.91, 208.03, 208.66, 211.77}) {
        auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluids));
        AS->set_mole_fractions(z);
        try { AS->update(CoolProp::PT_INPUTS, 60e5, T); std::printf("T=%.2f phase=%d Q=%.6g rho=%.6g\n", T, (int)AS->phase(), AS->Q(), AS->rhomolar()); }
        catch (std::exception& e) { std::printf("T=%.2f threw %s\n", T, e.what()); }
    }
}
