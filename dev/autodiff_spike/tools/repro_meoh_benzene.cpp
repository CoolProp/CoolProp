#include "CoolProp/AbstractState.h"
#include <cstdio>
#include <memory>
int main() {
    for (double x : {0.50, 0.52, 0.54, 0.56, 0.58, 0.60, 0.70, 0.76, 0.78, 0.80}) {
        auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "methanol&benzene"));
        AS->set_mole_fractions({x, 1 - x});
        try {
            AS->update(CoolProp::PT_INPUTS, 101325, 308.15);
            std::printf("x=%.2f phase=%d Q=%.6g rho=%.8g\n", x, (int)AS->phase(), AS->Q(), AS->rhomolar());
        } catch (std::exception& e) {
            std::printf("x=%.2f threw %s\n", x, e.what());
        }
    }
}
