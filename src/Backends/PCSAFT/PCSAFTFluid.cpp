#include <string>
#include <vector>
#include <cmath>

#include "CoolProp/fluids/PCSAFTFluid.h"

namespace CoolProp {

void PCSAFTFluid::calc_water_sigma(double t) {
    if (t > 473.16) {
        throw ValueError("The current function for sigma for water is only valid for temperatures below 473.15 K.");
    } else if (t < 273) {
        throw ValueError("The current function for sigma for water is only valid for temperatures above 273.15 K.");
    }

    params.sigma = 3.8395 + 1.2828 * exp(-0.0074944 * t) - 1.3939 * exp(-0.00056029 * t);
}

} /* namespace CoolProp */
