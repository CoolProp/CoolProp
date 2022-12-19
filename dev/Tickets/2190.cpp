
#include <vector>
#include <string>
#include <chrono>
#include <memory>

#include "CoolProp.h"
#include "HumidAirProp.h"

int main(int argc, const char* argv[]) {

    CoolProp::set_debug_level(1000);

    double T = 393.15;
    double p = 101325;
    double R = 0.1;

    double h = HumidAir::HAPropsSI("H", "T", T, "P", p, "R", R);
    double R_test = HumidAir::HAPropsSI("R", "T", T, "P", p, "H", h);

    return 0;
}
