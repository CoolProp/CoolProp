#include <CoolProp/CoolPropLib.h>

#include <cmath>
#include <cstdio>

int main() {
    const double temperature = PropsSI("T", "P", 101325.0, "Q", 0.0, "Water");

    if (!std::isfinite(temperature) || temperature < 350.0 || temperature > 400.0) {
        std::fprintf(stderr, "CoolProp static C wrapper returned an invalid boiling temperature: %.17g K\n", temperature);
        return 1;
    }

    std::printf("CoolProp static C wrapper water boiling temperature: %.12g K\n", temperature);
    return 0;
}
