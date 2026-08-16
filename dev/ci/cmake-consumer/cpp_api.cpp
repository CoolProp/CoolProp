#include <CoolProp/CoolProp.h>
// Intentionally unused: verifies nested installed headers and transitive Eigen includes.
#include <CoolProp/numerics/PolyMath.h>

#include <cmath>
#include <iostream>
#include <string>
#include <string_view>

int main() {
    constexpr std::string_view output = "T";
    const double temperature = CoolProp::PropsSI(std::string(output), "P", 101325.0, "Q", 0.0, "Water");

    if (!std::isfinite(temperature) || temperature < 350.0 || temperature > 400.0) {
        std::cerr << "CoolProp C++ API returned an invalid boiling temperature: " << temperature << " K\n";
        return 1;
    }

    std::cout << "CoolProp C++ API water boiling temperature: " << temperature << " K\n";
    return 0;
}
