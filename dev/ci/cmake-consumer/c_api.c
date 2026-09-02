#include <CoolProp/CoolPropLib.h>

#include <math.h>
#include <stdio.h>

int main(void) {
    const double temperature = PropsSI("T", "P", 101325.0, "Q", 0.0, "Water");

    if (!isfinite(temperature) || temperature < 350.0 || temperature > 400.0) {
        fprintf(stderr, "CoolProp C API returned an invalid boiling temperature: %.17g K\n", temperature);
        return 1;
    }

    printf("CoolProp C API water boiling temperature: %.12g K\n", temperature);
    return 0;
}
