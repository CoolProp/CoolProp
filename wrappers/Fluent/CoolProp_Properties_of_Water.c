#include "udf.h"
#include <string.h>
struct string;
const char FLUID[] = "Water";
const real gauge = 101325; /*operating pressure in pascal (as defined in fluent) */
double PropsSI(char*, char*, double, char*, double, char*);

/*Density of the FLUID[]*/
DEFINE_PROPERTY(water_density, c, t) {
    real temperature = C_T(c, t);
    real pressure = C_P(c, t) + gauge;
    real density;
    density = PropsSI((char*)"D", (char*)"T", temperature, (char*)"P", pressure, (char*)FLUID);
    return density;
}
/*Thermal Conductivity of the FLUID[]*/
DEFINE_PROPERTY(water_thermalConductivity, c, t) {
    real thermalConductivity;
    real temperature = C_T(c, t);
    real pressure = C_P(c, t) + gauge;
    thermalConductivity = PropsSI((char*)"L", (char*)"T", temperature, (char*)"P", pressure, (char*)FLUID);
    return thermalConductivity;
}
/*Viscosity of the FLUID[]*/
DEFINE_PROPERTY(water_viscosity, c, t) {
    real viscosity;
    real temperature = C_T(c, t);
    real pressure = C_P(c, t) + gauge;
    viscosity = PropsSI((char*)"V", (char*)"T", temperature, (char*)"P", pressure, (char*)FLUID);
    return viscosity;
}
/*Specific heat of the FLUID[]*/
DEFINE_SPECIFIC_HEAT(water_specificHeat, temperature, Tref, enthalpy, yi) {
    real pressure;
    /*
        density = 1.7730;
	*/
    /* The following commented code is supposed to get the pressure
	 from the cell to use with Coolprop. Left commented because
	 specific heat depends very little on pressure. Will increase
	 computational time significantly. */

    Domain* domain = Get_Domain(1);
    Thread* t;
    cell_t c;
    thread_loop_c(t, domain) {
        begin_c_loop(c, t) {
            real pressure = C_P(c, t) + gauge;
            temperature = C_T(c, t);
            real specificHeat;
            specificHeat = PropsSI((char*)"C", (char*)"T", temperature, (char*)"P", pressure, (char*)FLUID);
            *enthalpy = specificHeat * (temperature - Tref);
            return specificHeat;
        }
        end_c_loop(c, t)
    }
}
/* test coolprop integration */
DEFINE_ON_DEMAND(call_coolprop_water) {
    real p, t, density, specificHeat, thermalConductivity, enthalpy, viscosity;

    p = 101325.0;
    t = 298.15;

    density = PropsSI((char*)"D", (char*)"T", t, (char*)"P", p, (char*)FLUID);
    specificHeat = PropsSI((char*)"C", (char*)"T", t, (char*)"P", p, (char*)FLUID);
    viscosity = PropsSI((char*)"V", (char*)"T", t, (char*)"P", p, (char*)FLUID);
    thermalConductivity = PropsSI((char*)"L", (char*)"T", t, (char*)"P", p, (char*)FLUID);
    enthalpy = PropsSI((char*)"H", (char*)"T", t, (char*)"P", p, (char*)FLUID);

    Message("p = %lf, T = %lf => density = %lf\n", p, t, density);
    Message("p = %lf, T = %lf => specific heat = %lf\n", p, t, specificHeat);
    Message("p = %lf, T = %lf => viscosity = %lf\n", p, t, viscosity);
    Message("p = %lf, T = %lf => thermal conductivity = %lf\n", p, t, thermalConductivity);
    Message("p = %lf, T = %lf => enthalpy = %lf\n", p, t, enthalpy);
}
