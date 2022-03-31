/*
* UDF TO CALCULATE FLUID PROPERTIES BASED ON THE OPEN-SOURCE
* THERMODYNAMIC LIBRARY COOLPROP
*/

/* PropsSI() is faster if you use temperature and density as
 * independent state variables */

#include "udf.h"

const char* FLUID = "CarbonDioxide";
const real gauge = 101325;  //operating pressure in pascal (as defined in fluent)

double PropsSI(char*, char*, double, char*, double, char*);

/* REAL GAS VISCOSITY */
DEFINE_PROPERTY(cell_viscosity, c, t) {
    real viscosity;
    real temperature = C_T(c, t);
    real density = C_R(c, t);
    viscosity = PropsSI((char*)"V", (char*)"T", temperature, (char*)"D", density, (char*)FLUID);
    return viscosity;
}

/* REAL GAS THERMAL CONDUCTIVITY */
DEFINE_PROPERTY(cell_thermalConductivity, c, t) {
    real thermalConductivity;
    real temperature = C_T(c, t);
    real density = C_R(c, t);
    thermalConductivity = PropsSI((char*)"L", (char*)"T", temperature, (char*)"D", density, (char*)FLUID);
    return thermalConductivity;
}

/* REAL GAS SPECIFIC MASS */
DEFINE_PROPERTY(cell_density, c, t) {
    real density;
    real temperature = C_T(c, t);
    real pressure = C_P(c, t) + gauge;
    density = PropsSI((char*)"D", (char*)"T", temperature, (char*)"P", pressure, (char*)FLUID);
    return density;
}

/* REAL GAS SPECIFIC HEAT */
DEFINE_SPECIFIC_HEAT(cell_specificHeat, temperature, Tref, enthalpy, yi) {
    real density;
    density = 1.7730;

    /* The following commented code is supposed to get the pressure
	 from the cell to use with Coolprop. Left commented because
	 specific heat depends very little on pressure. Will increase
	 computational time significantly. */

    /*
	* Domain *domain = Get_Domain(1);
	* Thread *t;
	* cell_t c;
	* thread_loop_c(t, domain)
	* {
	* 	begin_c_loop(c, t)
	* 	{
	* 		pressure = C_P(c, t);
	* 	}end_c_loop(c, t)
	* }
	*/

    real specificHeat;

    specificHeat = PropsSI((char*)"C", (char*)"T", temperature, (char*)"D", density, (char*)FLUID);
    *enthalpy = specificHeat * (temperature - Tref);
    return specificHeat;
}

/* Execute on demand UDF to test if the library was built correctly */
DEFINE_ON_DEMAND(call_coolprop) {
    real p, t, density, specificHeat, thermalConductivity, enthalpy, viscosity;

    p = 100000.0;
    t = 300.0;

    density = PropsSI((char*)"D", (char*)"T", t, (char*)"P", p, (char*)FLUID);
    specificHeat = PropsSI((char*)"C", (char*)"T", t, (char*)"D", density, (char*)FLUID);
    viscosity = PropsSI((char*)"V", (char*)"T", t, (char*)"D", density, (char*)FLUID);
    thermalConductivity = PropsSI((char*)"L", (char*)"T", t, (char*)"D", density, (char*)FLUID);
    enthalpy = PropsSI((char*)"H", (char*)"T", t, (char*)"D", density, (char*)FLUID);

    Message("p = %lf, T = %lf => density = %lf\n", p, t, density);
    Message("p = %lf, T = %lf => specific heat = %lf\n", p, t, specificHeat);
    Message("p = %lf, T = %lf => viscosity = %lf\n", p, t, viscosity);
    Message("p = %lf, T = %lf => thermal conductivity = %lf\n", p, t, thermalConductivity);
    Message("p = %lf, T = %lf => enthalpy = %lf\n", p, t, enthalpy);
}
