#include "CoolProp.h"
#include "HumidAirProp.h"
#include "AbstractState.h"
#include "crossplatform_shared_ptr.h"
#include <iostream>
#include <stdlib.h>

using namespace CoolProp; // For normal CoolProp calls
using namespace HumidAir; // For HAPropsSI
int main() {
    double T, h, p, D;

    printf("\n************ Information *************\n");
    printf("CoolProp version:\t%s\n",     get_global_param_string("version").c_str());
    printf("CoolProp gitrevision:\t%s\n", get_global_param_string("gitrevision").c_str());
    printf("CoolProp fluids:\t%s\n",      get_global_param_string("FluidsList").c_str());

    printf("\n************ USING HIGH LEVEL INTERFACE *************\n");
    printf("Critical temperature of water: %f K\n", Props1SI("Water", "Tcrit"));
    printf("Boiling temperature of water at 101325 Pa: %f K\n",   PropsSI("T", "P", 101325, "Q", 0, "Water"));
    printf("Phase of water at 101325 Pa and 300 K: %s\n",  PhaseSI("P", 101325, "Q", 0, "Water").c_str());
    printf("cp of water at 101325 Pa and 300 K: %g J/kg/K\n",  PropsSI("C", "P", 101325, "T", 300, "Water"));
    printf("cp of water (using derivatives) at 101325 Pa and 300 K: %g J/kg/K\n", PropsSI("d(H)/d(T)|P", "P", 101325, "T", 300, "Water"));

    printf("\n************ HUMID AIR PROPERTIES *************\n");
    printf("Humidity ratio of 50%% rel. hum. air at 300 K, 101325 Pa: %f kg_w/kg_da\n", HAPropsSI("W", "T", 300, "P", 101325, "R", 0.5));
    printf("Relative humidity from last calculation: %f (fractional)\n",  HAPropsSI("R", "T", 300, "P", 101325, "W", HAPropsSI("W", "T", 300, "P", 101325, "R", 0.5)));
    
    printf("\n************ BRINES AND SECONDARY WORKING FLUIDS *************\n");
    printf("Density of 50%% (mass) ethylene glycol/water at 300 K, 101325 Pa: %f kg/m^3\n", PropsSI("D", "T", 300, "P", 101325, "INCOMP::MEG-50%"));
    printf("Viscosity of Therminol D12 at 350 K, 101325 Pa: %f Pa-s\n",  PropsSI("V", "T", 350, "P", 101325, "INCOMP::TD12"));

	printf("\n************ LOW LEVEL INTERFACE *************\n");
	shared_ptr<AbstractState> Water(AbstractState::factory("HEOS", "Water"));
	printf("Critical temperature of water: %f K\n", Water->T_critical());
	Water->update(PQ_INPUTS, 101325, 1);
    printf("Boiling temperature of water at 101325 Pa: %f K\n", Water->T());
    Water->update(PT_INPUTS, 101325, 300);
	printf("Phase index of water at 101325 Pa and 300 K: %d\n", Water->phase());
    printf("cp of water at 101325 Pa and 300 K: %g J/kg/K\n", Water->cpmass());
	printf("cp of water (using derivatives) at 101325 Pa and 300 K: %g J/kg/K\n", Water->first_partial_deriv(iHmass, iT, iP));

    return EXIT_SUCCESS;

}
