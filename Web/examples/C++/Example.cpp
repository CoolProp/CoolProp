#include "CoolProp.h"
#include "HumidAirProp.h"
#include "CPState.h"
#include <iostream>
#include <stdlib.h>

#pragma GCC diagnostic ignored "-Wwrite-strings" //Ignore char* warnings

int main() {
    double T, h, p, D;

    printf("CoolProp version:\t%s\n",     get_global_param_string("version").c_str());
    printf("CoolProp gitrevision:\t%s\n", get_global_param_string("gitrevision").c_str());
    printf("CoolProp fluids:\t%s\n",      get_global_param_string("FluidsList").c_str());

    printf("\n************ USING EOS *************\n");

    printf("FLUID STATE INDEPENDENT INPUTS\n");
    printf("Critical Density Propane: %f kg/m^3\n", Props1("Propane", "rhocrit"));

    printf("\nTWO PHASE INPUTS (Pressure)\n");
    printf("Density of saturated liquid Propane at 101.325 kPa: %f kg/m^3\n",   Props("D", 'P', 101.325, 'Q', 0, "Propane"));
    printf("Density of saturated vapor R290 at 101.325 kPa:     %f kg/m^3\n",   Props("D", 'P', 101.325, 'Q', 1, "R290"));

    printf("\nTWO PHASE INPUTS (Temperature)\n");
    printf("Density of saturated liquid Propane at 300 K: %f kg/m^3\n", Props("D", 'T', 300, 'Q', 0, "Propane"));
    printf("Density of saturated vapor R290 at 300 K:     %f kg/m^3\n", Props("D", 'T', 300, 'Q', 1, "R290"));

    printf("\nSINGLE PHASE CYCLE (Propane)\n");
    p = Props("P", 'T', 300, 'D', 1, "Propane");
    h = Props("H", 'T', 300, 'D', 1, "Propane");
    printf("T,D -> P,H : 300,1 -> %f,%f\n", p, h);

    T = Props("T", 'P', p, 'H', h, "Propane");
    D = Props("D", 'P', p, 'H', h, "Propane");
    printf("P,H -> T,D : %f, %f -> %f, %f\n", p, h, T, D);

    printf("\n************ USING TTSE ***************\n");
    enable_TTSE_LUT("Propane");

    printf("TWO PHASE INPUTS (Pressure)\n");
    printf("Density of saturated liquid Propane at 101.325 kPa: %f kg/m^3\n", Props("D", 'P', 101.325, 'Q', 0, "Propane"));
    printf("Density of saturated vapor R290 at 101.325 kPa:     %f kg/m^3\n", Props("D", 'P', 101.325, 'Q', 1, "R290"));

    printf("\nTWO PHASE INPUTS (Temperature)");
    printf("Density of saturated liquid Propane at 300 K: %f kg/m^3\n", Props("D", 'T', 300, 'Q', 0, "Propane"));
    printf("Density of saturated vapor R290 at 300 K:     %f kg/m^3\n", Props("D", 'T', 300, 'Q', 1, "R290"));

    printf("\nSINGLE PHASE CYCLE (propane)\n");
    p = Props("P", 'T', 300, 'D', 1, "Propane");
    h = Props("H", 'T', 300, 'D', 1, "Propane");
    printf("T,D -> P,H : 300,1 -> %f,%f", p, h);

    T = Props("T", 'P', p, 'H', h, "Propane");
    D = Props("D", 'P', p, 'H', h, "Propane");
    printf("P,H -> T,D : %f, %f -> %f, %f\n", p, h, T, D);

    disable_TTSE_LUT("Propane");

    try
    {
        printf("\n************ USING REFPROP ***************\n");
        std::string RPName = std::string("REFPROP-")+get_fluid_param_string("Propane","REFPROPname");
        printf("TWO PHASE INPUTS (Pressure)");
        printf("Density of saturated liquid Propane at 101.325 kPa: %f kg/m^3\n", Props("D", 'P', 101.325, 'Q', 0, RPName));
        printf("Density of saturated vapor R290 at 101.325 kPa:     %f kg/m^3\n", Props("D", 'P', 101.325, 'Q', 1, RPName));

        printf("\nTWO PHASE INPUTS (Temperature)");
        printf("Density of saturated liquid Propane at 300 K: %f kg/m^3\n", Props("D", 'T', 300, 'Q', 0, RPName));
        printf("Density of saturated vapor R290 at 300 K:     %f kg/m^3\n", Props("D", 'T', 300, 'Q', 1, RPName));

        printf("\nSINGLE PHASE CYCLE (propane)\n");
        p = Props("P", 'T', 300, 'D', 1, RPName);
        h = Props("H", 'T', 300, 'D', 1, RPName);
        printf("T,D -> P,H : 300,1 -> %f,%f\n", p, h);

        T = Props("T", 'P', p, 'H', h, RPName);
        D = Props("D", 'P', p, 'H', h, RPName);
        printf("P,H -> T,D : %f, %f -> %f, %f\n", p, h, T, D);
    }
    catch (std::exception &e)
    {
        printf("\n************ CANT USE REFPROP ************\n");
    }

    printf("\n************ CHANGE UNIT SYSTEM (default is kSI) *************\n");
    set_standard_unit_system(UNIT_SYSTEM_SI);
    printf("Vapor pressure of water at 373.15 K in SI units (Pa):   %f\n", Props("P", 'T', 373.15, 'Q', 0, "Water"));

    set_standard_unit_system(UNIT_SYSTEM_KSI);
    printf("Vapor pressure of water at 373.15 K in kSI units (kPa): %f\n", Props("P", 'T', 373.15, 'Q', 0, "Water"));

    printf("\n************ BRINES AND SECONDARY WORKING FLUIDS *************\n");
    printf("Density of 50%% (mass) ethylene glycol/water at 300 K, 101.325 kPa: %f kg/m^3\n", Props("D", 'T', 300, 'P', 101.325, "EG-50%"));
    printf("Viscosity of Therminol D12 at 350 K, 101.325 kPa: %f Pa-s\n",  Props("V", 'T', 350, 'P', 101.325, "TD12"));

    printf("\n************ HUMID AIR PROPERTIES *************\n");
    printf("Humidity ratio of 50%% rel. hum. air at 300 K, 101.325 kPa: %f kg_w/kg_da\n", HAProps("W", "T", 300, "P", 101.325, "R", 0.5));
    printf("Relative humidity from last calculation: %f (fractional)\n",   HAProps("R", "T", 300, "P", 101.325, "W", HAProps("W", "T", 300, "P", 101.325, "R", 0.5)));
    return 0;
}
