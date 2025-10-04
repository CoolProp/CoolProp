#include "udf.h"
#include "CoolPropLib.h"
#include <string>

#pragma comment(lib, "C:\\XX\\XX\\CoolProp.lib") 

// User settings
static const std::string FLUID_NAME = "R134a";

DEFINE_ON_DEMAND(test)
{
    double PropResult = PropsSI("T", "P", 101325, "Q", 0, "Water");
    Message0("Saturation temperature of Water at 1 atm = %g K\n", PropResult);    
}

// Density via CoolProp
DEFINE_PROPERTY(rho_coolprop, cell, thread)
{
    real P = C_P(cell, thread);
    real H = C_H(cell, thread);
    double rho = PropsSI("VISCOSITY",
                         "P",    static_cast<double>(P),
                         "Hmass",static_cast<double>(H),
                         FLUID_NAME.c_str());
    return static_cast<real>(rho);
}

// Viscosity via CoolProp
DEFINE_PROPERTY(mu_coolprop, cell, thread)
{
    real P = C_P(cell, thread);
    real H = C_H(cell, thread);
    double mu = PropsSI("VISCOSITY",
                        "P",    static_cast<double>(P),
                        "Hmass",static_cast<double>(H),
                        FLUID_NAME.c_str());
    return static_cast<real>(mu);
}

// Thermal conductivity via CoolProp
DEFINE_PROPERTY(k_coolprop, cell, thread)
{
    real P = C_P(cell, thread);
    real H = C_H(cell, thread);
    double k = PropsSI("CONDUCTIVITY",
                       "P",    static_cast<double>(P),
                       "Hmass",static_cast<double>(H),
                       FLUID_NAME.c_str());
    return static_cast<real>(k);
}
