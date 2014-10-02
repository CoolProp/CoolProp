#ifndef COOLPROP_CONFIGURATION
#define COOLPROP_CONFIGURATION

#include "Exceptions.h"
#include "CoolPropTools.h"

enum configuration_keys {NORMALIZE_GAS_CONSTANTS, CRITICAL_SPLINES_ENABLED};

namespace CoolProp
{
//
//class Configuration
//{
//public:
//    Configuration();
//    ~Configuration();
//    
//};

/// Return the value of a boolean key from the configuration
bool get_config_bool(configuration_keys key);
double get_config_double(configuration_keys key);
std::string get_config_string(configuration_keys key);

}

#endif // COOLPROP_CONFIGURATION
