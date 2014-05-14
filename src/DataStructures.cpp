

#include "DataStructures.h"
#include "Exceptions.h"
#include "CoolPropTools.h"

namespace CoolProp{

struct parameter_info
{
    int key;
    std::string short_desc, IO, units, description;
public:
    parameter_info(int key, std::string short_desc, std::string IO, std::string units, std::string description): key(key), short_desc(short_desc), IO(IO), units(units), description(description){};
};

parameter_info parameter_info_list[] = {
    /// Input/Output parameters
    parameter_info(iT, "T", "IO", "K", "Temperature"),
    parameter_info(iP, "P", "IO", "Pa", "Pressure"),
    parameter_info(iDmolar, "Dmolar","IO","mol/m^3","Molar density"),
    parameter_info(iHmolar, "Hmolar","IO","J/mol","Molar specific enthalpy"),
    parameter_info(iSmolar, "Smolar","IO","J/mol/K","Molar specific entropy"),
    parameter_info(iUmolar, "Umolar","IO","J/mol","Molar specific internal energy"),
    parameter_info(iDmass, "Dmass","IO","kg/m^3","Mass density"),
    parameter_info(iHmass, "Hmass","IO","J/kg","Mass specific enthalpy"),
    parameter_info(iSmass, "Smass","IO","J/kg/K","Mass specific entropy"),
    parameter_info(iUmass, "Umass","IO","J/kg","Mass specific internal energy"),
    parameter_info(iQ, "Q","IO","mol/mol","Mass vapor quality"),
    parameter_info(iDelta, "Delta","IO","-","Reduced density (rho/rhoc)"),
    parameter_info(iTau, "Tau","IO","-","Reciprocal reduced temperature (Tc/T)"),
    /// Output only
    parameter_info(iCpmolar, "Cpmolar","O","J/mol/K","Molar specific constant presssure specific heat"),
    parameter_info(iCpmass, "Cpmass","O","J/kg/K","Mass specific constant presssure specific heat"),
    parameter_info(iCvmolar, "Cvmolar","O","J/mol/K","Molar specific constant volume specific heat"),
    parameter_info(iCvmass, "Cvmass","O","J/kg/K","Mass specific constant volume specific heat"),
    parameter_info(iGWP20, "GWP20","O","-","20-year gobal warming potential"),
    parameter_info(iGWP100, "GWP100","O","-","100-year gobal warming potential"),
    parameter_info(iGWP500, "GWP500","O","-","500-year gobal warming potential"),
    parameter_info(iFH, "FH","O","-","Flammability hazard"),
    parameter_info(iHH, "HH","O","-","Health hazard"),
    parameter_info(iPH, "PH","O","-","Physical hazard"),
    parameter_info(iODP, "ODP","O","-","Ozone depletion potential"),
    parameter_info(iBvirial, "Bvirial","O","-","Second virial coefficient"),
    parameter_info(iCvirial, "Cvirial","O","-","Third virial coefficient"),
    parameter_info(idBvirial_dT, "dBvirial_dT","O","-","Derivative of second virial coefficient with respect to T"),
    parameter_info(idCvirial_dT, "dCvirial_dT","O","-","Derivative of third virial coefficient with respect to T"),
    parameter_info(imolar_mass, "molar_mass","O","kg/mol","Molar mass"),
    parameter_info(irhomolar_reducing, "rhomolar_reducing","O","mol/m^3","Molar density at reducing point"),
    parameter_info(irhomolar_critical, "rhomolar_critical","O","mol/m^3","Molar density at critical point"),
    parameter_info(iT_reducing, "T_reducing","O","K","Temperature at the reducing point"),
    parameter_info(iT_critical, "T_critical","O","K","Temperature at the critical point"),
    parameter_info(iisothermal_compressibility, "isothermal_compressibility","O","1/Pa","Isothermal compressibility")
};

class ParameterInformation
{
public:
    std::map<int, std::string> short_desc_map, description_map, IO_map, units_map;
    std::map<std::string, int> index_map;
    ParameterInformation()
    {
        int N = sizeof(parameter_info_list)/sizeof(parameter_info_list[0]);
        for (int i = 0; i < N; ++i)
        {
            parameter_info &el = parameter_info_list[i];
            short_desc_map.insert(std::pair<int, std::string>(el.key, el.short_desc));
            IO_map.insert(std::pair<int, std::string>(el.key, el.IO));
            units_map.insert(std::pair<int, std::string>(el.key, el.units));
            description_map.insert(std::pair<int, std::string>(el.key, el.description));
            index_map.insert(std::pair<std::string, int>(el.short_desc, el.key));
        }
        // Backward compatibility aliases
        index_map.insert(std::pair<std::string, int>("D", iDmass));
        index_map.insert(std::pair<std::string, int>("H", iHmass));
        index_map.insert(std::pair<std::string, int>("S", iSmass));
        index_map.insert(std::pair<std::string, int>("U", iUmass));
        index_map.insert(std::pair<std::string, int>("C", iCpmass));
        index_map.insert(std::pair<std::string, int>("O", iCvmass));
    }
};

static ParameterInformation parameter_information;

std::string get_parameter_information(int key, std::string info)
{
    std::map<int, std::string> *M;
    std::map<int, std::string>::iterator it;

    // Hook up the right map (since they are all of the same type)
    if (!info.compare("IO")){
        M = &(parameter_information.IO_map);
    }
    else if (!info.compare("short")){
        M = &(parameter_information.short_desc_map);
    }
    else if (!info.compare("long")){
        M = &(parameter_information.description_map);
    }
    else if (!info.compare("units")){
        M = &(parameter_information.units_map);
    }
    else
        throw ValueError(format("Bad info string [%s] to get_parameter_information",info.c_str()));

    // Try to find it
    it = (*M).find(key);
    // If equal to end, not found
    if (it != (*M).end())
    {
        // Found it, return it
        return it->second;
    }
    else
    {
        throw ValueError(format("Unable to match the key [%d] in get_parameter_information for info [%s]",key, info.c_str()));
    }
}

int get_parameter_index(const std::string &param_name)
{
    std::map<std::string, int>::iterator it;
    // Try to find it
    it = parameter_information.index_map.find(param_name);
    // If equal to end, not found
    if (it != parameter_information.index_map.end())
    {
        // Found, return it
        return it->second;
    }
    else
    {
        throw ValueError(format("Your input name [%s] is not valid in get_parameter_index (names are case sensitive)",param_name.c_str()));
    }
}

struct input_pair_info
{
    int key;
    std::string short_desc, long_desc;
public:
    input_pair_info(int key, std::string short_desc, std::string long_desc): key(key), short_desc(short_desc), long_desc(long_desc){};
};

input_pair_info input_pair_list[] = {
    input_pair_info(QT_INPUTS,"QT_INPUTS","Molar quality, Temperature in K"),
    input_pair_info(PQ_INPUTS,"PQ_INPUTS","Pressure in Pa, Molar quality"),
    input_pair_info(PT_INPUTS, "PT_INPUTS","Pressure in Pa, Temperature in K"),

    input_pair_info(DmassT_INPUTS, "DmassT_INPUTS", "Mass density in kg/m^3, Temperature in K"),
    input_pair_info(DmolarT_INPUTS, "DmolarT_INPUTS", "Molar density in mol/m^3, Temperature in K"),
    input_pair_info(HmassT_INPUTS, "HmassT_INPUTS", "Enthalpy in J/kg, Temperature in K"),
    input_pair_info(HmolarT_INPUTS, "HmolarT_INPUTS", "Enthalpy in J/mol, Temperature in K"),
    input_pair_info(SmassT_INPUTS, "SmassT_INPUTS", "Entropy in J/kg/K, Temperature in K"),
    input_pair_info(SmolarT_INPUTS, "SmolarT_INPUTS", "Entropy in J/mol/K, Temperature in K"),
    input_pair_info(TUmass_INPUTS, "TUmass_INPUTS", "Temperature in K, Internal energy in J/kg"),
    input_pair_info(TUmolar_INPUTS, "TUmolar_INPUTS", "Temperature in K, Internal energy in J/mol"),
                 
    input_pair_info(DmassP_INPUTS, "DmassP_INPUTS", "Mass density in kg/m^3, Pressure in Pa"),
    input_pair_info(DmolarP_INPUTS, "DmolarP_INPUTS", "Molar density in mol/m^3, Pressure in Pa"),
    input_pair_info(HmassP_INPUTS, "HmassP_INPUTS", "Enthalpy in J/kg, Pressure in Pa"),
    input_pair_info(HmolarP_INPUTS, "HmolarP_INPUTS", "Enthalpy in J/mol, Pressure in Pa"),
    input_pair_info(PSmass_INPUTS, "PSmass_INPUTS", "Pressure in Pa, Entropy in J/kg/K"),
    input_pair_info(PSmolar_INPUTS, "PSmolar_INPUTS", "Pressure in Pa, Entropy in J/mol/K "),
    input_pair_info(PUmass_INPUTS, "PUmass_INPUTS", "Pressure in Pa, Internal energy in J/kg"),
    input_pair_info(PUmolar_INPUTS, "PUmolar_INPUTS", "Pressure in Pa, Internal energy in J/mol"),

    input_pair_info(DmassHmass_INPUTS, "DmassHmass_INPUTS","Mass density in kg/m^3, Enthalpy in J/kg"),
    input_pair_info(DmolarHmolar_INPUTS, "DmolarHmolar_INPUTS","Molar density in mol/m^3, Enthalpy in J/mol"),
    input_pair_info(DmassSmass_INPUTS, "DmassSmass_INPUTS","Mass density in kg/m^3, Entropy in J/kg/K"),
    input_pair_info(DmolarSmolar_INPUTS, "DmolarSmolar_INPUTS","Molar density in mol/m^3, Entropy in J/mol/K"),
    input_pair_info(DmassUmass_INPUTS, "DmassUmass_INPUTS","Mass density in kg/m^3, Internal energy in J/kg"),
    input_pair_info(DmolarUmolar_INPUTS, "DmolarUmolar_INPUTS","Molar density in mol/m^3, Internal energy in J/mol"),

    input_pair_info(HmassSmass_INPUTS, "HmassSmass_INPUTS", "Enthalpy in J/kg, Entropy in J/kg/K"),
    input_pair_info(HmolarSmolar_INPUTS, "HmolarSmolar_INPUTS", "Enthalpy in J/mol, Entropy in J/mol/K"),
    input_pair_info(SmassUmass_INPUTS, "SmassUmass_INPUTS", "Entropy in J/kg/K, Internal energy in J/kg"),
    input_pair_info(SmolarUmolar_INPUTS, "SmolarUmolar_INPUTS", "Entropy in J/mol/K, Internal energy in J/mol"),
};

class InputPairInformation
{
public:
    std::map<int, std::string> short_desc_map, long_desc_map;
    InputPairInformation()
    {
        int N = sizeof(input_pair_list)/sizeof(input_pair_list[0]);
        for (int i = 0; i < N; ++i)
        {
            short_desc_map.insert(std::pair<int, std::string>(input_pair_list[i].key, input_pair_list[i].short_desc));
            long_desc_map.insert(std::pair<int, std::string>(input_pair_list[i].key, input_pair_list[i].long_desc));
        }
    }
};

static InputPairInformation input_pair_information;

std::string get_input_pair_short_desc(int pair)
{
    return input_pair_information.short_desc_map[pair];
}
std::string get_input_pair_long_desc(int pair)
{
    return input_pair_information.long_desc_map[pair];
}

} /* namespace CoolProp */