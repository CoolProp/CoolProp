

#include "DataStructures.h"
#include "Exceptions.h"
#include "CoolPropTools.h"
#include "CoolProp.h"

namespace CoolProp{

struct parameter_info
{    
    int key;
    std::string short_desc, IO, units, description;
    bool trivial; ///< True if the input is trivial, and can be directly calculated (constants like critical properties, etc.)
public:
    parameter_info(int key, std::string short_desc, std::string IO, std::string units, std::string description, bool trivial): key(key), short_desc(short_desc), IO(IO), units(units), description(description), trivial(trivial){};
};

parameter_info parameter_info_list[] = {
    /// Input/Output parameters
    parameter_info(iT, "T", "IO", "K", "Temperature",false),
    parameter_info(iP, "P", "IO", "Pa", "Pressure",false),
    parameter_info(iDmolar, "Dmolar","IO","mol/m^3","Molar density",false),
    parameter_info(iHmolar, "Hmolar","IO","J/mol","Molar specific enthalpy",false),
    parameter_info(iSmolar, "Smolar","IO","J/mol/K","Molar specific entropy",false),
    parameter_info(iUmolar, "Umolar","IO","J/mol","Molar specific internal energy",false),
    parameter_info(iGmolar, "Gmolar","O","J/mol","Molar specific Gibbs energy",false),
    parameter_info(iDmass, "Dmass","IO","kg/m^3","Mass density",false),
    parameter_info(iHmass, "Hmass","IO","J/kg","Mass specific enthalpy",false),
    parameter_info(iSmass, "Smass","IO","J/kg/K","Mass specific entropy",false),
    parameter_info(iUmass, "Umass","IO","J/kg","Mass specific internal energy",false),
    parameter_info(iGmass, "Gmass","O","J/kg","Mass specific Gibbs energy",false),
    parameter_info(iQ, "Q","IO","mol/mol","Mass vapor quality",false),
    parameter_info(iDelta, "Delta","IO","-","Reduced density (rho/rhoc)",false),
    parameter_info(iTau, "Tau","IO","-","Reciprocal reduced temperature (Tc/T)",false),
    /// Output only
    parameter_info(iCpmolar, "Cpmolar","O","J/mol/K","Molar specific constant presssure specific heat",false),
    parameter_info(iCpmass, "Cpmass","O","J/kg/K","Mass specific constant presssure specific heat",false),
    parameter_info(iCvmolar, "Cvmolar","O","J/mol/K","Molar specific constant volume specific heat",false),
    parameter_info(iCvmass, "Cvmass","O","J/kg/K","Mass specific constant volume specific heat",false),
    parameter_info(iCp0molar, "Cp0molar","O","J/mol/K","Ideal gas molar specific constant presssure specific heat",false),
    parameter_info(iCp0mass, "Cp0mass","O","J/kg/K","Ideal gas mass specific constant presssure specific heat",false),
    parameter_info(iGWP20, "GWP20","O","-","20-year gobal warming potential",true),
    parameter_info(iGWP100, "GWP100","O","-","100-year gobal warming potential",true),
    parameter_info(iGWP500, "GWP500","O","-","500-year gobal warming potential",true),
    parameter_info(iFH, "FH","O","-","Flammability hazard",true),
    parameter_info(iHH, "HH","O","-","Health hazard",true),
    parameter_info(iPH, "PH","O","-","Physical hazard",true),
    parameter_info(iODP, "ODP","O","-","Ozone depletion potential",true),
    parameter_info(iBvirial, "Bvirial","O","-","Second virial coefficient",false),
    parameter_info(iCvirial, "Cvirial","O","-","Third virial coefficient",false),
    parameter_info(idBvirial_dT, "dBvirial_dT","O","-","Derivative of second virial coefficient with respect to T",false),
    parameter_info(idCvirial_dT, "dCvirial_dT","O","-","Derivative of third virial coefficient with respect to T",false),
    parameter_info(igas_constant, "gas_constant","O","J/mol/K","Molar gas constant",true),
	parameter_info(imolar_mass, "molar_mass","O","kg/mol","Molar mass",true),
    parameter_info(irhomass_reducing, "rhomass_reducing","O","kg/m^3","Mass density at reducing point",true),
    parameter_info(irhomolar_reducing, "rhomolar_reducing","O","mol/m^3","Molar density at reducing point",true),
    parameter_info(irhomolar_critical, "rhomolar_critical","O","mol/m^3","Molar density at critical point",true),
    parameter_info(irhomass_critical, "rhomass_critical","O","kg/m^3","Mass density at critical point",true),
    parameter_info(iT_reducing, "T_reducing","O","K","Temperature at the reducing point",true),
    parameter_info(iT_critical, "T_critical","O","K","Temperature at the critical point",true),
    parameter_info(iT_triple, "T_triple","O","K","Temperature at the triple point",true),
    parameter_info(iT_max, "T_max","O","K","Maximum temperature limit",true),
    parameter_info(iT_min, "T_min","O","K","Minimum temperature limit",true),
    parameter_info(iP_min, "P_min","O","Pa","Minimum pressure limit",true),
    parameter_info(iP_max, "P_max","O","Pa","Maximum pressure limit",true),
    parameter_info(iP_critical, "p_critical","O","Pa","Pressure at the critical point",true),
	parameter_info(iP_reducing, "p_reducing","O","Pa","Pressure at the reducing point",true),
    parameter_info(iP_triple, "p_triple","O","Pa","Pressure at the triple point (pure only)",true),
    parameter_info(ifraction_min, "fraction_min","O","-","Fraction (mole, mass, volume) minimum value for incompressible solutions",true),
    parameter_info(ifraction_max, "fraction_max","O","-","Fraction (mole, mass, volume) maximum value for incompressible solutions",true),
    parameter_info(iT_freeze, "T_freeze","O","-","Freezing temperature incompressible solutions",true),
    
    parameter_info(ispeed_sound, "speed_of_sound","O","m/s","Speed of sound",false),
    parameter_info(iviscosity, "viscosity","O","Pa-s","Viscosity",false),
    parameter_info(iconductivity, "conductivity","O","W/m/K","Thermal conductivity",false),
    parameter_info(isurface_tension, "surface_tension","O","N/m","Surface tension",false),
    parameter_info(iPrandtl, "Prandtl","O","-","Prandtl number",false),
    
    parameter_info(iisothermal_compressibility, "isothermal_compressibility","O","1/Pa","Isothermal compressibility",false),
    parameter_info(iisobaric_expansion_coefficient, "isobaric_expansion_coefficient","O","1/K","Isobaric expansion coefficient",false),
    parameter_info(iZ, "Z","O","-","Compressibility factor",false),
    parameter_info(ifundamental_derivative_of_gas_dynamics, "fundamental_derivative_of_gas_dynamics","O","-","Fundamental_derivative_of_gas_dynamics",false),
    
    parameter_info(ialphar, "alphar","O","-","Residual Helmholtz energy",false),
    parameter_info(idalphar_dtau_constdelta, "dalphar_dtau_constdelta","O","-","Derivative of residual Helmholtz energy with tau",false),
    parameter_info(idalphar_ddelta_consttau, "dalphar_ddelta_consttau","O","-","Derivative of residual Helmholtz energy with delta",false),
    
    parameter_info(ialpha0, "alpha0","O","-","Ideal Helmholtz energy",false),
    parameter_info(idalpha0_dtau_constdelta, "dalpha0_dtau_constdelta","O","-","Derivative of ideal Helmholtz energy with tau",false),
    parameter_info(idalpha0_ddelta_consttau, "dalpha0_ddelta_consttau","O","-","Derivative of ideal Helmholtz energy with delta",false),
    
    parameter_info(iPhase, "Phase","O","-","Phase index as a float",false),
    
};

class ParameterInformation
{
public:
    std::map<int, bool> trivial_map;
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
            trivial_map.insert(std::pair<int, bool>(el.key, el.trivial));
        }
        // Backward compatibility aliases
        index_map.insert(std::pair<std::string, int>("D", iDmass));
        index_map.insert(std::pair<std::string, int>("H", iHmass));
        index_map.insert(std::pair<std::string, int>("M", imolar_mass));
        index_map.insert(std::pair<std::string, int>("S", iSmass));
        index_map.insert(std::pair<std::string, int>("U", iUmass));
        index_map.insert(std::pair<std::string, int>("C", iCpmass));
        index_map.insert(std::pair<std::string, int>("O", iCvmass));
        index_map.insert(std::pair<std::string, int>("V", iviscosity));
        index_map.insert(std::pair<std::string, int>("L", iconductivity));
        index_map.insert(std::pair<std::string, int>("pcrit", iP_critical));
        index_map.insert(std::pair<std::string, int>("Pcrit", iP_critical));
        index_map.insert(std::pair<std::string, int>("Tcrit", iT_critical));
        index_map.insert(std::pair<std::string, int>("Ttriple", iT_triple));
        index_map.insert(std::pair<std::string, int>("ptriple", iP_triple));
        index_map.insert(std::pair<std::string, int>("rhocrit", irhomass_critical));
        index_map.insert(std::pair<std::string, int>("Tmin", iT_min));
        index_map.insert(std::pair<std::string, int>("Tmax", iT_max));
        index_map.insert(std::pair<std::string, int>("pmax", iP_max));
        index_map.insert(std::pair<std::string, int>("pmin", iP_min));
        index_map.insert(std::pair<std::string, int>("molemass", imolar_mass));
        index_map.insert(std::pair<std::string, int>("molarmass", imolar_mass));
        index_map.insert(std::pair<std::string, int>("A", ispeed_sound));
		index_map.insert(std::pair<std::string, int>("I", isurface_tension));

        std::map<std::string,int>::iterator it;
        for(it = index_map.begin(); it != index_map.end(); ++it )
        {
            // Add all upper-case aliases for EES support (fine to just do it if is already there)
            index_map.insert(std::pair<std::string, int>(upper(it->first), it->second));
        }
    }
};

static ParameterInformation parameter_information;

bool is_trivial_parameter(int key)
{
    std::map<int, bool>::iterator it;

    // Try to find it
    it = parameter_information.trivial_map.find(key);
    // If equal to end, not found
    if (it != parameter_information.trivial_map.end())
    {
        // Found it, return it
        return it->second;
    }
    else
    {
        throw ValueError(format("Unable to match the key [%d: %s] in is_trivial_parameter",key, get_parameter_information(key, "short").c_str()));
    }
}

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

/// Return a list of parameters
std::string get_csv_parameter_list()
{
    std::vector<std::string> strings;
    std::map<std::string,int>::iterator it;
    for(it = parameter_information.index_map.begin(); it != parameter_information.index_map.end(); ++it )
    {
        strings.push_back(it->first);
    }
    return strjoin(strings, ",");
}
bool is_valid_parameter(const std::string &param_name, parameters &iOutput)
{
    std::map<std::string, int>::iterator it;
    
    // Try to find it
    it = parameter_information.index_map.find(param_name);
    // If equal to end, not found
    if (it != parameter_information.index_map.end()){
        // Found, return it
        iOutput = static_cast<parameters>(it->second);
        return true;
    }
    else{
        return false;
    }
}

bool is_valid_first_derivative(const std::string & name, parameters &iOf, parameters &iWrt, parameters &iConstant)
{
    std::size_t iN0, iN1, iD0, iD1;
    parameters Of, Wrt, Constant;
    if (get_debug_level() > 5){std::cout << format("is_valid_first_derivative(%s)",name.c_str());}
    // There should be exactly one /
    // There should be exactly one |
    
    // Suppose we start with "d(P)/d(T)|Dmolar"
    std::vector<std::string> split_at_bar = strsplit(name, '|'); // "d(P)/d(T)"  and "Dmolar"
    if (split_at_bar.size() != 2){return false;}
    
    std::vector<std::string> split_at_slash = strsplit(split_at_bar[0], '/'); // "d(P)" and "d(T)"
    if (split_at_slash.size() != 2){return false;}
    
    iN0 = split_at_slash[0].find("(");
    iN1 = split_at_slash[0].find(")", iN0);
    if (!(iN0 > 0 && iN1 > 0 && iN1 > iN0)){return false;}
    std::string num = split_at_slash[0].substr(iN0+1, iN1-2);
    
    iD0 = split_at_slash[1].find("(");
    iD1 = split_at_slash[1].find(")", iD0);
    if (!(iD0 > 0 && iD1 > 0 && iD1 > iD0)){return false;}
    std::string den = split_at_slash[1].substr(iD0+1, iD1-2);
    
    if (is_valid_parameter(num, Of) && is_valid_parameter(den, Wrt) && is_valid_parameter(split_at_bar[1], Constant)){
        iOf = Of; iWrt = Wrt; iConstant = Constant; return true;
    }
    else{
        return false;
    }
}

bool is_valid_second_derivative(const std::string & name, parameters &iOf1, parameters &iWrt1, parameters &iConstant1, parameters &iWrt2, parameters &iConstant2)
{
    if (get_debug_level() > 5){std::cout << format("is_valid_second_derivative(%s)",name.c_str());}
    
    // Suppose we start with "d(d(P)/d(Dmolar)|T)/d(Dmolar)|T"
    std::size_t i_bar = name.rfind('|');
    if (i_bar == std::string::npos){return false;}
    std::string constant2 = name.substr(i_bar+1); // "T"
    if (!is_valid_parameter(constant2, iConstant2)){return false;};
    std::string left_of_bar = name.substr(0, i_bar); // "d(d(P)/d(Dmolar)|T)/d(Dmolar)"
    
    std::size_t i_slash = left_of_bar.rfind('/');
    if (i_slash == std::string::npos){return false;}
    
    std::string left_of_slash = left_of_bar.substr(0, i_slash); // "d(d(P)/d(Dmolar)|T)"
    std::size_t iN0 = left_of_slash.find("(");
    std::size_t iN1 = left_of_slash.rfind(")");
    if (!(iN0 > 0 && iN1 > 0 && iN1 > iN0)){return false;}
    std::string num = left_of_slash.substr(iN0+1, iN1-2); // "d(P)/d(Dmolar)|T"
    if (!is_valid_first_derivative(num, iOf1, iWrt1, iConstant1)){return false;}
    
    std::string right_of_slash = left_of_bar.substr(i_slash+1); // "d(Dmolar)"
    std::size_t iD0 = right_of_slash.find("(");
    std::size_t iD1 = right_of_slash.rfind(")");
    if (!(iD0 > 0 && iD1 > 0 && iD1 > iD0)){return false;}
    std::string den = right_of_slash.substr(iD0+1, iD1-2); // "Dmolar"
    if (!is_valid_parameter(den, iWrt2)){return false;}
    
    // If we haven't quit yet, all is well
    return true;
}

struct phase_info
{    
    phases key;
    std::string short_desc, long_desc;
public:
    phase_info(phases key, std::string short_desc, std::string long_desc): key(key), short_desc(short_desc), long_desc(long_desc){};
};
            
phase_info phase_info_list[] = {
    phase_info(iphase_liquid, "phase_liquid",""),
    phase_info(iphase_gas, "phase_gas",""),
    phase_info(iphase_twophase, "phase_twophase",""),
    phase_info(iphase_supercritical, "phase_supercritical", ""),
    phase_info(iphase_supercritical_gas, "phase_supercritical_gas", "p < pc, T > Tc"),
    phase_info(iphase_supercritical_liquid, "phase_supercritical_liquid", "p > pc, T < Tc"),
    phase_info(iphase_critical_point, "phase_critical_point", "p = pc, T = Tc"),
    phase_info(iphase_unknown, "phase_unknown", ""),
    phase_info(iphase_not_imposed, "phase_not_imposed", ""),
};

class PhaseInformation
{
public:
    std::map<phases, std::string> short_desc_map, long_desc_map;
    std::map<std::string, phases> index_map;
    PhaseInformation()
    {
        int N = sizeof(phase_info_list)/sizeof(phase_info_list[0]);
        for (int i = 0; i < N; ++i)
        {
            phase_info &el = phase_info_list[i];
            short_desc_map.insert(std::pair<phases, std::string>(el.key, el.short_desc));
            long_desc_map.insert(std::pair<phases, std::string>(el.key, el.long_desc));
            index_map.insert(std::pair<std::string, phases>(el.short_desc, el.key));
        }
    }
};
static PhaseInformation phase_information;

std::string get_phase_short_desc(phases phase)
{
    return phase_information.short_desc_map[phase];
}
bool is_valid_phase(const std::string &phase_name, phases &iOutput)
{
    std::map<std::string, phases>::iterator it;
    
    // Try to find it
    it = phase_information.index_map.find(phase_name);
    // If equal to end, not found
    if (it != phase_information.index_map.end()){
        // Found, return it
        iOutput = static_cast<phases>(it->second);
        return true;
    }
    else{
        return false;
    }
}

phases get_phase_index(const std::string &param_name)
{
    phases iPhase;
    if (is_valid_phase(param_name, iPhase)){
        return iPhase;
    }
    else{
        throw ValueError(format("Your input name [%s] is not valid in get_phase_index (names are case sensitive)",param_name.c_str()));
    }
}
parameters get_parameter_index(const std::string &param_name)
{
    parameters iOutput;
    if (is_valid_parameter(param_name, iOutput)){
        return iOutput;
    }
    else{
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
    input_pair_info(QSmolar_INPUTS,"QS_INPUTS","Molar quality, Entropy in J/mol/K"),
    input_pair_info(QSmass_INPUTS,"QS_INPUTS","Molar quality, Entropy in J/kg/K"),
    input_pair_info(HmolarQ_INPUTS,"HQ_INPUTS","Enthalpy in J/mol, Molar quality"),
    input_pair_info(HmassQ_INPUTS,"HQ_INPUTS","Enthalpy in J/kg, Molar quality"),
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



#ifdef ENABLE_CATCH
#include "catch.hpp"

TEST_CASE("Check that all parameters are descibed","")
{
    for (int i = 1; i < CoolProp::iundefined_parameter; ++i){
        std::ostringstream ss;
        ss << "Parameter index," << i << "last index:" << CoolProp::iundefined_parameter;
        SECTION(ss.str(), "")
        {   
            std::string prior;
            if (i > 1){
                CHECK_NOTHROW(prior = CoolProp::get_parameter_information(i-1,"short"));
                CAPTURE(prior);
            }
            CHECK_NOTHROW(CoolProp::get_parameter_information(i,"short"));
        }
    }
}

TEST_CASE("Check that all phases are descibed","[phase_index]")
{
    for (int i = 0; i < CoolProp::iphase_not_imposed; ++i){
        std::ostringstream ss;
        ss << "Parameter index," << i << "last index:" << CoolProp::iundefined_parameter;
        SECTION(ss.str(), "")
        {   
            std::string stringrepr;
            int key;
            CHECK_NOTHROW(stringrepr = CoolProp::get_phase_short_desc(static_cast<CoolProp::phases>(i)));
            CAPTURE(stringrepr);
            CHECK_NOTHROW(key = CoolProp::get_phase_index(stringrepr));
            CAPTURE(key);
            CHECK(key == i);
        }
    }
}

#endif

