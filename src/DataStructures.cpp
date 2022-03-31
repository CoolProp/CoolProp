

#include "DataStructures.h"
#include "Exceptions.h"
#include "CoolPropTools.h"
#include "CoolProp.h"

namespace CoolProp {

struct parameter_info
{
    int key;
    const char *short_desc, *IO, *units, *description;
    bool trivial;  ///< True if the input is trivial, and can be directly calculated (constants like critical properties, etc.)
};

const parameter_info parameter_info_list[] = {
  /// Input/Output parameters
  {iT, "T", "IO", "K", "Temperature", false},
  {iP, "P", "IO", "Pa", "Pressure", false},
  {iDmolar, "Dmolar", "IO", "mol/m^3", "Molar density", false},
  {iHmolar, "Hmolar", "IO", "J/mol", "Molar specific enthalpy", false},
  {iSmolar, "Smolar", "IO", "J/mol/K", "Molar specific entropy", false},
  {iUmolar, "Umolar", "IO", "J/mol", "Molar specific internal energy", false},
  {iGmolar, "Gmolar", "O", "J/mol", "Molar specific Gibbs energy", false},
  {iHelmholtzmolar, "Helmholtzmolar", "O", "J/mol", "Molar specific Helmholtz energy", false},
  {iDmass, "Dmass", "IO", "kg/m^3", "Mass density", false},
  {iHmass, "Hmass", "IO", "J/kg", "Mass specific enthalpy", false},
  {iSmass, "Smass", "IO", "J/kg/K", "Mass specific entropy", false},
  {iUmass, "Umass", "IO", "J/kg", "Mass specific internal energy", false},
  {iGmass, "Gmass", "O", "J/kg", "Mass specific Gibbs energy", false},
  {iHelmholtzmass, "Helmholtzmass", "O", "J/kg", "Mass specific Helmholtz energy", false},
  {iQ, "Q", "IO", "mol/mol", "Molar vapor quality", false},
  {iDelta, "Delta", "IO", "-", "Reduced density (rho/rhoc)", false},
  {iTau, "Tau", "IO", "-", "Reciprocal reduced temperature (Tc/T)", false},
  /// Output only
  {iCpmolar, "Cpmolar", "O", "J/mol/K", "Molar specific constant pressure specific heat", false},
  {iCpmass, "Cpmass", "O", "J/kg/K", "Mass specific constant pressure specific heat", false},
  {iCvmolar, "Cvmolar", "O", "J/mol/K", "Molar specific constant volume specific heat", false},
  {iCvmass, "Cvmass", "O", "J/kg/K", "Mass specific constant volume specific heat", false},
  {iCp0molar, "Cp0molar", "O", "J/mol/K", "Ideal gas molar specific constant pressure specific heat", false},
  {iCp0mass, "Cp0mass", "O", "J/kg/K", "Ideal gas mass specific constant pressure specific heat", false},
  {iHmolar_residual, "Hmolar_residual", "O", "J/mol/K", "Residual molar enthalpy", false},
  {iSmolar_residual, "Smolar_residual", "O", "J/mol/K", "Residual molar entropy (sr/R = s(T,rho) - s^0(T,rho))", false},
  {iGmolar_residual, "Gmolar_residual", "O", "J/mol/K", "Residual molar Gibbs energy", false},
  {iGWP20, "GWP20", "O", "-", "20-year global warming potential", true},
  {iGWP100, "GWP100", "O", "-", "100-year global warming potential", true},
  {iGWP500, "GWP500", "O", "-", "500-year global warming potential", true},
  {iFH, "FH", "O", "-", "Flammability hazard", true},
  {iHH, "HH", "O", "-", "Health hazard", true},
  {iPH, "PH", "O", "-", "Physical hazard", true},
  {iODP, "ODP", "O", "-", "Ozone depletion potential", true},
  {iBvirial, "Bvirial", "O", "-", "Second virial coefficient", false},
  {iCvirial, "Cvirial", "O", "-", "Third virial coefficient", false},
  {idBvirial_dT, "dBvirial_dT", "O", "-", "Derivative of second virial coefficient with respect to T", false},
  {idCvirial_dT, "dCvirial_dT", "O", "-", "Derivative of third virial coefficient with respect to T", false},
  {igas_constant, "gas_constant", "O", "J/mol/K", "Molar gas constant", true},
  {imolar_mass, "molar_mass", "O", "kg/mol", "Molar mass", true},
  {iacentric_factor, "acentric", "O", "-", "Acentric factor", true},
  {idipole_moment, "dipole_moment", "O", "C-m", "Dipole moment", true},
  {irhomass_reducing, "rhomass_reducing", "O", "kg/m^3", "Mass density at reducing point", true},
  {irhomolar_reducing, "rhomolar_reducing", "O", "mol/m^3", "Molar density at reducing point", true},
  {irhomolar_critical, "rhomolar_critical", "O", "mol/m^3", "Molar density at critical point", true},
  {irhomass_critical, "rhomass_critical", "O", "kg/m^3", "Mass density at critical point", true},
  {iT_reducing, "T_reducing", "O", "K", "Temperature at the reducing point", true},
  {iT_critical, "T_critical", "O", "K", "Temperature at the critical point", true},
  {iT_triple, "T_triple", "O", "K", "Temperature at the triple point", true},
  {iT_max, "T_max", "O", "K", "Maximum temperature limit", true},
  {iT_min, "T_min", "O", "K", "Minimum temperature limit", true},
  {iP_min, "P_min", "O", "Pa", "Minimum pressure limit", true},
  {iP_max, "P_max", "O", "Pa", "Maximum pressure limit", true},
  {iP_critical, "p_critical", "O", "Pa", "Pressure at the critical point", true},
  {iP_reducing, "p_reducing", "O", "Pa", "Pressure at the reducing point", true},
  {iP_triple, "p_triple", "O", "Pa", "Pressure at the triple point (pure only)", true},
  {ifraction_min, "fraction_min", "O", "-", "Fraction (mole, mass, volume) minimum value for incompressible solutions", true},
  {ifraction_max, "fraction_max", "O", "-", "Fraction (mole, mass, volume) maximum value for incompressible solutions", true},
  {iT_freeze, "T_freeze", "O", "K", "Freezing temperature for incompressible solutions", true},

  {ispeed_sound, "speed_of_sound", "O", "m/s", "Speed of sound", false},
  {iviscosity, "viscosity", "O", "Pa-s", "Viscosity", false},
  {iconductivity, "conductivity", "O", "W/m/K", "Thermal conductivity", false},
  {isurface_tension, "surface_tension", "O", "N/m", "Surface tension", false},
  {iPrandtl, "Prandtl", "O", "-", "Prandtl number", false},

  {iisothermal_compressibility, "isothermal_compressibility", "O", "1/Pa", "Isothermal compressibility", false},
  {iisobaric_expansion_coefficient, "isobaric_expansion_coefficient", "O", "1/K", "Isobaric expansion coefficient", false},
  {iisentropic_expansion_coefficient, "isentropic_expansion_coefficient", "O", "-", "Isentropic expansion coefficient", false},
  {iZ, "Z", "O", "-", "Compressibility factor", false},
  {ifundamental_derivative_of_gas_dynamics, "fundamental_derivative_of_gas_dynamics", "O", "-", "Fundamental derivative of gas dynamics", false},
  {iPIP, "PIP", "O", "-", "Phase identification parameter", false},

  {ialphar, "alphar", "O", "-", "Residual Helmholtz energy", false},
  {idalphar_dtau_constdelta, "dalphar_dtau_constdelta", "O", "-", "Derivative of residual Helmholtz energy with tau", false},
  {idalphar_ddelta_consttau, "dalphar_ddelta_consttau", "O", "-", "Derivative of residual Helmholtz energy with delta", false},

  {ialpha0, "alpha0", "O", "-", "Ideal Helmholtz energy", false},
  {idalpha0_dtau_constdelta, "dalpha0_dtau_constdelta", "O", "-", "Derivative of ideal Helmholtz energy with tau", false},
  {idalpha0_ddelta_consttau, "dalpha0_ddelta_consttau", "O", "-", "Derivative of ideal Helmholtz energy with delta", false},
  {id2alpha0_ddelta2_consttau, "d2alpha0_ddelta2_consttau", "O", "-", "Second derivative of ideal Helmholtz energy with delta", false},
  {id3alpha0_ddelta3_consttau, "d3alpha0_ddelta3_consttau", "O", "-", "Third derivative of ideal Helmholtz energy with delta", false},

  {iPhase, "Phase", "O", "-", "Phase index as a float", false},

};

class ParameterInformation
{
   public:
    std::map<int, bool> trivial_map;
    std::map<int, std::string> short_desc_map, description_map, IO_map, units_map;
    std::map<std::string, int> index_map;
    ParameterInformation() {
        const parameter_info* const end = parameter_info_list + sizeof(parameter_info_list) / sizeof(parameter_info_list[0]);
        for (const parameter_info* el = parameter_info_list; el != end; ++el) {
            short_desc_map.insert(std::pair<int, std::string>(el->key, el->short_desc));
            IO_map.insert(std::pair<int, std::string>(el->key, el->IO));
            units_map.insert(std::pair<int, std::string>(el->key, el->units));
            description_map.insert(std::pair<int, std::string>(el->key, el->description));
            index_map_insert(el->short_desc, el->key);
            trivial_map.insert(std::pair<int, bool>(el->key, el->trivial));
        }
        // Backward compatibility aliases
        index_map_insert("D", iDmass);
        index_map_insert("H", iHmass);
        index_map_insert("M", imolar_mass);
        index_map_insert("S", iSmass);
        index_map_insert("U", iUmass);
        index_map_insert("C", iCpmass);
        index_map_insert("O", iCvmass);
        index_map_insert("G", iGmass);
        index_map_insert("V", iviscosity);
        index_map_insert("L", iconductivity);
        index_map_insert("pcrit", iP_critical);
        index_map_insert("Pcrit", iP_critical);
        index_map_insert("Tcrit", iT_critical);
        index_map_insert("Ttriple", iT_triple);
        index_map_insert("ptriple", iP_triple);
        index_map_insert("rhocrit", irhomass_critical);
        index_map_insert("Tmin", iT_min);
        index_map_insert("Tmax", iT_max);
        index_map_insert("pmax", iP_max);
        index_map_insert("pmin", iP_min);
        index_map_insert("molemass", imolar_mass);
        index_map_insert("molarmass", imolar_mass);
        index_map_insert("A", ispeed_sound);
        index_map_insert("I", isurface_tension);
    }

   private:
    void index_map_insert(const std::string& desc, int key) {
        index_map.insert(std::pair<std::string, int>(desc, key));
        index_map.insert(std::pair<std::string, int>(upper(desc), key));
    }
};

static ParameterInformation parameter_information;

bool is_trivial_parameter(int key) {
    // Try to find it
    std::map<int, bool>::const_iterator it = parameter_information.trivial_map.find(key);
    // If equal to end, not found
    if (it != parameter_information.trivial_map.end()) {
        // Found it, return it
        return it->second;
    } else {
        throw ValueError(format("Unable to match the key [%d: %s] in is_trivial_parameter", key, get_parameter_information(key, "short").c_str()));
    }
}

std::string get_parameter_information(int key, const std::string& info) {
    std::map<int, std::string>* M;

    // Hook up the right map (since they are all of the same type)
    if (!info.compare("IO")) {
        M = &(parameter_information.IO_map);
    } else if (!info.compare("short")) {
        M = &(parameter_information.short_desc_map);
    } else if (!info.compare("long")) {
        M = &(parameter_information.description_map);
    } else if (!info.compare("units")) {
        M = &(parameter_information.units_map);
    } else
        throw ValueError(format("Bad info string [%s] to get_parameter_information", info.c_str()));

    // Try to find it
    std::map<int, std::string>::const_iterator it = M->find(key);
    // If equal to end, not found
    if (it != M->end()) {
        // Found it, return it
        return it->second;
    } else {
        throw ValueError(format("Unable to match the key [%d] in get_parameter_information for info [%s]", key, info.c_str()));
    }
}

/// Return a list of parameters
std::string get_csv_parameter_list() {
    std::vector<std::string> strings;
    for (std::map<std::string, int>::const_iterator it = parameter_information.index_map.begin(); it != parameter_information.index_map.end(); ++it) {
        strings.push_back(it->first);
    }
    return strjoin(strings, ",");
}
bool is_valid_parameter(const std::string& param_name, parameters& iOutput) {
    // Try to find it
    std::map<std::string, int>::const_iterator it = parameter_information.index_map.find(param_name);
    // If equal to end, not found
    if (it != parameter_information.index_map.end()) {
        // Found, return it
        iOutput = static_cast<parameters>(it->second);
        return true;
    } else {
        return false;
    }
}

bool is_valid_first_derivative(const std::string& name, parameters& iOf, parameters& iWrt, parameters& iConstant) {
    if (get_debug_level() > 5) {
        std::cout << format("is_valid_first_derivative(%s)", name.c_str());
    }
    // There should be exactly one /
    // There should be exactly one |

    // Suppose we start with "d(P)/d(T)|Dmolar"
    std::vector<std::string> split_at_bar = strsplit(name, '|');  // "d(P)/d(T)"  and "Dmolar"
    if (split_at_bar.size() != 2) {
        return false;
    }

    std::vector<std::string> split_at_slash = strsplit(split_at_bar[0], '/');  // "d(P)" and "d(T)"
    if (split_at_slash.size() != 2) {
        return false;
    }

    std::size_t i0 = split_at_slash[0].find("(");
    std::size_t i1 = split_at_slash[0].find(")", i0);
    if (!((i0 > 0) && (i0 != std::string::npos) && (i1 > (i0 + 1)) && (i1 != std::string::npos))) {
        return false;
    }
    std::string num = split_at_slash[0].substr(i0 + 1, i1 - i0 - 1);

    i0 = split_at_slash[1].find("(");
    i1 = split_at_slash[1].find(")", i0);
    if (!((i0 > 0) && (i0 != std::string::npos) && (i1 > (i0 + 1)) && (i1 != std::string::npos))) {
        return false;
    }
    std::string den = split_at_slash[1].substr(i0 + 1, i1 - i0 - 1);

    parameters Of, Wrt, Constant;
    if (is_valid_parameter(num, Of) && is_valid_parameter(den, Wrt) && is_valid_parameter(split_at_bar[1], Constant)) {
        iOf = Of;
        iWrt = Wrt;
        iConstant = Constant;
        return true;
    } else {
        return false;
    }
}

bool is_valid_first_saturation_derivative(const std::string& name, parameters& iOf, parameters& iWrt) {
    if (get_debug_level() > 5) {
        std::cout << format("is_valid_first_saturation_derivative(%s)", name.c_str());
    }
    // There should be exactly one /
    // There should be exactly one |

    // Suppose we start with "d(P)/d(T)|sigma"
    std::vector<std::string> split_at_bar = strsplit(name, '|');  // "d(P)/d(T)"  and "sigma"
    if (split_at_bar.size() != 2) {
        return false;
    }

    std::vector<std::string> split_at_slash = strsplit(split_at_bar[0], '/');  // "d(P)" and "d(T)"
    if (split_at_slash.size() != 2) {
        return false;
    }

    std::size_t i0 = split_at_slash[0].find("(");
    std::size_t i1 = split_at_slash[0].find(")", i0);
    if (!((i0 > 0) && (i0 != std::string::npos) && (i1 > (i0 + 1)) && (i1 != std::string::npos))) {
        return false;
    }
    std::string num = split_at_slash[0].substr(i0 + 1, i1 - i0 - 1);

    i0 = split_at_slash[1].find("(");
    i1 = split_at_slash[1].find(")", i0);
    if (!((i0 > 0) && (i0 != std::string::npos) && (i1 > (i0 + 1)) && (i1 != std::string::npos))) {
        return false;
    }
    std::string den = split_at_slash[1].substr(i0 + 1, i1 - i0 - 1);

    parameters Of, Wrt;
    if (is_valid_parameter(num, Of) && is_valid_parameter(den, Wrt) && upper(split_at_bar[1]) == "SIGMA") {
        iOf = Of;
        iWrt = Wrt;
        return true;
    } else {
        return false;
    }
}

bool is_valid_second_derivative(const std::string& name, parameters& iOf1, parameters& iWrt1, parameters& iConstant1, parameters& iWrt2,
                                parameters& iConstant2) {
    if (get_debug_level() > 5) {
        std::cout << format("is_valid_second_derivative(%s)", name.c_str());
    }

    // Suppose we start with "d(d(P)/d(Dmolar)|T)/d(Dmolar)|T"
    std::size_t i = name.rfind('|');
    if ((i == 0) || (i == std::string::npos)) {
        return false;
    }
    std::string constant2 = name.substr(i + 1);  // "T"
    if (!is_valid_parameter(constant2, iConstant2)) {
        return false;
    };
    std::string left_of_bar = name.substr(0, i);  // "d(d(P)/d(Dmolar)|T)/d(Dmolar)"

    i = left_of_bar.rfind('/');
    if ((i == 0) || (i == std::string::npos)) {
        return false;
    }
    std::string left_of_slash = left_of_bar.substr(0, i);    // "d(d(P)/d(Dmolar)|T)"
    std::string right_of_slash = left_of_bar.substr(i + 1);  // "d(Dmolar)"

    i = left_of_slash.find("(");
    std::size_t i1 = left_of_slash.rfind(")");
    if (!((i > 0) && (i != std::string::npos) && (i1 > (i + 1)) && (i1 != std::string::npos))) {
        return false;
    }
    std::string num = left_of_slash.substr(i + 1, i1 - i - 1);  // "d(P)/d(Dmolar)|T"
    if (!is_valid_first_derivative(num, iOf1, iWrt1, iConstant1)) {
        return false;
    }

    i = right_of_slash.find("(");
    i1 = right_of_slash.rfind(")");
    if (!((i > 0) && (i != std::string::npos) && (i1 > (i + 1)) && (i1 != std::string::npos))) {
        return false;
    }
    std::string den = right_of_slash.substr(i + 1, i1 - i - 1);  // "Dmolar"
    if (!is_valid_parameter(den, iWrt2)) {
        return false;
    }

    // If we haven't quit yet, all is well
    return true;
}

struct phase_info
{
    phases key;
    const char *short_desc, *long_desc;
};

const phase_info phase_info_list[] = {
  {iphase_liquid, "phase_liquid", ""},
  {iphase_gas, "phase_gas", ""},
  {iphase_twophase, "phase_twophase", ""},
  {iphase_supercritical, "phase_supercritical", ""},
  {iphase_supercritical_gas, "phase_supercritical_gas", "p < pc, T > Tc"},
  {iphase_supercritical_liquid, "phase_supercritical_liquid", "p > pc, T < Tc"},
  {iphase_critical_point, "phase_critical_point", "p = pc, T = Tc"},
  {iphase_unknown, "phase_unknown", ""},
  {iphase_not_imposed, "phase_not_imposed", ""},
};

class PhaseInformation
{
   public:
    std::map<phases, std::string> short_desc_map, long_desc_map;
    std::map<std::string, phases> index_map;
    PhaseInformation() {
        const phase_info* const end = phase_info_list + sizeof(phase_info_list) / sizeof(phase_info_list[0]);
        for (const phase_info* el = phase_info_list; el != end; ++el) {
            short_desc_map.insert(std::pair<phases, std::string>(el->key, el->short_desc));
            long_desc_map.insert(std::pair<phases, std::string>(el->key, el->long_desc));
            index_map.insert(std::pair<std::string, phases>(el->short_desc, el->key));
        }
    }
};
static PhaseInformation phase_information;

const std::string& get_phase_short_desc(phases phase) {
    return phase_information.short_desc_map[phase];
}
bool is_valid_phase(const std::string& phase_name, phases& iOutput) {
    // Try to find it
    std::map<std::string, phases>::const_iterator it = phase_information.index_map.find(phase_name);
    // If equal to end, not found
    if (it != phase_information.index_map.end()) {
        // Found, return it
        iOutput = static_cast<phases>(it->second);
        return true;
    } else {
        return false;
    }
}

phases get_phase_index(const std::string& param_name) {
    phases iPhase;
    if (is_valid_phase(param_name, iPhase)) {
        return iPhase;
    } else {
        throw ValueError(format("Your input name [%s] is not valid in get_phase_index (names are case sensitive)", param_name.c_str()));
    }
}
parameters get_parameter_index(const std::string& param_name) {
    parameters iOutput;
    if (is_valid_parameter(param_name, iOutput)) {
        return iOutput;
    } else {
        throw ValueError(format("Your input name [%s] is not valid in get_parameter_index (names are case sensitive)", param_name.c_str()));
    }
}

struct input_pair_info
{
    input_pairs key;
    const char *short_desc, *long_desc;
};

const input_pair_info input_pair_list[] = {
  {QT_INPUTS, "QT_INPUTS", "Molar quality, Temperature in K"},
  {QSmolar_INPUTS, "QS_INPUTS", "Molar quality, Entropy in J/mol/K"},
  {QSmass_INPUTS, "QS_INPUTS", "Molar quality, Entropy in J/kg/K"},
  {HmolarQ_INPUTS, "HQ_INPUTS", "Enthalpy in J/mol, Molar quality"},
  {HmassQ_INPUTS, "HQ_INPUTS", "Enthalpy in J/kg, Molar quality"},
  {DmassQ_INPUTS, "DmassQ_INPUTS", "Molar density kg/m^3, Molar quality"},
  {DmolarQ_INPUTS, "DmolarQ_INPUTS", "Molar density in mol/m^3, Molar quality"},

  {PQ_INPUTS, "PQ_INPUTS", "Pressure in Pa, Molar quality"},

  {PT_INPUTS, "PT_INPUTS", "Pressure in Pa, Temperature in K"},

  {DmassT_INPUTS, "DmassT_INPUTS", "Mass density in kg/m^3, Temperature in K"},
  {DmolarT_INPUTS, "DmolarT_INPUTS", "Molar density in mol/m^3, Temperature in K"},
  {HmassT_INPUTS, "HmassT_INPUTS", "Enthalpy in J/kg, Temperature in K"},
  {HmolarT_INPUTS, "HmolarT_INPUTS", "Enthalpy in J/mol, Temperature in K"},
  {SmassT_INPUTS, "SmassT_INPUTS", "Entropy in J/kg/K, Temperature in K"},
  {SmolarT_INPUTS, "SmolarT_INPUTS", "Entropy in J/mol/K, Temperature in K"},
  {TUmass_INPUTS, "TUmass_INPUTS", "Temperature in K, Internal energy in J/kg"},
  {TUmolar_INPUTS, "TUmolar_INPUTS", "Temperature in K, Internal energy in J/mol"},

  {DmassP_INPUTS, "DmassP_INPUTS", "Mass density in kg/m^3, Pressure in Pa"},
  {DmolarP_INPUTS, "DmolarP_INPUTS", "Molar density in mol/m^3, Pressure in Pa"},
  {HmassP_INPUTS, "HmassP_INPUTS", "Enthalpy in J/kg, Pressure in Pa"},
  {HmolarP_INPUTS, "HmolarP_INPUTS", "Enthalpy in J/mol, Pressure in Pa"},
  {PSmass_INPUTS, "PSmass_INPUTS", "Pressure in Pa, Entropy in J/kg/K"},
  {PSmolar_INPUTS, "PSmolar_INPUTS", "Pressure in Pa, Entropy in J/mol/K "},
  {PUmass_INPUTS, "PUmass_INPUTS", "Pressure in Pa, Internal energy in J/kg"},
  {PUmolar_INPUTS, "PUmolar_INPUTS", "Pressure in Pa, Internal energy in J/mol"},

  {DmassHmass_INPUTS, "DmassHmass_INPUTS", "Mass density in kg/m^3, Enthalpy in J/kg"},
  {DmolarHmolar_INPUTS, "DmolarHmolar_INPUTS", "Molar density in mol/m^3, Enthalpy in J/mol"},
  {DmassSmass_INPUTS, "DmassSmass_INPUTS", "Mass density in kg/m^3, Entropy in J/kg/K"},
  {DmolarSmolar_INPUTS, "DmolarSmolar_INPUTS", "Molar density in mol/m^3, Entropy in J/mol/K"},
  {DmassUmass_INPUTS, "DmassUmass_INPUTS", "Mass density in kg/m^3, Internal energy in J/kg"},
  {DmolarUmolar_INPUTS, "DmolarUmolar_INPUTS", "Molar density in mol/m^3, Internal energy in J/mol"},

  {HmassSmass_INPUTS, "HmassSmass_INPUTS", "Enthalpy in J/kg, Entropy in J/kg/K"},
  {HmolarSmolar_INPUTS, "HmolarSmolar_INPUTS", "Enthalpy in J/mol, Entropy in J/mol/K"},
  {SmassUmass_INPUTS, "SmassUmass_INPUTS", "Entropy in J/kg/K, Internal energy in J/kg"},
  {SmolarUmolar_INPUTS, "SmolarUmolar_INPUTS", "Entropy in J/mol/K, Internal energy in J/mol"},
};

class InputPairInformation
{
   public:
    std::map<input_pairs, std::string> short_desc_map, long_desc_map;
    std::map<std::string, input_pairs> index_map;
    InputPairInformation() {
        const input_pair_info* const end = input_pair_list + sizeof(input_pair_list) / sizeof(input_pair_list[0]);
        for (const input_pair_info* el = input_pair_list; el != end; ++el) {
            short_desc_map.insert(std::pair<input_pairs, std::string>(el->key, el->short_desc));
            long_desc_map.insert(std::pair<input_pairs, std::string>(el->key, el->long_desc));
            index_map.insert(std::pair<std::string, input_pairs>(el->short_desc, el->key));
        }
    }
};

static InputPairInformation input_pair_information;

input_pairs get_input_pair_index(const std::string& input_pair_name) {
    std::map<std::string, input_pairs>::iterator it = input_pair_information.index_map.find(input_pair_name);
    if (it != input_pair_information.index_map.end()) {
        return it->second;
    } else {
        throw ValueError(format("Your input name [%s] is not valid in get_input_pair_index (names are case sensitive)", input_pair_name.c_str()));
    }
}

const std::string& get_input_pair_short_desc(input_pairs pair) {
    return input_pair_information.short_desc_map[pair];
}
const std::string& get_input_pair_long_desc(input_pairs pair) {
    return input_pair_information.long_desc_map[pair];
}
void split_input_pair(input_pairs pair, parameters& p1, parameters& p2) {
    switch (pair) {
        case QT_INPUTS:
            p1 = iQ;
            p2 = iT;
            break;
        case QSmolar_INPUTS:
            p1 = iQ;
            p2 = iSmolar;
            break;
        case QSmass_INPUTS:
            p1 = iQ;
            p2 = iSmass;
            break;
        case HmolarQ_INPUTS:
            p1 = iHmolar;
            p2 = iQ;
            break;
        case HmassQ_INPUTS:
            p1 = iHmass;
            p2 = iQ;
            break;
        case PQ_INPUTS:
            p1 = iP;
            p2 = iQ;
            break;
        case PT_INPUTS:
            p1 = iP;
            p2 = iT;
            break;
        case DmassT_INPUTS:
            p1 = iDmass;
            p2 = iT;
            break;
        case DmolarT_INPUTS:
            p1 = iDmolar;
            p2 = iT;
            break;
        case HmassT_INPUTS:
            p1 = iHmass;
            p2 = iT;
            break;
        case HmolarT_INPUTS:
            p1 = iHmolar;
            p2 = iT;
            break;
        case SmassT_INPUTS:
            p1 = iSmass;
            p2 = iT;
            break;
        case SmolarT_INPUTS:
            p1 = iSmolar;
            p2 = iT;
            break;
        case TUmass_INPUTS:
            p1 = iT;
            p2 = iUmass;
            break;
        case TUmolar_INPUTS:
            p1 = iT;
            p2 = iUmolar;
            break;
        case DmassP_INPUTS:
            p1 = iDmass;
            p2 = iP;
            break;
        case DmolarP_INPUTS:
            p1 = iDmolar;
            p2 = iP;
            break;
        case DmassQ_INPUTS:
            p1 = iDmass;
            p2 = iQ;
            break;
        case DmolarQ_INPUTS:
            p1 = iDmolar;
            p2 = iQ;
            break;
        case HmassP_INPUTS:
            p1 = iHmass;
            p2 = iP;
            break;
        case HmolarP_INPUTS:
            p1 = iHmolar;
            p2 = iP;
            break;
        case PSmass_INPUTS:
            p1 = iP;
            p2 = iSmass;
            break;
        case PSmolar_INPUTS:
            p1 = iP;
            p2 = iSmolar;
            break;
        case PUmass_INPUTS:
            p1 = iP;
            p2 = iUmass;
            break;
        case PUmolar_INPUTS:
            p1 = iP;
            p2 = iUmolar;
            break;
        case DmassHmass_INPUTS:
            p1 = iDmass;
            p2 = iHmass;
            break;
        case DmolarHmolar_INPUTS:
            p1 = iDmolar;
            p2 = iHmolar;
            break;
        case DmassSmass_INPUTS:
            p1 = iDmass;
            p2 = iSmass;
            break;
        case DmolarSmolar_INPUTS:
            p1 = iDmolar;
            p2 = iSmolar;
            break;
        case DmassUmass_INPUTS:
            p1 = iDmass;
            p2 = iUmass;
            break;
        case DmolarUmolar_INPUTS:
            p1 = iDmolar;
            p2 = iUmolar;
            break;
        case HmassSmass_INPUTS:
            p1 = iHmass;
            p2 = iSmass;
            break;
        case HmolarSmolar_INPUTS:
            p1 = iHmolar;
            p2 = iSmolar;
            break;
        case SmassUmass_INPUTS:
            p1 = iSmass;
            p2 = iUmass;
            break;
        case SmolarUmolar_INPUTS:
            p1 = iSmolar;
            p2 = iUmolar;
            break;
        default:
            throw ValueError(format("Invalid input pair"));
    }
}

struct backend_family_info
{
    backend_families family;
    const char* name;
};

struct backend_info
{
    backends backend;
    const char* name;
    backend_families family;
};

const backend_family_info backend_family_list[] = {
  {HEOS_BACKEND_FAMILY, "HEOS"},   {REFPROP_BACKEND_FAMILY, "REFPROP"}, {INCOMP_BACKEND_FAMILY, "INCOMP"},   {IF97_BACKEND_FAMILY, "IF97"},
  {TREND_BACKEND_FAMILY, "TREND"}, {TTSE_BACKEND_FAMILY, "TTSE"},       {BICUBIC_BACKEND_FAMILY, "BICUBIC"}, {SRK_BACKEND_FAMILY, "SRK"},
  {PR_BACKEND_FAMILY, "PR"},       {VTPR_BACKEND_FAMILY, "VTPR"},       {PCSAFT_BACKEND_FAMILY, "PCSAFT"}};

const backend_info backend_list[] = {{HEOS_BACKEND_PURE, "HelmholtzEOSBackend", HEOS_BACKEND_FAMILY},
                                     {HEOS_BACKEND_MIX, "HelmholtzEOSMixtureBackend", HEOS_BACKEND_FAMILY},
                                     {REFPROP_BACKEND_PURE, "REFPROPBackend", REFPROP_BACKEND_FAMILY},
                                     {REFPROP_BACKEND_MIX, "REFPROPMixtureBackend", REFPROP_BACKEND_FAMILY},
                                     {INCOMP_BACKEND, "IncompressibleBackend", INCOMP_BACKEND_FAMILY},
                                     {IF97_BACKEND, "IF97Backend", IF97_BACKEND_FAMILY},
                                     {TREND_BACKEND, "TRENDBackend", TREND_BACKEND_FAMILY},
                                     {TTSE_BACKEND, "TTSEBackend", TTSE_BACKEND_FAMILY},
                                     {BICUBIC_BACKEND, "BicubicBackend", BICUBIC_BACKEND_FAMILY},
                                     {SRK_BACKEND, "SRKBackend", SRK_BACKEND_FAMILY},
                                     {PR_BACKEND, "PengRobinsonBackend", PR_BACKEND_FAMILY},
                                     {VTPR_BACKEND, "VTPRBackend", VTPR_BACKEND_FAMILY},
                                     {PCSAFT_BACKEND, "PCSAFTBackend", PCSAFT_BACKEND_FAMILY}};

class BackendInformation
{
   public:
    std::map<backend_families, std::string> family_name_map;  /// < from family to family name
    std::map<backends, backend_families> backend_family_map;  /// < from backend to family
    std::map<backends, std::string> backend_name_map;         /// < from backend to backend name

    std::map<std::string, backend_families> family_name_map_r;  /// < from backend name **or** family name to family
    std::map<std::string, backends> backend_name_map_r;         /// < from backend name to backend

    BackendInformation() {
        const backend_family_info* const family_end = backend_family_list + sizeof(backend_family_list) / sizeof(backend_family_list[0]);
        for (const backend_family_info* el = backend_family_list; el != family_end; ++el) {
            family_name_map.insert(std::pair<backend_families, std::string>(el->family, el->name));
            family_name_map_r.insert(std::pair<std::string, backend_families>(el->name, el->family));
        }
        const backend_info* const backend_end = backend_list + sizeof(backend_list) / sizeof(backend_list[0]);
        for (const backend_info* el = backend_list; el != backend_end; ++el) {
            backend_family_map.insert(std::pair<backends, backend_families>(el->backend, el->family));
            backend_name_map.insert(std::pair<backends, std::string>(el->backend, el->name));
            backend_name_map_r.insert(std::pair<std::string, backends>(el->name, el->backend));
            family_name_map_r.insert(std::pair<std::string, backend_families>(el->name, el->family));
        }
    }
};

static BackendInformation backend_information;

/// Convert a string into the enum values
void extract_backend_families(std::string backend_string, backend_families& f1, backend_families& f2) {
    f1 = INVALID_BACKEND_FAMILY;
    f2 = INVALID_BACKEND_FAMILY;
    std::size_t i = backend_string.find("&");
    std::map<std::string, backend_families>::const_iterator it;
    if (i != std::string::npos) {
        it = backend_information.family_name_map_r.find(backend_string.substr(0, i));  // Before "&"
        if (it != backend_information.family_name_map_r.end()) f1 = it->second;
        it = backend_information.family_name_map_r.find(backend_string.substr(i + 1));  // After "&"
        if (it != backend_information.family_name_map_r.end()) f2 = it->second;
    } else {
        it = backend_information.family_name_map_r.find(backend_string);
        if (it != backend_information.family_name_map_r.end()) f1 = it->second;
    }
}

void extract_backend_families_string(std::string backend_string, backend_families& f1, std::string& f2) {
    backend_families f2_enum;
    extract_backend_families(backend_string, f1, f2_enum);
    std::map<backend_families, std::string>::const_iterator it;
    it = backend_information.family_name_map.find(f2_enum);
    if (it != backend_information.family_name_map.end())
        f2 = it->second;
    else
        f2.clear();
}

std::string get_backend_string(backends backend) {
    std::map<backends, std::string>::const_iterator it;
    it = backend_information.backend_name_map.find(backend);
    if (it != backend_information.backend_name_map.end())
        return it->second;
    else
        return std::string("");
}

} /* namespace CoolProp */

#ifdef ENABLE_CATCH
#    include <catch2/catch_all.hpp>
#    include <sstream>

TEST_CASE("Check that all parameters are described", "") {
    for (int i = 1; i < CoolProp::iundefined_parameter; ++i) {
        std::ostringstream ss;
        ss << "Parameter index," << i << "last index:" << CoolProp::iundefined_parameter;
        SECTION(ss.str(), "") {
            std::string prior;
            if (i > 1) {
                CHECK_NOTHROW(prior = CoolProp::get_parameter_information(i - 1, "short"));
                CAPTURE(prior);
            }
            CHECK_NOTHROW(CoolProp::get_parameter_information(i, "short"));
        }
    }
}

TEST_CASE("Check that all phases are described", "[phase_index]") {
    for (int i = 0; i < CoolProp::iphase_not_imposed; ++i) {
        std::ostringstream ss;
        ss << "Parameter index," << i << "last index:" << CoolProp::iundefined_parameter;
        SECTION(ss.str(), "") {
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
