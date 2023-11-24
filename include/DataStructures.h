/*
 * DataStructures.h
 *
 *  Created on: 21 Dec 2013
 *      Author: jowr
 */

#ifndef DATASTRUCTURES_H_
#define DATASTRUCTURES_H_

#include "CPnumerics.h"
#include "Exceptions.h"
#include <map>
namespace CoolProp {

struct SimpleState
{
    double rhomolar, T, p, hmolar, smolar, umolar, Q;
    SimpleState() {
        fill(_HUGE);
    }
    void fill(double v) {
        rhomolar = v;
        T = v;
        p = v;
        hmolar = v;
        smolar = v;
        umolar = v;
        Q = v;
    }
    bool is_valid() {
        return ValidNumber(rhomolar) && ValidNumber(T) && ValidNumber(hmolar) && ValidNumber(p);
    }
};

struct CriticalState : SimpleState
{
    bool stable;
    CriticalState() : stable(false) {
        fill(_HUGE);
    }
};

/// A modified class for the state point at the maximum saturation entropy on the vapor curve
struct SsatSimpleState : public SimpleState
{
    enum SsatSimpleStateEnum
    {
        SSAT_MAX_NOT_SET = 0,
        SSAT_MAX_DOESNT_EXIST,
        SSAT_MAX_DOES_EXIST
    };
    SsatSimpleStateEnum exists;
    SsatSimpleState() : exists(SSAT_MAX_NOT_SET) {}
};

/// --------------------------------------------------
/// Define some constants that will be used throughout
/// --------------------------------------------------
/// These are constants for the input and output parameters
/// The structure is taken directly from the AbstractState class.
//
// !! If you add a parameter, update the map in the corresponding CPP file !!
enum parameters
{
    INVALID_PARAMETER = 0,

    // General parameters
    igas_constant,       ///< Ideal-gas constant
    imolar_mass,         ///< Molar mass
    iacentric_factor,    ///< Acentric factor
    irhomolar_reducing,  ///< Molar density used for the reducing state
    irhomolar_critical,  ///< Molar density used for the critical point
    iT_reducing,         ///< Temperature at the reducing state
    iT_critical,         ///< Temperature at the critical point
    irhomass_reducing,   ///< Mass density at the reducing state
    irhomass_critical,   ///< Mass density at the critical point
    iP_critical,         ///< Pressure at the critical point
    iP_reducing,         ///< Pressure at the reducing point
    iT_triple,           ///< Triple point temperature
    iP_triple,           ///< Triple point pressure
    iT_min,              ///< Minimum temperature
    iT_max,              ///< Maximum temperature
    iP_max,              ///< Maximum pressure
    iP_min,              ///< Minimum pressure
    idipole_moment,      ///< Dipole moment

    // Bulk properties
    iT,      ///< Temperature
    iP,      ///< Pressure
    iQ,      ///< Vapor quality
    iTau,    ///< Reciprocal reduced temperature
    iDelta,  ///< Reduced density

    // Molar specific thermodynamic properties
    iDmolar,           ///< Mole-based density
    iHmolar,           ///< Mole-based enthalpy
    iSmolar,           ///< Mole-based entropy
    iCpmolar,          ///< Mole-based constant-pressure specific heat
    iCp0molar,         ///< Mole-based ideal-gas constant-pressure specific heat
    iCvmolar,          ///< Mole-based constant-volume specific heat
    iUmolar,           ///< Mole-based internal energy
    iGmolar,           ///< Mole-based Gibbs energy
    iHelmholtzmolar,   ///< Mole-based Helmholtz energy
    iHmolar_residual,  ///< The residual molar enthalpy
    iSmolar_residual,  ///< The residual molar entropy (as a function of temperature and density)
    iGmolar_residual,  ///< The residual molar Gibbs energy

    // Mass specific thermodynamic properties
    iDmass,          ///< Mass-based density
    iHmass,          ///< Mass-based enthalpy
    iSmass,          ///< Mass-based entropy
    iCpmass,         ///< Mass-based constant-pressure specific heat
    iCp0mass,        ///< Mass-based ideal-gas specific heat
    iCvmass,         ///< Mass-based constant-volume specific heat
    iUmass,          ///< Mass-based internal energy
    iGmass,          ///< Mass-based Gibbs energy
    iHelmholtzmass,  ///< Mass-based Helmholtz energy

    // Transport properties
    iviscosity,        ///< Viscosity
    iconductivity,     ///< Thermal conductivity
    isurface_tension,  ///< Surface tension
    iPrandtl,          ///< The Prandtl number

    // Derivative-based terms
    ispeed_sound,                       ///< Speed of sound
    iisothermal_compressibility,        ///< Isothermal compressibility
    iisobaric_expansion_coefficient,    ///< Isobaric expansion coefficient
    iisentropic_expansion_coefficient,  ///< Isentropic expansion coefficient

    // Fundamental derivative of gas dynamics
    ifundamental_derivative_of_gas_dynamics,  ///< The fundamental derivative of gas dynamics

    // Derivatives of the residual non-dimensionalized Helmholtz energy with respect to the EOS variables
    ialphar,
    idalphar_dtau_constdelta,
    idalphar_ddelta_consttau,

    // Derivatives of the ideal-gas non-dimensionalized Helmholtz energy with respect to the EOS variables
    ialpha0,
    idalpha0_dtau_constdelta,
    idalpha0_ddelta_consttau,
    id2alpha0_ddelta2_consttau,
    id3alpha0_ddelta3_consttau,

    // Other functions and derivatives
    iBvirial,      ///< Second virial coefficient
    iCvirial,      ///< Third virial coefficient
    idBvirial_dT,  ///< Derivative of second virial coefficient with temperature
    idCvirial_dT,  ///< Derivative of third virial coefficient with temperature
    iZ,            ///< The compressibility factor Z = p*v/(R*T)
    iPIP,          ///< The phase identification parameter of Venkatarathnam and Oellrich

    // Accessors for incompressibles
    ifraction_min,  ///< The minimum fraction (mole, mass, volume) for incompressibles
    ifraction_max,  ///< The maximum fraction (mole,mass,volume) for incompressibles
    iT_freeze,      ///< The freezing temperature for incompressibles

    // Environmental parameters
    iGWP20,               ///< The 20-year global warming potential
    iGWP100,              ///< The 100-year global warming potential
    iGWP500,              ///< The 500-year global warming potential
    iFH,                  ///< Fire hazard index
    iHH,                  ///< Health hazard index
    iPH,                  ///< Physical hazard index
    iODP,                 ///< Ozone depletion potential (R-11 = 1.0)
    iPhase,               ///< The phase index of the given state
    iundefined_parameter  ///< The last parameter, so we can check that all parameters are described in DataStructures.cpp

};
// !! If you add a parameter, update the map in the corresponding CPP file !!
// !! Also update phase_lookup_string() in CoolProp.cpp                    !!

/// These are constants for the phases of the fluid
enum phases
{
    iphase_liquid,                ///< Subcritical liquid
    iphase_supercritical,         ///< Supercritical (p > pc, T > Tc)
    iphase_supercritical_gas,     ///< Supercritical gas (p < pc, T > Tc)
    iphase_supercritical_liquid,  ///< Supercritical liquid (p > pc, T < Tc)
    iphase_critical_point,        ///< At the critical point
    iphase_gas,                   ///< Subcritical gas
    iphase_twophase,              ///< Twophase
    iphase_unknown,               ///< Unknown phase
    iphase_not_imposed
};  ///< Phase is not imposed

/// Constants for the different PC-SAFT association schemes (see Huang and Radosz 1990)
enum schemes
{
    i1,
    i2a,
    i2b,
    i3a,
    i3b,
    i4a,
    i4b,
    i4c
};

/// Return information about the parameter
/// @param key The key, one of iT, iP, etc.
/// @param info The thing you want, one of "IO" ("IO" if input/output, "O" if output only), "short" (very short description), "long" (a longer description), "units"
std::string get_parameter_information(int key, const std::string& info);

/// Return the enum key corresponding to the parameter name ("Dmolar" for instance)
parameters get_parameter_index(const std::string& param_name);

/// Return true if passed phase name is valid, otherwise false
/// @param phase_name The phase name string to be checked ("phase_liquid" for instance)
/// @param iOutput Gets updated with the phases enum value if phase_name is found
bool is_valid_phase(const std::string& phase_name, phases& iOutput);

/// Return the enum key corresponding to the phase name ("phase_liquid" for instance)
phases get_phase_index(const std::string& param_name);

/// Return true if passed PC-SAFT association scheme name is valid, otherwise false
/// @param scheme_name The association scheme string to be checked ("2B" for instance)
/// @param iOutput Gets updated with the schemes enum value if scheme_name is found
bool is_valid_scheme(const std::string &scheme_name, schemes &iOutput);

/// Return the enum key corresponding to the association scheme name ("2B" for instance)
schemes get_scheme_index(const std::string &scheme_name);

/// Returns true if the input is trivial (constants, critical parameters, etc.)
bool is_trivial_parameter(int key);

/// Returns true if a valid parameter, and sets value in the variable iOutput
bool is_valid_parameter(const std::string& name, parameters& iOutput);

/// Returns true if the string corresponds to a valid first derivative
///
/// If it is a value derivative, the variables are set to the parts of the derivative
bool is_valid_first_derivative(const std::string& name, parameters& iOf, parameters& iWrt, parameters& iConstant);

/// Returns true if the string corresponds to a valid first saturation derivative - e.g. "d(P)/d(T)|sigma" for instance
///
/// If it is a valid derivative, the variables are set to the parts of the derivative
bool is_valid_first_saturation_derivative(const std::string& name, parameters& iOf, parameters& iWrt);

/// Returns true if the string corresponds to a valid second derivative
///
/// If it is a value derivative, the variables are set to the parts of the derivative
bool is_valid_second_derivative(const std::string& name, parameters& iOf1, parameters& iWrt1, parameters& iConstant1, parameters& iWrt2,
                                parameters& iConstant2);

/// Get a comma separated list of parameters
std::string get_csv_parameter_list();

/// These are constants for the compositions
enum composition_types
{
    IFRAC_MASS,
    IFRAC_MOLE,
    IFRAC_VOLUME,
    IFRAC_UNDEFINED,
    IFRAC_PURE
};

/// These are unit types for the fluid
enum fluid_types
{
    FLUID_TYPE_PURE,
    FLUID_TYPE_PSEUDOPURE,
    FLUID_TYPE_REFPROP,
    FLUID_TYPE_INCOMPRESSIBLE_LIQUID,
    FLUID_TYPE_INCOMPRESSIBLE_SOLUTION,
    FLUID_TYPE_UNDEFINED
};

// !! If you add a parameter, update the map in the corresponding CPP file !!
/// These are input pairs that can be used for the update function (in each pair, input keys are sorted alphabetically)
enum input_pairs
{
    INPUT_PAIR_INVALID = 0,  // Default (invalid) value
    QT_INPUTS,               ///< Molar quality, Temperature in K
    PQ_INPUTS,               ///< Pressure in Pa, Molar quality
    QSmolar_INPUTS,          ///< Molar quality, Entropy in J/mol/K
    QSmass_INPUTS,           ///< Molar quality, Entropy in J/kg/K
    HmolarQ_INPUTS,          ///< Enthalpy in J/mol, Molar quality
    HmassQ_INPUTS,           ///< Enthalpy in J/kg, Molar quality
    DmolarQ_INPUTS,          ///< Density in mol/m^3, Molar quality
    DmassQ_INPUTS,           ///< Density in kg/m^3, Molar quality

    PT_INPUTS,  ///< Pressure in Pa, Temperature in K

    DmassT_INPUTS,   ///< Mass density in kg/m^3, Temperature in K
    DmolarT_INPUTS,  ///< Molar density in mol/m^3, Temperature in K
    HmolarT_INPUTS,  ///< Enthalpy in J/mol, Temperature in K
    HmassT_INPUTS,   ///< Enthalpy in J/kg, Temperature in K
    SmolarT_INPUTS,  ///< Entropy in J/mol/K, Temperature in K
    SmassT_INPUTS,   ///< Entropy in J/kg/K, Temperature in K
    TUmolar_INPUTS,  ///< Temperature in K, Internal energy in J/mol
    TUmass_INPUTS,   ///< Temperature in K, Internal energy in J/kg

    DmassP_INPUTS,   ///< Mass density in kg/m^3, Pressure in Pa
    DmolarP_INPUTS,  ///< Molar density in mol/m^3, Pressure in Pa
    HmassP_INPUTS,   ///< Enthalpy in J/kg, Pressure in Pa
    HmolarP_INPUTS,  ///< Enthalpy in J/mol, Pressure in Pa
    PSmass_INPUTS,   ///< Pressure in Pa, Entropy in J/kg/K
    PSmolar_INPUTS,  ///< Pressure in Pa, Entropy in J/mol/K
    PUmass_INPUTS,   ///< Pressure in Pa, Internal energy in J/kg
    PUmolar_INPUTS,  ///< Pressure in Pa, Internal energy in J/mol

    HmassSmass_INPUTS,    ///< Enthalpy in J/kg, Entropy in J/kg/K
    HmolarSmolar_INPUTS,  ///< Enthalpy in J/mol, Entropy in J/mol/K
    SmassUmass_INPUTS,    ///< Entropy in J/kg/K, Internal energy in J/kg
    SmolarUmolar_INPUTS,  ///< Entropy in J/mol/K, Internal energy in J/mol

    DmassHmass_INPUTS,    ///< Mass density in kg/m^3, Enthalpy in J/kg
    DmolarHmolar_INPUTS,  ///< Molar density in mol/m^3, Enthalpy in J/mol
    DmassSmass_INPUTS,    ///< Mass density in kg/m^3, Entropy in J/kg/K
    DmolarSmolar_INPUTS,  ///< Molar density in mol/m^3, Entropy in J/mol/K
    DmassUmass_INPUTS,    ///< Mass density in kg/m^3, Internal energy in J/kg
    DmolarUmolar_INPUTS,  ///< Molar density in mol/m^3, Internal energy in J/mol
};
// !! If you add or remove a parameter, update the map in the corresponding CPP file !!

inline bool match_pair(parameters key1, parameters key2, parameters x1, parameters x2, bool& swap) {
    swap = !(key1 == x1);
    return ((key1 == x1 && key2 == x2) || (key2 == x1 && key1 == x2));
};
/**
 * @brief Generate an update pair from key, value pairs
 *
 * If the input pair is valid, v1 and v2 will correspond to the returned output pair
 *
 * @param key1 The first input key
 * @param value1 The first input value
 * @param key2 The second input key
 * @param value2 The second input value
 * @param out1 The first output value
 * @param out2 The second output value
 * @return pair, or INPUT_PAIR_INVALID if not valid
 */
template <class T>
CoolProp::input_pairs generate_update_pair(parameters key1, T value1, parameters key2, T value2, T& out1, T& out2) throw() {
    CoolProp::input_pairs pair;
    bool swap;

    if (match_pair(key1, key2, iQ, iT, swap)) {
        pair = QT_INPUTS;  ///< Molar quality, Temperature in K
    } else if (match_pair(key1, key2, iP, iQ, swap)) {
        pair = PQ_INPUTS;  ///< Pressure in Pa, Molar quality
    } else if (match_pair(key1, key2, iP, iT, swap)) {
        pair = PT_INPUTS;  ///< Pressure in Pa, Temperature in K
    } else if (match_pair(key1, key2, iDmolar, iT, swap)) {
        pair = DmolarT_INPUTS;  // Molar density in mol/m^3, Temperature in K
    } else if (match_pair(key1, key2, iDmass, iT, swap)) {
        pair = DmassT_INPUTS;  // Mass density in kg/m^3, Temperature in K
    } else if (match_pair(key1, key2, iHmolar, iT, swap)) {
        pair = HmolarT_INPUTS;  // Enthalpy in J/mol, Temperature in K
    } else if (match_pair(key1, key2, iHmass, iT, swap)) {
        pair = HmassT_INPUTS;  // Enthalpy in J/kg, Temperature in K
    } else if (match_pair(key1, key2, iSmolar, iT, swap)) {
        pair = SmolarT_INPUTS;  // Entropy in J/mol/K, Temperature in K
    } else if (match_pair(key1, key2, iSmass, iT, swap)) {
        pair = SmassT_INPUTS;  // Entropy in J/kg/K, Temperature in K
    } else if (match_pair(key1, key2, iT, iUmolar, swap)) {
        pair = TUmolar_INPUTS;  // Temperature in K, Internal energy in J/mol
    } else if (match_pair(key1, key2, iT, iUmass, swap)) {
        pair = TUmass_INPUTS;  // Temperature in K, Internal energy in J/kg
    } else if (match_pair(key1, key2, iDmass, iHmass, swap)) {
        pair = DmassHmass_INPUTS;  // Mass density in kg/m^3, Enthalpy in J/kg
    } else if (match_pair(key1, key2, iDmolar, iHmolar, swap)) {
        pair = DmolarHmolar_INPUTS;  // Molar density in mol/m^3, Enthalpy in J/mol
    } else if (match_pair(key1, key2, iDmass, iSmass, swap)) {
        pair = DmassSmass_INPUTS;  // Mass density in kg/m^3, Entropy in J/kg/K
    } else if (match_pair(key1, key2, iDmolar, iSmolar, swap)) {
        pair = DmolarSmolar_INPUTS;  // Molar density in mol/m^3, Entropy in J/mol/K
    } else if (match_pair(key1, key2, iDmass, iUmass, swap)) {
        pair = DmassUmass_INPUTS;  // Mass density in kg/m^3, Internal energy in J/kg
    } else if (match_pair(key1, key2, iDmolar, iUmolar, swap)) {
        pair = DmolarUmolar_INPUTS;  // Molar density in mol/m^3, Internal energy in J/mol
    } else if (match_pair(key1, key2, iDmass, iP, swap)) {
        pair = DmassP_INPUTS;  // Mass density in kg/m^3, Pressure in Pa
    } else if (match_pair(key1, key2, iDmolar, iP, swap)) {
        pair = DmolarP_INPUTS;  // Molar density in mol/m^3, Pressure in Pa
    } else if (match_pair(key1, key2, iDmass, iQ, swap)) {
        pair = DmassQ_INPUTS;  // Mass density in kg/m^3, molar vapor quality
    } else if (match_pair(key1, key2, iDmolar, iQ, swap)) {
        pair = DmolarQ_INPUTS;  // Molar density in mol/m^3, molar vapor quality
    } else if (match_pair(key1, key2, iHmass, iP, swap)) {
        pair = HmassP_INPUTS;  // Enthalpy in J/kg, Pressure in Pa
    } else if (match_pair(key1, key2, iHmolar, iP, swap)) {
        pair = HmolarP_INPUTS;  // Enthalpy in J/mol, Pressure in Pa
    } else if (match_pair(key1, key2, iP, iSmass, swap)) {
        pair = PSmass_INPUTS;  // Pressure in Pa, Entropy in J/kg/K
    } else if (match_pair(key1, key2, iP, iSmolar, swap)) {
        pair = PSmolar_INPUTS;  // Pressure in Pa, Entropy in J/mol/K
    } else if (match_pair(key1, key2, iP, iUmass, swap)) {
        pair = PUmass_INPUTS;  // Pressure in Pa, Internal energy in J/kg
    } else if (match_pair(key1, key2, iP, iUmolar, swap)) {
        pair = PUmolar_INPUTS;  // Pressure in Pa, Internal energy in J/mol
    } else if (match_pair(key1, key2, iHmass, iSmass, swap)) {
        pair = HmassSmass_INPUTS;  // Enthalpy in J/kg, Entropy in J/kg/K
    } else if (match_pair(key1, key2, iHmolar, iSmolar, swap)) {
        pair = HmolarSmolar_INPUTS;  // Enthalpy in J/mol, Entropy in J/mol/K
    } else if (match_pair(key1, key2, iSmass, iUmass, swap)) {
        pair = SmassUmass_INPUTS;  ///< Entropy in J/kg/K, Internal energy in J/kg
    } else if (match_pair(key1, key2, iSmolar, iUmolar, swap)) {
        pair = SmolarUmolar_INPUTS;  ///< Entropy in J/mol/K, Internal energy in J/mol
    } else {
        pair = INPUT_PAIR_INVALID;
        return pair;
    }

    if (!swap) {
        out1 = value1;
        out2 = value2;
    } else {
        out1 = value2;
        out2 = value1;
    }
    return pair;
};

/// Get the input pair index associated with its string representation
input_pairs get_input_pair_index(const std::string& input_pair_name);

/// Return the short description of an input pair key ("DmolarT_INPUTS" for instance)
const std::string& get_input_pair_short_desc(input_pairs pair);

/// Return the long description of an input pair key ("Molar density in mol/m^3, Temperature in K" for instance)
const std::string& get_input_pair_long_desc(input_pairs pair);

/// Split an input pair into parameters for the two parts that form the pair
void split_input_pair(input_pairs pair, parameters& p1, parameters& p2);

extern std::string get_mixture_binary_pair_data(const std::string& CAS1, const std::string& CAS2, const std::string& param);
extern void set_mixture_binary_pair_data(const std::string& CAS1, const std::string& CAS2, const std::string& param, const double val);
extern std::string get_mixture_binary_pair_pcsaft(const std::string& CAS1, const std::string& CAS2, const std::string& param);
extern void set_mixture_binary_pair_pcsaft(const std::string& CAS1, const std::string& CAS2, const std::string& param, const double val);

/// The structure is taken directly from the AbstractState class.
// !! If you add a parameter, update the map in the corresponding CPP file !!
enum backend_families
{
    INVALID_BACKEND_FAMILY = 0,
    HEOS_BACKEND_FAMILY,
    REFPROP_BACKEND_FAMILY,
    INCOMP_BACKEND_FAMILY,
    IF97_BACKEND_FAMILY,
    TREND_BACKEND_FAMILY,
    TTSE_BACKEND_FAMILY,
    BICUBIC_BACKEND_FAMILY,
    SRK_BACKEND_FAMILY,
    PR_BACKEND_FAMILY,
    VTPR_BACKEND_FAMILY,
    PCSAFT_BACKEND_FAMILY
};
enum backends
{
    INVALID_BACKEND = 0,
    HEOS_BACKEND_PURE,
    HEOS_BACKEND_MIX,
    REFPROP_BACKEND_PURE,
    REFPROP_BACKEND_MIX,
    INCOMP_BACKEND,
    IF97_BACKEND,
    TREND_BACKEND,
    TTSE_BACKEND,
    BICUBIC_BACKEND,
    SRK_BACKEND,
    PR_BACKEND,
    VTPR_BACKEND,
    PCSAFT_BACKEND
};

/// Convert a string into the enum values
void extract_backend_families(std::string backend_string, backend_families& f1, backend_families& f2);
void extract_backend_families_string(std::string backend_string, backend_families& f1, std::string& f2);
std::string get_backend_string(backends backend);

#if !defined(NO_FMTLIB) && FMT_VERSION >= 90000
/// Allows enums to be formatted
inline int format_as(parameters parameter) {
    return fmt::underlying(parameter);
}

inline int format_as(phases phase) {
    return fmt::underlying(phase);
}

inline int format_as(schemes scheme) {
    return fmt::underlying(scheme);
}

inline int format_as(composition_types type) {
    return fmt::underlying(type);
}

inline int format_as(fluid_types type) {
    return fmt::underlying(type);
}

inline int format_as(input_pairs pair) {
    return fmt::underlying(pair);
}

inline int format_as(backend_families family) {
    return fmt::underlying(family);
}

inline int format_as(backends backend) {
    return fmt::underlying(backend);
}
#endif

} /* namespace CoolProp */
#endif /* DATASTRUCTURES_H_ */
