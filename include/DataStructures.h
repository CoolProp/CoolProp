/*
 * DataStructures.h
 *
 *  Created on: 21 Dec 2013
 *      Author: jowr
 */

#ifndef DATASTRUCTURES_H_
#define DATASTRUCTURES_H_

#include "Exceptions.h"
#include <map>
namespace CoolProp {

struct SimpleState
{
    double rhomolar, T, p, hmolar, smolar, umolar;
};

/// --------------------------------------------------
/// Define some constants that will be used throughout
/// --------------------------------------------------
/// These are constants for the input and output parameters
/// The structure is taken directly from the AbstractState class.
//
// !! If you add a parameter, update the map in the corresponding CPP file !!
enum parameters{

    // General parameters
    imolar_mass, irhomolar_reducing, irhomolar_critical, iT_reducing, iT_critical,

    // Bulk properties
    iT,  iP, iQ, iTau, iDelta,

    // Molar specific thermodynamic properties
    iDmolar, iHmolar, iSmolar, iCpmolar, iCp0molar, iCvmolar, iUmolar, iGmolar, 

    // Mass specific thermodynamic properties
    iDmass, iHmass, iSmass, iCpmass, iCp0mass, iCvmass, iUmass, iGmass,

    // Smoothing functions for density
    idrhodh_constp_smoothed, idrhodp_consth_smoothed, irho_smoothed,

    // Transport properties
    iviscosity, iconductivity, isurface_tension, iPrandtl,

    // Derivative-based terms
    ispeed_sound, iisothermal_compressibility, iisobaric_expansion_coefficient,
    
    // Fundamental derivative of gas dynamics
    ifundamental_derivative_of_gas_dynamics, id2pdv2_consts,

    // Derivatives of the residual non-dimensionalized Helmholtz energy with respect to the EOS variables
    ialphar, idalphar_dtau_constdelta, idalphar_ddelta_consttau, 
    
    // Derivatives of the ideal-gas non-dimensionalized Helmholtz energy with respect to the EOS variables
    ialpha0, idalpha0_dtau_constdelta, idalpha0_ddelta_consttau,

    // Other functions and derivatives
    iBvirial, iCvirial, iZ, idBvirial_dT, idCvirial_dT, idZdDelta, idZdTau,

    // Environmental parameters
    iGWP20, iGWP100, iGWP500, iFH, iHH, iPH, iODP
};
// !! If you add a parameter, update the map in the corresponding CPP file !!

/// Return information about the parameter
/// @param key The key, one of iT, iP, etc.
/// @param info The thing you want, one of "IO" ("IO" if input/output, "O" if output only), "short" (very short description), "long" (a longer description), "units"
std::string get_parameter_information(int key, std::string info);
/// Return the integer key corresponding to the parameter name ("Dmolar" for instance)
int get_parameter_index(const std::string &param_name);

/// These are constants for the phases of the fluid
enum phases {iphase_liquid, iphase_supercritical, iphase_gas, iphase_twophase, iphase_unknown};

/// These are unit types for the fluid
enum fluid_types{FLUID_TYPE_PURE, FLUID_TYPE_PSEUDOPURE, FLUID_TYPE_REFPROP, FLUID_TYPE_INCOMPRESSIBLE_LIQUID, FLUID_TYPE_INCOMPRESSIBLE_SOLUTION, FLUID_TYPE_UNDEFINED};

// !! If you add a parameter, update the map in the corresponding CPP file !!
/// These are input pairs that can be used (in each pair, input keys are sorted alphabetically)
enum input_pairs{
    QT_INPUTS, ///< Molar quality, Temperature in K
    PQ_INPUTS, ///< Pressure in Pa, Molar quality
                 
    PT_INPUTS, ///< Pressure in Pa, Temperature in K

    DmassT_INPUTS, ///< Mass density in kg/m^3, Temperature in K
    DmolarT_INPUTS, ///< Molar density in mol/m^3, Temperature in K
    HmolarT_INPUTS, ///< Enthalpy in J/mol, Temperature in K
    HmassT_INPUTS, ///< Enthalpy in J/kg, Temperature in K
    SmolarT_INPUTS, ///< Entropy in J/mol/K, Temperature in K
    SmassT_INPUTS, ///< Entropy in J/kg/K, Temperature in K
    TUmolar_INPUTS, ///< Temperature in K, Internal energy in J/mol
    TUmass_INPUTS, ///< Temperature in K, Internal energy in J/kg
                 
    DmassP_INPUTS, ///< Mass density in kg/m^3, Pressure in Pa
    DmolarP_INPUTS, ///< Molar density in mol/m^3, Pressure in Pa
    HmassP_INPUTS, ///< Enthalpy in J/kg, Pressure in Pa
    HmolarP_INPUTS, ///< Enthalpy in J/mol, Pressure in Pa
    PSmass_INPUTS, ///< Pressure in Pa, Entropy in J/kg/K
    PSmolar_INPUTS, ///< Pressure in Pa, Entropy in J/mol/K 
    PUmass_INPUTS, ///< Pressure in Pa, Internal energy in J/kg
    PUmolar_INPUTS, ///< Pressure in Pa, Internal energy in J/mol
                 
    HmassSmass_INPUTS, ///< Enthalpy in J/kg, Entropy in J/kg/K
    HmolarSmolar_INPUTS, ///< Enthalpy in J/mol, Entropy in J/mol/K
    SmassUmass_INPUTS, ///< Entropy in J/kg/K, Internal energy in J/kg
    SmolarUmolar_INPUTS, ///< Entropy in J/mol/K, Internal energy in J/mol
                 
    DmassHmass_INPUTS, ///< Mass density in kg/m^3, Enthalpy in J/kg
    DmolarHmolar_INPUTS, ///< Molar density in mol/m^3, Enthalpy in J/mol
    DmassSmass_INPUTS, ///< Mass density in kg/m^3, Entropy in J/kg/K
    DmolarSmolar_INPUTS, ///< Molar density in mol/m^3, Entropy in J/mol/K
    DmassUmass_INPUTS, ///< Mass density in kg/m^3, Internal energy in J/kg
    DmolarUmolar_INPUTS, ///< Molar density in mol/m^3, Internal energy in J/mol
};
// !! If you add or remove a parameter, update the map in the corresponding CPP file !!

inline bool match_pair(long key1, long key2, long x1, long x2, bool &swap)
{
    swap = !(key1 == x1);
    return ((key1 == x1 && key2 == x2) || (key2 == x1 && key1 == x2));
};
template<class T> long generate_update_pair(long key1, T value1, long key2, T value2, T &out1, T&out2)
    {
        long pair;
        bool swap;

        if (match_pair(key1, key2, iQ, iT, swap)){
            pair = QT_INPUTS; ///< Molar quality, Temperature in K
        }
        else if (match_pair(key1, key2, iP, iQ, swap)){
            pair = PQ_INPUTS; ///< Pressure in Pa, Molar quality
        }
        else if (match_pair(key1, key2, iP, iT, swap)){
            pair = PT_INPUTS; ///< Pressure in Pa, Temperature in K
        }
        else if (match_pair(key1, key2, iDmolar, iT, swap)){
            pair = DmolarT_INPUTS; // Molar density in mol/m^3, Temperature in K
        }
        else if (match_pair(key1, key2, iDmass, iT, swap)){
            pair = DmassT_INPUTS; // Mass density in kg/m^3, Temperature in K
        }
        else if (match_pair(key1, key2, iHmolar, iT, swap)){
            pair = HmolarT_INPUTS; // Enthalpy in J/mol, Temperature in K
        }
        else if (match_pair(key1, key2, iHmass, iT, swap)){
            pair = HmassT_INPUTS; // Enthalpy in J/kg, Temperature in K
        }
        else if (match_pair(key1, key2, iSmolar, iT, swap)){
            pair = SmolarT_INPUTS; // Entropy in J/mol/K, Temperature in K
        }
        else if (match_pair(key1, key2, iSmass, iT, swap)){
            pair = SmassT_INPUTS; // Entropy in J/kg/K, Temperature in K
        }
        else if (match_pair(key1, key2, iT, iUmolar, swap)){
            pair = TUmolar_INPUTS; // Temperature in K, Internal energy in J/mol
        }
        else if (match_pair(key1, key2, iT, iUmass, swap)){
            pair = TUmass_INPUTS; // Temperature in K, Internal energy in J/kg
        }
        else if (match_pair(key1, key2, iDmass, iHmass, swap)){
            pair = DmassHmass_INPUTS; // Mass density in kg/m^3, Enthalpy in J/kg
        }
        else if (match_pair(key1, key2, iDmolar, iHmolar, swap)){
            pair = DmolarHmolar_INPUTS; // Molar density in mol/m^3, Enthalpy in J/mol
        }
        else if (match_pair(key1, key2, iDmass, iSmass, swap)){
            pair = DmassSmass_INPUTS; // Mass density in kg/m^3, Entropy in J/kg/K
        }
        else if (match_pair(key1, key2, iDmolar, iSmolar, swap)){
            pair = DmolarSmolar_INPUTS; // Molar density in mol/m^3, Entropy in J/mol/K
        }
        else if (match_pair(key1, key2, iDmass, iUmass, swap)){
            pair = DmassUmass_INPUTS; // Mass density in kg/m^3, Internal energy in J/kg
        }
        else if (match_pair(key1, key2, iDmolar, iUmolar, swap)){
            pair = DmolarUmolar_INPUTS; // Molar density in mol/m^3, Internal energy in J/mol
        }
        else if (match_pair(key1, key2, iDmass, iP, swap)){
            pair = DmassP_INPUTS; // Mass density in kg/m^3, Pressure in Pa
        }
        else if (match_pair(key1, key2, iDmolar, iP, swap)){
            pair = DmolarP_INPUTS; // Molar density in mol/m^3, Pressure in Pa
        }
        else if (match_pair(key1, key2, iHmass, iP, swap)){
            pair = HmassP_INPUTS; // Enthalpy in J/kg, Pressure in Pa
        }
        else if (match_pair(key1, key2, iHmolar, iP, swap)){
            pair = HmolarP_INPUTS; // Enthalpy in J/mol, Pressure in Pa
        }
        else if (match_pair(key1, key2, iP, iSmass, swap)){
            pair = PSmass_INPUTS; // Pressure in Pa, Entropy in J/kg/K
        }
        else if (match_pair(key1, key2, iP, iSmolar, swap)){
            pair = PSmolar_INPUTS; // Pressure in Pa, Entropy in J/mol/K 
        }
        else if (match_pair(key1, key2, iP, iUmass, swap)){
            pair = PUmass_INPUTS; // Pressure in Pa, Internal energy in J/kg
        }
        else if (match_pair(key1, key2, iP, iUmolar, swap)){
            pair = PUmolar_INPUTS; // Pressure in Pa, Internal energy in J/mol
        }
        else
            throw ValueError("Invalid set of inputs to generate_update_pair");

        /*
        HmassSmass_INPUTS, ///< Enthalpy in J/kg, Entropy in J/kg/K
        HmolarSmolar_INPUTS, ///< Enthalpy in J/mol, Entropy in J/mol/K
        SmassUmass_INPUTS, ///< Entropy in J/kg/K, Internal energy in J/kg
        SmolarUmolar_INPUTS, ///< Entropy in J/mol/K, Internal energy in J/mol   
        */

        if (!swap)
        {
            out1 = value1;
            out2 = value2;
        }
        else
        {
            out1 = value2;
            out2 = value1;
        }
        return pair;
    };

/// Return the short description of an input pair key ("DmolarT_INPUTS" for instance)
std::string get_input_pair_short_desc(int pair);
/// Return the long description of an input pair key ("Molar density in mol/m^3, Temperature in K" for instance)
std::string get_input_pair_long_desc(int pair);



} /* namespace CoolProp */
#endif /* DATASTRUCTURES_H_ */
