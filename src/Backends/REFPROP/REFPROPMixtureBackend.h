/*
 * AbstractBackend.h
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#ifndef REFPROPMIXTUREBACKEND_H_
#define REFPROPMIXTUREBACKEND_H_

#include "AbstractState.h"

#include <vector>

namespace CoolProp {


class REFPROPMixtureBackend : public AbstractState  {
protected:
    std::size_t Ncomp;
    bool _mole_fractions_set;
    static bool _REFPROP_supported;
    std::vector<double> mole_fractions, mass_fractions;
    std::vector<double> mole_fractions_liq, mole_fractions_vap;
public:
    REFPROPMixtureBackend(){};

    /// The instantiator
    /// @param fluid_names The vector of strings of the fluid components, without file ending
    REFPROPMixtureBackend(const std::vector<std::string>& fluid_names);
    virtual ~REFPROPMixtureBackend();

    // REFPROP backend uses mole fractions
    bool using_mole_fractions(){return true;}
    bool using_mass_fractions(){return false;}
    bool using_volu_fractions(){return false;}

    /// Updating function for REFPROP
    /**
    In this function we take a pair of thermodynamic states, those defined in the input_pairs
    enumeration and update all the internal variables that we can.  REFPROP calculates
    a lot of other state variables each time you use a flash routine so we cache all the
    outputs that we can, which saves on computational time.

    @param input_pair Integer key from CoolProp::input_pairs to the two inputs that will be passed to the function
    @param value1 First input value
    @param value2 Second input value
    */
    void update(CoolProp::input_pairs,
                double value1,
                double value2
                );

    long double calc_molar_mass(void);

    /// Returns true if REFPROP is supported on this platform
    bool REFPROP_supported(void);

    long double calc_cpmolar_idealgas(void);

    long double calc_first_partial_deriv(parameters Of, parameters Wrt, parameters Constant);

    /// Set the fluids in REFPROP DLL by calling the SETUPdll function
    /**
    @param fluid_names The vector of strings of the fluid components, without file ending
    */
    void set_REFPROP_fluids(const std::vector<std::string> &fluid_names);

    /// Set the mole fractions
    /**
    @param mole_fractions The vector of mole fractions of the components
    */
    void set_mole_fractions(const std::vector<long double> &mole_fractions);

    /// Set the mass fractions
    /**
    @param mass_fractions The vector of mass fractions of the components
    */
    void set_mass_fractions(const std::vector<long double> &mass_fractions);

    void calc_phase_envelope(const std::string &type);

    std::vector<long double> calc_mole_fractions_liquid(void){return std::vector<long double>(mole_fractions_liq.begin(), mole_fractions_liq.end());}
    std::vector<long double> calc_mole_fractions_vapor(void){return std::vector<long double>(mole_fractions_vap.begin(), mole_fractions_vap.end());}

    /// Check if the mole fractions have been set, etc.
    void check_status();

    /// Get the viscosity [Pa-s] (based on the temperature and density in the state class)
    long double calc_viscosity(void);
    /// Get the thermal conductivity [W/m/K] (based on the temperature and density in the state class)
    long double calc_conductivity(void);
    /// Get the surface tension [N/m] (based on the temperature in the state class).  Invalid for temperatures above critical point or below triple point temperature
    long double calc_surface_tension(void);

    long double calc_fugacity_coefficient(int i);
    long double calc_melting_line(int param, int given, long double value);
    bool has_melting_line(){return true;};
    double calc_melt_Tmax();
    long double calc_T_critical(void);
    long double calc_p_critical(void);
    long double calc_rhomolar_critical(void);
    long double calc_Ttriple(void);

    /// A wrapper function to calculate the limits for the EOS
    void limits(double &Tmin, double &Tmax, double &rhomolarmax, double &pmax);
    /// Calculate the maximum pressure
    long double calc_pmax(void);
    /// Calculate the maximum temperature
    long double calc_Tmax(void);
};

} /* namespace CoolProp */
#endif /* REFPROPMIXTUREBACKEND_H_ */
