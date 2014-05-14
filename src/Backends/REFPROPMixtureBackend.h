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
    int Ncomp;
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
	void update(long input_pair,
				double value1,
				double value2
				);

	/// Returns true if REFPROP is supported on this platform
	bool REFPROP_supported(void);

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

	/// Check if the mole fractions have been set, etc.
	void check_status();

	/// Get the viscosity [Pa-s] (based on the temperature and density in the state class)
	long double calc_viscosity(void);
	/// Get the thermal conductivity [W/m/K] (based on the temperature and density in the state class)
	long double calc_conductivity(void);
	/// Get the surface tension [N/m] (based on the temperature in the state class).  Invalid for temperatures above critical point or below triple point temperature
	long double calc_surface_tension(void);

    long double calc_fugacity_coefficient(int i);
    
    long double calc_melt_p_T(long double T);
    long double calc_melt_T_p(long double p);
    long double calc_melt_rho_T(long double T);
    double calc_melt_Tmax();
};

} /* namespace CoolProp */
#endif /* REFPROPMIXTUREBACKEND_H_ */
