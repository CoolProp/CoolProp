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
    private:
    std::string cached_component_string;
protected:
    std::size_t Ncomp;
    bool _mole_fractions_set;
    
    static std::size_t instance_counter;
    static bool _REFPROP_supported;
    std::vector<CoolPropDbl> mole_fractions_long_double; // read-only
    std::vector<double> mole_fractions, mass_fractions;
    std::vector<double> mole_fractions_liq, mole_fractions_vap;
	std::vector<std::string> fluid_names;

	
	/// Call the PHIXdll function in the dll
	CoolPropDbl call_phixdll(long itau, long idelta);
	/// Call the PHI0dll function in the dll
	CoolPropDbl call_phi0dll(long itau, long idelta);
	
public:
    REFPROPMixtureBackend(){cached_component_string = ""; instance_counter++;};

    /// The instantiator
    /// @param fluid_names The vector of strings of the fluid components, without file ending
    REFPROPMixtureBackend(const std::vector<std::string>& fluid_names);
    virtual ~REFPROPMixtureBackend();
	
	std::vector<std::string> calc_fluid_names(){return fluid_names;};

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

    CoolPropDbl calc_molar_mass(void);
    
    void check_loaded_fluid(void);

    /// Returns true if REFPROP is supported on this platform
    bool REFPROP_supported(void);

    CoolPropDbl calc_cpmolar_idealgas(void);

    /// Set the fluids in REFPROP DLL by calling the SETUPdll function
    /**
    @param fluid_names The vector of strings of the fluid components, without file ending
    */
    void set_REFPROP_fluids(const std::vector<std::string> &fluid_names);

    /// Set the mole fractions
    /**
    @param mole_fractions The vector of mole fractions of the components
    */
    void set_mole_fractions(const std::vector<CoolPropDbl> &mole_fractions);

    /// Set the mass fractions
    /**
    @param mass_fractions The vector of mass fractions of the components
    */
    void set_mass_fractions(const std::vector<CoolPropDbl> &mass_fractions);
    
    const std::vector<CoolPropDbl> &get_mole_fractions(){return mole_fractions_long_double;};

    void calc_phase_envelope(const std::string &type);

    std::vector<CoolPropDbl> calc_mole_fractions_liquid(void){return std::vector<CoolPropDbl>(mole_fractions_liq.begin(), mole_fractions_liq.end());}
    std::vector<CoolPropDbl> calc_mole_fractions_vapor(void){return std::vector<CoolPropDbl>(mole_fractions_vap.begin(), mole_fractions_vap.end());}

    /// Check if the mole fractions have been set, etc.
    void check_status();

    /// Get the viscosity [Pa-s] (based on the temperature and density in the state class)
    CoolPropDbl calc_viscosity(void);
    /// Get the thermal conductivity [W/m/K] (based on the temperature and density in the state class)
    CoolPropDbl calc_conductivity(void);
    /// Get the surface tension [N/m] (based on the temperature in the state class).  Invalid for temperatures above critical point or below triple point temperature
    CoolPropDbl calc_surface_tension(void);

    CoolPropDbl calc_fugacity_coefficient(int i);
    CoolPropDbl calc_melting_line(int param, int given, CoolPropDbl value);
    bool has_melting_line(){return true;};
    double calc_melt_Tmax();
    CoolPropDbl calc_T_critical(void);
	CoolPropDbl calc_T_reducing(void);
    CoolPropDbl calc_p_critical(void);
    CoolPropDbl calc_rhomolar_critical(void);
	CoolPropDbl calc_rhomolar_reducing(void);
    CoolPropDbl calc_Ttriple(void);
	CoolPropDbl calc_gas_constant(void);

    /// A wrapper function to calculate the limits for the EOS
    void limits(double &Tmin, double &Tmax, double &rhomolarmax, double &pmax);
    /// Calculate the maximum pressure
    CoolPropDbl calc_pmax(void);
    /// Calculate the maximum temperature
    CoolPropDbl calc_Tmax(void);
	
	/// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r\f$ (dimensionless)
	CoolPropDbl calc_alphar(void){return call_phixdll(0,0);};
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta}\f$ (dimensionless)
	CoolPropDbl calc_dalphar_dDelta(void){ return call_phixdll(0,1); };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\tau}\f$ (dimensionless)
	CoolPropDbl calc_dalphar_dTau(void){ return call_phixdll(1,0); };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta}\f$ (dimensionless)
	CoolPropDbl calc_d2alphar_dDelta2(void){ return call_phixdll(0,2); };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\tau}\f$ (dimensionless)
	CoolPropDbl calc_d2alphar_dDelta_dTau(void){ return call_phixdll(1,1); };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\tau\tau}\f$ (dimensionless)
	CoolPropDbl calc_d2alphar_dTau2(void){ return call_phixdll(2,0); };
	/// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta\delta}\f$ (dimensionless)
    CoolPropDbl calc_d3alphar_dDelta3(void){ return call_phixdll(0,3); };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta\tau}\f$ (dimensionless)
	CoolPropDbl calc_d3alphar_dDelta2_dTau(void){ return call_phixdll(1,2); };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\tau\tau}\f$ (dimensionless)
	CoolPropDbl calc_d3alphar_dDelta_dTau2(void){ return call_phixdll(2,1); };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\tau\tau\tau}\f$ (dimensionless)
	CoolPropDbl calc_d3alphar_dTau3(void){ return call_phixdll(3,0); };

    CoolPropDbl calc_alpha0(void){ return call_phi0dll(0,0); };
    CoolPropDbl calc_dalpha0_dDelta(void){ return call_phi0dll(0,1); };
    CoolPropDbl calc_dalpha0_dTau(void){ return call_phi0dll(1,0); };
    CoolPropDbl calc_d2alpha0_dDelta2(void){ return call_phi0dll(0,2); };
    CoolPropDbl calc_d2alpha0_dDelta_dTau(void){ return call_phi0dll(1,1); };
    CoolPropDbl calc_d2alpha0_dTau2(void){ return call_phi0dll(2,0); };
    CoolPropDbl calc_d3alpha0_dDelta3(void){ return call_phi0dll(0,3); };
    CoolPropDbl calc_d3alpha0_dDelta2_dTau(void){ return call_phi0dll(1,2);	};
    CoolPropDbl calc_d3alpha0_dDelta_dTau2(void){ return call_phi0dll(2,1); };
    CoolPropDbl calc_d3alpha0_dTau3(void){ return call_phi0dll(3,0); };
};

} /* namespace CoolProp */
#endif /* REFPROPMIXTUREBACKEND_H_ */
