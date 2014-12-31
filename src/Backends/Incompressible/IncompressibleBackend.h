
#ifndef INCOMPRESSIBLEBACKEND_H_
#define INCOMPRESSIBLEBACKEND_H_

#include "DataStructures.h"
#include "IncompressibleFluid.h"
#include "AbstractState.h"
#include "Exceptions.h"

#include <vector>

namespace CoolProp {

class IncompressibleBackend : public AbstractState  {
protected:
    //int Ncomp;
    //static bool _REFPROP_supported;
    //int _fractions_id;
    std::vector<long double> _fractions;
    IncompressibleFluid *fluid;

    /// Set the fractions
    /**
    @param fractions The vector of fractions of the components converted to the correct input
    */
    void set_fractions(const std::vector<long double> &fractions);

public:
    IncompressibleBackend();
    virtual ~IncompressibleBackend(){};

    /// The instantiator
    /// @param fluid object, mostly for testing purposes
    IncompressibleBackend(IncompressibleFluid* fluid);
    /// The instantiator
    /// @param fluid_name the string with the fluid name
    IncompressibleBackend(const std::string &fluid_name);
    /// The instantiator
    /// @param component_names The vector of strings of the fluid components, without file ending
    IncompressibleBackend(const std::vector<std::string> &component_names);

    // Incompressible backend uses different compositions
    bool using_mole_fractions(void){return  this->fluid->getxid()==IFRAC_MOLE;};
    bool using_mass_fractions(void){return (this->fluid->getxid()==IFRAC_MASS || this->fluid->getxid()==IFRAC_PURE);};
    bool using_volu_fractions(void){return  this->fluid->getxid()==IFRAC_VOLUME;};

    /// Updating function for incompressible fluid
    /**
    In this function we take a pair of thermodynamic states, those defined in the input_pairs
    enumeration and update all the internal variables that we can.

    @param input_pair Integer key from CoolProp::input_pairs to the two inputs that will be passed to the function
    @param value1 First input value
    @param value2 Second input value
    */
    void update(CoolProp::input_pairs input_pair, double value1, double value2);

    /// Set the mole fractions
    /**
    @param mole_fractions The vector of mole fractions of the components
    */
    void set_mole_fractions(const std::vector<long double> &mole_fractions);
    const std::vector<long double> & get_mole_fractions(void){throw NotImplementedError("get_mole_fractions not implemented for this backend");};

    /// Set the mass fractions
    /**
    @param mass_fractions The vector of mass fractions of the components
    */
    void set_mass_fractions(const std::vector<long double> &mass_fractions);

    /// Set the volume fractions
    /**
    @param volu_fractions The vector of volume fractions of the components
    */
    void set_volu_fractions(const std::vector<long double> &volu_fractions);

    /// Check if the mole fractions have been set, etc.
    void check_status();

    /// Calculate T given pressure and density
    /**
    @param rhomass The mass density in kg/m^3
    @param p The pressure in Pa
    @returns T The temperature in K
    */
    long double DmassP_flash(long double rhomass, long double p);
    /// Calculate T given pressure and enthalpy
    /**
    @param hmass The mass enthalpy in J/kg
    @param p The pressure in Pa
    @returns T The temperature in K
    */
    long double HmassP_flash(long double hmass, long double p);
    /// Calculate T given pressure and entropy
    /**
    @param smass The mass entropy in J/kg/K
    @param p The pressure in Pa
    @returns T The temperature in K
    */
    long double PSmass_flash(long double p, long double smass);

    /// Calculate T given pressure and internal energy
    /**
    @param umass The mass internal energy in J/kg
    @param p The pressure in Pa
    @returns T The temperature in K
    */
    long double PUmass_flash(long double p, long double umass);

    /// Get the viscosity [Pa-s]
    long double calc_viscosity(void){return fluid->visc(_T, _p, _fractions[0]);};
    /// Get the thermal conductivity [W/m/K] (based on the temperature and pressure in the state class)
    long double calc_conductivity(void){return fluid->cond(_T, _p, _fractions[0]);};

    long double calc_rhomass(void){return fluid->rho(_T, _p, _fractions[0]);};
    long double calc_hmass(void){return fluid->h(_T, _p, _fractions[0]);};
    long double calc_smass(void){return fluid->s(_T, _p, _fractions[0]);};
    long double calc_umass(void){return fluid->u(_T, _p, _fractions[0]);};
    long double calc_cpmass(void){return fluid->cp(_T, _p, _fractions[0]);};
    long double calc_cvmass(void){return fluid->cv(_T, _p, _fractions[0]);};
    long double calc_Tmax(void){return fluid->getTmax();};
    long double calc_Tmin(void){return fluid->getTmin();};
    
    long double calc_fraction_min(void){return fluid->getxmin();};
    long double calc_fraction_max(void){return fluid->getxmax();};
    long double calc_T_freeze(void){
        return fluid->Tfreeze(_p, _fractions[0]);};
		
	long double calc_melting_line(int param, int given, long double value){
		if (param == iT && given == iP){
			return fluid->Tfreeze(value, _fractions[0]);
		}
		else{
			throw ValueError("For incompressibles, the only valid inputs to calc_melting_line are T(p)");
		}
	};
        
    std::string calc_name(void){return fluid->getDescription();};
};

} /* namespace CoolProp */
#endif /* INCOMPRESSIBLEBACKEND_H_ */
