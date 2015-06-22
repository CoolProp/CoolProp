
#ifndef IF97BACKEND_H_
#define IF97BACKEND_H_

#include "DataStructures.h"
#include "IF97.h"
#include "AbstractState.h"
#include "Exceptions.h"
#include <vector>

namespace CoolProp {

class IF97Backend : public AbstractState  {

protected:

	/// Additional cached elements used for the partial derivatives
	CachedElement  _cpmass, _cvmass, _hmass, _rhomass, _smass, _umass;

public:
    std::string backend_name(void){return "IF97Backend";}

    // REQUIRED BUT NOT USED FUNCTIONS
    bool using_mole_fractions(void){return false;};
    bool using_mass_fractions(void){return true;}; // But actually it doesn't matter since it is only a pure fluid
    bool using_volu_fractions(void){return false;};
    void set_mole_fractions(const std::vector<CoolPropDbl> &mole_fractions){throw NotImplementedError("Mole composition has not been implemented.");};
    void set_mass_fractions(const std::vector<CoolPropDbl> &mass_fractions){}; // Not implemented, but don't throw any errors
    void set_volu_fractions(const std::vector<CoolPropDbl> &volu_fractions){throw NotImplementedError("Volume composition has not been implemented.");};
    const std::vector<CoolPropDbl> & get_mole_fractions(void){throw NotImplementedError("get_mole_fractions composition has not been implemented.");};

    /// Updating function for incompressible fluid
    /**
    In this function we take a pair of thermodynamic states, those defined in the input_pairs
    enumeration and update all the internal variables that we can.

    @param input_pair Integer key from CoolProp::input_pairs to the two inputs that will be passed to the function
    @param value1 First input value
    @param value2 Second input value
    */
    void update(CoolProp::input_pairs input_pair, double value1, double value2){
        switch(input_pair){
            case PT_INPUTS: _p = value1; _T = value2; break;
            default:
                throw ValueError("bad input_pair");
        }
    };

    /// Clear all the cached values
    bool clear();

    /// Check if the mole fractions have been set, etc.
    void check_status();

    /** We have to override some of the functions from the AbstractState.
	 *  The incompressibles are only mass-based and do not support conversion
	 *  from molar to specific quantities.
	 *  We also have a few new cached variables that we need.
	 */
	/// Return the mass density in kg/m^3
    double rhomass(void){ return calc_rhomass(); }
    CoolPropDbl calc_rhomass(void){ return IF97::rhomass_Tp(_T, _p); }
	/// Return the mass enthalpy in J/kg
	double hmass(void){return calc_hmass();}
    CoolPropDbl calc_hmass(void){ return IF97::hmass_Tp(_T, _p); }
	/// Return the molar entropy in J/mol/K
	double smass(void){return calc_smass();}
    CoolPropDbl calc_smass(void){ return IF97::smass_Tp(_T, _p); }
	/// Return the molar internal energy in J/mol
	double umass(void){return calc_umass();}
    CoolPropDbl calc_umass(void){ return IF97::umass_Tp(_T, _p); }
	/// Return the mass constant pressure specific heat in J/kg/K
	double cpmass(void){return calc_cpmass();}
    CoolPropDbl calc_cpmass(void){ return IF97::cpmass_Tp(_T, _p); }
    /// Return the mass constant volume specific heat in J/kg/K
	double cvmass(void){return calc_cvmass();}
    CoolPropDbl calc_cvmass(void){ return IF97::cvmass_Tp(_T, _p); }
};

} /* namespace CoolProp */
#endif /* IF97BACKEND_H_ */
