/*

The goal of this backend is to allow the Helmholtz-based structure for cubics and to replace the entire
multi-fluid model with a one-fluid model.  The other changes are relatively trivial.  The primary 
change is to replace the core residual Helmholtz energy derivatives from HelmholtzEOSMixtureBackend
with the derivatives from this class. 

The core code for the Helmholtz translations is from the publication 
"Helmholtz energy translations for common cubic equations of state for use in one-fluid and multi-fluid mixture models"
by Ian H. Bell and Andreas Jaeger, J. Res. NIST, 2016

*/

#ifndef CUBICBACKEND_H_
#define CUBICBACKEND_H_

#include "DataStructures.h"
#include "GeneralizedCubic.h"
#include "AbstractState.h"
#include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#include "Exceptions.h"
#include <vector>

namespace CoolProp {

// Forward declaration for use in initalization of AbstractCubicBackend
class CubicResidualHelmholtz;

class AbstractCubicBackend : public HelmholtzEOSMixtureBackend  {

protected:
    shared_ptr<AbstractCubic> cubic;
public:
	
	/// Set the pointer to the residual helmholtz class, etc.
	void setup();

	/// Get a reference to the shared pointer managing the generalized cubic class
	shared_ptr<AbstractCubic> &get_cubic(){ return cubic; };
	
    bool using_mole_fractions(void){return true;};
    bool using_mass_fractions(void){return false;}; 
    bool using_volu_fractions(void){return false;};

    void set_mole_fractions(const std::vector<CoolPropDbl> &mole_fractions){throw NotImplementedError("Mole composition has not been implemented.");};
    void set_mass_fractions(const std::vector<CoolPropDbl> &mass_fractions){}; // Not implemented, but don't throw any errors
    void set_volu_fractions(const std::vector<CoolPropDbl> &volu_fractions){throw NotImplementedError("Volume composition has not been implemented.");};
    const std::vector<CoolPropDbl> & get_mole_fractions(void){throw NotImplementedError("get_mole_fractions composition has not been implemented.");};

	/// Calculate the gas constant in J/mol/K
	CoolPropDbl calc_gas_constant(void){
		return 8.314498; //TODO: get from class
	};
	/// Get the reducing state to be used
	SimpleState calc_reducing_state_nocache(const std::vector<CoolPropDbl> & mole_fractions)
	{
		SimpleState reducing;
		reducing.T = cubic->T_r;
		reducing.rhomolar = cubic->rho_r;
		return reducing;
	};
};

class SRKBackend : public AbstractCubicBackend  {

public:
	SRKBackend(const std::vector<double> &Tc, 
		       const std::vector<double> &pc, 
		       const std::vector<double> &acentric,
               double R_u) {
        cubic.reset(new SRK(Tc, pc, acentric, R_u));
		setup();
		
    };
	SRKBackend(double Tc, 
		       double pc, 
		       double acentric,
               double R_u) {
        cubic.reset(new SRK(Tc, pc, acentric, R_u));
		is_pure_or_pseudopure = true;
		setup();
    }
};

class PengRobinsonBackend : public AbstractCubicBackend  {

public:
	PengRobinsonBackend(const std::vector<double> &Tc, 
		       const std::vector<double> &pc, 
		       const std::vector<double> &acentric,
               double R_u) {
        cubic.reset(new PengRobinson(Tc, pc, acentric, R_u));
		setup();
    };
	PengRobinsonBackend(double Tc, 
		       double pc, 
		       double acentric,
               double R_u) {
        cubic.reset(new PengRobinson(Tc, pc, acentric, R_u));
		is_pure_or_pseudopure = true;
		setup();
    }
};

/**
 * This class implements all the derivatives of the Helmholtz energy (as well as composition derivatives) that are required for the
 * cubic backends
 */
class CubicResidualHelmholtz : public ResidualHelmholtz
{

protected:
	AbstractCubicBackend *ACB;
public:
	CubicResidualHelmholtz(){ ACB = NULL; };
	CubicResidualHelmholtz(AbstractCubicBackend * ACB) : ACB(ACB) {};

    /// All the derivatives of the residual Helmholtz energy w.r.t. tau and delta that do not involve composition derivative
    virtual HelmholtzDerivatives all(HelmholtzEOSMixtureBackend &HEOS, const std::vector<CoolPropDbl> &mole_fractions, bool cache_values = false)
    {
		HelmholtzDerivatives a;
		std::vector<double> z = std::vector<double>(mole_fractions.begin(), mole_fractions.end());
        shared_ptr<AbstractCubic> &cubic = ACB->get_cubic();
		double tau = HEOS.tau(), delta = HEOS.delta();
		a.alphar = cubic->alphar(tau, delta, z, 0, 0);
		a.dalphar_dtau = cubic->alphar(tau, delta, z, 1, 0);
		a.dalphar_ddelta = cubic->alphar(tau, delta, z, 0, 1);
        return a;
    }
    virtual CoolPropDbl dalphar_dxi(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag){
        return 0;
    }
	virtual CoolPropDbl d2alphar_dxi_dTau(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag){
        return 0;
    }
    virtual CoolPropDbl d2alphar_dxi_dDelta(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag){
        return 0;
    }
	virtual CoolPropDbl d3alphar_dxi_dTau2(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag){
        return 0;
    }
    virtual CoolPropDbl d3alphar_dxi_dDelta_dTau(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag){
        return 0;
    }
    virtual CoolPropDbl d3alphar_dxi_dDelta2(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag){
        return 0;
    }

    virtual CoolPropDbl d2alphardxidxj(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag){
        return 0;
    }
    virtual CoolPropDbl d3alphar_dxi_dxj_dTau(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag  xN_flag){
        return 0;
    }
    virtual CoolPropDbl d3alphar_dxi_dxj_dDelta(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag){
        return 0;
    }
    

    virtual CoolPropDbl d4alphar_dxi_dTau3(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag){
        return 0;
    }
    virtual CoolPropDbl d4alphar_dxi_dDelta2_dTau(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag){
        return 0;
    }
    virtual CoolPropDbl d4alphar_dxi_dDelta_dTau2(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag){
        return 0;
    }
    virtual CoolPropDbl d4alphar_dxi_dDelta3(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag){
        return 0;
    }
    virtual CoolPropDbl d4alphar_dxi_dxj_dTau2(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag  xN_flag){
        return 0;
    }
    virtual CoolPropDbl d4alphar_dxi_dxj_dDelta_dTau(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag  xN_flag){
        return 0;
    }
    virtual CoolPropDbl d4alphar_dxi_dxj_dDelta2(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag){
        return 0;
    }
	virtual CoolPropDbl d3alphardxidxjdxk(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, std::size_t k, x_N_dependency_flag xN_flag){
        return 0;
    }
};

} /* namespace CoolProp */
#endif /* IF97BACKEND_H_ */
