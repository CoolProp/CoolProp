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

#include "CoolPropTools.h"
#include "DataStructures.h"
#include "GeneralizedCubic.h"
#include "CubicsLibrary.h"
#include "Configuration.h"
#include "AbstractState.h"
#include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#include "Exceptions.h"
#include <vector>

namespace CoolProp {

// Forward declaration for use in initialization of AbstractCubicBackend
class CubicResidualHelmholtz;

class AbstractCubicBackend : public HelmholtzEOSMixtureBackend
{
   protected:
    shared_ptr<AbstractCubic> cubic;
    std::vector<CubicLibrary::CubicsValues> components;  ///< The components that are in use
   public:
    /// Set the pointer to the residual helmholtz class, etc.
    void setup(bool generate_SatL_and_SatV = true);

    /// Set the alpha function based on the alpha function defined in the components vector;
    void set_alpha_from_components();

    /// Set the non-dimensionalized Helmholtz energy based on the fluids defined in the components vector
    void set_alpha0_from_components();

    /// Get a reference to the shared pointer managing the generalized cubic class
    shared_ptr<AbstractCubic>& get_cubic() {
        return cubic;
    };

    std::vector<std::string> calc_fluid_names(void);

    bool using_mole_fractions(void) {
        return true;
    };
    bool using_mass_fractions(void) {
        return false;
    };
    bool using_volu_fractions(void) {
        return false;
    };

    void set_mass_fractions(const std::vector<CoolPropDbl>& mass_fractions) {
        throw NotImplementedError("Mass composition has not been implemented.");
    };
    void set_volu_fractions(const std::vector<CoolPropDbl>& volu_fractions) {
        throw NotImplementedError("Volume composition has not been implemented.");
    };
    const std::vector<CoolPropDbl>& get_mole_fractions(void) {
        return this->mole_fractions;
    };

    const double get_fluid_constant(std::size_t i, parameters param) const {
        switch (param) {
            case iP_critical:
                return cubic->get_pc()[i];
            case iT_reducing:
            case iT_critical:
                return cubic->get_Tc()[i];
            case iacentric_factor:
                return cubic->get_acentric()[i];
            case imolar_mass:
                return components[i].molemass;
            case iT_triple:
                return HelmholtzEOSMixtureBackend::get_components()[i].EOS().sat_min_liquid.T;  // From the base class data structure
            case iP_triple:
                return HelmholtzEOSMixtureBackend::get_components()[i].EOS().sat_min_liquid.p;  // From the base class data structure
            case irhomolar_reducing:
            case irhomolar_critical:
                return components[i].rhomolarc;
            case igas_constant:
                return get_config_double(R_U_CODATA);
            default:
                throw ValueError(format("I don't know what to do with this fluid constant: %s", get_parameter_information(param, "short").c_str()));
        }
    }

    /// Calculate the gas constant in J/mol/K
    CoolPropDbl calc_gas_constant(void) {
        return cubic->get_R_u();
    };
    /// Get the reducing state to be used
    SimpleState calc_reducing_state_nocache(const std::vector<CoolPropDbl>& mole_fractions) {
        SimpleState reducing;
        reducing.T = cubic->get_Tr();
        reducing.rhomolar = cubic->get_rhor();
        return reducing;
    };
    CoolPropDbl calc_reduced_density(void) {
        return _rhomolar / get_cubic()->get_rhor();
    };
    CoolPropDbl calc_reciprocal_reduced_temperature(void) {
        return get_cubic()->get_Tr() / _T;
    };
    std::vector<double> spinodal_densities();

    CoolPropDbl calc_T_critical(void) {
        if (is_pure_or_pseudopure) {
            return cubic->get_Tc()[0];
        } else {
            return HelmholtzEOSMixtureBackend::calc_T_critical();
        }
    };
    CoolPropDbl calc_p_critical(void) {
        if (is_pure_or_pseudopure) {
            return cubic->get_pc()[0];
        } else {
            return HelmholtzEOSMixtureBackend::calc_p_critical();
        }
    };
    CoolPropDbl calc_acentric_factor(void) {
        if (is_pure_or_pseudopure) {
            return cubic->get_acentric()[0];
        } else {
            throw ValueError("acentric factor cannot be calculated for mixtures");
        }
    }
    CoolPropDbl calc_rhomolar_critical(void) {
        if (is_pure_or_pseudopure) {
            // Curve fit from all the pure fluids in CoolProp (thanks to recommendation of A. Kazakov)
            double v_c_Lmol = 2.14107171795 * (cubic->get_Tc()[0] / cubic->get_pc()[0] * 1000) + 0.00773144012514;  // [L/mol]
            return 1 / (v_c_Lmol / 1000.0);
        } else {
            return HelmholtzEOSMixtureBackend::calc_rhomolar_critical();
        }
    };

    /// \brief Get linear mole fraction weighting of the critical molar volumes and temperatures
    /// these are used in te
    void get_linear_reducing_parameters(double& rhomolar, double& T);

    /// Get the the starting values for the critical point evaluation routines
    void get_critical_point_starting_values(double& delta0, double& tau0);

    /// Get the search radius in delta and tau for the tracer, scaled appropriately for cubic
    void get_critical_point_search_radii(double& R_delta, double& R_tau);

    /// Checking function to see if we should stop the tracing of the critical contour
    bool get_critical_is_terminated(double& delta, double& tau);

    CoolPropDbl calc_alphar_deriv_nocache(const int nTau, const int nDelta, const std::vector<CoolPropDbl>& mole_fractions, const CoolPropDbl& tau,
                                          const CoolPropDbl& delta);

    /// Calculate the pressure in most computationally efficient manner
    CoolPropDbl calc_pressure_nocache(CoolPropDbl T, CoolPropDbl rhomolar);

    /// Update the state for DT inputs if phase is imposed. Otherwise delegate to base class
    virtual void update_DmolarT();

    virtual void update(CoolProp::input_pairs input_pair, double value1, double value2);

    /** Use the cubic EOS to solve for density
     *
     * \f[
     * \left[
     * \begin{array}{l}
     * Z^3 \\
     * +Z^2[B(\Delta_1+\Delta_2-1)-1]
     * +Z[A-B(\Delta_1+\Delta_2)-B^2(\Delta_1+\Delta_2-\Delta_1\Delta_2)]
     * +[-AB-\Delta_1\Delta2(-B^2-B^3)]
     * \right]
     * \f]
     * with
     * \f[ A = \frac{ap}{R^2T^2} \f]
     * \f[ B = \frac{bp}{RT} \f]
     *
     * Sympy code:
     * R,T,v,b,a,Delta_1,Delta_2,p,Z,A,B = symbols('R,T,v,b,a,Delta_1,Delta_2,p,Z,A,B')
     * eqn = (R*T/(v-b)-a/(v+Delta_1*b)/(v+Delta_2*b)-p).subs(v,Z*R*T/p)
     * eqn2 = eqn*(R*T*Z/p-b)*(R*T*Z/p+Delta_2*b)*(R*T*Z/p+Delta_1*b)/(R**3*T**3/p**2)
     * collect(expand(factor(-eqn2)),Z).subs(b*p/(R*T),B).subs(a*p/(R**2*T**2),A)
     *
     */
    void rho_Tp_cubic(CoolPropDbl T, CoolPropDbl p, int& Nsolns, double& rho0, double& rho1, double& rho2);

    /// In this class, we are already doing cubic evaluation, just delegate to our function
    CoolPropDbl solver_rho_Tp_SRK(CoolPropDbl T, CoolPropDbl p, phases phase) {
        return solver_rho_Tp(T, p);
    };
    /**
     * /brief Solve for rho = f(T,p)
     *
     * You can often get three solutions, to overcome this problem you must either specify the phase, or provide a reasonable guess value for rho_guess, but not both
     */
    CoolPropDbl solver_rho_Tp(CoolPropDbl T, CoolPropDbl p, CoolPropDbl rho_guess = -1);

    CoolPropDbl solver_rho_Tp_global(CoolPropDbl T, CoolPropDbl p, CoolPropDbl rhomax);

    /// Update the state used to calculate the tangent-plane-distance
    void update_TPD_state() {
        TPD_state.reset(get_copy());
    };

    /// Cubic backend flashes for PQ, and QT
    void saturation(CoolProp::input_pairs inputs);

    CoolPropDbl calc_molar_mass(void);

    void set_binary_interaction_double(const std::size_t i1, const std::size_t i2, const std::string& parameter, const double value);
    double get_binary_interaction_double(const std::size_t i1, const std::size_t i2, const std::string& parameter);

    void set_binary_interaction_double(const std::string& CAS1, const std::string& CAS2, const std::string& parameter, const double value) {
        throw ValueError("set_binary_interaction_double not defined for AbstractCubic not defined for CAS #");
    }
    double get_binary_interaction_double(const std::string& CAS1, const std::string& CAS2, const std::string& parameter) {
        throw ValueError("get_binary_interaction_double not defined for AbstractCubic not defined for CAS #");
    };

    // Return a 1-1 copy of this class
    virtual HelmholtzEOSMixtureBackend* get_copy(bool generate_SatL_and_SatV = true) = 0;

    // Copy the entire kij matrix from another instance in one shot
    void copy_k(AbstractCubicBackend* donor);

    //
    void copy_all_alpha_functions(AbstractCubicBackend* donor);

    /// Copy the internals from another class into this one (kij, alpha functions, cp0 functions, etc.)
    void copy_internals(AbstractCubicBackend& donor);

    // Set the cubic alpha function's constants:
    void set_cubic_alpha_C(const size_t i, const std::string& parameter, const double c1, const double c2, const double c3);

    // Set fluid parameter (currently the volume translation parameter)
    void set_fluid_parameter_double(const size_t i, const std::string& parameter, const double value);

    // Get fluid parameter (currently the volume translation parameter)
    double get_fluid_parameter_double(const size_t i, const std::string& parameter);
};

class SRKBackend : public AbstractCubicBackend
{

   public:
    SRKBackend(const std::vector<double>& Tc, const std::vector<double>& pc, const std::vector<double>& acentric, double R_u,
               bool generate_SatL_and_SatV = true) {
        cubic.reset(new SRK(Tc, pc, acentric, R_u));
        setup(generate_SatL_and_SatV);
    };
    SRKBackend(double Tc, double pc, double acentric, double R_u, bool generate_SatL_and_SatV = true) {
        cubic.reset(new SRK(Tc, pc, acentric, R_u));
        setup(generate_SatL_and_SatV);
    }
    SRKBackend(const std::vector<std::string> fluid_identifiers, const double R_u = get_config_double(R_U_CODATA),
               bool generate_SatL_and_SatV = true) {
        std::vector<double> Tc, pc, acentric;
        N = fluid_identifiers.size();
        components.resize(N);
        for (std::size_t i = 0; i < fluid_identifiers.size(); ++i) {
            components[i] = CubicLibrary::get_cubic_values(fluid_identifiers[i]);
            Tc.push_back(components[i].Tc);
            pc.push_back(components[i].pc);
            acentric.push_back(components[i].acentric);
        }
        cubic.reset(new SRK(Tc, pc, acentric, R_u));
        setup(generate_SatL_and_SatV);
    }
    HelmholtzEOSMixtureBackend* get_copy(bool generate_SatL_and_SatV = true) {
        AbstractCubicBackend* ACB = new SRKBackend(cubic->get_Tc(), cubic->get_pc(), cubic->get_acentric(), cubic->get_R_u(), generate_SatL_and_SatV);
        ACB->copy_internals(*this);
        return static_cast<HelmholtzEOSMixtureBackend*>(ACB);
    }
    std::string backend_name(void) {
        return get_backend_string(SRK_BACKEND);
    }
};

class PengRobinsonBackend : public AbstractCubicBackend
{

   public:
    PengRobinsonBackend(){};  // Default constructor (make sure you know what you are doing)
    PengRobinsonBackend(const std::vector<double>& Tc, const std::vector<double>& pc, const std::vector<double>& acentric, double R_u,
                        bool generate_SatL_and_SatV = true) {
        cubic.reset(new PengRobinson(Tc, pc, acentric, R_u));
        setup(generate_SatL_and_SatV);
    };
    PengRobinsonBackend(double Tc, double pc, double acentric, double R_u, bool generate_SatL_and_SatV = true) {
        cubic.reset(new PengRobinson(Tc, pc, acentric, R_u));
        setup(generate_SatL_and_SatV);
    };
    PengRobinsonBackend(const std::vector<std::string> fluid_identifiers, const double R_u = get_config_double(R_U_CODATA),
                        bool generate_SatL_and_SatV = true) {
        std::vector<double> Tc, pc, acentric;
        N = fluid_identifiers.size();
        components.resize(N);
        for (std::size_t i = 0; i < fluid_identifiers.size(); ++i) {
            components[i] = CubicLibrary::get_cubic_values(fluid_identifiers[i]);
            Tc.push_back(components[i].Tc);
            pc.push_back(components[i].pc);
            acentric.push_back(components[i].acentric);
        }
        cubic.reset(new PengRobinson(Tc, pc, acentric, R_u));
        setup(generate_SatL_and_SatV);
    };
    HelmholtzEOSMixtureBackend* get_copy(bool generate_SatL_and_SatV = true) {
        AbstractCubicBackend* ACB =
          new PengRobinsonBackend(cubic->get_Tc(), cubic->get_pc(), cubic->get_acentric(), cubic->get_R_u(), generate_SatL_and_SatV);
        ACB->copy_internals(*this);
        return static_cast<HelmholtzEOSMixtureBackend*>(ACB);
    }
    std::string backend_name(void) {
        return get_backend_string(PR_BACKEND);
    }
};

/**
 * This class implements all the derivatives of the Helmholtz energy (as well as composition derivatives) that are required for the
 * cubic backends
 */
class CubicResidualHelmholtz : public ResidualHelmholtz
{

   protected:
    AbstractCubicBackend* ACB;

   public:
    CubicResidualHelmholtz() {
        ACB = NULL;
    };
    CubicResidualHelmholtz(AbstractCubicBackend* ACB) : ACB(ACB){};

    // copy assignment
    CubicResidualHelmholtz& operator=(CubicResidualHelmholtz& other) {
        ACB = other.ACB;
        return *this;
    }

    /// All the derivatives of the residual Helmholtz energy w.r.t. tau and delta that do not involve composition derivative
    virtual HelmholtzDerivatives all(HelmholtzEOSMixtureBackend& HEOS, const std::vector<CoolPropDbl>& mole_fractions, double tau, double delta,
                                     bool cache_values = false) {
        HelmholtzDerivatives a;
        std::vector<double> z = std::vector<double>(mole_fractions.begin(), mole_fractions.end());
        shared_ptr<AbstractCubic>& cubic = ACB->get_cubic();
        a.alphar = cubic->alphar(tau, delta, z, 0, 0);
        a.dalphar_dtau = cubic->alphar(tau, delta, z, 1, 0);
        a.dalphar_ddelta = cubic->alphar(tau, delta, z, 0, 1);
        a.d2alphar_dtau2 = cubic->alphar(tau, delta, z, 2, 0);
        a.d2alphar_ddelta_dtau = cubic->alphar(tau, delta, z, 1, 1);
        a.d2alphar_ddelta2 = cubic->alphar(tau, delta, z, 0, 2);
        a.d3alphar_dtau3 = cubic->alphar(tau, delta, z, 3, 0);
        a.d3alphar_ddelta_dtau2 = cubic->alphar(tau, delta, z, 2, 1);
        a.d3alphar_ddelta2_dtau = cubic->alphar(tau, delta, z, 1, 2);
        a.d3alphar_ddelta3 = cubic->alphar(tau, delta, z, 0, 3);
        a.d4alphar_dtau4 = cubic->alphar(tau, delta, z, 4, 0);
        a.d4alphar_ddelta_dtau3 = cubic->alphar(tau, delta, z, 3, 1);
        a.d4alphar_ddelta2_dtau2 = cubic->alphar(tau, delta, z, 2, 2);
        a.d4alphar_ddelta3_dtau = cubic->alphar(tau, delta, z, 1, 3);
        a.d4alphar_ddelta4 = cubic->alphar(tau, delta, z, 0, 4);
        return a;
    }
    virtual CoolPropDbl dalphar_dxi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
        return ACB->get_cubic()->d_alphar_dxi(HEOS.tau(), HEOS.delta(), HEOS.get_mole_fractions_doubleref(), 0, 0, i, xN_flag == XN_INDEPENDENT);
    }
    virtual CoolPropDbl d2alphar_dxi_dTau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
        return ACB->get_cubic()->d_alphar_dxi(HEOS.tau(), HEOS.delta(), HEOS.get_mole_fractions_doubleref(), 1, 0, i, xN_flag == XN_INDEPENDENT);
    }
    virtual CoolPropDbl d2alphar_dxi_dDelta(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
        return ACB->get_cubic()->d_alphar_dxi(HEOS.tau(), HEOS.delta(), HEOS.get_mole_fractions_doubleref(), 0, 1, i, xN_flag == XN_INDEPENDENT);
    }
    virtual CoolPropDbl d3alphar_dxi_dTau2(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
        return ACB->get_cubic()->d_alphar_dxi(HEOS.tau(), HEOS.delta(), HEOS.get_mole_fractions_doubleref(), 2, 0, i, xN_flag == XN_INDEPENDENT);
    }
    virtual CoolPropDbl d3alphar_dxi_dDelta_dTau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
        return ACB->get_cubic()->d_alphar_dxi(HEOS.tau(), HEOS.delta(), HEOS.get_mole_fractions_doubleref(), 1, 1, i, xN_flag == XN_INDEPENDENT);
    }
    virtual CoolPropDbl d3alphar_dxi_dDelta2(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
        return ACB->get_cubic()->d_alphar_dxi(HEOS.tau(), HEOS.delta(), HEOS.get_mole_fractions_doubleref(), 0, 2, i, xN_flag == XN_INDEPENDENT);
    }

    virtual CoolPropDbl d2alphardxidxj(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) {
        return ACB->get_cubic()->d2_alphar_dxidxj(HEOS.tau(), HEOS.delta(), HEOS.get_mole_fractions_doubleref(), 0, 0, i, j,
                                                  xN_flag == XN_INDEPENDENT);
    }
    virtual CoolPropDbl d3alphar_dxi_dxj_dTau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) {
        return ACB->get_cubic()->d2_alphar_dxidxj(HEOS.tau(), HEOS.delta(), HEOS.get_mole_fractions_doubleref(), 1, 0, i, j,
                                                  xN_flag == XN_INDEPENDENT);
    }
    virtual CoolPropDbl d3alphar_dxi_dxj_dDelta(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) {
        return ACB->get_cubic()->d2_alphar_dxidxj(HEOS.tau(), HEOS.delta(), HEOS.get_mole_fractions_doubleref(), 0, 1, i, j,
                                                  xN_flag == XN_INDEPENDENT);
    }

    virtual CoolPropDbl d4alphar_dxi_dTau3(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
        return ACB->get_cubic()->d_alphar_dxi(HEOS.tau(), HEOS.delta(), HEOS.get_mole_fractions_doubleref(), 3, 0, i, xN_flag == XN_INDEPENDENT);
    }
    virtual CoolPropDbl d4alphar_dxi_dDelta2_dTau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
        return ACB->get_cubic()->d_alphar_dxi(HEOS.tau(), HEOS.delta(), HEOS.get_mole_fractions_doubleref(), 1, 2, i, xN_flag == XN_INDEPENDENT);
    }
    virtual CoolPropDbl d4alphar_dxi_dDelta_dTau2(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
        return ACB->get_cubic()->d_alphar_dxi(HEOS.tau(), HEOS.delta(), HEOS.get_mole_fractions_doubleref(), 2, 1, i, xN_flag == XN_INDEPENDENT);
    }
    virtual CoolPropDbl d4alphar_dxi_dDelta3(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
        return ACB->get_cubic()->d_alphar_dxi(HEOS.tau(), HEOS.delta(), HEOS.get_mole_fractions_doubleref(), 0, 3, i, xN_flag == XN_INDEPENDENT);
    }
    virtual CoolPropDbl d4alphar_dxi_dxj_dTau2(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) {
        return ACB->get_cubic()->d2_alphar_dxidxj(HEOS.tau(), HEOS.delta(), HEOS.get_mole_fractions_doubleref(), 2, 0, i, j,
                                                  xN_flag == XN_INDEPENDENT);
    }
    virtual CoolPropDbl d4alphar_dxi_dxj_dDelta_dTau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) {
        return ACB->get_cubic()->d2_alphar_dxidxj(HEOS.tau(), HEOS.delta(), HEOS.get_mole_fractions_doubleref(), 1, 1, i, j,
                                                  xN_flag == XN_INDEPENDENT);
    }
    virtual CoolPropDbl d4alphar_dxi_dxj_dDelta2(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) {
        return ACB->get_cubic()->d2_alphar_dxidxj(HEOS.tau(), HEOS.delta(), HEOS.get_mole_fractions_doubleref(), 0, 2, i, j,
                                                  xN_flag == XN_INDEPENDENT);
    }
    virtual CoolPropDbl d3alphardxidxjdxk(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, std::size_t k,
                                          x_N_dependency_flag xN_flag) {
        return ACB->get_cubic()->d3_alphar_dxidxjdxk(HEOS.tau(), HEOS.delta(), HEOS.get_mole_fractions_doubleref(), 0, 0, i, j, k,
                                                     xN_flag == XN_INDEPENDENT);
    }
    virtual CoolPropDbl d4alphar_dxi_dxj_dxk_dTau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, std::size_t k,
                                                  x_N_dependency_flag xN_flag) {
        return ACB->get_cubic()->d3_alphar_dxidxjdxk(HEOS.tau(), HEOS.delta(), HEOS.get_mole_fractions_doubleref(), 1, 0, i, j, k,
                                                     xN_flag == XN_INDEPENDENT);
    }
    virtual CoolPropDbl d4alphar_dxi_dxj_dxk_dDelta(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, std::size_t k,
                                                    x_N_dependency_flag xN_flag) {
        return ACB->get_cubic()->d3_alphar_dxidxjdxk(HEOS.tau(), HEOS.delta(), HEOS.get_mole_fractions_doubleref(), 0, 1, i, j, k,
                                                     xN_flag == XN_INDEPENDENT);
    }
};

} /* namespace CoolProp */
#endif /* CUBICBACKEND_H_ */
