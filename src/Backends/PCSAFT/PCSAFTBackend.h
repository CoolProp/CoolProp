#ifndef PCSAFTBACKEND_H_
#define PCSAFTBACKEND_H_

#include "CoolPropTools.h"
#include "DataStructures.h"
#include "PCSAFTLibrary.h"
#include "Configuration.h"
#include "AbstractState.h"
#include "Exceptions.h"
#include <vector>

using namespace std;

namespace CoolProp {

const static double kb = 1.380648465952442093e-23; // Boltzmann constant, J K^-1
const static double PI = 3.141592653589793;
const static double N_AV = 6.022140857e23; // Avagadro's number
const static double E_CHRG = 1.6021766208e-19; // elementary charge, units of coulomb
const static double perm_vac = 8.854187817e-22; //permittivity in vacuum, C V^-1 Angstrom^-1

class PCSAFTBackend : public AbstractState  {

protected:
    std::vector<PCSAFTFluid> components; ///< The components that are in use
    std::vector<double> k_ij; ///< binary interaction parameters
    bool is_pure_or_pseudopure; ///< A flag for whether the substance is a pure or pseudo-pure fluid (true) or a mixture (false)
    std::vector<CoolPropDbl> mole_fractions; ///< The bulk mole fractions of the mixture
    std::vector<double> mole_fractions_double; ///< A copy of the bulk mole fractions of the mixture stored as doubles
    std::vector<CoolPropDbl> K, ///< The K factors for the components
                             lnK; ///< The natural logarithms of the K factors of the components
    double dielc; ///< The dielectric constant of the solvent, if ion term is used

    std::size_t N; ///< Number of components

    bool ion_term; ///< Whether or not the ion term should be included
    bool polar_term; ///< Whether or not the dipole term should be included
    bool assoc_term; ///< Whether or not the association term should be included

    void post_update(bool optional_checks = true);

    // these functions are used internally to solve for association parameters
    vector<double> XA_find(vector<double> XA_guess, int ncomp, vector<double> delta_ij, double den,
        vector<double> x);
    vector<double> dXA_find(int ncA, int ncomp, vector<int> iA, vector<double> delta_ij,
        double den, vector<double> XA, vector<double> ddelta_dd, vector<double> x, int n_sites);
    vector<double> dXAdt_find(int ncA, vector<double> delta_ij, double den,
        vector<double> XA, vector<double> ddelta_dt, vector<double> x, int n_sites);
    double dielc_water(double t);

public:
    PCSAFTBackend(const std::vector<std::string> &component_names);
    PCSAFTBackend(const std::vector<PCSAFTFluid> &components);

    /// The name of the backend being used
    std::string backend_name(void) { return get_backend_string(PCSAFT_BACKEND); }

    bool using_mole_fractions(void){return true;};
    bool using_mass_fractions(void){return false;};
    bool using_volu_fractions(void){return false;};

    void set_mass_fractions(const std::vector<CoolPropDbl> &mass_fractions);
    void set_volu_fractions(const std::vector<CoolPropDbl> &volu_fractions){throw NotImplementedError("Volume composition has not been implemented.");};
    void set_mole_fractions(const std::vector<CoolPropDbl> &mole_fractions);
    const std::vector<CoolPropDbl> & get_mole_fractions(void){ return this->mole_fractions; };

    void resize(std::size_t N);

    virtual void update(CoolProp::input_pairs input_pair, double value1, double value2); // %%checked

    const double get_fluid_constant(std::size_t i, parameters param) const{
        // const PCSAFTFluid &fld = components[i];
        // switch(param){
        //     case im: return fld.m;
        //     case isigma: return fld.sigma;
        //     case iu: return fld.u;
        //     case iuAB: return fld.uAB;
        //     case ikappa: return fld.kappa;
        //     case idipm: return fld.dipm;
        //     case idipnum: return fld.dipnum;
        //     case iacentric_factor: return fld.EOS().acentric;
        //     case imolar_mass: return fld.EOS().molar_mass;
        //     default:
        //         throw ValueError(format("I don't know what to do with this fluid constant: %s", get_parameter_information(param,"short").c_str()));
        // }
    }

    // void set_binary_interaction_double(const std::size_t i, const std::size_t j, const std::string &parameter, const double value);
    // double get_binary_interaction_double(const std::size_t i, const std::size_t j, const std::string &parameter);

    /*  We have to override some of the functions from the AbstractState.
     *  IF97 is only mass-based and does not support conversion
     *  from mass- to molar-specific quantities.
     */
    // ************************************************************************* //
    //                   Basic Thermodynamic Functions                           //
    // ************************************************************************* //
    //
    /// Calculate the pressure in most computationally efficient manner
    CoolPropDbl calc_pressure_nocache(CoolPropDbl t, CoolPropDbl rho);

    /// Update the state for DT inputs if phase is imposed. Otherwise delegate to base class
    virtual void update_DmolarT();

    CoolPropDbl calc_alpha0(void); // ideal gas helmholtz energy term
    CoolPropDbl calc_alphar(void); // residual helmholtz energy
    CoolPropDbl calc_hmolar_nocache(CoolPropDbl T, CoolPropDbl rhomolar);
    CoolPropDbl calc_hmolar(void);
    CoolPropDbl calc_smolar_nocache(CoolPropDbl T, CoolPropDbl rhomolar);
    CoolPropDbl calc_smolar(void);
    CoolPropDbl calc_fugacity_coefficient(void);
    CoolPropDbl calc_gibbsmolar(void);
    CoolPropDbl calc_cpmolar(void);
    // CoolPropDbl calc_cp0molar(void);
    CoolPropDbl calc_compressibility_factor(double t, double rho);

    CoolPropDbl solver_rho_Tp(CoolPropDbl T, CoolPropDbl p, CoolPropDbl rhomolar_guess);
    double calc_SatLiquid(bool is_temperature_input, CoolPropDbl value1);
    double calc_SatVapor(bool is_temperature_input, CoolPropDbl value1);
    void flash_QT(PCSAFTBackend &PCSAFT);
    void flash_PQ(PCSAFTBackend &PCSAFT);

    phases calc_phase(void){return _phase;};
    /** \brief Specify the phase - this phase will always be used in calculations
     *
     * @param phase_index The index from CoolProp::phases
     */
    void calc_specify_phase(phases phase_index){ imposed_phase_index = phase_index; _phase = phase_index; }
    /**\brief Unspecify the phase - the phase is no longer imposed, different solvers can do as they like
     */
    void calc_unspecify_phase(){ imposed_phase_index = iphase_not_imposed;}
    //
    // ************************************************************************* //
    //                         Trivial Functions                                 //
    // ************************************************************************* //
    //
    double calc_molar_mass(void);
    //
    // ************************************************************************* //
    //                      Saturation Functions                                 //
    // ************************************************************************* //
    //
    double calc_pressure(void){ return _p; };
};
} /* namespace CoolProp */
#endif /* PCSAFTBACKEND_H_ */
