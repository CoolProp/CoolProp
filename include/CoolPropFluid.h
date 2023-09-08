/*
 * CoolPropFluid.h
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#ifndef COOLPROPFLUID_H_
#define COOLPROPFLUID_H_

#include "DataStructures.h"
#include "Helmholtz.h"
#include "Solvers.h"

#include <numeric>
#include <string>
#include <vector>
#include <map>
#include <cassert>
#include <iterator>
#include "Eigen/Core"
#include "PolyMath.h"
#include "Ancillaries.h"

namespace CoolProp {

struct BibTeXKeysStruct
{
    std::string EOS, CP0, VISCOSITY, CONDUCTIVITY, ECS_LENNARD_JONES, ECS_FITS, SURFACE_TENSION;
};

struct EnvironmentalFactorsStruct
{
    double GWP20, GWP100, GWP500, ODP, HH, PH, FH;
    std::string ASHRAE34;
};
struct CriticalRegionSplines
{
    double T_min, T_max, rhomolar_min, rhomolar_max;
    std::vector<double> cL, cV;
    bool enabled;
    CriticalRegionSplines() : T_min(_HUGE), T_max(_HUGE), rhomolar_min(_HUGE), rhomolar_max(_HUGE), enabled(false){};

    const void get_densities(double T, double rho_min, double rho_crit, double rho_max, double& rhoL, double& rhoV) const {
        int Nsoln = -1, Ngood = 0;
        double rho1 = 0, rho2 = 0, rho3 = 0;

        // -----------
        // Liquid part
        // -----------
        Ngood = 0;
        solve_cubic(cL[0], cL[1], cL[2], cL[3] - T, Nsoln, rho1, rho2, rho3);
        if (Nsoln == 1 && rho1 < rho_max && rho1 > rho_crit) {
            rhoL = rho1;
        } else {
            if (rho1 < rho_max && rho1 > rho_crit) {
                Ngood++;
                rhoL = rho1;
            }
            if (rho2 < rho_max && rho2 > rho_crit) {
                Ngood++;
                rhoL = rho2;
            }
            if (Nsoln > 2 && rho3 < rho_max && rho3 > rho_crit) {
                Ngood++;
                rhoL = rho3;
            }
            if (Ngood > 1) {
                throw ValueError(format("More than one liquid solution found for critical spline for T=%0.12g", T));
            };
            if (Ngood < 1) {
                throw ValueError(format("No liquid solution found for critical spline for T=%0.12g", T));
            };
        }

        // ----------
        // Vapor part
        // ----------
        Ngood = 0;
        Nsoln = 0;
        solve_cubic(cV[0], cV[1], cV[2], cV[3] - T, Nsoln, rho1, rho2, rho3);
        if (Nsoln == 1 && rho1 > rho_min && rho1 < rho_crit) {
            rhoV = rho1;
        } else {
            if (rho1 > rho_min && rho1 < rho_crit) {
                Ngood++;
                rhoV = rho1;
            }
            if (rho2 > rho_min && rho2 < rho_crit) {
                Ngood++;
                rhoV = rho2;
            }
            if (Nsoln > 2 && rho3 > rho_min && rho3 < rho_crit) {
                Ngood++;
                rhoV = rho3;
            }
            if (Ngood > 1) {
                throw ValueError(format("More than one vapor solution found for critical spline for T=%0.12g", T));
            };
            if (Ngood < 1) {
                throw ValueError(format("No vapor solution found for critical spline for T=%0.12g", T));
            };
        }
    };
};

/// A set of limits for the eos parameters
struct EOSLimits
{
    double Tmin, Tmax, rhomax, pmax;
};

struct ConductivityECSVariables
{
    std::string reference_fluid;
    CoolPropDbl psi_rhomolar_reducing, f_int_T_reducing;
    std::vector<CoolPropDbl> psi_a, psi_t, f_int_a, f_int_t;
};

struct ConductivityDiluteEta0AndPolyData
{
    std::vector<CoolPropDbl> A, t;
};

struct ConductivityDiluteRatioPolynomialsData
{
    CoolPropDbl T_reducing, p_reducing;
    std::vector<CoolPropDbl> A, B, n, m;
};
struct ConductivityDiluteVariables
{
    enum ConductivityDiluteEnum
    {
        CONDUCTIVITY_DILUTE_RATIO_POLYNOMIALS,
        CONDUCTIVITY_DILUTE_ETA0_AND_POLY,
        CONDUCTIVITY_DILUTE_CO2,
        CONDUCTIVITY_DILUTE_CO2_HUBER_JPCRD_2016,
        CONDUCTIVITY_DILUTE_ETHANE,
        CONDUCTIVITY_DILUTE_NONE,
        CONDUCTIVITY_DILUTE_NOT_SET
    };
    int type;
    ConductivityDiluteRatioPolynomialsData ratio_polynomials;
    ConductivityDiluteEta0AndPolyData eta0_and_poly;

    ConductivityDiluteVariables() {
        type = CONDUCTIVITY_DILUTE_NOT_SET;
    }
};

struct ConductivityResidualPolynomialAndExponentialData
{
    CoolPropDbl T_reducing, rhomass_reducing;
    std::vector<CoolPropDbl> A, t, d, gamma, l;
};

struct ConductivityResidualPolynomialData
{
    CoolPropDbl T_reducing, rhomass_reducing;
    std::vector<CoolPropDbl> B, t, d;
};
struct ConductivityResidualVariables
{
    enum ConductivityResidualEnum
    {
        CONDUCTIVITY_RESIDUAL_POLYNOMIAL,
        CONDUCTIVITY_RESIDUAL_POLYNOMIAL_AND_EXPONENTIAL,
        CONDUCTIVITY_RESIDUAL_R123,
        CONDUCTIVITY_RESIDUAL_CO2,
        CONDUCTIVITY_RESIDUAL_NOT_SET
    };
    int type;
    ConductivityResidualPolynomialData polynomials;
    ConductivityResidualPolynomialAndExponentialData polynomial_and_exponential;

    ConductivityResidualVariables() {
        type = CONDUCTIVITY_RESIDUAL_NOT_SET;
    }
};

struct ConductivityCriticalSimplifiedOlchowySengersData
{
    CoolPropDbl k, R0, gamma, nu, GAMMA, zeta0, qD, T_reducing, p_reducing, T_ref;
    ConductivityCriticalSimplifiedOlchowySengersData()
      :                    // Universal constants - can still be adjusted if need be
        k(1.3806488e-23),  //[J/K]
        R0(1.03),          //[-]
        gamma(1.239),      //[-]
        nu(0.63),          //[-]
        // Suggested default values - can be over-written
        GAMMA(0.0496),    //[-]
        zeta0(1.94e-10),  //[m]
        qD(2e9),          //[m]
        // Set to invalid number, can be provided in the JSON file
        // T_ref default is 1.5*Tc
        T_reducing(_HUGE),
        p_reducing(_HUGE),
        T_ref(_HUGE) {}
};
struct ConductivityCriticalVariables
{
    enum ConductivityResidualEnum
    {
        CONDUCTIVITY_CRITICAL_SIMPLIFIED_OLCHOWY_SENGERS,
        CONDUCTIVITY_CRITICAL_R123,
        CONDUCTIVITY_CRITICAL_AMMONIA,
        CONDUCTIVITY_CRITICAL_NONE,
        CONDUCTIVITY_CRITICAL_CARBONDIOXIDE_SCALABRIN_JPCRD_2006,
        CONDUCTIVITY_CRITICAL_NOT_SET
    };
    int type;
    ConductivityCriticalSimplifiedOlchowySengersData Olchowy_Sengers;

    ConductivityCriticalVariables() {
        type = CONDUCTIVITY_CRITICAL_NOT_SET;
    }
};

/// Variables for the dilute gas part
struct ViscosityDiluteGasCollisionIntegralData
{
    CoolPropDbl molar_mass, C;
    std::vector<CoolPropDbl> a, t;
};
struct ViscosityDiluteCollisionIntegralPowersOfTstarData
{
    CoolPropDbl T_reducing,  ///< Reducing temperature [K[
      C;                     ///< Leading constant
    std::vector<CoolPropDbl> a, t;
};
struct ViscosityDiluteGasPowersOfT
{
    std::vector<CoolPropDbl> a, t;
};
struct ViscosityDiluteGasPowersOfTr
{
    std::vector<CoolPropDbl> a, t;
    CoolPropDbl T_reducing;
};
struct ViscosityDiluteVariables
{
    enum ViscosityDiluteType
    {
        VISCOSITY_DILUTE_COLLISION_INTEGRAL,                  ///< Use \ref TransportRoutines::viscosity_dilute_collision_integral
        VISCOSITY_DILUTE_COLLISION_INTEGRAL_POWERS_OF_TSTAR,  ///< Use \ref TransportRoutines::viscosity_dilute_collision_integral_powers_of_T
        VISCOSITY_DILUTE_KINETIC_THEORY,                      ///< Use \ref TransportRoutines::viscosity_dilute_kinetic_theory
        VISCOSITY_DILUTE_ETHANE,                              ///< Use \ref TransportRoutines::viscosity_dilute_ethane
        VISCOSITY_DILUTE_CYCLOHEXANE,                         ///< Use \ref TransportRoutines::viscosity_dilute_cyclohexane
        VISCOSITY_DILUTE_CO2_LAESECKE_JPCRD_2017,             ///< Use \ref TransportRoutines::viscosity_dilute_CO2_LaeseckeJPCRD2017
        VISCOSITY_DILUTE_POWERS_OF_T,                         ///< Use \ref TransportRoutines::viscosity_dilute_powers_of_T
        VISCOSITY_DILUTE_POWERS_OF_TR,                        ///< Use \ref TransportRoutines::viscosity_dilute_powers_of_Tr
        VISCOSITY_DILUTE_NOT_SET
    };
    ViscosityDiluteType type;
    ViscosityDiluteGasCollisionIntegralData collision_integral;  ///< Data for \ref TransportRoutines::viscosity_dilute_collision_integral
    ViscosityDiluteCollisionIntegralPowersOfTstarData
      collision_integral_powers_of_Tstar;       ///< Data for \ref TransportRoutines::viscosity_dilute_collision_integral_powers_of_T
    ViscosityDiluteGasPowersOfT powers_of_T;    ///< Data for \ref TransportRoutines::viscosity_dilute_powers_of_T
    ViscosityDiluteGasPowersOfTr powers_of_Tr;  ///< Data for \ref TransportRoutines::viscosity_dilute_powers_of_Tr
    ViscosityDiluteVariables() {
        type = VISCOSITY_DILUTE_NOT_SET;
    }
};

struct ViscosityRainWaterFriendData
{
    std::vector<CoolPropDbl> b, t;
};
struct ViscosityInitialDensityEmpiricalData
{
    std::vector<CoolPropDbl> n, d, t;
    CoolPropDbl T_reducing, rhomolar_reducing;
};

struct ViscosityInitialDensityVariables
{
    enum ViscosityInitialDensityEnum
    {
        VISCOSITY_INITIAL_DENSITY_RAINWATER_FRIEND,  ///< Use \ref TransportRoutines::viscosity_initial_density_dependence_Rainwater_Friend
        VISCOSITY_INITIAL_DENSITY_EMPIRICAL,         ///< Use \ref TransportRoutines::viscosity_initial_density_dependence_empirical
        VISCOSITY_INITIAL_DENSITY_NOT_SET
    };
    ViscosityInitialDensityEnum type;
    ViscosityRainWaterFriendData rainwater_friend;   ///< Data for \ref TransportRoutines::viscosity_initial_density_dependence_Rainwater_Friend
    ViscosityInitialDensityEmpiricalData empirical;  ///< Data for \ref TransportRoutines::viscosity_initial_density_dependence_empirical
    ViscosityInitialDensityVariables() {
        type = VISCOSITY_INITIAL_DENSITY_NOT_SET;
    }
};

struct ViscosityModifiedBatschinskiHildebrandData
{
    std::vector<CoolPropDbl> a, d1, d2, t1, t2, f, g, h, p, q, gamma, l;
    CoolPropDbl T_reduce, rhomolar_reduce;
};
struct ViscosityFrictionTheoryData
{
    std::vector<CoolPropDbl> Aa, Aaa, Aaaa, Ar, Arr, Adrdr, Arrr, Ai, Aii, AdrAdr;
    int Na, Naa, Naaa, Nr, Nrr, Nrrr, Nii;
    CoolPropDbl c1, c2, T_reduce, rhomolar_reduce;
};
struct ViscosityHigherOrderVariables
{
    enum ViscosityHigherOrderEnum
    {
        VISCOSITY_HIGHER_ORDER_BATSCHINKI_HILDEBRAND,  ///< Use \ref TransportRoutines::viscosity_higher_order_modified_Batschinski_Hildebrand
        VISCOSITY_HIGHER_ORDER_HYDROGEN,               ///< Use \ref TransportRoutines::viscosity_hydrogen_higher_order_hardcoded
        VISCOSITY_HIGHER_ORDER_HEXANE,                 ///< Use \ref TransportRoutines::viscosity_hexane_higher_order_hardcoded
        VISCOSITY_HIGHER_ORDER_HEPTANE,                ///< Use \ref TransportRoutines::viscosity_heptane_higher_order_hardcoded
        VISCOSITY_HIGHER_ORDER_ETHANE,                 ///< Use \ref TransportRoutines::viscosity_ethane_higher_order_hardcoded
        VISCOSITY_HIGHER_ORDER_BENZENE,                ///< Use \ref TransportRoutines::viscosity_benzene_higher_order_hardcoded
        VISCOSITY_HIGHER_ORDER_TOLUENE,                ///< Use \ref TransportRoutines::viscosity_toluene_higher_order_hardcoded
        VISCOSITY_HIGHER_ORDER_CO2_LAESECKE_JPCRD_2017,///< Use \ref TransportRoutines::viscosity_CO2_higher_order_hardcoded_LaeseckeJPCRD2017
        VISCOSITY_HIGHER_ORDER_FRICTION_THEORY,        ///< Use \ref TransportRoutines::viscosity_higher_order_friction_theory
        VISCOSITY_HIGHER_ORDER_NOT_SET
    };
    ViscosityHigherOrderEnum type;
    ViscosityModifiedBatschinskiHildebrandData
      modified_Batschinski_Hildebrand;            ///< Data for \ref TransportRoutines::viscosity_higher_order_modified_Batschinski_Hildebrand
    ViscosityFrictionTheoryData friction_theory;  ///< Data for \ref TransportRoutines::viscosity_higher_order_friction_theory
    ViscosityHigherOrderVariables() {
        type = VISCOSITY_HIGHER_ORDER_NOT_SET;
    };
};

struct ViscosityRhoSrVariables
{
    std::vector<double> c_liq, c_vap;
    double C, x_crossover, rhosr_critical;
};
struct ViscosityECSVariables
{
    std::string reference_fluid;
    CoolPropDbl psi_rhomolar_reducing;
    std::vector<CoolPropDbl> psi_a, psi_t;
};
struct ViscosityChungData
{
    CoolPropDbl rhomolar_critical, acentric, T_critical, molar_mass, dipole_moment_D;
};

class TransportPropertyData
{
   public:
    enum ViscosityHardcodedEnum
    {
        VISCOSITY_HARDCODED_WATER,       ///< Use \ref TransportRoutines::viscosity_water_hardcoded
        VISCOSITY_HARDCODED_HEAVYWATER,  ///< Use \ref TransportRoutines::viscosity_heavywater_hardcoded
        VISCOSITY_HARDCODED_HELIUM,      ///< Use \ref TransportRoutines::viscosity_helium_hardcoded
        VISCOSITY_HARDCODED_R23,         ///< Use \ref TransportRoutines::viscosity_R23_hardcoded
        VISCOSITY_HARDCODED_METHANOL,    ///< Use \ref TransportRoutines::viscosity_methanol_hardcoded
        VISCOSITY_HARDCODED_M_XYLENE,    ///< Use \ref TransportRoutines::viscosity_m_xylene_hardcoded
        VISCOSITY_HARDCODED_O_XYLENE,    ///< Use \ref TransportRoutines::viscosity_o_xylene_hardcoded
        VISCOSITY_HARDCODED_P_XYLENE,    ///< Use \ref TransportRoutines::viscosity_p_xylene_hardcoded
        VISCOSITY_NOT_HARDCODED
    };
    enum ConductivityHardcodedEnum
    {
        CONDUCTIVITY_HARDCODED_WATER,       ///< Use \ref TransportRoutines::conductivity_hardcoded_water
        CONDUCTIVITY_HARDCODED_HEAVYWATER,  ///< Use \ref TransportRoutines::conductivity_hardcoded_heavywater
        CONDUCTIVITY_HARDCODED_R23,         ///< Use \ref TransportRoutines::conductivity_hardcoded_R23
        CONDUCTIVITY_HARDCODED_HELIUM,      ///< Use \ref TransportRoutines::conductivity_hardcoded_helium
        CONDUCTIVITY_HARDCODED_METHANE,     ///< Use \ref TransportRoutines::conductivity_hardcoded_methane
        CONDUCTIVITY_NOT_HARDCODED
    };
    ViscosityDiluteVariables viscosity_dilute;
    ViscosityInitialDensityVariables viscosity_initial;
    ViscosityHigherOrderVariables viscosity_higher_order;
    ViscosityECSVariables viscosity_ecs;
    ViscosityRhoSrVariables viscosity_rhosr;
    ViscosityChungData viscosity_Chung;

    ConductivityDiluteVariables conductivity_dilute;
    ConductivityResidualVariables conductivity_residual;
    ConductivityCriticalVariables conductivity_critical;
    ConductivityECSVariables conductivity_ecs;

    std::string BibTeX_viscosity,                      ///< The BibTeX key for the viscosity model
      BibTeX_conductivity;                             ///< The BibTeX key for the conductivity model
    bool viscosity_using_ECS;                          ///< A flag for whether to use extended corresponding states for viscosity.  False for no
    bool conductivity_using_ECS;                       ///< A flag for whether to use extended corresponding states for conductivity.  False for no
    bool viscosity_using_Chung;                        ///< A flag for whether to use Chung model. False for no
    bool viscosity_using_rhosr;                        ///< A flag for whether to use rho*sr CS model of Bell. False for no
    bool viscosity_model_provided;                     ///< A flag for whether viscosity model is provided.  False for no
    bool conductivity_model_provided;                  ///< A flag for whether thermal conductivity model is provided.  False for no
    CoolPropDbl sigma_eta,                             ///< The Lennard-Jones 12-6 \f$ \sigma \f$ parameter
      epsilon_over_k;                                  ///< The Lennard-Jones 12-6 \f$ \varepsilon/k \f$ parameter
    ViscosityHardcodedEnum hardcoded_viscosity;        ///< Hardcoded flags for the viscosity
    ConductivityHardcodedEnum hardcoded_conductivity;  ///< Hardcoded flags for the conductivity
    TransportPropertyData()
      : viscosity_using_ECS(false),
        conductivity_using_ECS(false),
        viscosity_using_Chung(false),
        viscosity_using_rhosr(false),
        viscosity_model_provided(false),
        conductivity_model_provided(false),
        sigma_eta(_HUGE),
        epsilon_over_k(_HUGE),
        hardcoded_viscosity(VISCOSITY_NOT_HARDCODED),
        hardcoded_conductivity(CONDUCTIVITY_NOT_HARDCODED) {}
};

struct Ancillaries
{
    SaturationAncillaryFunction pL, pV, rhoL, rhoV, hL, hLV, sL, sLV;
    MeltingLineVariables melting_line;
    SurfaceTensionCorrelation surface_tension;
};

/// The core class for an equation of state
/**
 This class holds the absolute minimum information to evaluate the equation
 of state.  This includes the reducing state, limits on the equation of state,
 the coefficients for the Helmholtz derivative terms.

 It does NOT include derived parameters like specific heat, enthalpy, etc.
*/
class EquationOfState
{
   public:
    EquationOfState(){};
    ~EquationOfState(){};
    SimpleState reduce,  ///< Reducing state used for the EOS (usually, but not always, the critical point)
      sat_min_liquid,    ///< The saturated liquid state at the minimum saturation temperature
      sat_min_vapor,     ///< The saturated vapor state at the minimum saturation temperature
      hs_anchor,         ///< A fixed anchor state at Tc*1.1 and rhoc*0.9 used as a reference state for enthalpy and entropy ancillary curves
      max_sat_T,         ///< The state at the maximum saturation temperature for pseudo-pure
      max_sat_p;         ///< The state at the maximum saturation pressure for pseudo-pure
    EOSLimits limits;    ///< Limits on the EOS
    double R_u,          ///< The universal gas constant used for this EOS (usually, but not always, 8.314472 J/mol/K)
      molar_mass,        ///< The molar mass in kg/mol (note NOT kg/kmol)
      acentric,          ///< The acentric factor \f$ \omega = -log_{10}\left(\frac{p_s(T/T_c=0.7)}{p_c}\right)-1\f$
      Ttriple,           ///< Triple point temperature (K)
      ptriple;           ///< Triple point pressure (Pa)
    bool pseudo_pure;    ///< Is a pseudo-pure fluid (true) or pure fluid (false)
    ResidualHelmholtzContainer alphar;  ///< The residual Helmholtz energy
    IdealHelmholtzContainer alpha0;     ///< The ideal Helmholtz energy
    std::string BibTeX_EOS,             ///< The bibtex key for the equation of state
      BibTeX_CP0;                       ///< The bibtex key for the ideal gas specific heat correlation
    CriticalRegionSplines
      critical_region_splines;  ///< A cubic spline in the form T = f(rho) for saturated liquid and saturated vapor curves in the near-critical region

    /// Validate the EOS that was just constructed
    void validate() {
        assert(R_u < 9 && R_u > 8);
        assert(molar_mass > 0.001 && molar_mass < 1);
    };
    CoolPropDbl baser(const CoolPropDbl& tau, const CoolPropDbl& delta) {
        return alphar.base(tau, delta);
    };
    // First partials
    CoolPropDbl dalphar_dDelta(const CoolPropDbl& tau, const CoolPropDbl& delta) {
        return alphar.dDelta(tau, delta);
    };
    CoolPropDbl dalphar_dTau(const CoolPropDbl& tau, const CoolPropDbl& delta) {
        return alphar.dTau(tau, delta);
    };
    // Second partials
    CoolPropDbl d2alphar_dDelta2(const CoolPropDbl& tau, const CoolPropDbl& delta) {
        return alphar.dDelta2(tau, delta);
    };
    CoolPropDbl d2alphar_dDelta_dTau(const CoolPropDbl& tau, const CoolPropDbl& delta) {
        return alphar.dDelta_dTau(tau, delta);
    };
    CoolPropDbl d2alphar_dTau2(const CoolPropDbl& tau, const CoolPropDbl& delta) {
        return alphar.dTau2(tau, delta);
    };
    // Third partials
    CoolPropDbl d3alphar_dDelta3(const CoolPropDbl& tau, const CoolPropDbl& delta) {
        return alphar.dDelta3(tau, delta);
    };
    CoolPropDbl d3alphar_dDelta2_dTau(const CoolPropDbl& tau, const CoolPropDbl& delta) {
        return alphar.dDelta2_dTau(tau, delta);
    };
    CoolPropDbl d3alphar_dDelta_dTau2(const CoolPropDbl& tau, const CoolPropDbl& delta) {
        return alphar.dDelta_dTau2(tau, delta);
    };
    CoolPropDbl d3alphar_dTau3(const CoolPropDbl& tau, const CoolPropDbl& delta) {
        return alphar.dTau3(tau, delta);
    };

    CoolPropDbl base0(const CoolPropDbl& tau, const CoolPropDbl& delta) {
        return alpha0.base(tau, delta);
    };
    // First partials
    CoolPropDbl dalpha0_dDelta(const CoolPropDbl& tau, const CoolPropDbl& delta) {
        return alpha0.dDelta(tau, delta);
    };
    CoolPropDbl dalpha0_dTau(const CoolPropDbl& tau, const CoolPropDbl& delta) {
        return alpha0.dTau(tau, delta);
    };
    // Second partials
    CoolPropDbl d2alpha0_dDelta2(const CoolPropDbl& tau, const CoolPropDbl& delta) {
        return alpha0.dDelta2(tau, delta);
    };
    CoolPropDbl d2alpha0_dDelta_dTau(const CoolPropDbl& tau, const CoolPropDbl& delta) {
        return alpha0.dDelta_dTau(tau, delta);
    };
    CoolPropDbl d2alpha0_dTau2(const CoolPropDbl& tau, const CoolPropDbl& delta) {
        return alpha0.dTau2(tau, delta);
    };
    // Third partials
    CoolPropDbl d3alpha0_dDelta3(const CoolPropDbl& tau, const CoolPropDbl& delta) {
        return alpha0.dDelta3(tau, delta);
    };
    CoolPropDbl d3alpha0_dDelta2_dTau(const CoolPropDbl& tau, const CoolPropDbl& delta) {
        return alpha0.dDelta2_dTau(tau, delta);
    };
    CoolPropDbl d3alpha0_dDelta_dTau2(const CoolPropDbl& tau, const CoolPropDbl& delta) {
        return alpha0.dDelta_dTau2(tau, delta);
    };
    CoolPropDbl d3alpha0_dTau3(const CoolPropDbl& tau, const CoolPropDbl& delta) {
        return alpha0.dTau3(tau, delta);
    };
};

/// A thermophysical property provider for critical and reducing values as well as derivatives of Helmholtz energy
/**
This fluid instance is populated using an entry from a JSON file
*/
class CoolPropFluid
{
   protected:
    // Transport property data
    std::string ECSReferenceFluid;  ///< A string that gives the name of the fluids that should be used for the ECS method for transport properties
    double ECS_qd;                  ///< The critical qd parameter for the Olchowy-Sengers cross-over term
   public:
    CoolPropFluid() : ECS_qd(-_HUGE) {
        this->ChemSpider_id = -1;
    };
    ~CoolPropFluid(){};
    const EquationOfState& EOS() const {
        return EOSVector[0];
    }  ///< Get a reference to the equation of state
    EquationOfState& EOS() {
        return EOSVector[0];
    }                                        ///< Get a reference to the equation of state
    std::vector<EquationOfState> EOSVector;  ///< The equations of state that could be used for this fluid

    std::string name;  ///< The name of the fluid
    std::string
      REFPROPname;  ///< The REFPROP-compliant name if REFPROP-"name" is not a compatible fluid name.  If not included, "name" is assumed to be a valid name for REFPROP
    std::string CAS;                   ///< The CAS number of the fluid
    std::string formula;               ///< The chemical formula, in LaTeX form
    std::vector<std::string> aliases;  ///< A vector of aliases of names for the fluid
    std::string InChI;                 ///< The InChI string for the fluid
    std::string InChIKey;              ///< The InChI key for the fluid
    std::string smiles;                ///< The SMILES identifier for the fluid
    int ChemSpider_id;                 ///< The ChemSpider identifier for the fluid
    std::string TwoDPNG_URL;           ///< The URL to a 2D representation of the molecule (from ChemSpider)

    BibTeXKeysStruct BibTeXKeys;             ///< The BibTeX keys associated
    EnvironmentalFactorsStruct environment;  ///< The environmental variables for global warming potential, ODP, etc.
    Ancillaries ancillaries;                 ///< The set of ancillary equations for dewpoint, bubblepoint, surface tension, etc.
    TransportPropertyData transport;
    SimpleState crit,  ///< The state at the critical point
      triple_liquid,   ///< The saturated liquid state at the triple point temperature
      triple_vapor;    ///< The saturated vapor state at the triple point temperature

    double gas_constant() {
        return EOS().R_u;
    };
    double molar_mass() {
        return EOS().molar_mass;
    };
};

#if !defined(NO_FMTLIB) && FMT_VERSION >= 90000
inline int format_as(ViscosityDiluteVariables::ViscosityDiluteType type) {
    return fmt::underlying(type);
}

inline int format_as(TransportPropertyData::ViscosityHardcodedEnum viscosity) {
    return fmt::underlying(viscosity);
}

inline int format_as(TransportPropertyData::ConductivityHardcodedEnum conductivity) {
    return fmt::underlying(conductivity);
}
#endif

} /* namespace CoolProp */
#endif /* COOLPROPFLUID_H_ */
