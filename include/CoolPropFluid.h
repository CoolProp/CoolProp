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
#include <assert.h>
#include <iterator>

namespace CoolProp {

struct BibTeXKeysStruct
{
    std::string EOS,
                CP0,
                VISCOSITY,
                CONDUCTIVITY,
                ECS_LENNARD_JONES,
                ECS_FITS,
                SURFACE_TENSION;
};

struct EnvironmentalFactorsStruct
{
    double GWP20, GWP100, GWP500, ODP, HH, PH, FH;
    std::string ASHRAE34;
};

/// A set of limits for the eos parameters
struct EOSLimits
{
    double Tmin, Tmax, rhomax, pmax;
};

struct ConductivityDiluteEta0AndPolyData{
    std::vector<long double> A, t;
};

struct ConductivityDiluteRatioPolynomialsData{
    long double T_reducing, p_reducing;
    std::vector<long double> A, B, n, m;
};
struct ConductivityDiluteVariables
{
    enum ConductivityDiluteEnum {CONDUCTIVITY_DILUTE_RATIO_POLYNOMIALS, 
                                 CONDUCTIVITY_DILUTE_ETA0_AND_POLY, 
                                 CONDUCTIVITY_DILUTE_CO2, 
                                 CONDUCTIVITY_DILUTE_ETHANE, 
                                 CONDUCTIVITY_DILUTE_NOT_SET
                                 };
    int type;
    ConductivityDiluteRatioPolynomialsData ratio_polynomials;
    ConductivityDiluteEta0AndPolyData eta0_and_poly;

    ConductivityDiluteVariables(){type = CONDUCTIVITY_DILUTE_NOT_SET;}
};

struct ConductivityResidualPolynomialAndExponentialData{
    long double T_reducing, rhomass_reducing;
    std::vector<long double> A, t, d, gamma, l;
};

struct ConductivityResidualPolynomialData{
    long double T_reducing, rhomass_reducing;
    std::vector<long double> B, t, d;
};
struct ConductivityResidualVariables
{
    enum ConductivityResidualEnum {CONDUCTIVITY_RESIDUAL_POLYNOMIAL,
                                   CONDUCTIVITY_RESIDUAL_POLYNOMIAL_AND_EXPONENTIAL,
                                   CONDUCTIVITY_RESIDUAL_R123,
                                   CONDUCTIVITY_RESIDUAL_CO2,
                                   CONDUCTIVITY_RESIDUAL_NOT_SET
                                   };
    int type;
    ConductivityResidualPolynomialData polynomials;
    ConductivityResidualPolynomialAndExponentialData polynomial_and_exponential;

    ConductivityResidualVariables(){type = CONDUCTIVITY_RESIDUAL_NOT_SET;}
};

struct ConductivityCriticalSimplifiedOlchowySengersData{
    long double T_reducing, p_reducing, k, R0, gamma, nu, qD, zeta0, GAMMA, T_ref;
    ConductivityCriticalSimplifiedOlchowySengersData(){
        // Universal constants - can still be adjusted if need be
        k = 1.3806488e-23; //[J/K]
        R0 = 1.03; //[-]
		gamma = 1.239; //[-]
		nu = 0.63; //[-]
        // Suggested default values - can be over-written
        GAMMA = 0.0496; //[-]
        zeta0 = 1.94e-10; //[m]
        qD = 1e9; //[m]

        // Set to invalid number, can be provided in the JSON file
        T_ref = _HUGE;
    }
};
struct ConductivityCriticalVariables
{
    enum ConductivityResidualEnum {CONDUCTIVITY_CRITICAL_SIMPLIFIED_OLCHOWY_SENGERS, 
                                   CONDUCTIVITY_CRITICAL_R123, 
                                   CONDUCTIVITY_CRITICAL_AMMONIA, 
                                   CONDUCTIVITY_CRITICAL_NONE, 
                                   CONDUCTIVITY_CRITICAL_NOT_SET
                                   };
    int type;
    ConductivityCriticalSimplifiedOlchowySengersData Olchowy_Sengers;

    ConductivityCriticalVariables(){type = CONDUCTIVITY_CRITICAL_NOT_SET; }
};

/// Variables for the dilute gas part 
struct ViscosityDiluteGasCollisionIntegralData
{
    long double molar_mass, C;
    std::vector<long double> a, t;
};
struct ViscosityDiluteCollisionIntegralPowersOfTstarData
{
    long double T_reducing, C;
    std::vector<long double> a, t;
};
struct ViscosityDiluteGasPowersOfT
{
    std::vector<long double> a, t;
};
struct ViscosityDiluteVariables
{
    enum ViscosityDiluteEnum {VISCOSITY_DILUTE_COLLISION_INTEGRAL, 
                              VISCOSITY_DILUTE_COLLISION_INTEGRAL_POWERS_OF_TSTAR, 
                              VISCOSITY_DILUTE_KINETIC_THEORY, 
                              VISCOSITY_DILUTE_ETHANE,
                              VISCOSITY_DILUTE_POWERS_OF_T,
                              VISCOSITY_DILUTE_NOT_SET
                              };
    int type;
    ViscosityDiluteGasCollisionIntegralData collision_integral;
    ViscosityDiluteCollisionIntegralPowersOfTstarData collision_integral_powers_of_Tstar;
    ViscosityDiluteGasPowersOfT powers_of_T;

    ViscosityDiluteVariables(){type = VISCOSITY_DILUTE_NOT_SET;}
};

struct ViscosityRainWaterFriendData
{
    std::vector<long double> b, t;
};
struct ViscosityInitialDensityVariables
{
    ViscosityRainWaterFriendData rainwater_friend;
};

struct ViscosityModifiedBatschinskiHildebrandData
{
    std::vector<long double> a,d1,d2,t1,t2,f,g,h,p,q,gamma, l;
    long double T_reduce, rhomolar_reduce;
};
struct ViscosityFrictionTheoryData
{
    std::vector<long double> Aa, Aaa, Aaaa, Ar, Arr, Adrdr, Arrr, Ai, Aii, AdrAdr;
    int Na, Naa, Naaa, Nr, Nrr, Nrrr, Nii;
    long double c1, c2, T_reduce, rhomolar_reduce;
};
struct ViscosityHigherOrderVariables
{
    enum ViscosityDiluteEnum {VISCOSITY_HIGHER_ORDER_BATSCHINKI_HILDEBRAND, 
                              VISCOSITY_HIGHER_ORDER_HYDROGEN,
                              VISCOSITY_HIGHER_ORDER_HEXANE,
                              VISCOSITY_HIGHER_ORDER_ETHANE,
                              VISCOSITY_HIGHER_ORDER_FRICTION_THEORY,
                              VISCOSITY_HIGHER_ORDER_NOT_SET
                              };
    int type;
    ViscosityModifiedBatschinskiHildebrandData modified_Batschinski_Hildebrand;
    ViscosityFrictionTheoryData friction_theory;
    ViscosityHigherOrderVariables(){type = VISCOSITY_HIGHER_ORDER_NOT_SET;};
};

struct ViscosityECSVariables{
    std::string reference_fluid;
    long double psi_rhomolar_reducing;
    std::vector<long double> psi_a, psi_t;
};

class TransportPropertyData
{
public:
    enum ViscosityDiluteEnum {VISCOSITY_HARDCODED_WATER, 
                              VISCOSITY_HARDCODED_HELIUM,
                              VISCOSITY_HARDCODED_R23,
                              VISCOSITY_NOT_HARDCODED
                              };
    enum ConductivityDiluteEnum {
                                 CONDUCTIVITY_HARDCODED_WATER,
                                 CONDUCTIVITY_HARDCODED_R23,
                                 CONDUCTIVITY_HARDCODED_HELIUM,
                                 CONDUCTIVITY_NOT_HARDCODED
                                 };
    ViscosityDiluteVariables viscosity_dilute;
    ViscosityInitialDensityVariables viscosity_initial;
    ViscosityHigherOrderVariables viscosity_higher_order;
    ViscosityECSVariables viscosity_ecs;
    ConductivityDiluteVariables conductivity_dilute;
    ConductivityResidualVariables conductivity_residual;
    ConductivityCriticalVariables conductivity_critical;

    std::string BibTeX_viscosity, BibTeX_conductivity;
    bool using_ECS; ///< A flag for whether to use extended corresponding states.  False for no
    long double sigma_eta, epsilon_over_k;
    int hardcoded_viscosity, hardcoded_conductivity;
    TransportPropertyData(){hardcoded_viscosity = VISCOSITY_NOT_HARDCODED; 
                            hardcoded_conductivity = CONDUCTIVITY_NOT_HARDCODED;
                            using_ECS = false;
    };
};

/**
The surface tension correlation class uses correlations for the surface tension that are all
of the form

\f[
\sigma = \sum_{i=0}^{k-1}a_i\left(1-\frac{T}{\tilde T_c}\right)^{n_i}
\f]

where \f$ \tilde T_c \f$ is the critical temperature used for the correlation which is 
almost always equal to the critical temperature of the equation of state.  Result for 
surface tension is in N/m
*/
class SurfaceTensionCorrelation
{
public:
    std::vector<long double> a, n, s;
    long double Tc;
    
    std::size_t N;

    std::string BibTeX;
    SurfaceTensionCorrelation(){};
    SurfaceTensionCorrelation(rapidjson::Value &json_code)
    {
        a = cpjson::get_long_double_array(json_code["a"]);
        n = cpjson::get_long_double_array(json_code["n"]);
        
        Tc = cpjson::get_double(json_code,"Tc");
        BibTeX = cpjson::get_string(json_code,"BibTeX");
        
        this->N = n.size();
        s = n;
    };
    long double evaluate(long double T)
    {
        if (a.empty()){ throw NotImplementedError(format("surface tension curve not provided"));}
        long double THETA = 1-T/Tc;
        for (std::size_t i = 0; i < N; ++i)
        {
            s[i] = a[i]*pow(THETA, n[i]);
        }
        return std::accumulate(s.begin(), s.end(), 0.0);
    }
};
/**
*/
class SaturationAncillaryFunction
{
private:
    std::vector<double> n, t, s;
    bool using_tau_r;
    double Tmax, Tmin, reducing_value, T_r;
    int type;
    enum ancillaryfunctiontypes{TYPE_NOT_EXPONENTIAL = 0, TYPE_EXPONENTIAL = 1};
    std::size_t N;
public:
    
    SaturationAncillaryFunction(){};
    SaturationAncillaryFunction(rapidjson::Value &json_code)
    {
        n = cpjson::get_double_array(json_code["n"]);
        t = cpjson::get_double_array(json_code["t"]);
        Tmin = cpjson::get_double(json_code,"Tmin");
        Tmax = cpjson::get_double(json_code,"Tmax");
        reducing_value = cpjson::get_double(json_code,"reducing_value");
        using_tau_r = cpjson::get_bool(json_code,"using_tau_r");
        T_r = cpjson::get_double(json_code,"T_r");
        std::string type = cpjson::get_string(json_code,"type");

        if (!type.compare("rhoLnoexp"))
            this->type = TYPE_NOT_EXPONENTIAL;
        else
            this->type = TYPE_EXPONENTIAL;
        this->N = n.size();
        s = n;
    };
    double evaluate(double T)
    {
        double THETA = 1-T/T_r;

        for (std::size_t i = 0; i < N; ++i)
        {
            s[i] = n[i]*pow(THETA, t[i]);
        }
        double summer = std::accumulate(s.begin(), s.end(), 0.0);
        
        if (type == TYPE_NOT_EXPONENTIAL)
        {    
            return reducing_value*(1+summer);
        }
        else
        {
            double tau_r_value;
            if (using_tau_r)
                tau_r_value = T_r/T;
            else
                tau_r_value = 1.0;
            return reducing_value*exp(tau_r_value*summer);
        }
    }
    double invert(double value)
    {
        // Invert the ancillary curve to get the temperature as a function of the output variable
        // Define the residual to be driven to zero
        class solver_resid : public FuncWrapper1D
        {
        public:
            int other;
            SaturationAncillaryFunction *anc;
            long double T, value, r, current_value;

            solver_resid(SaturationAncillaryFunction *anc, long double value) : anc(anc), value(value){};
                
            double call(double T){ 
                this->T = T;
                current_value = anc->evaluate(T);
                r = current_value - value;
                return r;
            };
        };
        solver_resid resid(this, value);
        std::string errstring;

        return Brent(resid,Tmin,Tmax,DBL_EPSILON,1e-12,100,errstring);
    }
};


class MeltingLine
{
};

struct Ancillaries
{
    SaturationAncillaryFunction pL, pV, rhoL, rhoV;
    MeltingLine melting_line;
    SurfaceTensionCorrelation surface_tension;
};

/// The core class for an equation of state
/** 
 This class holds the absolute minimum information to evaluate the equation 
 of state.  This includes the reducing state, limits on the equation of state,
 the coefficients for the Helmholtz derivative terms.

 It does NOT include derived parameters like specific heat, enthalpy, etc.
*/
class EquationOfState{
public:
    EquationOfState(){};
    ~EquationOfState(){};
    SimpleState reduce; ///< Reducing state used for the EOS (usually, but not always, the critical point)
    EOSLimits limits; ///< Limits on the EOS
    double R_u, ///< The universal gas constant used for this EOS (usually, but not always, 8.314472 J/mol/K)
           molar_mass, ///< The molar mass in kg/mol (note NOT kg/kmol)
           accentric, ///< The accentric factor \f$ \omega = -log_{10}\left(\frac{p_s(T/T_c=0.7)}{p_c}\right)-1\f$
           Ttriple, ///< Triple point temperature (K)
           ptriple, ///< Triple point pressure (Pa)
           rhoLtriple, ///< Density of liquid at triple point pressure (mol/m^3)
           rhoVtriple; ///< Density of vapor at triple point pressure (mol/m^3)
    bool pseudo_pure; ///< Is a pseudo-pure fluid (true) or pure fluid (false)
    ResidualHelmholtzContainer alphar; ///< The residual Helmholtz energy
    IdealHelmholtzContainer alpha0; ///< The ideal Helmholtz energy
    std::string BibTeX_EOS, ///< The bibtex key for the equation of state
                BibTeX_CP0; ///< The bibtex key for the ideal gas specific heat correlation

    /// Validate the EOS that was just constructed
    void validate()
    {
        assert(R_u < 9 && R_u > 8);
        assert(molar_mass > 0.001 && molar_mass < 1);
    };
    long double baser(const long double &tau, const long double &delta) throw()
    {
        return alphar.base(tau, delta);
    };
    // First partials
    long double dalphar_dDelta(const long double &tau, const long double &delta) throw()
    {
        return alphar.dDelta(tau, delta);
    };
    long double dalphar_dTau(const long double &tau, const long double &delta) throw()
    {
        return alphar.dTau(tau, delta);
    };
    // Second partials
    long double d2alphar_dDelta2(const long double &tau, const long double &delta) throw()
    {
        return alphar.dDelta2(tau, delta);
    };
    long double d2alphar_dDelta_dTau(const long double &tau, const long double &delta) throw()
    {
        return alphar.dDelta_dTau(tau, delta);
    };
    long double d2alphar_dTau2(const long double &tau, const long double &delta) throw()
    {
        return alphar.dTau2(tau, delta);
    };
    // Third partials
    long double d3alphar_dDelta3(const long double &tau, const long double &delta) throw()
    {
        return alphar.dDelta3(tau, delta);
    };
    long double d3alphar_dDelta2_dTau(const long double &tau, const long double &delta) throw()
    {
        return alphar.dDelta2_dTau(tau, delta);
    };
    long double d3alphar_dDelta_dTau2(const long double &tau, const long double &delta) throw()
    {
        return alphar.dDelta_dTau2(tau, delta);
    };
    long double d3alphar_dTau3(const long double &tau, const long double &delta) throw()
    {
        return alphar.dTau3(tau, delta);
    };


    long double base0(const long double &tau, const long double &delta) throw()
    {
        return alpha0.base(tau, delta);
    };
    // First partials
    long double dalpha0_dDelta(const long double &tau, const long double &delta) throw()
    {
        return alpha0.dDelta(tau, delta);
    };
    long double dalpha0_dTau(const long double &tau, const long double &delta) throw()
    {
        return alpha0.dTau(tau, delta);
    };
    // Second partials
    long double d2alpha0_dDelta2(const long double &tau, const long double &delta) throw()
    {
        return alpha0.dDelta2(tau, delta);
    };
    long double d2alpha0_dDelta_dTau(const long double &tau, const long double &delta) throw()
    {
        return alpha0.dDelta_dTau(tau, delta);
    };
    long double d2alpha0_dTau2(const long double &tau, const long double &delta) throw()
    {
        return alpha0.dTau2(tau, delta);
    };
    // Third partials
    long double d3alpha0_dDelta3(const long double &tau, const long double &delta) throw()
    {
        return alpha0.dDelta3(tau, delta);
    };
    long double d3alpha0_dDelta2_dTau(const long double &tau, const long double &delta) throw()
    {
        return alpha0.dDelta2_dTau(tau, delta);
    };
    long double d3alpha0_dDelta_dTau2(const long double &tau, const long double &delta) throw()
    {
        return alpha0.dDelta_dTau2(tau, delta);
    };
    long double d3alpha0_dTau3(const long double &tau, const long double &delta) throw()
    {
        return alpha0.dTau3(tau, delta);
    };
};

/// A thermophysical property provider for critical and reducing values as well as derivatives of Helmholtz energy
/**
This fluid instance is populated using an entry from a JSON file
*/
class CoolPropFluid {
    protected:
        // Transport property data
        std::string ECSReferenceFluid; ///< A string that gives the name of the fluids that should be used for the ECS method for transport properties
        double ECS_qd; ///< The critical qd parameter for the Olchowy-Sengers cross-over term
    public:
        CoolPropFluid(){};
        ~CoolPropFluid(){};
        EquationOfState *pEOS; ///< A pointer to the currently used EOS
        std::vector<EquationOfState> EOSVector; ///< The equations of state that could be used for this fluid

        std::string name; ///< The name of the fluid
        std::string REFPROPname; ///< The REFPROP-compliant name if REFPROP-"name" is not a compatible fluid name.  If not included, "name" is assumed to be a valid name for REFPROP
        std::string CAS; ///< The CAS number of the fluid
        std::vector <std::string> aliases; ///< A vector of aliases of names for the fluid

        /*std::vector<ViscosityCorrelation*> viscosity_vector; ///< The viscosity correlations that could be used for this fluid
        std::vector<ThermalConductivityCorrelation*> thermal_conductivity_vector; ///< The thermal conductivity correlations that could be used for this fluid
        std::vector<SurfaceTensionCorrelation*> surface_tension_vector; ///< The surface tension correlations that could be used for this fluid*/

        BibTeXKeysStruct BibTeXKeys;
        EnvironmentalFactorsStruct environment;
        Ancillaries ancillaries;
        TransportPropertyData transport;
        SimpleState crit;

        double gas_constant(){ return pEOS->R_u; };
        double molar_mass(){ return pEOS->molar_mass; };
};


} /* namespace CoolProp */
#endif /* COOLPROPFLUID_H_ */
