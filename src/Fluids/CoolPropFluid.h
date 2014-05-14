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

class ViscosityCorrelation
{
public:
    double dilute(double T, double rhomolar);
    double residual(double T, double rhomolar);
    double critical(double T, double rhomolar);
};
class ThermalConductivityCorrelation
{
public:
    double dilute(double T, double rhomolar);
    double residual(double T, double rhomolar);
    double critical(double T, double rhomolar);
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
    std::vector<long double> a, n, s;
    long double Tc;
    std::string BibTeX;
    std::size_t N;
public:
    
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
        virtual ~CoolPropFluid(){};
        EquationOfState *pEOS; ///< A pointer to the currently used EOS
        std::vector<EquationOfState> EOSVector; ///< The equations of state that could be used for this fluid

        std::string name; ///< The name of the fluid
        std::string REFPROPname; ///< The REFPROP-compliant name if REFPROP-"name" is not a compatible fluid name.  If not included, "name" is assumed to be a valid name for REFPROP
        std::string CAS; ///< The CAS number of the fluid
        std::vector <std::string> aliases; ///< A vector of aliases of names for the fluid

        std::vector<ViscosityCorrelation*> viscosity_vector; ///< The viscosity correlations that could be used for this fluid
        std::vector<ThermalConductivityCorrelation*> thermal_conductivity_vector; ///< The thermal conductivity correlations that could be used for this fluid
        std::vector<SurfaceTensionCorrelation*> surface_tension_vector; ///< The surface tension correlations that could be used for this fluid

        BibTeXKeysStruct BibTeXKeys;
        EnvironmentalFactorsStruct environment;
        Ancillaries ancillaries;

        double gas_constant(){ return pEOS->R_u; };
        double molar_mass(){ return pEOS->molar_mass; };
};


} /* namespace CoolProp */
#endif /* COOLPROPFLUID_H_ */
