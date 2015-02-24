#ifndef ANCILLARIES_H
#define ANCILLARIES_H

#include "Exceptions.h"
#include <vector>
#include "rapidjson/rapidjson_include.h"
#include "Eigen/Core"
#include "PolyMath.h"

namespace CoolProp{

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
    std::vector<CoolPropDbl> a, n, s;
    CoolPropDbl Tc;

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
    CoolPropDbl evaluate(CoolPropDbl T)
    {
        if (a.empty()){ throw NotImplementedError(format("surface tension curve not provided"));}
        CoolPropDbl THETA = 1-T/Tc;
        for (std::size_t i = 0; i < N; ++i)
        {
            s[i] = a[i]*pow(THETA, n[i]);
        }
        return std::accumulate(s.begin(), s.end(), 0.0);
    }
};
/**
 * 
 * This is generalized class that can be used to manage an ancillary curve, 
 * here they are ancillary curves for saturation pressure, density, enthalpy, entropy.
 * 
 * The form of the ancillary equation can take one of a number of forms:
 * 
 * a) So-called "exponential" form (type = TYPE_EXPONENTIAL) that has a form like
 * 
 * \f[ y = y_c\exp\left(\frac{T_c}{T}\sum(n_i \theta^{t_i})\right) \f] 
 * or 
 * \f[ y = y_c\exp\left(\sum(n_i \theta^{t_i})\right) \f]
 * 
 * b) A non-exponential form (type = TYPE_NOT_EXPONENTIAL) that has a form of
 *
 * \f[ y = y_c\left(1+\sum_i(n_i\theta^t_i)\right) \f]
 * with
 * \f[ \theta = \left(1-\frac{T}{T_c}\right) \f]
 * which is conveniently equal to zero at the critical point
 * 
 * c) Rational polynomial form (type = TYPE_RATIONAL_POLYNOMIAL) that has a form of
 * \f[ y = \frac{\sum_iA_iT^i}{\sum_iB_iT^i}\f]
 * where i is an integer, and the coefficients are in increasing order in both numerator and denominator
*/
class SaturationAncillaryFunction
{
private:
    Eigen::MatrixXd num_coeffs, ///< Coefficients for numerator in rational polynomial 
                    den_coeffs; ///< Coefficients for denominator in rational polynomial
    std::vector<double> n, t, s;
    bool using_tau_r;
    CoolPropDbl Tmax, Tmin, reducing_value, T_r, max_abs_error;
    enum ancillaryfunctiontypes{TYPE_NOT_SET = 0, 
                                TYPE_NOT_EXPONENTIAL, 
                                TYPE_EXPONENTIAL, 
                                TYPE_RATIONAL_POLYNOMIAL};
    ancillaryfunctiontypes type;
    std::size_t N;
public:

    SaturationAncillaryFunction(){type = TYPE_NOT_SET;};
    SaturationAncillaryFunction(rapidjson::Value &json_code);
    
    /// Return true if the ancillary is enabled
    bool enabled(void){return type != TYPE_NOT_SET;}
    
    /// Get the maximum absolute error for this fit
    /// @returns max_abs_error the maximum absolute error for ancillaries that are characterized by maximum absolute error
    CoolPropDbl get_max_abs_error(){return max_abs_error;};
    
    /// Evaluate this ancillary function, yielding for instance the saturated liquid density
    /// @param T The temperature in K
    /// @returns y the value of the ancillary function at temperature T
    double evaluate(double T);
    
    /// Invert this ancillary function, and calculate the temperature given the output the value of the function
    /// @param value The value of the output
    /// @param min_bound (optional) The minimum value for T; ignored if < 0
    /// @param max_bound (optional) The maximum value for T; ignored if < 0
    /// @returns T The temperature in K
    double invert(double value, double min_bound = -1, double max_bound = -1);
    
    /// Get the minimum temperature
    /// @returns T The minimum temperature in K
    double get_Tmin(void){return Tmin;};
    
    /// Get the maximum temperature
    /// @returns T The maximum temperature in K
    double get_Tmax(void){return Tmax;};
    
};

struct MeltingLinePiecewiseSimonSegment
{
    CoolPropDbl T_0, a, c, p_0, T_max, T_min, p_min, p_max;
};
struct MeltingLinePiecewiseSimonData
{
    std::vector<MeltingLinePiecewiseSimonSegment> parts;
};
class MeltingLinePiecewisePolynomialInTrSegment
{
public:
    std::vector<CoolPropDbl> a, t;
    CoolPropDbl T_0, p_0, T_max, T_min, p_min, p_max;
    CoolPropDbl evaluate(CoolPropDbl T)
    {
        CoolPropDbl summer = 0;
        for (std::size_t i =0; i < a.size(); ++i){
            summer += a[i]*(pow(T/T_0,t[i])-1);
        }
        return p_0*(1+summer);
    }
};
struct MeltingLinePiecewisePolynomialInTrData
{
    std::vector<MeltingLinePiecewisePolynomialInTrSegment> parts;
};
class MeltingLinePiecewisePolynomialInThetaSegment
{
public:
    std::vector<CoolPropDbl> a, t;
    CoolPropDbl T_0, p_0, T_max, T_min, p_min, p_max;
    
    CoolPropDbl evaluate(CoolPropDbl T)
    {
        CoolPropDbl summer = 0;
        for (std::size_t i =0; i < a.size(); ++i){
            summer += a[i]*pow(T/T_0-1,t[i]);
        }
        return p_0*(1+summer);
    }
};
struct MeltingLinePiecewisePolynomialInThetaData
{
    std::vector<MeltingLinePiecewisePolynomialInThetaSegment> parts;
};
class MeltingLineVariables
{
public:
    enum MeltingLineVariablesEnum{
        MELTING_LINE_SIMON_TYPE,
        MELTING_LINE_POLYNOMIAL_IN_TR_TYPE,
        MELTING_LINE_POLYNOMIAL_IN_THETA_TYPE,
        MELTING_LINE_NOT_SET
    };
    CoolPropDbl Tmin, Tmax, pmin, pmax;
    
    CoolPropDbl evaluate(int OF, int GIVEN, CoolPropDbl value);
    
    /// Evaluate the melting line to calculate the limits of the curve (Tmin/Tmax and pmin/pmax)
    void set_limits();
    
    bool enabled(){return type != MELTING_LINE_NOT_SET;};
    std::string BibTeX;
    CoolPropDbl T_m; ///< Melting temperature at 1 atmosphere
    MeltingLinePiecewiseSimonData simon;
    MeltingLinePiecewisePolynomialInTrData polynomial_in_Tr;
    MeltingLinePiecewisePolynomialInThetaData polynomial_in_Theta;
    int type;
    MeltingLineVariables(){type = MELTING_LINE_NOT_SET;};
};

} /* namespace CoolProp */
#endif