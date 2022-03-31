#ifndef ANCILLARIES_H
#define ANCILLARIES_H

#include "Exceptions.h"
#include <vector>
#include "rapidjson_include.h"
#include "Eigen/Core"
#include "PolyMath.h"

namespace CoolProp {

/**
The surface tension correlation class uses correlations for the surface tension that are all
of the form

\f[
\sigma = \sum_i a_i\left(1-\frac{T}{\tilde T_c}\right)^{n_i}
\f]

where \f$ \tilde T_c \f$ is the critical temperature used for the correlation which is
almost always equal to the critical temperature of the equation of state.  Result for
surface tension is in N/m
*/
class SurfaceTensionCorrelation
{
   public:
    std::vector<CoolPropDbl> a,  ///< the leading coefficients a_i
      n,                         ///< the powers n_i
      s;                         ///< a summation buffer
    CoolPropDbl Tc;              ///< critical temperature in K
    std::size_t N;               ///< number of a_i, n_i pairs
    std::string BibTeX;          ///< The BiBTeX key for the surface tension curve in use

    SurfaceTensionCorrelation() : Tc(_HUGE), N(0) {}
    SurfaceTensionCorrelation(rapidjson::Value& json_code) {
        a = cpjson::get_long_double_array(json_code["a"]);
        n = cpjson::get_long_double_array(json_code["n"]);

        Tc = cpjson::get_double(json_code, "Tc");
        BibTeX = cpjson::get_string(json_code, "BibTeX");

        this->N = n.size();
        s = n;
    };
    /// Actually evaluate the surface tension equation
    CoolPropDbl evaluate(CoolPropDbl T) {
        if (a.empty()) {
            throw NotImplementedError(format("surface tension curve not provided"));
        }
        if (T > Tc) {
            throw ValueError(format("Must be saturated state : T <= Tc"));
        }
        CoolPropDbl THETA = 1 - T / Tc;
        for (std::size_t i = 0; i < N; ++i) {
            s[i] = a[i] * pow(THETA, n[i]);
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
    Eigen::MatrixXd num_coeffs,   ///< Coefficients for numerator in rational polynomial
      den_coeffs;                 ///< Coefficients for denominator in rational polynomial
    std::vector<double> n, t, s;  // For TYPE_NOT_EXPONENTIAL & TYPE_EXPONENTIAL
    union
    {
        CoolPropDbl max_abs_error;  ///< For TYPE_RATIONAL_POLYNOMIAL
        struct
        {                                // For TYPE_NOT_EXPONENTIAL & TYPE_EXPONENTIAL
            bool using_tau_r;            ///< Whether the term \f$ \frac{T_c}{T} \f$ is included in the
            CoolPropDbl reducing_value,  ///< The value used to reduce the output variable
              T_r;                       ///< The temperature in K used to reduce the temperature (usually the critical temperature)
            std::size_t N;               ///< The number of values in the arrays
        };
    };
    CoolPropDbl Tmax,  ///< The maximum temperature in K
      Tmin;            ///< The minimum temperature in K
    enum ancillaryfunctiontypes
    {
        TYPE_NOT_SET = 0,
        TYPE_NOT_EXPONENTIAL,     ///< It is a non-exponential type of equation
        TYPE_EXPONENTIAL,         ///< It is an exponential type equation, with or without the T_c/T term
        TYPE_RATIONAL_POLYNOMIAL  ///< It is a rational polynomial equation
    };
    ancillaryfunctiontypes type;  ///< The type of ancillary curve being used
   public:
    SaturationAncillaryFunction() {
        type = TYPE_NOT_SET;
        Tmin = _HUGE;
        Tmax = _HUGE;
    };
    SaturationAncillaryFunction(rapidjson::Value& json_code);

    /// Return true if the ancillary is enabled (type is not TYPE_NOT_SET)
    bool enabled(void) {
        return type != TYPE_NOT_SET;
    }

    /// Get the maximum absolute error for this fit
    /// @returns max_abs_error the maximum absolute error for ancillaries that are characterized by maximum absolute error
    CoolPropDbl get_max_abs_error() {
        return max_abs_error;
    };

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

    /// Get the minimum temperature in K
    double get_Tmin(void) {
        return Tmin;
    };

    /// Get the maximum temperature in K
    double get_Tmax(void) {
        return Tmax;
    };
};

// ****************************************************************************
// ****************************************************************************
//                                 MELTING LINE
// ****************************************************************************
// ****************************************************************************

struct MeltingLinePiecewiseSimonSegment
{
    CoolPropDbl T_0, a, c, p_0, T_max, T_min, p_min, p_max;
};
struct MeltingLinePiecewiseSimonData
{
    std::vector<MeltingLinePiecewiseSimonSegment> parts;
};

/**
\brief The evaluator class for a melting curve formed of segments in the form

\f[
a_i((\frac{T}{T_0})^{t_i}-1)
\f]
*/
class MeltingLinePiecewisePolynomialInTrSegment
{
   public:
    std::vector<CoolPropDbl> a, t;
    CoolPropDbl T_0, p_0, T_max, T_min, p_min, p_max;
    CoolPropDbl evaluate(CoolPropDbl T) {
        CoolPropDbl summer = 0;
        for (std::size_t i = 0; i < a.size(); ++i) {
            summer += a[i] * (pow(T / T_0, t[i]) - 1);
        }
        return p_0 * (1 + summer);
    }
};
struct MeltingLinePiecewisePolynomialInTrData
{
    std::vector<MeltingLinePiecewisePolynomialInTrSegment> parts;
};

/**
 \brief The evaluator class for a melting curve formed of segments in the form

 \f[
 a_i(\frac{T}{T_0}-1)^t_i
 \f]
*/
class MeltingLinePiecewisePolynomialInThetaSegment
{
   public:
    std::vector<CoolPropDbl> a, t;
    CoolPropDbl T_0, p_0, T_max, T_min, p_min, p_max;

    CoolPropDbl evaluate(CoolPropDbl T) {
        CoolPropDbl summer = 0;
        for (std::size_t i = 0; i < a.size(); ++i) {
            summer += a[i] * pow(T / T_0 - 1, t[i]);
        }
        return p_0 * (1 + summer);
    }
};
struct MeltingLinePiecewisePolynomialInThetaData
{
    std::vector<MeltingLinePiecewisePolynomialInThetaSegment> parts;
};

class MeltingLineVariables
{
   public:
    enum MeltingLineVariablesEnum
    {
        MELTING_LINE_NOT_SET = 0,
        MELTING_LINE_SIMON_TYPE,                ///< A simon-type curve is in use
        MELTING_LINE_POLYNOMIAL_IN_TR_TYPE,     ///< a polynomial in \f$ T/T_c \f$ is in use
        MELTING_LINE_POLYNOMIAL_IN_THETA_TYPE,  ///< a polynomial in \f$ \theta \f$ is in use
    };
    CoolPropDbl Tmin,  ///< Minimum temperature in K
      Tmax,            ///< Maximum temperature in K
      pmin,            ///< Minimum pressure in Pa
      pmax;            ///< Maximum pressure in Pa

    std::string BibTeX;                   ///< BibTeX key for the melting curve in use
    CoolPropDbl T_m;                      ///< Melting temperature at 1 atmosphere
    MeltingLinePiecewiseSimonData simon;  /// The data used for a Simon-style curve
    MeltingLinePiecewisePolynomialInTrData
      polynomial_in_Tr;  /// The data needed for a melting curve formed of segments that are polynomials in \f$ T/T_c \f$
    MeltingLinePiecewisePolynomialInThetaData
      polynomial_in_Theta;  /// The data needed for a melting curve formed of segments that are polynomials in \f$ \theta \f$
    int type;

    MeltingLineVariables() : Tmin(_HUGE), Tmax(_HUGE), pmin(_HUGE), pmax(_HUGE), T_m(_HUGE), type(MELTING_LINE_NOT_SET){};

    /**
     * \brief Evaluate the melting line
     * @param OF The output variable
     * @param GIVEN The given variable
     * @param value The value of the given variable
     */
    CoolPropDbl evaluate(int OF, int GIVEN, CoolPropDbl value);

    /// Evaluate the melting line to calculate the limits of the curve (Tmin/Tmax and pmin/pmax)
    void set_limits();

    /// Return true if the ancillary is enabled (type is not the default value of MELTING_LINE_NOT_SET)
    bool enabled() {
        return type != MELTING_LINE_NOT_SET;
    };
};

} /* namespace CoolProp */
#endif