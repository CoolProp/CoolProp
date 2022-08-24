#ifndef POLYMATH_H
#define POLYMATH_H

#include "CoolProp.h"
#include "CoolPropTools.h"
#include "Exceptions.h"

#include <vector>
#include <string>
#include "Solvers.h"
#include "MatrixMath.h"
#include "unsupported/Eigen/Polynomials"

namespace CoolProp {

// Just a forward declaration
class Poly2DResidual;
class Poly2DFracResidual;

/// The base class for all Polynomials
/** A clear and straight forward implementation of polynomial operations. Still
 *  very basic, but serves its purpose.
 */
class Polynomial2D
{

   public:
    /// Constructors
    Polynomial2D(){};

    /// Destructor.  No implementation
    virtual ~Polynomial2D(){};

   public:
    /// Convert the coefficient vector.
    /// @param coefficients vector containing the ordered coefficients
    Eigen::MatrixXd convertCoefficients(const std::vector<double>& coefficients) {
        return vec_to_eigen(coefficients);
    }
    /// Convert the coefficient matrix.
    /// @param coefficients matrix containing the ordered coefficients
    Eigen::MatrixXd convertCoefficients(const std::vector<std::vector<double>>& coefficients) {
        return vec_to_eigen(coefficients);
    }

    /// Basic checks for coefficient vectors.
    /** Starts with only the first coefficient dimension
     *  and checks the matrix size against the parameters rows and columns. */
    /// @param coefficients matrix containing the ordered coefficients
    /// @param rows unsigned integer value that represents the desired degree of the polynomial in the 1st dimension
    /// @param columns unsigned integer value that represents the desired degree of the polynomial in the 2nd dimension
    bool checkCoefficients(const Eigen::MatrixXd& coefficients, const unsigned int rows, const unsigned int columns);

   public:
    /// Integration functions
    /** Integrating coefficients for polynomials is done by dividing the
     *  original coefficients by (i+1) and elevating the order by 1
     *  through adding a zero as first coefficient.
     *  Some reslicing needs to be applied to integrate along the x-axis.
     *  In the brine/solution equations, reordering of the parameters
     *  avoids this expensive operation. However, it is included for the
     *  sake of completeness.
     */
    /// @param coefficients matrix containing the ordered coefficients
    /// @param axis unsigned integer value that represents the desired direction of integration
    /// @param times integer value that represents the desired order of integration
    Eigen::MatrixXd integrateCoeffs(const Eigen::MatrixXd& coefficients, const int& axis, const int& times);

    /// Derivative coefficients calculation
    /** Deriving coefficients for polynomials is done by multiplying the
     *  original coefficients with i and lowering the order by 1.
     */
    /// @param coefficients matrix containing the ordered coefficients
    /// @param axis unsigned integer value that represents the desired direction of derivation
    /// @param times integer value that represents the desired order of derivation
    Eigen::MatrixXd deriveCoeffs(const Eigen::MatrixXd& coefficients, const int& axis = -1, const int& times = 1);

   public:
    /// The core functions to evaluate the polynomial
    /** It is here we implement the different special
     *  functions that allow us to specify certain
     *  types of polynomials.
     *
     *  Try to avoid many calls to the derivative and integral functions.
     *  Both of them have to calculate the new coefficients internally,
     *  which slows things down. Instead, you should use the deriveCoeffs
     *  and integrateCoeffs functions and store the coefficient matrix
     *  you need for future calls to evaluate derivative and integral.
     */
    /// @param coefficients vector containing the ordered coefficients
    /// @param x_in double value that represents the current input in the 1st dimension
    double evaluate(const Eigen::MatrixXd& coefficients, const double& x_in);

    /// @param coefficients vector containing the ordered coefficients
    /// @param x_in double value that represents the current input in the 1st dimension
    /// @param y_in double value that represents the current input in the 2nd dimension
    double evaluate(const Eigen::MatrixXd& coefficients, const double& x_in, const double& y_in);

    /// @param coefficients vector containing the ordered coefficients
    /// @param x_in double value that represents the current input in the 1st dimension
    /// @param y_in double value that represents the current input in the 2nd dimension
    /// @param axis unsigned integer value that represents the axis to derive for (0=x, 1=y)
    double derivative(const Eigen::MatrixXd& coefficients, const double& x_in, const double& y_in, const int& axis);

    /// @param coefficients vector containing the ordered coefficients
    /// @param x_in double value that represents the current input in the 1st dimension
    /// @param y_in double value that represents the current input in the 2nd dimension
    /// @param axis unsigned integer value that represents the axis to integrate for (0=x, 1=y)
    double integral(const Eigen::MatrixXd& coefficients, const double& x_in, const double& y_in, const int& axis);

   protected:
    // TODO: Why doe these base definitions not work with derived classes?
    /// Uses the Brent solver to find the roots of p(x_in,y_in)-z_in
    /// @param res Poly2DResidual object to calculate residuals and derivatives
    /// @param min double value that represents the minimum value
    /// @param max double value that represents the maximum value
    double solve_limits(Poly2DResidual* res, const double& min, const double& max);

    // TODO: Why doe these base definitions not work with derived classes?
    /// Uses the Newton solver to find the roots of p(x_in,y_in)-z_in
    /// @param res Poly2DResidual object to calculate residuals and derivatives
    /// @param guess double value that represents the start value
    double solve_guess(Poly2DResidual* res, const double& guess);

   public:
    /// Returns a vector with ALL the real roots of p(x_in,y_in)-z_in
    /// @param coefficients vector containing the ordered coefficients
    /// @param in double value that represents the current input in x (1st dimension) or y (2nd dimension)
    /// @param z_in double value that represents the current output in the 3rd dimension
    /// @param axis unsigned integer value that represents the axis to solve for (0=x, 1=y)
    Eigen::VectorXd solve(const Eigen::MatrixXd& coefficients, const double& in, const double& z_in, const int& axis);

    /// Uses the Brent solver to find the roots of p(x_in,y_in)-z_in
    /// @param coefficients vector containing the ordered coefficients
    /// @param in double value that represents the current input in x (1st dimension) or y (2nd dimension)
    /// @param z_in double value that represents the current output in the 3rd dimension
    /// @param min double value that represents the minimum value
    /// @param max double value that represents the maximum value
    /// @param axis unsigned integer value that represents the axis to solve for (0=x, 1=y)
    double solve_limits(const Eigen::MatrixXd& coefficients, const double& in, const double& z_in, const double& min, const double& max,
                        const int& axis);

    /// Uses the Newton solver to find the roots of p(x_in,y_in)-z_in
    /// @param coefficients vector containing the ordered coefficients
    /// @param in double value that represents the current input in x (1st dimension) or y (2nd dimension)
    /// @param z_in double value that represents the current output in the 3rd dimension
    /// @param guess double value that represents the start value
    /// @param axis unsigned integer value that represents the axis to solve for (0=x, 1=y)
    double solve_guess(const Eigen::MatrixXd& coefficients, const double& in, const double& z_in, const double& guess, const int& axis);

   protected:
    /// Simple polynomial function generator. <- Deprecated due to poor performance, use Horner-scheme instead
    /** Base function to produce n-th order polynomials
     *  based on the length of the coefficient vector.
     *  Starts with only the first coefficient at x^0. */
    double simplePolynomial(const std::vector<double>& coefficients, double x);
    DEPRECATED(double simplePolynomial(const std::vector<std::vector<double>>& coefficients, double x, double y));
    /// Horner function generator implementations
    /** Represent polynomials according to Horner's scheme.
     *  This avoids unnecessary multiplication and thus
     *  speeds up calculation.
     *  Deprecated since we moved everything to the Eigen framework.
     */
    double baseHorner(const std::vector<double>& coefficients, double x);
    DEPRECATED(double baseHorner(const std::vector<std::vector<double>>& coefficients, double x, double y));

    bool do_debug(void) {
        return get_debug_level() >= 500;
    }
};

class Poly2DResidual : public FuncWrapper1DWithDeriv
{
   protected:
    enum dims
    {
        iX,
        iY
    };
    Eigen::MatrixXd coefficients;
    bool derIsSet;
    Eigen::MatrixXd coefficientsDer;
    int axis;
    /// the fixed input != targetDim
    double in;
    /// Object that evaluates the equation
    Polynomial2D poly;
    /// Current output value
    double z_in;

   protected:
    Poly2DResidual();

   public:
    /// Residual of a polynomial
    /// @param poly polynomial object used to evaluate the calls
    /// @param coefficients
    /// @param in double value that represents the current input in x (1st dimension) or y (2nd dimension)
    /// @param z_in double value that represents the current output in the 3rd dimension
    /// @param axis unsigned integer value that represents the axis to solve for (0=x, 1=y)
    Poly2DResidual(Polynomial2D& poly, const Eigen::MatrixXd& coefficients, const double& in, const double& z_in, const int& axis);
    virtual ~Poly2DResidual(){};

    double call(double target);
    double deriv(double target);
};

/// A class for polynomials starting at an arbitrary degree.
/** It is implemented for the incompressibles and is a little messy, but seems to
 *  work fine for now. Besides handling arbitrary starting exponents for the
 *  polynomials, it can also calculate polynomials with a base value. This means
 *  that the independent variable no longer is x, but (x-x_base). For fitted
 *  functions, we often see such a design to enhance the fit quality/stability.
 */
class Polynomial2DFrac : public Polynomial2D
{

   public:
    /// Constructors
    Polynomial2DFrac(){};

    /// Destructor.  No implementation
    virtual ~Polynomial2DFrac(){};

   public:
    //    /// Integration functions
    //    /** Integrating coefficients for polynomials is done by dividing the
    //     *  original coefficients by (i+1) and elevating the order by 1
    //     *  through adding a zero as first coefficient.
    //     *  Some reslicing needs to be applied to integrate along the x-axis.
    //     *  In the brine/solution equations, reordering of the parameters
    //     *  avoids this expensive operation. However, it is included for the
    //     *  sake of completeness.
    //     */
    //    /// @param coefficients matrix containing the ordered coefficients
    //    /// @param axis unsigned integer value that represents the desired direction of integration
    //    /// @param times integer value that represents the desired order of integration
    //    /// @param firstExponent integer value that represents the first exponent of the polynomial in axis direction
    //    Eigen::MatrixXd integrateCoeffs(const Eigen::MatrixXd &coefficients, const int &axis, const int &times, const int &firstExponent);
    //
    /// Derivative coefficients calculation
    /** Deriving coefficients for polynomials is done by multiplying the
     *  original coefficients with i and lowering the order by 1.
     *
     *  Remember that the first exponent might need to be adjusted after derivation.
     *  It has to be lowered by times:
     *  derCoeffs = deriveCoeffs(coefficients, axis, times, firstExponent);
     *  firstExponent -= times;
     */
    /// @param coefficients matrix containing the ordered coefficients
    /// @param axis unsigned integer value that represents the desired direction of derivation
    /// @param times integer value that represents the desired order of derivation
    /// @param firstExponent integer value that represents the lowest exponent of the polynomial in axis direction
    Eigen::MatrixXd deriveCoeffs(const Eigen::MatrixXd& coefficients, const int& axis, const int& times, const int& firstExponent);

   public:
    /// The core functions to evaluate the polynomial
    /** It is here we implement the different special
     *  functions that allow us to specify certain
     *  types of polynomials.
     *
     *  Try to avoid many calls to the derivative and integral functions.
     *  Both of them have to calculate the new coefficients internally,
     *  which slows things down. Instead, you should use the deriveCoeffs
     *  and integrateCoeffs functions and store the coefficient matrix
     *  you need for future calls to evaluate derivative and integral.
     */
    /// @param coefficients vector containing the ordered coefficients
    /// @param x_in double value that represents the current input in the 1st dimension
    /// @param firstExponent integer value that represents the lowest exponent of the polynomial
    /// @param x_base double value that represents the base value for a centered fit in the 1st dimension
    double evaluate(const Eigen::MatrixXd& coefficients, const double& x_in, const int& firstExponent = 0, const double& x_base = 0.0);

    /// @param coefficients matrix containing the ordered coefficients
    /// @param x_in double value that represents the current input in the 1st dimension
    /// @param y_in double value that represents the current input in the 2nd dimension
    /// @param x_exp integer value that represents the lowest exponent of the polynomial in the 1st dimension
    /// @param y_exp integer value that represents the lowest exponent of the polynomial in the 2nd dimension
    /// @param x_base double value that represents the base value for a centered fit in the 1st dimension
    /// @param y_base double value that represents the base value for a centered fit in the 2nd dimension
    double evaluate(const Eigen::MatrixXd& coefficients, const double& x_in, const double& y_in, const int& x_exp, const int& y_exp,
                    const double& x_base = 0.0, const double& y_base = 0.0);

    /// @param coefficients vector containing the ordered coefficients
    /// @param x_in double value that represents the current input in the 1st dimension
    /// @param y_in double value that represents the current input in the 2nd dimension
    /// @param axis unsigned integer value that represents the axis to derive for (0=x, 1=y)
    /// @param x_exp integer value that represents the lowest exponent of the polynomial in the 1st dimension
    /// @param y_exp integer value that represents the lowest exponent of the polynomial in the 2nd dimension
    /// @param x_base double value that represents the base value for a centred fit in the 1st dimension
    /// @param y_base double value that represents the base value for a centred fit in the 2nd dimension
    double derivative(const Eigen::MatrixXd& coefficients, const double& x_in, const double& y_in, const int& axis, const int& x_exp,
                      const int& y_exp, const double& x_base = 0.0, const double& y_base = 0.0);

    /// @param coefficients vector containing the ordered coefficients
    /// @param x_in double value that represents the current input in the 1st dimension
    /// @param y_in double value that represents the current input in the 2nd dimension
    /// @param axis unsigned integer value that represents the axis to integrate for (0=x, 1=y)
    /// @param x_exp integer value that represents the lowest exponent of the polynomial in the 1st dimension
    /// @param y_exp integer value that represents the lowest exponent of the polynomial in the 2nd dimension
    /// @param x_base double value that represents the base value for a centred fit in the 1st dimension
    /// @param y_base double value that represents the base value for a centred fit in the 2nd dimension
    /// @param ax_val double value that represents the base value for the definite integral on the chosen axis.
    double integral(const Eigen::MatrixXd& coefficients, const double& x_in, const double& y_in, const int& axis, const int& x_exp, const int& y_exp,
                    const double& x_base = 0.0, const double& y_base = 0.0, const double& ax_val = 0.0);

   public:
    /// Returns a vector with ALL the real roots of p(x_in,y_in)-z_in
    /// @param coefficients vector containing the ordered coefficients
    /// @param in double value that represents the current input in x (1st dimension) or y (2nd dimension)
    /// @param z_in double value that represents the current output in the 3rd dimension
    /// @param axis unsigned integer value that represents the axis to solve for (0=x, 1=y)
    /// @param x_exp integer value that represents the lowest exponent of the polynomial in the 1st dimension
    /// @param y_exp integer value that represents the lowest exponent of the polynomial in the 2nd dimension
    /// @param x_base double value that represents the base value for a centred fit in the 1st dimension
    /// @param y_base double value that represents the base value for a centred fit in the 2nd dimension
    Eigen::VectorXd solve(const Eigen::MatrixXd& coefficients, const double& in, const double& z_in, const int& axis, const int& x_exp,
                          const int& y_exp, const double& x_base = 0.0, const double& y_base = 0.0);

    /// Uses the Brent solver to find the roots of p(x_in,y_in)-z_in
    /// @param coefficients vector containing the ordered coefficients
    /// @param in double value that represents the current input in x (1st dimension) or y (2nd dimension)
    /// @param z_in double value that represents the current output in the 3rd dimension
    /// @param min double value that represents the minimum value
    /// @param max double value that represents the maximum value
    /// @param axis unsigned integer value that represents the axis to solve for (0=x, 1=y)
    /// @param x_exp integer value that represents the lowest exponent of the polynomial in the 1st dimension
    /// @param y_exp integer value that represents the lowest exponent of the polynomial in the 2nd dimension
    /// @param x_base double value that represents the base value for a centred fit in the 1st dimension
    /// @param y_base double value that represents the base value for a centred fit in the 2nd dimension
    double solve_limits(const Eigen::MatrixXd& coefficients, const double& in, const double& z_in, const double& min, const double& max,
                        const int& axis, const int& x_exp, const int& y_exp, const double& x_base = 0.0, const double& y_base = 0.0);

    /// Uses the Newton solver to find the roots of p(x_in,y_in)-z_in
    /// @param coefficients vector containing the ordered coefficients
    /// @param in double value that represents the current input in x (1st dimension) or y (2nd dimension)
    /// @param z_in double value that represents the current output in the 3rd dimension
    /// @param guess double value that represents the start value
    /// @param axis unsigned integer value that represents the axis to solve for (0=x, 1=y)
    /// @param x_exp integer value that represents the lowest exponent of the polynomial in the 1st dimension
    /// @param y_exp integer value that represents the lowest exponent of the polynomial in the 2nd dimension
    /// @param x_base double value that represents the base value for a centred fit in the 1st dimension
    /// @param y_base double value that represents the base value for a centred fit in the 2nd dimension
    double solve_guess(const Eigen::MatrixXd& coefficients, const double& in, const double& z_in, const double& guess, const int& axis,
                       const int& x_exp, const int& y_exp, const double& x_base = 0.0, const double& y_base = 0.0);

    /// Uses the Brent solver to find the roots of Int(p(x_in,y_in))-z_in
    /// @param coefficients vector containing the ordered coefficients
    /// @param in double value that represents the current input in x (1st dimension) or y (2nd dimension)
    /// @param z_in double value that represents the current output in the 3rd dimension
    /// @param min double value that represents the minimum value
    /// @param max double value that represents the maximum value
    /// @param axis unsigned integer value that represents the axis to solve for (0=x, 1=y)
    /// @param x_exp integer value that represents the lowest exponent of the polynomial in the 1st dimension
    /// @param y_exp integer value that represents the lowest exponent of the polynomial in the 2nd dimension
    /// @param x_base double value that represents the base value for a centred fit in the 1st dimension
    /// @param y_base double value that represents the base value for a centred fit in the 2nd dimension
    /// @param int_axis axis for the integration (0=x, 1=y)
    double solve_limitsInt(const Eigen::MatrixXd& coefficients, const double& in, const double& z_in, const double& min, const double& max,
                           const int& axis, const int& x_exp, const int& y_exp, const double& x_base = 0.0, const double& y_base = 0.0,
                           const int& int_axis = 0);

    /// Uses the Newton solver to find the roots of Int(p(x_in,y_in))-z_in
    /// @param coefficients vector containing the ordered coefficients
    /// @param in double value that represents the current input in x (1st dimension) or y (2nd dimension)
    /// @param z_in double value that represents the current output in the 3rd dimension
    /// @param guess double value that represents the start value
    /// @param axis unsigned integer value that represents the axis to solve for (0=x, 1=y)
    /// @param x_exp integer value that represents the lowest exponent of the polynomial in the 1st dimension
    /// @param y_exp integer value that represents the lowest exponent of the polynomial in the 2nd dimension
    /// @param x_base double value that represents the base value for a centred fit in the 1st dimension
    /// @param y_base double value that represents the base value for a centred fit in the 2nd dimension
    /// @param int_axis axis for the integration (0=x, 1=y)
    double solve_guessInt(const Eigen::MatrixXd& coefficients, const double& in, const double& z_in, const double& guess, const int& axis,
                          const int& x_exp, const int& y_exp, const double& x_base = 0.0, const double& y_base = 0.0, const int& int_axis = 0);

   protected:
    /// @param nValue integer value that represents the order of the factorial
    double factorial(const int& nValue);

    /// @param nValue integer value that represents the upper part of the factorial
    /// @param nValue2 integer value that represents the lower part of the factorial
    double binom(const int& nValue, const int& nValue2);

    ///Helper function to calculate the D vector:
    /// @param m integer value that represents order
    /// @param x_in double value that represents the current input
    /// @param x_base double value that represents the basis for the fit
    Eigen::MatrixXd fracIntCentralDvector(const int& m, const double& x_in, const double& x_base);

    ///Indefinite integral of a centred polynomial divided by its independent variable
    /// @param coefficients vector containing the ordered coefficients
    /// @param x_in double value that represents the current input
    /// @param x_base double value that represents the basis for the fit
    double fracIntCentral(const Eigen::MatrixXd& coefficients, const double& x_in, const double& x_base);
};

class Poly2DFracResidual : public Poly2DResidual
{
   protected:
    int x_exp, y_exp;
    double x_base, y_base;
    /// Object that evaluates the equation
    Polynomial2DFrac poly;

   protected:
    Poly2DFracResidual();

   public:
    /// Residual of a polynomial divided by the independent variable
    /// @param poly polynomial object used to evaluate the calls
    /// @param coefficients
    /// @param in double value that represents the current input in x (1st dimension) or y (2nd dimension)
    /// @param z_in double value that represents the current output in the 3rd dimension
    /// @param axis unsigned integer value that represents the axis to solve for (0=x, 1=y)
    /// @param x_exp first exponent in x-direction
    /// @param y_exp first exponent in y-direction
    /// @param x_base base value for x (x = x_in - x_base)
    /// @param y_base base value for y (y = y_in - y_base)
    Poly2DFracResidual(Polynomial2DFrac& poly, const Eigen::MatrixXd& coefficients, const double& in, const double& z_in, const int& axis,
                       const int& x_exp, const int& y_exp, const double& x_base, const double& y_base);
    virtual ~Poly2DFracResidual(){};
    double call(double target);
    double deriv(double target);
};

class Poly2DFracIntResidual : public Poly2DFracResidual
{

   protected:
    int int_axis;
    Poly2DFracIntResidual();

   public:
    /// Residual of an integrated polynomial divided by the independent variable
    /// @param poly polynomial object used to evaluate the calls
    /// @param coefficients vector of coefficients
    /// @param in double value that represents the current input in x (1st dimension) or y (2nd dimension)
    /// @param z_in double value that represents the current output in the 3rd dimension
    /// @param axis unsigned integer value that represents the axis to solve for (0=x, 1=y)
    /// @param x_exp first exponent in x-direction
    /// @param y_exp first exponent in y-direction
    /// @param x_base base value for x (x = x_in - x_base)
    /// @param y_base base value for y (y = y_in - y_base)
    /// @param int_axis axis for the integration (0=x, 1=y)
    Poly2DFracIntResidual(Polynomial2DFrac& poly, const Eigen::MatrixXd& coefficients, const double& in, const double& z_in, const int& axis,
                          const int& x_exp, const int& y_exp, const double& x_base, const double& y_base, const int& int_axis);
    virtual ~Poly2DFracIntResidual(){};
    double call(double target);
    double deriv(double target);
};

//
//
//
//
//
//
//
//
//
///// The base class for Polynomials
//class BasePolynomial{
//
//public:
//    // Constructor
//    BasePolynomial();
//    // Destructor.  No implementation
//    virtual ~BasePolynomial(){};
//
//public:
//    /// Basic checks for coefficient vectors.
//    /** Starts with only the first coefficient dimension
//     *  and checks the vector length against parameter n. */
//    bool checkCoefficients(const Eigen::VectorXd &coefficients, const unsigned int n);
//    bool checkCoefficients(const Eigen::MatrixXd &coefficients, const unsigned int rows, const unsigned int columns);
//    bool checkCoefficients(const std::vector<double> &coefficients, const unsigned int n);
//    bool checkCoefficients(const std::vector< std::vector<double> > &coefficients, const unsigned int rows, const unsigned int columns);
//
//    /** Integrating coefficients for polynomials is done by dividing the
//     *  original coefficients by (i+1) and elevating the order by 1
//     *  through adding a zero as first coefficient.
//     *  Some reslicing needs to be applied to integrate along the x-axis.
//     *  In the brine/solution equations, reordering of the parameters
//     *  avoids this expensive operation. However, it is included for the
//     *  sake of completeness.
//     */
//    std::vector<double> integrateCoeffs(const std::vector<double> &coefficients);
//    std::vector< std::vector<double> > integrateCoeffs(const std::vector< std::vector<double> > &coefficients, bool axis);
//
//    /** Deriving coefficients for polynomials is done by multiplying the
//     *  original coefficients with i and lowering the order by 1.
//     *
//     *  It is not really deprecated, but untested and therefore a warning
//     *  is issued. Please check this method before you use it.
//     */
//    std::vector<double> deriveCoeffs(const std::vector<double> &coefficients);
//    std::vector< std::vector<double> > deriveCoeffs(const std::vector< std::vector<double> > &coefficients, unsigned int axis);
//
//private:
//    /** The core of the polynomial wrappers are the different
//     *  implementations that follow below. In case there are
//     *  new calculation schemes available, please do not delete
//     *  the implementations, but mark them as deprecated.
//     *  The old functions are good for debugging since the
//     *  structure is easier to read than the backward Horner-scheme
//     *  or the recursive Horner-scheme.
//     */
//
//    /// Simple polynomial function generator. <- Deprecated due to poor performance, use Horner-scheme instead
//    /** Base function to produce n-th order polynomials
//     *  based on the length of the coefficient vector.
//     *  Starts with only the first coefficient at x^0. */
//    DEPRECATED(double simplePolynomial(const std::vector<double> &coefficients, double x));
//    DEPRECATED(double simplePolynomial(const std::vector<std::vector<double> > &coefficients, double x, double y));
//
//    /// Simple integrated polynomial function generator.
//    /** Base function to produce integrals of n-th order polynomials based on
//     *  the length of the coefficient vector.
//     *  Starts with only the first coefficient at x^0 */
//    ///Indefinite integral in x-direction
//    double simplePolynomialInt(const std::vector<double> &coefficients, double x);
//    ///Indefinite integral in y-direction only
//    double simplePolynomialInt(const std::vector<std::vector<double> > &coefficients, double x, double y);
//
//    /// Simple integrated polynomial function generator divided by independent variable.
//    /** Base function to produce integrals of n-th order
//     *  polynomials based on the length of the coefficient
//     *  vector. Starts with only the first coefficient at x^0 */
//    ///Indefinite integral of a polynomial divided by its independent variable
//    double simpleFracInt(const std::vector<double> &coefficients, double x);
//    ///Indefinite integral of a polynomial divided by its 2nd independent variable
//    double simpleFracInt(const std::vector<std::vector<double> > &coefficients, double x, double y);
//
//    /** Simple integrated centred(!) polynomial function generator divided by independent variable.
//     *  We need to rewrite some of the functions in order to
//     *  use central fit. Having a central temperature xbase
//     *  allows for a better fit, but requires a different
//     *  formulation of the fracInt function group. Other
//     *  functions are not affected.
//     *  Starts with only the first coefficient at x^0 */
//    ///Helper function to calculate the D vector:
//    double factorial(double nValue);
//    double binom(double nValue, double nValue2);
//    std::vector<double> fracIntCentralDvector(int m, double x, double xbase);
//    ///Indefinite integral of a centred polynomial divided by its independent variable
//    double fracIntCentral(const std::vector<double> &coefficients, double x, double xbase);
//
//    /// Horner function generator implementations
//    /** Represent polynomials according to Horner's scheme.
//     *  This avoids unnecessary multiplication and thus
//     *  speeds up calculation.
//     */
//    double baseHorner(const std::vector<double> &coefficients, double x);
//    double baseHorner(const std::vector< std::vector<double> > &coefficients, double x, double y);
//    ///Indefinite integral in x-direction
//    double baseHornerInt(const std::vector<double> &coefficients, double x);
//    ///Indefinite integral in y-direction only
//    double baseHornerInt(const std::vector<std::vector<double> > &coefficients, double x, double y);
//    ///Indefinite integral of a polynomial divided by its independent variable
//    double baseHornerFracInt(const std::vector<double> &coefficients, double x);
//    ///Indefinite integral of a polynomial divided by its 2nd independent variable
//    double baseHornerFracInt(const std::vector<std::vector<double> > &coefficients, double x, double y);
//
//    /** Alternatives
//     *  Simple functions that heavily rely on other parts of this file.
//     *  We still need to check which combinations yield the best
//     *  performance.
//     */
//    ///Derivative in x-direction
//    double deriveIn2Steps(const std::vector<double> &coefficients, double x); // TODO: Check results!
//    ///Derivative in terms of x(axis=true) or y(axis=false).
//    double deriveIn2Steps(const std::vector< std::vector<double> > &coefficients, double x, double y, bool axis); // TODO: Check results!
//    ///Indefinite integral in x-direction
//    double integrateIn2Steps(const std::vector<double> &coefficients, double x);
//    ///Indefinite integral in terms of x(axis=true) or y(axis=false).
//    double integrateIn2Steps(const std::vector< std::vector<double> > &coefficients, double x, double y, bool axis);
//    ///Indefinite integral in x-direction of a polynomial divided by its independent variable
//    double fracIntIn2Steps(const std::vector<double> &coefficients, double x);
//    ///Indefinite integral in y-direction of a polynomial divided by its 2nd independent variable
//    double fracIntIn2Steps(const std::vector<std::vector<double> > &coefficients, double x, double y);
//    ///Indefinite integral of a centred polynomial divided by its 2nd independent variable
//    double fracIntCentral2Steps(const std::vector<std::vector<double> > &coefficients, double x, double y, double ybase);
//
//public:
//    /** Here we define the functions that should be used by the
//     *  respective implementations. Please do no use any other
//     *  method since this would break the purpose of this interface.
//     *  Note that the functions below are supposed to be aliases
//     *  to implementations declared elsewhere in this file.
//     */
//
//    /** Everything related to the normal polynomials goes in this
//     *  section, holds all the functions for evaluating polynomials.
//     */
//    /// Evaluates a one-dimensional polynomial for the given coefficients
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param x double value that represents the current input
//    virtual inline double polyval(const std::vector<double> &coefficients, double x){
//        return baseHorner(coefficients,x);
//    }
//
//    /// Evaluates a two-dimensional polynomial for the given coefficients
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param x double value that represents the current input in the 1st dimension
//    /// @param y double value that represents the current input in the 2nd dimension
//    virtual inline double polyval(const std::vector< std::vector<double> > &coefficients, double x, double y){
//        return baseHorner(coefficients,x,y);
//    }
//
//
//    /** Everything related to the integrated polynomials goes in this
//     *  section, holds all the functions for evaluating polynomials.
//     */
//    /// Evaluates the indefinite integral of a one-dimensional polynomial
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param x double value that represents the current input
//    virtual inline double polyint(const std::vector<double> &coefficients, double x){
//        return baseHornerInt(coefficients,x);
//    }
//
//    /// Evaluates the indefinite integral of a two-dimensional polynomial along the 2nd axis (y)
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param x double value that represents the current input in the 1st dimension
//    /// @param y double value that represents the current input in the 2nd dimension
//    virtual inline double polyint(const std::vector< std::vector<double> > &coefficients, double x, double y){
//        return baseHornerInt(coefficients,x,y);
//    }
//
//
//    /** Everything related to the derived polynomials goes in this
//     *  section, holds all the functions for evaluating polynomials.
//     */
//    /// Evaluates the derivative of a one-dimensional polynomial
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param x double value that represents the current input
//    virtual inline double polyder(const std::vector<double> &coefficients, double x){
//        return deriveIn2Steps(coefficients,x);
//    }
//
//    /// Evaluates the derivative of a two-dimensional polynomial along the 2nd axis (y)
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param x double value that represents the current input in the 1st dimension
//    /// @param y double value that represents the current input in the 2nd dimension
//    virtual inline double polyder(const std::vector< std::vector<double> > &coefficients, double x, double y){
//        return deriveIn2Steps(coefficients,x,y,false);
//    }
//
//
//    /** Everything related to the polynomials divided by one variable goes in this
//     *  section, holds all the functions for evaluating polynomials.
//     */
//    /// Evaluates the indefinite integral of a one-dimensional polynomial divided by its independent variable
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param x double value that represents the current position
//    virtual inline double polyfracval(const std::vector<double> &coefficients, double x){
//        return baseHorner(coefficients,x)/x;
//    }
//
//    /// Evaluates the indefinite integral of a two-dimensional polynomial divided by its 2nd independent variable
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param x double value that represents the current input in the 1st dimension
//    /// @param y double value that represents the current input in the 2nd dimension
//    virtual inline double polyfracval(const std::vector< std::vector<double> > &coefficients, double x, double y){
//        return baseHorner(coefficients,x,y)/y;
//    }
//
//
//    /** Everything related to the integrated polynomials divided by one variable goes in this
//     *  section, holds all the functions for solving polynomials.
//     */
//    /// Evaluates the indefinite integral of a one-dimensional polynomial divided by its independent variable
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param x double value that represents the current position
//    virtual inline double polyfracint(const std::vector<double> &coefficients, double x){
//        return baseHornerFracInt(coefficients,x);
//    }
//
//    /// Evaluates the indefinite integral of a two-dimensional polynomial divided by its 2nd independent variable
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param x double value that represents the current input in the 1st dimension
//    /// @param y double value that represents the current input in the 2nd dimension
//    virtual inline double polyfracint(const std::vector< std::vector<double> > &coefficients, double x, double y){
//        return baseHornerFracInt(coefficients,x,y);
//    }
//
//    /// Evaluates the indefinite integral of a centred one-dimensional polynomial divided by its independent variable
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param x double value that represents the current position
//    /// @param xbase central temperature for fitted function
//    virtual inline double polyfracintcentral(const std::vector<double> &coefficients, double x, double xbase){
//        return fracIntCentral(coefficients,x,xbase);
//    }
//
//    /// Evaluates the indefinite integral of a centred two-dimensional polynomial divided by its 2nd independent variable
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param x double value that represents the current input in the 1st dimension
//    /// @param y double value that represents the current input in the 2nd dimension
//    /// @param ybase central temperature for fitted function
//    virtual inline double polyfracintcentral(const std::vector< std::vector<double> > &coefficients, double x, double y, double ybase){
//        return fracIntCentral2Steps(coefficients,x,y,ybase);
//    }
//
//
//    /** Everything related to the derived polynomials divided by one variable goes in this
//     *  section, holds all the functions for solving polynomials.
//     */
//    /// Evaluates the derivative of a one-dimensional polynomial divided by its independent variable
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param x double value that represents the current position
//    virtual inline double polyfracder(const std::vector<double> &coefficients, double x){
//        throw CoolProp::NotImplementedError("Derivatives of polynomials divided by their independent variable have not been implemented."); // TODO: Implement polyfracder1D
//    }
//
//    /// Evaluates the derivative of a two-dimensional polynomial divided by its 2nd independent variable
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param x double value that represents the current input in the 1st dimension
//    /// @param y double value that represents the current input in the 2nd dimension
//    virtual inline double polyfracder(const std::vector< std::vector<double> > &coefficients, double x, double y){
//        throw CoolProp::NotImplementedError("Derivatives of polynomials divided by their independent variable have not been implemented."); // TODO: Implement polyfracder2D
//    }
//
//    /// Evaluates the derivative of a centred one-dimensional polynomial divided by its independent variable
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param x double value that represents the current position
//    /// @param xbase central temperature for fitted function
//    virtual inline double polyfracdercentral(const std::vector<double> &coefficients, double x, double xbase){
//        throw CoolProp::NotImplementedError("Derivatives of polynomials divided by their independent variable have not been implemented."); // TODO: Implement polyfracdercentral1D
//    }
//
//    /// Evaluates the derivative of a centred two-dimensional polynomial divided by its 2nd independent variable
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param x double value that represents the current input in the 1st dimension
//    /// @param y double value that represents the current input in the 2nd dimension
//    /// @param ybase central temperature for fitted function
//    virtual inline double polyfracdercentral(const std::vector< std::vector<double> > &coefficients, double x, double y, double ybase){
//        throw CoolProp::NotImplementedError("Derivatives of polynomials divided by their independent variable have not been implemented."); // TODO: Implement polyfracdercentral2D
//    }
//};
//
//
//
//
///** Implements the function wrapper interface and can be
// *  used by the solvers.
// *  TODO: Make multidimensional
// */
//class PolyResidual : public FuncWrapper1D {
//protected:
//    enum dims {i1D, i2D};
//    /// Object that evaluates the equation
//    BasePolynomial poly;
//    /// Current output value
//    double output, firstDim;
//    int dim;
//    std::vector< std::vector<double> > coefficients;
//private:
//    PolyResidual();
//public:
//    PolyResidual(const std::vector<double> &coefficients, double y);
//    PolyResidual(const std::vector< std::vector<double> > &coefficients, double x, double z);
//    virtual ~PolyResidual(){};
//    bool is2D(){return (this->dim==i2D);};
//    virtual double call(double x);
//    virtual double deriv(double x);
//};
//class PolyIntResidual : public PolyResidual {
//public:
//    PolyIntResidual(const std::vector<double> &coefficients, double y):PolyResidual(coefficients, y){};
//    PolyIntResidual(const std::vector< std::vector<double> > &coefficients, double x, double z):PolyResidual(coefficients, x, z){};
//    virtual double call(double x);
//    virtual double deriv(double x);
//};
//class PolyFracIntResidual : public PolyResidual {
//public:
//    PolyFracIntResidual(const std::vector<double> &coefficients, double y):PolyResidual(coefficients, y){};
//    PolyFracIntResidual(const std::vector< std::vector<double> > &coefficients, double x, double z):PolyResidual(coefficients, x, z){};
//    virtual double call(double x);
//    virtual double deriv(double x);
//};
//class PolyFracIntCentralResidual : public PolyResidual {
//protected:
//    double baseVal;
//public:
//    PolyFracIntCentralResidual(const std::vector<double> &coefficients, double y, double xBase):PolyResidual(coefficients, y){this->baseVal = xBase;};
//    PolyFracIntCentralResidual(const std::vector< std::vector<double> > &coefficients, double x, double z, double yBase): PolyResidual(coefficients, x, z){this->baseVal = yBase;};
//    virtual double call(double x);
//    virtual double deriv(double x);
//};
//class PolyDerResidual : public PolyResidual {
//public:
//    PolyDerResidual(const std::vector<double> &coefficients, double y):PolyResidual(coefficients, y){};
//    PolyDerResidual(const std::vector< std::vector<double> > &coefficients, double x, double z):PolyResidual(coefficients, x, z){};
//    virtual double call(double x);
//    virtual double deriv(double x);
//};
//
//
//
//
///** Implements the same public functions as the
// *  but solves the polynomial for the given value
// *  instead of evaluating it.
// *  TODO: This class does not check for bijective
// *        polynomials and is therefore a little
// *        fragile.
// */
//class PolynomialSolver : public BasePolynomial{
//private:
//    enum solvers {iNewton, iBrent};
//    int uses;
//    double guess, min, max;
//    double macheps, tol;
//    int maxiter;
//
//public:
//    // Constructor
//    PolynomialSolver();
//    // Destructor.  No implementation
//    virtual ~PolynomialSolver(){};
//
//public:
//    /** Here we redefine the functions that solve the polynomials.
//     *  These implementations all use the base class to evaluate
//     *  the polynomial during the solution process.
//     */
//
//    /** Everything related to the normal polynomials goes in this
//     *  section, holds all the functions for solving polynomials.
//     */
//    /// Solves a one-dimensional polynomial for the given coefficients
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param y double value that represents the current input
//    virtual double polyval(const std::vector<double> &coefficients, double y);
//
//    /// Solves a two-dimensional polynomial for the given coefficients
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param x double value that represents the current input in the 1st dimension
//    /// @param z double value that represents the current output
//    virtual double polyval(const std::vector< std::vector<double> > &coefficients, double x, double z);
//
//
//    /** Everything related to the integrated polynomials goes in this
//     *  section, holds all the functions for solving polynomials.
//     */
//    /// Solves the indefinite integral of a one-dimensional polynomial
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param y double value that represents the current output
//    virtual double polyint(const std::vector<double> &coefficients, double y);
//
//    /// Solves the indefinite integral of a two-dimensional polynomial along the 2nd axis (y)
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param x double value that represents the current input in the 1st dimension
//    /// @param z double value that represents the current output
//    virtual double polyint(const std::vector< std::vector<double> > &coefficients, double x, double z);
//
//
//    /** Everything related to the derived polynomials goes in this
//     *  section, holds all the functions for solving polynomials.
//     */
//    /// Solves the derivative of a one-dimensional polynomial
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param y double value that represents the current output
//    virtual double polyder(const std::vector<double> &coefficients, double y);
//
//    /// Solves the derivative of a two-dimensional polynomial along the 2nd axis (y)
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param x double value that represents the current input in the 1st dimension
//    /// @param z double value that represents the current output
//    virtual double polyder(const std::vector< std::vector<double> > &coefficients, double x, double z);
//
//
//    /** Everything related to the polynomials divided by one variable goes in this
//     *  section, holds all the functions for solving polynomials.
//     */
//    /// Solves the indefinite integral of a one-dimensional polynomial divided by its independent variable
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param y double value that represents the current output
//    virtual double polyfracval(const std::vector<double> &coefficients, double y);
//
//    /// Solves the indefinite integral of a two-dimensional polynomial divided by its 2nd independent variable
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param x double value that represents the current input in the 1st dimension
//    /// @param z double value that represents the current output
//    virtual double polyfracval(const std::vector< std::vector<double> > &coefficients, double x, double z);
//
//
//    /** Everything related to the integrated polynomials divided by one variable goes in this
//     *  section, holds all the functions for solving polynomials.
//     */
//    /// Solves the indefinite integral of a one-dimensional polynomial divided by its independent variable
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param y double value that represents the current output
//    virtual double polyfracint(const std::vector<double> &coefficients, double y);
//
//    /// Solves the indefinite integral of a two-dimensional polynomial divided by its 2nd independent variable
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param x double value that represents the current input in the 1st dimension
//    /// @param z double value that represents the current output
//    virtual double polyfracint(const std::vector< std::vector<double> > &coefficients, double x, double z);
//
//    /// Solves the indefinite integral of a centred one-dimensional polynomial divided by its independent variable
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param y double value that represents the current output
//    /// @param xbase central x-value for fitted function
//    virtual double polyfracintcentral(const std::vector<double> &coefficients, double y, double xbase);
//
//    /// Solves the indefinite integral of a centred two-dimensional polynomial divided by its 2nd independent variable
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param x double value that represents the current input in the 1st dimension
//    /// @param z double value that represents the current output
//    /// @param ybase central y-value for fitted function
//    virtual double polyfracintcentral(const std::vector< std::vector<double> > &coefficients, double x, double z, double ybase);
//
//
//    /** Everything related to the derived polynomials divided by one variable goes in this
//     *  section, holds all the functions for solving polynomials.
//     */
//    /// Solves the derivative of a one-dimensional polynomial divided by its independent variable
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param y double value that represents the current output
//    virtual double polyfracder(const std::vector<double> &coefficients, double y);
//
//    /// Solves the derivative of a two-dimensional polynomial divided by its 2nd independent variable
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param x double value that represents the current input in the 1st dimension
//    /// @param z double value that represents the current output
//    virtual double polyfracder(const std::vector< std::vector<double> > &coefficients, double x, double z);
//
//    /// Solves the derivative of a centred one-dimensional polynomial divided by its independent variable
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param y double value that represents the current output
//    /// @param xbase central x-value for fitted function
//    virtual double polyfracdercentral(const std::vector<double> &coefficients, double y, double xbase);
//
//    /// Solves the derivative of a centred two-dimensional polynomial divided by its 2nd independent variable
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param x double value that represents the current input in the 1st dimension
//    /// @param z double value that represents the current output
//    /// @param ybase central y-value for fitted function
//    virtual double polyfracdercentral(const std::vector< std::vector<double> > &coefficients, double x, double z, double ybase);
//
//
//    /** Set the solvers and updates either the guess values or the
//     *  boundaries for the variable to solve for.
//     */
//    /// Sets the guess value for the Newton solver and enables it.
//    /// @param guess double value that represents the guess value
//    virtual void setGuess(double guess);
//    /// Sets the limits for the Brent solver and enables it.
//    /// @param min double value that represents the lower boundary
//    /// @param max double value that represents the upper boundary
//    virtual void setLimits(double min, double max);
//    /// Solves the equations based on previously defined parameters.
//    /// @param min double value that represents the lower boundary
//    /// @param max double value that represents the upper boundary
//    virtual double solve(PolyResidual &res);
//};
//
//
///// The base class for exponential functions
//class BaseExponential{
//
//protected:
//    BasePolynomial poly;
//    bool POLYMATH_DEBUG;
//
//public:
//    BaseExponential();
//    virtual ~BaseExponential(){};
//
//public:
//    /// Evaluates an exponential function for the given coefficients
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param x double value that represents the current input
//    /// @param n int value that determines the kind of exponential function
//    double expval(const std::vector<double> &coefficients, double x, int n);
//
//    /// Evaluates an exponential function for the given coefficients
//    /// @param coefficients vector containing the ordered coefficients
//    /// @param x double value that represents the current input in the 1st dimension
//    /// @param y double value that represents the current input in the 2nd dimension
//    /// @param n int value that determines the kind of exponential function
//    double expval(const std::vector< std::vector<double> > &coefficients, double x, double y, int n);
//};

}; /* namespace CoolProp */
#endif
