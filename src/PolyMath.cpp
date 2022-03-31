
#include "PolyMath.h"

#include "Exceptions.h"
#include "MatrixMath.h"

#include <vector>
#include <string>
//#include <sstream>
//#include <numeric>
#include <math.h>
#include <iostream>

#include "Solvers.h"

#include "unsupported/Eigen/Polynomials"

namespace CoolProp {

/// Basic checks for coefficient vectors.
/** Starts with only the first coefficient dimension
 *  and checks the matrix size against the parameters rows and columns.
 */
/// @param coefficients matrix containing the ordered coefficients
/// @param rows unsigned integer value that represents the desired degree of the polynomial in the 1st dimension
/// @param columns unsigned integer value that represents the desired degree of the polynomial in the 2nd dimension
bool Polynomial2D::checkCoefficients(const Eigen::MatrixXd& coefficients, const unsigned int rows, const unsigned int columns) {
    if (static_cast<size_t>(coefficients.rows()) == rows) {
        if (static_cast<size_t>(coefficients.cols()) == columns) {
            return true;
        } else {
            throw ValueError(format("%s (%d): The number of columns %d does not match with %d. ", __FILE__, __LINE__, coefficients.cols(), columns));
        }
    } else {
        throw ValueError(format("%s (%d): The number of rows %d does not match with %d. ", __FILE__, __LINE__, coefficients.rows(), rows));
    }
    return false;
}

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
/// @param axis integer value that represents the desired direction of integration
/// @param times integer value that represents the desired order of integration
Eigen::MatrixXd Polynomial2D::integrateCoeffs(const Eigen::MatrixXd& coefficients, const int& axis = -1, const int& times = 1) {
    if (times < 0)
        throw ValueError(format("%s (%d): You have to provide a positive order for integration, %d is not valid. ", __FILE__, __LINE__, times));
    if (times == 0) return Eigen::MatrixXd(coefficients);
    Eigen::MatrixXd oldCoefficients;
    Eigen::MatrixXd newCoefficients(coefficients);

    switch (axis) {
        case 0:
            newCoefficients = Eigen::MatrixXd(coefficients);
            break;
        case 1:
            newCoefficients = Eigen::MatrixXd(coefficients.transpose());
            break;
        default:
            throw ValueError(
              format("%s (%d): You have to provide a dimension, 0 or 1, for integration, %d is not valid. ", __FILE__, __LINE__, axis));
            break;
    }

    std::size_t r, c;
    for (int k = 0; k < times; k++) {
        oldCoefficients = Eigen::MatrixXd(newCoefficients);
        r = oldCoefficients.rows(), c = oldCoefficients.cols();
        newCoefficients = Eigen::MatrixXd::Zero(r + 1, c);
        newCoefficients.block(1, 0, r, c) = oldCoefficients.block(0, 0, r, c);
        for (size_t i = 0; i < r; ++i) {
            for (size_t j = 0; j < c; ++j)
                newCoefficients(i + 1, j) /= (i + 1.);
        }
    }

    switch (axis) {
        case 0:
            break;
        case 1:
            newCoefficients.transposeInPlace();
            break;
        default:
            throw ValueError(
              format("%s (%d): You have to provide a dimension, 0 or 1, for integration, %d is not valid. ", __FILE__, __LINE__, axis));
            break;
    }

    return newCoefficients;
}

/// Derivative coefficients calculation
/** Deriving coefficients for polynomials is done by multiplying the
 *  original coefficients with i and lowering the order by 1.
 */
/// @param coefficients matrix containing the ordered coefficients
/// @param axis integer value that represents the desired direction of derivation
/// @param times integer value that represents the desired order of integration
Eigen::MatrixXd Polynomial2D::deriveCoeffs(const Eigen::MatrixXd& coefficients, const int& axis, const int& times) {
    if (times < 0)
        throw ValueError(format("%s (%d): You have to provide a positive order for derivation, %d is not valid. ", __FILE__, __LINE__, times));
    if (times == 0) return Eigen::MatrixXd(coefficients);
    // Recursion is also possible, but not recommended
    //Eigen::MatrixXd newCoefficients;
    //if (times > 1) newCoefficients = deriveCoeffs(coefficients, axis, times-1);
    //else newCoefficients = Eigen::MatrixXd(coefficients);
    Eigen::MatrixXd newCoefficients;

    switch (axis) {
        case 0:
            newCoefficients = Eigen::MatrixXd(coefficients);
            break;
        case 1:
            newCoefficients = Eigen::MatrixXd(coefficients.transpose());
            break;
        default:
            throw ValueError(
              format("%s (%d): You have to provide a dimension, 0 or 1, for integration, %d is not valid. ", __FILE__, __LINE__, axis));
            break;
    }

    std::size_t r, c, i, j;
    for (int k = 0; k < times; k++) {
        r = newCoefficients.rows(), c = newCoefficients.cols();
        for (i = 1; i < r; ++i) {
            for (j = 0; j < c; ++j)
                newCoefficients(i, j) *= i;
        }
        removeRow(newCoefficients, 0);
    }

    switch (axis) {
        case 0:
            break;
        case 1:
            newCoefficients.transposeInPlace();
            break;
        default:
            throw ValueError(
              format("%s (%d): You have to provide a dimension, 0 or 1, for integration, %d is not valid. ", __FILE__, __LINE__, axis));
            break;
    }

    return newCoefficients;
}

/// The core functions to evaluate the polynomial
/** It is here we implement the different special
 *  functions that allow us to specify certain
 *  types of polynomials.
 *  The derivative might bee needed during the
 *  solution process of the solver. It could also
 *  be a protected function...
 */
/// @param coefficients vector containing the ordered coefficients
/// @param x_in double value that represents the current input
double Polynomial2D::evaluate(const Eigen::MatrixXd& coefficients, const double& x_in) {
    double result = Eigen::poly_eval(makeVector(coefficients), x_in);
    if (this->do_debug())
        std::cout << "Running      1D evaluate(" << mat_to_string(coefficients) << ", x_in:" << vec_to_string(x_in) << "): " << result << std::endl;
    return result;
}
/// @param coefficients vector containing the ordered coefficients
/// @param x_in double value that represents the current input in the 1st dimension
/// @param y_in double value that represents the current input in the 2nd dimension
double Polynomial2D::evaluate(const Eigen::MatrixXd& coefficients, const double& x_in, const double& y_in) {
    size_t r = coefficients.rows();
    double result = evaluate(coefficients.row(r - 1), y_in);
    for (int i = static_cast<int>(r) - 2; i >= 0; i--) {
        result *= x_in;
        result += evaluate(coefficients.row(i), y_in);
    }
    if (this->do_debug())
        std::cout << "Running      2D evaluate(" << mat_to_string(coefficients) << ", x_in:" << vec_to_string(x_in)
                  << ", y_in:" << vec_to_string(y_in) << "): " << result << std::endl;
    return result;
}

/// @param coefficients vector containing the ordered coefficients
/// @param x_in double value that represents the current input in the 1st dimension
/// @param y_in double value that represents the current input in the 2nd dimension
/// @param axis unsigned integer value that represents the axis to derive for (0=x, 1=y)
double Polynomial2D::derivative(const Eigen::MatrixXd& coefficients, const double& x_in, const double& y_in, const int& axis = -1) {
    return this->evaluate(this->deriveCoeffs(coefficients, axis, 1), x_in, y_in);
}

/// @param coefficients vector containing the ordered coefficients
/// @param x_in double value that represents the current input in the 1st dimension
/// @param y_in double value that represents the current input in the 2nd dimension
/// @param axis unsigned integer value that represents the axis to integrate for (0=x, 1=y)
double Polynomial2D::integral(const Eigen::MatrixXd& coefficients, const double& x_in, const double& y_in, const int& axis = -1) {
    return this->evaluate(this->integrateCoeffs(coefficients, axis, 1), x_in, y_in);
}

/// Uses the Brent solver to find the roots of p(x_in,y_in)-z_in
/// @param res Poly2DResidual object to calculate residuals and derivatives
/// @param min double value that represents the minimum value
/// @param max double value that represents the maximum value
double Polynomial2D::solve_limits(Poly2DResidual* res, const double& min, const double& max) {
    if (do_debug()) std::cout << format("Called solve_limits with: min=%f and max=%f", min, max) << std::endl;
    double macheps = DBL_EPSILON;
    double tol = DBL_EPSILON * 1e3;
    int maxiter = 10;
    double result = Brent(res, min, max, macheps, tol, maxiter);
    if (this->do_debug()) std::cout << "Brent solver message: " << res->errstring << std::endl;
    return result;
}

/// Uses the Newton solver to find the roots of p(x_in,y_in)-z_in
/// @param res Poly2DResidual object to calculate residuals and derivatives
/// @param guess double value that represents the start value
double Polynomial2D::solve_guess(Poly2DResidual* res, const double& guess) {
    if (do_debug()) std::cout << format("Called solve_guess with: guess=%f ", guess) << std::endl;
    //set_debug_level(1000);
    double tol = DBL_EPSILON * 1e3;
    int maxiter = 10;
    double result = Newton(res, guess, tol, maxiter);
    if (this->do_debug()) std::cout << "Newton solver message: " << res->errstring << std::endl;
    return result;
}

/// @param coefficients vector containing the ordered coefficients
/// @param in double value that represents the current input in x (1st dimension) or y (2nd dimension)
/// @param z_in double value that represents the current output in the 3rd dimension
/// @param axis unsigned integer value that represents the axis to solve for (0=x, 1=y)
Eigen::VectorXd Polynomial2D::solve(const Eigen::MatrixXd& coefficients, const double& in, const double& z_in, const int& axis = -1) {
    std::size_t r = coefficients.rows(), c = coefficients.cols();
    Eigen::VectorXd tmpCoefficients;
    switch (axis) {
        case 0:
            tmpCoefficients = Eigen::VectorXd::Zero(r);
            for (size_t i = 0; i < r; i++) {
                tmpCoefficients(i, 0) = evaluate(coefficients.row(i), in);
            }
            break;
        case 1:
            tmpCoefficients = Eigen::VectorXd::Zero(c);
            for (size_t i = 0; i < c; i++) {
                tmpCoefficients(i, 0) = evaluate(coefficients.col(i), in);
            }
            break;
        default:
            throw ValueError(format("%s (%d): You have to provide a dimension, 0 or 1, for the solver, %d is not valid. ", __FILE__, __LINE__, axis));
            break;
    }
    tmpCoefficients(0, 0) -= z_in;
    if (this->do_debug()) std::cout << "Coefficients: " << mat_to_string(Eigen::MatrixXd(tmpCoefficients)) << std::endl;
    Eigen::PolynomialSolver<double, Eigen::Dynamic> polySolver(tmpCoefficients);
    std::vector<double> rootsVec;
    polySolver.realRoots(rootsVec);
    if (this->do_debug()) std::cout << "Real roots: " << vec_to_string(rootsVec) << std::endl;
    return vec_to_eigen(rootsVec);
    //return rootsVec[0]; // TODO: implement root selection algorithm
}

/// @param in double value that represents the current input in x (1st dimension) or y (2nd dimension)
/// @param z_in double value that represents the current output in the 3rd dimension
/// @param axis unsigned integer value that represents the axis to solve for (0=x, 1=y)
double Polynomial2D::solve_limits(const Eigen::MatrixXd& coefficients, const double& in, const double& z_in, const double& min, const double& max,
                                  const int& axis) {
    Poly2DResidual res = Poly2DResidual(*this, coefficients, in, z_in, axis);
    return solve_limits(&res, min, max);
}

/// @param in double value that represents the current input in x (1st dimension) or y (2nd dimension)
/// @param z_in double value that represents the current output in the 3rd dimension
/// @param guess double value that represents the start value
/// @param axis unsigned integer value that represents the axis to solve for (0=x, 1=y)
double Polynomial2D::solve_guess(const Eigen::MatrixXd& coefficients, const double& in, const double& z_in, const double& guess, const int& axis) {
    Poly2DResidual res = Poly2DResidual(*this, coefficients, in, z_in, axis);
    return solve_guess(&res, guess);
}

/// Simple polynomial function generator. <- Deprecated due to poor performance, use Horner-scheme instead
/** Base function to produce n-th order polynomials
 *  based on the length of the coefficient vector.
 *  Starts with only the first coefficient at x^0. */
double Polynomial2D::simplePolynomial(std::vector<double> const& coefficients, double x) {
    double result = 0.;
    for (unsigned int i = 0; i < coefficients.size(); i++) {
        result += coefficients[i] * pow(x, (int)i);
    }
    if (this->do_debug())
        std::cout << "Running simplePolynomial(" << vec_to_string(coefficients) << ", " << vec_to_string(x) << "): " << result << std::endl;
    return result;
}
double Polynomial2D::simplePolynomial(std::vector<std::vector<double>> const& coefficients, double x, double y) {
    double result = 0;
    for (unsigned int i = 0; i < coefficients.size(); i++) {
        result += pow(x, (int)i) * simplePolynomial(coefficients[i], y);
    }
    if (this->do_debug())
        std::cout << "Running simplePolynomial(" << vec_to_string(coefficients) << ", " << vec_to_string(x) << ", " << vec_to_string(y)
                  << "): " << result << std::endl;
    return result;
}

/// Horner function generator implementations
/** Represent polynomials according to Horner's scheme.
 *  This avoids unnecessary multiplication and thus
 *  speeds up calculation.
 */
double Polynomial2D::baseHorner(std::vector<double> const& coefficients, double x) {
    double result = 0;
    for (int i = static_cast<int>(coefficients.size()) - 1; i >= 0; i--) {
        result *= x;
        result += coefficients[i];
    }
    if (this->do_debug())
        std::cout << "Running       baseHorner(" << vec_to_string(coefficients) << ", " << vec_to_string(x) << "): " << result << std::endl;
    return result;
}

double Polynomial2D::baseHorner(std::vector<std::vector<double>> const& coefficients, double x, double y) {
    double result = 0;
    for (int i = static_cast<int>(coefficients.size() - 1); i >= 0; i--) {
        result *= x;
        result += baseHorner(coefficients[i], y);
    }
    if (this->do_debug())
        std::cout << "Running       baseHorner(" << vec_to_string(coefficients) << ", " << vec_to_string(x) << ", " << vec_to_string(y)
                  << "): " << result << std::endl;
    return result;
}

Poly2DResidual::Poly2DResidual(Polynomial2D& poly, const Eigen::MatrixXd& coefficients, const double& in, const double& z_in, const int& axis) {
    switch (axis) {
        case iX:
        case iY:
            this->axis = axis;
            break;
        default:
            throw ValueError(format("%s (%d): You have to provide a dimension to the solver, %d is not valid. ", __FILE__, __LINE__, axis));
            break;
    }

    this->poly = poly;
    this->coefficients = coefficients;
    this->derIsSet = false;
    this->in = in;
    this->z_in = z_in;
}

double Poly2DResidual::call(double target) {
    if (axis == iX) return poly.evaluate(coefficients, target, in) - z_in;
    if (axis == iY) return poly.evaluate(coefficients, in, target) - z_in;
    return -_HUGE;
}

double Poly2DResidual::deriv(double target) {
    if (!this->derIsSet) {
        this->coefficientsDer = poly.deriveCoeffs(coefficients, axis);
        this->derIsSet = true;
    }
    return poly.evaluate(coefficientsDer, target, in);
}

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
 *
 */
/// @param coefficients matrix containing the ordered coefficients
/// @param axis unsigned integer value that represents the desired direction of derivation
/// @param times integer value that represents the desired order of derivation
/// @param firstExponent integer value that represents the lowest exponent of the polynomial in axis direction
Eigen::MatrixXd Polynomial2DFrac::deriveCoeffs(const Eigen::MatrixXd& coefficients, const int& axis, const int& times, const int& firstExponent) {
    if (times < 0)
        throw ValueError(format("%s (%d): You have to provide a positive order for derivation, %d is not valid. ", __FILE__, __LINE__, times));
    if (times == 0) return Eigen::MatrixXd(coefficients);
    // Recursion is also possible, but not recommended
    //Eigen::MatrixXd newCoefficients;
    //if (times > 1) newCoefficients = deriveCoeffs(coefficients, axis, times-1);
    //else newCoefficients = Eigen::MatrixXd(coefficients);
    Eigen::MatrixXd newCoefficients;

    switch (axis) {
        case 0:
            newCoefficients = Eigen::MatrixXd(coefficients);
            break;
        case 1:
            newCoefficients = Eigen::MatrixXd(coefficients.transpose());
            break;
        default:
            throw ValueError(
              format("%s (%d): You have to provide a dimension, 0 or 1, for integration, %d is not valid. ", __FILE__, __LINE__, axis));
            break;
    }

    std::size_t r = newCoefficients.rows(), c = newCoefficients.cols();
    std::size_t i, j;
    for (int k = 0; k < times; k++) {
        for (i = 0; i < r; ++i) {
            for (j = 0; j < c; ++j) {
                newCoefficients(i, j) *= double(i) + double(firstExponent);
            }
        }
    }

    switch (axis) {
        case 0:
            break;
        case 1:
            newCoefficients.transposeInPlace();
            break;
        default:
            throw ValueError(
              format("%s (%d): You have to provide a dimension, 0 or 1, for integration, %d is not valid. ", __FILE__, __LINE__, axis));
            break;
    }

    return newCoefficients;
}

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
/// @param x_base double value that represents the base value for a centred fit in the 1st dimension
double Polynomial2DFrac::evaluate(const Eigen::MatrixXd& coefficients, const double& x_in, const int& firstExponent, const double& x_base) {
    size_t r = coefficients.rows();
    size_t c = coefficients.cols();

    if ((r != 1) && (c != 1)) {
        throw ValueError(format("%s (%d): You have a 2D coefficient matrix (%d,%d), please use the 2D functions. ", __FILE__, __LINE__,
                                coefficients.rows(), coefficients.cols()));
    }
    if ((firstExponent < 0) && (std::abs(x_in - x_base) < DBL_EPSILON)) {
        throw ValueError(
          format("%s (%d): A fraction cannot be evaluated with zero as denominator, x_in-x_base=%f ", __FILE__, __LINE__, x_in - x_base));
    }

    Eigen::MatrixXd tmpCoeffs(coefficients);
    if (c == 1) {
        tmpCoeffs.transposeInPlace();
        c = r;
        r = 1;
    }
    Eigen::MatrixXd newCoeffs;
    double negExp = 0;  // First we treat the negative exponents
    double posExp = 0;  // then the positive exponents

    for (int i = 0; i > firstExponent; i--) {  // only for firstExponent<0
        if (c > 0) {
            negExp += tmpCoeffs(0, 0);
            removeColumn(tmpCoeffs, 0);
        }
        negExp /= x_in - x_base;
        c = tmpCoeffs.cols();
    }

    for (int i = 0; i < firstExponent; i++) {  // only for firstExponent>0
        newCoeffs = Eigen::MatrixXd::Zero(r, c + 1);
        newCoeffs.block(0, 1, r, c) = tmpCoeffs.block(0, 0, r, c);
        tmpCoeffs = Eigen::MatrixXd(newCoeffs);
        c = tmpCoeffs.cols();
    }

    if (c > 0) posExp += Eigen::poly_eval(Eigen::RowVectorXd(tmpCoeffs), x_in - x_base);
    return negExp + posExp;
}

/// @param coefficients vector containing the ordered coefficients
/// @param x_in double value that represents the current input in the 1st dimension
/// @param y_in double value that represents the current input in the 2nd dimension
/// @param x_exp integer value that represents the lowest exponent of the polynomial in the 1st dimension
/// @param y_exp integer value that represents the lowest exponent of the polynomial in the 2nd dimension
/// @param x_base double value that represents the base value for a centred fit in the 1st dimension
/// @param y_base double value that represents the base value for a centred fit in the 2nd dimension
double Polynomial2DFrac::evaluate(const Eigen::MatrixXd& coefficients, const double& x_in, const double& y_in, const int& x_exp, const int& y_exp,
                                  const double& x_base, const double& y_base) {
    if ((x_exp < 0) && (std::abs(x_in - x_base) < DBL_EPSILON)) {
        throw ValueError(
          format("%s (%d): A fraction cannot be evaluated with zero as denominator, x_in-x_base=%f ", __FILE__, __LINE__, x_in - x_base));
    }
    if ((y_exp < 0) && (std::abs(y_in - y_base) < DBL_EPSILON)) {
        throw ValueError(
          format("%s (%d): A fraction cannot be evaluated with zero as denominator, y_in-y_base=%f ", __FILE__, __LINE__, y_in - y_base));
    }

    Eigen::MatrixXd tmpCoeffs(coefficients);
    Eigen::MatrixXd newCoeffs;
    size_t r = tmpCoeffs.rows();
    size_t c = tmpCoeffs.cols();
    double negExp = 0;  // First we treat the negative exponents
    double posExp = 0;  // then the positive exponents

    for (int i = 0; i > x_exp; i--) {  // only for x_exp<0
        r = tmpCoeffs.rows();
        if (r > 0) {
            negExp += evaluate(tmpCoeffs.row(0), y_in, y_exp, y_base);
            removeRow(tmpCoeffs, 0);
        }
        negExp /= x_in - x_base;
    }

    r = tmpCoeffs.rows();
    for (int i = 0; i < x_exp; i++) {  // only for x_exp>0
        newCoeffs = Eigen::MatrixXd::Zero(r + 1, c);
        newCoeffs.block(1, 0, r, c) = tmpCoeffs.block(0, 0, r, c);
        tmpCoeffs = Eigen::MatrixXd(newCoeffs);
        r += 1;  // r = tmpCoeffs.rows();
    }

    //r = tmpCoeffs.rows();
    if (r > 0) posExp += evaluate(tmpCoeffs.row(r - 1), y_in, y_exp, y_base);
    for (int i = static_cast<int>(r) - 2; i >= 0; i--) {
        posExp *= x_in - x_base;
        posExp += evaluate(tmpCoeffs.row(i), y_in, y_exp, y_base);
    }
    if (this->do_debug()) std::cout << "Running      2D evaluate(" << mat_to_string(coefficients) << ", " << std::endl;
    if (this->do_debug())
        std::cout << "x_in:" << vec_to_string(x_in) << ", y_in:" << vec_to_string(y_in) << ", x_exp:" << vec_to_string(x_exp)
                  << ", y_exp:" << vec_to_string(y_exp) << ", x_base:" << vec_to_string(x_base) << ", y_base:" << vec_to_string(y_base)
                  << "): " << negExp + posExp << std::endl;
    return negExp + posExp;
}

/// @param coefficients vector containing the ordered coefficients
/// @param x_in double value that represents the current input in the 1st dimension
/// @param y_in double value that represents the current input in the 2nd dimension
/// @param axis integer value that represents the axis to derive for (0=x, 1=y)
/// @param x_exp integer value that represents the lowest exponent of the polynomial in the 1st dimension
/// @param y_exp integer value that represents the lowest exponent of the polynomial in the 2nd dimension
/// @param x_base double value that represents the base value for a centred fit in the 1st dimension
/// @param y_base double value that represents the base value for a centred fit in the 2nd dimension
double Polynomial2DFrac::derivative(const Eigen::MatrixXd& coefficients, const double& x_in, const double& y_in, const int& axis, const int& x_exp,
                                    const int& y_exp, const double& x_base, const double& y_base) {
    Eigen::MatrixXd newCoefficients;
    int der_exp, other_exp;
    double der_val, other_val;
    double int_base, other_base;

    switch (axis) {
        case 0:
            newCoefficients = Eigen::MatrixXd(coefficients);
            der_exp = x_exp;
            other_exp = y_exp;
            der_val = x_in;
            other_val = y_in;
            int_base = x_base;
            other_base = y_base;
            break;
        case 1:
            newCoefficients = Eigen::MatrixXd(coefficients.transpose());
            der_exp = y_exp;
            other_exp = x_exp;
            der_val = y_in;
            other_val = x_in;
            int_base = y_base;
            other_base = x_base;
            break;
        default:
            throw ValueError(
              format("%s (%d): You have to provide a dimension, 0 or 1, for integration, %d is not valid. ", __FILE__, __LINE__, axis));
            break;
    }

    const int times = 1;
    newCoefficients = deriveCoeffs(newCoefficients, 0, times, der_exp);
    der_exp -= times;

    return evaluate(newCoefficients, der_val, other_val, der_exp, other_exp, int_base, other_base);
}

/// @param coefficients vector containing the ordered coefficients
/// @param x_in double value that represents the current input in the 1st dimension
/// @param y_in double value that represents the current input in the 2nd dimension
/// @param axis integer value that represents the axis to integrate for (0=x, 1=y)
/// @param x_exp integer value that represents the lowest exponent of the polynomial in the 1st dimension
/// @param y_exp integer value that represents the lowest exponent of the polynomial in the 2nd dimension
/// @param x_base double value that represents the base value for a centred fit in the 1st dimension
/// @param y_base double value that represents the base value for a centred fit in the 2nd dimension
/// @param ax_val double value that represents the base value for the definite integral on the chosen axis.
double Polynomial2DFrac::integral(const Eigen::MatrixXd& coefficients, const double& x_in, const double& y_in, const int& axis, const int& x_exp,
                                  const int& y_exp, const double& x_base, const double& y_base, const double& ax_val) {

    Eigen::MatrixXd newCoefficients;
    int int_exp, other_exp;
    double int_val, other_val;
    double int_base, other_base;

    switch (axis) {
        case 0:
            newCoefficients = Eigen::MatrixXd(coefficients);
            int_exp = x_exp;
            other_exp = y_exp;
            int_val = x_in;
            other_val = y_in;
            int_base = x_base;
            other_base = y_base;
            break;
        case 1:
            newCoefficients = Eigen::MatrixXd(coefficients.transpose());
            int_exp = y_exp;
            other_exp = x_exp;
            int_val = y_in;
            other_val = x_in;
            int_base = y_base;
            other_base = x_base;
            break;
        default:
            throw ValueError(
              format("%s (%d): You have to provide a dimension, 0 or 1, for integration, %d is not valid. ", __FILE__, __LINE__, axis));
            break;
    }

    if (int_exp < -1)
        throw NotImplementedError(
          format("%s (%d): This function is only implemented for lowest exponents >= -1, int_exp=%d ", __FILE__, __LINE__, int_exp));
    // TODO: Fix this and allow the direct calculation of integrals
    if (std::abs(ax_val) > DBL_EPSILON)
        throw NotImplementedError(
          format("%s (%d): This function is only implemented for indefinite integrals, ax_val=%d ", __FILE__, __LINE__, ax_val));

    size_t r = newCoefficients.rows();
    size_t c = newCoefficients.cols();

    if (int_exp == -1) {
        if (std::abs(int_base) < DBL_EPSILON) {
            Eigen::MatrixXd tmpCoefficients = newCoefficients.row(0) * log(int_val - int_base);
            newCoefficients = integrateCoeffs(newCoefficients.block(1, 0, r - 1, c), 0, 1);
            newCoefficients.row(0) = tmpCoefficients;
            return evaluate(newCoefficients, int_val, other_val, 0, other_exp, int_base, other_base);
        } else {
            // Reduce the coefficients to the integration dimension:
            newCoefficients = Eigen::MatrixXd(r, 1);
            for (std::size_t i = 0; i < r; i++) {
                newCoefficients(i, 0) = evaluate(coefficients.row(i), other_val, other_exp, other_base);
            }
            return fracIntCentral(newCoefficients.transpose(), int_val, int_base);
        }
    }

    Eigen::MatrixXd tmpCoeffs;
    r = newCoefficients.rows();
    for (int i = 0; i < int_exp; i++) {  // only for x_exp>0
        tmpCoeffs = Eigen::MatrixXd::Zero(r + 1, c);
        tmpCoeffs.block(1, 0, r, c) = newCoefficients.block(0, 0, r, c);
        newCoefficients = Eigen::MatrixXd(tmpCoeffs);
        r += 1;  // r = newCoefficients.rows();
    }

    return evaluate(integrateCoeffs(newCoefficients, 0, 1), int_val, other_val, 0, other_exp, int_base, other_base);
}

/// Returns a vector with ALL the real roots of p(x_in,y_in)-z_in
/// @param coefficients vector containing the ordered coefficients
/// @param in double value that represents the current input in x (1st dimension) or y (2nd dimension)
/// @param z_in double value that represents the current output in the 3rd dimension
/// @param axis integer value that represents the axis to solve for (0=x, 1=y)
/// @param x_exp integer value that represents the lowest exponent of the polynomial in the 1st dimension
/// @param y_exp integer value that represents the lowest exponent of the polynomial in the 2nd dimension
/// @param x_base double value that represents the base value for a centred fit in the 1st dimension
/// @param y_base double value that represents the base value for a centred fit in the 2nd dimension
Eigen::VectorXd Polynomial2DFrac::solve(const Eigen::MatrixXd& coefficients, const double& in, const double& z_in, const int& axis, const int& x_exp,
                                        const int& y_exp, const double& x_base, const double& y_base) {

    Eigen::MatrixXd newCoefficients;
    Eigen::VectorXd tmpCoefficients;
    int solve_exp, other_exp;
    double input;

    switch (axis) {
        case 0:
            newCoefficients = Eigen::MatrixXd(coefficients);
            solve_exp = x_exp;
            other_exp = y_exp;
            input = in - y_base;
            break;
        case 1:
            newCoefficients = Eigen::MatrixXd(coefficients.transpose());
            solve_exp = y_exp;
            other_exp = x_exp;
            input = in - x_base;
            break;
        default:
            throw ValueError(format("%s (%d): You have to provide a dimension, 0 or 1, for the solver, %d is not valid. ", __FILE__, __LINE__, axis));
            break;
    }

    if (this->do_debug()) std::cout << "Coefficients: " << mat_to_string(Eigen::MatrixXd(newCoefficients)) << std::endl;

    const size_t r = newCoefficients.rows();
    for (size_t i = 0; i < r; i++) {
        newCoefficients(i, 0) = evaluate(newCoefficients.row(i), input, other_exp);
    }

    //Eigen::VectorXd tmpCoefficients;
    if (solve_exp >= 0) {  // extend vector to zero exponent
        tmpCoefficients = Eigen::VectorXd::Zero(r + solve_exp);
        tmpCoefficients.block(solve_exp, 0, r, 1) = newCoefficients.block(0, 0, r, 1);
        tmpCoefficients(0, 0) -= z_in;
    } else {                                             // check if vector reaches to zero exponent
        int diff = 1 - static_cast<int>(r) - solve_exp;  // How many entries have to be added
        tmpCoefficients = Eigen::VectorXd::Zero(r + std::max(diff, 0));
        tmpCoefficients.block(0, 0, r, 1) = newCoefficients.block(0, 0, r, 1);
        tmpCoefficients(r + diff - 1, 0) -= z_in;
        throw NotImplementedError(format("%s (%d): Currently, there is no solver for the generalised polynomial, an exponent of %d is not valid. ",
                                         __FILE__, __LINE__, solve_exp));
    }

    if (this->do_debug()) std::cout << "Coefficients: " << mat_to_string(Eigen::MatrixXd(tmpCoefficients)) << std::endl;
    Eigen::PolynomialSolver<double, -1> polySolver(tmpCoefficients);
    std::vector<double> rootsVec;
    polySolver.realRoots(rootsVec);
    if (this->do_debug()) std::cout << "Real roots: " << vec_to_string(rootsVec) << std::endl;
    return vec_to_eigen(rootsVec);
    //return rootsVec[0]; // TODO: implement root selection algorithm
}

/// Uses the Brent solver to find the roots of p(x_in,y_in)-z_in
/// @param coefficients vector containing the ordered coefficients
/// @param in double value that represents the current input in x (1st dimension) or y (2nd dimension)
/// @param z_in double value that represents the current output in the 3rd dimension
/// @param min double value that represents the minimum value
/// @param max double value that represents the maximum value
/// @param axis integer value that represents the axis to solve for (0=x, 1=y)
/// @param x_exp integer value that represents the lowest exponent of the polynomial in the 1st dimension
/// @param y_exp integer value that represents the lowest exponent of the polynomial in the 2nd dimension
/// @param x_base double value that represents the base value for a centred fit in the 1st dimension
/// @param y_base double value that represents the base value for a centred fit in the 2nd dimension
double Polynomial2DFrac::solve_limits(const Eigen::MatrixXd& coefficients, const double& in, const double& z_in, const double& min, const double& max,
                                      const int& axis, const int& x_exp, const int& y_exp, const double& x_base, const double& y_base) {
    if (do_debug())
        std::cout << format("Called solve_limits with: %f, %f, %f, %f, %d, %d, %d, %f, %f", in, z_in, min, max, axis, x_exp, y_exp, x_base, y_base)
                  << std::endl;
    Poly2DFracResidual res = Poly2DFracResidual(*this, coefficients, in, z_in, axis, x_exp, y_exp, x_base, y_base);
    return Polynomial2D::solve_limits(&res, min, max);
}  //TODO: Implement tests for this solver

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
double Polynomial2DFrac::solve_guess(const Eigen::MatrixXd& coefficients, const double& in, const double& z_in, const double& guess, const int& axis,
                                     const int& x_exp, const int& y_exp, const double& x_base, const double& y_base) {
    if (do_debug())
        std::cout << format("Called solve_guess with: %f, %f, %f, %d, %d, %d, %f, %f", in, z_in, guess, axis, x_exp, y_exp, x_base, y_base)
                  << std::endl;
    Poly2DFracResidual res = Poly2DFracResidual(*this, coefficients, in, z_in, axis, x_exp, y_exp, x_base, y_base);
    return Polynomial2D::solve_guess(&res, guess);
}  //TODO: Implement tests for this solver

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
double Polynomial2DFrac::solve_limitsInt(const Eigen::MatrixXd& coefficients, const double& in, const double& z_in, const double& min,
                                         const double& max, const int& axis, const int& x_exp, const int& y_exp, const double& x_base,
                                         const double& y_base, const int& int_axis) {
    Poly2DFracIntResidual res = Poly2DFracIntResidual(*this, coefficients, in, z_in, axis, x_exp, y_exp, x_base, y_base, int_axis);
    return Polynomial2D::solve_limits(&res, min, max);
}  //TODO: Implement tests for this solver

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
double Polynomial2DFrac::solve_guessInt(const Eigen::MatrixXd& coefficients, const double& in, const double& z_in, const double& guess,
                                        const int& axis, const int& x_exp, const int& y_exp, const double& x_base, const double& y_base,
                                        const int& int_axis) {
    Poly2DFracIntResidual res = Poly2DFracIntResidual(*this, coefficients, in, z_in, axis, x_exp, y_exp, x_base, y_base, int_axis);
    return Polynomial2D::solve_guess(&res, guess);
}  //TODO: Implement tests for this solver

/** Simple integrated centred(!) polynomial function generator divided by independent variable.
 *  We need to rewrite some of the functions in order to
 *  use central fit. Having a central temperature xbase
 *  allows for a better fit, but requires a different
 *  formulation of the fracInt function group. Other
 *  functions are not affected.
 *  Starts with only the first coefficient at x^0 */

//Helper functions to calculate binomial coefficients:
//http://rosettacode.org/wiki/Evaluate_binomial_coefficients#C.2B.2B
/// @param nValue integer value that represents the order of the factorial
double Polynomial2DFrac::factorial(const int& nValue) {
    double value = 1;
    for (int i = 2; i <= nValue; i++)
        value = value * i;
    return value;
}
/// @param nValue integer value that represents the upper part of the factorial
/// @param nValue2 integer value that represents the lower part of the factorial
double Polynomial2DFrac::binom(const int& nValue, const int& nValue2) {
    if (nValue2 == 1) return nValue * 1.0;
    return (factorial(nValue)) / (factorial(nValue2) * factorial((nValue - nValue2)));
}
///Helper function to calculate the D vector:
/// @param m integer value that represents order
/// @param x_in double value that represents the current input
/// @param x_base double value that represents the basis for the fit
Eigen::MatrixXd Polynomial2DFrac::fracIntCentralDvector(const int& m, const double& x_in, const double& x_base) {
    if (m < 1) throw ValueError(format("%s (%d): You have to provide coefficients, a vector length of %d is not a valid. ", __FILE__, __LINE__, m));

    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(1, m);
    double tmp;
    // TODO: This can be optimized using the Horner scheme!
    for (int j = 0; j < m; j++) {  // loop through row
        tmp = pow(-1.0, j) * log(x_in) * pow(x_base, j);
        for (int k = 0; k < j; k++) {  // internal loop for every entry
            tmp += binom(j, k) * pow(-1.0, k) / (j - k) * pow(x_in, j - k) * pow(x_base, k);
        }
        D(0, j) = tmp;
    }
    return D;
}
///Indefinite integral of a centred polynomial divided by its independent variable
/// @param coefficients vector containing the ordered coefficients
/// @param x_in double value that represents the current input
/// @param x_base double value that represents the basis for the fit
double Polynomial2DFrac::fracIntCentral(const Eigen::MatrixXd& coefficients, const double& x_in, const double& x_base) {
    if (coefficients.rows() != 1) {
        throw ValueError(format("%s (%d): You have a 2D coefficient matrix (%d,%d), please use the 2D functions. ", __FILE__, __LINE__,
                                coefficients.rows(), coefficients.cols()));
    }
    int m = static_cast<int>(coefficients.cols());
    Eigen::MatrixXd D = fracIntCentralDvector(m, x_in, x_base);
    double result = 0;
    for (int j = 0; j < m; j++) {
        result += coefficients(0, j) * D(0, j);
    }
    if (this->do_debug())
        std::cout << "Running   fracIntCentral(" << mat_to_string(coefficients) << ", " << vec_to_string(x_in) << ", " << vec_to_string(x_base)
                  << "): " << result << std::endl;
    return result;
}

Poly2DFracResidual::Poly2DFracResidual(Polynomial2DFrac& poly, const Eigen::MatrixXd& coefficients, const double& in, const double& z_in,
                                       const int& axis, const int& x_exp, const int& y_exp, const double& x_base, const double& y_base)
  : Poly2DResidual(poly, coefficients, in, z_in, axis) {
    this->x_exp = x_exp;
    this->y_exp = y_exp;
    this->x_base = x_base;
    this->y_base = y_base;
}

double Poly2DFracResidual::call(double target) {
    if (axis == iX) return poly.evaluate(coefficients, target, in, x_exp, y_exp, x_base, y_base) - z_in;
    if (axis == iY) return poly.evaluate(coefficients, in, target, x_exp, y_exp, x_base, y_base) - z_in;
    return _HUGE;
}

double Poly2DFracResidual::deriv(double target) {
    if (axis == iX) return poly.derivative(coefficients, target, in, axis, x_exp, y_exp, x_base, y_base);
    if (axis == iY) return poly.derivative(coefficients, in, target, axis, x_exp, y_exp, x_base, y_base);
    return _HUGE;
}

Poly2DFracIntResidual::Poly2DFracIntResidual(Polynomial2DFrac& poly, const Eigen::MatrixXd& coefficients, const double& in, const double& z_in,
                                             const int& axis, const int& x_exp, const int& y_exp, const double& x_base, const double& y_base,
                                             const int& int_axis)
  : Poly2DFracResidual(poly, coefficients, in, z_in, axis, x_exp, y_exp, x_base, y_base) {
    this->int_axis = int_axis;
}

double Poly2DFracIntResidual::call(double target) {
    if (axis == iX) return poly.integral(coefficients, target, in, int_axis, x_exp, y_exp, x_base, y_base) - z_in;
    if (axis == iY) return poly.integral(coefficients, in, target, int_axis, x_exp, y_exp, x_base, y_base) - z_in;
    return _HUGE;
}

double Poly2DFracIntResidual::deriv(double target) {
    if (axis == iX) return poly.evaluate(coefficients, target, in, x_exp, y_exp, x_base, y_base);
    if (axis == iY) return poly.evaluate(coefficients, in, target, x_exp, y_exp, x_base, y_base);
    return _HUGE;
}

}  // namespace CoolProp

#ifdef ENABLE_CATCH
#    include <math.h>
#    include <iostream>
#    include <catch2/catch_all.hpp>

TEST_CASE("Internal consistency checks and example use cases for PolyMath.cpp", "[PolyMath]") {
    bool PRINT = false;
    std::string tmpStr;

    /// Test case for "SylthermXLT" by "Dow Chemicals"
    std::vector<double> cHeat;
    cHeat.clear();
    cHeat.push_back(+1.1562261074E+03);
    cHeat.push_back(+2.0994549103E+00);
    cHeat.push_back(+7.7175381057E-07);
    cHeat.push_back(-3.7008444051E-20);

    double deltaT = 0.1;
    double Tmin = 273.15 - 50;
    double Tmax = 273.15 + 250;
    double Tinc = 200;

    std::vector<std::vector<double>> cHeat2D;
    cHeat2D.push_back(cHeat);
    cHeat2D.push_back(cHeat);
    cHeat2D.push_back(cHeat);

    Eigen::MatrixXd matrix2D = CoolProp::vec_to_eigen(cHeat2D);

    Eigen::MatrixXd matrix2Dtmp;
    std::vector<std::vector<double>> vec2Dtmp;

    SECTION("Coefficient parsing") {
        CoolProp::Polynomial2D poly;
        CHECK_THROWS(poly.checkCoefficients(matrix2D, 4, 5));
        CHECK(poly.checkCoefficients(matrix2D, 3, 4));
    }

    SECTION("Coefficient operations") {
        Eigen::MatrixXd matrix;
        CoolProp::Polynomial2D poly;

        CHECK_THROWS(poly.integrateCoeffs(matrix2D));

        CHECK_NOTHROW(matrix = poly.integrateCoeffs(matrix2D, 0));
        tmpStr = CoolProp::mat_to_string(matrix2D);
        if (PRINT) std::cout << tmpStr << std::endl;
        tmpStr = CoolProp::mat_to_string(matrix);
        if (PRINT) std::cout << tmpStr << std::endl << std::endl;

        CHECK_NOTHROW(matrix = poly.integrateCoeffs(matrix2D, 1));
        tmpStr = CoolProp::mat_to_string(matrix2D);
        if (PRINT) std::cout << tmpStr << std::endl;
        tmpStr = CoolProp::mat_to_string(matrix);
        if (PRINT) std::cout << tmpStr << std::endl << std::endl;

        CHECK_THROWS(poly.deriveCoeffs(matrix2D));

        CHECK_NOTHROW(matrix = poly.deriveCoeffs(matrix2D, 0));
        tmpStr = CoolProp::mat_to_string(matrix2D);
        if (PRINT) std::cout << tmpStr << std::endl;
        tmpStr = CoolProp::mat_to_string(matrix);
        if (PRINT) std::cout << tmpStr << std::endl << std::endl;

        CHECK_NOTHROW(matrix = poly.deriveCoeffs(matrix2D, 1));
        tmpStr = CoolProp::mat_to_string(matrix2D);
        if (PRINT) std::cout << tmpStr << std::endl;
        tmpStr = CoolProp::mat_to_string(matrix);
        if (PRINT) std::cout << tmpStr << std::endl << std::endl;
    }

    SECTION("Evaluation and test values") {

        Eigen::MatrixXd matrix = CoolProp::vec_to_eigen(cHeat);
        CoolProp::Polynomial2D poly;

        double acc = 0.0001;

        double T = 273.15 + 50;
        double c = poly.evaluate(matrix, T, 0.0);
        double d = 1834.746;

        {
            CAPTURE(T);
            CAPTURE(c);
            CAPTURE(d);
            tmpStr = CoolProp::mat_to_string(matrix);
            CAPTURE(tmpStr);
            CHECK(check_abs(c, d, acc));
        }

        c = 2.0;
        c = poly.solve(matrix, 0.0, d, 0)[0];
        {
            CAPTURE(T);
            CAPTURE(c);
            CAPTURE(d);
            CHECK(check_abs(c, T, acc));
        }

        c = 2.0;
        c = poly.solve_limits(matrix, 0.0, d, -50, 750, 0);
        {
            CAPTURE(T);
            CAPTURE(c);
            CAPTURE(d);
            CHECK(check_abs(c, T, acc));
        }

        c = 2.0;
        c = poly.solve_guess(matrix, 0.0, d, 350, 0);
        {
            CAPTURE(T);
            CAPTURE(c);
            CAPTURE(d);
            CHECK(check_abs(c, T, acc));
        }

        //        T = 0.0;
        //        solve.setGuess(75+273.15);
        //        T = solve.polyval(cHeat,c);
        //        printf("Should be  :      T = %3.3f \t K    \n",273.15+50.0);
        //        printf("From object:      T = %3.3f \t K    \n",T);
        //
        //        T = 0.0;
        //        solve.setLimits(273.15+10,273.15+100);
        //        T = solve.polyval(cHeat,c);
        //        printf("Should be  :      T = %3.3f \t K    \n",273.15+50.0);
        //        printf("From object:      T = %3.3f \t K    \n",T);
    }

    SECTION("Integration and derivation tests") {

        CoolProp::Polynomial2D poly;

        Eigen::MatrixXd matrix(matrix2D);
        Eigen::MatrixXd matrixInt = poly.integrateCoeffs(matrix, 1);
        Eigen::MatrixXd matrixDer = poly.deriveCoeffs(matrix, 1);
        Eigen::MatrixXd matrixInt2 = poly.integrateCoeffs(matrix, 1, 2);
        Eigen::MatrixXd matrixDer2 = poly.deriveCoeffs(matrix, 1, 2);

        CHECK_THROWS(poly.evaluate(matrix, 0.0));

        double x = 0.3, y = 255.3, val1, val2, val3, val4;

        //CHECK( std::abs( polyInt.derivative(x,y,0)-poly2D.evaluate(x,y) ) <= 1e-10 );

        std::string tmpStr;

        double acc = 0.001;

        for (double T = Tmin; T < Tmax; T += Tinc) {
            val1 = poly.evaluate(matrix, x, T - deltaT);
            val2 = poly.evaluate(matrix, x, T + deltaT);
            val3 = (val2 - val1) / 2 / deltaT;
            val4 = poly.evaluate(matrixDer, x, T);
            CAPTURE(T);
            CAPTURE(val3);
            CAPTURE(val4);
            tmpStr = CoolProp::mat_to_string(matrix);
            CAPTURE(tmpStr);
            tmpStr = CoolProp::mat_to_string(matrixDer);
            CAPTURE(tmpStr);
            CHECK(check_abs(val3, val4, acc));
        }

        for (double T = Tmin; T < Tmax; T += Tinc) {
            val1 = poly.evaluate(matrixDer, x, T - deltaT);
            val2 = poly.evaluate(matrixDer, x, T + deltaT);
            val3 = (val2 - val1) / 2 / deltaT;
            val4 = poly.evaluate(matrixDer2, x, T);
            CAPTURE(T);
            CAPTURE(val3);
            CAPTURE(val4);
            tmpStr = CoolProp::mat_to_string(matrixDer);
            CAPTURE(tmpStr);
            tmpStr = CoolProp::mat_to_string(matrixDer2);
            CAPTURE(tmpStr);
            CHECK(check_abs(val3, val4, acc));
        }

        for (double T = Tmin; T < Tmax; T += Tinc) {
            val1 = poly.evaluate(matrixInt, x, T - deltaT);
            val2 = poly.evaluate(matrixInt, x, T + deltaT);
            val3 = (val2 - val1) / 2 / deltaT;
            val4 = poly.evaluate(matrix, x, T);
            CAPTURE(T);
            CAPTURE(val3);
            CAPTURE(val4);
            tmpStr = CoolProp::mat_to_string(matrixInt);
            CAPTURE(tmpStr);
            tmpStr = CoolProp::mat_to_string(matrix);
            CAPTURE(tmpStr);
            CHECK(check_abs(val3, val4, acc));
        }

        for (double T = Tmin; T < Tmax; T += Tinc) {
            val1 = poly.evaluate(matrixInt2, x, T - deltaT);
            val2 = poly.evaluate(matrixInt2, x, T + deltaT);
            val3 = (val2 - val1) / 2 / deltaT;
            val4 = poly.evaluate(matrixInt, x, T);
            CAPTURE(T);
            CAPTURE(val3);
            CAPTURE(val4);
            tmpStr = CoolProp::mat_to_string(matrixInt2);
            CAPTURE(tmpStr);
            tmpStr = CoolProp::mat_to_string(matrixInt);
            CAPTURE(tmpStr);
            CHECK(check_abs(val3, val4, acc));
        }

        for (double T = Tmin; T < Tmax; T += Tinc) {
            val1 = poly.evaluate(matrix, x, T);
            val2 = poly.derivative(matrixInt, x, T, 1);
            CAPTURE(T);
            CAPTURE(val1);
            CAPTURE(val2);
            CHECK(check_abs(val1, val2, acc));
        }

        for (double T = Tmin; T < Tmax; T += Tinc) {
            val1 = poly.derivative(matrix, x, T, 1);
            val2 = poly.evaluate(matrixDer, x, T);
            CAPTURE(T);
            CAPTURE(val1);
            CAPTURE(val2);
            CHECK(check_abs(val1, val2, acc));
        }
    }

    SECTION("Testing Polynomial2DFrac") {

        Eigen::MatrixXd matrix = CoolProp::vec_to_eigen(cHeat);
        CoolProp::Polynomial2D poly;
        CoolProp::Polynomial2DFrac frac;

        double acc = 0.0001;

        double T = 273.15 + 50;
        double a, b;
        double c = poly.evaluate(matrix, T, 0.0);
        double d = frac.evaluate(matrix, T, 0.0, 0, 0);

        {
            CAPTURE(T);
            CAPTURE(c);
            CAPTURE(d);
            tmpStr = CoolProp::mat_to_string(matrix);
            CAPTURE(tmpStr);
            CHECK(check_abs(c, d, acc));
        }

        c = poly.evaluate(matrix, T, 0.0) / T / T;
        d = frac.evaluate(matrix, T, 0.0, -2, 0);

        {
            CAPTURE(T);
            CAPTURE(c);
            CAPTURE(d);
            tmpStr = CoolProp::mat_to_string(matrix);
            CAPTURE(tmpStr);
            CHECK(check_abs(c, d, acc));
        }

        matrix = CoolProp::vec_to_eigen(cHeat2D);
        double y = 0.1;
        c = poly.evaluate(matrix, T, y) / T / T;
        d = frac.evaluate(matrix, T, y, -2, 0);

        {
            CAPTURE(T);
            CAPTURE(c);
            CAPTURE(d);
            tmpStr = CoolProp::mat_to_string(matrix);
            CAPTURE(tmpStr);
            CHECK(check_abs(c, d, acc));
        }

        c = poly.evaluate(matrix, T, y) / y / y;
        d = frac.evaluate(matrix, T, y, 0, -2);

        {
            CAPTURE(T);
            CAPTURE(c);
            CAPTURE(d);
            tmpStr = CoolProp::mat_to_string(matrix);
            CAPTURE(tmpStr);
            CHECK(check_abs(c, d, acc));
        }

        c = poly.evaluate(matrix, T, y) / T / T / y / y;
        d = frac.evaluate(matrix, T, y, -2, -2);

        {
            CAPTURE(T);
            CAPTURE(c);
            CAPTURE(d);
            tmpStr = CoolProp::mat_to_string(matrix);
            CAPTURE(tmpStr);
            CHECK(check_abs(c, d, acc));
        }

        c = poly.evaluate(matrix, T, y) / T / T * y * y;
        d = frac.evaluate(matrix, T, y, -2, 2);

        {
            CAPTURE(T);
            CAPTURE(c);
            CAPTURE(d);
            tmpStr = CoolProp::mat_to_string(matrix);
            CAPTURE(tmpStr);
            CHECK(check_abs(c, d, acc));
        }

        matrix = CoolProp::vec_to_eigen(cHeat);
        T = 273.15 + 50;
        c = 145.59157247249246;
        d = frac.integral(matrix, T, 0.0, 0, -1, 0) - frac.integral(matrix, 273.15 + 25, 0.0, 0, -1, 0);

        {
            CAPTURE(T);
            CAPTURE(c);
            CAPTURE(d);
            tmpStr = CoolProp::mat_to_string(matrix);
            CAPTURE(tmpStr);
            CHECK(check_abs(c, d, acc));
        }

        T = 423.15;
        c = 3460.895272;
        d = frac.integral(matrix, T, 0.0, 0, -1, 0, 348.15, 0.0);

        {
            CAPTURE(T);
            CAPTURE(c);
            CAPTURE(d);
            tmpStr = CoolProp::mat_to_string(matrix);
            CAPTURE(tmpStr);
            CHECK(check_abs(c, d, acc));
        }

        deltaT = 0.01;
        for (T = Tmin; T < Tmax; T += Tinc) {
            a = poly.evaluate(matrix, T - deltaT, y);
            b = poly.evaluate(matrix, T + deltaT, y);
            c = (b - a) / 2.0 / deltaT;
            d = frac.derivative(matrix, T, y, 0, 0, 0);
            CAPTURE(a);
            CAPTURE(b);
            CAPTURE(T);
            CAPTURE(c);
            CAPTURE(d);
            tmpStr = CoolProp::mat_to_string(matrix);
            CAPTURE(tmpStr);
            CHECK(check_abs(c, d, acc));
        }

        T = 273.15 + 150;
        c = -2.100108045;
        d = frac.derivative(matrix, T, 0.0, 0, 0, 0);
        {
            CAPTURE(T);
            CAPTURE(c);
            CAPTURE(d);
            tmpStr = CoolProp::mat_to_string(matrix);
            CAPTURE(tmpStr);
            CHECK(check_abs(c, d, acc));
        }

        c = -0.006456574589;
        d = frac.derivative(matrix, T, 0.0, 0, -1, 0);
        {
            CAPTURE(T);
            CAPTURE(c);
            CAPTURE(d);
            tmpStr = CoolProp::mat_to_string(matrix);
            CAPTURE(tmpStr);
            CHECK(check_abs(c, d, acc));
        }

        c = frac.evaluate(matrix, T, 0.0, 2, 0);
        d = frac.solve(matrix, 0.0, c, 0, 2, 0)[0];
        {
            CAPTURE(T);
            CAPTURE(c);
            CAPTURE(d);
            tmpStr = CoolProp::mat_to_string(matrix);
            CAPTURE(tmpStr);
            CHECK(check_abs(T, d, acc));
        }

        c = frac.evaluate(matrix, T, 0.0, 0, 0);
        d = frac.solve(matrix, 0.0, c, 0, 0, 0)[0];
        {
            CAPTURE(T);
            CAPTURE(c);
            CAPTURE(d);
            tmpStr = CoolProp::mat_to_string(matrix);
            CAPTURE(tmpStr);
            CHECK(check_abs(T, d, acc));
        }

        c = frac.evaluate(matrix, T, 0.0, -1, 0);
        CHECK_THROWS(d = frac.solve(matrix, 0.0, c, 0, -1, 0)[0]);
        //        {
        //        CAPTURE(T);
        //        CAPTURE(c);
        //        CAPTURE(d);
        //        tmpStr = CoolProp::mat_to_string(matrix);
        //        CAPTURE(tmpStr);
        //        CHECK( check_abs(T,d,acc) );
        //        }

        d = frac.solve_limits(matrix, 0.0, c, T - 10, T + 10, 0, -1, 0);
        {
            CAPTURE(T);
            CAPTURE(c);
            CAPTURE(d);
            tmpStr = CoolProp::mat_to_string(matrix);
            CAPTURE(tmpStr);
            CHECK(check_abs(T, d, acc));
        }

        d = frac.solve_guess(matrix, 0.0, c, T - 10, 0, -1, 0);
        {
            CAPTURE(T);
            CAPTURE(c);
            CAPTURE(d);
            tmpStr = CoolProp::mat_to_string(matrix);
            CAPTURE(tmpStr);
            CHECK(check_abs(T, d, acc));
        }

        c = -0.00004224550082;
        d = frac.derivative(matrix, T, 0.0, 0, -2, 0);
        {
            CAPTURE(T);
            CAPTURE(c);
            CAPTURE(d);
            tmpStr = CoolProp::mat_to_string(matrix);
            CAPTURE(tmpStr);
            CHECK(check_abs(c, d, acc));
        }

        c = frac.evaluate(matrix, T, 0.0, 0, 0, 0.0, 0.0);
        d = frac.solve(matrix, 0.0, c, 0, 0, 0, 0.0, 0.0)[0];
        {
            CAPTURE(T);
            CAPTURE(c);
            CAPTURE(d);
            tmpStr = CoolProp::mat_to_string(matrix);
            CAPTURE(tmpStr);
            tmpStr = CoolProp::mat_to_string(Eigen::MatrixXd(frac.solve(matrix, 0.0, c, 0, 0, 0, 250, 0.0)));
            CAPTURE(tmpStr);
            CHECK(check_abs(T, d, acc));
        }
    }
}

#endif /* ENABLE_CATCH */
