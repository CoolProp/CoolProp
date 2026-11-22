#ifndef SOLVERS_H
#define SOLVERS_H

#include <vector>
#include <string>
#include "CoolProp/Exceptions.h"
#include "CoolProp/detail/tools.h"

namespace CoolProp {

// *****************************************************************************
// *****************************************************************************
//                              SOLVER WRAPPER CLASSES
// *****************************************************************************
// *****************************************************************************

class FuncWrapper1D
{
   public:
    int errcode;
    std::string errstring;
    Dictionary options;
    int iter;
    int verbosity = 0;
    FuncWrapper1D() : errcode(0), errstring(""), iter(0) {};
    virtual ~FuncWrapper1D() = default;
    virtual double call(double) = 0;
    /**
     * /brief A function for checking whether the input is in range;
     *
     * Meant to be implemented by derived classes; return true if input is out of range
     */
    virtual bool input_not_in_range(double x) {
        return false;
    };
};

class FuncWrapper1DWithDeriv : public FuncWrapper1D
{
   public:
    virtual double deriv(double) = 0;
};

class FuncWrapper1DWithTwoDerivs : public FuncWrapper1DWithDeriv
{
   public:
    virtual double second_deriv(double) = 0;
};

class FuncWrapper1DWithThreeDerivs : public FuncWrapper1DWithTwoDerivs
{
   public:
    virtual double third_deriv(double) = 0;
};

class FuncWrapperND
{
   public:
    int errcode;
    std::string errstring;
    FuncWrapperND() : errcode(0), errstring("") {};
    virtual ~FuncWrapperND() = default;
    virtual std::vector<double> call(const std::vector<double>&) = 0;  // must be provided
    virtual std::vector<std::vector<double>> Jacobian(const std::vector<double>&);
};

// *****************************************************************************
// *****************************************************************************
//                              SOLVER ROUTINES
// *****************************************************************************
// *****************************************************************************

// Which solver to use
// -------------------
// Pick by what you know about the problem, not by habit:
//
// - Brent: you have an interval where f changes sign and no derivative. Robust,
//   finds a root inside the interval or raises an error, but needs the most
//   function calls. The safe default.
// - TOMS748 (boost::math::tools::toms748_solve, not in this file): the same
//   situation as Brent, usually with fewer function calls. Used in the HEOS
//   P+{H,S,U} flash.
// - BoundedNewton: an interval with a sign change AND a derivative that is
//   cheap to get. Same guarantees as Brent, a few steps less. Only pays off if
//   the derivative costs little next to f itself.
// - Newton, Halley, Householder4: a good first guess and derivatives, but no
//   interval. Fastest when they work, but they can step out of the range where
//   your equations are valid. Use with a guess you trust, and keep a bracketed
//   method as the fallback.
// - Secant, ExtrapolatingSecant, BoundedSecant: a good first guess and no
//   derivative. BoundedSecant clamps the steps to limits but does not need,
//   or check, a sign change.
// - NDNewtonRaphson_Jacobian: several unknowns at once, with a Jacobian.
//
// If in doubt, try more than one and measure. Time them on your own function
// and count the calls. In the incompressible flashes BoundedNewton was 1.2 to
// 1.9 times faster than TOMS748, but only because the derivative came almost
// for free; with a naively computed derivative TOMS748 won. The convergence
// order on paper does NOT tell you which one is faster for your problem.

// Single-Dimensional solvers, pointer versions
double Brent(FuncWrapper1D* f, double a, double b, double macheps, double t, int maxiter);
double Secant(FuncWrapper1D* f, double x0, double dx, double ftol, int maxiter);
double BoundedSecant(FuncWrapper1D* f, double x0, double xmin, double xmax, double dx, double ftol, int maxiter);
double ExtrapolatingSecant(FuncWrapper1D* f, double x0, double dx, double ftol, int maxiter);
double Newton(FuncWrapper1DWithDeriv* f, double x0, double ftol, int maxiter);
double BoundedNewton(FuncWrapper1DWithDeriv* f, double x0, double xmin, double xmax, double xtol, int maxiter);
double Halley(FuncWrapper1DWithTwoDerivs* f, double x0, double ftol, int maxiter, double xtol_rel = 1e-12);
double Householder4(FuncWrapper1DWithThreeDerivs* f, double x0, double ftol, int maxiter, double xtol_rel = 1e-12);

// Single-Dimensional solvers, refere
inline double Brent(FuncWrapper1D& f, double a, double b, double macheps, double t, int maxiter) {
    return Brent(&f, a, b, macheps, t, maxiter);
}
inline double Secant(FuncWrapper1D& f, double x0, double dx, double ftol, int maxiter) {
    return Secant(&f, x0, dx, ftol, maxiter);
}

inline double ExtrapolatingSecant(FuncWrapper1D& f, double x0, double dx, double ftol, int maxiter) {
    return ExtrapolatingSecant(&f, x0, dx, ftol, maxiter);
}
inline double BoundedSecant(FuncWrapper1D& f, double x0, double xmin, double xmax, double dx, double ftol, int maxiter) {
    return BoundedSecant(&f, x0, xmin, xmax, dx, ftol, maxiter);
}
inline double Newton(FuncWrapper1DWithDeriv& f, double x0, double ftol, int maxiter) {
    return Newton(&f, x0, ftol, maxiter);
}
inline double BoundedNewton(FuncWrapper1DWithDeriv& f, double x0, double xmin, double xmax, double xtol, int maxiter) {
    return BoundedNewton(&f, x0, xmin, xmax, xtol, maxiter);
}
inline double Halley(FuncWrapper1DWithTwoDerivs& f, double x0, double ftol, int maxiter, double xtol_rel = 1e-12) {
    return Halley(&f, x0, ftol, maxiter, xtol_rel);
}
inline double Householder4(FuncWrapper1DWithThreeDerivs& f, double x0, double ftol, int maxiter, double xtol_rel = 1e-12) {
    return Householder4(&f, x0, ftol, maxiter, xtol_rel);
}

// Multi-Dimensional solvers
std::vector<double> NDNewtonRaphson_Jacobian(FuncWrapperND* f, const std::vector<double>& x0, double tol, int maxiter, double w = 1.0);

}; /*namespace CoolProp*/
#endif
