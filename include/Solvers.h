#ifndef SOLVERS_H
#define SOLVERS_H

#include <vector>
#include <string>
#include "Exceptions.h"

namespace CoolProp
{

class FuncWrapper1D
{
public:
    int errcode;
    std::string errstring;    
    FuncWrapper1D() : errcode(0), errstring("") {};
    virtual ~FuncWrapper1D(){};
    virtual double call(double) = 0;
    /**
     * /brief A function for checking whether the input is in range;
     *
     * Meant to be implemented by derived classes; return true if input is out of range
     */
    virtual bool input_not_in_range(double x){ return false; };
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

class FuncWrapperND
{
public:
    int errcode;
    std::string errstring;    
    FuncWrapperND() : errcode(0), errstring("") {};
    virtual ~FuncWrapperND(){};
    virtual std::vector<double> call(const std::vector<double>&) = 0;// must be provided
    virtual std::vector<std::vector<double> > Jacobian(const std::vector<double>&);
};

// Single-Dimensional solvers, pointer versions
double Brent(FuncWrapper1D* f, double a, double b, double macheps, double t, int maxiter);
double Secant(FuncWrapper1D* f, double x0, double dx, double ftol, int maxiter);
double BoundedSecant(FuncWrapper1D* f, double x0, double xmin, double xmax, double dx, double ftol, int maxiter);
double Newton(FuncWrapper1DWithDeriv* f, double x0, double ftol, int maxiter);
double Halley(FuncWrapper1DWithTwoDerivs* f, double x0, double ftol, int maxiter);

// Single-Dimensional solvers
inline double Brent(FuncWrapper1D &f, double a, double b, double macheps, double t, int maxiter){
    return Brent(&f, a, b, macheps, t, maxiter);
}
inline double Secant(FuncWrapper1D &f, double x0, double dx, double ftol, int maxiter){
    return Secant(&f, x0, dx, ftol, maxiter);
}
inline double BoundedSecant(FuncWrapper1D &f, double x0, double xmin, double xmax, double dx, double ftol, int maxiter){
    return BoundedSecant(&f, x0, xmin, xmax, dx, ftol, maxiter);
}
inline double Newton(FuncWrapper1DWithDeriv &f, double x0, double ftol, int maxiter){
    return Newton(&f, x0, ftol, maxiter);
}
inline double Halley(FuncWrapper1DWithTwoDerivs &f, double x0, double ftol, int maxiter){
    return Halley(&f, x0, ftol, maxiter);
}

// Multi-Dimensional solvers
std::vector<double> NDNewtonRaphson_Jacobian(FuncWrapperND *f, std::vector<double> &x0, double tol, int maxiter);

}; /*namespace CoolProp*/
#endif
