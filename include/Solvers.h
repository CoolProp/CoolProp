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
    virtual ~FuncWrapper1D(){};
    virtual double call(double) = 0;
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
    FuncWrapperND(){};
    virtual ~FuncWrapperND(){};
    virtual std::vector<double> call(const std::vector<double>&) = 0;// must be provided
    virtual std::vector<std::vector<double> > Jacobian(const std::vector<double>&);
};

// Single-Dimensional solvers, pointer versions
double Brent(FuncWrapper1D* f, double a, double b, double macheps, double t, int maxiter, std::string &errstr);
double Secant(FuncWrapper1D* f, double x0, double dx, double ftol, int maxiter, std::string &errstring);
double BoundedSecant(FuncWrapper1D* f, double x0, double xmin, double xmax, double dx, double ftol, int maxiter, std::string &errstring);
double Newton(FuncWrapper1DWithDeriv* f, double x0, double ftol, int maxiter, std::string &errstring);
double Halley(FuncWrapper1DWithTwoDerivs* f, double x0, double ftol, int maxiter, std::string &errstring);

// Single-Dimensional solvers
inline double Brent(FuncWrapper1D &f, double a, double b, double macheps, double t, int maxiter, std::string &errstr){
    return Brent(&f, a, b, macheps, t, maxiter, errstr);
}
inline double Secant(FuncWrapper1D &f, double x0, double dx, double ftol, int maxiter, std::string &errstring){
    return Secant(&f, x0, dx, ftol, maxiter, errstring);
}
inline double BoundedSecant(FuncWrapper1D &f, double x0, double xmin, double xmax, double dx, double ftol, int maxiter, std::string &errstring){
    return BoundedSecant(&f, x0, xmin, xmax, dx, ftol, maxiter, errstring);
}
inline double Newton(FuncWrapper1DWithDeriv &f, double x0, double ftol, int maxiter, std::string &errstring){
    return Newton(&f, x0, ftol, maxiter, errstring);
}
inline double Halley(FuncWrapper1DWithTwoDerivs &f, double x0, double ftol, int maxiter, std::string &errstring){
    return Halley(&f, x0, ftol, maxiter, errstring);
}

// Multi-Dimensional solvers
std::vector<double> NDNewtonRaphson_Jacobian(FuncWrapperND *f, const std::vector<double> &x0, double tol, int maxiter, std::string *errstring);

}; /*namespace CoolProp*/
#endif
