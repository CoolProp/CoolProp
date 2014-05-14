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
	FuncWrapper1D(){};
	virtual ~FuncWrapper1D(){};
	virtual double call(double) = 0;
	virtual double deriv(double){throw NotImplementedError("deriv function not implemented");};
};

class FuncWrapperND
{
public:
	FuncWrapperND(){};
	virtual ~FuncWrapperND(){};
	virtual std::vector<double> call(std::vector<double>) = 0;// must be provided
	virtual std::vector<std::vector<double> > Jacobian(std::vector<double>){std::vector<std::vector<double> > J; return J;}; // optional
};

// Single-Dimensional solvers
double Brent(FuncWrapper1D &f, double a, double b, double macheps, double t, int maxiter, std::string &errstr);
double Secant(FuncWrapper1D &f, double x0, double dx, double ftol, int maxiter, std::string &errstring);
double BoundedSecant(FuncWrapper1D &f, double x0, double xmin, double xmax, double dx, double ftol, int maxiter, std::string &errstring);
double Newton(FuncWrapper1D &f, double x0, double ftol, int maxiter, std::string &errstring);

// Multi-Dimensional solvers
std::vector<double> NDNewtonRaphson_Jacobian(FuncWrapperND *f, std::vector<double> x0, double tol, int maxiter, std::string *errstring);

}; /*namespace CoolProp*/
#endif
