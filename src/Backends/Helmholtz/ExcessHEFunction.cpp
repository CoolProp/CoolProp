#include <memory>

#include "ReducingFunctions.h"
#include "ExcessHEFunction.h"

namespace CoolProp{

double ExcessTerm::alphar(double tau, double delta, const std::vector<long double> &x)
{
	double summer = 0;
	for (unsigned int i = 0; i < N-1; i++)
	{
		for (unsigned int j = i + 1; j < N; j++)
		{
			summer += x[i]*x[j]*F[i][j]*DepartureFunctionMatrix[i][j]->alphar(tau,delta);
		}
	}
	return summer;
}
double ExcessTerm::dalphar_dTau(double tau, double delta, const std::vector<long double> &x)
{
	double summer = 0;
	for (unsigned int i = 0; i < N-1; i++)
	{
		for (unsigned int j = i + 1; j < N; j++)
		{
			summer += x[i]*x[j]*F[i][j]*DepartureFunctionMatrix[i][j]->dalphar_dTau(tau,delta);
		}
	}
	return summer;
}
double ExcessTerm::dalphar_dDelta(double tau, double delta, const std::vector<long double> &x)
{
	double summer = 0;
	for (unsigned int i = 0; i < N-1; i++)
	{
		for (unsigned int j = i + 1; j < N; j++)
		{
			summer += x[i]*x[j]*F[i][j]*DepartureFunctionMatrix[i][j]->dalphar_dDelta(tau,delta);
		}
	}
	return summer;
}
double ExcessTerm::d2alphar_dDelta2(double tau, double delta, const std::vector<long double> &x)
{
	double summer = 0;
	for (unsigned int i = 0; i < N-1; i++)
	{
		for (unsigned int j = i + 1; j < N; j++)
		{
			summer += x[i]*x[j]*F[i][j]*DepartureFunctionMatrix[i][j]->d2alphar_dDelta2(tau,delta);
		}
	}
	return summer;
}
double ExcessTerm::d2alphar_dTau2(double tau, double delta, const std::vector<long double> &x)
{
	double summer = 0;
	for (unsigned int i = 0; i < N-1; i++)
	{
		for (unsigned int j = i + 1; j < N; j++)
		{
			summer += x[i]*x[j]*F[i][j]*DepartureFunctionMatrix[i][j]->d2alphar_dTau2(tau,delta);
		}
	}
	return summer;
}
double ExcessTerm::d2alphar_dDelta_dTau(double tau, double delta, const std::vector<long double> &x)
{
	double summer = 0;
	for (unsigned int i = 0; i < N-1; i++)
	{
		for (unsigned int j = i + 1; j < N; j++)
		{
			summer += x[i]*x[j]*F[i][j]*DepartureFunctionMatrix[i][j]->d2alphar_dDelta_dTau(tau,delta);
		}
	}
	return summer;
}
double ExcessTerm::dalphar_dxi(double tau, double delta, const std::vector<long double> &x, unsigned int i)
{
	double summer = 0;
	for (unsigned int k = 0; k < N; k++)
	{
		if (i != k)
		{
			summer += x[k]*F[i][k]*DepartureFunctionMatrix[i][k]->alphar(tau,delta);
		}
	}
	return summer;
}
double ExcessTerm::d2alphardxidxj(double tau, double delta, const std::vector<long double> &x, unsigned int i, unsigned int j)
{
	if (i != j)
	{
		return F[i][j]*DepartureFunctionMatrix[i][j]->alphar(tau,delta);
	}
	else
	{
		return 0;
	}
}
double ExcessTerm::d2alphar_dxi_dTau(double tau, double delta, const std::vector<long double> &x, unsigned int i)
{
	double summer = 0;
	for (unsigned int k = 0; k < N; k++)
	{
		if (i != k)
		{
			summer += x[k]*F[i][k]*DepartureFunctionMatrix[i][k]->dalphar_dTau(tau,delta);
		}
	}
	return summer;
}
double ExcessTerm::d2alphar_dxi_dDelta(double tau, double delta, const std::vector<long double> &x, unsigned int i)
{
	double summer = 0;
	for (unsigned int k = 0; k < N; k++)
	{
		if (i != k)
		{
			summer += x[k]*F[i][k]*DepartureFunctionMatrix[i][k]->dalphar_dDelta(tau,delta);
		}
	}
	return summer;
}

GERG2008DepartureFunction::GERG2008DepartureFunction(const std::vector<double> &n,const std::vector<double> &d,const std::vector<double> &t,
                                                     const std::vector<double> &eta,const std::vector<double> &epsilon,const std::vector<double> &beta,
                                                     const std::vector<double> &gamma, unsigned int Npower)
{
    /// Break up into power and gaussian terms
    {
        std::vector<long double> _n(n.begin(), n.begin()+Npower);
        std::vector<long double> _d(d.begin(), d.begin()+Npower);
        std::vector<long double> _t(t.begin(), t.begin()+Npower);
        std::vector<long double> _l(Npower, 0.0);
        phi.add_Power(_n, _d, _t, _l);
    }
    if (n.size() == Npower)
    {
        using_gaussian = false;
    }
    else
    {
        using_gaussian = true;
        std::vector<long double> _n(n.begin()+Npower,                   n.end());
        std::vector<long double> _d(d.begin()+Npower,                   d.end());
        std::vector<long double> _t(t.begin()+Npower,                   t.end());
        std::vector<long double> _eta(eta.begin()+Npower,             eta.end());
        std::vector<long double> _epsilon(epsilon.begin()+Npower, epsilon.end());
        std::vector<long double> _beta(beta.begin()+Npower,          beta.end());
        std::vector<long double> _gamma(gamma.begin()+Npower,       gamma.end());
        phi.add_GERG2008Gaussian(_n, _d, _t, _eta, _epsilon, _beta, _gamma);
    }
}

} /* namespace CoolProp */
