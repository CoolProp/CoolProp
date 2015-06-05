#ifndef EXCESSHE_FUNCTIONS_H
#define EXCESSHE_FUNCTIONS_H

#include <memory>
#include <vector>
#include "CoolPropFluid.h"
#include "crossplatform_shared_ptr.h"

namespace CoolProp{

typedef std::vector<std::vector<CoolPropDbl> > STLMatrix;

/** \brief The abstract base class for departure functions used in the excess part of the Helmholtz energy
 * 
 * The only code included in the ABC is the structure for the derivatives of the Helmholtz energy with 
 * the reduced density and reciprocal reduced temperature
 */
class DepartureFunction
{
public:
    DepartureFunction(){};
    virtual ~DepartureFunction(){};
    ResidualHelmholtzGeneralizedExponential phi;
    HelmholtzDerivatives derivs;

    void update(double tau, double delta){
        derivs.reset(0.0);
        phi.all(tau, delta, derivs);
    };

    double alphar(){ return derivs.alphar;};
    double dalphar_dDelta(){ return derivs.dalphar_ddelta;};
    double dalphar_dTau(){ return derivs.dalphar_dtau;};
    
    double d2alphar_dDelta2(){return derivs.d2alphar_ddelta2;};
    double d2alphar_dDelta_dTau(){return derivs.d2alphar_ddelta_dtau;};
    double d2alphar_dTau2(){return derivs.d2alphar_dtau2;};

	double d3alphar_dTau3(){ return derivs.d3alphar_dtau3; };
	double d3alphar_dDelta_dTau2(){ return derivs.d3alphar_ddelta_dtau2; }; 
	double d3alphar_dDelta2_dTau(){ return derivs.d3alphar_ddelta2_dtau; }; 
	double d3alphar_dDelta3(){ return derivs.d3alphar_ddelta3; };
};

/** \brief The departure function used by the GERG-2008 formulation
 * 
 * This departure function has a form like
 * \f[
 * \alphar^r_{ij} = \sum_k n_{ij,k}\delta^{d_{ij,k}}\tau^{t_{ij,k}} + \sum_k n_{ij,k}\delta^{d_{ij,k}}\tau^{t_{ij,k}}\exp[-\eta_{ij,k}(\delta-\varepsilon_{ij,k})^2-\beta_{ij,k}(\delta-\gamma_{ij,k})]
 * \f]
 * It is symmetric so \f$\alphar^r_{ij} = \alphar^r_{ji}\f$
 */
class GERG2008DepartureFunction : public DepartureFunction
{    
public:
    GERG2008DepartureFunction(){};
    GERG2008DepartureFunction(const std::vector<double> &n,const std::vector<double> &d,const std::vector<double> &t,
                              const std::vector<double> &eta,const std::vector<double> &epsilon,const std::vector<double> &beta,
                              const std::vector<double> &gamma, std::size_t Npower)
    {
        /// Break up into power and gaussian terms
        {
            std::vector<CoolPropDbl> _n(n.begin(), n.begin()+Npower);
            std::vector<CoolPropDbl> _d(d.begin(), d.begin()+Npower);
            std::vector<CoolPropDbl> _t(t.begin(), t.begin()+Npower);
            std::vector<CoolPropDbl> _l(Npower, 0.0);
            phi.add_Power(_n, _d, _t, _l);
        }
        if (n.size() == Npower)
        {
        }
        else
        {
            std::vector<CoolPropDbl> _n(n.begin()+Npower,                   n.end());
            std::vector<CoolPropDbl> _d(d.begin()+Npower,                   d.end());
            std::vector<CoolPropDbl> _t(t.begin()+Npower,                   t.end());
            std::vector<CoolPropDbl> _eta(eta.begin()+Npower,             eta.end());
            std::vector<CoolPropDbl> _epsilon(epsilon.begin()+Npower, epsilon.end());
            std::vector<CoolPropDbl> _beta(beta.begin()+Npower,          beta.end());
            std::vector<CoolPropDbl> _gamma(gamma.begin()+Npower,       gamma.end());
            phi.add_GERG2008Gaussian(_n, _d, _t, _eta, _epsilon, _beta, _gamma);
        }
    };
    ~GERG2008DepartureFunction(){};
};

/** \brief A polynomial/exponential departure function
 * 
 * This departure function has a form like
 * \f[
 * \alpha^r_{ij} = \sum_k n_{ij,k}\delta^{d_{ij,k}}\tau^{t_{ij,k}}\exp(-\delta^{l_{ij,k}})
 * \f]
 * It is symmetric so \f$\alphar^r_{ij} = \alphar^r_{ji}\f$
 */
class ExponentialDepartureFunction : public DepartureFunction
{
public:
    ExponentialDepartureFunction(){};
    ExponentialDepartureFunction(const std::vector<double> &n, const std::vector<double> &d,
                                 const std::vector<double> &t, const std::vector<double> &l)
                                 {
                                     std::vector<CoolPropDbl> _n(n.begin(), n.begin()+n.size());
                                     std::vector<CoolPropDbl> _d(d.begin(), d.begin()+d.size());
                                     std::vector<CoolPropDbl> _t(t.begin(), t.begin()+t.size());
                                     std::vector<CoolPropDbl> _l(l.begin(), l.begin()+l.size());
                                     phi.add_Power(_n, _d, _t, _l);
                                 };
    ~ExponentialDepartureFunction(){};
};

typedef shared_ptr<DepartureFunction> DepartureFunctionPointer;

class ExcessTerm
{
public:
    std::size_t N;
    std::vector<std::vector<DepartureFunctionPointer> > DepartureFunctionMatrix;
    STLMatrix F;

    ExcessTerm():N(0){};

    /// Resize the parts of this term
    void resize(std::size_t N){
        this->N = N;
        F.resize(N, std::vector<CoolPropDbl>(N, 0));
        DepartureFunctionMatrix.resize(N);
        for (std::size_t i = 0; i < N; ++i){
            DepartureFunctionMatrix[i].resize(N);
        }
    };
    /// Update the internal cached derivatives in each departure function
    void update(double tau, double delta){
        for (std::size_t i = 0; i < N; i++){
            for (std::size_t j = i + 1; j < N; j++){
                DepartureFunctionMatrix[i][j]->update(tau, delta);
            }
            for (std::size_t j = 0; j < i; j++){
                DepartureFunctionMatrix[i][j]->update(tau, delta);
            }
        }
    }

    double alphar(const std::vector<CoolPropDbl> &x)
    {
        double summer = 0;
        for (std::size_t i = 0; i < N-1; i++)
        {
            for (std::size_t j = i + 1; j < N; j++)
            {
                summer += x[i]*x[j]*F[i][j]*DepartureFunctionMatrix[i][j]->alphar();
            }
        }
        return summer;
    }
    double dalphar_dDelta(const std::vector<CoolPropDbl> &x)
    {
        double summer = 0;
        for (std::size_t i = 0; i < N-1; i++)
        {
            for (std::size_t j = i + 1; j < N; j++)
            {
                summer += x[i]*x[j]*F[i][j]*DepartureFunctionMatrix[i][j]->dalphar_dDelta();
            }
        }
        return summer;
    }
    double d2alphar_dDelta2(const std::vector<CoolPropDbl> &x)
    {
        double summer = 0;
        for (std::size_t i = 0; i < N-1; i++)
        {
            for (std::size_t j = i + 1; j < N; j++)
            {
                summer += x[i]*x[j]*F[i][j]*DepartureFunctionMatrix[i][j]->d2alphar_dDelta2();
            }
        }
        return summer;
    };
    double d2alphar_dDelta_dTau(const std::vector<CoolPropDbl> &x)
    {
        double summer = 0;
        for (std::size_t i = 0; i < N-1; i++)
        {
            for (std::size_t j = i + 1; j < N; j++)
            {
                summer += x[i]*x[j]*F[i][j]*DepartureFunctionMatrix[i][j]->d2alphar_dDelta_dTau();
            }
        }
        return summer;
    }
    double dalphar_dTau(const std::vector<CoolPropDbl> &x)
    {
        double summer = 0;
        for (std::size_t i = 0; i < N-1; i++)
        {
            for (std::size_t j = i + 1; j < N; j++)
            {
                summer += x[i]*x[j]*F[i][j]*DepartureFunctionMatrix[i][j]->dalphar_dTau();
            }
        }
        return summer;
    };
    double d2alphar_dTau2(const std::vector<CoolPropDbl> &x)
    {
        double summer = 0;
        for (std::size_t i = 0; i < N-1; i++)
        {
            for (std::size_t j = i + 1; j < N; j++)
            {
                summer += x[i]*x[j]*F[i][j]*DepartureFunctionMatrix[i][j]->d2alphar_dTau2();
            }
        }
        return summer;
    };
	double d3alphar_dTau3(const std::vector<CoolPropDbl> &x)
	{
		double summer = 0;
		for (std::size_t i = 0; i < N - 1; i++)
		{
			for (std::size_t j = i + 1; j < N; j++)
			{
				summer += x[i] * x[j] * F[i][j] * DepartureFunctionMatrix[i][j]->d3alphar_dTau3();
			}
		}
		return summer;
	};
	double d3alphar_dDelta_dTau2(const std::vector<CoolPropDbl> &x)
	{
		double summer = 0;
		for (std::size_t i = 0; i < N - 1; i++)
		{
			for (std::size_t j = i + 1; j < N; j++)
			{
				summer += x[i] * x[j] * F[i][j] * DepartureFunctionMatrix[i][j]->d3alphar_dDelta_dTau2();
			}
		}
		return summer;
	};
	double d3alphar_dDelta2_dTau(const std::vector<CoolPropDbl> &x)
	{
		double summer = 0;
		for (std::size_t i = 0; i < N - 1; i++)
		{
			for (std::size_t j = i + 1; j < N; j++)
			{
				summer += x[i] * x[j] * F[i][j] * DepartureFunctionMatrix[i][j]->d3alphar_dDelta2_dTau();
			}
		}
		return summer;
	};
	double d3alphar_dDelta3(const std::vector<CoolPropDbl> &x)
	{
		double summer = 0;
		for (std::size_t i = 0; i < N - 1; i++)
		{
			for (std::size_t j = i + 1; j < N; j++)
			{
				summer += x[i] * x[j] * F[i][j] * DepartureFunctionMatrix[i][j]->d3alphar_dDelta3();
			}
		}
		return summer;
	};

    double dalphar_dxi(const std::vector<CoolPropDbl> &x, std::size_t i)
    {
        double summer = 0;
        for (std::size_t k = 0; k < N; k++)
        {
            if (i != k)
            {
                summer += x[k]*F[i][k]*DepartureFunctionMatrix[i][k]->alphar();
            }
        }
        return summer;
    };
    double d2alphardxidxj(const std::vector<CoolPropDbl> &x, std::size_t i, std::size_t j)
    {
        if (i != j)
        {
            return F[i][j]*DepartureFunctionMatrix[i][j]->alphar();
        }
        else
        {
            return 0;
        }
    };
    double d3alphar_dxi_dxj_dDelta(const std::vector<CoolPropDbl> &x, std::size_t i, std::size_t j)
    {
        if (i != j)
        {
            return F[i][j]*DepartureFunctionMatrix[i][j]->dalphar_dDelta();
        }
        else
        {
            return 0;
        }
    };
    double d3alphar_dxi_dxj_dTau(const std::vector<CoolPropDbl> &x, std::size_t i, std::size_t j)
    {
        if (i != j)
        {
            return F[i][j]*DepartureFunctionMatrix[i][j]->dalphar_dTau();
        }
        else
        {
            return 0;
        }
    };
    double d3alphardxidxjdxk(const std::vector<CoolPropDbl> &x, std::size_t i, std::size_t j, std::size_t k)
    {
        return 0;
    };
    double d2alphar_dxi_dTau(const std::vector<CoolPropDbl> &x, std::size_t i)
    {
        double summer = 0;
        for (std::size_t k = 0; k < N; k++)
        {
            if (i != k)
            {
                summer += x[k]*F[i][k]*DepartureFunctionMatrix[i][k]->dalphar_dTau();
            }
        }
        return summer;
    };
    double d2alphar_dxi_dDelta(const std::vector<CoolPropDbl> &x, std::size_t i)
    {
        double summer = 0;
        for (std::size_t k = 0; k < N; k++)
        {
            if (i != k)
            {
                summer += x[k]*F[i][k]*DepartureFunctionMatrix[i][k]->dalphar_dDelta();
            }
        }
        return summer;
    };
	double d3alphar_dxi_dDelta2(const std::vector<CoolPropDbl> &x, std::size_t i)
	{
		double summer = 0;
		for (std::size_t k = 0; k < N; k++)
		{
			if (i != k)
			{
				summer += x[k] * F[i][k] * DepartureFunctionMatrix[i][k]->d2alphar_dDelta2();
			}
		}
		return summer;
	};
	double d3alphar_dxi_dTau2(const std::vector<CoolPropDbl> &x, std::size_t i)
	{
		double summer = 0;
		for (std::size_t k = 0; k < N; k++)
		{
			if (i != k)
			{
				summer += x[k] * F[i][k] * DepartureFunctionMatrix[i][k]->d2alphar_dTau2();
			}
		}
		return summer;
	};
	double d3alphar_dxi_dDelta_dTau(const std::vector<CoolPropDbl> &x, std::size_t i)
	{
		double summer = 0;
		for (std::size_t k = 0; k < N; k++)
		{
			if (i != k)
			{
				summer += x[k] * F[i][k] * DepartureFunctionMatrix[i][k]->d2alphar_dDelta_dTau();
			}
		}
		return summer;
	};
};

} /* namespace CoolProp */
#endif
