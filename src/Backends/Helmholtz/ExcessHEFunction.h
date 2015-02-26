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

    /// The excess Helmholtz energy of the binary pair
    /// Pure-virtual function (must be implemented in derived class
    virtual double alphar(double tau, double delta) = 0;
    virtual double dalphar_dDelta(double tau, double delta) = 0;
    virtual double d2alphar_dDelta2(double tau, double delta) = 0;
    virtual double d2alphar_dDelta_dTau(double tau, double delta) = 0;
    virtual double dalphar_dTau(double tau, double delta) = 0;
    virtual double d2alphar_dTau2(double tau, double delta) = 0;
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
protected:
    bool using_gaussian;
    ResidualHelmholtzGeneralizedExponential phi;
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
            using_gaussian = false;
        }
        else
        {
            using_gaussian = true;
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

    double alphar(double tau, double delta){return phi.base(tau, delta);};
    double dalphar_dDelta(double tau, double delta){return phi.dDelta(tau, delta);};
    double d2alphar_dDelta_dTau(double tau, double delta){return phi.dDelta_dTau(tau, delta);};
    double dalphar_dTau(double tau, double delta){return phi.dTau(tau, delta);};
    double d2alphar_dDelta2(double tau, double delta){return phi.dDelta2(tau, delta);};
    double d2alphar_dTau2(double tau, double delta){return phi.dTau2(tau, delta);};
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
protected:
    ResidualHelmholtzGeneralizedExponential phi;
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

    double alphar(double tau, double delta){return phi.base(tau, delta);};
    double dalphar_dDelta(double tau, double delta){return phi.dDelta(tau, delta);};
    double d2alphar_dDelta_dTau(double tau, double delta){return phi.dDelta_dTau(tau, delta);};
    double dalphar_dTau(double tau, double delta){return phi.dTau(tau, delta);};
    double d2alphar_dDelta2(double tau, double delta){return phi.dDelta2(tau, delta);};
    double d2alphar_dTau2(double tau, double delta){return phi.dTau2(tau, delta);};
};

typedef shared_ptr<DepartureFunction> DepartureFunctionPointer;

class ExcessTerm
{
public:
    std::size_t N;
    std::vector<std::vector<DepartureFunctionPointer> > DepartureFunctionMatrix;
    STLMatrix F;

    ExcessTerm(){};
    ~ExcessTerm(){};

    /// Resize the parts of this term
    void resize(std::size_t N){
        this->N = N;
        F.resize(N, std::vector<CoolPropDbl>(N, 0));
        DepartureFunctionMatrix.resize(N);
        for (std::size_t i = 0; i < N; ++i){
            DepartureFunctionMatrix[i].resize(N);
        }
    };

    double alphar(double tau, double delta, const std::vector<CoolPropDbl> &x)
    {
        double summer = 0;
        for (std::size_t i = 0; i < N-1; i++)
        {
            for (std::size_t j = i + 1; j < N; j++)
            {
                summer += x[i]*x[j]*F[i][j]*DepartureFunctionMatrix[i][j]->alphar(tau,delta);
            }
        }
        return summer;
    }
    double dalphar_dDelta(double tau, double delta, const std::vector<CoolPropDbl> &x)
    {
        double summer = 0;
        for (std::size_t i = 0; i < N-1; i++)
        {
            for (std::size_t j = i + 1; j < N; j++)
            {
                summer += x[i]*x[j]*F[i][j]*DepartureFunctionMatrix[i][j]->dalphar_dDelta(tau,delta);
            }
        }
        return summer;
    }
    double d2alphar_dDelta2(double tau, double delta, const std::vector<CoolPropDbl> &x)
    {
        double summer = 0;
        for (std::size_t i = 0; i < N-1; i++)
        {
            for (std::size_t j = i + 1; j < N; j++)
            {
                summer += x[i]*x[j]*F[i][j]*DepartureFunctionMatrix[i][j]->d2alphar_dDelta2(tau,delta);
            }
        }
        return summer;
    };
    double d2alphar_dDelta_dTau(double tau, double delta, const std::vector<CoolPropDbl> &x)
    {
        double summer = 0;
        for (std::size_t i = 0; i < N-1; i++)
        {
            for (std::size_t j = i + 1; j < N; j++)
            {
                summer += x[i]*x[j]*F[i][j]*DepartureFunctionMatrix[i][j]->d2alphar_dDelta_dTau(tau,delta);
            }
        }
        return summer;
    }
    double dalphar_dTau(double tau, double delta, const std::vector<CoolPropDbl> &x)
    {
        double summer = 0;
        for (std::size_t i = 0; i < N-1; i++)
        {
            for (std::size_t j = i + 1; j < N; j++)
            {
                summer += x[i]*x[j]*F[i][j]*DepartureFunctionMatrix[i][j]->dalphar_dTau(tau,delta);
            }
        }
        return summer;
    };
    double d2alphar_dTau2(double tau, double delta, const std::vector<CoolPropDbl> &x)
    {
        double summer = 0;
        for (std::size_t i = 0; i < N-1; i++)
        {
            for (std::size_t j = i + 1; j < N; j++)
            {
                summer += x[i]*x[j]*F[i][j]*DepartureFunctionMatrix[i][j]->d2alphar_dTau2(tau,delta);
            }
        }
        return summer;
    };
    double dalphar_dxi(double tau, double delta, const std::vector<CoolPropDbl> &x, std::size_t i)
    {
        double summer = 0;
        for (std::size_t k = 0; k < N; k++)
        {
            if (i != k)
            {
                summer += x[k]*F[i][k]*DepartureFunctionMatrix[i][k]->alphar(tau,delta);
            }
        }
        return summer;
    };
    double d2alphardxidxj(double tau, double delta, const std::vector<CoolPropDbl> &x, std::size_t i, std::size_t j)
    {
        if (i != j)
        {
            return F[i][j]*DepartureFunctionMatrix[i][j]->alphar(tau,delta);
        }
        else
        {
            return 0;
        }
    };
    double d2alphar_dxi_dTau(double tau, double delta, const std::vector<CoolPropDbl> &x, std::size_t i)
    {
        double summer = 0;
        for (std::size_t k = 0; k < N; k++)
        {
            if (i != k)
            {
                summer += x[k]*F[i][k]*DepartureFunctionMatrix[i][k]->dalphar_dTau(tau,delta);
            }
        }
        return summer;
    };
    double d2alphar_dxi_dDelta(double tau, double delta, const std::vector<CoolPropDbl> &x, std::size_t i)
    {
        double summer = 0;
        for (std::size_t k = 0; k < N; k++)
        {
            if (i != k)
            {
                summer += x[k]*F[i][k]*DepartureFunctionMatrix[i][k]->dalphar_dDelta(tau,delta);
            }
        }
        return summer;
    };
};

} /* namespace CoolProp */
#endif
