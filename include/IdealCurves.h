#include "AbstractState.h"
#include "crossplatform_shared_ptr.h"
#include "Solvers.h"
#include "CoolPropTools.h"
#include <string>

namespace CoolProp {

class CurveTracer : public FuncWrapper1D
{
   public:
    AbstractState* AS;
    double p0, T0, lnT, lnp, rho_guess;
    std::vector<double> T, p;
    enum OBJECTIVE_TYPE
    {
        OBJECTIVE_INVALID = 0,
        OBJECTIVE_CIRCLE,
        OBJECTIVE_T
    };
    OBJECTIVE_TYPE obj;
    CurveTracer(AbstractState* AS, double p0, double T0) : AS(AS), p0(p0), T0(T0), lnT(_HUGE), lnp(_HUGE), rho_guess(_HUGE), obj(OBJECTIVE_INVALID) {
        this->p.push_back(p0);
    };
    void init() {
        // Solve for Temperature for first point
        this->obj = OBJECTIVE_T;
        this->rho_guess = -1;
        this->T.push_back(Secant(this, T0, 0.001 * T0, 1e-10, 100));
    }

    virtual double objective(void) = 0;

    virtual double starting_direction() {
        return M_PI / 2.0;
    }

    double call(double t) {
        if (this->obj == OBJECTIVE_CIRCLE) {
            double T2, P2;
            this->TPcoords(t, lnT, lnp, T2, P2);
            this->AS->update(PT_INPUTS, P2, T2);
        } else {
            if (this->rho_guess < 0)
                this->AS->update(PT_INPUTS, this->p[this->p.size() - 1], t);
            else {
                GuessesStructure guesses;
                guesses.rhomolar = this->rho_guess;
                this->AS->update_with_guesses(PT_INPUTS, this->p[this->p.size() - 1], t, guesses);
            }
        }
        double r = this->objective();
        return r;
    }

    void TPcoords(double t, double lnT, double lnp, double& T, double& p) {
        double rlnT = 0.1, rlnp = 0.1;
        T = exp(lnT + rlnT * cos(t));
        p = exp(lnp + rlnp * sin(t));
    }

    void trace(std::vector<double>& T, std::vector<double>& p) {
        double t = this->starting_direction();
        for (int i = 0; i < 1000; ++i) {
            try {
                this->lnT = log(this->T[this->T.size() - 1]);
                this->lnp = log(this->p[this->p.size() - 1]);
                this->obj = OBJECTIVE_CIRCLE;
                t = Brent(this, t - M_PI / 2.0, t + M_PI / 2.0, DBL_EPSILON, 1e-10, 100);
                double T2, P2;
                this->TPcoords(t, this->lnT, this->lnp, T2, P2);
                this->T.push_back(T2);
                this->p.push_back(P2);
                if (this->T[this->T.size() - 1] < this->AS->keyed_output(iT_triple)
                    || this->p[this->p.size() - 1] > 1000 * this->AS->keyed_output(iP_critical)) {
                    break;
                }
            } catch (std::exception&) {
                break;
            }
        }
        T = this->T;
        p = this->p;
    }
};

class IdealCurveTracer : public CurveTracer
{
   public:
    IdealCurveTracer(AbstractState* AS, double p0, double T0) : CurveTracer(AS, p0, T0) {
        init();
    };
    /// Z = 1
    double objective(void) {
        return this->AS->keyed_output(iZ) - 1;
    };
};

class BoyleCurveTracer : public CurveTracer
{
   public:
    BoyleCurveTracer(AbstractState* AS, double p0, double T0) : CurveTracer(AS, p0, T0) {
        init();
    };
    /// dZ/dv|T = 0
    double objective(void) {
        double r =
          (this->AS->p() - this->AS->rhomolar() * this->AS->first_partial_deriv(iP, iDmolar, iT)) / (this->AS->gas_constant() * this->AS->T());
        return r;
    };
};
class JouleInversionCurveTracer : public CurveTracer
{
   public:
    JouleInversionCurveTracer(AbstractState* AS, double p0, double T0) : CurveTracer(AS, p0, T0) {
        init();
    };
    /// dZ/dT|v = 0
    double objective(void) {
        double r = (this->AS->gas_constant() * this->AS->T() * 1 / this->AS->rhomolar() * this->AS->first_partial_deriv(iP, iT, iDmolar)
                    - this->AS->p() * this->AS->gas_constant() / this->AS->rhomolar())
                   / POW2(this->AS->gas_constant() * this->AS->T());
        return r;
    };
};
class JouleThomsonCurveTracer : public CurveTracer
{
   public:
    JouleThomsonCurveTracer(AbstractState* AS, double p0, double T0) : CurveTracer(AS, p0, T0) {
        init();
    };
    /// dZ/dT|p = 0
    double objective(void) {
        double dvdT__constp = -this->AS->first_partial_deriv(iDmolar, iT, iP) / POW2(this->AS->rhomolar());
        double r = this->AS->p() / (this->AS->gas_constant() * POW2(this->AS->T())) * (this->AS->T() * dvdT__constp - 1 / this->AS->rhomolar());
        return r;
    };
};

} /* namespace CoolProp */