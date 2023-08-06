#ifndef ODEINTEGRATORS_H
#define ODEINTEGRATORS_H

#include <vector>

namespace ODEIntegrators {

/// The abstract class defining the interface for the integrator routines
class AbstractODEIntegrator
{
   public:
    virtual std::vector<double> get_initial_array() const = 0;

    virtual void pre_step_callback() = 0;

    virtual void post_deriv_callback() = 0;

    virtual void post_step_callback(double t, double h, std::vector<double>& x) = 0;

    virtual bool premature_termination() = 0;

    virtual void derivs(double t, std::vector<double>& x, std::vector<double>& f) = 0;
};

/**
     @brief Use the adaptive Runge-Kutta integrator to integrate a system of differential equations

     @param tmin Starting value of the independent variable.  ``t`` is in the closed range [``tmin``, ``tmax``]
     @param tmax Ending value for the independent variable.  ``t`` is in the closed range [``tmin``, ``tmax``]
     @param hmin Minimum step size, something like 1e-5 usually is good.  Don't make this too big or you may not be able to get a stable solution
     @param hmax Maximum step size
     @param eps_allowed Maximum absolute error of any CV per step allowed.  Don't make this parameter too big or you may not be able to get a stable solution.  Also don't make it too small because then you are going to run into truncation error.
     @param step_relax The relaxation factor that is used in the step resizing algorithm.  Should be less than 1.0; you can play with this parameter to improve the adaptive resizing, but should not be necessary.

     */
bool AdaptiveRK54(AbstractODEIntegrator& ode, double tmin, double tmax, double hmin, double hmax, double eps_allowed, double step_relax);

}  // namespace ODEIntegrators

#endif