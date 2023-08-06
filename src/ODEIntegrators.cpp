
#include "ODEIntegrators.h"
#include "Eigen/Core"
#include "CPstrings.h"
#include "Exceptions.h"
#include <algorithm>

bool ODEIntegrators::AdaptiveRK54(AbstractODEIntegrator& ode, double tstart, double tend, double hmin, double hmax, double eps_allowed,
                                  double step_relax) {
    // Get the starting array of variables of integration
    std::vector<double> xold = ode.get_initial_array();
    const long N = static_cast<long>(xold.size());

    // Start at an index of 0
    int Itheta = 0;
    double t0 = tstart;
    double h = hmin;

    // Figure out if t is increasing or decreasing in the integration and set a flag
    bool forwards_integration = ((tend - tstart) > 0);
    // If backwards integration, flip the sign of the step
    if (!forwards_integration) {
        h *= -1;
    }

    double max_error;

    std::vector<double> xnew1(N), xnew2(N), xnew3(N), xnew4(N), xnew5(N), f1(N), f2(N), f3(N), f4(N), f5(N), f6(N), error(N), xnew(N);

    // t is the independent variable here, where t takes on values in the bounded range [tmin,tmax]
    do {

        // Check for termination
        bool abort = ode.premature_termination();
        if (abort) {
            return abort;
        }

        bool stepAccepted = false, disableAdaptive = false;

        while (!stepAccepted) {

            // reset the flag
            disableAdaptive = false;

            // If the step would go beyond the end of the region of integration,
            // just take a step to the end of the region of integration
            if (forwards_integration && (t0 + h > tend)) {
                disableAdaptive = true;
                h = tend - t0;
            }
            if (!forwards_integration && (t0 + h < tend)) {
                disableAdaptive = true;
                h = tend - t0;
            }

            ode.pre_step_callback();

            // We check stepAccepted again because if the derived class
            // sets the variable stepAccepted, we should not actually do the evaluation
            if (!stepAccepted) {

                Eigen::Map<Eigen::VectorXd> xold_w(&(xold[0]), N);

                if (std::abs(h) < hmin && !disableAdaptive) {
                    // Step is too small, just use the minimum step size
                    h = (forwards_integration) ? hmin : -hmin;
                    disableAdaptive = true;
                }

                // Step 1: derivatives evaluated at old values
                ode.derivs(t0, xold, f1);

                // Call post derivative callback after the first derivative evaluation (which might cache values)
                ode.post_deriv_callback();

                Eigen::Map<Eigen::VectorXd> xnew1_w(&(xnew1[0]), N), f1_w(&(f1[0]), N);
                xnew1_w = xold_w + h * (1.0 / 5.0) * f1_w;

                ode.derivs(t0 + 1.0 / 5.0 * h, xnew1, f2);
                Eigen::Map<Eigen::VectorXd> xnew2_w(&(xnew2[0]), N), f2_w(&(f2[0]), N);
                xnew2_w = xold_w + h * (+3.0 / 40.0 * f1_w + 9.0 / 40.0 * f2_w);

                ode.derivs(t0 + 3.0 / 10.0 * h, xnew2, f3);
                Eigen::Map<Eigen::VectorXd> xnew3_w(&(xnew3[0]), N), f3_w(&(f3[0]), N);
                xnew3_w = xold_w + h * (3.0 / 10.0 * f1_w - 9.0 / 10.0 * f2_w + 6.0 / 5.0 * f3_w);

                ode.derivs(t0 + 3.0 / 5.0 * h, xnew3, f4);
                Eigen::Map<Eigen::VectorXd> xnew4_w(&(xnew4[0]), N), f4_w(&(f4[0]), N);
                xnew4_w = xold_w + h * (-11.0 / 54.0 * f1_w + 5.0 / 2.0 * f2_w - 70.0 / 27.0 * f3_w + 35.0 / 27.0 * f4_w);

                ode.derivs(t0 + h, xnew4, f5);
                Eigen::Map<Eigen::VectorXd> xnew5_w(&(xnew5[0]), N), f5_w(&(f5[0]), N);
                xnew5_w =
                  xold_w
                  + h * (1631.0 / 55296 * f1_w + 175.0 / 512.0 * f2_w + 575.0 / 13824.0 * f3_w + 44275.0 / 110592.0 * f4_w + 253.0 / 4096.0 * f5_w);

                // Updated values at the next step using 5-th order
                ode.derivs(t0 + 7.0 / 8.0 * h, xnew5, f6);
                Eigen::Map<Eigen::VectorXd> xnew_w(&(xnew[0]), N), f6_w(&(f6[0]), N);
                xnew_w = xold_w + h * (37.0 / 378.0 * f1_w + 250.0 / 621.0 * f3_w + 125.0 / 594.0 * f4_w + 512.0 / 1771.0 * f6_w);

                Eigen::Map<Eigen::VectorXd> error_w(&(error[0]), N);
                error_w =
                  h
                  * (-277.0 / 64512.0 * f1_w + 6925.0 / 370944.0 * f3_w - 6925.0 / 202752.0 * f4_w - 277.0 / 14336.0 * f5_w + 277.0 / 7084.0 * f6_w);

                max_error = error_w.norm();

                // If the error is too large, make the step size smaller and try
                // the step again
                if (disableAdaptive) {
                    // Accept the step regardless of whether the error
                    // is too large or not
                    stepAccepted = true;
                } else {
                    if (max_error > eps_allowed) {
                        // Take a smaller step next time, try again on this step
                        // But only if adaptive mode is on
                        // If eps_allowed == max_error (approximately), force the step to change to avoid infinite loop
                        h *= std::min(step_relax * pow(eps_allowed / max_error, 0.3), 0.999);
                        stepAccepted = false;
                    } else {
                        stepAccepted = true;
                    }
                }

            } else {
                std::cout << format("accepted");
            }
        }

        // Step has been accepted, update variables
        t0 += h;
        Itheta += 1;
        xold = xnew;

        ode.post_step_callback(t0, h, xnew);

        // The error is already below the threshold
        if (max_error < eps_allowed && disableAdaptive == false && max_error > 0) {
            // Take a bigger step next time, since eps_allowed>max_error, but don't
            // let the steps get much larger too quickly
            h *= step_relax * pow(eps_allowed / max_error, 0.2);
        }

        // Constrain the step to not be too large
        if (forwards_integration) {
            h = std::min(h, hmax);
        } else {
            h = -std::min(std::abs(h), hmax);
        }

        // Overshot the end, oops...  That's an error
        if (forwards_integration && (t0 - tend > +1e-3)) {
            throw CoolProp::ValueError(format("t0 - tend [%g] > 1e-3", t0 - tend));
        }
        if (!forwards_integration && (t0 - tend < -1e-3)) {
            throw CoolProp::ValueError(format("t0 - tend [%g] < -1e-3", t0 - tend));
        }
    } while (((forwards_integration) && t0 < tend - 1e-10) || ((!forwards_integration) && t0 > tend + 1e-10));

    // No termination was requested
    return false;
}

#if defined(ENABLE_CATCH)
#    include <catch2/catch_all.hpp>

TEST_CASE("Integrate y'=y", "[ODEIntegrator]") {
    class SimpleODEIntegrator : public ODEIntegrators::AbstractODEIntegrator
    {
       public:
        std::vector<double> t, h, y;

        virtual std::vector<double> get_initial_array() const {
            return std::vector<double>(1, 1);
        }

        virtual void pre_step_callback(){};

        virtual void post_deriv_callback(){};

        virtual void post_step_callback(double t, double h, std::vector<double>& y) {
            this->t.push_back(t);
            this->h.push_back(h);
            this->y.push_back(y[0]);
        };

        virtual bool premature_termination() {
            return false;
        };

        virtual void derivs(double t, std::vector<double>& y, std::vector<double>& yprime) {
            yprime[0] = y[0];
        };
    };

    SimpleODEIntegrator simple;
    ODEIntegrators::AdaptiveRK54(simple, 0, 4, 1e-4, 0.5, 1e-7, 0.9);
    double yfinal_integration = simple.y[simple.y.size() - 1];
    double tfinal_integration = simple.t[simple.t.size() - 1];

    double yfinal_analytic = exp(4.0);
    double error = yfinal_integration / yfinal_analytic - 1;

    CAPTURE(yfinal_analytic);
    CAPTURE(yfinal_integration);
    CAPTURE(tfinal_integration);
    CHECK(std::abs(error) < 1e-6);
    CHECK(std::abs(tfinal_integration - 4) < 1e-10);
    int rr = 0;
}
#endif
