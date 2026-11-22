#include <cmath>
#include <vector>
#include "CoolProp/numerics/Solvers.h"
#include "CoolProp/numerics/MatrixMath.h"
#include <iostream>
#include "CoolProp/detail/tools.h"
#include <Eigen/Dense>

namespace CoolProp {

/** \brief Calculate the Jacobian using numerical differentiation by column
 */
std::vector<std::vector<double>> FuncWrapperND::Jacobian(const std::vector<double>& x) {
    double epsilon = NAN;
    std::size_t N = x.size();
    std::vector<double> r, xp;
    std::vector<std::vector<double>> J(N, std::vector<double>(N, 0));
    std::vector<double> r0 = call(x);
    // Build the Jacobian by column
    for (std::size_t i = 0; i < N; ++i) {
        xp = x;
        epsilon = 0.001 * x[i];
        xp[i] += epsilon;
        r = call(xp);

        for (std::size_t j = 0; j < N; ++j) {
            J[j][i] = (r[j] - r0[j]) / epsilon;
        }
    }
    return J;
}

/**
In this formulation of the Multi-Dimensional Newton-Raphson solver the Jacobian matrix is known.
Therefore, the dx vector can be obtained from

J(x)dx=-f(x)

for a given value of x.  The pointer to the class FuncWrapperND that is passed in must implement the call() and Jacobian()
functions, each of which take the vector x. The data is managed using std::vector<double> vectors

@param f A pointer to an subclass of the FuncWrapperND class that implements the call() and Jacobian() functions
@param x The initial guess value for the solution
@param tol The root-sum-square of the errors from each of the components
@param maxiter The maximum number of iterations
@param w A relaxation multiplier on the step size, multiplying the normal step size
@returns If no errors are found, the solution.  Otherwise, _HUGE, the value for infinity
*/
std::vector<double> NDNewtonRaphson_Jacobian(FuncWrapperND* f, const std::vector<double>& x, double tol, int maxiter, double w) {
    int iter = 0;
    f->errstring.clear();
    std::vector<double> f0, v;
    std::vector<std::vector<double>> JJ;
    std::vector<double> x0 = x;
    Eigen::VectorXd r(x0.size());
    Eigen::MatrixXd J(x0.size(), x0.size());
    double error = 999;
    while (iter == 0 || std::abs(error) > tol) {
        f0 = f->call(x0);
        JJ = f->Jacobian(x0);

        for (std::size_t i = 0; i < x0.size(); ++i) {
            r(i) = f0[i];
            for (std::size_t j = 0; j < x0.size(); ++j) {
                J(i, j) = JJ[i][j];
            }
        }

        Eigen::VectorXd v = J.colPivHouseholderQr().solve(-r);

        // Update the guess
        double max_relchange = -1;
        for (std::size_t i = 0; i < x0.size(); i++) {
            x0[i] += w * v(i);
            double relchange = std::abs(v(i) / x0[i]);
            if (std::abs(x0[i]) > 1e-16 && relchange > max_relchange) {
                max_relchange = relchange;
            }
        }

        // Stop if the solution is not changing by more than numerical precision
        double max_abschange = v.cwiseAbs().maxCoeff();
        if (max_abschange < DBL_EPSILON * 100) {
            return x0;
        }
        if (max_relchange < 1e-12) {
            return x0;
        }
        error = root_sum_square(f0);
        if (iter > maxiter) {
            f->errstring = "reached maximum number of iterations";
            x0[0] = _HUGE;
        }
        iter++;
    }
    return x0;
}

/**
In the newton function, a 1-D Newton-Raphson solver is implemented using exact solutions.  An initial guess for the solution is provided.

@param f A pointer to an instance of the FuncWrapper1D class that implements the call() function
@param x0 The initial guess for the solution
@param ftol The absolute value of the tolerance accepted for the objective function
@param maxiter Maximum number of iterations
@returns If no errors are found, the solution, otherwise the value _HUGE, the value for infinity
*/
double Newton(FuncWrapper1DWithDeriv* f, double x0, double ftol, int maxiter) {
    double x = NAN, dx = NAN, dfdx = NAN, fval = 999;

    // Initialize
    f->iter = 0;
    f->errstring.clear();
    x = x0;
    while (f->iter < maxiter || std::abs(fval) > ftol) {
        fval = f->call(x);
        dfdx = f->deriv(x);
        dx = -fval / dfdx;

        if (!ValidNumber(fval)) {
            f->errstring = "Residual function in newton returned invalid number";
            throw ValueError("Residual function in newton returned invalid number");
        };

        if (f->verbosity > 0) {
            std::cout << format("i: %d, x: %0.15g, dx: %g, f: %g, dfdx: %g", f->iter, x, dx, fval, dfdx) << '\n';
        }

        x += dx;

        if (std::abs(dx / x) < 1e-11) {
            return x;
        }

        if (f->iter > maxiter) {
            f->errstring = "Newton reached maximum number of iterations";
            throw SolutionError(format("Newton reached maximum number of iterations"));
        }
        f->iter++;
    }
    return x;
}
/**
Newton-Raphson iteration safeguarded by a bracket, for a root of f in [xmin, xmax].

f is evaluated at both ends first and must change sign between them, as for
Brent's method. The iteration then keeps an interval on which f changes sign:
every Newton step that would leave that interval (a poor derivative, a turning
point of f) is replaced by a bisection step. So the result is always a root
inside [xmin, xmax], f is never evaluated outside it, and the iteration
converges like Newton near the root. A Newton step that is not at most half
the step before last is also replaced by a bisection, so steps that bounce
around inside the bracket cannot keep it from shrinking. This is the "rtsafe"
scheme of Numerical Recipes. With an approximate derivative the iteration
still converges, but the final step, and so the error of the result, can be
larger than xtol by the factor the derivative is off.

The stopping tolerance is xtol plus four units in the last place of x, so a
tolerance below what a double can resolve at large |x| cannot stall the
iteration between two neighbouring doubles.

A Newton step smaller than the tolerance ends the iteration, as in rtsafe.
That trusts the derivative: one that is far too steep where f is not small
makes the step small away from the root and stops the iteration there. Give
it an exact or at least a reasonable derivative.

@param f A pointer to an instance of the FuncWrapper1DWithDeriv class that implements call() and deriv()
@param x0 The initial guess. If it is not strictly inside (xmin, xmax), NAN for instance, the
iteration starts from the false-position point of the two end values instead.
@param xmin The lower bound for the solution
@param xmax The upper bound for the solution
@param xtol Absolute tolerance on x: the iteration stops once a Newton step, or the
bracketing interval, is smaller than this (plus four ulp of x, see above)
@param maxiter Maximum number of iterations after the two end evaluations
@returns The solution. Throws ValueError if f does not change sign on [xmin, xmax] or returns a
value that is not finite, and SolutionError if it does not converge within maxiter.
*/
double BoundedNewton(FuncWrapper1DWithDeriv* f, double x0, double xmin, double xmax, double xtol, int maxiter) {
    f->iter = 0;
    f->errstring.clear();
    if (!(xmin < xmax)) {
        throw ValueError(format("BoundedNewton needs xmin < xmax, got [%g, %g]", xmin, xmax));
    }
    const double fmin = f->call(xmin);
    const double fmax = f->call(xmax);
    if (!ValidNumber(fmin) || !ValidNumber(fmax)) {
        throw ValueError(format("BoundedNewton: f is not a valid number at the ends of [%g, %g]", xmin, xmax));
    }
    if (fmin == 0) return xmin;
    if (fmax == 0) return xmax;
    if ((fmin > 0) == (fmax > 0)) {
        throw ValueError(format("Inputs in BoundedNewton [%g, %g] do not bracket the root.  Function values are [%g, %g]", xmin, xmax, fmin, fmax));
    }
    // The ends of the bracket on the negative and on the positive side of f
    double x_neg = (fmin < 0) ? xmin : xmax;
    double x_pos = (fmin < 0) ? xmax : xmin;
    double x = (x0 > xmin && x0 < xmax) ? x0 : xmin - fmin * (xmax - xmin) / (fmax - fmin);
    // The false-position point can overflow (inf / inf = NaN) for huge end
    // values, or round onto an end; take the midpoint then
    if (!(x > xmin && x < xmax)) x = 0.5 * xmin + 0.5 * xmax;
    // Size of the step before last. A Newton step larger than half of it is
    // not converging fast enough (it can bounce around a turning point while
    // staying inside the bracket), so a bisection is taken instead.
    double dx_old = xmax - xmin, dx_last = xmax - xmin;
    for (f->iter = 1; f->iter <= maxiter; f->iter++) {
        const double fval = f->call(x);
        if (!ValidNumber(fval)) {
            f->errstring = "Residual function in BoundedNewton returned invalid number";
            throw ValueError(format("Residual function in BoundedNewton returned invalid number at x = %g", x));
        }
        if (fval == 0) return x;
        if (fval < 0) {
            x_neg = x;
        } else {
            x_pos = x;
        }
        const double lo = std::min(x_neg, x_pos), hi = std::max(x_neg, x_pos);
        // Absolute tolerance, but never below what doubles can resolve at x
        const double tol = xtol + 4 * DBL_EPSILON * std::abs(x);
        const double dfdx = f->deriv(x);
        // Bisection, unless the Newton step is usable. Halved before adding,
        // so that lo + hi cannot overflow.
        double x_new = 0.5 * lo + 0.5 * hi;
        if (ValidNumber(dfdx) && dfdx != 0) {
            const double dx = fval / dfdx;
            if (f->verbosity > 0) {
                std::cout << format("i: %d, x: %0.15g, dx: %g, f: %g, dfdx: %g, [%0.15g, %0.15g]", f->iter, x, dx, fval, dfdx, lo, hi) << '\n';
            }
            // Converged. Checked before the interval test: near the root the
            // step can be smaller than one ulp of x, so x - dx == x, which
            // lies on the interval edge and would read as a step outside it.
            // Only a step that stays inside the bracket counts: with a derivative
            // of the wrong sign a tiny step can point away from the root, and
            // the bisection below takes over instead.
            if (std::abs(dx) < tol && x - dx >= lo && x - dx <= hi) return x - dx;
            if (x - dx > lo && x - dx < hi && std::abs(dx) <= 0.5 * std::abs(dx_old)) x_new = x - dx;
        }
        // The bracket is as small as asked for, or as small as doubles allow
        // (the midpoint rounds onto x itself)
        if (hi - lo < tol || x_new == x) return x_new;
        dx_old = dx_last;
        dx_last = x_new - x;
        x = x_new;
    }
    f->errstring = "reached maximum number of iterations";
    throw SolutionError(format("BoundedNewton reached maximum number of iterations (%d)", maxiter));
}

/**
In the Halley's method solver, two derivatives of the input variable are needed, it yields the following method:

\f[
x_{n+1} = x_n - \frac {2 f(x_n) f'(x_n)} {2 {[f'(x_n)]}^2 - f(x_n) f''(x_n)}
\f]

http://en.wikipedia.org/wiki/Halley%27s_method

@param f A pointer to an instance of the FuncWrapper1DWithTwoDerivs class that implements the call() and two derivatives
@param x0 The initial guess for the solution
@param ftol The absolute value of the tolerance accepted for the objective function
@param maxiter Maximum number of iterations
@param xtol_rel The minimum allowable (relative) step size
@returns If no errors are found, the solution, otherwise the value _HUGE, the value for infinity
*/
double Halley(FuncWrapper1DWithTwoDerivs* f, double x0, double ftol, int maxiter, double xtol_rel) {
    double x = NAN, dx = NAN, fval = 999, dfdx = NAN, d2fdx2 = NAN;

    // Initialize
    f->iter = 0;
    f->errstring.clear();
    x = x0;

    // The relaxation factor (less than 1 for smaller steps)
    double omega = f->options.get_double("omega", 1.0);

    while (f->iter < 2 || std::abs(fval) > ftol) {
        if (f->input_not_in_range(x)) {
            std::string msg = format("Input [%0.15g] is out of range", x);
            f->errstring = msg;
            throw ValueError(msg);
        }

        fval = f->call(x);
        dfdx = f->deriv(x);
        d2fdx2 = f->second_deriv(x);

        if (!ValidNumber(fval)) {
            f->errstring = "Residual function in Halley returned invalid number";
            throw ValueError("Residual function in Halley returned invalid number");
        };
        if (!ValidNumber(dfdx)) {
            f->errstring = "Derivative function in Halley returned invalid number";
            throw ValueError("Derivative function in Halley returned invalid number");
        };

        dx = -omega * (2 * fval * dfdx) / (2 * POW2(dfdx) - fval * d2fdx2);

        if (f->verbosity > 0) {
            std::cout << format("i: %d, x: %0.15g, dx: %g, f: %g, dfdx: %g, d2fdx2: %g", f->iter, x, dx, fval, dfdx, d2fdx2) << '\n';
        }

        x += dx;

        if (std::abs(dx / x) < xtol_rel) {
            return x;
        }

        if (f->iter > maxiter) {
            f->errstring = "Halley reached maximum number of iterations";
            throw SolutionError(format("Halley reached maximum number of iterations"));
        }
        f->iter += 1;
    }
    return x;
}

/**
 In the 4-th order Householder method, three derivatives of the input variable are needed, it yields the following method:

 \f[
 x_{n+1} = x_n - f(x_n)\left( \frac {[f'(x_n)]^2 - f(x_n)f''(x_n)/2  } {[f'(x_n)]^3-f(x_n)f'(x_n)f''(x_n)+f'''(x_n)*[f(x_n)]^2/6 } \right)
 \f]

http://numbers.computation.free.fr/Constants/Algorithms/newton.ps

 @param f A pointer to an instance of the FuncWrapper1DWithThreeDerivs class that implements the call() and three derivatives
 @param x0 The initial guess for the solution
 @param ftol The absolute value of the tolerance accepted for the objective function
 @param maxiter Maximum number of iterations
 @param xtol_rel The minimum allowable (relative) step size
 @returns If no errors are found, the solution, otherwise the value _HUGE, the value for infinity
 */
double Householder4(FuncWrapper1DWithThreeDerivs* f, double x0, double ftol, int maxiter, double xtol_rel) {
    double x = NAN, dx = NAN, fval = 999, dfdx = NAN, d2fdx2 = NAN, d3fdx3 = NAN;

    // Initialization
    f->iter = 1;
    f->errstring.clear();
    x = x0;

    // The relaxation factor (less than 1 for smaller steps)
    double omega = f->options.get_double("omega", 1.0);

    while (f->iter < 2 || std::abs(fval) > ftol) {
        if (f->input_not_in_range(x)) {
            throw ValueError(format("Input [%g] is out of range", x));
        }

        fval = f->call(x);
        dfdx = f->deriv(x);
        d2fdx2 = f->second_deriv(x);
        d3fdx3 = f->third_deriv(x);

        if (!ValidNumber(fval)) {
            throw ValueError("Residual function in Householder4 returned invalid number");
        };
        if (!ValidNumber(dfdx)) {
            throw ValueError("Derivative function in Householder4 returned invalid number");
        };
        if (!ValidNumber(d2fdx2)) {
            throw ValueError("Second derivative function in Householder4 returned invalid number");
        };
        if (!ValidNumber(d3fdx3)) {
            throw ValueError("Third derivative function in Householder4 returned invalid number");
        };

        dx = -omega * fval * (POW2(dfdx) - fval * d2fdx2 / 2.0) / (POW3(dfdx) - fval * dfdx * d2fdx2 + d3fdx3 * POW2(fval) / 6.0);

        x += dx;

        if (std::abs(dx / x) < xtol_rel) {
            return x;
        }

        if (f->iter > maxiter) {
            f->errstring = "reached maximum number of iterations";
            throw SolutionError(format("Householder4 reached maximum number of iterations"));
        }
        f->iter += 1;
    }
    return x;
}

/**
In the secant function, a 1-D Newton-Raphson solver is implemented.  An initial guess for the solution is provided.

@param f A pointer to an instance of the FuncWrapper1D class that implements the call() function
@param x0 The initial guess for the solutionh
@param dx The initial amount that is added to x in order to build the numerical derivative
@param tol The absolute value of the tolerance accepted for the objective function
@param maxiter Maximum number of iterations
@returns If no errors are found, the solution, otherwise the value _HUGE, the value for infinity
*/
double Secant(FuncWrapper1D* f, double x0, double dx, double tol, int maxiter) {
#if defined(COOLPROP_DEEP_DEBUG)
    static std::vector<double> xlog, flog;
    xlog.clear();
    flog.clear();
#endif

    // Initialization
    double x1 = 0, x2 = 0, x3 = 0, y1 = 0, y2 = 0, x = x0, fval = 999;
    f->iter = 1;
    f->errstring.clear();

    // The relaxation factor (less than 1 for smaller steps)
    double omega = f->options.get_double("omega", 1.0);

    if (std::abs(dx) == 0) {
        f->errstring = "dx cannot be zero";
        return _HUGE;
    }
    while (f->iter <= 2 || std::abs(fval) > tol) {
        if (f->iter == 1) {
            x1 = x0;
            x = x1;
        }
        if (f->iter == 2) {
            x2 = x0 + dx;
            x = x2;
        }
        if (f->iter > 2) {
            x = x2;
        }

        if (f->input_not_in_range(x)) {
            throw ValueError(format("Input [%g] is out of range", x));
        }

        fval = f->call(x);

#if defined(COOLPROP_DEEP_DEBUG)
        xlog.push_back(x);
        flog.push_back(fval);
#endif

        if (!ValidNumber(fval)) {
            throw ValueError("Residual function in secant returned invalid number");
        };
        if (f->iter == 1) {
            y1 = fval;
        }
        if (f->iter > 1) {
            double deltax = x2 - x1;
            if (std::abs(deltax) < 1e-14) {
                return x;
            }
            y2 = fval;
            double deltay = y2 - y1;
            if (f->iter > 2 && std::abs(deltay) < 1e-14) {
                return x;
            }
            x3 = x2 - omega * y2 / (y2 - y1) * (x2 - x1);
            y1 = y2;
            x1 = x2;
            x2 = x3;
        }
        if (f->iter > maxiter) {
            f->errstring = std::string("reached maximum number of iterations");
            throw SolutionError(format("Secant reached maximum number of iterations"));
        }
        f->iter += 1;
    }
    return x3;
}

/**
In the secant function, a 1-D Newton-Raphson solver is implemented.  An initial guess for the solution is provided.

@param f A pointer to an instance of the FuncWrapper1D class that implements the call() function
@param x0 The initial guess for the solution
@param xmax The upper bound for the solution
@param xmin The lower bound for the solution
@param dx The initial amount that is added to x in order to build the numerical derivative
@param tol The absolute value of the tolerance accepted for the objective function
@param maxiter Maximum number of iterations
@returns If no errors are found, the solution, otherwise the value _HUGE, the value for infinity
*/
double BoundedSecant(FuncWrapper1D* f, double x0, double xmin, double xmax, double dx, double tol, int maxiter) {
    double x1 = 0, x2 = 0, x3 = 0, y1 = 0, y2 = 0, x = NAN, fval = 999;
    int iter = 1;
    f->errstring.clear();
    if (std::abs(dx) == 0) {
        f->errstring = "dx cannot be zero";
        return _HUGE;
    }
    while (iter <= 3 || std::abs(fval) > tol) {
        if (iter == 1) {
            x1 = x0;
            x = x1;
        } else if (iter == 2) {
            x2 = x0 + dx;
            x = x2;
        } else {
            x = x2;
        }
        fval = f->call(x);
        if (iter == 1) {
            y1 = fval;
        } else {
            y2 = fval;
            x3 = x2 - y2 / (y2 - y1) * (x2 - x1);
            // Check bounds, go half the way to the limit if limit is exceeded
            if (x3 < xmin) {
                x3 = (xmin + x2) / 2;
            }
            if (x3 > xmax) {
                x3 = (xmax + x2) / 2;
            }
            y1 = y2;
            x1 = x2;
            x2 = x3;
        }
        if (iter > maxiter) {
            f->errstring = "reached maximum number of iterations";
            throw SolutionError(format("BoundedSecant reached maximum number of iterations"));
        }
        iter = iter + 1;
    }
    f->errcode = 0;
    return x3;
}

/**
In the secant function, a 1-D Newton-Raphson solver is implemented.  An initial guess for the solution is provided.
Note that this is different than the Secant function because if something goes out of bounds, it will just make its best guess.

@param f A pointer to an instance of the FuncWrapper1D class that implements the call() function
@param x0 The initial guess for the solutionh
@param dx The initial amount that is added to x in order to build the numerical derivative
@param tol The absolute value of the tolerance accepted for the objective function
@param maxiter Maximum number of iterations
@returns If no errors are found, the solution, otherwise the value _HUGE, the value for infinity
*/
double ExtrapolatingSecant(FuncWrapper1D* f, double x0, double dx, double tol, int maxiter) {
#if defined(COOLPROP_DEEP_DEBUG)
    static std::vector<double> xlog, flog;
    xlog.clear();
    flog.clear();
#endif

    // Initialization
    double x1 = 0, x2 = 0, x3 = 0, y0 = 0, y1 = 0, y2 = 0, x = x0, fval = 999;
    f->iter = 1;
    f->errstring.clear();

    // The relaxation factor (less than 1 for smaller steps)
    double omega = f->options.get_double("omega", 1.0);

    if (std::abs(dx) == 0) {
        f->errstring = "dx cannot be zero";
        return _HUGE;
    }
    while (f->iter <= 2 || std::abs(fval) > tol) {
        if (f->iter == 1) {
            x1 = x0;
            x = x1;
        }
        if (f->iter == 2) {
            x2 = x0 + dx;
            x = x2;
        }
        if (f->iter > 2) {
            x = x2;
        }

        if (f->input_not_in_range(x)) {
            throw ValueError(format("Input [%g] is out of range", x));
        }

        fval = f->call(x);

#if defined(COOLPROP_DEEP_DEBUG)
        xlog.push_back(x);
        flog.push_back(fval);
#endif

        if (!ValidNumber(fval)) {
            if (f->iter == 1) {
                return x;
            } else {
                return x2 - omega * y1 / (y1 - y0) * (x2 - x1);
            }
        };
        if (f->iter == 1) {
            y1 = fval;
        }
        if (f->iter > 1) {
            double deltax = x2 - x1;
            if (std::abs(deltax) < 1e-14) {
                return x;
            }
            y2 = fval;
            double deltay = y2 - y1;
            if (f->iter > 2 && std::abs(deltay) < 1e-14) {
                return x;
            }
            x3 = x2 - omega * y2 / (y2 - y1) * (x2 - x1);
            y0 = y1;
            y1 = y2;
            x1 = x2;
            x2 = x3;
        }
        if (f->iter > maxiter) {
            f->errstring = std::string("reached maximum number of iterations");
            throw SolutionError(format("Secant reached maximum number of iterations"));
        }
        f->iter += 1;
    }
    return x3;
}

/**

This function implements a 1-D bounded solver using the algorithm from Brent, R. P., Algorithms for Minimization Without Derivatives.
Englewood Cliffs, NJ: Prentice-Hall, 1973. Ch. 3-4.

a and b must bound the solution of interest and f(a) and f(b) must have opposite signs.  If the function is continuous, there must be
at least one solution in the interval [a,b].

@param f A pointer to an instance of the FuncWrapper1D class that must implement the class() function
@param a The minimum bound for the solution of f=0
@param b The maximum bound for the solution of f=0
@param macheps The machine precision
@param t Tolerance (absolute)
@param maxiter Maximum number of steps allowed.  Will throw a SolutionError if the solution cannot be found
*/
double Brent(FuncWrapper1D* f, double a, double b, double macheps, double t, int maxiter) {
    int iter = 0;
    f->errstring.clear();
    double fa = NAN, fb = NAN, c = NAN, fc = NAN, m = NAN, tol = NAN, d = NAN, e = NAN, p = NAN, q = NAN, s = NAN, r = NAN;
    fa = f->call(a);
    fb = f->call(b);

    // If one of the boundaries is to within tolerance, just stop
    if (std::abs(fb) < t) {
        return b;
    }
    if (!ValidNumber(fb)) {
        throw ValueError(format("Brent's method f(b) is NAN for b = %g, other input was a = %g", b, a).c_str());
    }
    if (std::abs(fa) < t) {
        return a;
    }
    if (!ValidNumber(fa)) {
        throw ValueError(format("Brent's method f(a) is NAN for a = %g, other input was b = %g", a, b).c_str());
    }
    if (fa * fb > 0) {
        throw ValueError(format("Inputs in Brent [%f,%f] do not bracket the root.  Function values are [%f,%f]", a, b, fa, fb));
    }

    c = a;
    fc = fa;
    iter = 1;
    if (std::abs(fc) < std::abs(fb)) {
        // Goto ext: from Brent ALGOL code
        a = b;
        b = c;
        c = a;
        fa = fb;
        fb = fc;
        fc = fa;
    }
    d = b - a;
    e = b - a;
    m = 0.5 * (c - b);
    tol = 2 * macheps * std::abs(b) + t;
    while (std::abs(m) > tol && fb != 0) {
        // See if a bisection is forced
        if (std::abs(e) < tol || std::abs(fa) <= std::abs(fb)) {
            m = 0.5 * (c - b);
            d = e = m;
        } else {
            s = fb / fa;
            if (a == c) {
                //Linear interpolation
                p = 2 * m * s;
                q = 1 - s;
            } else {
                //Inverse quadratic interpolation
                q = fa / fc;
                r = fb / fc;
                m = 0.5 * (c - b);
                p = s * (2 * m * q * (q - r) - (b - a) * (r - 1));
                q = (q - 1) * (r - 1) * (s - 1);
            }
            if (p > 0) {
                q = -q;
            } else {
                p = -p;
            }
            s = e;
            e = d;
            m = 0.5 * (c - b);
            if (2 * p < 3 * m * q - std::abs(tol * q) || p < std::abs(0.5 * s * q)) {
                d = p / q;
            } else {
                m = 0.5 * (c - b);
                d = e = m;
            }
        }
        a = b;
        fa = fb;
        if (std::abs(d) > tol) {
            b += d;
        } else if (m > 0) {
            b += tol;
        } else {
            b += -tol;
        }
        fb = f->call(b);
        if (!ValidNumber(fb)) {
            throw ValueError(format("Brent's method f(t) is NAN for t = %g", b).c_str());
        }
        if (std::abs(fb) < macheps) {
            return b;
        }
        if (fb * fc > 0) {
            // Goto int: from Brent ALGOL code
            c = a;
            fc = fa;
            d = e = b - a;
        }
        if (std::abs(fc) < std::abs(fb)) {
            // Goto ext: from Brent ALGOL code
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        m = 0.5 * (c - b);
        tol = 2 * macheps * std::abs(b) + t;
        iter += 1;
        if (!ValidNumber(a)) {
            throw ValueError(format("Brent's method a is NAN").c_str());
        }
        if (!ValidNumber(b)) {
            throw ValueError(format("Brent's method b is NAN").c_str());
        }
        if (!ValidNumber(c)) {
            throw ValueError(format("Brent's method c is NAN").c_str());
        }
        if (iter > maxiter) {
            throw SolutionError(format("Brent's method reached maximum number of steps of %d ", maxiter));
        }
        if (std::abs(fb) < 2 * macheps * std::abs(b)) {
            return b;
        }
    }
    return b;
}

}; /* namespace CoolProp */

#ifdef ENABLE_CATCH
#    include <algorithm>
#    include <functional>
#    include <catch2/catch_all.hpp>

namespace {
/// Test residual: f, its derivative, and the range of x it was evaluated at
class TestResidual : public CoolProp::FuncWrapper1DWithDeriv
{
   public:
    std::function<double(double)> f, dfdx;
    double x_lowest = HUGE_VAL, x_highest = -HUGE_VAL;
    int n_calls = 0;
    TestResidual(std::function<double(double)> f, std::function<double(double)> dfdx) : f(std::move(f)), dfdx(std::move(dfdx)) {}
    double call(double x) override {
        x_lowest = std::min(x_lowest, x);
        x_highest = std::max(x_highest, x);
        ++n_calls;
        return f(x);
    }
    double deriv(double x) override {
        return dfdx(x);
    }
};
}  // namespace

TEST_CASE("BoundedNewton finds the root inside its bracket", "[Solvers]") {
    SECTION("Increasing function with the exact derivative") {
        TestResidual r([](double x) { return x * x * x - 2.0; }, [](double x) { return 3 * x * x; });
        const double x = CoolProp::BoundedNewton(&r, NAN, 0.5, 3.0, 1e-14, 100);
        CHECK(x == Catch::Approx(std::cbrt(2.0)).epsilon(1e-14));
        CHECK(r.n_calls < 12);
    }
    SECTION("Decreasing function, root found from either orientation") {
        TestResidual r([](double x) { return std::exp(-x) - 0.1; }, [](double x) { return -std::exp(-x); });
        CHECK(CoolProp::BoundedNewton(&r, 1.0, 0.0, 10.0, 1e-14, 100) == Catch::Approx(std::log(10.0)).epsilon(1e-14));
    }
    SECTION("A poor derivative still converges, and never leaves the interval") {
        TestResidual r([](double x) { return x * x * x - 2.0; }, [](double x) { return 30 * x * x; });  // 10 times too steep
        const double x = CoolProp::BoundedNewton(&r, NAN, 0.5, 3.0, 1e-13, 200);
        CHECK(x == Catch::Approx(std::cbrt(2.0)).epsilon(1e-12));
        CHECK(r.x_lowest >= 0.5);
        CHECK(r.x_highest <= 3.0);
    }
    SECTION("A turning point inside the interval does not lead it out of the bracket") {
        // f = (x - 1)^3 - (x - 1) + 0.5 has one root in [-1, 3] near x = -0.19 and
        // a local maximum and minimum in between, where Newton steps go astray.
        TestResidual r([](double x) { return std::pow(x - 1, 3) - (x - 1) + 0.5; }, [](double x) { return 3 * (x - 1) * (x - 1) - 1; });
        const double x = CoolProp::BoundedNewton(&r, 1.2, -1.0, 3.0, 1e-14, 200);
        CHECK(std::abs(std::pow(x - 1, 3) - (x - 1) + 0.5) < 1e-12);
        CHECK(r.x_lowest >= -1.0);
        CHECK(r.x_highest <= 3.0);
    }
    SECTION("No sign change on the interval throws instead of guessing") {
        TestResidual r([](double x) { return (x - 1) * (x - 1) - 0.25; }, [](double x) { return 2 * (x - 1); });  // roots 0.5 and 1.5
        CHECK_THROWS_AS(CoolProp::BoundedNewton(&r, NAN, 0.0, 2.0, 1e-14, 100), CoolProp::ValueError);
    }
    SECTION("An absolute tolerance below the spacing of doubles does not stall it") {
        // Above 2048 neighbouring doubles are 4.5e-13 apart, twice the
        // 2.2e-13 the incompressible flashes ask for. f jumps by 2e-12 at the
        // root, so |f| never drops below 1e-12, as when the rounding noise of
        // f (h at a few MJ/kg) is larger than the tolerance. The bracket then
        // ends as two neighbouring doubles and cannot shrink any further; the
        // iteration must stop there instead of running out of iterations.
        const double root = 2100.123456789;
        auto noisy = [root](double x) { return (x - root) + std::copysign(1e-12, x - root); };
        TestResidual r(noisy, [](double) { return 1.0; });
        double x = NAN;
        CHECK_NOTHROW(x = CoolProp::BoundedNewton(&r, NAN, 400.0, 2500.0, 1e3 * DBL_EPSILON, 100));
        CHECK(std::abs(x - root) < 1e-11);
    }
    SECTION("A root at an end of the interval is returned") {
        TestResidual r([](double x) { return x - 2.0; }, [](double x) { return 1.0; });
        CHECK(CoolProp::BoundedNewton(&r, NAN, 2.0, 5.0, 1e-14, 100) == 2.0);
    }
    SECTION("Huge end values do not turn the first point into NaN") {
        // fmax - fmin overflows to infinity, and the false-position point to NaN
        TestResidual r([](double x) { return 1e308 * ((x - 5.0) / 5.0); }, [](double) { return 1e308 / 5.0; });
        CHECK(CoolProp::BoundedNewton(&r, NAN, 0.0, 10.0, 1e-14, 100) == Catch::Approx(5.0).epsilon(1e-14));
        CHECK(r.x_lowest >= 0.0);
        CHECK(r.x_highest <= 10.0);
    }
    SECTION("Bisection between huge bounds does not overflow") {
        // A zero derivative forces bisection; lo + hi would overflow to infinity
        TestResidual r([](double x) { return x / 1.2e308 - 1.0; }, [](double) { return 0.0; });
        double x = NAN;
        CHECK_NOTHROW(x = CoolProp::BoundedNewton(&r, 1.1e308, 1e308, 1.5e308, 1e-14, 2000));
        CHECK(x == Catch::Approx(1.2e308).epsilon(1e-14));
        CHECK(r.x_highest <= 1.5e308);
    }
}
#endif /* ENABLE_CATCH */
