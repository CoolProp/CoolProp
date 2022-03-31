#include <vector>
#include "Solvers.h"
#include "math.h"
#include "MatrixMath.h"
#include <iostream>
#include "CoolPropTools.h"
#include <Eigen/Dense>

namespace CoolProp {

/** \brief Calculate the Jacobian using numerical differentiation by column
 */
std::vector<std::vector<double>> FuncWrapperND::Jacobian(const std::vector<double>& x) {
    double epsilon;
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
@param x0 The initial guess value for the solution
@param tol The root-sum-square of the errors from each of the components
@param maxiter The maximum number of iterations
@param errstring  A string with the returned error.  If the length of errstring is zero, no errors were found
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
    double x, dx, fval = 999;
    int iter = 1;
    f->errstring.clear();
    x = x0;
    while (iter < 2 || std::abs(fval) > ftol) {
        fval = f->call(x);
        dx = -fval / f->deriv(x);

        if (!ValidNumber(fval)) {
            throw ValueError("Residual function in newton returned invalid number");
        };

        x += dx;

        if (std::abs(dx / x) < 1e-11) {
            return x;
        }

        if (iter > maxiter) {
            f->errstring = "reached maximum number of iterations";
            throw SolutionError(format("Newton reached maximum number of iterations"));
        }
        iter = iter + 1;
    }
    return x;
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
    double x, dx, fval = 999, dfdx, d2fdx2;

    // Initialize
    f->iter = 0;
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

        if (!ValidNumber(fval)) {
            throw ValueError("Residual function in Halley returned invalid number");
        };
        if (!ValidNumber(dfdx)) {
            throw ValueError("Derivative function in Halley returned invalid number");
        };

        dx = -omega * (2 * fval * dfdx) / (2 * POW2(dfdx) - fval * d2fdx2);

        x += dx;

        if (std::abs(dx / x) < xtol_rel) {
            return x;
        }

        if (f->iter > maxiter) {
            f->errstring = "reached maximum number of iterations";
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
    double x, dx, fval = 999, dfdx, d2fdx2, d3fdx3;

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
    double x1 = 0, x2 = 0, x3 = 0, y1 = 0, y2 = 0, x, fval = 999;
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
    int iter;
    f->errstring.clear();
    double fa, fb, c, fc, m, tol, d, e, p, q, s, r;
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
