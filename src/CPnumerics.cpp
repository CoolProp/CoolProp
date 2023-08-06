#include "CPnumerics.h"
#include "MatrixMath.h"
#include <unsupported/Eigen/Polynomials>

double root_sum_square(const std::vector<double>& x) {
    double sum = 0;
    for (unsigned int i = 0; i < x.size(); i++) {
        sum += pow(x[i], 2);
    }
    return sqrt(sum);
}
double interp1d(const std::vector<double>* x, const std::vector<double>* y, double x0) {
    std::size_t i, L, R, M;
    L = 0;
    R = (*x).size() - 1;
    M = (L + R) / 2;
    // Use interval halving to find the indices which bracket the density of interest
    while (R - L > 1) {
        if (x0 >= (*x)[M]) {
            L = M;
            M = (L + R) / 2;
            continue;
        }
        if (x0 < (*x)[M]) {
            R = M;
            M = (L + R) / 2;
            continue;
        }
    }
    i = L;
    if (i < (*x).size() - 2) {
        // Go "forwards" with the interpolation range
        return QuadInterp((*x)[i], (*x)[i + 1], (*x)[i + 2], (*y)[i], (*y)[i + 1], (*y)[i + 2], x0);
    } else {
        // Go "backwards" with the interpolation range
        return QuadInterp((*x)[i], (*x)[i - 1], (*x)[i - 2], (*y)[i], (*y)[i - 1], (*y)[i - 2], x0);
    }
}
double powInt(double x, int y) {
    // Raise a double to an integer power
    // Overload not provided in math.h
    int i;
    double product = 1.0;
    double x_in;
    int y_in;

    if (y == 0) {
        return 1.0;
    }

    if (y < 0) {
        x_in = 1 / x;
        y_in = -y;
    } else {
        x_in = x;
        y_in = y;
    }

    if (y_in == 1) {
        return x_in;
    }

    product = x_in;
    for (i = 1; i < y_in; i++) {
        product = product * x_in;
    }
    return product;
}

void MatInv_2(double A[2][2], double B[2][2]) {
    double Det;
    //Using Cramer's Rule to solve

    Det = A[0][0] * A[1][1] - A[1][0] * A[0][1];
    B[0][0] = 1.0 / Det * A[1][1];
    B[1][1] = 1.0 / Det * A[0][0];
    B[1][0] = -1.0 / Det * A[1][0];
    B[0][1] = -1.0 / Det * A[0][1];
}

void solve_cubic(double a, double b, double c, double d, int& N, double& x0, double& x1, double& x2) {
    // 0 = ax^3 + b*x^2 + c*x + d

    // First check if the "cubic" is actually a second order or first order curve
    if (std::abs(a) < 10 * DBL_EPSILON) {
        if (std::abs(b) < 10 * DBL_EPSILON) {
            // Linear solution if a = 0 and b = 0
            x0 = -d / c;
            N = 1;
            return;
        } else {
            // Quadratic solution(s) if a = 0 and b != 0
            x0 = (-c + sqrt(c * c - 4 * b * d)) / (2 * b);
            x1 = (-c - sqrt(c * c - 4 * b * d)) / (2 * b);
            N = 2;
            return;
        }
    }

    // Ok, it is really a cubic

    // Discriminant
    double DELTA = 18 * a * b * c * d - 4 * b * b * b * d + b * b * c * c - 4 * a * c * c * c - 27 * a * a * d * d;
    // Coefficients for the depressed cubic t^3+p*t+q = 0
    double p = (3 * a * c - b * b) / (3 * a * a);
    double q = (2 * b * b * b - 9 * a * b * c + 27 * a * a * d) / (27 * a * a * a);

    if (DELTA < 0) {
        // One real root
        double t0;
        if (4 * p * p * p + 27 * q * q > 0 && p < 0) {
            t0 = -2.0 * std::abs(q) / q * sqrt(-p / 3.0) * cosh(1.0 / 3.0 * acosh(-3.0 * std::abs(q) / (2.0 * p) * sqrt(-3.0 / p)));
        } else {
            t0 = -2.0 * sqrt(p / 3.0) * sinh(1.0 / 3.0 * asinh(3.0 * q / (2.0 * p) * sqrt(3.0 / p)));
        }
        N = 1;
        x0 = t0 - b / (3 * a);
        x1 = t0 - b / (3 * a);
        x2 = t0 - b / (3 * a);
    } else  //(DELTA>0)
    {
        // Three real roots
        double t0 = 2.0 * sqrt(-p / 3.0) * cos(1.0 / 3.0 * acos(3.0 * q / (2.0 * p) * sqrt(-3.0 / p)) - 0 * 2.0 * M_PI / 3.0);
        double t1 = 2.0 * sqrt(-p / 3.0) * cos(1.0 / 3.0 * acos(3.0 * q / (2.0 * p) * sqrt(-3.0 / p)) - 1 * 2.0 * M_PI / 3.0);
        double t2 = 2.0 * sqrt(-p / 3.0) * cos(1.0 / 3.0 * acos(3.0 * q / (2.0 * p) * sqrt(-3.0 / p)) - 2 * 2.0 * M_PI / 3.0);

        N = 3;
        x0 = t0 - b / (3 * a);
        x1 = t1 - b / (3 * a);
        x2 = t2 - b / (3 * a);
    }
}
void solve_quartic(double a, double b, double c, double d, double e, int& N, double& x0, double& x1, double& x2, double& x3) {

    // 0 = ax^4 + b*x^3 + c*x^2 + d*x + e

    Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
    Eigen::VectorXd coeff(5);
    coeff << e, d, c, b, a;
    solver.compute(coeff);

    std::vector<double> realRoots;
    solver.realRoots(realRoots);
    N = static_cast<int>(realRoots.size());

    if (N > 0) {
        x0 = realRoots[0];
    }
    if (N > 1) {
        x1 = realRoots[1];
    }
    if (N > 2) {
        x2 = realRoots[2];
    }
    if (N > 3) {
        x3 = realRoots[3];
    }
}

bool SplineClass::build() {
    if (Nconstraints == 4) {
        std::vector<double> abcd = CoolProp::linsolve(A, B);
        a = abcd[0];
        b = abcd[1];
        c = abcd[2];
        d = abcd[3];
        return true;
    } else {
        throw CoolProp::ValueError(format("Number of constraints[%d] is not equal to 4", Nconstraints));
    }
}
bool SplineClass::add_value_constraint(double x, double y) {
    int i = Nconstraints;
    if (i == 4) return false;
    A[i][0] = x * x * x;
    A[i][1] = x * x;
    A[i][2] = x;
    A[i][3] = 1;
    B[i] = y;
    Nconstraints++;
    return true;
}
void SplineClass::add_4value_constraints(double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4) {
    add_value_constraint(x1, y1);
    add_value_constraint(x2, y2);
    add_value_constraint(x3, y3);
    add_value_constraint(x4, y4);
}
bool SplineClass::add_derivative_constraint(double x, double dydx) {
    int i = Nconstraints;
    if (i == 4) return false;
    A[i][0] = 3 * x * x;
    A[i][1] = 2 * x;
    A[i][2] = 1;
    A[i][3] = 0;
    B[i] = dydx;
    Nconstraints++;
    return true;
}
double SplineClass::evaluate(double x) {
    return a * x * x * x + b * x * x + c * x + d;
}
