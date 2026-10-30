
#include "CoolProp/CoolProp.h"
#include "CoolProp/Configuration.h"
#include "CoolProp/DataStructures.h"
#include "CoolProp/fluids/IncompressibleFluid.h"
#include "IncompressibleLibrary.h"
#include <cmath>
#include "CoolProp/numerics/MatrixMath.h"
#include "CoolProp/numerics/PolyMath.h"
#include "CoolProp/numerics/Solvers.h"
// For the Chebyshev-Lobatto interpolation matrix used at load time (validate)
// to re-expand cp(T)/T; deliberately included in this .cpp only, the public
// header carries plain Eigen members.
#include "CoolProp/superancillary/superancillary.h"
#include <Eigen/Core>

constexpr double INCOMP_EPSILON = DBL_EPSILON * 100.0;
constexpr double INCOMP_DELTA = INCOMP_EPSILON * 10.0;

namespace {
/// The exponential-type fits (baseExponential, baseLogexponential) have a real,
/// non-removable pole where their denominator x_den crosses zero: the
/// function genuinely diverges there, unlike the removable 0/0 that shows up
/// when differentiating an ordinary polynomial (see Polynomial2DFrac::evaluate
/// in PolyMath.cpp, which can resolve that case exactly instead of
/// approximating). Away from the pole, evaluate fnc directly; within a tiny
/// epsilon band around it, linearly interpolate between two points just
/// outside the band rather than evaluating at (or near) the singularity.
template <typename Func>
double evaluateAwayFromPole(double x_den, Func fnc) {
    double x_lo = -INCOMP_EPSILON;
    double x_hi = +INCOMP_EPSILON;
    if (x_den < x_lo || x_den > x_hi) {
        return fnc(x_den);
    }
    x_lo -= INCOMP_DELTA;
    x_hi += INCOMP_DELTA;
    const double f_lo = fnc(x_lo);
    const double f_hi = fnc(x_hi);
    return (f_hi - f_lo) / (x_hi - x_lo) * (x_den - x_lo) + f_lo;
}

// ---- Polynomial caloric fits (type == INCOMPRESSIBLE_POLYNOMIAL) -----------
//
// Rows multiply powers of (T - Tbase), columns multiply powers of
// (x - xbase), as in Polynomial2DFrac::evaluate. The derivative and the two
// integrals that h and s need are built once at load time (validate) as
// coefficient matrices in this same basis, so evaluating them costs no more
// than evaluating the fit itself.

/// Evaluate a coefficient matrix in the polynomial basis with Horner's
/// scheme: first along each row in (x - xbase), then down the rows in
/// (T - Tbase). The highest coefficient is the starting value rather than
/// 0 * dx, so a single-column (pure fluid) matrix never touches xbase.
double polyEvalMatrix(const Eigen::MatrixXd& M, double T, double x, double Tbase, double xbase) {
    if (M.rows() < 1 || M.cols() < 1) {
        // Same reasoning as chebEvalMatrix below: an empty derived matrix
        // must be loud, not a silent 0.0.
        throw CoolProp::ValueError("Polynomial coefficient matrix is empty; the fluid was not loaded through validate().");
    }
    const double dT = T - Tbase;
    const double dx = x - xbase;
    auto row_value = [&](Eigen::Index i) {
        double value = M(i, M.cols() - 1);
        for (Eigen::Index j = M.cols() - 2; j >= 0; --j) {
            value = value * dx + M(i, j);
        }
        return value;
    };
    double result = row_value(M.rows() - 1);
    for (Eigen::Index i = M.rows() - 2; i >= 0; --i) {
        result = result * dT + row_value(i);
    }
    return result;
}

/// Coefficients of d/dT: the term a_i (T - Tbase)^i becomes
/// i a_i (T - Tbase)^(i-1), so row i moves to row i-1 times i.
Eigen::MatrixXd polyDerivativeMatrix(const Eigen::MatrixXd& C) {
    const Eigen::Index n = C.rows();
    if (n <= 1) return Eigen::MatrixXd::Zero(1, C.cols());
    Eigen::MatrixXd D(n - 1, C.cols());
    for (Eigen::Index i = 1; i < n; ++i) {
        D.row(i - 1) = static_cast<double>(i) * C.row(i);
    }
    return D;
}

/// Coefficients of the indefinite integral over T: a_i (T - Tbase)^i
/// becomes a_i / (i+1) (T - Tbase)^(i+1). Row 0, the integration constant,
/// is left at zero.
Eigen::MatrixXd polyAntiderivativeMatrix(const Eigen::MatrixXd& C) {
    const Eigen::Index n = C.rows();
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(n + 1, C.cols());
    for (Eigen::Index i = 0; i < n; ++i) {
        B.row(i + 1) = C.row(i) / static_cast<double>(i + 1);
    }
    return B;
}

/// Coefficients of the indefinite integral Int f/T dT for a polynomial f.
/// With y = T - Tbase, each term y^i / T = y^i / (y + Tbase) splits by
/// polynomial division into
///     y^i / (y + Tbase) = sum_{m=0}^{i-1} (-Tbase)^m y^(i-1-m) + (-Tbase)^i / T
/// and so integrates to
///     sum_{m=0}^{i-1} (-Tbase)^m y^(i-m) / (i-m)  +  (-Tbase)^i ln(T).
/// The first part goes into P (same basis as C), the factor of ln(T) into
/// the single row L. This holds for Tbase == 0 as well, where it reduces to
/// a_0 ln(T) + sum_{i>0} a_i T^i / i. Unlike the per-call evaluation in
/// Polynomial2DFrac::integral, which expands in powers of T itself, the
/// polynomial part stays in powers of (T - Tbase) and has no pole anywhere.
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> polyDivideByTAntiderivative(const Eigen::MatrixXd& C, double Tbase) {
    const Eigen::Index n = C.rows();
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(std::max<Eigen::Index>(n, 1), C.cols());
    Eigen::MatrixXd L = Eigen::MatrixXd::Zero(1, C.cols());
    for (Eigen::Index i = 0; i < n; ++i) {
        double minus_Tbase_power = 1.0;  // (-Tbase)^m, built up step by step
        for (Eigen::Index m = 0; m < i; ++m) {
            P.row(i - m) += C.row(i) * (minus_Tbase_power / static_cast<double>(i - m));
            minus_Tbase_power *= -Tbase;
        }
        L.row(0) += C.row(i) * minus_Tbase_power;  // now (-Tbase)^i
    }
    return {P, L};
}

// ---- Chebyshev caloric fits (type == INCOMPRESSIBLE_CHEBYSHEV) -------------
//
// The 2D layout mirrors the polynomial fits: rows are Chebyshev-in-T
// coefficients on [cheb_Tmin, cheb_Tmax], columns multiply (x - cheb_xbase)^j.
// Evaluation collapses the composition direction to an effective 1D
// coefficient vector (a linear combination of the columns), then runs
// Clenshaw -- the same collapse-then-evaluate structure Polynomial2DFrac uses.

/// Effective 1D Chebyshev coefficients at composition x.
Eigen::VectorXd chebCollapse(const Eigen::MatrixXd& C, double x, double xbase) {
    Eigen::VectorXd dx(C.cols());
    double power = 1.0;
    for (Eigen::Index j = 0; j < C.cols(); ++j) {
        dx(j) = power;
        power *= (x - xbase);
    }
    return C * dx;
}

/// Clenshaw evaluation at the mapped argument u in [-1, 1]
/// (mirrors superancillary::ChebyshevExpansion::eval).
double chebClenshaw(const Eigen::VectorXd& a, double u) {
    const int Norder = static_cast<int>(a.size()) - 1;
    if (Norder < 1) return (Norder == 0) ? a(0) : 0.0;
    double u_k = 0, u_kp1 = a(Norder), u_kp2 = 0;
    for (int k = Norder - 1; k > 0; k--) {
        u_k = 2.0 * u * u_kp1 - u_kp2 + a(k);
        u_kp2 = u_kp1;
        u_kp1 = u_k;
    }
    return a(0) + u * u_kp1 - u_kp2;
}

/// Evaluate a coefficient matrix that shares data's domain and xbase
/// (data.coeffs itself, or one of the derived matrices).
double chebEvalMatrix(const Eigen::MatrixXd& M, const CoolProp::IncompressibleData& data, double T, double x) {
    if (M.rows() < 1 || M.cols() < 1) {
        // An empty matrix would silently evaluate to 0.0 (Clenshaw over no
        // coefficients) -- rho == 0 / a missing cp integral must be loud.
        // Derived matrices are only built by IncompressibleFluid::validate().
        throw CoolProp::ValueError("Chebyshev coefficient matrix is empty; the fluid was not loaded through validate().");
    }
    const double u = (2.0 * T - (data.cheb_Tmax + data.cheb_Tmin)) / (data.cheb_Tmax - data.cheb_Tmin);
    return chebClenshaw(chebCollapse(M, x, data.cheb_xbase), u);
}

/// Coefficients of d/dT, one column at a time. Standard descending
/// recurrence d_{k} = d_{k+2} + 2(k+1) a_{k+1}, d_0 halved, times du/dT.
Eigen::MatrixXd chebDerivativeMatrix(const Eigen::MatrixXd& C, double Tmin, double Tmax) {
    const Eigen::Index n = C.rows() - 1;  // degree
    const double scale = 2.0 / (Tmax - Tmin);
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(std::max<Eigen::Index>(n, 1), C.cols());
    for (Eigen::Index j = 0; j < C.cols(); ++j) {
        for (Eigen::Index k = n - 1; k >= 0; --k) {
            const double higher = (k + 2 <= n - 1) ? D(k + 2, j) : 0.0;
            D(k, j) = higher + 2.0 * static_cast<double>(k + 1) * C(k + 1, j);
        }
        if (n >= 1) D(0, j) *= 0.5;
        D.col(j) *= scale;
    }
    return D;
}

/// Coefficients of the indefinite integral over T, one column at a time.
/// Standard antiderivative recurrence (the inverse of the derivative one),
/// times dT/du = (Tmax - Tmin)/2. The b_0 coefficient is left at zero: the
/// integration constant always cancels against the reference-state value.
Eigen::MatrixXd chebAntiderivativeMatrix(const Eigen::MatrixXd& C, double Tmin, double Tmax) {
    const Eigen::Index n = C.rows() - 1;  // degree
    const double scale = 0.5 * (Tmax - Tmin);
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(C.rows() + 1, C.cols());
    for (Eigen::Index j = 0; j < C.cols(); ++j) {
        B(1, j) = C(0, j);
        if (n >= 1) B(2, j) = 0.25 * C(1, j);
        for (Eigen::Index k = 2; k <= n; ++k) {
            B(k + 1, j) += C(k, j) / (2.0 * static_cast<double>(k + 1));
            B(k - 1, j) -= C(k, j) / (2.0 * static_cast<double>(k - 1));
        }
        B.col(j) *= scale;
    }
    return B;
}

/// Re-expand f(T)/T column-wise on the same domain by Chebyshev-Lobatto
/// interpolation. f/T is analytic on [Tmin, Tmax] (Tmin > 0 is enforced at
/// load), so the coefficients decay geometrically; the degree grows until
/// the relative tail is at rounding level. Division by T is linear, so
/// treating each (x - xbase)^j column independently is exact for every x.
Eigen::MatrixXd chebDivideByTMatrix(const Eigen::MatrixXd& C, double Tmin, double Tmax) {
    const double tail_target = 1e-13;
    Eigen::Index degree = C.rows() + 8;
    const Eigen::Index max_degree = 96;
    while (true) {
        const auto LU = CoolProp::superancillary::detail::get_LU_matrices(static_cast<std::size_t>(degree));
        const Eigen::MatrixXd& L = std::get<0>(LU);
        // Chebyshev-Lobatto nodes mapped to [Tmin, Tmax]
        Eigen::VectorXd T_nodes(degree + 1);
        for (Eigen::Index i = 0; i <= degree; ++i) {
            const double u = std::cos(EIGEN_PI * static_cast<double>(i) / static_cast<double>(degree));
            T_nodes(i) = 0.5 * ((Tmax - Tmin) * u + (Tmax + Tmin));
        }
        Eigen::MatrixXd out(degree + 1, C.cols());
        double worst_tail = 0.0;
        for (Eigen::Index j = 0; j < C.cols(); ++j) {
            Eigen::VectorXd samples(degree + 1);
            for (Eigen::Index i = 0; i <= degree; ++i) {
                const double u = std::cos(EIGEN_PI * static_cast<double>(i) / static_cast<double>(degree));
                samples(i) = chebClenshaw(C.col(j), u) / T_nodes(i);
            }
            out.col(j) = L * samples;
            const double scale = out.col(j).cwiseAbs().maxCoeff();
            if (scale > 0.0) {
                worst_tail = std::max(worst_tail, out.col(j).tail(3).cwiseAbs().maxCoeff() / scale);
            }
        }
        if (worst_tail < tail_target || degree >= max_degree) {
            return out;
        }
        degree += 16;
    }
}
}  // namespace

namespace CoolProp {

/// A thermophysical property provider for all properties
/**
This fluid instance is populated using an entry from a JSON file and uses
simplified polynomial and exponential functions to calculate thermophysical
and transport properties.
*/
//IncompressibleFluid::IncompressibleFluid();

void IncompressibleFluid::validate() {
    // Caloric fits: check them and build the derived coefficient matrices
    // exactly once, at load time. Everything h and s need afterwards is a
    // plain evaluation of a stored matrix (see drhodTatPx, dhdTatPxdT and
    // dsdTatPxdT). Both the classic and the Chebyshev fits are prepared, so
    // that either can be selected per state without reloading anything.
    struct CaloricFit
    {
        IncompressibleData* data;
        const char* what;
        bool is_specific_heat;
    };
    for (const CaloricFit& fit : {CaloricFit{&density, "density", false}, CaloricFit{&density_cheb, "density", false},
                                  CaloricFit{&specific_heat, "specific heat", true}, CaloricFit{&specific_heat_cheb, "specific heat", true}}) {
        IncompressibleData* data = fit.data;
        if (data->type == IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL) {
            if (data->coeffs.rows() < 1 || data->coeffs.cols() < 1) {
                throw ValueError(
                  format("%s (%d): Polynomial %s fit for %s has an empty coefficient matrix.", __FILE__, __LINE__, fit.what, name.c_str()));
            }
            if (fit.is_specific_heat) {
                data->int_dT = polyAntiderivativeMatrix(data->coeffs);
                std::tie(data->int_dTdivT, data->int_dTdivT_log) = polyDivideByTAntiderivative(data->coeffs, Tbase);
            } else {
                data->ddT = polyDerivativeMatrix(data->coeffs);
            }
            continue;
        }
        if (data->type != IncompressibleData::INCOMPRESSIBLE_CHEBYSHEV) continue;
        if (!(data->cheb_Tmin > 0.0) || !(data->cheb_Tmax > data->cheb_Tmin)) {
            throw ValueError(format("%s (%d): Chebyshev %s fit for %s needs 0 < Tmin < Tmax, got [%g, %g].", __FILE__, __LINE__, fit.what,
                                    name.c_str(), data->cheb_Tmin, data->cheb_Tmax));
        }
        if (data->coeffs.rows() < 1 || data->coeffs.cols() < 1) {
            throw ValueError(format("%s (%d): Chebyshev %s fit for %s has an empty coefficient matrix.", __FILE__, __LINE__, fit.what, name.c_str()));
        }
        // The fit domain must cover the fluid's advertised range. The scaled
        // variable u = (2T - (cheb_Tmax + cheb_Tmin)) / (cheb_Tmax - cheb_Tmin)
        // is not clamped, and Trange is deliberately independent of the
        // fluid-level Tmin/Tmax (it is the data range), so a fit narrower than
        // the fluid range would let checkT() admit a temperature and then have
        // Clenshaw extrapolate outside [-1, 1], where a Chebyshev series
        // diverges fast and silently. A fit WIDER than the fluid range is fine
        // and is shipped deliberately (ZS25/ZS40 cp cover ~2 K below Tmin,
        // which is real measured data).
        //
        // Compared with a tolerance, not exactly: Trange is emitted as a
        // full-precision double from the Python fitter, so an endpoint that is
        // conceptually equal to Tmax can differ in the last bits. HFE2 ships
        // cheb_Tmax = 337.41999999999996 against Tmax = 337.42, and a bit-exact
        // compare would reject it -- taking the whole incompressible library
        // down at static-init time over 4e-14 K.
        const double cover_tol = 1e-6 * std::max(1.0, Tmax - Tmin);
        if (data->cheb_Tmin > Tmin + cover_tol || data->cheb_Tmax < Tmax - cover_tol) {
            throw ValueError(format("%s (%d): Chebyshev fit domain [%g, %g] for %s does not cover the fluid range [%g, %g].", __FILE__, __LINE__,
                                    data->cheb_Tmin, data->cheb_Tmax, name.c_str(), Tmin, Tmax));
        }
        if (fit.is_specific_heat) {
            data->int_dT = chebAntiderivativeMatrix(data->coeffs, data->cheb_Tmin, data->cheb_Tmax);
            data->int_dTdivT =
              chebAntiderivativeMatrix(chebDivideByTMatrix(data->coeffs, data->cheb_Tmin, data->cheb_Tmax), data->cheb_Tmin, data->cheb_Tmax);
        } else {
            data->ddT = chebDerivativeMatrix(data->coeffs, data->cheb_Tmin, data->cheb_Tmax);
        }
    }
    // The pure/solution classification comes from the number of density
    // columns (is_pure). If the classic and the Chebyshev density disagreed,
    // the fluid would change classification with the configuration, while
    // the library files it under one name list for good.
    if (density.type != IncompressibleData::INCOMPRESSIBLE_NOT_SET && density_cheb.type != IncompressibleData::INCOMPRESSIBLE_NOT_SET
        && (density.coeffs.cols() == 1) != (density_cheb.coeffs.cols() == 1)) {
        throw ValueError(format("%s (%d): The classic and Chebyshev density fits for %s disagree on whether the fluid is pure (%d vs %d columns).",
                                __FILE__, __LINE__, name.c_str(), static_cast<int>(density.coeffs.cols()),
                                static_cast<int>(density_cheb.coeffs.cols())));
    }
}

bool IncompressibleFluid::is_pure() {
    // validate() ensures both density fits agree, so either one will do
    if (densityFit(true).coeffs.cols() == 1) return true;
    return false;
}

bool IncompressibleFluid::chebyshevCaloricPreferred() {
    return get_config_bool(INCOMP_PREFER_CHEBYSHEV_CALORIC);
}

/// Base exponential function: exp(coeffs[0] / ((y - ybase) + coeffs[1]) - coeffs[2])
double IncompressibleFluid::baseExponential(const IncompressibleData& data, double y, double ybase) {
    Eigen::VectorXd coeffs = makeVector(data.coeffs);
    size_t r = coeffs.rows(), c = coeffs.cols();
    if (strict && (r != 3 || c != 1)) {
        throw ValueError(format("%s (%d): You have to provide a 3,1 matrix of coefficients, not  (%d,%d).", __FILE__, __LINE__, r, c));
    }
    const double x_den = (y - ybase) + coeffs[1];
    return evaluateAwayFromPole(x_den, [&](double x) { return exp(coeffs[0] / x - coeffs[2]); });
}

/// Base exponential function with logarithmic term:
/// exp(log(1/x + 1/x^2) * coeffs[1] + coeffs[2]) with x = (y - ybase) + coeffs[0]
double IncompressibleFluid::baseLogexponential(const IncompressibleData& data, double y, double ybase) {
    Eigen::VectorXd coeffs = makeVector(data.coeffs);
    size_t r = coeffs.rows(), c = coeffs.cols();
    if (strict && (r != 3 || c != 1)) {
        throw ValueError(format("%s (%d): You have to provide a 3,1 matrix of coefficients, not  (%d,%d).", __FILE__, __LINE__, r, c));
    }
    const double x_den = (y - ybase) + coeffs[0];
    return evaluateAwayFromPole(x_den, [&](double x) { return exp(log(1.0 / x + 1.0 / x / x) * coeffs[1] + coeffs[2]); });
}

double IncompressibleFluid::basePolyOffset(IncompressibleData data, double y, double z) {
    size_t r = data.coeffs.rows(), c = data.coeffs.cols();
    double offset = 0.0;
    double in = 0.0;
    Eigen::MatrixXd coeffs;
    if (r > 0 && c > 0) {
        offset = data.coeffs(0, 0);
        if (r == 1 && c > 1) {  // row vector -> function of z
            coeffs = Eigen::MatrixXd(data.coeffs.block(0, 1, r, c - 1));
            in = z;
        } else if (r > 1 && c == 1) {  // column vector -> function of y
            coeffs = Eigen::MatrixXd(data.coeffs.block(1, 0, r - 1, c));
            in = y;
        } else {
            throw ValueError(format("%s (%d): You have to provide a vector (1D matrix) of coefficients, not  (%d,%d).", __FILE__, __LINE__, r, c));
        }
        return poly.evaluate(coeffs, in, 0, offset);
    }
    throw ValueError(format("%s (%d): You have to provide a vector (1D matrix) of coefficients, not  (%d,%d).", __FILE__, __LINE__, r, c));
}

/// Build the error for a property this fluid carries no fit for.
/** The message this replaces named neither the fluid nor the property.  It
 *  read "IncompressibleFluid.cpp (388): The function type is not specified
 *  ("[0]"), are you sure the coefficients have been set?", which gave a user
 *  asking for something the reference data simply does not cover no way to
 *  tell that from a bug in CoolProp, and no way to see which fluid or which
 *  property was at fault when several were in play.  See issue #2567, where
 *  T_freeze on the ice slurries surfaced exactly this way.
 *
 *  A missing fit is a gap in the data, not a bad input, so the message says
 *  so plainly and leaves out the file and line: those point into CoolProp's
 *  own source, which is not where the answer is.
 *
 *  Returned rather than thrown so that every call site keeps its own throw
 *  statement.  The compiler then still sees each switch case as terminating;
 *  a helper that threw internally would need [[noreturn]] for that, and would
 *  produce "control reaches end of non-void function" warnings on any
 *  compiler that ignored the attribute.
 */
ValueError IncompressibleFluid::notDefinedError(const std::string& property) const {
    return ValueError(format("INCOMP::%s does not define %s: no fit for this property is shipped in its fluid "
                             "data, so there is nothing to evaluate. This is a gap in the reference data for "
                             "this fluid, not a problem with your inputs.",
                             this->name.c_str(), property.c_str()));
}

/// Density as a function of temperature, pressure and composition.
double IncompressibleFluid::rho(double T, double p, double x, bool chebyshev) {
    const IncompressibleData& fit = densityFit(chebyshev);
    switch (fit.type) {
        case IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL:
            return poly.evaluate(fit.coeffs, T, x, 0, 0, Tbase, xbase);
        case IncompressibleData::INCOMPRESSIBLE_CHEBYSHEV:
            return chebEvalMatrix(fit.coeffs, fit, T, x);
        case IncompressibleData::INCOMPRESSIBLE_EXPONENTIAL:
            return baseExponential(fit, T, 0.0);
        case IncompressibleData::INCOMPRESSIBLE_LOGEXPONENTIAL:
            return baseLogexponential(fit, T, 0.0);
        case IncompressibleData::INCOMPRESSIBLE_EXPPOLYNOMIAL:
            return exp(poly.evaluate(fit.coeffs, T, x, 0, 0, Tbase, xbase));
        case IncompressibleData::INCOMPRESSIBLE_POLYOFFSET:
            return basePolyOffset(fit, T, x);
        case IncompressibleData::INCOMPRESSIBLE_NOT_SET:
            throw notDefinedError("density");
        default:
            throw ValueError(format("%s (%d): Your function type \"[%d]\" is unknown.", __FILE__, __LINE__, fit.type));
    }
}

/// Heat capacities as a function of temperature, pressure and composition.
double IncompressibleFluid::c(double T, double p, double x, bool chebyshev) {
    const IncompressibleData& fit = specificHeatFit(chebyshev);
    switch (fit.type) {
        case IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL:
            return poly.evaluate(fit.coeffs, T, x, 0, 0, Tbase, xbase);
        case IncompressibleData::INCOMPRESSIBLE_CHEBYSHEV:
            return chebEvalMatrix(fit.coeffs, fit, T, x);
        case IncompressibleData::INCOMPRESSIBLE_NOT_SET:
            throw notDefinedError("a specific heat");
        default:
            throw ValueError(format("%s (%d): There is no predefined way to use this function type \"[%d]\" for specific heat.", __FILE__, __LINE__,
                                    fit.type));
    }
}

/// Viscosity as a function of temperature, pressure and composition.
double IncompressibleFluid::visc(double T, double p, double x) {
    switch (viscosity.type) {
        case IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL:
            return poly.evaluate(viscosity.coeffs, T, x, 0, 0, Tbase, xbase);
        case IncompressibleData::INCOMPRESSIBLE_EXPONENTIAL:
            return baseExponential(viscosity, T, 0.0);
        case IncompressibleData::INCOMPRESSIBLE_LOGEXPONENTIAL:
            return baseLogexponential(viscosity, T, 0.0);
        case IncompressibleData::INCOMPRESSIBLE_EXPPOLYNOMIAL:
            return exp(poly.evaluate(viscosity.coeffs, T, x, 0, 0, Tbase, xbase));
        case IncompressibleData::INCOMPRESSIBLE_POLYOFFSET:
            return basePolyOffset(viscosity, T, x);
        case IncompressibleData::INCOMPRESSIBLE_NOT_SET:
            throw notDefinedError("a viscosity");
        default:
            throw ValueError(format("%s (%d): Your function type \"[%d]\" is unknown.", __FILE__, __LINE__, viscosity.type));
    }
}
/// Thermal conductivity as a function of temperature, pressure and composition.
double IncompressibleFluid::cond(double T, double p, double x) {
    switch (conductivity.type) {
        case IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL:
            return poly.evaluate(conductivity.coeffs, T, x, 0, 0, Tbase, xbase);
        case IncompressibleData::INCOMPRESSIBLE_EXPONENTIAL:
            return baseExponential(conductivity, T, 0.0);
        case IncompressibleData::INCOMPRESSIBLE_LOGEXPONENTIAL:
            return baseLogexponential(conductivity, T, 0.0);
        case IncompressibleData::INCOMPRESSIBLE_EXPPOLYNOMIAL:
            return exp(poly.evaluate(conductivity.coeffs, T, x, 0, 0, Tbase, xbase));
        case IncompressibleData::INCOMPRESSIBLE_POLYOFFSET:
            return basePolyOffset(conductivity, T, x);
        case IncompressibleData::INCOMPRESSIBLE_NOT_SET:
            throw notDefinedError("a thermal conductivity");
        default:
            throw ValueError(format("%s (%d): Your function type \"[%d]\" is unknown.", __FILE__, __LINE__, conductivity.type));
    }
}
/// Saturation pressure as a function of temperature and composition.
double IncompressibleFluid::psat(double T, double x) {
    // Below TminPsat the polynomial fit is not valid. Returning 0.0
    // silently here surfaced as "Psat = 0 for the entire valid T range"
    // for fluids whose TminPsat is at or above Tmax (e.g. MEG, where
    // TminPsat = Tmax = 373.15 K leaves no valid range): #2209. Throw a
    // controlled exception instead so callers see the limitation rather
    // than a misleading zero.
    if (T <= this->TminPsat) {
        throw ValueError(format("Saturation pressure is not available below TminPsat=%g K (T=%g K)", this->TminPsat, T));
    }
    switch (p_sat.type) {
        case IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL:
            return poly.evaluate(p_sat.coeffs, T, x, 0, 0, Tbase, xbase);
        case IncompressibleData::INCOMPRESSIBLE_EXPONENTIAL:
            return baseExponential(p_sat, T, 0.0);
        case IncompressibleData::INCOMPRESSIBLE_LOGEXPONENTIAL:
            return baseLogexponential(p_sat, T, 0.0);
        case IncompressibleData::INCOMPRESSIBLE_EXPPOLYNOMIAL:
            return exp(poly.evaluate(p_sat.coeffs, T, x, 0, 0, Tbase, xbase));
        case IncompressibleData::INCOMPRESSIBLE_POLYOFFSET:
            return basePolyOffset(p_sat, T, x);
        case IncompressibleData::INCOMPRESSIBLE_NOT_SET:
            throw notDefinedError("a saturation pressure");
        default:
            throw ValueError(format("%s (%d): Your function type \"[%d]\" is unknown.", __FILE__, __LINE__, p_sat.type));
    }
}
/// Freezing temperature as a function of pressure and composition.
double IncompressibleFluid::Tfreeze(double p, double x) {
    switch (T_freeze.type) {
        case IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL:
            return poly.evaluate(T_freeze.coeffs, p, x, 0, 0, 0.0, xbase);
        case IncompressibleData::INCOMPRESSIBLE_EXPONENTIAL:
            return baseExponential(T_freeze, x, 0.0);
        case IncompressibleData::INCOMPRESSIBLE_LOGEXPONENTIAL:
            return baseLogexponential(T_freeze, x, 0.0);
        case IncompressibleData::INCOMPRESSIBLE_EXPPOLYNOMIAL:
            return exp(poly.evaluate(T_freeze.coeffs, p, x, 0, 0, 0.0, xbase));
        case IncompressibleData::INCOMPRESSIBLE_POLYOFFSET:
            return basePolyOffset(T_freeze, p, x);
        case IncompressibleData::INCOMPRESSIBLE_NOT_SET:
            throw notDefinedError("a freezing temperature");
        default:
            throw ValueError(format("%s (%d): Your function type \"[%d]\" is unknown.", __FILE__, __LINE__, T_freeze.type));
    }
}

/* Below are direct calculations of the derivatives. Nothing
 * special is going on, we simply use the polynomial class to
 * derive the different functions with respect to temperature.
 */
/// Partial derivative of density with respect to temperature at constant pressure and composition
double IncompressibleFluid::drhodTatPx(double T, double p, double x, bool chebyshev) {
    const IncompressibleData& fit = densityFit(chebyshev);
    switch (fit.type) {
        case IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL:
            // d/dT, derived once at load time in validate()
            return polyEvalMatrix(fit.ddT, T, x, Tbase, xbase);
        case IncompressibleData::INCOMPRESSIBLE_CHEBYSHEV:
            return chebEvalMatrix(fit.ddT, fit, T, x);
        case IncompressibleData::INCOMPRESSIBLE_NOT_SET:
            throw notDefinedError("density");
        default:
            throw ValueError(
              format("%s (%d): There is no predefined way to use this function type \"[%d]\" for density.", __FILE__, __LINE__, fit.type));
    }
}
/// Partial derivative of entropy
//  with respect to temperature at constant pressure and composition
//  integrated in temperature
double IncompressibleFluid::dsdTatPxdT(double T, double p, double x, bool chebyshev) {
    const IncompressibleData& fit = specificHeatFit(chebyshev);
    switch (fit.type) {
        case IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL:
            // Indefinite Int c/T dT, derived once at load time in validate():
            // a polynomial part plus a ln(T) term, see polyDivideByTAntiderivative
            return polyEvalMatrix(fit.int_dTdivT, T, x, Tbase, xbase) + std::log(T) * polyEvalMatrix(fit.int_dTdivT_log, T, x, Tbase, xbase);
        case IncompressibleData::INCOMPRESSIBLE_CHEBYSHEV:
            // Indefinite Int c/T dT, derived once at load time in validate()
            return chebEvalMatrix(fit.int_dTdivT, fit, T, x);
        case IncompressibleData::INCOMPRESSIBLE_NOT_SET:
            throw notDefinedError("a specific heat");
        default:
            throw ValueError(
              format("%s (%d): There is no predefined way to use this function type \"[%d]\" for entropy.", __FILE__, __LINE__, fit.type));
    }
}
/// Partial derivative of enthalpy
//  with respect to temperature at constant pressure and composition
//  integrated in temperature
double IncompressibleFluid::dhdTatPxdT(double T, double p, double x, bool chebyshev) {
    const IncompressibleData& fit = specificHeatFit(chebyshev);
    switch (fit.type) {
        case IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL:
            // Indefinite Int c dT, derived once at load time in validate()
            return polyEvalMatrix(fit.int_dT, T, x, Tbase, xbase);
        case IncompressibleData::INCOMPRESSIBLE_CHEBYSHEV:
            // Indefinite Int c dT, derived once at load time in validate()
            return chebEvalMatrix(fit.int_dT, fit, T, x);
        case IncompressibleData::INCOMPRESSIBLE_NOT_SET:
            throw notDefinedError("a specific heat");
        default:
            throw ValueError(
              format("%s (%d): There is no predefined way to use this function type \"[%d]\" for entropy.", __FILE__, __LINE__, fit.type));
    }
}

/// Mass fraction conversion function
/** If the fluid type is mass-based, it does not do anything. Otherwise,
 *  it converts the mass fraction to the required input. */
double IncompressibleFluid::inputFromMass(double T, double x) {
    if (this->xid == IFRAC_PURE) {
        return _HUGE;
    } else if (this->xid == IFRAC_MASS) {
        return x;
    } else {
        throw NotImplementedError("Mass composition conversion has not been implemented.");
        //switch (mass2input.type) {
        //    case IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL:
        //        return poly.evaluate(mass2input.coeffs, T, x, 0, 0, 0.0, 0.0); // TODO: make sure Tbase and xbase are defined in the correct way
        //        break;
        //    case IncompressibleData::INCOMPRESSIBLE_EXPONENTIAL:
        //        return baseExponential(mass2input, x, 0.0);
        //        break;
        //    case IncompressibleData::INCOMPRESSIBLE_LOGEXPONENTIAL:
        //        return baseLogexponential(mass2input, x, 0.0);
        //        break;
        //    case IncompressibleData::INCOMPRESSIBLE_EXPPOLYNOMIAL:
        //        return exp(poly.evaluate(mass2input.coeffs, T, x, 0, 0, 0.0, 0.0)); // TODO: make sure Tbase and xbase are defined in the correct way
        //        break;
        //    case IncompressibleData::INCOMPRESSIBLE_POLYOFFSET:
        //        return basePolyOffset(mass2input, T, x);
        //        break;
        //    case IncompressibleData::INCOMPRESSIBLE_NOT_SET:
        //        throw ValueError(format("%s (%d): The function type is not specified (\"[%d]\"), are you sure the coefficients have been set?",__FILE__,__LINE__,mass2input.type));
        //        break;
        //    default:
        //        throw ValueError(format("%s (%d): Your function type \"[%d]\" is unknown.",__FILE__,__LINE__,mass2input.type));
        //        break;
        //}
        //return _HUGE;
    }
}

/// Volume fraction conversion function
/** If the fluid type is volume-based, it does not do anything. Otherwise,
 *  it converts the volume fraction to the required input. */
double IncompressibleFluid::inputFromVolume(double T, double x) {
    if (this->xid == IFRAC_PURE) {
        return _HUGE;
    } else if (this->xid == IFRAC_VOLUME) {
        return x;
    } else {
        throw NotImplementedError("Volume composition conversion has not been implemented.");
        //switch (volume2input.type) {
        //    case IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL:
        //        return poly.evaluate(volume2input.coeffs, T, x, 0, 0, 0.0, 0.0); // TODO: make sure Tbase and xbase are defined in the correct way
        //        break;
        //    case IncompressibleData::INCOMPRESSIBLE_EXPONENTIAL:
        //        return baseExponential(volume2input, x, 0.0);
        //        break;
        //    case IncompressibleData::INCOMPRESSIBLE_LOGEXPONENTIAL:
        //        return baseLogexponential(volume2input, x, 0.0);
        //        break;
        //    case IncompressibleData::INCOMPRESSIBLE_EXPPOLYNOMIAL:
        //        return exp(poly.evaluate(volume2input.coeffs, T, x, 0, 0, 0.0, 0.0)); // TODO: make sure Tbase and xbase are defined in the correct way
        //        break;
        //    case IncompressibleData::INCOMPRESSIBLE_POLYOFFSET:
        //        return basePolyOffset(volume2input, T, x);
        //        break;
        //    case IncompressibleData::INCOMPRESSIBLE_NOT_SET:
        //        throw ValueError(format("%s (%d): The function type is not specified (\"[%d]\"), are you sure the coefficients have been set?",__FILE__,__LINE__,volume2input.type));
        //        break;
        //    default:
        //        throw ValueError(format("%s (%d): Your function type \"[%d]\" is unknown.",__FILE__,__LINE__,volume2input.type));
        //        break;
        //}
        //return _HUGE;
    }
}

/// Mole fraction conversion function
/** If the fluid type is mole-based, it does not do anything. Otherwise,
 *  it converts the mole fraction to the required input. */
double IncompressibleFluid::inputFromMole(double T, double x) {
    if (this->xid == IFRAC_PURE) {
        return _HUGE;
    } else if (this->xid == IFRAC_MOLE) {
        return x;
    } else {
        throw NotImplementedError("Mole composition conversion has not been implemented.");
        /*
        switch (mole2input.type) {
            case IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL:
                return poly.evaluate(mole2input.coeffs, T, x, 0, 0, 0.0, 0.0); // TODO: make sure Tbase and xbase are defined in the correct way
                break;
            case IncompressibleData::INCOMPRESSIBLE_EXPONENTIAL:
                return baseExponential(mole2input, x, 0.0);
                break;
            case IncompressibleData::INCOMPRESSIBLE_LOGEXPONENTIAL:
                return baseLogexponential(mole2input, x, 0.0);
                break;
            case IncompressibleData::INCOMPRESSIBLE_EXPPOLYNOMIAL:
                return exp(poly.evaluate(mole2input.coeffs, T, x, 0, 0, 0.0, 0.0)); // TODO: make sure Tbase and xbase are defined in the correct way
                break;
            case IncompressibleData::INCOMPRESSIBLE_POLYOFFSET:
                return basePolyOffset(mole2input, T, x);
                break;
            case IncompressibleData::INCOMPRESSIBLE_NOT_SET:
                throw notDefinedError("a mole fraction conversion");
            default:
                throw ValueError(format("%s (%d): Your function type \"[%d]\" is unknown.",__FILE__,__LINE__,mole2input.type));
                break;
        }
        return _HUGE;
        */
    }
}

/* Some functions can be inverted directly, those are listed
 * here. It is also possible to solve for other quantities, but
 * that involves some more sophisticated processing and is not
 * done here, but in the backend, T(h,p) for example.
 */
/// Temperature as a function of density, pressure and composition.
double IncompressibleFluid::T_rho(double Dmass, double p, double x, bool chebyshev) {
    const IncompressibleData& fit = densityFit(chebyshev);
    double d_raw = Dmass;  // No changes needed, no reference values...
    switch (fit.type) {
        case IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL:
            return poly.solve_limits(fit.coeffs, x, d_raw, Tmin, Tmax, 0, 0, 0, Tbase, xbase);
        case IncompressibleData::INCOMPRESSIBLE_CHEBYSHEV: {
            // Brent on the forward evaluation over the fluid's range,
            // mirroring the backend's HmassP/PSmass flashes; the Chebyshev
            // basis has no (and needs no) direct polynomial root-solve.
            class Residual : public FuncWrapper1D
            {
               public:
                IncompressibleFluid* fluid;
                double p, x, target;
                bool chebyshev;
                Residual(IncompressibleFluid* fluid, double p, double x, double target, bool chebyshev)
                  : fluid(fluid), p(p), x(x), target(target), chebyshev(chebyshev) {}
                double call(double T) override {
                    return fluid->rho(T, p, x, chebyshev) - target;
                }
            } res(this, p, x, d_raw, chebyshev);
            return Brent(&res, Tmin, Tmax, DBL_EPSILON, 1e3 * DBL_EPSILON, 100);
        }
        case IncompressibleData::INCOMPRESSIBLE_NOT_SET:
            throw notDefinedError("density");
        default:
            // The density fit's type, not the specific heat's: T_rho
            // switches on the density fit, so reporting the heat capacity's
            // type code here named the wrong property.
            throw ValueError(format("%s (%d): There is no predefined way to use this function type \"[%d]\" for inverse density.", __FILE__, __LINE__,
                                    fit.type));
    }
}
/// Temperature as a function of heat capacities as a function of temperature, pressure and composition.
double IncompressibleFluid::T_c(double Cmass, double p, double x, bool chebyshev) {
    const IncompressibleData& fit = specificHeatFit(chebyshev);
    double c_raw = Cmass;  // No changes needed, no reference values...
    switch (fit.type) {
        case IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL:
            return poly.solve_limits(fit.coeffs, x, c_raw, Tmin, Tmax, 0, 0, 0, Tbase, xbase);
        case IncompressibleData::INCOMPRESSIBLE_CHEBYSHEV: {
            // See the matching comment in T_rho above.
            class Residual : public FuncWrapper1D
            {
               public:
                IncompressibleFluid* fluid;
                double p, x, target;
                bool chebyshev;
                Residual(IncompressibleFluid* fluid, double p, double x, double target, bool chebyshev)
                  : fluid(fluid), p(p), x(x), target(target), chebyshev(chebyshev) {}
                double call(double T) override {
                    return fluid->c(T, p, x, chebyshev) - target;
                }
            } res(this, p, x, c_raw, chebyshev);
            return Brent(&res, Tmin, Tmax, DBL_EPSILON, 1e3 * DBL_EPSILON, 100);
        }
        case IncompressibleData::INCOMPRESSIBLE_NOT_SET:
            throw notDefinedError("a specific heat");
        default:
            throw ValueError(format("%s (%d): There is no predefined way to use this function type \"[%d]\" for inverse specific heat.", __FILE__,
                                    __LINE__, fit.type));
    }
}

/*
 * Some more functions to provide a single implementation
 * of important routines.
 * We start with the check functions that can validate input
 * in terms of pressure p, temperature T and composition x.
 */
/// Check validity of temperature input.
/** Compares the given temperature T to the result of a
 *  freezing point calculation. This is not necessarily
 *  defined for all fluids, default values do not cause errors. */
bool IncompressibleFluid::checkT(double T, double p, double x) {
    if (Tmin <= 0.) throw ValueError("Please specify the minimum temperature.");
    if (Tmax <= 0.) throw ValueError("Please specify the maximum temperature.");
    if ((Tmin > T) || (T > Tmax)) throw ValueError(format("Your temperature %f is not between %f and %f.", T, Tmin, Tmax));
    double TF = 0.0;
    if (T_freeze.type != IncompressibleData::INCOMPRESSIBLE_NOT_SET) TF = Tfreeze(p, x);
    if (T < TF) throw ValueError(format("Your temperature %f is below the freezing point of %f.", T, TF));
    return true;
}

/// Check validity of pressure input.
/** Compares the given pressure p to the saturation pressure at
 *  temperature T and throws and exception if p is lower than
 *  the saturation conditions.
 *  The default value for psat is -1 yielding true if psat
 *  is not redefined in the subclass.
 *  */
bool IncompressibleFluid::checkP(double T, double p, double x) {
    double ps = 0.0;
    if (p_sat.type != IncompressibleData::INCOMPRESSIBLE_NOT_SET) {
        // psat() now throws below TminPsat (#2209). For this validation
        // path, T below TminPsat means the saturation curve is not
        // available — there is nothing to compare p against, so skip
        // the p >= psat check rather than propagating the throw to
        // callers that didn't ask for psat in the first place
        // (e.g. PropsSI("D","P",101325,"T",300,"INCOMP::DowQ") where
        // DowQ has TminPsat > 300 K).
        try {
            ps = psat(T, x);
        } catch (const ValueError&) {
            ps = 0.0;
        }
    }
    if (p < 0.0) throw ValueError(format("You cannot use negative pressures: %f < %f. ", p, 0.0));
    if (ps > 0.0 && p < ps) throw ValueError(format("Equations are valid for liquid phase only: %f < %f (psat). ", p, ps));
    return true;
}

/// Check validity of composition input.
/** Compares the given composition x to a stored minimum and
 *  maximum value. Enforces the redefinition of xmin and
 *  xmax since the default values cause an error. */
bool IncompressibleFluid::checkX(double x) {
    // is_in_closed_range, not (x < lo || x > hi): the latter is false for NaN.
    if (!is_in_closed_range(0.0, 1.0, xmin)) throw ValueError("Please specify the minimum concentration between 0 and 1.");
    if (!is_in_closed_range(0.0, 1.0, xmax)) throw ValueError("Please specify the maximum concentration between 0 and 1.");
    // is_in_closed_range swaps its bounds, so inverted data would otherwise be
    // accepted as the swapped range.
    if (!(xmin <= xmax)) throw ValueError(format("The minimum concentration %g exceeds the maximum concentration %g.", xmin, xmax));
    if (!is_in_closed_range(xmin * (1 - INCOMP_EPSILON), xmax * (1 + INCOMP_EPSILON), x)) {
        throw ValueError(format("Your composition %g is not between %g and %g.", x, xmin, xmax));
    }
    return true;
}

} /* namespace CoolProp */

// Testing still needs to be enhanced.
/* Below, I try to carry out some basic tests for both 2D and 1D
 * polynomials as well as the exponential functions for vapour
 * pressure etc.
 */
#ifdef ENABLE_CATCH
#    include <math.h>
#    include <iostream>
#    include <catch2/catch_all.hpp>
#    include "Tests/TestObjects.h"

Eigen::MatrixXd makeMatrix(const std::vector<double>& coefficients) {
    //IncompressibleClass::checkCoefficients(coefficients,18);
    std::vector<std::vector<double>> matrix;
    std::vector<double> tmpVector;

    tmpVector.clear();
    tmpVector.push_back(coefficients[0]);
    tmpVector.push_back(coefficients[6]);
    tmpVector.push_back(coefficients[11]);
    tmpVector.push_back(coefficients[15]);
    matrix.push_back(tmpVector);

    tmpVector.clear();
    tmpVector.push_back(coefficients[1] * 100.0);
    tmpVector.push_back(coefficients[7] * 100.0);
    tmpVector.push_back(coefficients[12] * 100.0);
    tmpVector.push_back(coefficients[16] * 100.0);
    matrix.push_back(tmpVector);

    tmpVector.clear();
    tmpVector.push_back(coefficients[2] * 100.0 * 100.0);
    tmpVector.push_back(coefficients[8] * 100.0 * 100.0);
    tmpVector.push_back(coefficients[13] * 100.0 * 100.0);
    tmpVector.push_back(coefficients[17] * 100.0 * 100.0);
    matrix.push_back(tmpVector);

    tmpVector.clear();
    tmpVector.push_back(coefficients[3] * 100.0 * 100.0 * 100.0);
    tmpVector.push_back(coefficients[9] * 100.0 * 100.0 * 100.0);
    tmpVector.push_back(coefficients[14] * 100.0 * 100.0 * 100.0);
    tmpVector.push_back(0.0);
    matrix.push_back(tmpVector);

    tmpVector.clear();
    tmpVector.push_back(coefficients[4] * 100.0 * 100.0 * 100.0 * 100.0);
    tmpVector.push_back(coefficients[10] * 100.0 * 100.0 * 100.0 * 100.0);
    tmpVector.push_back(0.0);
    tmpVector.push_back(0.0);
    matrix.push_back(tmpVector);

    tmpVector.clear();
    tmpVector.push_back(coefficients[5] * 100.0 * 100.0 * 100.0 * 100.0 * 100.0);
    tmpVector.push_back(0.0);
    tmpVector.push_back(0.0);
    tmpVector.push_back(0.0);
    matrix.push_back(tmpVector);

    tmpVector.clear();
    return CoolProp::vec_to_eigen(matrix).transpose();
}

TEST_CASE("Internal consistency checks and example use cases for the incompressible fluids", "[IncompressibleFluids]") {
    SECTION("Test case for \"SylthermXLT\" by Dow Chemicals") {

        std::vector<double> cRho;
        cRho.push_back(+1.1563685145E+03);
        cRho.push_back(-1.0269048032E+00);
        cRho.push_back(-9.3506079577E-07);
        cRho.push_back(+1.0368116627E-09);
        CoolProp::IncompressibleData density;
        density.type = CoolProp::IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL;
        density.coeffs = CoolProp::vec_to_eigen(cRho);

        std::vector<double> cHeat;
        cHeat.push_back(+1.1562261074E+03);
        cHeat.push_back(+2.0994549103E+00);
        cHeat.push_back(+7.7175381057E-07);
        cHeat.push_back(-3.7008444051E-20);
        CoolProp::IncompressibleData specific_heat;
        specific_heat.type = CoolProp::IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL;
        specific_heat.coeffs = CoolProp::vec_to_eigen(cHeat);

        std::vector<double> cCond;
        cCond.push_back(+1.6121957379E-01);
        cCond.push_back(-1.3023781944E-04);
        cCond.push_back(-1.4395238766E-07);
        CoolProp::IncompressibleData conductivity;
        conductivity.type = CoolProp::IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL;
        conductivity.coeffs = CoolProp::vec_to_eigen(cCond);

        std::vector<double> cVisc;
        cVisc.push_back(+1.0337654989E+03);
        cVisc.push_back(-4.3322764383E+01);
        cVisc.push_back(+1.0715062356E+01);
        CoolProp::IncompressibleData viscosity;
        viscosity.type = CoolProp::IncompressibleData::INCOMPRESSIBLE_EXPONENTIAL;
        viscosity.coeffs = CoolProp::vec_to_eigen(cVisc);

        CoolProp::IncompressibleFluid XLT;
        XLT.setName("XLT");
        XLT.setDescription("SylthermXLT");
        XLT.setReference("Dow Chemicals data sheet");
        XLT.setTmax(533.15);
        XLT.setTmin(173.15);
        XLT.setxmax(0.0);
        XLT.setxmin(0.0);
        XLT.setTminPsat(533.15);

        XLT.setTbase(0.0);
        XLT.setxbase(0.0);

        /// Setters for the coefficients
        XLT.setDensity(density);
        XLT.setSpecificHeat(specific_heat);
        XLT.setViscosity(viscosity);
        XLT.setConductivity(conductivity);
        //XLT.setPsat(parse_coefficients(fluid_json, "saturation_pressure", false));
        //XLT.setTfreeze(parse_coefficients(fluid_json, "T_freeze", false));
        //XLT.setVolToMass(parse_coefficients(fluid_json, "volume2mass", false));
        //XLT.setMassToMole(parse_coefficients(fluid_json, "mass2mole", false));

        /// A function to check coefficients and equation types.
        //XLT.validate();
        double acc = 0.0001;
        double val = 0;
        double res = 0;

        // Prepare the results and compare them to the calculated values
        double T = 273.15 + 50;
        double p = 10e5;
        double x = 0.0;

        // Compare density
        val = 824.4615702148608;
        res = XLT.rho(T, p, x);
        {
            CAPTURE(T);
            CAPTURE(val);
            CAPTURE(res);
            CHECK(check_abs(val, res, acc));
        }

        // Compare cp
        val = 1834.7455527670554;
        res = XLT.c(T, p, x);
        {
            CAPTURE(T);
            CAPTURE(val);
            CAPTURE(res);
            CHECK(check_abs(val, res, acc));
        }

        // Check property functions
        CHECK_THROWS(XLT.s(T, p, x));
        CHECK_THROWS(XLT.h(T, p, x));
        CHECK_THROWS(XLT.u(T, p, x));

        // Compare v
        val = 0.0008931435169681835;
        res = XLT.visc(T, p, x);
        {
            CAPTURE(T);
            CAPTURE(val);
            CAPTURE(res);
            CHECK(check_abs(val, res, acc));
        }

        // Compare l
        val = 0.10410086156049088;
        res = XLT.cond(T, p, x);
        {
            CAPTURE(T);
            CAPTURE(val);
            CAPTURE(res);
            CHECK(check_abs(val, res, acc));
        }
    }

    SECTION("Test case for Methanol from SecCool") {

        CoolProp::IncompressibleFluid CH3OH = CoolProp::get_incompressible_fluid("MMA2");

        // Prepare the results and compare them to the calculated values
        double acc = 0.0001;
        double T = 273.15 + 10;
        double p = 10e5;
        double x = 0.25;
        double expected = 0;
        double actual = 0;

        // Compare density
        expected = 963.2886528091547;
        actual = CH3OH.rho(T, p, x);
        {
            CAPTURE(T);
            CAPTURE(p);
            CAPTURE(x);
            CAPTURE(expected);
            CAPTURE(actual);
            CHECK(check_abs(expected, actual, acc));
        }

        // Compare cp
        expected = 3993.9748117022423;
        actual = CH3OH.c(T, p, x);
        {
            CAPTURE(T);
            CAPTURE(p);
            CAPTURE(x);
            CAPTURE(expected);
            CAPTURE(actual);
            CHECK(check_abs(expected, actual, acc));
        }

        // Check property functions
        CHECK_THROWS(CH3OH.s(T, p, x));
        CHECK_THROWS(CH3OH.h(T, p, x));
        CHECK_THROWS(CH3OH.u(T, p, x));

        // Compare viscosity
        expected = 0.0023970245009602097;  // Pa-s   // Melinder 3.48 mPas @ 24.87 wt%, 0C
        actual = CH3OH.visc(T, p, x);
        {
            CAPTURE(T);
            CAPTURE(p);
            CAPTURE(x);
            CAPTURE(expected);
            CAPTURE(actual);
            std::string errmsg = CoolProp::get_global_param_string("errstring");
            CAPTURE(errmsg);
            CHECK(check_abs(expected, actual, acc));
        }

        // Compare conductivity
        expected = 0.44791148414693727;
        actual = CH3OH.cond(T, p, x);
        {
            CAPTURE(T);
            CAPTURE(p);
            CAPTURE(x);
            CAPTURE(expected);
            CAPTURE(actual);
            std::string errmsg = CoolProp::get_global_param_string("errstring");
            CAPTURE(errmsg);
            CHECK(check_abs(expected, actual, acc));
        }

        // Compare Tfreeze
        expected = -20.02 + 273.15;  // 253.1293105454671;
        actual = CH3OH.Tfreeze(p, x);
        {
            CAPTURE(T);
            CAPTURE(p);
            CAPTURE(x);
            CAPTURE(expected);
            CAPTURE(actual);
            std::string errmsg = CoolProp::get_global_param_string("errstring");
            CAPTURE(errmsg);
            CHECK(check_abs(expected, actual, acc));
        }
    }
}

#endif /* ENABLE_CATCH */
