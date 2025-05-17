/*
# NIST Disclaimer of Copyright and Warranty

This software was developed by employees of the National Institute of
Standards and Technology (NIST), an agency of the Federal Government
and is being made available as a public service. Pursuant to title 17
United States Code Section 105, works of NIST employees are not
subject to copyright protection in the United States. This software
may be subject to foreign copyright. Permission in the United States
and in foreign countries, to the extent that NIST may hold copyright,
to use, copy, modify, create derivative works, and distribute this
software and its documentation without fee is hereby granted on a
non-exclusive basis, provided that this notice and disclaimer of
warranty appears in all copies.

THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND,
EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED
TO, ANY WARRANTY THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS,
ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE, AND FREEDOM FROM INFRINGEMENT, AND ANY WARRANTY THAT THE
DOCUMENTATION WILL CONFORM TO THE SOFTWARE, OR ANY WARRANTY THAT THE
SOFTWARE WILL BE ERROR FREE. IN NO EVENT SHALL NIST BE LIABLE FOR ANY
DAMAGES, INCLUDING, BUT NOT LIMITED TO, DIRECT, INDIRECT, SPECIAL OR
CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR IN ANY WAY
CONNECTED WITH THIS SOFTWARE, WHETHER OR NOT BASED UPON WARRANTY,
CONTRACT, TORT, OR OTHERWISE, WHETHER OR NOT INJURY WAS SUSTAINED BY
PERSONS OR PROPERTY OR OTHERWISE, AND WHETHER OR NOT LOSS WAS
SUSTAINED FROM, OR AROSE OUT OF THE RESULTS OF, OR USE OF, THE
SOFTWARE OR SERVICES PROVIDED HEREUNDER.

Subsequent edits by Ian Bell
*/

#pragma once

#include <iostream>
#include <utility>
#include <optional>
#include <Eigen/Dense>

#include "boost/math/tools/toms748_solve.hpp"
#include "nlohmann/json.hpp"

namespace CoolProp{
namespace superancillary{

namespace detail{

// From https://arxiv.org/pdf/1401.5766.pdf (Algorithm #3)
template<typename Matrix>
inline void balance_matrix(const Matrix &A, Matrix &Aprime, Matrix &D) {
    const int p = 2;
    double beta = 2; // Radix base (2?)
    int iter = 0;
    Aprime = A;
    D = Matrix::Identity(A.rows(), A.cols());
    bool converged = false;
    do {
        converged = true;
        for (Eigen::Index i = 0; i < A.rows(); ++i) {
            double c = Aprime.col(i).template lpNorm<p>();
            double r = Aprime.row(i).template lpNorm<p>();
            double s = pow(c, p) + pow(r, p);
            double f = 1;
            //if (!ValidNumber(c)){
            //    std::cout << A << std::endl;
            //    throw std::range_error("c is not a valid number in balance_matrix"); }
            //if (!ValidNumber(r)) { throw std::range_error("r is not a valid number in balance_matrix"); }
            while (c < r/beta) {
                c *= beta;
                r /= beta;
                f *= beta;
            }
            while (c >= r*beta) {
                c /= beta;
                r *= beta;
                f /= beta;
            }
            if (pow(c, p) + pow(r, p) < 0.95*s) {
                converged = false;
                D(i, i) *= f;
                Aprime.col(i) *= f;
                Aprime.row(i) /= f;
            }
        }
        iter++;
        if (iter > 50) {
            break;
        }
    } while (!converged);
}

inline void companion_matrix_transposed(const Eigen::ArrayXd &coeffs, Eigen::MatrixXd &A) {
    auto N = coeffs.size() - 1; // degree
    if (A.rows() != N){ throw std::invalid_argument("A.rows() != N"); }
    A.setZero();
    // First col
    A(1, 0) = 1;
    // Last col
    A.col(N-1) = -coeffs.head(N)/(2.0*coeffs(N));
    A(N - 2, N - 1) += 0.5;
    // All the other cols
    for (int j = 1; j < N - 1; ++j) {
        A(j - 1, j) = 0.5;
        A(j + 1, j) = 0.5;
    }
}
inline void companion_matrix_transposed(const std::vector<double> &coeffs, Eigen::MatrixXd &A) {
    Eigen::ArrayXd coeffs_ = Eigen::Map<const Eigen::ArrayXd>(&coeffs[0], coeffs.size());
    return companion_matrix_transposed(coeffs_, A);
}


/**
* @brief Get the L and U matrices needed for transformations between nodes and function values in a Chebyshev expansion
*
* The L matrix is used to convert from functional values to coefficients, as in \f[ \vec{c} = \mathbf{L}\vec{f} \f]
* The U matrix is used to convert from coefficients to functional values, as in \f[ \vec{f} = \mathbf{U}\vec{c} \f]
* \param N the degree of the expansion (one less than the number of coefficients)
*/

inline auto get_LU_matrices(std::size_t N){
    Eigen::MatrixXd L(N + 1, N + 1); ///< Matrix of coefficients
    Eigen::MatrixXd U(N + 1, N + 1); ///< Matrix of coefficients
    for (std::size_t j = 0; j <= N; ++j) {
        for (std::size_t k = j; k <= N; ++k) {
            double p_j = (j == 0 || j == N) ? 2 : 1;
            double p_k = (k == 0 || k == N) ? 2 : 1;
            double cosjPikN = cos((j*EIGEN_PI*k) / N);
            L(j, k) = 2.0/(p_j*p_k*N)*cosjPikN;
            // Exploit symmetry to fill in the symmetric elements in the matrix
            L(k, j) = L(j, k);
            
            U(j, k) = cosjPikN;
            // Exploit symmetry to fill in the symmetric elements in the matrix
            U(k, j) = U(j, k);
        }
    }
    return std::make_tuple(L, U);
}

inline double M_element_norm(const std::vector<double>& x, Eigen::Index M){
    Eigen::Map<const Eigen::ArrayXd> X(&x[0], x.size());
    return X.tail(M).matrix().norm() / X.head(M).matrix().norm();
}

inline double M_element_norm(const Eigen::ArrayXd& x, Eigen::Index M){
    return x.tail(M).matrix().norm() / x.head(M).matrix().norm();
} 

/**
 std::function<double(double)>
 */

template<typename Function, typename Container>
inline auto dyadic_splitting(const std::size_t N, const Function& func, const double xmin, const double xmax,
    const int M=3, const double tol=1e-12, const int max_refine_passes = 8,
    const std::optional<std::function<void(int, const Container&)>>& callback = std::nullopt) -> Container
{
    using CE_t = std::decay_t<decltype(Container().front())>;
    using ArrayType = std::decay_t<decltype(CE_t().coeff())>;
    
    Eigen::MatrixXd Lmat, Umat;
    std::tie(Lmat, Umat) = detail::get_LU_matrices(N);
    
    auto builder = [&](double xmin, double xmax) -> CE_t{
        
        auto get_nodes_realworld = [N, xmin, xmax]() -> Eigen::ArrayXd{
            Eigen::ArrayXd nodes  = (Eigen::ArrayXd::LinSpaced(N + 1, 0, static_cast<double>(N)).array()*EIGEN_PI / N).cos();
            return ((xmax - xmin)*nodes + (xmax + xmin))*0.5;
        };
        
        Eigen::ArrayXd x = get_nodes_realworld();

        // Apply the function to do the transformation
        // of the functional values at the nodes
        for (auto j = 0L; j < x.size(); ++j){
            x(j) = func(j, static_cast<long>(x.size()), x(j));
        }
        
        // And now rebuild the expansion by left-multiplying by the L matrix
        Eigen::ArrayXd c = Lmat*x.matrix();
//        std::cout << "c: " << c << std::endl;
        // Check if any coefficients are invalid, stop if so
        if (!c.allFinite() ) {
            throw std::invalid_argument("At least one coefficient is non-finite");
        }
        if constexpr (std::is_same_v<ArrayType, std::vector<double>>){
            return {xmin, xmax, std::vector<double>(&c[0], &c[0] + c.size())};
        }
        else{
            // Probably an Eigen::ArrayXd, just pass that into constructor
            return {xmin, xmax, c};
        }
    };
    
    // Start off with the full domain from xmin to xmax
    Container expansions;
    expansions.emplace_back(builder(xmin, xmax));

    // Now enter into refinement passes
    for (int refine_pass = 0; refine_pass < max_refine_passes; ++refine_pass) {
        bool all_converged = true;
        // Start at the right and move left because insertions will make the length increase
        for (int iexpansion = static_cast<int>(expansions.size())-1; iexpansion >= 0; --iexpansion) {
            auto& expan = expansions[iexpansion];
            auto err = M_element_norm(expan.coeff(), M);
//            std::cout << "err: " <<  err << std::endl;
            if (err > tol) {
                // Splitting is required, do a dyadic split
                auto xmid = (expan.xmin() + expan.xmax()) / 2;
                CE_t newleft{builder(expan.xmin(), xmid)};
                CE_t newright{builder(xmid, expan.xmax())};
                
                expansions.at(iexpansion) = std::move(newleft);
                expansions.insert(expansions.begin() + iexpansion+1, newright);
                all_converged = false;
            }
//            std::cout << expansions.size() << std::endl;
        }
        if (callback) {
            callback.value()(refine_pass, expansions);
        }
        if (all_converged) { break; }
    }
    return expansions;
}


}

/**
 A Chebyshev expansion of the form
 \f[
 y = \sum_i c_iT_i(x)
 \f]
 where c are the coefficients and \f$T_i\f$ are the Chebyshev basis functions of the first kind
 
 More advanced tools are possible in the ChebTools package in C++, but the essentials are available here
 
 */
template<typename ArrayType>
class ChebyshevExpansion{
private:
    double m_xmin, ///< The minimum value of the independent variable
           m_xmax; ///< The maximum value of the independent variable
    ArrayType m_coeff; ///< The coefficients of the expansion
public:
    ChebyshevExpansion() {};

    /// Constructor with bounds and coefficients
    ChebyshevExpansion(double xmin, double xmax, const ArrayType& coeff) : m_xmin(xmin), m_xmax(xmax), m_coeff(coeff){};
    
    /// Copy constructor
    ChebyshevExpansion(const ChebyshevExpansion& ex) = default;
    ChebyshevExpansion& operator=(ChebyshevExpansion&& ex) = default;
    ChebyshevExpansion& operator=(const ChebyshevExpansion& ex) = default;
    
    /// Get the minimum value of the independent variable
    const auto xmin() const { return m_xmin; }
    
    /// Get the maximum value of the independent variable
    const auto xmax() const { return m_xmax; }
    
    /// Get a const view on the expansion coefficients
    const auto& coeff() const { return m_coeff; }

    /** Evaluate the expansion with Clenshaw's method
     \param x The value of the independent variable
     */
    template<typename T>
    auto eval(const T& x) const{
        // Scale to (-1, 1)
        T xscaled = (2.0*x - (m_xmax + m_xmin)) / (m_xmax - m_xmin);
        int Norder = static_cast<int>(m_coeff.size()) - 1;
        T u_k = 0, u_kp1 = m_coeff[Norder], u_kp2 = 0;
        for (int k = Norder-1; k > 0; k--){ // k must be signed!
            // Do the recurrent calculation
            u_k = 2.0*xscaled*u_kp1 - u_kp2 + m_coeff[k];
            // Update the values
            u_kp2 = u_kp1; u_kp1 = u_k;
        }
        T retval = m_coeff[0] + xscaled*u_kp1 - u_kp2; // This seems inefficient but required to ensure that Eigen::Array are evaluated and an expression is not returned
        return retval;
    }
    
    /// A vectorized variant (for use with Python interface)
    template<typename T>
    auto eval_many(const T& x, T& y) const{
        if (x.size() != y.size()){ throw std::invalid_argument("x and y are not the same size"); }
        for (auto i = 0; i < x.size(); ++i){ y(i) = eval(x(i)); }
    }

    /// A vectorized variant in which arrays are C-style, assumed to be of the same length
    template<typename T>
    auto eval_manyC(const T x[], T y[], std::size_t N) const{
        for (std::size_t i = 0; i < N; ++i){ y[i] = eval(x[i]); }
    }
    
    /// A vectorized variant (for use with Python interface)
    template<typename T>
    auto eval_Eigen(const T& x, T& y) const{
        if (x.size() != y.size()){ throw std::invalid_argument("x and y are not the same size"); }
        y = eval(x);
    }
    
    /// Chebyshev-Lobatto nodes \f$\cos(\pi j/N), j = 0,..., N \f$ mapped to the range [xmin, xmax]
    Eigen::ArrayXd get_nodes_realworld() const {
        Eigen::Index N = m_coeff.size()-1;
        Eigen::ArrayXd nodes  = (Eigen::ArrayXd::LinSpaced(N + 1, 0, N).array()*EIGEN_PI / N).cos();
        return ((m_xmax - m_xmin)*nodes + (m_xmax + m_xmin))*0.5;
    }
    
    /**
    Solve for independent variable given a bracketing interval [a,b] with the use of the TOMS 748
    algorithm from boost which is an improvement over Brent's method (in all cases, asymptotically,
    acccording to the boost docs)
 
    \param y target value to be matched
    \param a left bound for interval
    \param b right bound for interval
    \param bits number of bits to be matched in TOMS748 solver
    \param max_iter maximum nimber of function calls
    \param boundsytol tolerance that is considered to be the right solution
    \return Tuple of value of x and the number of function evaluations required
     */
    auto solve_for_x_count(double y, double a, double b, unsigned int bits, std::size_t max_iter, double boundsytol) const{
        using namespace boost::math::tools;
        std::size_t counter = 0;
        auto f = [&](double x){ counter++; return eval(x) - y; };
        double fa = f(a);
        if (std::abs(fa) < boundsytol){ return std::make_tuple(a, std::size_t{1}); }
        double fb = f(b);
        if (std::abs(fb) < boundsytol){ return std::make_tuple(b, std::size_t{2}); }
        boost::math::uintmax_t max_iter_ = static_cast<boost::math::uintmax_t>(max_iter);
        auto [l, r] = toms748_solve(f, a, b, fa, fb, eps_tolerance<double>(bits), max_iter_);
        return std::make_tuple((r+l)/2.0, counter);
    }
    
    /// Return the value of x only for given value of y
    /// \sa solve_for_x_count
    auto solve_for_x(double y, double a, double b, unsigned int bits, std::size_t max_iter, double boundsytol) const{
        return std::get<0>(solve_for_x_count(y, a, b, bits, max_iter, boundsytol));
    }
    
    /// A vectorized variant (for use with Python interface)
    template<typename T>
    auto solve_for_x_many(const T& y, double a, double b, unsigned int bits, std::size_t max_iter, double boundsytol, T& x, T& counts) const{
        if (x.size() != y.size()){ throw std::invalid_argument("x and y are not the same size"); }
        for (auto i = 0; i < x.size(); ++i){ std::tie(x(i), counts(i)) = solve_for_x_count(y(i), a, b, bits, max_iter, boundsytol); }
    }

    /// A vectorized variant in which arrays are C-style, assumed to be of the same length
    template<typename T>
    auto solve_for_x_manyC(const T y[], std::size_t N, double a, double b, unsigned int bits, std::size_t max_iter, double boundsytol, T x[], T counts[]) const{
        for (std::size_t i = 0; i < N; ++i){
            std::tie(x[i], counts[i]) = solve_for_x_count(y[i], a, b, bits, max_iter, boundsytol);
        }
    }
    
    ArrayType do_derivs(std::size_t Nderiv) const{
        // See Mason and Handscomb, p. 34, Eq. 2.52
        // and example in https ://github.com/numpy/numpy/blob/master/numpy/polynomial/chebyshev.py#L868-L964
        ArrayType c = m_coeff;
        for (std::size_t deriv_counter = 0; deriv_counter < Nderiv; ++deriv_counter) {
            std::size_t N = c.size() - 1, ///< Degree of the expansion
                        Nd = N - 1; ///< Degree of the derivative expansion
            ArrayType cd(N);
            for (std::size_t r = 0; r <= Nd; ++r) {
                cd[r] = 0;
                for (std::size_t k = r + 1; k <= N; ++k) {
                    // Terms where k-r is odd have values, otherwise, they are zero
                    if ((k - r) % 2 == 1) {
                        cd[r] += 2*k*c[k];
                    }
                }
                // The first term with r = 0 is divided by 2 (the single prime in Mason and Handscomb, p. 34, Eq. 2.52)
                if (r == 0) {
                    cd[r] /= 2;
                }
                // Rescale the values if the range is not [-1,1].  Arrives from the derivative of d(xreal)/d(x_{-1,1})
                cd[r] /= (m_xmax-m_xmin)/2.0;
            }
            if (Nderiv == 1) {
                return cd;
            }
            else{
                c = cd;
            }
        }
        return c;
    }
};

static_assert(std::is_move_assignable_v<ChebyshevExpansion<std::vector<double>>>);
//static_assert(std::is_copy_assignable_v<ChebyshevExpansion<std::vector<double>>>);

/// Data associated with monotonic expansion
struct MonotonicExpansionMatch{
    std::size_t idx; ///< The index of the expansion that has been matched
    double ymin, ///< The minimum value of the dependent variable
           ymax, ///< The maximum value of the dependent variable
           xmin, ///< The minimum value of the independent variable
           xmax; ///< The maximum value of the independent variable
    /// Check if a value of the dependent variable is within this match
    bool contains_y (double y) const { return y >= ymin && y <= ymax; }
};

/// Data associated with a monotonic interval
struct IntervalMatch{
    std::vector<MonotonicExpansionMatch> expansioninfo; ///< The information about the expansions for this interval
    double xmin, ///< The minimum value of the independent variable
           xmax, ///< The maximum value of the independent variable
           ymin, ///< The minimum value of the dependent variable
           ymax; ///< The maximum value of the dependent variable
    /// Check if a value of the dependent variable is within this interval
    bool contains_y (double y) const { return y >= ymin && y <= ymax; }
};

/**
 A set of multiple Chebyshev expansions covering an interval [xmin, xmax]. This is a 1D approximation.
 Practically speaking the independent variable is temperature, but the code was left generic to highlight the generiticity
 of the approach
 
 At construction, the independent variable is subdivided into portions
 that are each monotonic in the dependent variable, to facilitate later rootfinding in domains that are therefore
 known to be monotonic, and therefore invertible
*/
template<typename ArrayType=Eigen::ArrayXd>
struct ChebyshevApproximation1D{
private:
    const double thresh_imag = 1e-15; ///< The threshold below which a complex number is considered to be imaginary
    std::vector<ChebyshevExpansion<ArrayType>> m_expansions; ///< The collection of expansions forming the approximation
    std::vector<double> m_x_at_extrema; ///< The values of the independent variable at the extrema of the expansions
    std::vector<IntervalMatch> m_monotonic_intervals; ///< The intervals that are monotonic

    Eigen::ArrayXd head(const Eigen::ArrayXd& c, Eigen::Index N) const{
        return c.head(N);
    }
    std::vector<double> head(const std::vector<double>& c, Eigen::Index N) const{
        return std::vector<double>(c.begin(), c.begin()+N);
    }
    
    /** Determine the values of x for the extrema where y'(x)=0 according to the expansions
     \param expansions The set of expansions that are to be traversed to identify extrema
     \param thresh_im The threshold on the imaginary part of an eigenvalue rootfinding solution to deem it to be "real enough"
    */
    std::vector<double> determine_extrema(const std::decay_t<decltype(m_expansions)>& expansions, double thresh_im) const{
        std::vector<double> x_at_extrema;
        Eigen::MatrixXd companion_matrix, cprime, D;
        for (auto& expan : expansions){
            // Get the coefficients of the derivative's expansion (for x in [-1, 1])
            auto cd = expan.do_derivs(1);
            // First, right-trim the coefficients that are equal to zero
            int ilastnonzero = static_cast<int>(cd.size())-1;
            for (int i = static_cast<int>(cd.size())-1; i >= 0; --i){
                if (std::abs(cd[i]) != 0){
                    ilastnonzero = i; break;
                }
            }
            if (ilastnonzero != static_cast<int>(cd.size()-1)){
                cd = head(cd, ilastnonzero);
            }
            // Then do eigenvalue rootfinding after balancing
            // Define working buffers here to avoid allocations all over the place
            if (companion_matrix.rows() != static_cast<int>(cd.size()-1)){
                companion_matrix.resize(cd.size()-1, cd.size()-1); companion_matrix.setZero();
                D.resizeLike(companion_matrix); D.setZero();
                cprime.resizeLike(companion_matrix); cprime.setZero();
            }
            detail::companion_matrix_transposed(cd, companion_matrix);
            cprime = companion_matrix;
            
            detail::balance_matrix(companion_matrix, cprime, D);
            for (auto &root : companion_matrix.eigenvalues()){
                // re_n11 is in [-1,1], need to rescale back to real units
                auto re_n11 = root.real(), im = root.imag();
                if (std::abs(im) < thresh_im){
                    if (re_n11 >= -1 && re_n11 <= 1){
                        x_at_extrema.push_back(((expan.xmax() - expan.xmin())*re_n11 + (expan.xmax() + expan.xmin())) / 2.0);
                    }
                }
            }
        }
        std::sort(x_at_extrema.begin(), x_at_extrema.end());
        return x_at_extrema;
    }
    
    /** Build information about intervals in which the function to be approximated is monotonic (invertible)
     \param x_at_extrema The values of x where extema exist inside the edges of the overall interval
     \note The edges of the intervals are added to the front and end of the interval before the comparisons begin
     */
    std::vector<IntervalMatch> build_monotonic_intervals(const std::vector<double>& x_at_extrema) const{
        std::vector<IntervalMatch> intervals;
        
        auto sort = [](double& x, double &y){ if (x > y){ std::swap(x, y);} };
        
        if (false){//x_at_extrema.empty()){
            auto xmin = m_expansions.front().xmin();
            auto xmax = m_expansions.back().xmax();
            auto ymin = eval(xmin), ymax = eval(xmax);
            sort(ymin, ymax);
            MonotonicExpansionMatch mem;
            mem.idx = 0;
            mem.xmin = xmin;
            mem.xmax = xmax;
            mem.ymin = ymin;
            mem.ymax = ymax;
            IntervalMatch im;
            im.expansioninfo = {mem};
            im.xmin = mem.xmin;
            im.xmax = mem.xmax;
            im.ymin = mem.ymin;
            im.ymax = mem.ymax;
            intervals.push_back(im);
        }
        else{
            auto newx = x_at_extrema;
            newx.insert(newx.begin(), m_expansions.front().xmin());
            newx.insert(newx.end(), m_expansions.back().xmax());
            
            /// See: https://stackoverflow.com/a/20023538
            auto interval_intersection = [](const auto&t1, const auto& t2){
                auto a = std::max(t1.xmin(), t2.xmin);
                auto b = std::min(t1.xmax(), t2.xmax);
                return std::make_tuple(a, b);
            };
            
            for (auto j = 0; j < static_cast<int>(newx.size()-1); ++j){
                double xmin = newx[j], xmax = newx[j+1];
                IntervalMatch im;
                // Loop over the expansions that contain one of the values of x
                // that intersect the interval defined by [xmin, xmax]
                for (auto i = 0UL; i < m_expansions.size(); ++i){
                    struct A {double xmin, xmax; };
                    auto [a,b] = interval_intersection(m_expansions[i], A{xmin, xmax} );
                    if (a < b){
                        const auto&e = m_expansions[i];
                        MonotonicExpansionMatch mem;
                        mem.idx = i;
                        mem.xmin = a;
                        mem.xmax = b;
                        double ymin = e.eval(a), ymax = e.eval(b); sort(ymin, ymax);
                        mem.ymin = ymin;
                        mem.ymax = ymax;
                        im.expansioninfo.push_back(mem);
                    }
                }
                // These are the limits for the entire interval
                im.xmin = xmin;
                im.xmax = xmax;
                double yminoverall = eval(xmin), ymaxoverall = eval(xmax); sort(yminoverall, ymaxoverall);
                im.ymin = yminoverall;
                im.ymax = ymaxoverall;
                intervals.push_back(im);
            }
        }
        return intervals;
    }
    double m_xmin, ///< The minimum value of the independent variable
           m_xmax; ///< The maximum value of the independent variable
    
public:
    
    // Move constructor given a vector of expansions
    ChebyshevApproximation1D(std::vector<ChebyshevExpansion<ArrayType>> && expansions) :
        m_expansions(std::move(expansions)),
        m_x_at_extrema(determine_extrema(m_expansions, thresh_imag)),
        m_monotonic_intervals(build_monotonic_intervals(m_x_at_extrema)),
        m_xmin(get_expansions().front().xmin()),
        m_xmax(get_expansions().back().xmax())
    {}
    
    ChebyshevApproximation1D(const ChebyshevApproximation1D & other) :
        m_expansions(other.m_expansions),
        m_x_at_extrema(other.m_x_at_extrema),
        m_monotonic_intervals(other.m_monotonic_intervals),
        m_xmin(other.xmin()),
        m_xmax(other.xmax())
    {}
    
    ChebyshevApproximation1D& operator=(ChebyshevApproximation1D && other) = default;
//    ChebyshevApproximation1D& operator=(const ChebyshevApproximation1D& other) = default;
    ChebyshevApproximation1D& operator=(ChebyshevApproximation1D other) {
        std::swap(m_expansions, other.m_expansions);
        std::swap(m_x_at_extrema, other.m_x_at_extrema);
        std::swap(m_monotonic_intervals, other.m_monotonic_intervals);
        std::swap(m_xmin, other.m_xmin);
        std::swap(m_xmax, other.m_xmax);
        return *this;
    }
    
    /// Get a const view on the expansions owned by the approximation instance
    const auto& get_expansions() const { return m_expansions; }
    
    /// Get a const view on values of x at the extrema
    const auto& get_x_at_extrema() const { return m_x_at_extrema; }
    
    /// Get a const view on the monotonic intervals identified
    const auto& get_monotonic_intervals() const { return m_monotonic_intervals; }
    
    /// Get the minimum x value
    const auto xmin() const { return m_xmin; }
    
    /// Get the maximum x value
    const auto xmax() const { return m_xmax; }
    
    /// Check whether the function is monotonic, if so some simplifications can be made to
    /// rootfinding in some cases
    bool is_monotonic() const {
        return m_monotonic_intervals.size() == 1;
    }
    
    /** Return the index of the expansion that is desired
     * \param x value of x
     * \returns The index of the expansion in the array of expansions that the point is within
     */
    auto get_index(double x) const {
        
        // https://proquest.safaribooksonline.com/9780321637413
        // https://web.stanford.edu/class/archive/cs/cs107/cs107.1202/lab1/
        auto midpoint_Knuth = [](int x, int y) {
            return (x & y) + ((x ^ y) >> 1);
        };
        
        int iL = 0U, iR = static_cast<int>(m_expansions.size()) - 1, iM;
        while (iR - iL > 1) {
            iM = midpoint_Knuth(iL, iR);
            if (x >= m_expansions[iM].xmin()) {
                iL = iM;
            }
            else {
                iR = iM;
            }
        }
        return (x < m_expansions[iL].xmax()) ? iL : iR;
    };
    
    /** Evaluate the value from the expansion
     \param x The value of the independent variable
     \returns y The evaluation of y(x) from the expansion
     */
    double eval(double x) const{
        return m_expansions[get_index(x)].eval(x);
    }
    
    /// A vectorized and templated getter (for calling from python)
    template<typename Container>
    const auto eval_many(const Container& x, Container& y) const {
        if (x.size() != y.size()){ throw std::invalid_argument("x and y are not the same size"); }
        for (auto i = 0U; i < x.size(); ++i){
            y(i) = eval(x(i));
        }
    }

    /// A vectorized and templated getter without any allocation or checking
    template<typename Container>
    const auto eval_manyC(const Container x[], Container y[], std::size_t N) const {
        for (auto i = 0U; i < N; ++i){
            y[i] = eval(x[i]);
        }
    }
    
    /// Find the intervals containing the value of y
    const std::vector<IntervalMatch> get_intervals_containing_y(double y) const{
        std::vector<IntervalMatch> matches;
        for (auto & interval : m_monotonic_intervals){
            if (y >= interval.ymin && y <= interval.ymax){
                matches.push_back(interval);
            }
        }
        return matches;
    }
    
    /** Solve for (possibly multiple) values of the independent variable x given a value of the dependent variable y
     */
    const auto get_x_for_y(double y, unsigned int bits, std::size_t max_iter, double boundsftol) const {
        std::vector<std::pair<double, int>> solns;
        for (const auto& interval: m_monotonic_intervals){
            // Each interval is required to be monotonic internally, so if the value of
            // y is within the y values at the endpoints it is a candidate
            // for rootfinding, otherwise continue
            if (interval.contains_y(y)){
                // Now look at the expansions that intersect the interval
                for (const auto& ei: interval.expansioninfo){
                    // Again, since the portions of the expansions are required
                    // to be monotonic, if it is contained then a solution must exist
                    if (ei.contains_y(y)){
                        const ChebyshevExpansion<ArrayType>& e = m_expansions[ei.idx];
                        auto [xvalue, num_steps] = e.solve_for_x_count(y, ei.xmin, ei.xmax, bits, max_iter, boundsftol);
                        solns.emplace_back(xvalue, static_cast<int>(num_steps));
                    }
                }
            }
        }
        return solns;
    }
    
    /// A vectorized and templated getter (for calling from python)
    template<typename Container>
    const auto count_x_for_y_many(const Container& y, unsigned int bits, std::size_t max_iter, double boundsftol, Container& x) const {
        if (x.size() != y.size()){ throw std::invalid_argument("x and y are not the same size"); }
        for (auto i = 0U; i < x.size(); ++i){
            x(i) = get_x_for_y(y(i), bits, max_iter, boundsftol).size();
        }
    }

    /// A vectorized and templated getter (for calling from python)
    template<typename Container>
    const auto count_x_for_y_manyC(const Container y[], size_t N, unsigned int bits, std::size_t max_iter, double boundsftol, Container x[]) const {
        for (auto i = 0U; i < N; ++i){
            x[i] = get_x_for_y(y[i], bits, max_iter, boundsftol).size();
        }
    }
};

static_assert(std::is_copy_constructible_v<ChebyshevApproximation1D<std::vector<double>>>);
static_assert(std::is_copy_assignable_v<ChebyshevApproximation1D<std::vector<double>>>);

struct SuperAncillaryTwoPhaseSolution{
    double T, ///< The temperature, in K
           q; ///< The vapor quality
    std::size_t counter; ///< Counter for how many steps have been taken
};

/**
 
A superancillary object is formed of a number of one dimensional Chebyshev approximations, one for each phase, property pair.
 
 Loaded from the file are density and pressure as functions of temperature, and a thermodynamic model can be used to build
 the
 
 */
template<typename ArrayType=Eigen::ArrayXd>
class SuperAncillary{
private:
    /// These ones must always be present
    ChebyshevApproximation1D<ArrayType> m_rhoL, ///< Approximation of \f$\rho'(T)\f$
                                        m_rhoV, ///< Approximation of \f$\rho''(T)\f$
                                        m_p; ///< Approximation of \f$p(T)\f$
    
    // These ones *may* be present
    std::optional<ChebyshevApproximation1D<ArrayType>> m_hL, ///< Approximation of \f$h'(T)\f$
                                        m_hV, ///< Approximation of \f$h''(T)\f$
                                        m_sL, ///< Approximation of \f$s'(T)\f$
                                        m_sV, ///< Approximation of \f$s''(T)\f$
                                        m_uL, ///< Approximation of \f$u'(T)\f$
                                        m_uV, ///< Approximation of \f$u''(T)\f$
                                        m_invlnp;///< Approximation of \f$T(ln(p))\f$
    
    double m_Tmin; ///< The minimum temperature, in K
    double m_Tcrit_num; ///< The numerical critical temperature, in K
    double m_rhocrit_num; ///< The numerical critical density, in mol/m^3
    double m_pmin; ///< The minimum pressure, in Pa
    double m_pmax; ///< The maximum pressure, in Pa
    
    

    /** A convenience function to load a ChebyshevExpansion from a JSON data structure
     \param j The JSON data
     \param key The key to be loaded from the superancillary block, probably one of "jexpansions_rhoL", "jexpansions_rhoV", or "jexpansions_p"
     */
    auto loader(const nlohmann::json &j, const std::string& key){
        std::vector<ChebyshevExpansion<ArrayType>> buf;
        // If you want to use Eigen...
        // auto toeig = [](const nlohmann::json &j) -> Eigen::ArrayXd{
        //     auto vec = j.get<std::vector<double>>();
        //     return Eigen::Map<const Eigen::ArrayXd>(&vec[0], vec.size());
        // };
        // for (auto& block : j[key]){
        //     buf.emplace_back(block.at("xmin"), block.at("xmax"), toeig(block.at("coef")));
        // }
        for (auto& block : j[key]){
            buf.emplace_back(block.at("xmin"), block.at("xmax"), block.at("coef"));
        }
        return buf;
    }
    
    /** Make an inverse ChebyshevApproximation1D for T(p)
     \param Ndegree The degree of the expansion in each interval. In double precision, 12 to 16 is a good plan if you will not be doing eigenvalue rootfinding. If you will be doing eigenvalue rootfinding, use a lower degree expansion, 8 is good.
     */
    auto make_invlnp(Eigen::Index Ndegree){
        
        auto pmin = m_p.eval(m_p.xmin());
        auto pmax = m_p.eval(m_p.xmax());
        // auto N = m_p.get_expansions().front().coeff().size()-1;
        const double EPSILON = std::numeric_limits<double>::epsilon();
        
        auto func = [&](long i, long Npts, double lnp) -> double{
            double p = exp(lnp);
            auto solns = m_p.get_x_for_y(p, 64, 100U, 1e-8);
            
            if (solns.size() != 1){
                if ((i == 0 || i == Npts-1) && ( p > pmin*(1-EPSILON*1000) && p < pmin*(1+EPSILON*1000))){
                    return m_p.get_monotonic_intervals()[0].xmin;
                }
                if ((i == 0 || i == Npts-1) && ( p > pmax*(1-EPSILON*1000) && p < pmax*(1+EPSILON*1000))){
                    return m_p.get_monotonic_intervals()[0].xmax;
                }
                std::stringstream ss;
                ss << std::setprecision(20) << "Must have one solution for T(p), got " << solns.size() << " for " << p << " Pa; limits are [" << pmin << + " Pa , " << pmax << " Pa]; i is " << i;
                throw std::invalid_argument(ss.str());
            }
            auto [T, iters] = solns[0];
            return T;
        };
        
        using CE_t = std::vector<ChebyshevExpansion<ArrayType>>;
        return detail::dyadic_splitting<decltype(func), CE_t>(Ndegree, func, log(pmin), log(pmax), 3, 1e-12, 26);
    }
    
    // using PropertyPairs = properties::PropertyPairs; ///< A convenience alias to save some typing
    
public:
    /// Reading in a data structure in the JSON format of https://pubs.aip.org/aip/jpr/article/53/1/013102/3270194
    /// which includes sets of Chebyshev expansions for rhoL, rhoV, and p
    SuperAncillary(const nlohmann::json &j) :
    m_rhoL(std::move(loader(j, "jexpansions_rhoL"))),
    m_rhoV(std::move(loader(j, "jexpansions_rhoV"))),
    m_p(std::move(loader(j, "jexpansions_p"))),
    m_Tmin(m_p.xmin()),
    m_Tcrit_num(j.at("meta").at("Tcrittrue / K")),
    m_rhocrit_num(j.at("meta").at("rhocrittrue / mol/m^3")),
    m_pmin(m_p.eval(m_p.xmin())),
    m_pmax(m_p.eval(m_p.xmax()))
    {};
    
    /** Load the superancillary with the data passed in as a string blob. This constructor delegates directly to the the one that consumes JSON
     * \param s The string-encoded JSON data for the superancillaries
     */
    SuperAncillary(const std::string& s) : SuperAncillary(nlohmann::json::parse(s)) {};
    
    /** Get a const reference to a ChebyshevApproximation1D
     
     \param k The key for the property (D,S,H,P,U)
     \param Q The vapor quality, either 0 or 1
     */
    const auto& get_approx1d(char k, short Q) const {
        auto get_or_throw = [&](const auto& v) -> const auto& {
            if (v){
                return v.value();
            }
            else{
                throw std::invalid_argument("unable to get the variable "+std::string(1, k)+", make sure it has been added to superancillary");
            }
        };
        switch(k){
            case 'P': return m_p;
            case 'D': return (Q == 0) ? m_rhoL : m_rhoV;
            case 'H': return (Q == 0) ? get_or_throw(m_hL) : get_or_throw(m_hV);
            case 'S': return (Q == 0) ? get_or_throw(m_sL) : get_or_throw(m_sV);
            case 'U': return (Q == 0) ? get_or_throw(m_uL) : get_or_throw(m_uV);
            default: throw std::invalid_argument("Bad key of '" + std::string(1, k) + "'");
        }
    }
    /// Get a const reference to the inverse approximation for T(ln(p))
    const auto& get_invlnp(){
        // Lazily construct on the first access
        if (!m_invlnp){
            // Degree of expansion is the same as 
            auto Ndegree = m_p.get_expansions()[0].coeff().size()-1;
            m_invlnp = make_invlnp(Ndegree);
        }
        return m_invlnp;
    }
    /// Get the minimum pressure in Pa
    const double get_pmin() const{ return m_pmin; }
    /// Get the maximum pressure in Pa
    const double get_pmax() const{ return m_pmax; }
    /// Get the minimum temperature in K
    const double get_Tmin() const{ return m_Tmin; }
    /// Get the numerical critical temperature in K
    const double get_Tcrit_num() const{ return m_Tcrit_num; }
    /// Get the numerical critical density in mol/m^3
    const double get_rhocrit_num() const{ return m_rhocrit_num; }
    
    /**
     Using the provided function that gives y(T, rho), build the ancillaries for this variable based on the ancillaries for rhoL and rhoV
     \param var The key for the property (H,S,U)
     \param caller A function that takes temperature and molar density and returns the property of interest, molar enthalpy in the case of H, etc.
     */
    void add_variable(char var, const std::function<double(double, double)> & caller){
        Eigen::MatrixXd Lmat, Umat;
        std::tie(Lmat, Umat) = detail::get_LU_matrices(12);
        
        auto builder = [&](char var, auto& variantL, auto& variantV){
            std::vector<ChebyshevExpansion<ArrayType>> newexpL, newexpV;
            const auto& expsL = get_approx1d('D', 0);
            const auto& expsV = get_approx1d('D', 1);
            if (expsL.get_expansions().size() != expsV.get_expansions().size()){
                throw std::invalid_argument("L&V are not the same size");
            }
            for (auto i = 0U; i < expsL.get_expansions().size(); ++i){
                const auto& expL = expsL.get_expansions()[i];
                const auto& expV = expsV.get_expansions()[i];
                const auto& T = expL.get_nodes_realworld();
                // Get the functional values at the Chebyshev nodes
                Eigen::ArrayXd funcL = Umat*expL.coeff().matrix();
                Eigen::ArrayXd funcV = Umat*expV.coeff().matrix();
                // Apply the function inplace to do the transformation
                // of the functional values at the nodes
                for (auto j = 0; j < funcL.size(); ++j){
                    funcL(j) = caller(T(j), funcL(j));
                    funcV(j) = caller(T(j), funcV(j));
                }
                // And now rebuild the expansions by left-multiplying by the L matrix
                newexpL.emplace_back(expL.xmin(), expL.xmax(), (Lmat*funcL.matrix()).eval());
                newexpV.emplace_back(expV.xmin(), expV.xmax(), (Lmat*funcV.matrix()).eval());
            }
            
            variantL.emplace(std::move(newexpL));
            variantV.emplace(std::move(newexpV));
        };
        
        switch(var){
            case 'H': builder(var, m_hL, m_hV); break;
            case 'S': builder(var, m_sL, m_sV); break;
            case 'U': builder(var, m_uL, m_uV); break;
            default: throw std::invalid_argument("nope");
        }
    }
    
    /** Given the value of Q in {0,1}, evaluate one of the the ChebyshevApproximation1D
     \param T Temperature, in K
     \param k Property key, in {D,P,H,S,U}
     \param Q Vapor quality, in {0,1}
     */
    double eval_sat(double T, char k, short Q) const {
        if (Q == 1 || Q == 0){
            return get_approx1d(k, Q).eval(T);
        }
        else{
            throw std::invalid_argument("invalid_value for Q");
        }
    }
    
    /**
    A vectorized version of eval_sat for wrapping in Python interface and profiling
     */
    template <typename Container>
    void eval_sat_many(const Container& T, char k, short Q, Container& y) const {
        if (T.size() != y.size()){ throw std::invalid_argument("T and y are not the same size"); }
        const auto& approx = get_approx1d(k, Q);
        for (auto i =  0; i < T.size(); ++i){
            y(i) = approx.eval(T(i));
        }
    }

    /**
    A vectorized version of eval_sat for wrapping in Python interface and profiling
     */
    template <typename Container>
    void eval_sat_manyC(const Container T[], std::size_t N, char k, short Q, Container y[]) const {
        // if (T.size() != y.size()){ throw std::invalid_argument("T and y are not the same size"); }
        const auto& approx = get_approx1d(k, Q);
        for (std::size_t i =  0; i < N; ++i){
            y[i] = approx.eval(T[i]);
        }
    }
    
    /** A convenience function to pass off to the ChebyshevApproximation1D and do an inversion calculation for a value of the variable for a saturated state
     
     \param propval The value of the property
     \param k Property key, in {D,P,H,S,U}
     \param Q Vapor quality, in {0,1}
     \param bits passed to toms748 algorithm
     \param max_iter Maximum allowed number of function calls
     \param boundsftol A functional value stopping condition to test on the endpoints
     */
    auto solve_for_T(double propval, char k, bool Q, unsigned int bits=64, unsigned int max_iter=100, double boundsftol=1e-13) const{
        return get_approx1d(k, Q).get_x_for_y(propval, bits, max_iter, boundsftol);
    }
    
    /** Get the non-iterative vapor quality q given the temperature T and the value of the thermodynamic variable
    \param T Temperature, in K
    \param propval The value of the given thermodynamic variable
    \param k Property key, in {D,P,H,S,U}
     \returns The non-iterative vapor quality based on the values from the superancillary functions
    */
    auto get_vaporquality(double T, double propval, char k) const {
        if (k == 'D'){
            // Need to special-case here because v = q*v_V + (1-q)*v_V but it is NOT(!!!!) the case that rho = q*rho_V + (1-q)*rho_L
            double v_L = 1/get_approx1d('D', 0).eval(T);
            double v_V = 1/get_approx1d('D', 1).eval(T);
            return (1/propval-v_L)/(v_V-v_L);
        }
        else{
            double L = get_approx1d(k, 0).eval(T);
            double V = get_approx1d(k, 1).eval(T);
            return (propval-L)/(V-L);
        }
    }
    
    /**
     \brief Use the inverted pressure superancillary to calculate temperature given pressure
     \param p The pressure (not its logarithm!), in Pa
     \returns T The temperature, in K
     */
    auto get_T_from_p(double p) {
        return get_invlnp().value().eval(log(p));
    }
    
    /**
     \brief Return the evaluated value of the thermodynamic variable, given the temperature and vapor quality.
     
     \param T Temperature, in K
     \param q Vapor quality, in [0,1]
     \param k Property key, in {D,P,H,S,U}
     */
    auto get_yval(double T, double q, char k) const{
            
        if (k == 'D'){
            // Need to special-case here because v = q*v_V + (1-q)*v_V but it is NOT(!!!!) the case that rho = q*rho_V + (1-q)*rho_L
            double v_L = 1/get_approx1d('D', 0).eval(T);
            double v_V = 1/get_approx1d('D', 1).eval(T);
            double v = q*v_V + (1-q)*v_L;
            return 1/v; // rho = 1/v
        }
        else{
            double L = get_approx1d(k, 0).eval(T);
            double V = get_approx1d(k, 1).eval(T);
            return q*V + (1-q)*L;
        }
    }
    
    /// A vectorized version of get_yval for profiling in Python
    template <typename Container>
    void get_yval_many(const Container& T, char k, const Container& q, Container& y) const{
        if (T.size() != y.size() || T.size() != q.size()){ throw std::invalid_argument("T, q, and y are not all the same size"); }
        
        const auto& L = get_approx1d(k, 0);
        const auto& V = get_approx1d(k, 1);
        
        if (k == 'D'){
            for (auto i = 0; i < T.size(); ++i){
                // Need to special-case here because v = q*v_V + (1-q)*v_V but it is NOT(!!!!) the case that rho = q*rho_V + (1-q)*rho_L
                double v_L = 1.0/L.eval(T(i));
                double v_V = 1.0/V.eval(T(i));
                double v = q(i)*v_V + (1-q(i))*v_L;
                y(i) = 1/v;
            }
        }
        else{
            for (auto i = 0; i < T.size(); ++i){
                y(i) = q(i)*V.eval(T(i)) + (1-q(i))*L.eval(T(i));
            }
        }
    }
    /** Determine all the values of temperature that correspond to intersections with the superancillary function, for both the vapor and liquid phases
     \param k Property key, in {D,P,H,S,U}
     \param val Value of the thermodynamic variable
     \param bits passed to toms748 algorithm
     \param max_iter Maximum allowed number of function calls
     \param boundsftol A functional value stopping condition to test on the endpoints
    */
    auto get_all_intersections(const char k, const double val, unsigned int bits, std::size_t max_iter, double boundsftol) const{
        const auto& L = get_approx1d(k, 0);
        const auto& V = get_approx1d(k, 1);
        auto TsatL = L.get_x_for_y(val, bits, max_iter, boundsftol);
        const auto TsatV = V.get_x_for_y(val, bits, max_iter, boundsftol);
        for (auto &&el : TsatV){
            TsatL.push_back(el);
        }
//        TsatL.insert(TsatL.end(),
//                     std::make_move_iterator(TsatV.begin() + TsatV.size()),
//                     std::make_move_iterator(TsatV.end()));
        return TsatL;
    }
    
    /**
    \brief Iterate to find a value of temperature and vapor quality corresponding to the two given thermodynamic variables, if such a solution exists. This is the lower-level function used by the solve_XX methods
     
     \note The temperature range must bound the solution, you might need to call get_all_intersections and parse its solutions to construct bounded intervals
     \param Tmin Minimum temperature, in K
     \param Tmax Maximum temperature, in K
     \param ch1 The key for the first variable, in {T,D,P,H,S,U}
     \param val1 The value for the first variable
     \param ch2 The key for the second variable, in {T,D,P,H,S,U}
     \param val2 The value for the second variable
     \param bits passed to toms748 algorithm
     \param max_iter Maximum allowed number of function calls
     \param boundsftol A functional value stopping condition to test on the endpoints
     */
    std::optional<SuperAncillaryTwoPhaseSolution> iterate_for_Tq_XY(double Tmin, double Tmax, char ch1, double val1, char ch2, double val2, unsigned int bits, std::size_t max_iter, double boundsftol) const {
        
        std::size_t counter = 0;
        auto f = [&](double T_){
            counter++;
            double q_fromv1 = get_vaporquality(T_, val1, ch1);
            double resid = get_yval(T_, q_fromv1, ch2) - val2;
            return resid;
        };
        using namespace boost::math::tools;
        double fTmin = f(Tmin), fTmax = f(Tmax);
        
        // First we check if we are converged enough (TOMS748 does not stop based on function value)
        double T;
        
        // A little lambda to make it easier to return
        // in different logical branches
        auto returner = [&](){
            SuperAncillaryTwoPhaseSolution soln;
            soln.T = T;
            soln.q = get_vaporquality(T, val1, ch1);
            soln.counter = counter;
            return soln;
        };
        
        if (std::abs(fTmin) < boundsftol){
            T = Tmin; return returner();
        }
        if (std::abs(fTmax) < boundsftol){
            T = Tmax; return returner();
        }
        if (fTmin*fTmax > 0){
            // No sign change, this means that the inputs are not within the two-phase region
            // and thus no iteration is needed
            return std::nullopt;
        }
        // Neither bound meets the convergence criterion, we need to iterate on temperature
        try{
            boost::math::uintmax_t max_iter_ = static_cast<boost::math::uintmax_t>(max_iter);
            auto [l, r] = toms748_solve(f, Tmin, Tmax, fTmin, fTmax, eps_tolerance<double>(bits), max_iter_);
            T = (r+l)/2.0;
            return returner();
        }
        catch(...){
            std::cout << "fTmin,fTmax: " << fTmin << "," << fTmax << std::endl;
            throw;
        }
    }
    
    /**
     Given a saturated density and another property other than T, solve for the temperature and vapor quality
     \param rho The molar density
     \param propval The value of the other property
     \param k Property key, in {D,P,H,S,U}
     \param bits passed to toms748 algorithm
     \param max_iter Maximum allowed number of function calls
     \param boundsftol A functional value stopping condition to test on the endpoints
     */
    std::optional<SuperAncillaryTwoPhaseSolution> solve_for_Tq_DX(const double rho, const double propval, const char k, unsigned int bits, std::size_t max_iter, double boundsftol) const {
        
        const auto& Lrho = get_approx1d('D', 0);
        const auto& Vrho = get_approx1d('D', 1);
        auto Tsat = get_all_intersections(k, propval, bits, max_iter, boundsftol);
        std::size_t rhosat_soln_count = Tsat.size();
        
        std::tuple<double, double> Tsearchrange;
        if (rhosat_soln_count == 1){
            // Search the complete range from Tmin to the intersection point where rhosat(T) = rho
            // obtained just above
            Tsearchrange = std::make_tuple(Lrho.xmin*0.999, std::get<0>(Tsat[0]));
        }else if (rhosat_soln_count == 2){
            double y1 = std::get<0>(Tsat[0]), y2 = std::get<0>(Tsat[1]);
            if (y2 < y1){ std::swap(y2, y2); } // Required to be sorted in increasing value
            Tsearchrange = std::make_tuple(y1, y2);
        }
        else{
            throw std::invalid_argument("cannot have number of solutions other than 1 or 2; got "+std::to_string(rhosat_soln_count)+" solutions");
        }
        auto [a, b] = Tsearchrange;
        return iterate_for_Tq_XY(a, b, 'D', rho, k, propval, bits, max_iter, boundsftol);
    }
    
    /// A vectorize version of solve_for_Tq_DX for use in the Python interface for profiling
    template <typename Container>
    void solve_for_Tq_DX_many(const Container& rho, const Container& propval, const char k, unsigned int bits, std::size_t max_iter, double boundsftol, Container& T, Container& q, Container& count){
        if (std::set<std::size_t>({rho.size(), propval.size(), T.size(), q.size(), count.size()}).size() != 1){
            throw std::invalid_argument("rho, propval, T, q, count are not all the same size");
        }
        for (auto i = 0U; i < T.size(); ++i){
            auto osoln = solve_for_Tq_DX(rho(i), propval(i), k, bits, max_iter, boundsftol);
            if (osoln){
                const auto& soln = osoln.value();
                T(i) = soln.T;
                q(i) = soln.q;
                count(i) = soln.counter;
            }
            else{
                T(i) = -1;
                q(i) = -1;
                count(i) = -1;
            }
        }
    }
    
//     /**
//     The high-level function used to carry out a solution. It handles all the different permutations of variables and delegates to lower-level functions to actually od the calculations
     
//      \param pair The enumerated pair of thermodynamic variables being provided
//      \param val1 The first value
//      \param val2 The second value
//      */
//     auto flash(PropertyPairs pair, double val1, double val2) const -> std::optional<SuperAncillaryTwoPhaseSolution>{
//         double T, q;
//         std::size_t counter = 0;
//         double boundsftol = 1e-12;
        
//         auto returner = [&](){
//             SuperAncillaryTwoPhaseSolution soln;
//             soln.T = T;
//             soln.q = q;
//             soln.counter = counter;
//             return soln;
//         };
        
//         auto get_T_other = [&](){
//             switch (pair){
//                 case PropertyPairs::DT: return std::make_tuple(val2, val1, 'D');
//                 case PropertyPairs::HT: return std::make_tuple(val2, val1, 'H');
//                 case PropertyPairs::ST: return std::make_tuple(val2, val1, 'S');
//                 case PropertyPairs::TU: return std::make_tuple(val1, val2, 'U');
//                 default:
//                     throw std::invalid_argument("not valid");
//             };
//         };
//         auto get_p_other = [&](){
//             switch (pair){
//                 case PropertyPairs::DP: return std::make_tuple(val2, val1, 'D');
//                 case PropertyPairs::HP: return std::make_tuple(val2, val1, 'H');
//                 case PropertyPairs::PS: return std::make_tuple(val1, val2, 'S');
//                 case PropertyPairs::PU: return std::make_tuple(val1, val2, 'U');
//                 default:
//                     throw std::invalid_argument("not valid");
//             };
//         };
//         auto get_rho_other = [&](){
//             switch (pair){
//                 case PropertyPairs::DH: return std::make_tuple(val1, val2, 'H');
//                 case PropertyPairs::DS: return std::make_tuple(val1, val2, 'S');
//                 case PropertyPairs::DU: return std::make_tuple(val1, val2, 'U');
//                 default:
//                     throw std::invalid_argument("not valid");
//             };
//         };
        
//         // Given the arguments, unpack them into (xchar, xval, ychar, yval) where x is the
//         // variable of convenience that will be used to do the interval intersection and then
//         // the iteration will be in the other variable
//         auto parse_HSU = [&](){
//             switch (pair){
//                 case PropertyPairs::HS: return std::make_tuple('S', val2, 'H', val1);
//                 case PropertyPairs::SU: return std::make_tuple('S', val1, 'U', val2);
//                 case PropertyPairs::HU: return std::make_tuple('H', val1, 'U', val2);
//                 default:
//                     throw std::invalid_argument("not valid");
//             };
//         };

//         switch(pair){
//             // Case 0, PT is always single-phase, by definition
//             case PropertyPairs::PT: throw std::invalid_argument("PT inputs are trivial");
                
//             // Case 1, T is a given variable
//             case PropertyPairs::DT:
//             case PropertyPairs::HT:
//             case PropertyPairs::ST:
//             case PropertyPairs::TU:
//             {
//                 auto [Tval, other, otherchar] = get_T_other();
//                 q = get_vaporquality(Tval, other, otherchar);
//                 T = Tval;
//                 if (0 <= q && q <= 1){
//                     return returner();
//                 }
//                 break;
//             }
//             // Case 2, p is a given variable, p gets converted to T, and then q calculated
//             case PropertyPairs::DP:
//             case PropertyPairs::HP:
//             case PropertyPairs::PS:
//             case PropertyPairs::PU:
//             {
//                 auto [p, other, otherchar] = get_p_other();
//                 T = get_T_from_p(p);
//                 q = get_vaporquality(T, other, otherchar); // Based on the T and the other variable
//                 if (0 <= q && q <= 1){
//                     return returner();
//                 }
//                 break;
//             }
//             // Case 3, rho is a given variable, special case that because it is relatively simple
//             // and only water and heavy water have more than one density solution along the saturation curve
//             case PropertyPairs::DH:
//             case PropertyPairs::DS:
//             case PropertyPairs::DU:
//             {
//                 auto [rho, otherval, otherchar] = get_rho_other();
//                 auto optsoln = solve_for_Tq_DX(rho, otherval, otherchar, 64, 100, boundsftol);
//                 if (optsoln){
//                     const auto& soln = optsoln.value();
//                     T = soln.T;
//                     q = soln.q;
//                     counter = soln.counter;
//                     return returner();
//                 }
//                 break;
//             }
//             // Case 4, handle all the cases where complicated interval checks are
//             // necessary
//             case PropertyPairs::HS:
//             case PropertyPairs::HU:
//             case PropertyPairs::SU:
//             {
//                 auto [xchar, xval, ychar, yval] = parse_HSU();
//                 auto Tsat = get_all_intersections(xchar, xval, 64, 100, boundsftol);
//                 auto sort_predicate = [](const auto& x, const auto& y) {
//                     return std::get<0>(x) < std::get<0>(y);
//                 };
//                 std::sort(Tsat.begin(), Tsat.end(), sort_predicate);
// //                std::cout << Tsat.size() << " T crossings for " << xval << std::endl;
// //                for (auto & el : Tsat){
// //                    std::cout << std::get<0>(el) << std::endl;
// //                }
//                 if (Tsat.size() == 1){
//                     // A single crossing, in which the other temperature bound for iteration
//                     // is assumed to be the minimum temperature of the superancillary
//                     double Tmin = get_approx1d('D', 0).xmin*0.9999;
//                     double Tmax = std::get<0>(Tsat[0]);
//                     auto o = iterate_for_Tq_XY(Tmin, Tmax, xchar, xval, ychar, yval, 64, 100, boundsftol);
//                     if (o){ return o.value(); }
//                 }
//                 else if (Tsat.size() % 2 == 0){ // Even number of crossings
//                     // neighboring intersections should bound the solution, or if not, the point is single-phase
//                     for (auto i = 0; i < Tsat.size(); i += 2){
//                         auto o = iterate_for_Tq_XY(std::get<0>(Tsat[i]), std::get<0>(Tsat[i+1]), xchar, xval, ychar, yval, 64, 100, boundsftol);
//                         if (o){ return o.value(); }
//                     }
//                 }
//                 else{
//                     throw std::invalid_argument("Don't know what to do about this number of T crossings ("+std::to_string(Tsat.size())+") yet");
//                 }
//                 break;
//             }
//             default:
//                 throw std::invalid_argument("These inputs are not yet supported in flash");
//         }
//         return std::nullopt;
//     }
    
//     /// A vectorized version of the flash function used in the Python interface for profiling
//     template <typename Container>
//     void flash_many(const PropertyPairs ppair, const Container& val1, const Container& val2, Container& T, Container& q, Container& count){
//         if (std::set<std::size_t>({val1.size(), val2.size(), T.size(), q.size(), count.size()}).size() != 1){
//             throw std::invalid_argument("val1, val2, T, q, count are not all the same size");
//         }
//         for (auto i = 0U; i < T.size(); ++i){
//             try{
//                 auto osoln = flash(ppair, val1(i), val2(i));
//                 if (osoln){
//                     const auto& soln = osoln.value();
//                     T(i) = soln.T;
//                     q(i) = soln.q;
//                     count(i) = soln.counter;
//                 }
//                 else{
//                     T(i) = -1;
//                     q(i) = -1;
//                     count(i) = -1;
//                 }
//             }
//             catch(std::exception&e){
// //                std::cout << e.what() << " for " << val1(i) << "," << val2(i) << std::endl;
//                 T(i) = -2;
//                 q(i) = -2;
//                 count(i) = -2;
//             }
            
//         }
//     }
};

static_assert(std::is_copy_constructible_v<SuperAncillary<std::vector<double>>>);
static_assert(std::is_copy_assignable_v<SuperAncillary<std::vector<double>>>);

}
}
