#ifndef COOLPROP_NUMERICS_H
#define COOLPROP_NUMERICS_H

#include <vector>
#include <set>
#include <cfloat>
#include <stdlib.h> // For abs
#include <algorithm> // For max
#include <numeric>
#include "CPstrings.h"
#include "Exceptions.h"

#if defined(HUGE_VAL) && !defined(_HUGE)
    # define _HUGE HUGE_VAL
#else
    // GCC Version of huge value macro
    #if defined(HUGE) && !defined(_HUGE)
    #  define _HUGE HUGE
    #endif
#endif

inline bool ValidNumber(double x)
{
    // Idea from http://www.johndcook.com/IEEE_exceptions_in_cpp.html
    return (x <= DBL_MAX && x >= -DBL_MAX);
};

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

#ifndef COOLPROP_OK
#define COOLPROP_OK 1
#endif

// Undefine these terrible macros defined in windows header
#undef min
#undef max

/// Make a linearly spaced vector of points
template <typename T> std::vector<T> linspace(T xmin, T xmax, std::size_t n) {
    std::vector<T> x(n, 0.0);
    
    for ( std::size_t i = 0;  i < n; ++i) {
        x[i] = (xmax-xmin)/(n-1)*i+xmin;
    }
    return x;
}
/// Make a base-10 logarithmically spaced vector of points
template <typename T> std::vector<T> log10space(T xmin, T xmax, std::size_t n) {
    std::vector<T> x(n, 0.0);
    T logxmin = log10(xmin), logxmax = log10(xmax);
    
    for ( std::size_t i = 0;  i < n; ++i) {
        x[i] = exp((logxmax-logxmin)/(n-1)*i+logxmin);
    }
    return x;
}
/// Make a base-e logarithmically spaced vector of points
template <typename T> std::vector<T> logspace(T xmin, T xmax, std::size_t n) {
    std::vector<T> x(n, 0.0);
    T logxmin = log(xmin), logxmax = log(xmax);
    
    for ( std::size_t i = 0;  i < n; ++i) {
        x[i] = exp((logxmax-logxmin)/(n-1)*i+logxmin);
    }
    return x;
}

/**
 * @brief Use bisection to find the inputs that bisect the value you want, the trick
 * here is that this function is allowed to have "holes" where parts of the the array are 
 * also filled with invalid numbers for which ValidNumber(x) is false
 * @param vec The vector to be bisected
 * @param val The value to be found
 * @param i The index to the left of the final point; i and i+1 bound the value
 */
template <typename T> void bisect_vector(const std::vector<T> &vec, T val, std::size_t &i)
{
    T rL, rM, rR;
    std::size_t N = vec.size(), L = 0, R = N-1, M = (L+R)/2;
    // Move the right limits in until they are good
    while (!ValidNumber(vec[R])){
        if (R == 1){ throw CoolProp::ValueError("All the values in bisection vector are invalid"); }
        R--;
    }
    // Move the left limits in until they are good
    while (!ValidNumber(vec[L])){
        if (L == vec.size()-1){ throw CoolProp::ValueError("All the values in bisection vector are invalid"); }
        L++;
    }
    rL = vec[L] - val; rR = vec[R] - val;
    while (R - L > 1){
        if (!ValidNumber(vec[M])){
            std::size_t MR = M, ML = M;
            // Move middle-right to the right until it is ok
            while (!ValidNumber(vec[MR])){
                if (MR == vec.size()-1){ throw CoolProp::ValueError("All the values in bisection vector are invalid"); }
                MR++;
            }
            // Move middle-left to the left until it is ok
            while (!ValidNumber(vec[ML])){
                if (ML == 1){ throw CoolProp::ValueError("All the values in bisection vector are invalid"); }
                ML--;
            }
            T rML = vec[ML] - val; 
            T rMR = vec[MR] - val;
            // Figure out which chunk is the good part
            if (rR*rML > 0 && rL*rML < 0){
                // solution is between L and ML
                R = ML; rR = vec[ML] - val;
            }
            else if (rR*rMR < 0 && rL*rMR > 0){
                // solution is between R and MR
                L = MR; rL = vec[MR] - val;
            }
            else{
                throw CoolProp::ValueError(format("Unable to bisect segmented vector; neither chunk contains the solution %g %g %g %g", rL, rML, rMR, rR));
            }
            M = (L+R)/2;
        }
        else{
            rM = vec[M] - val;
            if (rR*rM > 0 && rL*rM < 0){
                // solution is between L and M
                R = M; rR = vec[R] - val;
            }
            else{
                // solution is between R and M
                L = M; rL = vec[L] - val;
            }
            M = (L+R)/2;
        }
    }
    i = L;
}

/**
 * @brief Use bisection to find the inputs that bisect the value you want, the trick
 * here is that this function is allowed to have "holes" where parts of the the array are 
 * also filled with invalid numbers for which ValidNumber(x) is false
 * @param matrix The vector to be bisected
 * @param j The index of the matric in the off-grain dimension
 * @param val The value to be found
 * @param i The index to the left of the final point; i and i+1 bound the value
 */
template <typename T> void bisect_segmented_vector_slice(const std::vector<std::vector<T> > &mat, std::size_t j, T val, std::size_t &i)
{
    T rL, rM, rR;
    std::size_t N = mat[j].size(), L = 0, R = N-1, M = (L+R)/2;
    // Move the right limits in until they are good
    while (!ValidNumber(mat[R][j])){
        if (R == 1){ throw CoolProp::ValueError("All the values in bisection vector are invalid"); }
        R--;
    }
    rR = mat[R][j] - val;
    // Move the left limits in until they are good
    while (!ValidNumber(mat[L][j])){
        if (L == mat.size()-1){ throw CoolProp::ValueError("All the values in bisection vector are invalid"); }
        L++;
    }
    rL = mat[L][j] - val;
    while (R - L > 1){
        if (!ValidNumber(mat[M][j])){
            std::size_t MR = M, ML = M;
            // Move middle-right to the right until it is ok
            while (!ValidNumber(mat[MR][j])){
                if (MR == mat.size()-1){ throw CoolProp::ValueError("All the values in bisection vector are invalid"); }
                MR++;
            }
            // Move middle-left to the left until it is ok
            while (!ValidNumber(mat[ML][j])){
                if (ML == 1){ throw CoolProp::ValueError("All the values in bisection vector are invalid"); }
                ML--;
            }
            T rML = mat[ML][j] - val; 
            T rMR = mat[MR][j] - val;
            // Figure out which chunk is the good part
            if (rR*rMR > 0 && rL*rML < 0){
                // solution is between L and ML
                R = ML; rR = mat[ML][j] - val;
            }
            else if (rR*rMR < 0 && rL*rML > 0){
                // solution is between R and MR
                L = MR; rL = mat[MR][j] - val;
            }
            else{
                throw CoolProp::ValueError(format("Unable to bisect segmented vector slice; neither chunk contains the solution %g %g %g %g", rL,rML,rMR,rR));
            }
            M = (L+R)/2;
        }
        else{
            rM = mat[M][j] - val;
            if (rR*rM > 0 && rL*rM < 0){
                // solution is between L and M
                R = M; rR = mat[R][j] - val;
            }
            else{
                // solution is between R and M
                L = M; rL = mat[L][j] - val;
            }
            M = (L+R)/2;
        }
    }
    i = L;
}

// From http://rosettacode.org/wiki/Power_set#C.2B.2B
inline std::size_t powerset_dereference(std::set<std::size_t>::const_iterator v) { return *v; };
      
// From http://rosettacode.org/wiki/Power_set#C.2B.2B
inline std::set<std::set<std::size_t> > powerset(std::set<std::size_t> const& set)
{
  std::set<std::set<std::size_t> > result;
  std::vector<std::set<std::size_t>::const_iterator> elements;
  do
  {
    std::set<std::size_t> tmp;
    std::transform(elements.begin(), elements.end(),
                   std::inserter(tmp, tmp.end()),
                   powerset_dereference);
    result.insert(tmp);
    if (!elements.empty() && ++elements.back() == set.end())
    {
      elements.pop_back();
    }
    else
    {
      std::set<std::size_t>::const_iterator iter;
      if (elements.empty())
      {
        iter = set.begin();
      }
      else
      {
        iter = elements.back();
        ++iter;
      }
      for (; iter != set.end(); ++iter)
      {
        elements.push_back(iter);
      }
    }
  } while (!elements.empty());
 
  return result;
}

/// Some functions related to testing and comparison of values
bool inline check_abs(double A, double B, double D){
    double max = std::abs(A);
    double min = std::abs(B);
    if (min>max) {
        max = min;
        min = std::abs(A);
    }
    if (max>DBL_EPSILON*1e3) return ( ( 1.0-min/max*1e0 ) < D );
    else throw CoolProp::ValueError(format("Too small numbers: %f cannot be tested with an accepted error of %f for a machine precision of %f. ",max,D,DBL_EPSILON));
};
bool inline check_abs(double A, double B){
    return check_abs(A,B,1e5*DBL_EPSILON);
};

template<class T> void normalize_vector(std::vector<T> &x)
{
    // Sum up all the elements in the vector
    T sumx = std::accumulate( x.begin(), x.end(), static_cast<T>(0) );
    // Normalize the components by dividing each by the sum
    for (std::size_t i = 0; i < x.size(); ++i){
        x[i] /= sumx;
    }
};

/// A spline is a curve given by the form y = ax^3 + bx^2 + c*x + d
/// As there are 4 constants, 4 constraints are needed to create the spline.  These constraints could be the derivative or value at a point
/// Often, the value and derivative of the value are known at two points.
class SplineClass
{
protected:
    int Nconstraints;
    std::vector<std::vector<double> > A;
    std::vector<double> B;
public:
    double a,b,c,d;
    SplineClass():Nconstraints(0),A(4, std::vector<double>(4, 0)),B(4,0),a(_HUGE),b(_HUGE),c(_HUGE),d(_HUGE){}
    bool build(void);
    bool add_value_constraint(double x, double y);
    void add_4value_constraints(double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4);
    bool add_derivative_constraint(double x, double dydx);
    double evaluate(double x);
};

/// from http://stackoverflow.com/a/5721830/1360263
template<class T> T factorial(T n) 
{
    if (n == 0)
        return 1;
    return n * factorial(n - 1);
}
/// see https://proofwiki.org/wiki/Nth_Derivative_of_Mth_Power
/// and https://proofwiki.org/wiki/Definition:Falling_Factorial
template<class T1, class T2> T1 nth_derivative_of_x_to_m(T1 x, T2 n, T2 m)
{
    if (n > m){
        return 0;
    }
    else{
        return factorial(m)/factorial(m-n)*pow(x, m-n);
    }
}

void MatInv_2(double A[2][2] , double B[2][2]);

double root_sum_square(const std::vector<double> &x);
double interp1d(const std::vector<double> *x, const std::vector<double> *y, double x0);
double powInt(double x, int y);
    
template<class T> T POW2(T x) { return x*x; }
template<class T> T POW3(T x) { return POW2(x)*x; }
template<class T> T POW4(T x) { return POW2(x)*POW2(x); }
#define POW5(x) ((x)*(x)*(x)*(x)*(x))
#define POW6(x) ((x)*(x)*(x)*(x)*(x)*(x))
#define POW7(x) ((x)*(x)*(x)*(x)*(x)*(x)*(x))

template<class T> T LinearInterp(T x0, T x1, T y0, T y1, T x)
{
    return (y1-y0)/(x1-x0)*(x-x0)+y0;
};
template<class T1, class T2> T2 LinearInterp(const std::vector<T1> &x, const std::vector<T1> &y, std::size_t i0, std::size_t i1, T2 val)
{
    return LinearInterp(x[i0],x[i1],y[i0],y[i1], static_cast<T1>(val));
};

template<class T> T QuadInterp(T x0, T x1, T x2, T f0, T f1, T f2, T x)
{
    /* Quadratic interpolation.
    Based on method from Kreyszig,
    Advanced Engineering Mathematics, 9th Edition
    */
    T L0, L1, L2;
    L0=((x-x1)*(x-x2))/((x0-x1)*(x0-x2));
    L1=((x-x0)*(x-x2))/((x1-x0)*(x1-x2));
    L2=((x-x0)*(x-x1))/((x2-x0)*(x2-x1));
    return L0*f0+L1*f1+L2*f2;
};
template<class T1, class T2> T2 QuadInterp(const std::vector<T1> &x, const std::vector<T1> &y, std::size_t i0, std::size_t i1, std::size_t i2, T2 val)
{
    return QuadInterp(x[i0],x[i1],x[i2],y[i0],y[i1],y[i2],static_cast<T1>(val));
};

template<class T> T CubicInterp( T x0, T x1, T x2, T x3, T f0, T f1, T f2, T f3, T x)
{
    /*
    Lagrange cubic interpolation as from
    http://nd.edu/~jjwteach/441/PdfNotes/lecture6.pdf
    */
    T L0,L1,L2,L3;
    L0=((x-x1)*(x-x2)*(x-x3))/((x0-x1)*(x0-x2)*(x0-x3));
    L1=((x-x0)*(x-x2)*(x-x3))/((x1-x0)*(x1-x2)*(x1-x3));
    L2=((x-x0)*(x-x1)*(x-x3))/((x2-x0)*(x2-x1)*(x2-x3));
    L3=((x-x0)*(x-x1)*(x-x2))/((x3-x0)*(x3-x1)*(x3-x2));
    return L0*f0 + L1*f1 + L2*f2 + L3*f3;
};
/** /brief Calculate the first derivative of the function using a cubic interpolation form
    */
template<class T> T CubicInterpFirstDeriv( T x0, T x1, T x2, T x3, T f0, T f1, T f2, T f3, T x)
{
    // Based on http://math.stackexchange.com/a/809946/66405
    T L0=((x-x1)*(x-x2)*(x-x3))/((x0-x1)*(x0-x2)*(x0-x3));
    T dL0_dx = (1/(x-x1) + 1/(x-x2) + 1/(x-x3) )*L0;
    T L1=((x-x0)*(x-x2)*(x-x3))/((x1-x0)*(x1-x2)*(x1-x3));
    T dL1_dx = (1/(x-x0) + 1/(x-x2) + 1/(x-x3) )*L1;
    T L2=((x-x0)*(x-x1)*(x-x3))/((x2-x0)*(x2-x1)*(x2-x3));
    T dL2_dx = (1/(x-x0) + 1/(x-x1) + 1/(x-x3) )*L2;
    T L3=((x-x0)*(x-x1)*(x-x2))/((x3-x0)*(x3-x1)*(x3-x2));
    T dL3_dx = (1/(x-x0) + 1/(x-x1) + 1/(x-x2) )*L3;
    return dL0_dx*f0 + dL1_dx*f1 + dL2_dx*f2 + dL3_dx*f3;
};
template<class T1, class T2> T2 CubicInterp(const std::vector<T1> &x, const std::vector<T1> &y, std::size_t i0, std::size_t i1, std::size_t i2, std::size_t i3, T2 val)
{
    return CubicInterp(x[i0],x[i1],x[i2],x[i3],y[i0],y[i1],y[i2],y[i3],static_cast<T1>(val));
};

template<class T> T is_in_closed_range( T x1, T x2, T x)
{
    return (x >= std::min(x1,x2) && x <= std::max(x1,x2));
};

/** \brief Solve a cubic with coefficients in decreasing order
    * 
    * 0 = ax^3 + b*x^2 + c*x + d
    * 
    * @param a The x^3 coefficient
    * @param b The x^2 coefficient
    * @param c The x^1 coefficient
    * @param d The x^0 coefficient
    * @param N The number of unique real solutions found
    * @param x0 The first solution found
    * @param x1 The second solution found
    * @param x2 The third solution found
    */
void solve_cubic(double a, double b, double c, double d, int &N, double &x0, double &x1, double &x2);

template<class T> inline T min3(T x1, T x2, T x3){return std::min(std::min(x1, x2), x3);};
template<class T> inline T max3(T x1, T x2, T x3){return std::max(std::max(x1, x2), x3);};
template<class T> inline T min4(T x1, T x2, T x3, T x4){return std::min(std::min(std::min(x1, x2), x3), x4);};
template<class T> inline T max4(T x1, T x2, T x3, T x4){return std::max(std::max(std::max(x1, x2), x3), x4);};

inline bool double_equal(double a, double b){return std::abs(a - b) <= 1 * DBL_EPSILON * std::max(std::abs(a), std::abs(b));};

template<class T> T max_abs_value(const std::vector<T> &x)
{
    T max = 0;
    std::size_t N = x.size();
    for (std::size_t i = 0; i < N; ++i)
    {
        T axi = std::abs(x[i]);
        if (axi > max){ max = axi; }
    }
    return max;
}
    
template<class T> T min_abs_value(const std::vector<T> &x)
{
    T min = 1e40;
    std::size_t N = x.size();
    for (std::size_t i = 0; i < N; ++i)
    {
        T axi = std::abs(x[i]);
        if (axi < min){ min = axi; }
    }
    return min;
}

inline int Kronecker_delta(std::size_t i, std::size_t j){ 
    if (i == j) { 
        return static_cast<int>(1); 
    } 
    else { 
        return static_cast<int>(0);
    } 
};
inline int Kronecker_delta(int i, int j){
    if (i == j) {
        return 1;
    }
    else {
        return 0;
    }
};

/// Sort three values in place; see http://codereview.stackexchange.com/a/64763
template<typename T>
void sort3(T &a, T &b, T &c){
    if (a > b) {
        std::swap(a, b);
    }
    if (a > c) {
        std::swap(a, c);
    }
    if (b > c) {
        std::swap(b, c);
    }
}


/**
* Due to the periodicity of angles, you need to handle the case where the
* angles wrap around - suppose theta_d is 6.28 and you are at an angles of 0.1 rad,
* the difference should be around 0.1, not -6.27
* 
* This brilliant method is from http://blog.lexique-du-net.com/index.php?post/Calculate-the-real-difference-between-two-angles-keeping-the-sign
* and the comment of user tk
* 
* Originally implemented in PDSim
*/
template<class T> T angle_difference(T angle1, T angle2){
    return fmod(angle1 - angle2 + M_PI, 2*M_PI) - M_PI;
}

/// A simple function for use in wrappers where macros cause problems
inline double get_HUGE(){ return _HUGE; }

#if defined(_MSC_VER)
// Microsoft version of math.h doesn't include acosh or asinh, so we just define them here.
// It was included from Visual Studio 2013
#if _MSC_VER < 1800
static double acosh(double x)
{
    return log(x + sqrt(x*x - 1.0));
}
static double asinh(double value)
{
    if(value>0){
        return log(value + sqrt(value * value + 1));
    }
    else{
        return -log(-value + sqrt(value * value + 1));
    }
}
#endif
#endif

#if defined(__powerpc__)
// PPC version of math.h doesn't include acosh or asinh, so we just define them here
static double acosh(double x)
{
        return log(x + sqrt(x*x - 1.0) );
}
static double asinh(double value)
{
    if(value>0){
        return log(value + sqrt(value * value + 1));
    }
    else{
        return -log(-value + sqrt(value * value + 1));
    }
}
#endif

#if defined(__powerpc__)
    #undef EOS
#endif

#endif