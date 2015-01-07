#ifndef COOLPROPTOOLS_H
#define COOLPROPTOOLS_H

    #define _CRT_SECURE_NO_WARNINGS

    #include "PlatformDetermination.h"
    #include "Exceptions.h"

    #include <string>
    #include <vector>
    #include <cmath>
    #include "float.h"

    #ifndef M_PI
    #  define M_PI 3.14159265358979323846
    #endif

    #ifndef COOLPROP_OK
    #define COOLPROP_OK 1
    #endif

    #if defined(HUGE_VAL) && !defined(_HUGE)
        # define _HUGE HUGE_VAL
    #else
        // GCC Version of huge value macro
        #if defined(HUGE) && !defined(_HUGE)
        #  define _HUGE HUGE
        #endif
    #endif

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
        #undef min
        #undef max
        #undef EOS
    #endif

    inline bool ValidNumber(double x)
    {
        // Idea from http://www.johndcook.com/IEEE_exceptions_in_cpp.html
        return (x <= DBL_MAX && x >= -DBL_MAX);
    };

    /// Define the deprecated macro to give compile-time warnings
    #ifdef __GNUC__
        #define DEPRECATED(func) func __attribute__ ((deprecated))
    #elif defined(_MSC_VER)
        #define DEPRECATED(func) __declspec(deprecated) func
    #else
        #pragma message("WARNING: You need to implement DEPRECATED for this compiler")
        #define DEPRECATED(func) func
    #endif

    #include <algorithm>
    #include <functional>
    #include <cctype>
    #include <map>
    #include <locale>
    #include <fstream>
    #include <cerrno>
    #include <numeric>

    /// The following code for the trim functions was taken from http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring
    // trim from start
    inline std::string &strlstrip(std::string &s) {
            s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
            return s;
    }
    // trim from end
    inline std::string &strrstrip(std::string &s) {
            s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
            return s;
    }
    // trim from both ends
    inline std::string &strstrip(std::string &s) {
            return strlstrip(strrstrip(s));
    }

    // Get all the contents of a file and dump into a STL string
    // Thanks to http://stackoverflow.com/questions/2602013/read-whole-ascii-file-into-c-stdstring
    std::string get_file_contents(const char *filename);

    // Missing string printf
    std::string format(const char* fmt, ...);
    // Missing string split - like in Python
    std::vector<std::string> strsplit(std::string s, char del);

    inline std::string upper(const std::string str_)
    {
        std::string str = str_;
        std::transform(str.begin(), str.end(), str.begin(), ::toupper);
        return str;
    }
	
	inline std::string lower(const std::string str_)
    {
        std::string str = str_;
        std::transform(str.begin(), str.end(), str.begin(), ::tolower);
        return str;
    }

    std::string strjoin(std::vector<std::string> strings, std::string delim);

    void MatInv_2(double A[2][2] , double B[2][2]);

    double root_sum_square(std::vector<double> x);
    double interp1d(std::vector<double> *x, std::vector<double> *y, double x0);
    double powInt(double x, int y);
    
    #define POW2(x) ((x)*(x))
    #define POW3(x) ((x)*(x)*(x))
    #define POW4(x) ((x)*(x)*(x)*(x))
    #define POW5(x) ((x)*(x)*(x)*(x)*(x))
    #define POW6(x) ((x)*(x)*(x)*(x)*(x)*(x))
    #define POW7(x) ((x)*(x)*(x)*(x)*(x)*(x)*(x))

    template<class T> T LinearInterp(T x0, T x1, T y0, T y1, T x)
    {
        return (y1-y0)/(x1-x0)*(x-x0)+y0;
    };
    template<class T> T LinearInterp(std::vector<T> x, std::vector<T> y, std::size_t i0, std::size_t i1, T val)
    {
        return LinearInterp(x[i0],x[i1],y[i0],y[i1],val);
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
    template<class T> T QuadInterp(std::vector<T> x, std::vector<T> y, std::size_t i0, std::size_t i1, std::size_t i2, T val)
    {
        return QuadInterp(x[i0],x[i1],x[i2],y[i0],y[i1],y[i2],val);
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
        return L0*f0+L1*f1+L2*f2+L3*f3;
    };
    template<class T> T CubicInterp(std::vector<T> x, std::vector<T> y, std::size_t i0, std::size_t i1, std::size_t i2, std::size_t i3, T val)
    {
        return CubicInterp(x[i0],x[i1],x[i2],x[i3],y[i0],y[i1],y[i2],y[i3],val);
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

    inline double min3(double x1, double x2, double x3){return std::min(std::min(x1, x2), x3);};
    inline double max3(double x1, double x2, double x3){return std::max(std::max(x1, x2), x3);};

    inline bool double_equal(double a, double b){return std::abs(a - b) <= 1 * DBL_EPSILON * std::max(std::abs(a), std::abs(b));};

    template<class T> T max_abs_value(std::vector<T> x)
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
    
    template<class T> T min_abs_value(std::vector<T> x)
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

    inline int Kronecker_delta(int i, int j){if (i == j) {return 1;} else {return 0;}};

    class Dictionary
    {
    private:
        std::map<std::string, double> numbers;
        std::map<std::string, std::string> strings;
        std::map<std::string, std::vector<double> > double_vectors;
        std::map<std::string, std::vector<std::string> > string_vectors;
    public:
        Dictionary(){};
        bool is_empty(void){return numbers.empty() && strings.empty() && double_vectors.empty() && string_vectors.empty();}
        void add_string(std::string s1, std::string s2){ strings.insert(std::pair<std::string, std::string>(s1, s2));}
        void add_number(std::string s1, double d){ numbers.insert(std::pair<std::string, double>(s1, d));}
        void add_double_vector(std::string s1, std::vector<double> d){ double_vectors.insert(std::pair<std::string, std::vector<double> >(s1, d));}
        void add_string_vector(std::string s1, std::vector<std::string> d){ string_vectors.insert(std::pair<std::string, std::vector<std::string> >(s1, d));}
        std::string get_string(std::string s)
        {
            if (strings.find(s) != strings.end()){ 
                return strings[s]; 
            } 
            else{ 
                throw CoolProp::ValueError(format("%s could not be matched in get_string",s.c_str())); 
            }
        };
        double get_double(std::string s)
        {
            if (numbers.find(s) != numbers.end()){ 
                return numbers[s]; 
            } 
            else{ 
                throw CoolProp::ValueError(format("%s could not be matched in get_number",s.c_str())); 
            }
        };
        double get_number(std::string s)
        {
            return get_double(s);
        };
        std::vector<double> get_double_vector(std::string s)
        {
            if (double_vectors.find(s) != double_vectors.end()){ 
                return double_vectors[s]; 
            } 
            else{ 
                throw CoolProp::ValueError(format("%s could not be matched in get_double_vector",s.c_str()));
            }
        };
        std::vector<std::string> get_string_vector(std::string s)
        {
            if (string_vectors.find(s) != string_vectors.end()){ 
                return string_vectors[s]; 
            } 
            else{ 
                throw CoolProp::ValueError(format("%s could not be matched in get_string_vector",s.c_str()));
            }
        };
    };
    /// Utility function to clear a std::map of pointers
    //http://stackoverflow.com/questions/569110/why-is-memory-still-accessible-after-stdmapclear-is-called
    template <typename M> void freeClear( M & amap ) {
        for ( typename M::iterator it = amap.begin(); it != amap.end(); ++it ) {
            delete it->second;
        }
        amap.clear();
    }

    /// Make a linearly spaced vector of points
    template <typename T> std::vector<T> linspace(T xmin, T xmax, int n) {
        std::vector<T> x(n, 0.0);
        
        for ( std::size_t i = 0;  i < n; ++i) {
            x[i] = (xmax-xmin)/(n-1)*i+xmin;
        }
        return x;
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

#define CATCH_ALL_ERRORS_RETURN_HUGE(x)    try{                                                                                \
                                                x                                                                              \
                                            }                                                                                  \
                                            catch(const std::exception& e){                                                    \
                                                return _HUGE;                                                                  \
                                            }                                                                                  \
                                            catch(...){                                                                        \
                                                return _HUGE;                                                                  \
                                            }

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
        SplineClass();
        bool build(void);
        bool add_value_constraint(double x, double y);
        void add_4value_constraints(double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4);
        bool add_derivative_constraint(double x, double dydx);
        double evaluate(double x);
    };
#endif
