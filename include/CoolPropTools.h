#ifndef COOLPROPTOOLS_H
#define COOLPROPTOOLS_H

    #define _CRT_SECURE_NO_WARNINGS

    #include "PlatformDetermination.h"

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

    #ifdef HUGE_VAL
    #  define _HUGE HUGE_VAL
    #else
    // GCC Version of huge value macro
    #ifdef HUGE 
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

	std::string strjoin(std::vector<std::string> strings, std::string delim);
	
	void MatInv_2(double A[2][2] , double B[2][2]);

	double root_sum_square(std::vector<double> x);
	double interp1d(std::vector<double> *x, std::vector<double> *y, double x0);
	double powInt(double x, int y);
    
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
	template<class T> T QuadInterp(std::vector<T> x, std::vector<T> y, int i0, int i1, int i2, T val)
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
    template<class T> T CubicInterp(std::vector<T> x, std::vector<T> y, int i0, int i1, int i2, int i3, T val)
    {
        return CubicInterp(x[i0],x[i1],x[i2],x[i3],y[i0],y[i1],y[i2],y[i3],val);
    };

    template<class T> T is_in_closed_range( T x1, T x2, T x)
    {
        if (x1 > x2)
        {
            std::swap(x1, x2);
        }

	    return (x >= x1 && x <= x2);
    };

	void solve_cubic(double a, double b, double c, double d, int &N, double &x0, double &x1, double &x2);
	
	inline double min3(double x1, double x2, double x3){return std::min(std::min(x1, x2), x3);};
	inline double max3(double x1, double x2, double x3){return std::max(std::max(x1, x2), x3);};

	inline bool double_equal(double a, double b){return fabs(a - b) <= 1 * DBL_EPSILON * std::max(fabs(a), fabs(b));};

    template<class T> T max_abs_value(std::vector<T> x)
    {
        T max = 0;
        std::size_t N = x.size();
        for (std::size_t i = 0; i < N; ++i)
        {
            T axi = fabs(x[i]);
            if (axi > max){ max = axi; }
        }
        return max;
    }

	inline int Kronecker_delta(int i, int j){if (i == j) {return 1;} else {return 0;}};

    class Dictionary
    {
    private:
        std::map<std::string, double> numbers;
        std::map<std::string, std::string> strings;
        std::map<std::string, std::vector<double> > double_vectors;
    public:
        Dictionary(){};
        void add_string(std::string s1, std::string s2){ strings.insert(std::pair<std::string, std::string>(s1, s2));}
        void add_number(std::string s1, double d){ numbers.insert(std::pair<std::string, double>(s1, d));}
        void add_double_vector(std::string s1, std::vector<double> d){ double_vectors.insert(std::pair<std::string, std::vector<double> >(s1, d));}
        std::string get_string(std::string s)
        {
            if (strings.find(s) != strings.end()){ return strings[s]; } else{ throw std::exception(); }
        };
        double get_number(std::string s)
        {
            if (numbers.find(s) != numbers.end()){ return numbers[s]; } else{ throw std::exception(); }
        };
        std::vector<double> get_double_vector(std::string s)
        {
            if (double_vectors.find(s) != double_vectors.end()){ return double_vectors[s]; } else{ throw std::exception(); }
        };
    };

#define CATCH_ALL_ERRORS_RETURN_HUGE(x)    try{                                                                                \
                                                x                                                                              \
                                            }                                                                                  \
                                            catch(const std::exception& e){                                                    \
                                                std::cout << e.what() << std::endl;                                            \
                                                return _HUGE;                                                                  \
                                            }                                                                                  \
                                            catch(...){                                                                        \
                                                return _HUGE;                                                                  \
                                            }                                                        

#endif
