#ifndef COOLPROPTOOLS_H
#define COOLPROPTOOLS_H

    #define _CRT_SECURE_NO_WARNINGS

    #include "PlatformDetermination.h"
    #include "Exceptions.h"
    #include <string>
    #include <vector>
    #include <cmath>
    #include "float.h"
    #include <iterator>
    #include <algorithm>
    #include <functional>
    #include <cctype>
    #include <map>
    #include <locale>
    #include <fstream>
    #include <cerrno>
    #include <numeric>
    #include <set>
    #include <sys/types.h>
    #include <sys/stat.h>    

    // This will kill the horrible min and max macros 
    #ifndef NOMINMAX
        #define NOMINMAX
    #endif
    
    #if defined(__ISWINDOWS__)
        #define UNICODE
        #define _UNICODE
        #include "Windows.h"
        #include <windows.h> // for the CreateDirectory function
    #else
        #include <unistd.h>
        #include <pwd.h>
    #endif

    #ifndef __has_feature         // Optional of course.
        #define __has_feature(x) 0  // Compatibility with non-clang compilers.
    #endif  
    
    // see http://stackoverflow.com/questions/18298280/how-to-declare-a-variable-as-thread-local-portably
    #ifndef thread_local
        #if __STDC_VERSION__ >= 201112 && !defined __STDC_NO_THREADS__
        # define thread_local _Thread_local
        #elif defined _WIN32 && ( \
              defined _MSC_VER || \
              defined __ICL || \
              defined __DMC__ || \
              defined __BORLANDC__ )
            #define thread_local __declspec(thread) 
        #elif defined(__ISAPPLE__) && (defined(__llvm__) || defined(__clang__)) && !__has_feature(cxx_thread_local)
            #define thread_local 
        /* note that ICC (linux) and Clang are covered by __GNUC__ */
        #elif defined __GNUC__ || \
              defined __SUNPRO_C || \
              defined __xlC__
            #define thread_local __thread
        #else
            #error "Cannot define thread_local"
//          #define thread_local
        #endif
    #endif

    typedef long double CoolPropDbl;

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
    
    /// Copy string to wstring
    /// Dangerous if the string has non-ASCII characters; from http://stackoverflow.com/a/8969776/1360263 
    inline void StringToWString(const std::string &s, std::wstring &ws)
    {
        ws = std::wstring(s.begin(), s.end());
    }

    #if defined(__ISWINDOWS__)
    /// From http://stackoverflow.com/a/17827724/1360263
    inline bool IsBrowsePath(const std::wstring& path)
    {
        return (path == L"." || path == L"..");
    }
    inline unsigned long long CalculateDirSize(const std::wstring &path, std::vector<std::wstring> *errVect = NULL)
    {
        unsigned long long size = 0;
        WIN32_FIND_DATAW data;
        HANDLE sh = NULL;
        sh = FindFirstFileW((path + L"\\*").c_str(), &data);

        if (sh == INVALID_HANDLE_VALUE )
        {
            //if we want, store all happened error  
            if (errVect != NULL)
                errVect ->push_back(path);
            return size;
        }

        do
        {
            // skip current and parent
            if (!IsBrowsePath(data.cFileName))
            {
                // if found object is ...
                if (data.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)
                    // directory, then search it recursievly
                    size += CalculateDirSize(path + L"\\" + data.cFileName, NULL);
                else
                    // otherwise get object size and add it to directory size
                    size += data.nFileSizeHigh * (unsigned long long)(MAXDWORD) + data.nFileSizeLow;
            }

        } while (FindNextFileW(sh, &data)); // do

        FindClose(sh);

        return size;
    } 
    #else
    /// Get the size of a directory in bytes
    unsigned long long CalculateDirSize(const std::string &path);
    #endif

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
    std::vector<std::string> strsplit(const std::string &s, char del);

    inline std::string upper(std::string str)
    {
        std::transform(str.begin(), str.end(), str.begin(), ::toupper);
        return str;
    }
	
	inline std::string lower(std::string str)
    {
        std::transform(str.begin(), str.end(), str.begin(), ::tolower);
        return str;
    }

    std::string strjoin(const std::vector<std::string> &strings, const std::string &delim);

    void MatInv_2(double A[2][2] , double B[2][2]);

    double root_sum_square(const std::vector<double> &x);
    double interp1d(const std::vector<double> *x, const std::vector<double> *y, double x0);
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
    template<class T> T LinearInterp(const std::vector<T> &x, const std::vector<T> &y, std::size_t i0, std::size_t i1, T val)
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
    template<class T> T QuadInterp(const std::vector<T> &x, const std::vector<T> &y, std::size_t i0, std::size_t i1, std::size_t i2, T val)
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
    template<class T> T CubicInterp(const std::vector<T> &x, const std::vector<T> &y, std::size_t i0, std::size_t i1, std::size_t i2, std::size_t i3, T val)
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

    inline int Kronecker_delta(int i, int j){if (i == j) {return 1;} else {return 0;}};

    class Dictionary
    {
    private:
        typedef std::map<std::string, double> numbers_map;
        numbers_map numbers;
        typedef std::map<std::string, std::string> strings_map;
        strings_map strings;
        typedef std::map<std::string, std::vector<double> > double_vectors_map;
        double_vectors_map double_vectors;
        typedef std::map<std::string, std::vector<std::string> > string_vectors_map;
        string_vectors_map string_vectors;
    public:
        Dictionary(){};
        bool is_empty(void) const {return numbers.empty() && strings.empty() && double_vectors.empty() && string_vectors.empty();}
        void add_string(const std::string &s1, const std::string &s2){ strings.insert(std::pair<std::string, std::string>(s1, s2));}
        void add_number(const std::string &s1, double d){ numbers.insert(std::pair<std::string, double>(s1, d));}
        void add_double_vector(const std::string &s1, const std::vector<double> &d){ double_vectors.insert(std::pair<std::string, std::vector<double> >(s1, d));}
        void add_string_vector(const std::string &s1, const std::vector<std::string> &d){ string_vectors.insert(std::pair<std::string, std::vector<std::string> >(s1, d));}
        std::string get_string(const std::string &s) const
        {
            strings_map::const_iterator i = strings.find(s);
            if (i != strings.end()){
                return i->second;
            }
            else{
                throw CoolProp::ValueError(format("%s could not be matched in get_string",s.c_str()));
            }
        };
        double get_double(const std::string &s) const
        {
            numbers_map::const_iterator i = numbers.find(s);
            if (i != numbers.end()){
                return i->second;
            }
            else{
                throw CoolProp::ValueError(format("%s could not be matched in get_number",s.c_str()));
            }
        };
        double get_number(const std::string &s) const
        {
            return get_double(s);
        };
        const std::vector<double>& get_double_vector(const std::string &s) const
        {
            double_vectors_map::const_iterator i = double_vectors.find(s);
            if (i != double_vectors.end()){
                return i->second;
            }
            else{
                throw CoolProp::ValueError(format("%s could not be matched in get_double_vector",s.c_str()));
            }
        };
        const std::vector<std::string>& get_string_vector(const std::string &s) const
        {
            string_vectors_map::const_iterator i = string_vectors.find(s);
            if (i != string_vectors.end()){
                return i->second;
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

template <typename T> void bisect_vector(const std::vector<T> &vec, T val, std::size_t &i)
{
    T rL, rM, rR;
    std::size_t N = vec.size(), L = 0, R = N-1, M = (L+R)/2;
    rL = vec[L] - val; rR = vec[R] - val;
    while (R - L > 1){
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
   
       /// Get the user's home directory;  It is believed that is is always a place that files can be written
    inline std::string get_home_dir(void)
    {
        // See http://stackoverflow.com/questions/2552416/how-can-i-find-the-users-home-dir-in-a-cross-platform-manner-using-c
        #if defined(__ISLINUX__)
            char *home = NULL;
            home = getenv("HOME");
            return std::string(home);
        #elif defined(__ISAPPLE__)
            char *home = NULL;
            home = getenv("HOME");
            if (home==NULL) {
              struct passwd* pwd = getpwuid(getuid());
              if (pwd) {
                home = pwd->pw_dir;
              }
            }
            if (home==NULL) {
              throw CoolProp::NotImplementedError("Could not detect home directory.");
            } 
            return std::string(home);
        #elif defined(__ISWINDOWS__)
			
            char * pUSERPROFILE = getenv("USERPROFILE");
            if (pUSERPROFILE != NULL) {
                return std::string(pUSERPROFILE);
            } else {
                char * pHOMEDRIVE = getenv("HOMEDRIVE");
                char * pHOMEPATH = getenv("HOMEPATH");
                if (pHOMEDRIVE != NULL && pHOMEPATH != NULL) {
                    return std::string(pHOMEDRIVE) + std::string(pHOMEPATH);
                } else {
                    return std::string("");
                }
            }
        #else
            throw CoolProp::NotImplementedError("This function is not defined for your platform.");
        #endif
    };
    
    inline bool path_exists(const std::string &path)
    {
        #if defined(__ISWINDOWS__) // Defined for 32-bit and 64-bit windows
            struct _stat buf;
            // Get data associated with path using the windows libraries, 
            // and if you can (result == 0), the path exists
            if ( _stat( path.c_str(), &buf) == 0)
                return true;
            else
                return false;
        #elif defined(__ISLINUX__) || defined(__ISAPPLE__)
            struct stat st;
            if(lstat(path.c_str(),&st) == 0) {
                if(S_ISDIR(st.st_mode)) return true;
                if(S_ISREG(st.st_mode)) return true;
                return false;
            } else {
                return false;
            }
        #else
            throw CoolProp::NotImplementedError("This function is not defined for your platform.");
        #endif
    };

    inline void make_dirs(std::string file_path)
    {
        std::replace( file_path.begin(), file_path.end(), '\\', '/'); // replace all '\' with '/'
        
        #if defined(__ISWINDOWS__)
        const char sep = '\\'; // well, Windows (and DOS) allows forward slash "/", too :)
        #else
        const char sep = '/';
        #endif
        
        std::vector<std::string> pathsplit = strsplit(file_path,'/');
        std::string path = pathsplit[0]; // will throw if pathsplit.size() == 0
        for (std::size_t i = 0, sz = pathsplit.size(); i < sz; i++)
        {
            if (!path_exists(path))
            {
                #if defined(__ISWINDOWS__) // Defined for 32-bit and 64-bit windows
                    int errcode = CreateDirectoryA((LPCSTR)path.c_str(),NULL);
                    if (errcode == 0){
                        switch(GetLastError()){
                            case ERROR_ALREADY_EXISTS:
                                break;
                            case ERROR_PATH_NOT_FOUND:
                                throw CoolProp::ValueError(format("Unable to make the directory %s",path.c_str()));
                            default:
                                break;
                        }
                        
                    }
                #else
                    #if defined(__powerpc__)
                    #else
                        mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                    #endif
                #endif
            }
            if (i < (sz-1))
                path += format("%c%s", sep, pathsplit[i+1].c_str());
        }
    };
#endif
