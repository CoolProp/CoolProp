#define _CRT_SECURE_NO_WARNINGS
#include <string>
#include <vector>
#include <cstdio>
#include <cstdarg>
#include <stdlib.h>
#include <memory>
#include "math.h"
#include "stdio.h"
#include "float.h"
#include <string.h>
#include "CoolPropTools.h"
#include "MatrixMath.h"
#include "Exceptions.h"
#include "crossplatform_shared_ptr.h"

#if !defined(__ISWINDOWS__)
#include <ftw.h>
#include <stdint.h>
#include <iostream>

static thread_local unsigned long long ftw_summer; // An evil global variable for the ftw function
int ftw_function(const char *fpath, const struct stat *sb, int tflag, struct FTW *ftwbuf){
   ftw_summer += sb->st_size;
   return 0;           /* To tell nftw() to continue */
}
unsigned long long CalculateDirSize(const std::string &path){
    ftw_summer = 0;
    int flags = 0 | FTW_DEPTH | FTW_PHYS;
    nftw(path.c_str(), ftw_function, 20, flags);
    double temp = ftw_summer;
    ftw_summer = 0;
    return temp;
}
#endif

double root_sum_square(const std::vector<double> &x)
{
    double sum = 0;
    for (unsigned int i=0; i<x.size(); i++)
    {
        sum += pow(x[i],2);
    }
    return sqrt(sum);
}
std::string format(const char* fmt, ...)
{
    const int size = 512;
    struct deleter{ static void delarray(void* p) { delete[] p; } }; // to use delete[]
    shared_ptr<char> buffer(new char[size], deleter::delarray); // I'd prefer unique_ptr, but it's only available since c++11
    va_list vl;
    va_start(vl,fmt);
    int nsize = vsnprintf(buffer.get(),size,fmt,vl);
    if(size<=nsize){//fail delete buffer and try again
        buffer.reset(new char[++nsize], deleter::delarray);//+1 for /0
        nsize = vsnprintf(buffer.get(),nsize,fmt,vl);
    }
    va_end(vl);
    return buffer.get();
}

std::vector<std::string> strsplit(const std::string &s, char del)
{
    std::vector<std::string> v;
    std::string::const_iterator i1 = s.begin(), i2;
    while (true){
        i2 = std::find(i1, s.end(), del);
        v.push_back(std::string(i1, i2));
        if (i2 == s.end())
            break;
        i1 = i2+1;
    }
    return v;
}
    
double interp1d(const std::vector<double> *x, const std::vector<double> *y, double x0)
{
    std::size_t i,L,R,M;
    L=0;
    R=(*x).size()-1;
    M=(L+R)/2;
    // Use interval halving to find the indices which bracket the density of interest
    while (R-L>1)
    {
        if (x0 >= (*x)[M])
        { L=M; M=(L+R)/2; continue;}
        if (x0 < (*x)[M])
        { R=M; M=(L+R)/2; continue;}
    }
    i=L;
    if (i<(*x).size()-2)
    {
        // Go "forwards" with the interpolation range
        return QuadInterp((*x)[i],(*x)[i+1],(*x)[i+2],(*y)[i],(*y)[i+1],(*y)[i+2],x0);
    }
    else
    {
        // Go "backwards" with the interpolation range
        return QuadInterp((*x)[i],(*x)[i-1],(*x)[i-2],(*y)[i],(*y)[i-1],(*y)[i-2],x0);
    }
}
double powInt(double x, int y)
{
    // Raise a double to an integer power
    // Overload not provided in math.h
    int i;
    double product=1.0;
    double x_in;
    int y_in;
    
    if (y==0)
    {
        return 1.0;
    }
    
    if (y<0)
    {
        x_in=1/x;
        y_in=-y;
    }
    else
    {
        x_in=x;
        y_in=y;
    }

    if (y_in==1)
    {
        return x_in;
    }    
    
    product=x_in;
    for (i=1;i<y_in;i++)
    {
        product=product*x_in;
    }
    return product;
}

void MatInv_2(double A[2][2] , double B[2][2])
{
    double Det;
    //Using Cramer's Rule to solve

    Det=A[0][0]*A[1][1]-A[1][0]*A[0][1];
    B[0][0]=1.0/Det*A[1][1];
    B[1][1]=1.0/Det*A[0][0];
    B[1][0]=-1.0/Det*A[1][0];
    B[0][1]=-1.0/Det*A[0][1];
}


std::string get_file_contents(const char *filename)
{
    std::ifstream in(filename, std::ios::in | std::ios::binary);
    if (in)
    {
        std::string contents;
        in.seekg(0, std::ios::end);
        contents.resize((unsigned int) in.tellg());
        in.seekg(0, std::ios::beg);
        in.read(&contents[0], contents.size());
        in.close();
        return(contents);
    }
    throw(errno);
}

void solve_cubic(double a, double b, double c, double d, int &N, double &x0, double &x1, double &x2)
{
    // 0 = ax^3 + b*x^2 + c*x + d
    
    // First check if the "cubic" is actually a second order or first order curve
    if (std::abs(a) < 10*DBL_EPSILON){
        if (std::abs(b) < 10*DBL_EPSILON){
            // Linear solution if a = 0 and b = 0
            x0 = -d/c;
            N = 1;
            return;
        }
        else{
            // Quadratic solution(s) if a = 0 and b != 0
            x0 = (-c+sqrt(c*c-4*b*d))/(2*b);
            x1 = (-c-sqrt(c*c-4*b*d))/(2*b);
            N = 2;
            return;
        }
    }
    
    // Ok, it is really a cubic

    // Discriminant
    double DELTA = 18*a*b*c*d-4*b*b*b*d+b*b*c*c-4*a*c*c*c-27*a*a*d*d;
    // Coefficients for the depressed cubic t^3+p*t+q = 0
    double p = (3*a*c-b*b)/(3*a*a);
    double q = (2*b*b*b-9*a*b*c+27*a*a*d)/(27*a*a*a);

    if (DELTA<0)
    {
        // One real root
        double t0;
        if (4*p*p*p+27*q*q>0 && p<0)
        {
            t0 = -2.0*std::abs(q)/q*sqrt(-p/3.0)*cosh(1.0/3.0*acosh(-3.0*std::abs(q)/(2.0*p)*sqrt(-3.0/p)));
        }
        else
        {
            t0 = -2.0*sqrt(p/3.0)*sinh(1.0/3.0*asinh(3.0*q/(2.0*p)*sqrt(3.0/p)));
        }
        N = 1;
        x0 = t0-b/(3*a);
        x1 = t0-b/(3*a);
        x2 = t0-b/(3*a);
    }
    else //(DELTA>0)
    {
        // Three real roots
        double t0 = 2.0*sqrt(-p/3.0)*cos(1.0/3.0*acos(3.0*q/(2.0*p)*sqrt(-3.0/p))-0*2.0*M_PI/3.0);
        double t1 = 2.0*sqrt(-p/3.0)*cos(1.0/3.0*acos(3.0*q/(2.0*p)*sqrt(-3.0/p))-1*2.0*M_PI/3.0);
        double t2 = 2.0*sqrt(-p/3.0)*cos(1.0/3.0*acos(3.0*q/(2.0*p)*sqrt(-3.0/p))-2*2.0*M_PI/3.0);

        N = 3;
        x0 = t0-b/(3*a);
        x1 = t1-b/(3*a);
        x2 = t2-b/(3*a);
    }
}

std::string strjoin(const std::vector<std::string> &strings, const std::string &delim)
{
    // Empty input vector
    if (strings.empty()){return "";}

    std::string output = strings[0];
    for (unsigned int i = 1; i < strings.size(); i++)
    {
        output += format("%s%s",delim.c_str(),strings[i].c_str());
    }
    return output;
}


SplineClass::SplineClass()
{
    Nconstraints = 0;
    A.resize(4, std::vector<double>(4, 0));
    B.resize(4,0);
}
bool SplineClass::build()
{
    if (Nconstraints == 4)
    {
        std::vector<double> abcd = CoolProp::linsolve(A,B);
        a = abcd[0];
        b = abcd[1];
        c = abcd[2];
        d = abcd[3];
        return true;
    }
    else
    {
        throw CoolProp::ValueError(format("Number of constraints[%d] is not equal to 4", Nconstraints));
    }
}
bool SplineClass::add_value_constraint(double x, double y)
{
    int i = Nconstraints;
    if (i == 4)
        return false;
    A[i][0] = x*x*x;
    A[i][1] = x*x;
    A[i][2] = x;
    A[i][3] = 1;
    B[i] = y;
    Nconstraints++;
    return true;
}
void SplineClass::add_4value_constraints(double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4)
{
    add_value_constraint(x1, y1);
    add_value_constraint(x2, y2);
    add_value_constraint(x3, y3);
    add_value_constraint(x4, y4);
}
bool SplineClass::add_derivative_constraint(double x, double dydx)
{
    int i = Nconstraints;
    if (i == 4)
        return false;
    A[i][0] = 3*x*x;
    A[i][1] = 2*x;
    A[i][2] = 1;
    A[i][3] = 0;
    B[i] = dydx;
    Nconstraints++;
    return true;
}
double SplineClass::evaluate(double x)
{
    return a*x*x*x+b*x*x+c*x+d;
}