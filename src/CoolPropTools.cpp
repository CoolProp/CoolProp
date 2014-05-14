#define _CRT_SECURE_NO_WARNINGS
#include <string>
#include <vector>
#include <cstdio>
#include <cstdarg>
#include <stdlib.h>
#include "math.h"
#include "stdio.h"
#include "float.h"
#include <string.h>
#include "CoolPropTools.h"

double root_sum_square(std::vector<double> x)
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
    int size = 512;
    char* buffer = 0;
    buffer = new char[size];
    va_list vl;
    va_start(vl,fmt);
    int nsize = vsnprintf(buffer,size,fmt,vl);
    if(size<=nsize){//fail delete buffer and try again
        delete buffer; buffer = 0;
        buffer = new char[nsize+1];//+1 for /0
        nsize = vsnprintf(buffer,size,fmt,vl);
    }
    std::string ret(buffer);
    va_end(vl);
    delete[] buffer;
    return ret;
}

std::vector<std::string> strsplit(std::string s, char del)
{
	int iL = 0, iR = 0, N;
	N = s.size();
	std::vector<std::string> v;
	// Find the first instance of the delimiter
	iR = s.find_first_of(del);
	// Delimiter not found, return the same string again
	if (iR < 0){
		v.push_back(s);
		return v;
	}
	while (iR != N-1)
	{
		v.push_back(s.substr(iL,iR-iL));
		iL = iR;
		iR = s.find_first_of(del,iR+1);
		// Move the iL to the right to avoid the delimiter
		iL += 1;
		if (iR == (int)std::string::npos)
		{
			v.push_back(s.substr(iL,N-iL));
			break;
		}
	}
	return v;
}
    
double interp1d(std::vector<double> *x, std::vector<double> *y, double x0)
{
	unsigned int i,L,R,M;
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
			t0 = -2.0*fabs(q)/q*sqrt(-p/3.0)*cosh(1.0/3.0*acosh(-3.0*fabs(q)/(2.0*p)*sqrt(-3.0/p)));
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

std::string strjoin(std::vector<std::string> strings, std::string delim)
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
