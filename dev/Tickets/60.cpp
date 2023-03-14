#include "../../CoolProp/CoolPropTools.h"
#include "../../CoolProp/IncompBase.h"
#include "../../CoolProp/IncompLiquid.h"
#include "../../CoolProp/IncompSolution.h"

#if defined(__ISWINDOWS__)
#    include <Windows.h>
#    include "stdafx.h"
#else
#    include <sys/time.h>
#    include <ctime>
#    include <stdint.h>
#endif

#include <vector>
#include <math.h>
#include <string.h>
#include <iostream>
#include <map>
#include <utility>

/** A class to test the different implementations,
 *  exposes more functions than the other ones.
 */
class IncompressibleTest : public IncompressibleClass
{
   public:
    IncompressibleTest(){};
    ~IncompressibleTest(){};
    // Some functions need top be overwritten!

   public:
    double testSi(std::vector<double> const& coefficients, double x) {
        return simplePolynomial(coefficients, x);
    }

    double testHo(std::vector<double> const& coefficients, double x) {
        return baseHorner(coefficients, x);
    }

    double testSiInt(std::vector<double> const& coefficients, double T) {
        return simplePolynomialInt(coefficients, T);
    }

    double testHoInt(std::vector<double> const& coefficients, double T) {
        return baseHornerInt(coefficients, T);
    }

    double test2sInt(std::vector<double> const& coefficients, double T) {
        return integrateIn2Steps(coefficients, T);
    }

    double testSiInt(std::vector<double> const& coefficients, double T1, double T0) {
        return simplePolynomialInt(coefficients, T1, T0);
    }

    double testHoInt(std::vector<double> const& coefficients, double T1, double T0) {
        return baseHornerInt(coefficients, T1, T0);
    }

    double test2sInt(std::vector<double> const& coefficients, double T1, double T0) {
        return integrateIn2Steps(coefficients, T1, T0);
    }

    double testSiFra(std::vector<double> const& coefficients, double T) {
        return simpleFracInt(coefficients, T);
    }

    double testHoFra(std::vector<double> const& coefficients, double T) {
        return baseHornerFra(coefficients, T);
    }

    double test2sFra(std::vector<double> const& coefficients, double T) {
        return fracIntIn2Steps(coefficients, T);
    }

    double testSiFra(std::vector<double> const& coefficients, double T1, double T0) {
        return simpleFracInt(coefficients, T1, T0);
    }

    double testHoFra(std::vector<double> const& coefficients, double T1, double T0) {
        return baseHornerFra(coefficients, T1, T0);
    }

    double test2sFra(std::vector<double> const& coefficients, double T1, double T0) {
        return fracIntIn2Steps(coefficients, T1, T0);
    }

    // And the same in 2D
    double testSi(std::vector<std::vector<double>> const& coefficients, double x, double T) {
        return simplePolynomial(coefficients, x, T);
    }

    double testHo(std::vector<std::vector<double>> const& coefficients, double x, double T) {
        return baseHorner(coefficients, x, T);
    }

    double testSiInt(std::vector<std::vector<double>> const& coefficients, double x, double T) {
        return simplePolynomialInt(coefficients, x, T);
    }

    double testHoInt(std::vector<std::vector<double>> const& coefficients, double x, double T) {
        return baseHornerInt(coefficients, x, T);
    }

    double test2sInt(std::vector<std::vector<double>> const& coefficients, double x, double T) {
        return integrateIn2Steps(coefficients, x, T, false);
    }

    double testSiInt(std::vector<std::vector<double>> const& coefficients, double x, double T1, double T0) {
        return simplePolynomialInt(coefficients, x, T1, T0);
    }

    double testHoInt(std::vector<std::vector<double>> const& coefficients, double x, double T1, double T0) {
        return baseHornerInt(coefficients, x, T1, T0);
    }

    double test2sInt(std::vector<std::vector<double>> const& coefficients, double x, double T1, double T0) {
        return integrateIn2Steps(coefficients, x, T1, T0);
    }

    double testSiFra(std::vector<std::vector<double>> const& coefficients, double x, double T) {
        return simpleFracInt(coefficients, x, T);
    }

    double testHoFra(std::vector<std::vector<double>> const& coefficients, double x, double T) {
        return baseHornerFra(coefficients, x, T);
    }

    double test2sFra(std::vector<std::vector<double>> const& coefficients, double x, double T) {
        return fracIntIn2Steps(coefficients, x, T);
    }

    double testSiFra(std::vector<std::vector<double>> const& coefficients, double x, double T1, double T0) {
        return simpleFracInt(coefficients, x, T1, T0);
    }

    double testHoFra(std::vector<std::vector<double>> const& coefficients, double x, double T1, double T0) {
        return baseHornerFra(coefficients, x, T1, T0);
    }

    double test2sFra(std::vector<std::vector<double>> const& coefficients, double x, double T1, double T0) {
        return fracIntIn2Steps(coefficients, x, T1, T0);
    }
};

/* Returns the amount of milliseconds elapsed since the UNIX epoch. Works on both
 * windows and linux.
 * Taken from http://stackoverflow.com/questions/1861294/how-to-calculate-execution-time-of-a-code-snippet-in-c
 * */
uint64_t getTimeValue() {
#if defined(__ISWINDOWS__)
    /* Windows */
    FILETIME ft;
    LARGE_INTEGER li;

    /* Get the amount of 100 nano seconds intervals elapsed since January 1, 1601 (UTC) and copy it
	 * to a LARGE_INTEGER structure. */
    GetSystemTimeAsFileTime(&ft);
    li.LowPart = ft.dwLowDateTime;
    li.HighPart = ft.dwHighDateTime;

    uint64_t ret = li.QuadPart;
    ret -= 116444736000000000LL; /* Convert from file time to UNIX epoch time. */
    ret /= 10;                   /* From 100 nano seconds (10^-7) to 1 microsecond (10^-6) intervals */

    return ret;
#else
    /* Linux */
    struct timeval tv;

    gettimeofday(&tv, NULL);

    uint64_t ret = tv.tv_usec;
    /* In micro seconds (10^-6), add the seconds (10^0)
	 * after converting them to microseconds (10^-6) */
    ret += (tv.tv_sec * 1000000);

    return ret;
#endif
}

std::vector<std::vector<double>> makeMatrix(std::vector<double> const& coefficients) {
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
    tmpVector.push_back(coefficients[1]);
    tmpVector.push_back(coefficients[7]);
    tmpVector.push_back(coefficients[12]);
    tmpVector.push_back(coefficients[16]);
    matrix.push_back(tmpVector);

    tmpVector.clear();
    tmpVector.push_back(coefficients[2]);
    tmpVector.push_back(coefficients[8]);
    tmpVector.push_back(coefficients[13]);
    tmpVector.push_back(coefficients[17]);
    matrix.push_back(tmpVector);

    tmpVector.clear();
    tmpVector.push_back(coefficients[3]);
    tmpVector.push_back(coefficients[9]);
    tmpVector.push_back(coefficients[14]);
    tmpVector.push_back(0.0);
    matrix.push_back(tmpVector);

    tmpVector.clear();
    tmpVector.push_back(coefficients[4]);
    tmpVector.push_back(coefficients[10]);
    tmpVector.push_back(0.0);
    tmpVector.push_back(0.0);
    matrix.push_back(tmpVector);

    tmpVector.clear();
    tmpVector.push_back(coefficients[5]);
    tmpVector.push_back(0.0);
    tmpVector.push_back(0.0);
    tmpVector.push_back(0.0);
    matrix.push_back(tmpVector);

    tmpVector.clear();
    return matrix;
}

std::map<std::string, double> testObject(IncompressibleTest* fluid, std::vector<std::vector<double>> coeffs, int exponent, bool print) {

    std::map<std::string, double> results;

    uint64_t runs = 10;
    for (int i = 1; i < exponent; i++) {
        runs *= 10;
    }

    if (print) std::cout << "Runs: " << runs << std::endl;

    std::vector<double> simCoeffs;
    if (coeffs.size() > 0) {
        simCoeffs = coeffs[1];
    } else {
        results["error"] = -_HUGE;
        return results;
    }

    if (print) {
        std::cout << "simCoeffs: " << simCoeffs[0];
        for (int i = 1; i < simCoeffs.size(); i++) {
            std::cout << ", " << simCoeffs[i];
        }
        std::cout << std::endl;
    }

    uint64_t start;
    uint64_t end;
    double time;

    double Tin = 280;
    double T0 = Tin - 10;
    double xin = 0.25;

    // Testing the 1D functions first
    fluid->setDebug(false);
    start = getTimeValue();
    for (int i = 0; i < runs; i++) {
        fluid->testSi(simCoeffs, Tin + i / runs);
    }
    end = getTimeValue();
    time = (end - start) * 1e3 / runs;
    results["1D std sim simple"] = time;
    fluid->setDebug(print);
    fluid->testSi(simCoeffs, Tin);

    fluid->setDebug(false);
    start = getTimeValue();
    for (int i = 0; i < runs; i++) {
        fluid->testHo(simCoeffs, Tin + i / runs);
    }
    end = getTimeValue();
    time = (end - start) * 1e3 / runs;
    results["1D std sim Horner"] = time;
    fluid->setDebug(print);
    fluid->testHo(simCoeffs, Tin);

    // Testing the 1D integrators
    fluid->setDebug(false);
    start = getTimeValue();
    for (int i = 0; i < runs; i++) {
        fluid->testSiInt(simCoeffs, Tin + i / runs);
    }
    end = getTimeValue();
    time = (end - start) * 1e3 / runs;
    results["1D int ind simple"] = time;
    fluid->setDebug(print);
    fluid->testSiInt(simCoeffs, Tin);

    fluid->setDebug(false);
    start = getTimeValue();
    for (int i = 0; i < runs; i++) {
        fluid->testHoInt(simCoeffs, Tin + i / runs);
    }
    end = getTimeValue();
    time = (end - start) * 1e3 / runs;
    results["1D int ind Horner"] = time;
    fluid->setDebug(print);
    fluid->testHoInt(simCoeffs, Tin);

    fluid->setDebug(false);
    start = getTimeValue();
    for (int i = 0; i < runs; i++) {
        fluid->test2sInt(simCoeffs, Tin + i / runs);
    }
    end = getTimeValue();
    time = (end - start) * 1e3 / runs;
    results["1D int ind 2Steps"] = time;
    fluid->setDebug(print);
    fluid->test2sInt(simCoeffs, Tin);

    // Testing the 1D definite integrators
    fluid->setDebug(false);
    start = getTimeValue();
    for (int i = 0; i < runs; i++) {
        fluid->testSiInt(simCoeffs, Tin + i / runs, T0);
    }
    end = getTimeValue();
    time = (end - start) * 1e3 / runs;
    results["1D int def simple"] = time;
    fluid->setDebug(print);
    fluid->testSiInt(simCoeffs, Tin, T0);

    fluid->setDebug(false);
    start = getTimeValue();
    for (int i = 0; i < runs; i++) {
        fluid->testHoInt(simCoeffs, Tin + i / runs, T0);
    }
    end = getTimeValue();
    time = (end - start) * 1e3 / runs;
    results["1D int def Horner"] = time;
    fluid->setDebug(print);
    fluid->testHoInt(simCoeffs, Tin, T0);

    fluid->setDebug(false);
    start = getTimeValue();
    for (int i = 0; i < runs; i++) {
        fluid->test2sInt(simCoeffs, Tin + i / runs, T0);
    }
    end = getTimeValue();
    time = (end - start) * 1e3 / runs;
    results["1D int def 2Steps"] = time;
    fluid->setDebug(print);
    fluid->test2sInt(simCoeffs, Tin, T0);

    // Testing the 1D fraction integrators
    fluid->setDebug(false);
    start = getTimeValue();
    for (int i = 0; i < runs; i++) {
        fluid->testSiFra(simCoeffs, Tin + i / runs);
    }
    end = getTimeValue();
    time = (end - start) * 1e3 / runs;
    results["1D fra ind simple"] = time;
    fluid->setDebug(print);
    fluid->testSiFra(simCoeffs, Tin);

    fluid->setDebug(false);
    start = getTimeValue();
    for (int i = 0; i < runs; i++) {
        fluid->testHoFra(simCoeffs, Tin + i / runs);
    }
    end = getTimeValue();
    time = (end - start) * 1e3 / runs;
    results["1D fra ind Horner"] = time;
    fluid->setDebug(print);
    fluid->testHoFra(simCoeffs, Tin);

    fluid->setDebug(false);
    start = getTimeValue();
    for (int i = 0; i < runs; i++) {
        fluid->test2sFra(simCoeffs, Tin + i / runs);
    }
    end = getTimeValue();
    time = (end - start) * 1e3 / runs;
    results["1D fra ind 2Steps"] = time;
    fluid->setDebug(print);
    fluid->test2sFra(simCoeffs, Tin);

    // Testing the 1D definite fraction integrators
    fluid->setDebug(false);
    start = getTimeValue();
    for (int i = 0; i < runs; i++) {
        fluid->testSiFra(simCoeffs, Tin + i / runs, T0);
    }
    end = getTimeValue();
    time = (end - start) * 1e3 / runs;
    results["1D fra def simple"] = time;
    fluid->setDebug(print);
    fluid->testSiFra(simCoeffs, Tin, T0);

    fluid->setDebug(false);
    start = getTimeValue();
    for (int i = 0; i < runs; i++) {
        fluid->testHoFra(simCoeffs, Tin + i / runs, T0);
    }
    end = getTimeValue();
    time = (end - start) * 1e3 / runs;
    results["1D fra def Horner"] = time;
    fluid->setDebug(print);
    fluid->testHoFra(simCoeffs, Tin, T0);

    fluid->setDebug(false);
    start = getTimeValue();
    for (int i = 0; i < runs; i++) {
        fluid->test2sFra(simCoeffs, Tin + i / runs, T0);
    }
    end = getTimeValue();
    time = (end - start) * 1e3 / runs;
    results["1D fra def 2Steps"] = time;
    fluid->setDebug(print);
    fluid->test2sFra(simCoeffs, Tin, T0);

    // Testing the 2D functions
    fluid->setDebug(false);
    start = getTimeValue();
    for (int i = 0; i < runs; i++) {
        fluid->testSi(coeffs, xin, Tin + i / runs);
    }
    end = getTimeValue();
    time = (end - start) * 1e3 / runs;
    results["2D std sim simple"] = time;
    fluid->setDebug(print);
    fluid->testSi(coeffs, xin, Tin);

    fluid->setDebug(false);
    start = getTimeValue();
    for (int i = 0; i < runs; i++) {
        fluid->testHo(coeffs, xin, Tin + i / runs);
    }
    end = getTimeValue();
    time = (end - start) * 1e3 / runs;
    results["2D std sim Horner"] = time;
    fluid->setDebug(print);
    fluid->testHo(coeffs, xin, Tin);

    // Testing the 2D integrators
    fluid->setDebug(false);
    start = getTimeValue();
    for (int i = 0; i < runs; i++) {
        fluid->testSiInt(coeffs, xin, Tin + i / runs);
    }
    end = getTimeValue();
    time = (end - start) * 1e3 / runs;
    results["2D int ind simple"] = time;
    fluid->setDebug(print);
    fluid->testSiInt(coeffs, xin, Tin);

    fluid->setDebug(false);
    start = getTimeValue();
    for (int i = 0; i < runs; i++) {
        fluid->testHoInt(coeffs, xin, Tin + i / runs);
    }
    end = getTimeValue();
    time = (end - start) * 1e3 / runs;
    results["2D int ind Horner"] = time;
    fluid->setDebug(print);
    fluid->testHoInt(coeffs, xin, Tin);

    fluid->setDebug(false);
    start = getTimeValue();
    for (int i = 0; i < runs; i++) {
        fluid->test2sInt(coeffs, xin, Tin + i / runs);
    }
    end = getTimeValue();
    time = (end - start) * 1e3 / runs;
    results["2D int ind 2Steps"] = time;
    fluid->setDebug(print);
    fluid->test2sInt(coeffs, xin, Tin);

    // Testing the 2D definite integrators
    fluid->setDebug(false);
    start = getTimeValue();
    for (int i = 0; i < runs; i++) {
        fluid->testSiInt(coeffs, xin, Tin + i / runs, T0);
    }
    end = getTimeValue();
    time = (end - start) * 1e3 / runs;
    results["2D int def simple"] = time;
    fluid->setDebug(print);
    fluid->testSiInt(coeffs, xin, Tin, T0);

    fluid->setDebug(false);
    start = getTimeValue();
    for (int i = 0; i < runs; i++) {
        fluid->testHoInt(coeffs, xin, Tin + i / runs, T0);
    }
    end = getTimeValue();
    time = (end - start) * 1e3 / runs;
    results["2D int def Horner"] = time;
    fluid->setDebug(print);
    fluid->testHoInt(coeffs, xin, Tin, T0);

    fluid->setDebug(false);
    start = getTimeValue();
    for (int i = 0; i < runs; i++) {
        fluid->test2sInt(coeffs, xin, Tin + i / runs, T0);
    }
    end = getTimeValue();
    time = (end - start) * 1e3 / runs;
    results["2D int def 2Steps"] = time;
    fluid->setDebug(print);
    fluid->test2sInt(coeffs, xin, Tin, T0);

    // Testing the 2D fraction integrators
    fluid->setDebug(false);
    start = getTimeValue();
    for (int i = 0; i < runs; i++) {
        fluid->testSiFra(coeffs, xin, Tin + i / runs);
    }
    end = getTimeValue();
    time = (end - start) * 1e3 / runs;
    results["2D fra ind simple"] = time;
    fluid->setDebug(print);
    fluid->testSiFra(coeffs, xin, Tin);

    fluid->setDebug(false);
    start = getTimeValue();
    for (int i = 0; i < runs; i++) {
        fluid->testHoFra(coeffs, xin, Tin + i / runs);
    }
    end = getTimeValue();
    time = (end - start) * 1e3 / runs;
    results["2D fra ind Horner"] = time;
    fluid->setDebug(print);
    fluid->testHoFra(coeffs, xin, Tin);

    fluid->setDebug(false);
    start = getTimeValue();
    for (int i = 0; i < runs; i++) {
        fluid->test2sFra(coeffs, xin, Tin + i / runs);
    }
    end = getTimeValue();
    time = (end - start) * 1e3 / runs;
    results["2D fra ind 2Steps"] = time;
    fluid->setDebug(print);
    fluid->test2sFra(coeffs, xin, Tin);

    // Testing the 2D definite fraction integrators
    fluid->setDebug(false);
    start = getTimeValue();
    for (int i = 0; i < runs; i++) {
        fluid->testSiFra(coeffs, xin, Tin + i / runs, T0);
    }
    end = getTimeValue();
    time = (end - start) * 1e3 / runs;
    results["2D fra def simple"] = time;
    fluid->setDebug(print);
    fluid->testSiFra(coeffs, xin, Tin, T0);

    fluid->setDebug(false);
    start = getTimeValue();
    for (int i = 0; i < runs; i++) {
        fluid->testHoFra(coeffs, xin, Tin + i / runs, T0);
    }
    end = getTimeValue();
    time = (end - start) * 1e3 / runs;
    results["2D fra def Horner"] = time;
    fluid->setDebug(print);
    fluid->testHoFra(coeffs, xin, Tin, T0);

    fluid->setDebug(false);
    start = getTimeValue();
    for (int i = 0; i < runs; i++) {
        fluid->test2sFra(coeffs, xin, Tin + i / runs, T0);
    }
    end = getTimeValue();
    time = (end - start) * 1e3 / runs;
    results["2D fra def 2Steps"] = time;
    fluid->setDebug(print);
    fluid->test2sFra(coeffs, xin, Tin, T0);

    return results;
}

int main(int argc, const char* argv[]) {

    std::vector<double> tmpVector;
    std::vector<std::vector<double>> coeffs;
    int exponent = 6;

    IncompressibleTest* fluid = new IncompressibleTest();

    tmpVector.clear();
    tmpVector.push_back(1081.6353100);
    tmpVector.push_back(-2.4559523700);
    tmpVector.push_back(0.0058152057);
    tmpVector.push_back(-7.500013E-05);
    tmpVector.push_back(-7.575759E-07);
    tmpVector.push_back(1.666671E-07);
    tmpVector.push_back(-5.6609963900);
    tmpVector.push_back(0.1002726190);
    tmpVector.push_back(-0.0004797330);
    tmpVector.push_back(1.333333E-06);
    tmpVector.push_back(3.636364E-08);
    tmpVector.push_back(-0.0852857143);
    tmpVector.push_back(0.0007904762);
    tmpVector.push_back(1.428571E-06);
    tmpVector.push_back(6.666668E-07);
    tmpVector.push_back(-0.0037650794);
    tmpVector.push_back(3.333333E-05);
    tmpVector.push_back(6.984127E-07);
    coeffs.clear();
    coeffs = makeMatrix(tmpVector);

    std::map<std::string, double> time = testObject(fluid, coeffs, exponent, true);

    char buff[100];
    std::map<std::string, double>::iterator iter;
    for (iter = time.begin(); iter != time.end(); ++iter) {
        sprintf(buff, "%8.3f", iter->second);
        std::cout << "Time consumption in " << iter->first << ": " << buff << " ns per call for 1e" << exponent << " calls." << std::endl;
    }

    //std::cout << "Time consumption liquid: "   << time_liq << " µs per call from 1e" << exponent << " calls." << std::endl;

    //std::cout << "Time consumption liquid: "   << time_liq << " µs per call from 1e" << exponent << " calls." << std::endl;
    //std::cout << "Time consumption solution: " << time_sol << " µs per call from 1e" << exponent << " calls." << std::endl;

    //	SecCoolSolution* obj = new MethanolSolution();
    //	double x = 0.25;
    //	double T = 5.0 + 273.15;
    //	double p = 3e5;
    //
    //	obj->testInputs(T + 00, p, x);
    //	obj->testInputs(T + 05, p, x);
    //	obj->testInputs(T + 10, p, x);
    //	obj->testInputs(T + 15, p, x);
}

//double result = coefficients[0] * log(T);
//if (coefficients.size() > 1) {
//	for (unsigned int i=1; i<coefficients.size(); i++){
//		result += 1/(i) * coefficients[i] * pow(T,(int)(i));
//	}
//}
//
//
//double result = coefficients[0] * log(T);
//if (coefficients.size() > 1) {
//	std::vector<double> newCoeffs(coefficients.begin() + 1, coefficients.end());
//	result += polyint(newCoeffs,T);
//}
