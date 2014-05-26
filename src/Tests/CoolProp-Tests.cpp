

#include "AbstractState.h"
#include "DataStructures.h"

// ############################################
//                      TESTS
// ############################################

#if defined(ENABLE_CATCH)

#include "catch.hpp"

namespace TransportValidation{

// A structure to hold the values for one validation call
struct vel
{
public:
    std::string in1, in2, out, fluid;
    double v1, v2, tol, expected;
    vel(std::string fluid, std::string in1, double v1, std::string in2, double v2, std::string out, double expected, double tol)
    {
        this->in1 = in1; this->in2 = in2; this->fluid = fluid;
        this->v1 = v1; this->v2 = v2; this->expected = expected; 
        this->tol = tol;
    };
};

vel viscosity_validation_data[] = {
// From Vogel, JPCRD, 1998
vel("Propane", "T", 90, "Dmolar", 16.52e3, "V", 7388e-6, 1e-3),
vel("Propane", "T", 150, "Dmolar", 15.14e3, "V", 656.9e-6, 5e-3),
vel("Propane", "T", 600, "Dmolar", 10.03e3, "V", 73.92e-6, 5e-3),
vel("Propane", "T", 280, "Dmolar", 11.78e3, "V", 117.4e-6,1e-3),

// Huber, FPE, 2004
vel("n-Octane", "T", 300, "Dmolar", 6177.2, "V", 553.60e-6, 1e-3),
vel("n-Nonane", "T", 300, "Dmolar", 5619.1, "V", 709.53e-6, 1e-3),
vel("n-Decane", "T", 300, "Dmolar", 5150.4, "V", 926.44e-6, 1e-3),

// Huber, Energy & Fuels, 2004
vel("n-Dodecane", "T", 300, "Dmolar", 4411.5, "V", 1484.8e-6, 1e-3),
vel("n-Dodecane", "T", 500, "Dmolar", 3444.7, "V", 183.76e-6, 1e-3),

// Huber, I&ECR, 2006
vel("R125", "T", 300, "Dmolar", 10596.9998, "V", 177.37e-6, 1e-3),
vel("R125", "T", 400, "Dmolar", 30.631, "V", 17.070e-6, 1e-3),

// From REFPROP 9.1 since Huber I&ECR 2003 does not provide validation data
vel("R134a", "T", 185, "Q", 0, "V", 0.0012698376398294414, 1e-3),
vel("R134a", "T", 185, "Q", 1, "V", 7.4290821400170869e-006, 1e-3),
vel("R134a", "T", 360, "Q", 0, "V", 7.8146319978982133e-005, 1e-3),
vel("R134a", "T", 360, "Q", 1, "V", 1.7140264998576107e-005, 1e-3),

// From REFPROP 9.1 since Kiselev, IECR, 2005 does not provide validation data
vel("Ethanol", "T", 300, "Q", 0, "V", 0.0010439017679191723, 1e-3),
vel("Ethanol", "T", 300, "Q", 1, "V", 8.8293820936046416e-006, 1e-3),
vel("Ethanol", "T", 500, "Q", 0, "V", 6.0979347125450671e-005, 1e-3),
vel("Ethanol", "T", 500, "Q", 1, "V", 1.7229157141572511e-005, 1e-3),

// From CoolProp v5 implementation of correlation - more or less agrees with REFPROP
// Errata in BibTeX File
vel("Hydrogen", "T", 35, "Dmass", 100, "V", 5.6658815444656859e-005, 1e-3),

// From Meng 2012 experimental data (note erratum in BibTeX file)
vel("DimethylEther", "T", 253.146, "Dmass", 734.28, "V", 0.20444e-3, 3e-3),
vel("DimethylEther", "T", 373.132, "Dmass", 613.78, "V", 0.09991e-3, 3e-3),

// From Fenghour, JPCRD, 1995
vel("Ammonia", "T", 200, "Dmolar", 3.9, "V", 6.95e-6, 1e-3),
vel("Ammonia", "T", 200, "Dmolar", 42754.4, "V", 507.28e-6, 1e-3),
vel("Ammonia", "T", 398, "Dmolar", 7044.7, "V", 17.67e-6, 1e-3),
vel("Ammonia", "T", 398, "Dmolar", 21066.7, "V", 43.95e-6, 1e-3),

// From Lemmon and Jacobsen, JPCRD, 2004
vel("Nitrogen", "T", 100, "Dmolar", 1e-14, "V", 6.90349e-6, 1e-3),
vel("Nitrogen", "T", 300, "Dmolar", 1e-14, "V", 17.8771e-6, 1e-3),
vel("Nitrogen", "T", 100, "Dmolar", 25000, "V", 79.7418e-6, 1e-3),
vel("Nitrogen", "T", 200, "Dmolar", 10000, "V", 21.0810e-6, 1e-3),
vel("Nitrogen", "T", 300, "Dmolar", 5000, "V", 20.7430e-6, 1e-3),
vel("Nitrogen", "T", 126.195, "Dmolar", 11180, "V", 18.2978e-6, 1e-3),
vel("Argon", "T", 100, "Dmolar", 1e-14, "V", 8.18940e-6, 1e-3),
vel("Argon", "T", 300, "Dmolar", 1e-14, "V", 22.7241e-6, 1e-3),
vel("Argon", "T", 100, "Dmolar", 33000, "V", 184.232e-6, 1e-3),
vel("Argon", "T", 200, "Dmolar", 10000, "V", 25.5662e-6, 1e-3),
vel("Argon", "T", 300, "Dmolar", 5000, "V", 26.3706e-6, 1e-3),
vel("Argon", "T", 150.69, "Dmolar", 13400, "V", 27.6101e-6, 1e-3),
vel("Oxygen", "T", 100, "Dmolar", 1e-14, "V", 7.70243e-6, 1e-3),
vel("Oxygen", "T", 300, "Dmolar", 1e-14, "V", 20.6307e-6, 1e-3),
vel("Oxygen", "T", 100, "Dmolar", 35000, "V", 172.136e-6, 1e-3),
vel("Oxygen", "T", 200, "Dmolar", 10000, "V", 22.4445e-6, 1e-3),
vel("Oxygen", "T", 300, "Dmolar", 5000, "V", 23.7577e-6, 1e-3),
vel("Oxygen", "T", 154.6, "Dmolar", 13600, "V", 24.7898e-6, 1e-3),
vel("Air", "T", 100, "Dmolar", 1e-14, "V", 7.09559e-6, 1e-3),
vel("Air", "T", 300, "Dmolar", 1e-14, "V", 18.5230e-6, 1e-3),
vel("Air", "T", 100, "Dmolar", 28000, "V", 107.923e-6, 1e-3),
vel("Air", "T", 200, "Dmolar", 10000, "V", 21.1392e-6, 1e-3),
vel("Air", "T", 300, "Dmolar", 5000, "V", 21.3241e-6, 1e-3),
vel("Air", "T", 132.64, "Dmolar", 10400, "V", 17.7623e-6, 1e-3),

// From Michailidou, JPCRD, 2013
vel("Hexane", "T", 250, "Dmass", 1e-14, "V", 5.2584e-6, 1e-3),
vel("Hexane", "T", 400, "Dmass", 1e-14, "V", 8.4149e-6, 1e-3),
vel("Hexane", "T", 550, "Dmass", 1e-14, "V", 11.442e-6, 1e-3),
vel("Hexane", "T", 250, "Dmass", 700, "V", 528.2e-6, 1e-3),
vel("Hexane", "T", 400, "Dmass", 600, "V", 177.62e-6, 1e-3),
vel("Hexane", "T", 550, "Dmass", 500, "V", 95.002e-6, 1e-3),

// From Fenghour, JPCRD, 1998
vel("CO2", "T", 220, "Dmass", 2.440, "V", 11.06e-6, 1e-3),
vel("CO2", "T", 300, "Dmass", 1.773, "V", 15.02e-6, 1e-3),
vel("CO2", "T", 800, "Dmass", 0.662, "V", 35.09e-6, 1e-3),
vel("CO2", "T", 304, "Dmass", 254.3205, "V", 20.99e-6, 1e-3),
vel("CO2", "T", 220, "Dmass", 1194.86, "V", 269.37e-6, 1e-3),
vel("CO2", "T", 300, "Dmass", 1029.27, "V", 132.55e-6, 1e-3),
vel("CO2", "T", 800, "Dmass", 407.828, "V", 48.74e-6, 1e-3),

// Tanaka, IJT, 1996
vel("R123", "T", 265, "Dmass", 1545.8, "V", 627.1e-6, 1e-3),
vel("R123", "T", 265, "Dmass", 1.614, "V", 9.534e-6, 1e-3),
vel("R123", "T", 415, "Dmass", 1079.4, "V", 121.3e-6, 1e-3),
vel("R123", "T", 415, "Dmass", 118.9, "V", 15.82e-6, 1e-3),
//
// Krauss, IJT, 1996
vel("R152A", "T", 242, "Dmass", 1025.5, "V", 347.3e-6, 1e-3),
vel("R152A", "T", 242, "Dmass", 2.4868, "V", 8.174e-6, 1e-3),
vel("R152A", "T", 384, "Dmass", 504.51, "V", 43.29e-6, 5e-3),
vel("R152A", "T", 384, "Dmass", 239.35, "V", 21.01e-6, 10e-3),

// Huber, JPCRD, 2008 and IAPWS
vel("Water", "T", 298.15, "Dmass", 998, "V", 889.735100e-6, 1e-3),
vel("Water", "T", 298.15, "Dmass", 1200, "V", 1437.649467e-6, 1e-3),
vel("Water", "T", 373.15, "Dmass", 1000, "V", 307.883622e-6, 1e-3),
vel("Water", "T", 433.15, "Dmass", 1, "V", 14.538324e-6, 1e-3),
vel("Water", "T", 433.15, "Dmass", 1000, "V", 217.685358e-6, 1e-3),
vel("Water", "T", 873.15, "Dmass", 1, "V", 32.619287e-6, 1e-3),
vel("Water", "T", 873.15, "Dmass", 100, "V", 35.802262e-6, 1e-3),
vel("Water", "T", 873.15, "Dmass", 600, "V", 77.430195e-6, 1e-3),
vel("Water", "T", 1173.15, "Dmass", 1, "V", 44.217245e-6, 1e-3),
vel("Water", "T", 1173.15, "Dmass", 100, "V", 47.640433e-6, 1e-3),
vel("Water", "T", 1173.15, "Dmass", 400, "V", 64.154608e-6, 1e-3),

// Quinones-Cisneros, JPCRD, 2012
vel("SF6", "T", 300, "Dmass", 1e-14, "V", 15.2887e-6, 1e-4),
vel("SF6", "T", 300, "Dmass", 5.92, "V", 15.3043e-6, 1e-4),
vel("SF6", "T", 300, "Dmass", 1345.1, "V", 117.417e-6, 1e-4),
vel("SF6", "T", 400, "Dmass", 1e-14, "V", 19.6796e-6, 1e-4),
vel("SF6", "T", 400, "Dmass", 278.47, "V", 24.4272e-6, 1e-4),
vel("SF6", "T", 400, "Dmass", 1123.8, "V", 84.7835e-6, 1e-4),

// Vogel, JCED, 2011
vel("H2S", "T", 405.76, "Dmass", 0.740, "V", 15.56e-6, 1e-2),
// Quinones-Cisneros, JCED, 2012
vel("H2S", "T", 373.15, "P", 1e6, "V", 15.9e-6, 1e-3),
vel("H2S", "T", 373.15, "P", 10e6, "V", 48.3e-6, 1e-3),
vel("H2S", "T", 423.15, "P", 1e6, "V", 16.5e-6, 1e-3),
vel("H2S", "T", 423.15, "P", 10e6, "V", 21.2e-6, 1e-3),

// Geller, Purdue Conference, 2000
vel("R410A", "T", 243.15, "Q", 0, "V", 238.61e-6, 5e-2),
vel("R410A", "T", 243.15, "Q", 1, "V", 10.37e-6, 5e-2),
vel("R410A", "T", 333.15, "Q", 0, "V", 70.71e-6, 5e-2),
vel("R410A", "T", 333.15, "Q", 1, "V", 19.19e-6, 5e-2),
vel("R407C", "T", 243.15, "Q", 0, "V", 304.18e-6, 1e-2),
vel("R407C", "T", 243.15, "Q", 1, "V", 9.83e-6, 1e-2),
vel("R407C", "T", 333.15, "Q", 0, "V", 95.96e-6, 1e-2),
vel("R407C", "T", 333.15, "Q", 1, "V", 16.38e-6, 1e-2),
vel("R404A", "T", 243.15, "Q", 0, "V", 264.67e-6, 1e-2),
vel("R404A", "T", 243.15, "Q", 1, "V", 10.13e-6, 1e-2),
vel("R404A", "T", 333.15, "Q", 0, "V", 73.92e-6, 1e-2),
vel("R404A", "T", 333.15, "Q", 1, "V", 18.56e-6, 1e-2),
vel("R507A", "T", 243.15, "Q", 0, "V", 284.59e-6, 3e-2),
vel("R507A", "T", 243.15, "Q", 1, "V", 9.83e-6, 1e-2),
vel("R507A", "T", 333.15, "Q", 0, "V", 74.37e-6, 1e-2),
vel("R507A", "T", 333.15, "Q", 1, "V", 19.35e-6, 1e-2),

// From Arp, NIST, 1998
vel("Helium", "T", 3.6, "P", 0.180e6, "V", 3.745e-6, 1e-2),
vel("Helium", "T", 50, "P", 0.180e6, "V", 6.376e-6, 1e-2),
vel("Helium", "T", 400, "P", 0.180e6, "V", 24.29e-6, 1e-2),

// From Shan, ASHRAE, 2000
vel("R23", "T", 180, "Dmolar", 21097, "V", 353.88e-6, 1e-4),
vel("R23", "T", 420, "Dmolar", 7564, "V", 39.459e-6, 1e-4),
vel("R23", "T", 370, "Dmolar", 32.62, "V", 18.213e-6, 1e-4),

// From Friend, JPCRD, 1991
vel("Ethane", "T", 100, "Dmolar", 21330, "V", 878.6e-6, 1e-2),
vel("Ethane", "T", 430, "Dmolar", 12780, "V", 58.70e-6, 1e-2),
vel("Ethane", "T", 500, "Dmolar", 11210, "V", 48.34e-6, 1e-2),

// From REFPROP 9.1 since no data provided
vel("n-Butane", "T", 150, "Q", 0, "V", 0.0013697657668, 1e-4),
vel("n-Butane", "T", 400, "Q", 1, "V", 1.2027464524762453e-005, 1e-4),
vel("IsoButane", "T", 120, "Q", 0, "V", 0.0060558450757844271, 1e-4),
vel("IsoButane", "T", 400, "Q", 1, "V", 1.4761041187617117e-005, 2e-4),
vel("R134a", "T", 175, "Q", 0, "V", 0.0017558494524138289, 1e-4),
vel("R134a", "T", 360, "Q", 1, "V", 1.7140264998576107e-005, 1e-4),

};

class TransportValidationFixture
{
protected:
    long double actual, x1, x2;
    CoolProp::AbstractState *pState;
    int pair;
public:
    TransportValidationFixture(){ pState = NULL; }
    ~TransportValidationFixture(){ delete pState; }
    void set_backend(std::string backend, std::string fluid_name){
        pState = CoolProp::AbstractState::factory(backend, fluid_name);
    }
    void set_pair(std::string &in1, double v1, std::string &in2, double v2){ 
        double o1, o2;
        long iin1 = CoolProp::get_parameter_index(in1);
        long iin2 = CoolProp::get_parameter_index(in2);
        long pair = CoolProp::generate_update_pair(iin1, v1, iin2, v2, o1, o2);
        pState->update(pair, o1, o2);
    }
    void get_value(long key)
    {
        actual = pState->keyed_output(key);
    }
};

TEST_CASE_METHOD(TransportValidationFixture, "Compare viscosities against published data", "[viscosity]")
{
    int inputsN = sizeof(viscosity_validation_data)/sizeof(viscosity_validation_data[0]);
    for (int i = 0; i < inputsN; ++i)
    {
        vel el = viscosity_validation_data[i];
        CHECK_NOTHROW(set_backend("HEOS", el.fluid));

        CAPTURE(el.fluid);
        CAPTURE(el.in1);
        CAPTURE(el.v1);
        CAPTURE(el.in2);
        CAPTURE(el.v2);
        CHECK_NOTHROW(set_pair(el.in1, el.v1, el.in2, el.v2));
        get_value(CoolProp::iviscosity);
        CAPTURE(el.expected);
        CAPTURE(actual);
        CHECK(fabs(actual/el.expected-1) < el.tol);
    }
}

vel conductivity_validation_data[] = {
///\todo Re-enable the conductivity tests that fail due to not having viscosity correlation

// From Assael, JPCRD, 2013
vel("Hexane", "T", 250, "Dmass", 700, "L", 137.62e-3, 1e-4),
vel("Hexane", "T", 400, "Dmass", 2, "L", 23.558e-3, 1e-4),
vel("Hexane", "T", 400, "Dmass", 650, "L", 129.28e-3, 2e-4),
vel("Hexane", "T", 510, "Dmass", 2, "L", 36.772e-3, 1e-4),

// From Assael, JPCRD, 2013
//vel("Heptane", "T", 250, "Dmass", 720, "L", 137.09e-3, 1e-4),
//vel("Heptane", "T", 400, "Dmass", 2, "L", 21.794e-3, 1e-4),
//vel("Heptane", "T", 400, "Dmass", 650, "L", 120.75e-3, 1e-4),
//vel("Heptane", "T", 535, "Dmass", 100, "L", 51.655e-3, 1e-4),

// From Assael, JPCRD, 2013
vel("Ethanol", "T", 300, "Dmass", 850, "L", 209.68e-3, 1e-4),
vel("Ethanol", "T", 400, "Dmass", 2, "L", 26.108e-3, 1e-4),
vel("Ethanol", "T", 400, "Dmass", 690, "L", 149.21e-3, 1e-4),
vel("Ethanol", "T", 500, "Dmass", 10, "L", 39.594e-3, 1e-4),

//// From Assael, JPCRD, 2012
//vel("Toluene", "T", 298.15, "Dmass", 1e-15, "L", 10.749e-3, 1e-4),
//vel("Toluene", "T", 298.15, "Dmass", 862.948, "L", 130.66e-3, 1e-4),
//vel("Toluene", "T", 298.15, "Dmass", 876.804, "L", 136.70e-3, 1e-4),
//vel("Toluene", "T", 595, "Dmass", 1e-15, "L", 40.538e-3, 1e-4),
//vel("Toluene", "T", 595, "Dmass", 46.512, "L", 41.549e-3, 1e-4),
//vel("Toluene", "T", 185, "Dmass", 1e-15, "L", 4.3758e-3, 1e-4),
//vel("Toluene", "T", 185, "Dmass", 968.821, "L", 158.24e-3, 1e-4),

// From Assael, JPCRD, 2012
vel("SF6", "T", 298.15, "Dmass", 1e-13, "L", 12.952e-3, 1e-4),
vel("SF6", "T", 298.15, "Dmass", 100, "L", 14.127e-3, 1e-4),
vel("SF6", "T", 298.15, "Dmass", 1600, "L", 70.747e-3, 1e-4),
vel("SF6", "T", 310, "Dmass", 1e-13, "L", 13.834e-3, 1e-4),
vel("SF6", "T", 310, "Dmass", 1200, "L", 49.173e-3, 1e-4),
vel("SF6", "T", 480, "Dmass", 100, "L", 28.863e-3, 1e-4),

//// From Assael, JPCRD, 2012
//vel("Benzene", "T", 290, "Dmass", 890, "L", 147.66e-3, 1e-4),
//vel("Benzene", "T", 500, "Dmass", 2, "L", 30.174e-3, 1e-4),
//vel("Benzene", "T", 500, "Dmass", 32, "L", 32.175e-3, 1e-4),
//vel("Benzene", "T", 500, "Dmass", 800, "L", 141.24e-3, 1e-4),
//vel("Benzene", "T", 575, "Dmass", 1.7, "L", 37.763e-3, 1e-4),

// From Assael, JPCRD, 2011
vel("Hydrogen", "T", 298.15, "Dmass", 1e-13, "L", 185.67e-3, 1e-4),
vel("Hydrogen", "T", 298.15, "Dmass", 0.80844, "L", 186.97e-3, 1e-4),
vel("Hydrogen", "T", 298.15, "Dmass", 14.4813, "L", 201.35e-3, 1e-4),
vel("Hydrogen", "T", 35, "Dmass", 1e-13, "L", 26.988e-3, 1e-4),
vel("Hydrogen", "T", 35, "Dmass", 30, "L", 75.594e-3, 1e-4),
vel("Hydrogen", "T", 18, "Dmass", 1e-13, "L", 13.875e-3, 1e-4),
vel("Hydrogen", "T", 18, "Dmass", 75, "L", 104.48e-3, 1e-4),
/*vel("ParaHydrogen", "T", 298.15, "Dmass", 1e-13, "L", 192.38e-3, 1e-4),
vel("ParaHydrogen", "T", 298.15, "Dmass", 0.80844, "L", 192.81e-3, 1e-4),
vel("ParaHydrogen", "T", 298.15, "Dmass", 14.4813, "L", 207.85e-3, 1e-4),
vel("ParaHydrogen", "T", 35, "Dmass", 1e-13, "L", 27.222e-3, 1e-4),
vel("ParaHydrogen", "T", 35, "Dmass", 30, "L", 70.335e-3, 1e-4),
vel("ParaHydrogen", "T", 18, "Dmass", 1e-13, "L", 13.643e-3, 1e-4),
vel("ParaHydrogen", "T", 18, "Dmass", 75, "L", 100.52e-3, 1e-4),*/


// Some of these don't work
vel("R125", "T", 341, "Dmass", 600, "L", 0.056564297849487, 2e-2),
vel("R125", "T", 200, "Dmass", 1e-13, "L", 0.007036843623086, 2e-2),
vel("IsoButane", "T", 390, "Dmass", 387.09520158645068, "L", 0.063039, 1e-2),
vel("IsoButane", "T", 390, "Dmass", 85.76703973869482, "L", 0.036603, 2e-2),
vel("n-Butane", "T", 415, "Dmass", 360.01895129934866, "L", 0.067045, 1e-2),
vel("n-Butane", "T", 415, "Dmass", 110.3113177144, "L", 0.044449, 1e-2),

// From Huber, FPE, 2005
vel("n-Octane", "T", 300, "Dmolar", 6177.2, "L", 0.12836, 1e-4),
vel("n-Nonane", "T", 300, "Dmolar", 5619.4, "L", 0.13031, 1e-4),
//vel("n-Decane", "T", 300, "Dmass", 5150.4, "L", 0.13280, 1e-4), // no viscosity

// From Huber, EF, 2004
vel("n-Dodecane", "T", 300, "Dmolar", 4411.5, "L", 0.13829, 1e-4),
vel("n-Dodecane", "T", 500, "Dmolar", 3444.7, "L", 0.09384, 1e-4),
vel("n-Dodecane", "T", 660, "Dmolar", 1500.98, "L", 0.090346, 1e-4),

// From REFPROP 9.1 since no data provided in Marsh, 2002
vel("n-Propane", "T", 368, "Q", 0, "L", 0.07282154952457,1e-3),
vel("n-Propane", "T", 368, "Dmolar", 1e-10, "L", 0.0266135388745317,1e-4),

// From Perkins, JCED, 2011
//vel("R1234yf", "T", 250, "Dmass", 2.80006, "L", 0.0098481, 1e-4),
//vel("R1234yf", "T", 300, "Dmass", 4.671556, "L", 0.013996, 1e-4),
//vel("R1234yf", "T", 250, "Dmass", 1299.50, "L", 0.088574, 1e-4),
//vel("R1234yf", "T", 300, "Dmass", 1182.05, "L", 0.075245, 1e-4),
//vel("R1234ze(E)", "T", 250, "Dmass", 2.80451, "L", 0.0098503, 1e-4),
//vel("R1234ze(E)", "T", 300, "Dmass", 4.67948, "L", 0.013933, 1e-4),
//vel("R1234ze(E)", "T", 250, "Dmass", 1349.37, "L", 0.10066, 1e-4),
//vel("R1234ze(E)", "T", 300, "Dmass", 1233.82, "L", 0.085389, 1e-4),

// From Laesecke, IJR 1995
vel("R123", "T", 180, "Dmass", 1739, "L", 110.9e-3, 1e-3),
vel("R123", "T", 180, "Dmass", 0.2873e-2, "L", 2.473e-3, 1e-3),
vel("R123", "T", 430, "Dmass", 996.35, "L", 45.62e-3, 1e-3),
vel("R123", "T", 430, "Dmass", 166.9,  "L", 21.03e-3, 1e-3),

// From Vesovic, JPCRD, 1990
vel("CO2", "T", 220, "Dmass", 2.440, "L", 10.90e-3, 1e-3),
vel("CO2", "T", 300, "Dmass", 1.773, "L", 16.77e-3, 1e-3),
vel("CO2", "T", 800, "Dmass", 0.662, "L", 56.65e-3, 1e-3),
vel("CO2", "T", 304, "Dmass", 254.3205, "L", 42.52e-3, 1e-3),
vel("CO2", "T", 220, "Dmass", 1194.86, "L", 187.50e-3, 1e-3),
vel("CO2", "T", 300, "Dmass", 1029.27, "L", 137.61e-3, 1e-3),
vel("CO2", "T", 800, "Dmass", 407.828, "L", 78.47e-3, 1e-3),

};

TEST_CASE_METHOD(TransportValidationFixture, "Compare thermal conductivities against published data", "[conductivity]")
{
    int inputsN = sizeof(conductivity_validation_data)/sizeof(conductivity_validation_data[0]);
    for (int i = 0; i < inputsN; ++i)
    {
        vel el = conductivity_validation_data[i];
        CHECK_NOTHROW(set_backend("HEOS", el.fluid));
        CAPTURE(el.fluid);
        CAPTURE(el.in1);
        CAPTURE(el.v1);
        CAPTURE(el.in2);
        CAPTURE(el.v2);
        CHECK_NOTHROW(set_pair(el.in1, el.v1, el.in2, el.v2));
        get_value(CoolProp::iconductivity);
        CAPTURE(el.expected);
        CAPTURE(actual);
        CHECK(fabs(actual/el.expected-1) < el.tol);
    }
}


}; /* namespace TransportValidation */

static int inputs[] = {
    CoolProp::DmolarT_INPUTS,
    //CoolProp::SmolarT_INPUTS,
    //CoolProp::HmolarT_INPUTS, 
    //CoolProp::TUmolar_INPUTS,

    CoolProp::DmolarP_INPUTS, 
    CoolProp::DmolarHmolar_INPUTS, 
    CoolProp::DmolarSmolar_INPUTS, 
    CoolProp::DmolarUmolar_INPUTS,
        
    /*
    CoolProp::HmolarP_INPUTS,
    CoolProp::PSmolar_INPUTS,
    CoolProp::PUmolar_INPUTS, 
    */

    /*
    CoolProp::HmolarSmolar_INPUTS, 
    CoolProp::HmolarUmolar_INPUTS, 
    CoolProp::SmolarUmolar_INPUTS 
    */
};

class ConsistencyFixture
{
protected:
    long double hmolar, pmolar, smolar, umolar, rhomolar, T, p, x1, x2;
    CoolProp::AbstractState *pState;
    int pair;
public:
    ConsistencyFixture(){
        pState = NULL;
    }
    ~ConsistencyFixture(){
        delete pState;
    }
    void set_backend(std::string backend, std::string fluid_name){
        pState = CoolProp::AbstractState::factory(backend, fluid_name);
    }
    void set_pair(int pair){ 
        this->pair = pair;
    }
    void set_TP(long double T, long double p)
    {
        this->T = T; this->p = p;
        CoolProp::AbstractState &State = *pState;

        // Start with T,P as inputs, cycle through all the other pairs that are supported
        State.update(CoolProp::PT_INPUTS, p, T);
            
        // Set the other state variables
        rhomolar = State.rhomolar(); hmolar = State.hmolar(); smolar = State.smolar(); umolar = State.umolar();
    }
    void get_variables()
    {
        CoolProp::AbstractState &State = *pState;
            
        switch (pair)
        {
        /// In this group, T is one of the known inputs, iterate for the other one (easy)
        case CoolProp::HmolarT_INPUTS:
            x1 = hmolar; x2 = T;  break;
        case CoolProp::SmolarT_INPUTS:
            x1 = smolar; x2 = T; break;
        case CoolProp::TUmolar_INPUTS:
            x1 = T; x2 = umolar; break;
        case CoolProp::DmolarT_INPUTS:
            x1 = rhomolar; x2 = T; break;

        /// In this group, D is one of the known inputs, iterate for the other one (a little bit harder)
        case CoolProp::DmolarHmolar_INPUTS:
            x1 = rhomolar; x2 = hmolar; break;
        case CoolProp::DmolarSmolar_INPUTS:
            x1 = rhomolar; x2 = smolar; break;
        case CoolProp::DmolarUmolar_INPUTS:
            x1 = rhomolar; x2 = umolar; break;
        case CoolProp::DmolarP_INPUTS:
            x1 = rhomolar; x2 = p; break;

        /// In this group, p is one of the known inputs (a little less easy)
        case CoolProp::HmolarP_INPUTS:
            x1 = hmolar; x2 = p; break;
        case CoolProp::PSmolar_INPUTS:
            x1 = p; x2 = smolar; break;
        case CoolProp::PUmolar_INPUTS:
            x1 = p; x2 = umolar; break;

        case CoolProp::HmolarSmolar_INPUTS:
            x1 = hmolar; x2 = smolar; break;
        case CoolProp::SmolarUmolar_INPUTS:
            x1 = smolar; x2 = umolar; break;
        }
    }
    void single_phase_consistency_check()
    {
        CoolProp::AbstractState &State = *pState;
        State.update(pair, x1, x2);

        // Make sure we end up back at the same temperature and pressure we started out with
        if(fabs(T-State.T()) > 1e-2) throw CoolProp::ValueError(format("Error on T [%g K] is greater than 1e-2",fabs(State.T()-T)));
        if(fabs(p-State.p())/p*100 > 1e-2)  throw CoolProp::ValueError(format("Error on p [%g %%] is greater than 1e-2 %%",fabs(p-State.p())/p ));
    }
};

TEST_CASE_METHOD(ConsistencyFixture, "Test all input pairs for CO2 using all valid backends", "[]")
{
    CHECK_NOTHROW(set_backend("HEOS", "CO2"));
        
    int inputsN = sizeof(inputs)/sizeof(inputs[0]);
    for (double p = 600000; p < 800000000.0; p *= 5)
    {
        for (double T = 220; T < pState->Tmax(); T += 5)
        {
            CHECK_NOTHROW(set_TP(T, p));

            for (int i = 0; i < inputsN; ++i)
            {
                int pair = inputs[i];
                std::string pair_desc = CoolProp::get_input_pair_short_desc(pair);
                set_pair(pair);
                CAPTURE(pair_desc);
                CAPTURE(T);
                CAPTURE(p);
                get_variables();
                CAPTURE(x1);
                CAPTURE(x2);
                CHECK_NOTHROW(single_phase_consistency_check());
            }
        }
    }
}

#endif
