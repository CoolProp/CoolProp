

#include "AbstractState.h"
#include "DataStructures.h"
#include "../Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#include "../Backends/Helmholtz/HelmholtzEOSBackend.h"
// ############################################
//                      TESTS
// ############################################

#if defined(ENABLE_CATCH)

#    include "crossplatform_shared_ptr.h"
#    include <catch2/catch_all.hpp>
#    include "CoolPropTools.h"
#    include "CoolProp.h"

using namespace CoolProp;

namespace TransportValidation {

// A structure to hold the values for one validation call
struct vel
{
   public:
    std::string in1, in2, out, fluid;
    double v1, v2, tol, expected;
    vel(std::string fluid, std::string in1, double v1, std::string in2, double v2, std::string out, double expected, double tol) {
        this->in1 = in1;
        this->in2 = in2;
        this->fluid = fluid;
        this->v1 = v1;
        this->v2 = v2;
        this->expected = expected;
        this->tol = tol;
    };
};

vel viscosity_validation_data[] = {
  // From Vogel, JPCRD, 1998
  vel("Propane", "T", 90, "Dmolar", 16.52e3, "V", 7388e-6, 1e-3),
  vel("Propane", "T", 150, "Dmolar", 15.14e3, "V", 656.9e-6, 5e-3),
  vel("Propane", "T", 600, "Dmolar", 10.03e3, "V", 73.92e-6, 5e-3),
  vel("Propane", "T", 280, "Dmolar", 11.78e3, "V", 117.4e-6, 1e-3),

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
  vel("Hydrogen", "T", 35, "Dmass", 100, "V", 5.47889e-005, 1e-3),

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

  // From Assael, JPCRD, 2014
  vel("Heptane", "T", 250, "Dmass", 1e-14, "V", 4.9717e-6, 1e-3),
  vel("Heptane", "T", 400, "Dmass", 1e-14, "V", 7.8361e-6, 1e-3),
  vel("Heptane", "T", 550, "Dmass", 1e-14, "V", 10.7394e-6, 1e-3),
  vel("Heptane", "T", 250, "Dmass", 720, "V", 725.69e-6, 1e-3),
  vel("Heptane", "T", 400, "Dmass", 600, "V", 175.94e-6, 1e-3),
  vel("Heptane", "T", 550, "Dmass", 500, "V", 95.105e-6, 1e-3),

  // From Fenghour, JPCRD, 1998
  vel("CO2", "T", 220, "Dmass", 2.440, "V", 11.06e-6, 1e-3),
  vel("CO2", "T", 300, "Dmass", 1.773, "V", 15.02e-6, 1e-3),
  vel("CO2", "T", 800, "Dmass", 0.662, "V", 35.09e-6, 1e-3),
  vel("CO2", "T", 304, "Dmass", 254.320, "V", 20.99e-6, 1e-2),  // no critical enhancement
  vel("CO2", "T", 220, "Dmass", 1194.86, "V", 269.37e-6, 1e-3),
  vel("CO2", "T", 300, "Dmass", 1029.27, "V", 132.55e-6, 1e-3),
  vel("CO2", "T", 800, "Dmass", 407.828, "V", 48.74e-6, 1e-3),

  // Tanaka, IJT, 1996
  vel("R123", "T", 265, "Dmass", 1545.8, "V", 627.1e-6, 1e-3),
  vel("R123", "T", 265, "Dmass", 1.614, "V", 9.534e-6, 1e-3),
  vel("R123", "T", 415, "Dmass", 1079.4, "V", 121.3e-6, 1e-3),
  vel("R123", "T", 415, "Dmass", 118.9, "V", 15.82e-6, 1e-3),

  // Krauss, IJT, 1996
  vel("R152A", "T", 242, "Dmass", 1025.5, "V", 347.3e-6, 1e-3),
  vel("R152A", "T", 242, "Dmass", 2.4868, "V", 8.174e-6, 1e-3),
  vel("R152A", "T", 384, "Dmass", 504.51, "V", 43.29e-6, 5e-3),
  vel("R152A", "T", 384, "Dmass", 239.35, "V", 21.01e-6, 10e-3),

  // Huber, JPCRD, 2008 and IAPWS
  vel("Water", "T", 298.15, "Dmass", 998, "V", 889.735100e-6, 1e-7),
  vel("Water", "T", 298.15, "Dmass", 1200, "V", 1437.649467e-6, 1e-7),
  vel("Water", "T", 373.15, "Dmass", 1000, "V", 307.883622e-6, 1e-7),
  vel("Water", "T", 433.15, "Dmass", 1, "V", 14.538324e-6, 1e-7),
  vel("Water", "T", 433.15, "Dmass", 1000, "V", 217.685358e-6, 1e-7),
  vel("Water", "T", 873.15, "Dmass", 1, "V", 32.619287e-6, 1e-7),
  vel("Water", "T", 873.15, "Dmass", 100, "V", 35.802262e-6, 1e-7),
  vel("Water", "T", 873.15, "Dmass", 600, "V", 77.430195e-6, 1e-7),
  vel("Water", "T", 1173.15, "Dmass", 1, "V", 44.217245e-6, 1e-7),
  vel("Water", "T", 1173.15, "Dmass", 100, "V", 47.640433e-6, 1e-7),
  vel("Water", "T", 1173.15, "Dmass", 400, "V", 64.154608e-6, 1e-7),
  vel("Water", "T", 647.35, "Dmass", 122, "V", 25.520677e-6, 1e-7),
  vel("Water", "T", 647.35, "Dmass", 222, "V", 31.337589e-6, 1e-7),
  vel("Water", "T", 647.35, "Dmass", 272, "V", 36.228143e-6, 1e-7),
  vel("Water", "T", 647.35, "Dmass", 322, "V", 42.961579e-6, 1e-7),
  vel("Water", "T", 647.35, "Dmass", 372, "V", 45.688204e-6, 1e-7),
  vel("Water", "T", 647.35, "Dmass", 422, "V", 49.436256e-6, 1e-7),

  // Quinones-Cisneros, JPCRD, 2012
  vel("SF6", "T", 300, "Dmass", 1e-14, "V", 15.2887e-6, 1e-4),
  vel("SF6", "T", 300, "Dmass", 5.92, "V", 15.3043e-6, 1e-4),
  vel("SF6", "T", 300, "Dmass", 1345.1, "V", 117.417e-6, 1e-4),
  vel("SF6", "T", 400, "Dmass", 1e-14, "V", 19.6796e-6, 1e-4),
  vel("SF6", "T", 400, "Dmass", 278.47, "V", 24.4272e-6, 1e-4),
  vel("SF6", "T", 400, "Dmass", 1123.8, "V", 84.7835e-6, 1e-4),

  // Quinones-Cisneros, JCED, 2012, data from validation
  vel("H2S", "T", 200, "P", 1000e5, "V", 0.000460287, 1e-3),
  vel("H2S", "T", 200, "P", 0.251702e5, "V", 8.02322E-06, 1e-3),
  vel("H2S", "T", 596.961, "P", 1000e5, "V", 6.94741E-05, 1e-3),
  vel("H2S", "T", 596.961, "P", 1e5, "V", 2.38654E-05, 1e-3),

  // Geller, Purdue Conference, 2000
  //vel("R410A", "T", 243.15, "Q", 0, "V", 238.61e-6, 5e-2),
  //vel("R410A", "T", 243.15, "Q", 1, "V", 10.37e-6, 5e-2),
  //vel("R410A", "T", 333.15, "Q", 0, "V", 70.71e-6, 5e-2),
  //vel("R410A", "T", 333.15, "Q", 1, "V", 19.19e-6, 5e-2),
  //vel("R407C", "T", 243.15, "Q", 0, "V", 304.18e-6, 1e-2),
  //vel("R407C", "T", 243.15, "Q", 1, "V", 9.83e-6, 1e-2),
  //vel("R407C", "T", 333.15, "Q", 0, "V", 95.96e-6, 1e-2),
  //vel("R407C", "T", 333.15, "Q", 1, "V", 16.38e-6, 1e-2),
  //vel("R404A", "T", 243.15, "Q", 0, "V", 264.67e-6, 1e-2),
  //vel("R404A", "T", 243.15, "Q", 1, "V", 10.13e-6, 1e-2),
  //vel("R404A", "T", 333.15, "Q", 0, "V", 73.92e-6, 1e-2),
  //vel("R404A", "T", 333.15, "Q", 1, "V", 18.56e-6, 1e-2),
  //vel("R507A", "T", 243.15, "Q", 0, "V", 284.59e-6, 3e-2),
  //vel("R507A", "T", 243.15, "Q", 1, "V", 9.83e-6, 1e-2),
  //vel("R507A", "T", 333.15, "Q", 0, "V", 74.37e-6, 1e-2),
  //vel("R507A", "T", 333.15, "Q", 1, "V", 19.35e-6, 1e-2),

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

  // From Xiang, JPCRD, 2006
  vel("Methanol", "T", 300, "Dmass", 0.12955, "V", 0.009696e-3, 1e-3),
  vel("Methanol", "T", 300, "Dmass", 788.41, "V", 0.5422e-3, 1e-3),
  vel("Methanol", "T", 630, "Dmass", 0.061183, "V", 0.02081e-3, 1e-3),
  vel("Methanol", "T", 630, "Dmass", 888.50, "V", 0.2405e-3, 1e-1),  // They use a different EOS in the high pressure region

  // From REFPROP 9.1 since no data provided
  vel("n-Butane", "T", 150, "Q", 0, "V", 0.0013697657668, 1e-4),
  vel("n-Butane", "T", 400, "Q", 1, "V", 1.2027464524762453e-005, 1e-4),
  vel("IsoButane", "T", 120, "Q", 0, "V", 0.0060558450757844271, 1e-4),
  vel("IsoButane", "T", 400, "Q", 1, "V", 1.4761041187617117e-005, 2e-4),
  vel("R134a", "T", 175, "Q", 0, "V", 0.0017558494524138289, 1e-4),
  vel("R134a", "T", 360, "Q", 1, "V", 1.7140264998576107e-005, 1e-4),

  // From Tariq, JPCRD, 2014
  vel("Cyclohexane", "T", 300, "Dmolar", 1e-10, "V", 7.058e-6, 1e-4),
  vel("Cyclohexane", "T", 300, "Dmolar", 0.0430e3, "V", 6.977e-6, 1e-4),
  vel("Cyclohexane", "T", 300, "Dmolar", 9.1756e3, "V", 863.66e-6, 1e-4),
  vel("Cyclohexane", "T", 300, "Dmolar", 9.9508e3, "V", 2850.18e-6, 1e-4),
  vel("Cyclohexane", "T", 500, "Dmolar", 1e-10, "V", 11.189e-6, 1e-4),
  vel("Cyclohexane", "T", 500, "Dmolar", 6.0213e3, "V", 94.842e-6, 1e-4),
  vel("Cyclohexane", "T", 500, "Dmolar", 8.5915e3, "V", 380.04e-6, 1e-4),
  vel("Cyclohexane", "T", 700, "Dmolar", 1e-10, "V", 15.093e-6, 1e-4),
  vel("Cyclohexane", "T", 700, "Dmolar", 7.4765e3, "V", 176.749e-6, 1e-4),

  // From Avgeri, JPCRD, 2014
  vel("Benzene", "T", 300, "Dmass", 1e-10, "V", 7.625e-6, 1e-4),
  vel("Benzene", "T", 400, "Dmass", 1e-10, "V", 10.102e-6, 1e-4),
  vel("Benzene", "T", 550, "Dmass", 1e-10, "V", 13.790e-6, 1e-4),
  vel("Benzene", "T", 300, "Dmass", 875, "V", 608.52e-6, 1e-4),
  vel("Benzene", "T", 400, "Dmass", 760, "V", 211.74e-6, 1e-4),
  vel("Benzene", "T", 550, "Dmass", 500, "V", 60.511e-6, 1e-4),

  // From Cao, JPCRD, 2016
  vel("m-Xylene", "T", 300, "Dmolar", 1e-10, "V", 6.637e-6, 1e-4),
  vel("m-Xylene", "T", 300, "Dmolar", 0.04 * 1e3, "V", 6.564e-6, 1e-4),
  vel("m-Xylene", "T", 300, "Dmolar", 8.0849 * 1e3, "V", 569.680e-6, 1e-4),
  vel("m-Xylene", "T", 300, "Dmolar", 8.9421 * 1e3, "V", 1898.841e-6, 1e-4),
  vel("m-Xylene", "T", 400, "Dmolar", 1e-10, "V", 8.616e-6, 1e-4),
  vel("m-Xylene", "T", 400, "Dmolar", 0.04 * 1e3, "V", 8.585e-6, 1e-4),
  vel("m-Xylene", "T", 400, "Dmolar", 7.2282 * 1e3, "V", 238.785e-6, 1e-4),
  vel("m-Xylene", "T", 400, "Dmolar", 8.4734 * 1e3, "V", 718.950e-6, 1e-4),
  vel("m-Xylene", "T", 600, "Dmolar", 1e-10, "V", 12.841e-6, 1e-4),
  vel("m-Xylene", "T", 600, "Dmolar", 0.04 * 1e3, "V", 12.936e-6, 1e-4),
  vel("m-Xylene", "T", 600, "Dmolar", 7.6591 * 1e3, "V", 299.164e-6, 1e-4),

  // From Cao, JPCRD, 2016
  vel("o-Xylene", "T", 300, "Dmolar", 1e-10, "V", 6.670e-6, 1e-4),
  vel("o-Xylene", "T", 300, "Dmolar", 0.04 * 1e3, "V", 6.598e-6, 1e-4),
  vel("o-Xylene", "T", 300, "Dmolar", 8.2369 * 1e3, "V", 738.286e-6, 1e-4),
  vel("o-Xylene", "T", 300, "Dmolar", 8.7845 * 1e3, "V", 1645.436e-6, 1e-4),
  vel("o-Xylene", "T", 400, "Dmolar", 1e-10, "V", 8.658e-6, 1e-4),
  vel("o-Xylene", "T", 400, "Dmolar", 0.04 * 1e3, "V", 8.634e-6, 1e-4),
  vel("o-Xylene", "T", 400, "Dmolar", 7.4060 * 1e3, "V", 279.954e-6, 1e-4),
  vel("o-Xylene", "T", 400, "Dmolar", 8.2291 * 1e3, "V", 595.652e-6, 1e-4),
  vel("o-Xylene", "T", 600, "Dmolar", 1e-10, "V", 12.904e-6, 1e-4),
  vel("o-Xylene", "T", 600, "Dmolar", 0.04 * 1e3, "V", 13.018e-6, 1e-4),
  vel("o-Xylene", "T", 600, "Dmolar", 7.2408 * 1e3, "V", 253.530e-6, 1e-4),

  // From Balogun, JPCRD, 2016
  vel("p-Xylene", "T", 300, "Dmolar", 1e-10, "V", 6.604e-6, 1e-4),
  vel("p-Xylene", "T", 300, "Dmolar", 0.049 * 1e3, "V", 6.405e-6, 1e-4),
  vel("p-Xylene", "T", 300, "Dmolar", 8.0548 * 1e3, "V", 593.272e-6, 1e-4),
  vel("p-Xylene", "T", 300, "Dmolar", 8.6309 * 1e3, "V", 1266.337e-6, 1e-4),
  vel("p-Xylene", "T", 400, "Dmolar", 1e-10, "V", 8.573e-6, 1e-4),
  vel("p-Xylene", "T", 400, "Dmolar", 7.1995 * 1e3, "V", 239.202e-6, 1e-4),
  vel("p-Xylene", "T", 400, "Dmolar", 8.0735 * 1e3, "V", 484.512e-6, 1e-4),
  vel("p-Xylene", "T", 600, "Dmolar", 1e-10, "V", 12.777e-6, 1e-4),
  vel("p-Xylene", "T", 600, "Dmolar", 7.0985 * 1e3, "V", 209.151e-6, 1e-4),

  // From Mylona, JPCRD, 2014
  vel("EthylBenzene", "T", 617, "Dmass", 316, "V", 33.22e-6, 1e-2),

  // Heavy Water, IAPWS formulation
  vel("HeavyWater", "T", 0.5000 * 643.847, "Dmass", 3.07 * 358, "V", 12.0604912273 * 55.2651e-6, 1e-5),
  vel("HeavyWater", "T", 0.9000 * 643.847, "Dmass", 2.16 * 358, "V", 1.6561616211 * 55.2651e-6, 1e-5),
  vel("HeavyWater", "T", 1.2000 * 643.847, "Dmass", 0.8 * 358, "V", 0.7651099154 * 55.2651e-6, 1e-5),

  // Toluene, Avgeri, JPCRD, 2015
  vel("Toluene", "T", 300, "Dmass", 1e-10, "V", 7.023e-6, 1e-4),
  vel("Toluene", "T", 400, "Dmass", 1e-10, "V", 9.243e-6, 1e-4),
  vel("Toluene", "T", 550, "Dmass", 1e-10, "V", 12.607e-6, 1e-4),
  vel("Toluene", "T", 300, "Dmass", 865, "V", 566.78e-6, 1e-4),
  vel("Toluene", "T", 400, "Dmass", 770, "V", 232.75e-6, 1e-4),
  vel("Toluene", "T", 550, "Dmass", 550, "V", 80.267e-6, 1e-4),

};

class TransportValidationFixture
{
   protected:
    CoolPropDbl actual, x1, x2;
    shared_ptr<CoolProp::AbstractState> pState;
    CoolProp::input_pairs pair;

   public:
    TransportValidationFixture() {}
    ~TransportValidationFixture() {}
    void set_backend(std::string backend, std::string fluid_name) {
        pState.reset(CoolProp::AbstractState::factory(backend, fluid_name));
    }
    void set_pair(std::string& in1, double v1, std::string& in2, double v2) {
        double o1, o2;
        parameters iin1 = CoolProp::get_parameter_index(in1);
        parameters iin2 = CoolProp::get_parameter_index(in2);
        CoolProp::input_pairs pair = CoolProp::generate_update_pair(iin1, v1, iin2, v2, o1, o2);
        pState->update(pair, o1, o2);
    }
    void get_value(parameters key) {
        actual = pState->keyed_output(key);
    }
};

TEST_CASE_METHOD(TransportValidationFixture, "Compare viscosities against published data", "[viscosity],[transport]") {
    int inputsN = sizeof(viscosity_validation_data) / sizeof(viscosity_validation_data[0]);
    for (int i = 0; i < inputsN; ++i) {
        vel el = viscosity_validation_data[i];
        CHECK_NOTHROW(set_backend("HEOS", el.fluid));

        CAPTURE(el.fluid);
        CAPTURE(el.in1);
        CAPTURE(el.v1);
        CAPTURE(el.in2);
        CAPTURE(el.v2);
        CHECK_NOTHROW(set_pair(el.in1, el.v1, el.in2, el.v2));
        CHECK_NOTHROW(get_value(CoolProp::iviscosity));
        CAPTURE(el.expected);
        CAPTURE(actual);
        CHECK(std::abs(actual / el.expected - 1) < el.tol);
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
  vel("Heptane", "T", 250, "Dmass", 720, "L", 137.09e-3, 1e-4),
  vel("Heptane", "T", 400, "Dmass", 2, "L", 21.794e-3, 1e-4),
  vel("Heptane", "T", 400, "Dmass", 650, "L", 120.75e-3, 1e-4),
  vel("Heptane", "T", 535, "Dmass", 100, "L", 51.655e-3, 3e-3),  // Relaxed tolerance because conductivity was fit using older viscosity correlation

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
  vel("SF6", "T", 298.15, "Dmass", 100, "L", 14.126e-3, 1e-4),
  vel("SF6", "T", 298.15, "Dmass", 1600, "L", 69.729e-3, 1e-4),
  vel("SF6", "T", 310, "Dmass", 1e-13, "L", 13.834e-3, 1e-4),
  vel("SF6", "T", 310, "Dmass", 1200, "L", 48.705e-3, 1e-4),
  vel("SF6", "T", 480, "Dmass", 100, "L", 28.847e-3, 1e-4),

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
  vel("Hydrogen", "T", 35, "Dmass", 30, "L", 0.0770177, 1e-4),  // Updated since Assael uses a different viscosity correlation
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
  vel("R125", "T", 341, "Dmass", 600, "L", 0.0565642978494, 2e-4),
  vel("R125", "T", 200, "Dmass", 1e-13, "L", 0.007036843623086, 2e-4),
  vel("IsoButane", "T", 390, "Dmass", 387.09520158645068, "L", 0.063039, 2e-4),
  vel("IsoButane", "T", 390, "Dmass", 85.76703973869482, "L", 0.036603, 2e-4),
  vel("n-Butane", "T", 415, "Dmass", 360.01895129934866, "L", 0.067045, 2e-4),
  vel("n-Butane", "T", 415, "Dmass", 110.3113177144, "L", 0.044449, 1e-4),

  // From Huber, FPE, 2005
  vel("n-Octane", "T", 300, "Dmolar", 6177.2, "L", 0.12836, 1e-4),
  vel("n-Nonane", "T", 300, "Dmolar", 5619.4, "L", 0.13031, 1e-4),
  //vel("n-Decane", "T", 300, "Dmass", 5150.4, "L", 0.13280, 1e-4), // no viscosity

  // From Huber, EF, 2004
  vel("n-Dodecane", "T", 300, "Dmolar", 4411.5, "L", 0.13829, 1e-4),
  vel("n-Dodecane", "T", 500, "Dmolar", 3444.7, "L", 0.09384, 1e-4),
  vel("n-Dodecane", "T", 660, "Dmolar", 1500.98, "L", 0.090346, 1e-4),

  // From REFPROP 9.1 since no data provided in Marsh, 2002
  vel("n-Propane", "T", 368, "Q", 0, "L", 0.07282154952457, 1e-3),
  vel("n-Propane", "T", 368, "Dmolar", 1e-10, "L", 0.0266135388745317, 1e-4),

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
  vel("R123", "T", 180, "Dmass", 1739, "L", 110.9e-3, 2e-4),
  vel("R123", "T", 180, "Dmass", 0.2873e-2, "L", 2.473e-3, 1e-3),
  vel("R123", "T", 430, "Dmass", 996.35, "L", 45.62e-3, 1e-3),
  vel("R123", "T", 430, "Dmass", 166.9, "L", 21.03e-3, 1e-3),

  // From Scalabrin, JPCRD, 2006
  vel("CO2", "T", 218, "Q", 0, "L", 181.09e-3, 1e-4),
  vel("CO2", "T", 218, "Q", 1, "L", 10.837e-3, 1e-4),
  vel("CO2", "T", 304, "Q", 0, "L", 140.3e-3, 1e-4),
  vel("CO2", "T", 304, "Q", 1, "L", 217.95e-3, 1e-4),
  vel("CO2", "T", 225, "Dmass", 0.23555, "L", 11.037e-3, 1e-4),
  vel("CO2", "T", 275, "Dmass", 1281.64, "L", 238.44e-3, 1e-4),

  // From Friend, JPCRD, 1991
  vel("Ethane", "T", 100, "Dmass", 1e-13, "L", 3.46e-3, 1e-2),
  vel("Ethane", "T", 230, "Dmolar", 16020, "L", 126.2e-3, 1e-2),
  vel("Ethane", "T", 440, "Dmolar", 1520, "L", 45.9e-3, 1e-2),
  vel("Ethane", "T", 310, "Dmolar", 4130, "L", 45.4e-3, 1e-2),

  // From Lemmon and Jacobsen, JPCRD, 2004
  vel("Nitrogen", "T", 100, "Dmolar", 1e-14, "L", 9.27749e-3, 1e-4),
  vel("Nitrogen", "T", 300, "Dmolar", 1e-14, "L", 25.9361e-3, 1e-4),
  vel("Nitrogen", "T", 100, "Dmolar", 25000, "L", 103.834e-3, 1e-4),
  vel("Nitrogen", "T", 200, "Dmolar", 10000, "L", 36.0099e-3, 1e-4),
  vel("Nitrogen", "T", 300, "Dmolar", 5000, "L", 32.7694e-3, 1e-4),
  vel("Nitrogen", "T", 126.195, "Dmolar", 11180, "L", 675.800e-3, 1e-4),
  vel("Argon", "T", 100, "Dmolar", 1e-14, "L", 6.36587e-3, 1e-4),
  vel("Argon", "T", 300, "Dmolar", 1e-14, "L", 17.8042e-3, 1e-4),
  vel("Argon", "T", 100, "Dmolar", 33000, "L", 111.266e-3, 1e-4),
  vel("Argon", "T", 200, "Dmolar", 10000, "L", 26.1377e-3, 1e-4),
  vel("Argon", "T", 300, "Dmolar", 5000, "L", 23.2302e-3, 1e-4),
  vel("Argon", "T", 150.69, "Dmolar", 13400, "L", 856.793e-3, 1e-4),
  vel("Oxygen", "T", 100, "Dmolar", 1e-14, "L", 8.94334e-3, 1e-4),
  vel("Oxygen", "T", 300, "Dmolar", 1e-14, "L", 26.4403e-3, 1e-4),
  vel("Oxygen", "T", 100, "Dmolar", 35000, "L", 146.044e-3, 1e-4),
  vel("Oxygen", "T", 200, "Dmolar", 10000, "L", 34.6124e-3, 1e-4),
  vel("Oxygen", "T", 300, "Dmolar", 5000, "L", 32.5491e-3, 1e-4),
  vel("Oxygen", "T", 154.6, "Dmolar", 13600, "L", 377.476e-3, 1e-4),
  vel("Air", "T", 100, "Dmolar", 1e-14, "L", 9.35902e-3, 1e-4),
  vel("Air", "T", 300, "Dmolar", 1e-14, "L", 26.3529e-3, 1e-4),
  vel("Air", "T", 100, "Dmolar", 28000, "L", 119.221e-3, 1e-4),
  vel("Air", "T", 200, "Dmolar", 10000, "L", 35.3185e-3, 1e-4),
  vel("Air", "T", 300, "Dmolar", 5000, "L", 32.6062e-3, 1e-4),
  vel("Air", "T", 132.64, "Dmolar", 10400, "L", 75.6231e-3, 1e-4),

  // Huber, JPCRD, 2012
  vel("Water", "T", 298.15, "Dmass", 1e-14, "L", 18.4341883e-3, 1e-6),
  vel("Water", "T", 298.15, "Dmass", 998, "L", 607.712868e-3, 1e-6),
  vel("Water", "T", 298.15, "Dmass", 1200, "L", 799.038144e-3, 1e-6),
  vel("Water", "T", 873.15, "Dmass", 1e-14, "L", 79.1034659e-3, 1e-6),
  vel("Water", "T", 647.35, "Dmass", 1, "L", 51.9298924e-3, 1e-6),
  vel("Water", "T", 647.35, "Dmass", 122, "L", 130.922885e-3, 2e-4),
  vel("Water", "T", 647.35, "Dmass", 222, "L", 367.787459e-3, 2e-4),
  vel("Water", "T", 647.35, "Dmass", 272, "L", 757.959776e-3, 2e-4),
  vel("Water", "T", 647.35, "Dmass", 322, "L", 1443.75556e-3, 2e-4),
  vel("Water", "T", 647.35, "Dmass", 372, "L", 650.319402e-3, 2e-4),
  vel("Water", "T", 647.35, "Dmass", 422, "L", 448.883487e-3, 2e-4),
  vel("Water", "T", 647.35, "Dmass", 750, "L", 600.961346e-3, 2e-4),

  // From Shan, ASHRAE, 2000
  vel("R23", "T", 180, "Dmolar", 21097, "L", 143.19e-3, 1e-4),
  vel("R23", "T", 420, "Dmolar", 7564, "L", 50.19e-3, 2e-4),
  vel("R23", "T", 370, "Dmolar", 32.62, "L", 17.455e-3, 1e-4),

  // From REFPROP 9.1 since no sample data provided in Tufeu
  vel("Ammonia", "T", 310, "Dmolar", 34320, "L", 0.45223303481784971, 1e-4),
  vel("Ammonia", "T", 395, "Q", 0, "L", 0.2264480769301, 1e-4),

  // From Hands, Cryogenics, 1981
  vel("Helium", "T", 800, "P", 1e5, "L", 0.3085, 1e-2),
  vel("Helium", "T", 300, "P", 1e5, "L", 0.1560, 1e-2),
  vel("Helium", "T", 20, "P", 1e5, "L", 0.0262, 1e-2),
  vel("Helium", "T", 8, "P", 1e5, "L", 0.0145, 1e-2),
  vel("Helium", "T", 4, "P", 20e5, "L", 0.0255, 1e-2),
  vel("Helium", "T", 8, "P", 20e5, "L", 0.0308, 1e-2),
  vel("Helium", "T", 20, "P", 20e5, "L", 0.0328, 1e-2),
  vel("Helium", "T", 4, "P", 100e5, "L", 0.0385, 3e-2),
  vel("Helium", "T", 8, "P", 100e5, "L", 0.0566, 3e-2),
  vel("Helium", "T", 20, "P", 100e5, "L", 0.0594, 1e-2),
  vel("Helium", "T", 4, "P", 1e5, "L", 0.0186, 1e-2),
  vel("Helium", "T", 4, "P", 2e5, "L", 0.0194, 1e-2),
  vel("Helium", "T", 5.180, "P", 2.3e5, "L", 0.0195, 1e-1),
  vel("Helium", "T", 5.2, "P", 2.3e5, "L", 0.0202, 1e-1),
  vel("Helium", "T", 5.230, "P", 2.3e5, "L", 0.0181, 1e-1),
  vel("Helium", "T", 5.260, "P", 2.3e5, "L", 0.0159, 1e-1),
  vel("Helium", "T", 5.3, "P", 2.3e5, "L", 0.0149, 1e-1),

  // Geller, IJT, 2001 - based on experimental data, no validation data provided
  //vel("R404A", "T", 253.03, "P", 0.101e6, "L", 0.00991, 0.03),
  //vel("R404A", "T", 334.38, "P", 2.176e6, "L", 19.93e-3, 0.03),
  //vel("R407C", "T", 253.45, "P", 0.101e6, "L", 0.00970, 0.03),
  //vel("R407C", "T", 314.39, "P", 0.458e6, "L", 14.87e-3, 0.03),
  //vel("R410A", "T", 260.32, "P", 0.101e6, "L", 0.01043, 0.03),
  //vel("R410A", "T", 332.09, "P", 3.690e6, "L", 22.76e-3, 0.03),
  //vel("R507A", "T", 254.85, "P", 0.101e6, "L", 0.01007, 0.03),
  //vel("R507A", "T", 333.18, "P", 2.644e6, "L", 21.31e-3, 0.03),

  // From REFPROP 9.1 since no data provided
  vel("R134a", "T", 240, "D", 1e-10, "L", 0.008698768, 1e-4),
  vel("R134a", "T", 330, "D", 1e-10, "L", 0.015907606, 1e-4),
  vel("R134a", "T", 330, "Q", 0, "L", 0.06746432253, 1e-4),
  vel("R134a", "T", 240, "Q", 1, "L", 0.00873242359, 1e-4),

  // Mylona, JPCRD, 2014
  vel("o-Xylene", "T", 635, "D", 270, "L", 96.4e-3, 1e-2),
  vel("m-Xylene", "T", 616, "D", 220, "L", 79.5232e-3, 1e-2),  // CoolProp is correct, paper is incorrect (it seems)
  vel("p-Xylene", "T", 620, "D", 287, "L", 107.7e-3, 1e-2),
  vel("EthylBenzene", "T", 617, "D", 316, "L", 140.2e-3, 1e-2),
  // dilute values
  vel("o-Xylene", "T", 300, "D", 1e-12, "L", 13.68e-3, 1e-3),
  vel("o-Xylene", "T", 600, "D", 1e-12, "L", 41.6e-3, 1e-3),
  vel("m-Xylene", "T", 300, "D", 1e-12, "L", 9.45e-3, 1e-3),
  vel("m-Xylene", "T", 600, "D", 1e-12, "L", 40.6e-3, 1e-3),
  vel("p-Xylene", "T", 300, "D", 1e-12, "L", 10.57e-3, 1e-3),
  vel("p-Xylene", "T", 600, "D", 1e-12, "L", 41.73e-3, 1e-3),
  vel("EthylBenzene", "T", 300, "D", 1e-12, "L", 9.71e-3, 1e-3),
  vel("EthylBenzene", "T", 600, "D", 1e-12, "L", 41.14e-3, 1e-3),

  // Friend, JPCRD, 1989
  vel("Methane", "T", 100, "D", 1e-12, "L", 9.83e-3, 1e-3),
  vel("Methane", "T", 400, "D", 1e-12, "L", 49.96e-3, 1e-3),
  vel("Methane", "T", 182, "Q", 0, "L", 82.5e-3, 5e-3),
  vel("Methane", "T", 100, "Dmolar", 28.8e3, "L", 234e-3, 1e-2),

  // Sykioti, JPCRD, 2013
  vel("Methanol", "T", 300, "Dmass", 850, "L", 241.48e-3, 1e-2),
  vel("Methanol", "T", 400, "Dmass", 2, "L", 25.803e-3, 1e-2),
  vel("Methanol", "T", 400, "Dmass", 690, "L", 183.59e-3, 1e-2),
  vel("Methanol", "T", 500, "Dmass", 10, "L", 40.495e-3, 1e-2),

  // Heavy Water, IAPWS formulation
  vel("HeavyWater", "T", 0.5000 * 643.847, "Dmass", 3.07 * 358, "V", 835.786416818 * 0.742128e-3, 1e-5),
  vel("HeavyWater", "T", 0.9000 * 643.847, "Dmass", 2.16 * 358, "V", 627.777590127 * 0.742128e-3, 1e-5),
  vel("HeavyWater", "T", 1.2000 * 643.847, "Dmass", 0.8 * 358, "V", 259.605241187 * 0.742128e-3, 1e-5),

  // Vassiliou, JPCRD, 2015
  vel("Cyclopentane", "T", 512, "Dmass", 1e-12, "L", 37.042e-3, 1e-5),
  vel("Cyclopentane", "T", 512, "Dmass", 400, "L", 69.698e-3, 1e-1),
  vel("Isopentane", "T", 460, "Dmass", 1e-12, "L", 35.883e-3, 1e-4),
  vel("Isopentane", "T", 460, "Dmass", 329.914, "L", 59.649e-3, 1e-1),
  vel("n-Pentane", "T", 460, "Dmass", 1e-12, "L", 34.048e-3, 1e-5),
  vel("n-Pentane", "T", 460, "Dmass", 377.687, "L", 71.300e-3, 1e-1),
};

TEST_CASE_METHOD(TransportValidationFixture, "Compare thermal conductivities against published data", "[conductivity],[transport]") {
    int inputsN = sizeof(conductivity_validation_data) / sizeof(conductivity_validation_data[0]);
    for (int i = 0; i < inputsN; ++i) {
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
        CHECK(std::abs(actual / el.expected - 1) < el.tol);
    }
}

}; /* namespace TransportValidation */

static CoolProp::input_pairs inputs[] = {
  CoolProp::DmolarT_INPUTS,
  //CoolProp::SmolarT_INPUTS,
  //CoolProp::HmolarT_INPUTS,
  //CoolProp::TUmolar_INPUTS,

  //    CoolProp::DmolarP_INPUTS,
  //    CoolProp::DmolarHmolar_INPUTS,
  //    CoolProp::DmolarSmolar_INPUTS,
  //    CoolProp::DmolarUmolar_INPUTS,
  //
  //    CoolProp::HmolarP_INPUTS,
  //    CoolProp::PSmolar_INPUTS,
  //    CoolProp::PUmolar_INPUTS,
  //
  /*
    CoolProp::HmolarSmolar_INPUTS,
    CoolProp::HmolarUmolar_INPUTS,
    CoolProp::SmolarUmolar_INPUTS
    */
};

class ConsistencyFixture
{
   protected:
    CoolPropDbl hmolar, pmolar, smolar, umolar, rhomolar, T, p, x1, x2;
    shared_ptr<CoolProp::AbstractState> pState;
    CoolProp::input_pairs pair;

   public:
    ConsistencyFixture() {}
    ~ConsistencyFixture() {}
    void set_backend(std::string backend, std::string fluid_name) {
        pState.reset(CoolProp::AbstractState::factory(backend, fluid_name));
    }
    void set_pair(CoolProp::input_pairs pair) {
        this->pair = pair;
    }
    void set_TP(CoolPropDbl T, CoolPropDbl p) {
        this->T = T;
        this->p = p;
        CoolProp::AbstractState& State = *pState;

        // Start with T,P as inputs, cycle through all the other pairs that are supported
        State.update(CoolProp::PT_INPUTS, p, T);

        // Set the other state variables
        rhomolar = State.rhomolar();
        hmolar = State.hmolar();
        smolar = State.smolar();
        umolar = State.umolar();
    }
    void get_variables() {

        switch (pair) {
            /// In this group, T is one of the known inputs, iterate for the other one (easy)
            case CoolProp::HmolarT_INPUTS:
                x1 = hmolar;
                x2 = T;
                break;
            case CoolProp::SmolarT_INPUTS:
                x1 = smolar;
                x2 = T;
                break;
            case CoolProp::TUmolar_INPUTS:
                x1 = T;
                x2 = umolar;
                break;
            case CoolProp::DmolarT_INPUTS:
                x1 = rhomolar;
                x2 = T;
                break;

            /// In this group, D is one of the known inputs, iterate for the other one (a little bit harder)
            case CoolProp::DmolarHmolar_INPUTS:
                x1 = rhomolar;
                x2 = hmolar;
                break;
            case CoolProp::DmolarSmolar_INPUTS:
                x1 = rhomolar;
                x2 = smolar;
                break;
            case CoolProp::DmolarUmolar_INPUTS:
                x1 = rhomolar;
                x2 = umolar;
                break;
            case CoolProp::DmolarP_INPUTS:
                x1 = rhomolar;
                x2 = p;
                break;

            /// In this group, p is one of the known inputs (a little less easy)
            case CoolProp::HmolarP_INPUTS:
                x1 = hmolar;
                x2 = p;
                break;
            case CoolProp::PSmolar_INPUTS:
                x1 = p;
                x2 = smolar;
                break;
            case CoolProp::PUmolar_INPUTS:
                x1 = p;
                x2 = umolar;
                break;

            case CoolProp::HmolarSmolar_INPUTS:
                x1 = hmolar;
                x2 = smolar;
                break;
            case CoolProp::SmolarUmolar_INPUTS:
                x1 = smolar;
                x2 = umolar;
                break;

            default:
                throw CoolProp::ValueError();
        }
    }
    void single_phase_consistency_check() {
        CoolProp::AbstractState& State = *pState;
        State.update(pair, x1, x2);

        // Make sure we end up back at the same temperature and pressure we started out with
        if (State.Q() < 1 && State.Q() > 0) throw CoolProp::ValueError(format("Q [%g] is between 0 and 1; two-phase solution", State.Q()));
        if (std::abs(T - State.T()) > 1e-2) throw CoolProp::ValueError(format("Error on T [%Lg K] is greater than 1e-2", std::abs(State.T() - T)));
        if (std::abs(p - State.p()) / p * 100 > 1e-2)
            throw CoolProp::ValueError(format("Error on p [%Lg %%] is greater than 1e-2 %%", std::abs(p - State.p()) / p * 100));
    }
    void subcritical_pressure_liquid() {
        // Subcritical pressure liquid
        int inputsN = sizeof(inputs) / sizeof(inputs[0]);
        for (double p = pState->p_triple() * 1.1; p < pState->p_critical(); p *= 3) {
            double Ts = PropsSI("T", "P", p, "Q", 0, "Water");
            double Tmelt = pState->melting_line(CoolProp::iT, CoolProp::iP, p);
            for (double T = Tmelt; T < Ts - 0.1; T += 0.1) {
                CHECK_NOTHROW(set_TP(T, p));

                for (int i = 0; i < inputsN; ++i) {
                    CoolProp::input_pairs pair = inputs[i];
                    std::string pair_desc = CoolProp::get_input_pair_short_desc(pair);
                    set_pair(pair);
                    CAPTURE(pair_desc);
                    CAPTURE(T);
                    CAPTURE(p);
                    get_variables();
                    CAPTURE(x1);
                    CAPTURE(x2);
                    CAPTURE(Ts);
                    CHECK_NOTHROW(single_phase_consistency_check());
                    double rhomolar_RP = PropsSI("Dmolar", "P", p, "T", T, "REFPROP::Water");
                    if (ValidNumber(rhomolar_RP)) {
                        CAPTURE(rhomolar_RP);
                        CAPTURE(rhomolar);
                        CHECK(std::abs((rhomolar_RP - rhomolar) / rhomolar) < 1e-3);
                    }
                }
            }
        }
    }
};

TEST_CASE_METHOD(ConsistencyFixture, "Test all input pairs for Water using all valid backends", "[consistency]") {
    CHECK_NOTHROW(set_backend("HEOS", "Water"));
    subcritical_pressure_liquid();

    //    int inputsN = sizeof(inputs)/sizeof(inputs[0]);
    //    for (double p = 600000; p < pState->pmax(); p *= 3)
    //    {
    //        for (double T = 220; T < pState->Tmax(); T += 1)
    //        {
    //            CHECK_NOTHROW(set_TP(T, p));
    //
    //            for (int i = 0; i < inputsN; ++i)
    //            {
    //                CoolProp::input_pairs pair = inputs[i];
    //                std::string pair_desc = CoolProp::get_input_pair_short_desc(pair);
    //                set_pair(pair);
    //                CAPTURE(pair_desc);
    //                CAPTURE(T);
    //                CAPTURE(p);
    //                get_variables();
    //                CAPTURE(x1);
    //                CAPTURE(x2);
    //                CHECK_NOTHROW(single_phase_consistency_check());
    //            }
    //        }
    //    }
}

TEST_CASE("Test saturation properties for a few fluids", "[saturation],[slow]") {
    SECTION("sat_p") {
        std::vector<double> pv = linspace(Props1SI("CO2", "ptriple"), Props1SI("CO2", "pcrit") - 1e-6, 5);

        SECTION("All pressures are ok")
        for (std::size_t i = 0; i < pv.size(); ++i) {
            CAPTURE(pv[i]);
            double T = CoolProp::PropsSI("T", "P", pv[i], "Q", 0, "CO2");
        }
    }
}

class HumidAirDewpointFixture
{
   public:
    shared_ptr<CoolProp::AbstractState> AS;
    std::vector<std::string> fluids;
    std::vector<double> z;
    void setup(double zH2O) {
        double z_Air[4] = {0.7810, 0.2095, 0.0092, 0.0003};  // N2, O2, Ar, CO2
        z.resize(5);
        z[0] = zH2O;
        for (int i = 0; i < 4; ++i) {
            z[i + 1] = (1 - zH2O) * z_Air[i];
        }
    }
    void run_p(double p) {
        CAPTURE(p);
        for (double zH2O = 0.999; zH2O > 0; zH2O -= 0.001) {
            setup(zH2O);
            AS->set_mole_fractions(z);
            CAPTURE(zH2O);
            CHECK_NOTHROW(AS->update(PQ_INPUTS, p, 1));
            if (AS->T() < 273.15) {
                break;
            }
        }
    }
    void run_checks() {
        fluids = strsplit("Water&Nitrogen&Oxygen&Argon&CO2", '&');
        AS.reset(AbstractState::factory("HEOS", fluids));
        run_p(1e5);
        run_p(1e6);
        run_p(1e7);
    }
};
TEST_CASE_METHOD(HumidAirDewpointFixture, "Humid air dewpoint calculations", "[humid_air_dewpoint]") {
    run_checks();
}

TEST_CASE("Test consistency between Gernert models in CoolProp and Gernert models in REFPROP", "[Gernert]") {
    // See https://groups.google.com/forum/?fromgroups#!topic/catch-forum/mRBKqtTrITU
    std::string mixes[] = {"CO2[0.7]&Argon[0.3]", "CO2[0.7]&Water[0.3]", "CO2[0.7]&Nitrogen[0.3]"};
    for (int i = 0; i < 3; ++i) {
        const char* ykey = mixes[i].c_str();
        std::ostringstream ss1;
        ss1 << mixes[i];
        SECTION(ss1.str(), "") {
            double Tnbp_CP, Tnbp_RP;
            CHECK_NOTHROW(Tnbp_CP = PropsSI("T", "P", 101325, "Q", 1, "HEOS::" + mixes[i]));
            CAPTURE(Tnbp_CP);
            CHECK_NOTHROW(Tnbp_RP = PropsSI("T", "P", 101325, "Q", 1, "REFPROP::" + mixes[i]));
            CAPTURE(Tnbp_RP);
            double diff = std::abs(Tnbp_CP / Tnbp_RP - 1);
            CHECK(diff < 1e-6);
        }
    }
}

TEST_CASE("Tests for solvers in P,T flash using Water", "[flash],[PT]") {
    SECTION("Check that T,P for saturated state yields error") {
        double Ts, ps, rho;
        CHECK_NOTHROW(Ts = PropsSI("T", "P", 101325, "Q", 0, "Water"));
        CHECK(ValidNumber(Ts));
        CHECK_NOTHROW(ps = PropsSI("P", "T", Ts, "Q", 0, "Water"));
        CHECK(ValidNumber(ps));
        CAPTURE(Ts);
        CAPTURE(ps);
        CHECK_NOTHROW(rho = PropsSI("D", "T", Ts, "P", ps, "Water"));
        CAPTURE(rho);
        CHECK(!ValidNumber(rho));
    }
    SECTION("Subcritical p slightly subcooled should be ok") {
        double Ts, rho, dT = 1e-4;
        CHECK_NOTHROW(Ts = PropsSI("T", "P", 101325, "Q", 0, "Water"));
        CAPTURE(Ts);
        CHECK(ValidNumber(Ts));
        CAPTURE(dT);
        CHECK_NOTHROW(rho = PropsSI("D", "T", Ts - dT, "P", 101325, "Water"));
        CAPTURE(rho);
        CHECK(ValidNumber(rho));
    }
    SECTION("Subcritical p slightly superheated should be ok") {
        double Ts, rho, dT = 1e-4;
        CHECK_NOTHROW(Ts = PropsSI("T", "P", 101325, "Q", 0, "Water"));
        CAPTURE(Ts);
        CHECK(ValidNumber(Ts));
        CAPTURE(dT);
        CHECK_NOTHROW(rho = PropsSI("D", "T", Ts + dT, "P", 101325, "Water"));
        CAPTURE(rho);
        CHECK(ValidNumber(rho));
    }
}

TEST_CASE("Tests for solvers in P,Y flash using Water", "[flash],[PH],[PS],[PU]") {
    double Ts, y, T2;
    // See https://groups.google.com/forum/?fromgroups#!topic/catch-forum/mRBKqtTrITU
    std::string Ykeys[] = {"H", "S", "U", "Hmass", "Smass", "Umass", "Hmolar", "Smolar", "Umolar"};
    for (int i = 0; i < 9; ++i) {
        const char* ykey = Ykeys[i].c_str();
        std::ostringstream ss1;
        ss1 << "Subcritical superheated P," << ykey;
        SECTION(ss1.str(), "") {
            double dT = 10;
            CHECK_NOTHROW(Ts = PropsSI("T", "P", 101325, "Q", 0, "Water"));
            CHECK(ValidNumber(Ts));
            CAPTURE(Ts);
            CHECK_NOTHROW(y = PropsSI(ykey, "T", Ts + dT, "P", 101325, "Water"));
            CAPTURE(dT);
            CAPTURE(y);
            CHECK(ValidNumber(y));
            CHECK_NOTHROW(T2 = PropsSI("T", ykey, y, "P", 101325, "Water"));
            CAPTURE(CoolProp::get_global_param_string("errstring"));
            CAPTURE(T2);
            CHECK(ValidNumber(T2));
        }
        std::ostringstream ss2;
        ss2 << "Subcritical barely superheated P," << ykey;
        SECTION(ss2.str(), "") {
            double dT = 1e-3;
            CHECK_NOTHROW(Ts = PropsSI("T", "P", 101325, "Q", 0, "Water"));
            CHECK(ValidNumber(Ts));
            CAPTURE(Ts);
            CHECK_NOTHROW(y = PropsSI(ykey, "T", Ts + dT, "P", 101325, "Water"));
            CAPTURE(dT);
            CAPTURE(y);
            CHECK(ValidNumber(y));
            CHECK_NOTHROW(T2 = PropsSI("T", ykey, y, "P", 101325, "Water"));
            CAPTURE(CoolProp::get_global_param_string("errstring"));
            CAPTURE(T2);
            CHECK(ValidNumber(T2));
        }
        std::ostringstream ss3;
        ss3 << "Subcritical subcooled P," << ykey;
        SECTION(ss3.str(), "") {
            double dT = -10;
            CHECK_NOTHROW(Ts = PropsSI("T", "P", 101325, "Q", 0, "Water"));
            CHECK(ValidNumber(Ts));
            CAPTURE(Ts);
            CHECK_NOTHROW(y = PropsSI(ykey, "T", Ts + dT, "P", 101325, "Water"));
            CAPTURE(dT);
            CAPTURE(y);
            CHECK(ValidNumber(y));
            CHECK_NOTHROW(T2 = PropsSI("T", ykey, y, "P", 101325, "Water"));
            CAPTURE(CoolProp::get_global_param_string("errstring"));
            CAPTURE(T2);
            CHECK(ValidNumber(T2));
        }
        std::ostringstream ss4;
        ss4 << "Subcritical barely subcooled P," << ykey;
        SECTION(ss4.str(), "") {
            double dT = -1e-3;
            CHECK_NOTHROW(Ts = PropsSI("T", "P", 101325, "Q", 0, "Water"));
            CHECK(ValidNumber(Ts));
            CAPTURE(Ts);
            CHECK_NOTHROW(y = PropsSI(ykey, "T", Ts + dT, "P", 101325, "Water"));
            CAPTURE(dT);
            CAPTURE(y);
            CHECK(ValidNumber(y));
            CHECK_NOTHROW(T2 = PropsSI("T", ykey, y, "P", 101325, "Water"));
            CAPTURE(CoolProp::get_global_param_string("errstring"));
            CAPTURE(T2);
            CHECK(ValidNumber(T2));
        }
        std::ostringstream ss5;
        ss5 << "Supercritical P," << ykey;
        SECTION(ss5.str(), "") {
            double Tc = Props1SI("Water", "Tcrit");
            double pc = Props1SI("Water", "pcrit");
            double p = pc * 1.3;
            double T = Tc * 1.3;
            CAPTURE(T);
            CAPTURE(p);
            CHECK(ValidNumber(T));
            CHECK(ValidNumber(p));
            CHECK_NOTHROW(y = PropsSI(ykey, "P", p, "T", T, "Water"));
            CAPTURE(y);
            CHECK(ValidNumber(y));
            CHECK_NOTHROW(T2 = PropsSI("T", ykey, y, "P", p, "Water"));
            CAPTURE(CoolProp::get_global_param_string("errstring"));
            CAPTURE(T2);
            CHECK(ValidNumber(T2));
        }
        std::ostringstream ss6;
        ss6 << "Supercritical \"gas\" P," << ykey;
        SECTION(ss6.str(), "") {
            double Tc = Props1SI("Water", "Tcrit");
            double pc = Props1SI("Water", "pcrit");
            double p = pc * 0.7;
            double T = Tc * 1.3;
            CAPTURE(T);
            CAPTURE(p);
            CHECK(ValidNumber(T));
            CHECK(ValidNumber(p));
            CHECK_NOTHROW(y = PropsSI(ykey, "P", p, "T", T, "Water"));
            CAPTURE(y);
            CHECK(ValidNumber(y));
            CHECK_NOTHROW(T2 = PropsSI("T", ykey, y, "P", p, "Water"));
            CAPTURE(CoolProp::get_global_param_string("errstring"));
            CAPTURE(T2);
            CHECK(ValidNumber(T2));
        }
        std::ostringstream ss7;
        ss7 << "Supercritical \"liquid\" P," << ykey;
        SECTION(ss7.str(), "") {
            double Tc = Props1SI("Water", "Tcrit");
            double pc = Props1SI("Water", "pcrit");
            double p = pc * 2;
            double T = Tc * 0.5;
            CAPTURE(T);
            CAPTURE(p);
            CHECK(ValidNumber(T));
            CHECK(ValidNumber(p));
            CHECK_NOTHROW(y = PropsSI(ykey, "P", p, "T", T, "Water"));
            CAPTURE(y);
            CHECK(ValidNumber(y));
            CHECK_NOTHROW(T2 = PropsSI("T", ykey, y, "P", p, "Water"));
            CAPTURE(CoolProp::get_global_param_string("errstring"));
            CAPTURE(T2);
            CHECK(ValidNumber(T2));
        }
    }
}

TEST_CASE("Tests for solvers in P,H flash using Propane", "[flashdups],[flash],[PH],[consistency]") {
    double hmolar, hmass;
    SECTION("5 times PH with HEOS AbstractState yields same results every time", "") {
        shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "n-Propane"));

        CHECK_NOTHROW(AS->update(CoolProp::PT_INPUTS, 101325, 300));
        hmolar = AS->hmolar();
        hmass = AS->hmass();
        CHECK_NOTHROW(AS->update(CoolProp::HmassP_INPUTS, hmass, 101325));
        CHECK_NOTHROW(AS->update(CoolProp::HmolarP_INPUTS, hmolar, 101325));
        hmolar = AS->hmolar();
        hmass = AS->hmass();
        CHECK_NOTHROW(AS->update(CoolProp::HmassP_INPUTS, hmass, 101325));
        CHECK_NOTHROW(AS->update(CoolProp::HmolarP_INPUTS, hmolar, 101325));
        hmolar = AS->hmolar();
        hmass = AS->hmass();
        CHECK_NOTHROW(AS->update(CoolProp::HmassP_INPUTS, hmass, 101325));
        CHECK_NOTHROW(AS->update(CoolProp::HmolarP_INPUTS, hmolar, 101325));
        hmolar = AS->hmolar();
        hmass = AS->hmass();
        CHECK_NOTHROW(AS->update(CoolProp::HmassP_INPUTS, hmass, 101325));
        CHECK_NOTHROW(AS->update(CoolProp::HmolarP_INPUTS, hmolar, 101325));
        hmolar = AS->hmolar();
        hmass = AS->hmass();
        CHECK_NOTHROW(AS->update(CoolProp::HmassP_INPUTS, hmass, 101325));
        CHECK_NOTHROW(AS->update(CoolProp::HmolarP_INPUTS, hmolar, 101325));
    }
}

TEST_CASE("Multiple calls to state class are consistent", "[flashdups],[flash],[PH],[consistency]") {
    double hmolar, hmass;
    SECTION("3 times PH with HEOS AbstractState yields same results every time", "") {
        shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "n-Propane"));

        CHECK_NOTHROW(AS->update(CoolProp::PT_INPUTS, 101325, 300));
        hmolar = AS->hmolar();
        hmass = AS->hmass();
        CHECK_NOTHROW(AS->update(CoolProp::HmassP_INPUTS, hmass, 101325));
        CHECK_NOTHROW(AS->update(CoolProp::HmolarP_INPUTS, hmolar, 101325));
        hmolar = AS->hmolar();
        hmass = AS->hmass();
        CHECK_NOTHROW(AS->update(CoolProp::HmassP_INPUTS, hmass, 101325));
        CHECK_NOTHROW(AS->update(CoolProp::HmolarP_INPUTS, hmolar, 101325));
        hmolar = AS->hmolar();
        hmass = AS->hmass();
        CHECK_NOTHROW(AS->update(CoolProp::HmassP_INPUTS, hmass, 101325));
        CHECK_NOTHROW(AS->update(CoolProp::HmolarP_INPUTS, hmolar, 101325));
    }
}

TEST_CASE("Test first partial derivatives using PropsSI", "[derivatives]") {
    double T = 300;
    SECTION("Check drhodp|T 3 ways", "") {
        shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "n-Propane"));
        AS->update(CoolProp::PT_INPUTS, 101325, T);

        double drhomolardp__T_AbstractState = AS->first_partial_deriv(CoolProp::iDmolar, CoolProp::iP, CoolProp::iT);
        double drhomolardp__T_PropsSI_num =
          (PropsSI("Dmolar", "T", T, "P", 101325 + 1e-3, "n-Propane") - PropsSI("Dmolar", "T", T, "P", 101325 - 1e-3, "n-Propane")) / (2 * 1e-3);
        double drhomolardp__T_PropsSI = PropsSI("d(Dmolar)/d(P)|T", "T", T, "P", 101325, "n-Propane");

        CAPTURE(drhomolardp__T_AbstractState);
        CAPTURE(drhomolardp__T_PropsSI_num);
        CAPTURE(drhomolardp__T_PropsSI);
        double rel_err_exact = std::abs((drhomolardp__T_AbstractState - drhomolardp__T_PropsSI) / drhomolardp__T_PropsSI);
        double rel_err_approx = std::abs((drhomolardp__T_PropsSI_num - drhomolardp__T_PropsSI) / drhomolardp__T_PropsSI);
        CHECK(rel_err_exact < 1e-7);
        CHECK(rel_err_approx < 1e-7);
    }
    SECTION("Check drhodp|T 3 ways for water", "") {
        T = 80 + 273.15;
        shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Water"));
        AS->update(CoolProp::PT_INPUTS, 101325, T);

        double drhomolardp__T_AbstractState = AS->first_partial_deriv(CoolProp::iDmolar, CoolProp::iP, CoolProp::iT);
        double drhomolardp__T_PropsSI_num =
          (PropsSI("Dmolar", "T", T, "P", 101325 + 1, "Water") - PropsSI("Dmolar", "T", T, "P", 101325 - 1, "Water")) / (2 * 1);
        double drhomolardp__T_PropsSI = PropsSI("d(Dmolar)/d(P)|T", "T", T, "P", 101325, "Water");

        CAPTURE(drhomolardp__T_AbstractState);
        CAPTURE(drhomolardp__T_PropsSI_num);
        CAPTURE(drhomolardp__T_PropsSI);
        double rel_err_exact = std::abs((drhomolardp__T_AbstractState - drhomolardp__T_PropsSI) / drhomolardp__T_PropsSI);
        double rel_err_approx = std::abs((drhomolardp__T_PropsSI_num - drhomolardp__T_PropsSI) / drhomolardp__T_PropsSI);
        CHECK(rel_err_exact < 1e-4);
        CHECK(rel_err_approx < 1e-4);
    }
    SECTION("Check dpdrho|T 3 ways for water", "") {
        T = 80 + 273.15;
        shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Water"));
        AS->update(CoolProp::PT_INPUTS, 101325, T);
        CoolPropDbl rhomolar = AS->rhomolar();
        double dpdrhomolar__T_AbstractState = AS->first_partial_deriv(CoolProp::iP, CoolProp::iDmolar, CoolProp::iT);
        double dpdrhomolar__T_PropsSI_num =
          (PropsSI("P", "T", T, "Dmolar", rhomolar + 1e-3, "Water") - PropsSI("P", "T", T, "Dmolar", rhomolar - 1e-3, "Water")) / (2 * 1e-3);
        double dpdrhomolar__T_PropsSI = PropsSI("d(P)/d(Dmolar)|T", "T", T, "P", 101325, "Water");
        CAPTURE(rhomolar);
        CAPTURE(dpdrhomolar__T_AbstractState);
        CAPTURE(dpdrhomolar__T_PropsSI_num);
        CAPTURE(dpdrhomolar__T_PropsSI);
        double rel_err_exact = std::abs((dpdrhomolar__T_AbstractState - dpdrhomolar__T_PropsSI) / dpdrhomolar__T_PropsSI);
        double rel_err_approx = std::abs((dpdrhomolar__T_PropsSI_num - dpdrhomolar__T_PropsSI) / dpdrhomolar__T_PropsSI);
        CHECK(rel_err_exact < 1e-6);
        CHECK(rel_err_approx < 1e-6);
    }
    SECTION("Check dpdrho|T 3 ways for water using mass based", "") {
        T = 80 + 273.15;
        shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Water"));
        AS->update(CoolProp::PT_INPUTS, 101325, T);
        CoolPropDbl rhomass = AS->rhomass();
        double dpdrhomass__T_AbstractState = AS->first_partial_deriv(CoolProp::iP, CoolProp::iDmass, CoolProp::iT);
        double dpdrhomass__T_PropsSI_num =
          (PropsSI("P", "T", T, "Dmass", rhomass + 1e-3, "Water") - PropsSI("P", "T", T, "Dmass", rhomass - 1e-3, "Water")) / (2 * 1e-3);
        double dpdrhomass__T_PropsSI = PropsSI("d(P)/d(Dmass)|T", "T", T, "P", 101325, "Water");
        CAPTURE(rhomass);
        CAPTURE(dpdrhomass__T_AbstractState);
        CAPTURE(dpdrhomass__T_PropsSI_num);
        CAPTURE(dpdrhomass__T_PropsSI);
        double rel_err_exact = std::abs((dpdrhomass__T_AbstractState - dpdrhomass__T_PropsSI) / dpdrhomass__T_PropsSI);
        double rel_err_approx = std::abs((dpdrhomass__T_PropsSI_num - dpdrhomass__T_PropsSI) / dpdrhomass__T_PropsSI);
        CHECK(rel_err_exact < 1e-7);
        CHECK(rel_err_approx < 1e-7);
    }
    SECTION("Invalid first partial derivatives", "") {
        CHECK(!ValidNumber(PropsSI("d()/d(P)|T", "T", 300, "P", 101325, "n-Propane")));
        CHECK(!ValidNumber(PropsSI("d(Dmolar)/d()|T", "T", 300, "P", 101325, "n-Propane")));
        CHECK(!ValidNumber(PropsSI("d(Dmolar)/d(P)|", "T", 300, "P", 101325, "n-Propane")));
        CHECK(!ValidNumber(PropsSI("d(XXXX)/d(P)|T", "T", 300, "P", 101325, "n-Propane")));
        CHECK(!ValidNumber(PropsSI("d(Dmolar)d(P)|T", "T", 300, "P", 101325, "n-Propane")));
        CHECK(!ValidNumber(PropsSI("d(Dmolar)/d(P)T", "T", 300, "P", 101325, "n-Propane")));
        CHECK(!ValidNumber(PropsSI("d(Bvirial)/d(P)T", "T", 300, "P", 101325, "n-Propane")));
        CHECK(!ValidNumber(PropsSI("d(Tcrit)/d(P)T", "T", 300, "P", 101325, "n-Propane")));
    }
}

TEST_CASE("Test second partial derivatives", "[derivatives]") {
    double T = 300;
    SECTION("Check d2pdrho2|T 3 ways", "") {
        shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Water"));
        double rhomolar = 60000;
        AS->update(CoolProp::DmolarT_INPUTS, rhomolar, T);
        double p = AS->p();

        double d2pdrhomolar2__T_AbstractState =
          AS->second_partial_deriv(CoolProp::iP, CoolProp::iDmolar, CoolProp::iT, CoolProp::iDmolar, CoolProp::iT);
        // Centered second derivative
        double del = 1e0;
        double d2pdrhomolar2__T_PropsSI_num =
          (PropsSI("P", "T", T, "Dmolar", rhomolar + del, "Water") - 2 * PropsSI("P", "T", T, "Dmolar", rhomolar, "Water")
           + PropsSI("P", "T", T, "Dmolar", rhomolar - del, "Water"))
          / pow(del, 2);
        double d2pdrhomolar2__T_PropsSI = PropsSI("d(d(P)/d(Dmolar)|T)/d(Dmolar)|T", "T", T, "Dmolar", rhomolar, "Water");

        CAPTURE(d2pdrhomolar2__T_AbstractState);
        CAPTURE(d2pdrhomolar2__T_PropsSI_num);
        double rel_err_exact = std::abs((d2pdrhomolar2__T_AbstractState - d2pdrhomolar2__T_PropsSI) / d2pdrhomolar2__T_PropsSI);
        double rel_err_approx = std::abs((d2pdrhomolar2__T_PropsSI_num - d2pdrhomolar2__T_AbstractState) / d2pdrhomolar2__T_AbstractState);
        CHECK(rel_err_exact < 1e-5);
        CHECK(rel_err_approx < 1e-5);
    }
    SECTION("Valid second partial derivatives", "") {
        CHECK(ValidNumber(PropsSI("d(d(Hmolar)/d(P)|T)/d(T)|Dmolar", "T", 300, "P", 101325, "n-Propane")));
    }
    SECTION("Invalid second partial derivatives", "") {
        CHECK(!ValidNumber(PropsSI("d(d()/d(P)|T)/d()|", "T", 300, "P", 101325, "n-Propane")));
        CHECK(!ValidNumber(PropsSI("dd(Dmolar)/d()|T)|T", "T", 300, "P", 101325, "n-Propane")));
    }
    SECTION("Check derivatives with respect to T", "") {
        shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Propane"));
        double rhomolar = 100, dT = 1e-1;
        AS->update(CoolProp::DmolarT_INPUTS, rhomolar, T);

        // base state
        CoolPropDbl T0 = AS->T(), rhomolar0 = AS->rhomolar(), hmolar0 = AS->hmolar(), smolar0 = AS->smolar(), umolar0 = AS->umolar(), p0 = AS->p();
        CoolPropDbl dhdT_rho_ana = AS->first_partial_deriv(CoolProp::iHmolar, CoolProp::iT, CoolProp::iDmolar);
        CoolPropDbl d2hdT2_rho_ana = AS->second_partial_deriv(CoolProp::iHmolar, CoolProp::iT, CoolProp::iDmolar, CoolProp::iT, CoolProp::iDmolar);
        CoolPropDbl dsdT_rho_ana = AS->first_partial_deriv(CoolProp::iSmolar, CoolProp::iT, CoolProp::iDmolar);
        CoolPropDbl d2sdT2_rho_ana = AS->second_partial_deriv(CoolProp::iSmolar, CoolProp::iT, CoolProp::iDmolar, CoolProp::iT, CoolProp::iDmolar);
        CoolPropDbl dudT_rho_ana = AS->first_partial_deriv(CoolProp::iUmolar, CoolProp::iT, CoolProp::iDmolar);
        CoolPropDbl d2udT2_rho_ana = AS->second_partial_deriv(CoolProp::iUmolar, CoolProp::iT, CoolProp::iDmolar, CoolProp::iT, CoolProp::iDmolar);
        CoolPropDbl dpdT_rho_ana = AS->first_partial_deriv(CoolProp::iP, CoolProp::iT, CoolProp::iDmolar);
        CoolPropDbl d2pdT2_rho_ana = AS->second_partial_deriv(CoolProp::iP, CoolProp::iT, CoolProp::iDmolar, CoolProp::iT, CoolProp::iDmolar);

        // increment T
        AS->update(CoolProp::DmolarT_INPUTS, rhomolar, T + dT);
        CoolPropDbl Tpt = AS->T(), rhomolarpt = AS->rhomolar(), hmolarpt = AS->hmolar(), smolarpt = AS->smolar(), umolarpt = AS->umolar(),
                    ppt = AS->p();
        // decrement T
        AS->update(CoolProp::DmolarT_INPUTS, rhomolar, T - dT);
        CoolPropDbl Tmt = AS->T(), rhomolarmt = AS->rhomolar(), hmolarmt = AS->hmolar(), smolarmt = AS->smolar(), umolarmt = AS->umolar(),
                    pmt = AS->p();

        CoolPropDbl dhdT_rho_num = (hmolarpt - hmolarmt) / (2 * dT);
        CoolPropDbl d2hdT2_rho_num = (hmolarpt - 2 * hmolar0 + hmolarmt) / pow(dT, 2);
        CoolPropDbl dsdT_rho_num = (smolarpt - smolarmt) / (2 * dT);
        CoolPropDbl d2sdT2_rho_num = (smolarpt - 2 * smolar0 + smolarmt) / pow(dT, 2);
        CoolPropDbl dudT_rho_num = (umolarpt - umolarmt) / (2 * dT);
        CoolPropDbl d2udT2_rho_num = (umolarpt - 2 * umolar0 + umolarmt) / pow(dT, 2);
        CoolPropDbl dpdT_rho_num = (ppt - pmt) / (2 * dT);
        CoolPropDbl d2pdT2_rho_num = (ppt - 2 * p0 + pmt) / pow(dT, 2);

        CAPTURE(format("%0.15Lg", d2pdT2_rho_ana).c_str());

        double tol = 1e-4;
        CHECK(std::abs((dhdT_rho_num - dhdT_rho_ana) / dhdT_rho_ana) < tol);
        CHECK(std::abs((d2hdT2_rho_num - d2hdT2_rho_ana) / d2hdT2_rho_ana) < tol);
        CHECK(std::abs((dpdT_rho_num - dpdT_rho_ana) / dpdT_rho_ana) < tol);
        CHECK(std::abs((d2pdT2_rho_num - d2pdT2_rho_ana) / d2pdT2_rho_ana) < tol);
        CHECK(std::abs((dsdT_rho_num - dsdT_rho_ana) / dsdT_rho_ana) < tol);
        CHECK(std::abs((d2sdT2_rho_num - d2sdT2_rho_ana) / d2sdT2_rho_ana) < tol);
        CHECK(std::abs((dudT_rho_num - dudT_rho_ana) / dudT_rho_ana) < tol);
        CHECK(std::abs((d2udT2_rho_num - d2udT2_rho_ana) / d2udT2_rho_ana) < tol);
    }

    SECTION("Check derivatives with respect to rho", "") {
        shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Propane"));
        double rhomolar = 100, drho = 1e-1;
        AS->update(CoolProp::DmolarT_INPUTS, rhomolar, T);

        // base state
        CoolPropDbl T0 = AS->T(), rhomolar0 = AS->rhomolar(), hmolar0 = AS->hmolar(), smolar0 = AS->smolar(), umolar0 = AS->umolar(), p0 = AS->p();
        CoolPropDbl dhdrho_T_ana = AS->first_partial_deriv(CoolProp::iHmolar, CoolProp::iDmolar, CoolProp::iT);
        CoolPropDbl d2hdrho2_T_ana = AS->second_partial_deriv(CoolProp::iHmolar, CoolProp::iDmolar, CoolProp::iT, CoolProp::iDmolar, CoolProp::iT);
        CoolPropDbl dsdrho_T_ana = AS->first_partial_deriv(CoolProp::iSmolar, CoolProp::iDmolar, CoolProp::iT);
        CoolPropDbl d2sdrho2_T_ana = AS->second_partial_deriv(CoolProp::iSmolar, CoolProp::iDmolar, CoolProp::iT, CoolProp::iDmolar, CoolProp::iT);
        CoolPropDbl dudrho_T_ana = AS->first_partial_deriv(CoolProp::iUmolar, CoolProp::iDmolar, CoolProp::iT);
        CoolPropDbl d2udrho2_T_ana = AS->second_partial_deriv(CoolProp::iUmolar, CoolProp::iDmolar, CoolProp::iT, CoolProp::iDmolar, CoolProp::iT);
        CoolPropDbl dpdrho_T_ana = AS->first_partial_deriv(CoolProp::iP, CoolProp::iDmolar, CoolProp::iT);
        CoolPropDbl d2pdrho2_T_ana = AS->second_partial_deriv(CoolProp::iP, CoolProp::iDmolar, CoolProp::iT, CoolProp::iDmolar, CoolProp::iT);

        // increment rho
        AS->update(CoolProp::DmolarT_INPUTS, rhomolar + drho, T);
        CoolPropDbl Tpr = AS->T(), rhomolarpr = AS->rhomolar(), hmolarpr = AS->hmolar(), smolarpr = AS->smolar(), umolarpr = AS->umolar(),
                    ppr = AS->p();
        // decrement rho
        AS->update(CoolProp::DmolarT_INPUTS, rhomolar - drho, T);
        CoolPropDbl Tmr = AS->T(), rhomolarmr = AS->rhomolar(), hmolarmr = AS->hmolar(), smolarmr = AS->smolar(), umolarmr = AS->umolar(),
                    pmr = AS->p();

        CoolPropDbl dhdrho_T_num = (hmolarpr - hmolarmr) / (2 * drho);
        CoolPropDbl d2hdrho2_T_num = (hmolarpr - 2 * hmolar0 + hmolarmr) / pow(drho, 2);
        CoolPropDbl dsdrho_T_num = (smolarpr - smolarmr) / (2 * drho);
        CoolPropDbl d2sdrho2_T_num = (smolarpr - 2 * smolar0 + smolarmr) / pow(drho, 2);
        CoolPropDbl dudrho_T_num = (umolarpr - umolarmr) / (2 * drho);
        CoolPropDbl d2udrho2_T_num = (umolarpr - 2 * umolar0 + umolarmr) / pow(drho, 2);
        CoolPropDbl dpdrho_T_num = (ppr - pmr) / (2 * drho);
        CoolPropDbl d2pdrho2_T_num = (ppr - 2 * p0 + pmr) / pow(drho, 2);

        CAPTURE(format("%0.15Lg", d2pdrho2_T_ana).c_str());

        double tol = 1e-4;
        CHECK(std::abs((dhdrho_T_num - dhdrho_T_ana) / dhdrho_T_ana) < tol);
        CHECK(std::abs((d2hdrho2_T_num - d2hdrho2_T_ana) / d2hdrho2_T_ana) < tol);
        CHECK(std::abs((dpdrho_T_num - dpdrho_T_ana) / dpdrho_T_ana) < tol);
        CHECK(std::abs((d2pdrho2_T_num - d2pdrho2_T_ana) / d2pdrho2_T_ana) < tol);
        CHECK(std::abs((dsdrho_T_num - dsdrho_T_ana) / dsdrho_T_ana) < tol);
        CHECK(std::abs((d2sdrho2_T_num - d2sdrho2_T_ana) / d2sdrho2_T_ana) < tol);
        CHECK(std::abs((dudrho_T_num - dudrho_T_ana) / dudrho_T_ana) < tol);
        CHECK(std::abs((d2udrho2_T_num - d2udrho2_T_ana) / d2udrho2_T_ana) < tol);
    }
    SECTION("Check second mixed partial(h,p) with respect to rho", "") {
        shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Propane"));
        double dhmass = 1.0, T = 300;
        AS->update(CoolProp::QT_INPUTS, 0.0, T);
        double deriv1 = AS->first_partial_deriv(iDmass, iP, iHmass);
        double deriv_analyt = AS->second_partial_deriv(iDmass, iP, iHmass, iHmass, iP);
        double deriv_analyt2 = AS->second_partial_deriv(iDmass, iHmass, iP, iP, iHmass);
        AS->update(CoolProp::HmassP_INPUTS, AS->hmass() - 1, AS->p());
        double deriv2 = AS->first_partial_deriv(iDmass, iP, iHmass);
        double deriv_num = (deriv1 - deriv2) / dhmass;
        CAPTURE(deriv_num);
        CAPTURE(deriv_analyt);

        double tol = 1e-4;
        CHECK(std::abs((deriv_num - deriv_analyt) / deriv_analyt) < tol);
    }
}

TEST_CASE("REFPROP names for coolprop fluids", "[REFPROPName]") {
    std::vector<std::string> fluids = strsplit(CoolProp::get_global_param_string("fluids_list"), ',');
    for (std::size_t i = 0; i < fluids.size(); ++i) {
        std::ostringstream ss1;
        ss1 << "Check that REFPROP fluid name for fluid " << fluids[i] << " is valid";
        SECTION(ss1.str(), "") {
            std::string RPName = get_fluid_param_string(fluids[i], "REFPROPName");
            CHECK(!RPName.empty());
            CAPTURE(RPName);
            if (!RPName.compare("N/A")) {
                break;
            }
            CHECK(ValidNumber(Props1SI("REFPROP::" + RPName, "molemass")));
            CHECK(ValidNumber(Props1SI(RPName, "molemass")));
        }
    }
}
TEST_CASE("Backwards compatibility for REFPROP v4 fluid name convention", "[REFPROP_backwards_compatibility]") {
    SECTION("REFPROP-", "") {
        double val = Props1SI("REFPROP-Water", "Tcrit");
        std::string err = get_global_param_string("errstring");
        CAPTURE(val);
        CAPTURE(err);
        CHECK(ValidNumber(val));
    }
    SECTION("REFPROP-MIX:", "") {
        double val = PropsSI("T", "P", 101325, "Q", 0, "REFPROP-MIX:Methane[0.5]&Ethane[0.5]");
        std::string err = get_global_param_string("errstring");
        CAPTURE(val);
        CAPTURE(err);
        CHECK(ValidNumber(val));
    }
}

class AncillaryFixture
{
   public:
    std::string name;
    void run_checks() {
        std::vector<std::string> fluids = strsplit(CoolProp::get_global_param_string("fluids_list"), ',');
        for (std::size_t i = 0; i < fluids.size(); ++i) {
            shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", fluids[i]));
            do_sat(AS);
        }
    }
    void do_sat(shared_ptr<CoolProp::AbstractState>& AS) {
        for (double f = 0.1; f < 1; f += 0.4) {
            double Tc = AS->T_critical();
            double Tt = AS->Ttriple();
            double T = f * Tc + (1 - f) * Tt;
            name = strjoin(AS->fluid_names(), "&");

            AS->update(CoolProp::QT_INPUTS, 0, T);
            check_rhoL(AS);
            check_pL(AS);

            AS->update(CoolProp::QT_INPUTS, 1, T);
            check_rhoV(AS);
            check_pV(AS);
        }
    }
    void check_pL(const shared_ptr<CoolProp::AbstractState>& AS) {
        double p_EOS = AS->saturated_liquid_keyed_output(iP);
        double p_anc = AS->saturation_ancillary(CoolProp::iP, 0, CoolProp::iT, AS->T());
        double err = std::abs(p_EOS - p_anc) / p_anc;
        CAPTURE(name);
        CAPTURE("pL");
        CAPTURE(p_EOS);
        CAPTURE(p_anc);
        CAPTURE(AS->T());
        CHECK(err < 0.02);
    }
    void check_pV(const shared_ptr<CoolProp::AbstractState>& AS) {
        double p_EOS = AS->saturated_liquid_keyed_output(iP);
        double p_anc = AS->saturation_ancillary(CoolProp::iP, 1, CoolProp::iT, AS->T());
        double err = std::abs(p_EOS - p_anc) / p_anc;
        CAPTURE(name);
        CAPTURE("pV");
        CAPTURE(p_EOS);
        CAPTURE(p_anc);
        CAPTURE(AS->T());
        CHECK(err < 0.02);
    }
    void check_rhoL(const shared_ptr<CoolProp::AbstractState>& AS) {
        double rho_EOS = AS->saturated_liquid_keyed_output(iDmolar);
        double rho_anc = AS->saturation_ancillary(CoolProp::iDmolar, 0, CoolProp::iT, AS->T());
        double err = std::abs(rho_EOS - rho_anc) / rho_anc;
        CAPTURE("rhoL");
        CAPTURE(name);
        CAPTURE(rho_EOS);
        CAPTURE(rho_anc);
        CAPTURE(AS->T());
        CHECK(err < 0.03);
    }
    void check_rhoV(const shared_ptr<CoolProp::AbstractState>& AS) {
        double rho_EOS = AS->saturated_vapor_keyed_output(iDmolar);
        double rho_anc = AS->saturation_ancillary(CoolProp::iDmolar, 1, CoolProp::iT, AS->T());
        double err = std::abs(rho_EOS - rho_anc) / rho_anc;
        CAPTURE("rhoV");
        CAPTURE(name);
        CAPTURE(rho_EOS);
        CAPTURE(rho_anc);
        CAPTURE(AS->T());
        CHECK(err < 0.03);
    }
};
TEST_CASE_METHOD(AncillaryFixture, "Ancillary functions", "[ancillary]") {
    run_checks();
};

TEST_CASE("Triple point checks", "[triple_point]") {
    std::vector<std::string> fluids = strsplit(CoolProp::get_global_param_string("fluids_list"), ',');
    for (std::size_t i = 0; i < fluids.size(); ++i) {
        std::vector<std::string> names(1, fluids[i]);
        shared_ptr<CoolProp::HelmholtzEOSMixtureBackend> HEOS(new CoolProp::HelmholtzEOSMixtureBackend(names));
        // Skip pseudo-pure
        if (!HEOS->is_pure()) {
            continue;
        }

        std::ostringstream ss1;
        ss1 << "Minimum saturation temperature state matches for liquid " << fluids[i];
        SECTION(ss1.str(), "") {
            REQUIRE_NOTHROW(HEOS->update(CoolProp::QT_INPUTS, 0, HEOS->Ttriple()));
            double p_EOS = HEOS->p();
            double p_sat_min_liquid = HEOS->get_components()[0].EOS().sat_min_liquid.p;
            double err_sat_min_liquid = std::abs(p_EOS - p_sat_min_liquid) / p_sat_min_liquid;
            CAPTURE(p_EOS);
            CAPTURE(p_sat_min_liquid);
            CAPTURE(err_sat_min_liquid);
            if (p_EOS < 1e-3) {
                continue;
            }  // Skip very low pressure below 1 mPa
            CHECK(err_sat_min_liquid < 1e-3);
        }
        std::ostringstream ss2;
        ss2 << "Minimum saturation temperature state matches for vapor " << fluids[i];
        SECTION(ss2.str(), "") {
            REQUIRE_NOTHROW(HEOS->update(CoolProp::QT_INPUTS, 1, HEOS->Ttriple()));

            double p_EOS = HEOS->p();
            double p_sat_min_vapor = HEOS->get_components()[0].EOS().sat_min_vapor.p;
            double err_sat_min_vapor = std::abs(p_EOS - p_sat_min_vapor) / p_sat_min_vapor;
            CAPTURE(p_EOS);
            CAPTURE(p_sat_min_vapor);
            CAPTURE(err_sat_min_vapor);
            if (p_EOS < 1e-3) {
                continue;
            }  // Skip very low pressure below 1 mPa
            CHECK(err_sat_min_vapor < 1e-3);
        }
        std::ostringstream ss3;
        ss3 << "Minimum saturation temperature state matches for vapor " << fluids[i];
        SECTION(ss3.str(), "") {
            REQUIRE_NOTHROW(HEOS->update(CoolProp::PQ_INPUTS, HEOS->p_triple(), 1));

            double T_EOS = HEOS->T();
            double T_sat_min_vapor = HEOS->get_components()[0].EOS().sat_min_vapor.T;
            double err_sat_min_vapor = std::abs(T_EOS - T_sat_min_vapor);
            CAPTURE(T_EOS);
            CAPTURE(T_sat_min_vapor);
            CAPTURE(err_sat_min_vapor);
            CHECK(err_sat_min_vapor < 1e-3);
        }
        std::ostringstream ss4;
        ss4 << "Minimum saturation temperature state matches for liquid " << fluids[i];
        SECTION(ss4.str(), "") {
            REQUIRE_NOTHROW(HEOS->update(CoolProp::PQ_INPUTS, HEOS->p_triple(), 0));
            double T_EOS = HEOS->T();
            double T_sat_min_vapor = HEOS->get_components()[0].EOS().sat_min_vapor.T;
            double err_sat_min_vapor = std::abs(T_EOS - T_sat_min_vapor);
            CAPTURE(T_EOS);
            CAPTURE(T_sat_min_vapor);
            CAPTURE(err_sat_min_vapor);
            CHECK(err_sat_min_vapor < 1e-3);
        }
        //        std::ostringstream ss2;
        //        ss2 << "Liquid density error < 3% for fluid " << fluids[i] << " at " << T << " K";
        //        SECTION(ss2.str(), "")
        //        {
        //            double rho_EOS = AS->rhomolar();
        //            double rho_anc = AS->saturation_ancillary(CoolProp::iDmolar, 0, CoolProp::iT, T);
        //            double err = std::abs(rho_EOS-rho_anc)/rho_anc;
        //            CAPTURE(rho_EOS);
        //            CAPTURE(rho_anc);
        //            CAPTURE(T);
        //            CHECK(err < 0.03);
        //        }
        //        std::ostringstream ss3;
        //        ss3 << "Vapor density error < 3% for fluid " << fluids[i] << " at " << T << " K";
        //        SECTION(ss3.str(), "")
        //        {
        //            double rho_EOS = AS->rhomolar();
        //            double rho_anc = AS->saturation_ancillary(CoolProp::iDmolar, 1, CoolProp::iT, T);
        //            double err = std::abs(rho_EOS-rho_anc)/rho_anc;
        //            CAPTURE(rho_EOS);
        //            CAPTURE(rho_anc);
        //            CAPTURE(T);
        //            CHECK(err < 0.03);
        //        }
    }
}

class SatTFixture
{
   public:
    std::string name;
    double Tc;
    void run_checks() {
        std::vector<std::string> fluids = strsplit(CoolProp::get_global_param_string("fluids_list"), ',');
        for (std::size_t i = 0; i < fluids.size(); ++i) {
            shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", fluids[i]));
            do_sat(AS);
        }
    }
    void do_sat(shared_ptr<CoolProp::AbstractState>& AS) {
        Tc = AS->T_critical();
        name = strjoin(AS->fluid_names(), "&");
        check_at_Tc(AS);
        double Tt = AS->Ttriple();
        if (AS->fluid_param_string("pure") == "true") {
            Tc = std::min(Tc, AS->T_reducing());
        }
        for (double j = 0.1; j > 1e-10; j /= 10) {
            check_QT(AS, Tc - j);
        }
    }
    void check_at_Tc(const shared_ptr<CoolProp::AbstractState>& AS) {
        CAPTURE("Check @ Tc");
        CAPTURE(name);
        CHECK_NOTHROW(AS->update(QT_INPUTS, 0, Tc));
    }
    void check_QT(const shared_ptr<CoolProp::AbstractState>& AS, double T) {
        std::string test_name = "Check --> Tc";
        CAPTURE(test_name);
        CAPTURE(name);
        CAPTURE(T);
        CHECK_NOTHROW(AS->update(QT_INPUTS, 0, T));
    }
};
TEST_CASE_METHOD(SatTFixture, "Test that saturation solvers solve all the way to T = Tc", "[sat_T_to_Tc]") {
    run_checks();
};

TEST_CASE("Check mixtures with fluid name aliases", "[mixture_name_aliasing]") {
    shared_ptr<CoolProp::AbstractState> AS1, AS2;
    AS1.reset(CoolProp::AbstractState::factory("HEOS", "EBENZENE&P-XYLENE"));
    AS2.reset(CoolProp::AbstractState::factory("HEOS", "EthylBenzene&P-XYLENE"));
    REQUIRE(AS1->fluid_names().size() == AS2->fluid_names().size());
    std::size_t N = AS1->fluid_names().size();
    for (std::size_t i = 0; i < N; ++i) {
        CAPTURE(i);
        CHECK(AS1->fluid_names()[i] == AS2->fluid_names()[i]);
    }
}

TEST_CASE("Predefined mixtures", "[predefined_mixtures]") {
    SECTION("PropsSI") {
        double val = PropsSI("Dmolar", "P", 101325, "T", 300, "Air.mix");
        std::string err = get_global_param_string("errstring");
        CAPTURE(val);
        CAPTURE(err);
        CHECK(ValidNumber(val));
    }
}
TEST_CASE("Test that reference states yield proper values using high-level interface", "[reference_states]") {
    struct ref_entry
    {
        std::string name;
        double hmass, smass;
        std::string in1;
        double val1;
        std::string in2;
        double val2;
    };
    std::string fluids[] = {"n-Propane", "R134a", "R124"};
    ref_entry entries[3] = {{"IIR", 200000, 1000, "T", 273.15, "Q", 0}, {"ASHRAE", 0, 0, "T", 233.15, "Q", 0}, {"NBP", 0, 0, "P", 101325, "Q", 0}};
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            std::ostringstream ss1;
            ss1 << "Check state for " << fluids[i] << " for " + entries[j].name + " reference state ";
            SECTION(ss1.str(), "") {
                // First reset the reference state
                set_reference_stateS(fluids[i], "DEF");
                // Then set to desired reference state
                set_reference_stateS(fluids[i], entries[j].name);
                // Calculate the values
                double hmass = PropsSI("Hmass", entries[j].in1, entries[j].val1, entries[j].in2, entries[j].val2, fluids[i]);
                double smass = PropsSI("Smass", entries[j].in1, entries[j].val1, entries[j].in2, entries[j].val2, fluids[i]);
                CHECK(std::abs(hmass - entries[j].hmass) < 1e-8);
                CHECK(std::abs(smass - entries[j].smass) < 1e-8);
                // Then reset the reference state
                set_reference_stateS(fluids[i], "DEF");
            }
        }
    }
}
TEST_CASE("Test that reference states yield proper values using low-level interface", "[reference_states]") {
    struct ref_entry
    {
        std::string name;
        double hmass, smass;
        parameters in1;
        double val1;
        parameters in2;
        double val2;
    };
    std::string fluids[] = {"n-Propane", "R134a", "R124"};
    ref_entry entries[3] = {{"IIR", 200000, 1000, iT, 273.15, iQ, 0}, {"ASHRAE", 0, 0, iT, 233.15, iQ, 0}, {"NBP", 0, 0, iP, 101325, iQ, 0}};
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            std::ostringstream ss1;
            ss1 << "Check state for " << fluids[i] << " for " + entries[j].name + " reference state ";
            SECTION(ss1.str(), "") {
                double val1, val2;
                input_pairs pair = generate_update_pair(entries[j].in1, entries[j].val1, entries[j].in2, entries[j].val2, val1, val2);
                // Generate a state instance
                shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", fluids[i]));
                AS->update(pair, val1, val2);
                double hmass0 = AS->hmass();
                double smass0 = AS->smass();
                // First reset the reference state
                set_reference_stateS(fluids[i], "DEF");
                AS->update(pair, val1, val2);
                double hmass00 = AS->hmass();
                double smass00 = AS->smass();
                CHECK(std::abs(hmass00 - hmass0) < 1e-10);
                CHECK(std::abs(smass00 - smass0) < 1e-10);

                // Then set to desired reference state
                set_reference_stateS(fluids[i], entries[j].name);

                // Should not change existing instance
                AS->clear();
                AS->update(pair, val1, val2);
                double hmass1 = AS->hmass();
                double smass1 = AS->smass();
                CHECK(std::abs(hmass1 - hmass0) < 1e-10);
                CHECK(std::abs(smass1 - smass0) < 1e-10);

                // New instance - should get updated reference state
                shared_ptr<CoolProp::AbstractState> AS2(CoolProp::AbstractState::factory("HEOS", fluids[i]));
                AS2->update(pair, val1, val2);
                double hmass2 = AS2->hmass();
                double smass2 = AS2->smass();
                CHECK(std::abs(hmass2 - entries[j].hmass) < 1e-8);
                CHECK(std::abs(smass2 - entries[j].smass) < 1e-8);

                // Then reset the reference state
                set_reference_stateS(fluids[i], "DEF");
            }
        }
    }
}

class FixedStateFixture
{
   public:
    void run_fluid(const std::string& fluid, const std::string& state, const std::string& ref_state) {

        // Skip impossible reference states
        if (Props1SI("Ttriple", fluid) > 233.15 && ref_state == "ASHRAE") {
            return;
        }
        if (Props1SI("Tcrit", fluid) < 233.15 && ref_state == "ASHRAE") {
            return;
        }
        if (Props1SI("Tcrit", fluid) < 273.15 && ref_state == "IIR") {
            return;
        }
        if (Props1SI("Ttriple", fluid) > 273.15 && ref_state == "IIR") {
            return;
        }
        if (Props1SI("ptriple", fluid) > 101325 && ref_state == "NBP") {
            return;
        }

        // First reset the reference state
        if (ref_state != "DEF") {
            set_reference_stateS(fluid, "DEF");
            try {
                // Then try to set to the specified reference state
                set_reference_stateS(fluid, ref_state);
            } catch (std::exception& e) {
                // Then set the reference state back to the default
                set_reference_stateS(fluid, "DEF");
                CAPTURE(e.what());
                REQUIRE(false);
            }
        }

        std::ostringstream name;
        name << "Check state for " << state << " for " << fluid << " for reference state " << ref_state;
        CAPTURE(name.str());

        std::vector<std::string> fl(1, fluid);
        shared_ptr<CoolProp::HelmholtzEOSMixtureBackend> HEOS(new CoolProp::HelmholtzEOSMixtureBackend(fl));

        // Skip the saturation maxima states for pure fluids
        if (HEOS->is_pure() && (state == "max_sat_T" || state == "max_sat_p")) {
            return;
        }

        // Get the state
        CoolProp::SimpleState _state = HEOS->calc_state(state);
        HEOS->specify_phase(iphase_gas);  // something homogenous
        // Bump a tiny bit for EOS with non-analytic parts
        double f = 1.0;
        if ((fluid == "Water" || fluid == "CarbonDioxide") && (state == "reducing" || state == "critical")) {
            f = 1.00001;
        }
        HEOS->update(CoolProp::DmolarT_INPUTS, _state.rhomolar * f, _state.T * f);
        CAPTURE(_state.hmolar);
        CAPTURE(_state.smolar);
        CHECK(ValidNumber(_state.hmolar));
        CHECK(ValidNumber(_state.smolar));
        double EOS_hmolar = HEOS->hmolar();
        double EOS_smolar = HEOS->smolar();
        CAPTURE(EOS_hmolar);
        CAPTURE(EOS_smolar);
        CHECK(std::abs(EOS_hmolar - _state.hmolar) < 1e-2);
        CHECK(std::abs(EOS_smolar - _state.smolar) < 1e-2);
        // Then set the reference state back to the default
        set_reference_stateS(fluid, "DEF");
    };
    void run_checks() {

        std::vector<std::string> fluids = strsplit(CoolProp::get_global_param_string("fluids_list"), ',');
        for (std::size_t i = 0; i < fluids.size(); ++i) {
            std::string ref_state[4] = {"DEF", "IIR", "ASHRAE", "NBP"};
            for (std::size_t j = 0; j < 4; ++j) {
                std::string states[] = {"hs_anchor", "reducing", "critical", "max_sat_T", "max_sat_p", "triple_liquid", "triple_vapor"};
                for (std::size_t k = 0; k < 7; ++k) {
                    run_fluid(fluids[i], states[k], ref_state[j]);
                }
            }
        }
    }
};
TEST_CASE_METHOD(FixedStateFixture, "Test that enthalpies and entropies are correct for fixed states for all reference states", "[fixed_states]") {
    run_checks();
};  // !!!! check this

TEST_CASE("Check the first partial derivatives", "[first_saturation_partial_deriv]") {
    const int number_of_pairs = 10;
    struct pair
    {
        parameters p1, p2;
    };
    pair pairs[number_of_pairs] = {{iP, iT}, {iDmolar, iT}, {iHmolar, iT}, {iSmolar, iT}, {iUmolar, iT},
                                   {iT, iP}, {iDmolar, iP}, {iHmolar, iP}, {iSmolar, iP}, {iUmolar, iP}};
    shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "n-Propane"));
    for (std::size_t i = 0; i < number_of_pairs; ++i) {
        // See https://groups.google.com/forum/?fromgroups#!topic/catch-forum/mRBKqtTrITU
        std::ostringstream ss1;
        ss1 << "Check first partial derivative for d(" << get_parameter_information(pairs[i].p1, "short") << ")/d("
            << get_parameter_information(pairs[i].p2, "short") << ")|sat";
        SECTION(ss1.str(), "") {
            AS->update(QT_INPUTS, 1, 300);
            CoolPropDbl p = AS->p();
            CoolPropDbl analytical = AS->first_saturation_deriv(pairs[i].p1, pairs[i].p2);
            CAPTURE(analytical);
            CoolPropDbl numerical;
            if (pairs[i].p2 == iT) {
                AS->update(QT_INPUTS, 1, 300 + 1e-5);
                CoolPropDbl v1 = AS->keyed_output(pairs[i].p1);
                AS->update(QT_INPUTS, 1, 300 - 1e-5);
                CoolPropDbl v2 = AS->keyed_output(pairs[i].p1);
                numerical = (v1 - v2) / (2e-5);
            } else if (pairs[i].p2 == iP) {
                AS->update(PQ_INPUTS, p + 1e-2, 1);
                CoolPropDbl v1 = AS->keyed_output(pairs[i].p1);
                AS->update(PQ_INPUTS, p - 1e-2, 1);
                CoolPropDbl v2 = AS->keyed_output(pairs[i].p1);
                numerical = (v1 - v2) / (2e-2);
            } else {
                throw ValueError();
            }
            CAPTURE(numerical);
            CHECK(std::abs(numerical / analytical - 1) < 1e-4);
        }
    }
}

TEST_CASE("Check the second saturation derivatives", "[second_saturation_partial_deriv]") {
    const int number_of_pairs = 5;
    struct pair
    {
        parameters p1, p2, p3;
    };
    pair pairs[number_of_pairs] = {{iT, iP, iP}, {iDmolar, iP, iP}, {iHmolar, iP, iP}, {iSmolar, iP, iP}, {iUmolar, iP, iP}};
    shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "n-Propane"));
    for (std::size_t i = 0; i < number_of_pairs; ++i) {
        // See https://groups.google.com/forum/?fromgroups#!topic/catch-forum/mRBKqtTrITU
        std::ostringstream ss1;
        ss1 << "Check second saturation derivative for d2(" << get_parameter_information(pairs[i].p1, "short") << ")/d("
            << get_parameter_information(pairs[i].p2, "short") << ")2|sat";
        SECTION(ss1.str(), "") {
            AS->update(QT_INPUTS, 1, 300);
            CoolPropDbl p = AS->p();
            CoolPropDbl analytical = AS->second_saturation_deriv(pairs[i].p1, pairs[i].p2, pairs[i].p3);
            CAPTURE(analytical);
            CoolPropDbl numerical;
            if (pairs[i].p2 == iT) {
                throw NotImplementedError();
            } else if (pairs[i].p2 == iP) {
                AS->update(PQ_INPUTS, p + 1e-2, 1);
                CoolPropDbl v1 = AS->first_saturation_deriv(pairs[i].p1, pairs[i].p2);
                AS->update(PQ_INPUTS, p - 1e-2, 1);
                CoolPropDbl v2 = AS->first_saturation_deriv(pairs[i].p1, pairs[i].p2);
                numerical = (v1 - v2) / (2e-2);
            } else {
                throw ValueError();
            }
            CAPTURE(numerical);
            CHECK(std::abs(numerical / analytical - 1) < 1e-4);
        }
    }
}

TEST_CASE("Check the first two-phase derivative", "[first_two_phase_deriv]") {
    const int number_of_pairs = 4;
    struct pair
    {
        parameters p1, p2, p3;
    };
    pair pairs[number_of_pairs] = {{iDmass, iP, iHmass}, {iDmolar, iP, iHmolar}, {iDmolar, iHmolar, iP}, {iDmass, iHmass, iP}};
    shared_ptr<CoolProp::HelmholtzEOSBackend> AS(new CoolProp::HelmholtzEOSBackend("n-Propane"));
    for (std::size_t i = 0; i < number_of_pairs; ++i) {
        // See https://groups.google.com/forum/?fromgroups#!topic/catch-forum/mRBKqtTrITU
        std::ostringstream ss1;
        ss1 << "for (" << get_parameter_information(pairs[i].p1, "short") << ", " << get_parameter_information(pairs[i].p2, "short") << ", "
            << get_parameter_information(pairs[i].p3, "short") << ")";
        SECTION(ss1.str(), "") {
            AS->update(QT_INPUTS, 0.3, 300);
            CoolPropDbl numerical;
            CoolPropDbl analytical = AS->first_two_phase_deriv(pairs[i].p1, pairs[i].p2, pairs[i].p3);
            CAPTURE(analytical);

            CoolPropDbl out1, out2;
            CoolPropDbl v2base, v3base;
            v2base = AS->keyed_output(pairs[i].p2);
            v3base = AS->keyed_output(pairs[i].p3);
            CoolPropDbl v2plus = v2base * 1.001;
            CoolPropDbl v2minus = v2base * 0.999;
            CoolProp::input_pairs input_pair1 = generate_update_pair(pairs[i].p2, v2plus, pairs[i].p3, v3base, out1, out2);
            AS->update(input_pair1, out1, out2);
            CoolPropDbl v1 = AS->keyed_output(pairs[i].p1);
            CoolProp::input_pairs input_pair2 = generate_update_pair(pairs[i].p2, v2minus, pairs[i].p3, v3base, out1, out2);
            AS->update(input_pair2, out1, out2);
            CoolPropDbl v2 = AS->keyed_output(pairs[i].p1);

            numerical = (v1 - v2) / (v2plus - v2minus);
            CAPTURE(numerical);
            CHECK(std::abs(numerical / analytical - 1) < 1e-4);
        }
    }
}

TEST_CASE("Check the second two-phase derivative", "[second_two_phase_deriv]") {
    SECTION("d2rhodhdp", "") {
        shared_ptr<CoolProp::HelmholtzEOSBackend> AS(new CoolProp::HelmholtzEOSBackend("n-Propane"));
        AS->update(QT_INPUTS, 0.3, 300);
        CoolPropDbl analytical = AS->second_two_phase_deriv(iDmolar, iHmolar, iP, iP, iHmolar);
        CAPTURE(analytical);
        CoolPropDbl pplus = AS->p() * 1.001, pminus = AS->p() * 0.999, h = AS->hmolar();
        AS->update(HmolarP_INPUTS, h, pplus);
        CoolPropDbl v1 = AS->first_two_phase_deriv(iDmolar, iHmolar, iP);
        AS->update(HmolarP_INPUTS, h, pminus);
        CoolPropDbl v2 = AS->first_two_phase_deriv(iDmolar, iHmolar, iP);
        CoolPropDbl numerical = (v1 - v2) / (pplus - pminus);
        CAPTURE(numerical);
        CHECK(std::abs(numerical / analytical - 1) < 1e-6);
    }
    SECTION("d2rhodhdp using mass", "") {
        shared_ptr<CoolProp::HelmholtzEOSBackend> AS(new CoolProp::HelmholtzEOSBackend("n-Propane"));
        AS->update(QT_INPUTS, 0.3, 300);
        CoolPropDbl analytical = AS->second_two_phase_deriv(iDmass, iHmass, iP, iP, iHmass);
        CAPTURE(analytical);
        CoolPropDbl pplus = AS->p() * 1.001, pminus = AS->p() * 0.999, h = AS->hmass();
        AS->update(HmassP_INPUTS, h, pplus);
        CoolPropDbl v1 = AS->first_two_phase_deriv(iDmass, iHmass, iP);
        AS->update(HmassP_INPUTS, h, pminus);
        CoolPropDbl v2 = AS->first_two_phase_deriv(iDmass, iHmass, iP);
        CoolPropDbl numerical = (v1 - v2) / (pplus - pminus);
        CAPTURE(numerical);
        CHECK(std::abs(numerical / analytical - 1) < 1e-6);
    }
}

TEST_CASE("Check the first two-phase derivative using splines", "[first_two_phase_deriv_splined]") {
    const int number_of_pairs = 4;
    struct pair
    {
        parameters p1, p2, p3;
    };
    pair pairs[number_of_pairs] = {{iDmass, iP, iHmass}, {iDmolar, iP, iHmolar}, {iDmolar, iHmolar, iP}, {iDmass, iHmass, iP}};
    shared_ptr<CoolProp::HelmholtzEOSBackend> AS(new CoolProp::HelmholtzEOSBackend("n-Propane"));
    for (std::size_t i = 0; i < number_of_pairs; ++i) {
        // See https://groups.google.com/forum/?fromgroups#!topic/catch-forum/mRBKqtTrITU
        std::ostringstream ss1;
        ss1 << "for (" << get_parameter_information(pairs[i].p1, "short") << ", " << get_parameter_information(pairs[i].p2, "short") << ", "
            << get_parameter_information(pairs[i].p3, "short") << ")";
        SECTION(ss1.str(), "") {
            AS->update(QT_INPUTS, 0.2, 300);
            CoolPropDbl numerical;
            CoolPropDbl analytical = AS->first_two_phase_deriv_splined(pairs[i].p1, pairs[i].p2, pairs[i].p3, 0.3);
            CAPTURE(analytical);

            CoolPropDbl out1, out2;
            CoolPropDbl v2base, v3base;
            v2base = AS->keyed_output(pairs[i].p2);
            v3base = AS->keyed_output(pairs[i].p3);
            CoolPropDbl v2plus = v2base * 1.00001;
            CoolPropDbl v2minus = v2base * 0.99999;

            CoolProp::input_pairs input_pair1 = generate_update_pair(pairs[i].p2, v2plus, pairs[i].p3, v3base, out1, out2);
            AS->update(input_pair1, out1, out2);
            CoolPropDbl v1 = AS->first_two_phase_deriv_splined(pairs[i].p1, pairs[i].p1, pairs[i].p1, 0.3);

            CoolProp::input_pairs input_pair2 = generate_update_pair(pairs[i].p2, v2minus, pairs[i].p3, v3base, out1, out2);
            AS->update(input_pair2, out1, out2);
            CoolPropDbl v2 = AS->first_two_phase_deriv_splined(pairs[i].p1, pairs[i].p1, pairs[i].p1, 0.3);

            numerical = (v1 - v2) / (v2plus - v2minus);
            CAPTURE(numerical);
            CHECK(std::abs(numerical / analytical - 1) < 1e-8);
        }
    }
}

TEST_CASE("Check the phase flags", "[phase]") {
    SECTION("subcooled liquid") {
        shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Water"));
        AS->update(PT_INPUTS, 101325, 300);
        CHECK(AS->phase() == iphase_liquid);
    }
    SECTION("superheated gas") {
        shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Water"));
        AS->update(PT_INPUTS, 101325, 400);
        CHECK(AS->phase() == iphase_gas);
    }
    SECTION("supercritical gas") {
        shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Water"));
        AS->update(PT_INPUTS, 1e5, 800);
        CHECK(AS->phase() == iphase_supercritical_gas);
    }
    SECTION("supercritical liquid") {
        shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Water"));
        AS->update(PT_INPUTS, 1e8, 500);
        CHECK(AS->phase() == iphase_supercritical_liquid);
    }
    SECTION("supercritical") {
        shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Water"));
        AS->update(PT_INPUTS, 1e8, 800);
        CHECK(AS->phase() == iphase_supercritical);
    }
}

TEST_CASE("Check the changing of reducing function constants", "[reducing]") {
    double z0 = 0.2;
    std::vector<double> z(2);
    z[0] = z0;
    z[1] = 1 - z[0];
    shared_ptr<CoolProp::AbstractState> AS1(CoolProp::AbstractState::factory("HEOS", "Methane&Ethane"));
    shared_ptr<CoolProp::AbstractState> AS2(CoolProp::AbstractState::factory("HEOS", "Methane&Ethane"));
    AS1->set_mole_fractions(z);
    AS2->set_mole_fractions(z);
    std::vector<CoolProp::CriticalState> pts1 = AS1->all_critical_points();
    double gammaT = AS2->get_binary_interaction_double(0, 1, "gammaT");
    AS2->set_binary_interaction_double(0, 1, "gammaT", gammaT * 0.7);
    std::vector<CoolProp::CriticalState> pts2 = AS2->all_critical_points();
    double Tdiff = abs(pts2[0].T - pts1[0].T);
    CHECK(Tdiff > 1e-3);  // Make sure that it actually got the change to the interaction parameters
}

TEST_CASE("Check the PC-SAFT pressure function", "[pcsaft_pressure]") {
    double p = 101325.;
    double p_calc = CoolProp::PropsSI("P", "T", 320., "Dmolar", 9033.11420899, "PCSAFT::TOLUENE");
    CHECK(abs((p_calc / p) - 1) < 1e-5);

    p_calc = CoolProp::PropsSI("P", "T", 274., "Dmolar", 55530.40512318346, "PCSAFT::WATER");
    CHECK(abs((p_calc / p) - 1) < 1e-5);

    p_calc = CoolProp::PropsSI("P", "T", 305., "Dmolar", 16965.43663595, "PCSAFT::ACETIC ACID");
    CHECK(abs((p_calc / p) - 1) < 1e-5);

    p_calc = CoolProp::PropsSI("P", "T", 240., "Dmolar", 15865.69021378, "PCSAFT::DIMETHYL ETHER");
    CHECK(abs((p_calc / p) - 1) < 1e-5);

    p_calc = CoolProp::PropsSI("P", "T", 298.15, "Dmolar", 9368.9036823, "PCSAFT::METHANOL[0.055]&CYCLOHEXANE[0.945]");
    CHECK(abs((p_calc / p) - 1) < 1e-5);

    p_calc = CoolProp::PropsSI("P", "T", 298.15, "Dmolar", 55740.157290833515,
                               "PCSAFT::Na+[0.010579869455908]&Cl-[0.010579869455908]&WATER[0.978840261088184]");
    CHECK(abs((p_calc / p) - 1) < 1e-5);

    p = CoolProp::PropsSI("P", "T", 100., "Q", 0, "PCSAFT::PROPANE");
    double rho = 300;
    double phase = CoolProp::PropsSI("Phase", "T", 100., "Dmolar", rho, "PCSAFT::PROPANE");
    CHECK(phase == get_phase_index("phase_twophase"));
    p_calc = CoolProp::PropsSI("P", "T", 100, "Dmolar", rho, "PCSAFT::PROPANE");
    CHECK(abs((p_calc / p) - 1) < 1e-4);
}

TEST_CASE("Check the PC-SAFT density function", "[pcsaft_density]") {
    double den = 9033.114209728405;
    double den_calc = CoolProp::PropsSI("Dmolar", "T|liquid", 320., "P", 101325., "PCSAFT::TOLUENE");
    CHECK(abs((den_calc / den) - 1) < 1e-5);

    den = 55530.40512318346;
    den_calc = CoolProp::PropsSI("Dmolar", "T|liquid", 274., "P", 101325, "PCSAFT::WATER");
    CHECK(abs((den_calc / den) - 1) < 1e-5);

    den = 16965.436637145376;
    den_calc = CoolProp::PropsSI("Dmolar", "T|liquid", 305., "P", 101325, "PCSAFT::ACETIC ACID");
    CHECK(abs((den_calc / den) - 1) < 1e-5);

    den = 15865.690215090615;
    den_calc = CoolProp::PropsSI("Dmolar", "T|liquid", 240., "P", 101325, "PCSAFT::DIMETHYL ETHER");
    CHECK(abs((den_calc / den) - 1) < 1e-5);

    den = 9368.90368306872;
    den_calc = CoolProp::PropsSI("Dmolar", "T|liquid", 298.15, "P", 101325, "PCSAFT::METHANOL[0.055]&CYCLOHEXANE[0.945]");
    CHECK(abs((den_calc / den) - 1) < 1e-5);

    den = 55740.157290833515;
    den_calc =
      CoolProp::PropsSI("Dmolar", "T|liquid", 298.15, "P", 101325, "PCSAFT::Na+[0.010579869455908]&Cl-[0.010579869455908]&WATER[0.978840261088184]");
    CHECK(abs((den_calc / den) - 1) < 1e-5);

    den = 16621.0;
    den_calc = CoolProp::PropsSI("Dmolar", "T|liquid", 85.525, "P", 1.7551e-4, "PCSAFT::PROPANE");
    CHECK(abs((den_calc / den) - 1) < 1e-2);

    den = 1.9547e-7;
    den_calc = CoolProp::PropsSI("Dmolar", "T|gas", 85.525, "P", 1.39e-4, "PCSAFT::PROPANE");
    CHECK(abs((den_calc / den) - 1) < 1e-2);

    den = 11346.0;
    den_calc = CoolProp::PropsSI("Dmolar", "T|liquid", 293, "P", 833240, "PCSAFT::PROPANE");
    CHECK(abs((den_calc / den) - 1) < 1e-2);

    den = 623.59;
    den_calc = CoolProp::PropsSI("Dmolar", "T", 430, "P", 2000000, "PCSAFT::PROPANE");
    CHECK(abs((den_calc / den) - 1) < 1e-2);

    den = 623.59;
    den_calc = CoolProp::PropsSI("Dmolar", "T", 430, "P", 2000000, "PCSAFT::PROPANE");
    CHECK(abs((den_calc / den) - 1) < 1e-2);
}

TEST_CASE("Check the PC-SAFT residual enthalpy function", "[pcsaft_enthalpy]") {
    double h = -36809.962122036086;
    double h_calc = CoolProp::PropsSI("Hmolar_residual", "T|liquid", 325., "Dmolar", 8983.377722763931, "PCSAFT::TOLUENE");
    CHECK(abs((h_calc / h) - 1) < 1e-5);

    h = -362.6832840695562;
    h_calc = CoolProp::PropsSI("Hmolar_residual", "T|gas", 325., "Dmolar", 39.44490805826904, "PCSAFT::TOLUENE");
    CHECK(abs((h_calc / h) - 1) < 1e-5);

    h = -38925.302571456035;
    h_calc = CoolProp::PropsSI("Hmolar_residual", "T|liquid", 325., "Dmolar", 16655.844528563375, "PCSAFT::ACETIC ACID");
    CHECK(abs((h_calc / h) - 1) < 1e-5);

    h = -15393.870073928741;
    h_calc = CoolProp::PropsSI("Hmolar_residual", "T|gas", 325., "Dmolar", 85.70199446609787, "PCSAFT::ACETIC ACID");
    CHECK(abs((h_calc / h) - 1) < 1e-5);

    h = -18037.24422056259;
    h_calc = CoolProp::PropsSI("Hmolar_residual", "T|liquid", 325., "Dmolar", 12963.391139983729, "PCSAFT::DIMETHYL ETHER");
    CHECK(abs((h_calc / h) - 1) < 1e-5);

    h = -92.66136745908202;
    h_calc = CoolProp::PropsSI("Hmolar_residual", "T|gas", 325., "Dmolar", 37.9473393419189, "PCSAFT::DIMETHYL ETHER");
    CHECK(abs((h_calc / h) - 1) < 1e-5);

    // checks based on values from the HEOS backend
    h = CoolProp::PropsSI("Hmolar_residual", "T|liquid", 325., "Dmolar", 8983.377722763931, "HEOS::TOLUENE");
    h_calc = CoolProp::PropsSI("Hmolar_residual", "T|liquid", 325., "Dmolar", 8983.377722763931, "PCSAFT::TOLUENE");
    CHECK(abs(h_calc - h) < 600.);

    h = CoolProp::PropsSI("Hmolar_residual", "T|gas", 325., "Dmolar", 39.44490805826904, "HEOS::TOLUENE");
    h_calc = CoolProp::PropsSI("Hmolar_residual", "T|gas", 325., "Dmolar", 39.44490805826904, "PCSAFT::TOLUENE");
    CHECK(abs(h_calc - h) < 600.);

    h = CoolProp::PropsSI("Hmolar_residual", "T|liquid", 325., "Dmolar", 54794.1, "HEOS::WATER");
    h_calc = CoolProp::PropsSI("Hmolar_residual", "T|liquid", 325., "Dmolar", 54794.1, "PCSAFT::WATER");
    CHECK(abs(h_calc - h) < 600.);

    h = CoolProp::PropsSI("Hmolar_residual", "T|gas", 325., "Dmolar", 0.370207, "HEOS::WATER");
    h_calc = CoolProp::PropsSI("Hmolar_residual", "T|gas", 325., "Dmolar", 0.370207, "PCSAFT::WATER");
    CHECK(abs(h_calc - h) < 600.);
}

TEST_CASE("Check the PC-SAFT residual entropy function", "[pcsaft_entropy]") {
    // checks based on values from working PC-SAFT code
    double s = -50.81694890352192;
    double s_calc = CoolProp::PropsSI("Smolar_residual", "T|liquid", 325., "Dmolar", 8983.377722763931, "PCSAFT::TOLUENE");
    CHECK(abs((s_calc / s) - 1) < 1e-5);

    s = -0.2929618646219797;
    s_calc = CoolProp::PropsSI("Smolar_residual", "T|gas", 325., "Dmolar", 39.44490805826904, "PCSAFT::TOLUENE");
    CHECK(abs((s_calc / s) - 1) < 1e-5);

    s = -47.42736805661422;
    s_calc = CoolProp::PropsSI("Smolar_residual", "T|liquid", 325., "Dmolar", 16655.844528563375, "PCSAFT::ACETIC ACID");
    CHECK(abs((s_calc / s) - 1) < 1e-5);

    s = -34.0021996393859;
    s_calc = CoolProp::PropsSI("Smolar_residual", "T|gas", 325., "Dmolar", 85.70199446609787, "PCSAFT::ACETIC ACID");
    CHECK(abs((s_calc / s) - 1) < 1e-5);

    s = -25.91216157948035;
    s_calc = CoolProp::PropsSI("Smolar_residual", "T|liquid", 325., "Dmolar", 12963.391139983729, "PCSAFT::DIMETHYL ETHER");
    CHECK(abs((s_calc / s) - 1) < 1e-5);

    s = -0.0842409121406476;
    s_calc = CoolProp::PropsSI("Smolar_residual", "T|gas", 325., "Dmolar", 37.9473393419189, "PCSAFT::DIMETHYL ETHER");
    CHECK(abs((s_calc / s) - 1) < 1e-5);

    // checks based on values from the HEOS backend
    s = CoolProp::PropsSI("Smolar_residual", "T|liquid", 325., "Dmolar", 8983.377722763931, "HEOS::TOLUENE");
    s_calc = CoolProp::PropsSI("Smolar_residual", "T|liquid", 325., "Dmolar", 8983.377722763931, "PCSAFT::TOLUENE");
    CHECK(abs(s_calc - s) < 3.);

    s = CoolProp::PropsSI("Smolar_residual", "T|gas", 325., "Dmolar", 39.44490805826904, "HEOS::TOLUENE");
    s_calc = CoolProp::PropsSI("Smolar_residual", "T|gas", 325., "Dmolar", 39.44490805826904, "PCSAFT::TOLUENE");
    CHECK(abs(s_calc - s) < 3.);

    s = CoolProp::PropsSI("Smolar_residual", "T|liquid", 325., "Dmolar", 54794.1, "HEOS::WATER");
    s_calc = CoolProp::PropsSI("Smolar_residual", "T|liquid", 325., "Dmolar", 54794.1, "PCSAFT::WATER");
    CHECK(abs(s_calc - s) < 3.);

    s = CoolProp::PropsSI("Smolar_residual", "T|gas", 325., "Dmolar", 0.370207, "HEOS::WATER");
    s_calc = CoolProp::PropsSI("Smolar_residual", "T|gas", 325., "Dmolar", 0.370207, "PCSAFT::WATER");
    CHECK(abs(s_calc - s) < 3.);
}

TEST_CASE("Check the PC-SAFT residual gibbs energy function", "[pcsaft_gibbs]") {
    double g = -5489.471870270737;
    double g_calc = CoolProp::PropsSI("Gmolar_residual", "T|liquid", 325., "Dmolar", 8983.377722763931, "PCSAFT::TOLUENE");
    CHECK(abs((g_calc / g) - 1) < 1e-5);

    g = -130.63592030187894;
    g_calc = CoolProp::PropsSI("Gmolar_residual", "T|gas", 325., "Dmolar", 39.44490805826904, "PCSAFT::TOLUENE");
    CHECK(abs((g_calc / g) - 1) < 1e-5);

    g = -7038.128334100866;
    g_calc = CoolProp::PropsSI("Gmolar_residual", "T|liquid", 325., "Dmolar", 16655.844528563375, "PCSAFT::ACETIC ACID");
    CHECK(abs((g_calc / g) - 1) < 1e-5);

    g = -2109.4916554917604;
    g_calc = CoolProp::PropsSI("Gmolar_residual", "T|gas", 325., "Dmolar", 85.70199446609787, "PCSAFT::ACETIC ACID");
    CHECK(abs((g_calc / g) - 1) < 1e-5);

    g = 6180.230281553767;
    g_calc = CoolProp::PropsSI("Gmolar_residual", "T|liquid", 325., "Dmolar", 12963.391139983729, "PCSAFT::DIMETHYL ETHER");
    CHECK(abs((g_calc / g) - 1) < 1e-5);

    g = -33.03853932580277;
    g_calc = CoolProp::PropsSI("Gmolar_residual", "T|gas", 325., "Dmolar", 37.9473393419189, "PCSAFT::DIMETHYL ETHER");
    CHECK(abs((g_calc / g) - 1) < 1e-5);
}

TEST_CASE("Check vapor pressures calculated using PC-SAFT", "[pcsaft_vapor_pressure]") {
    double vp = 3290651.18080112;
    double vp_calc = CoolProp::PropsSI("P", "T", 572.6667, "Q", 0, "PCSAFT::TOLUENE");
    CHECK(abs((vp_calc / vp) - 1) < 1e-3);

    vp = 66917.67387203;
    vp_calc = CoolProp::PropsSI("P", "T", 362, "Q", 0, "PCSAFT::WATER");
    CHECK(abs((vp_calc / vp) - 1) < 1e-3);

    vp = 190061.78088909;
    vp_calc = CoolProp::PropsSI("P", "T", 413.5385, "Q", 0, "PCSAFT::ACETIC ACID");
    CHECK(abs((vp_calc / vp) - 1) < 1e-3);

    vp = 623027.07850612;
    vp_calc = CoolProp::PropsSI("P", "T", 300., "Q", 0, "PCSAFT::DIMETHYL ETHER");
    CHECK(abs((vp_calc / vp) - 1) < 1e-3);

    vp = 1.7551e-4;
    vp_calc = CoolProp::PropsSI("P", "T", 85.525, "Q", 0, "PCSAFT::PROPANE");
    CHECK(abs((vp_calc / vp) - 1) < 0.1);

    vp = 8.3324e5;
    vp_calc = CoolProp::PropsSI("P", "T", 293, "Q", 0, "PCSAFT::PROPANE");
    CHECK(abs((vp_calc / vp) - 1) < 0.01);

    vp = 42.477e5;
    vp_calc = CoolProp::PropsSI("P", "T", 369.82, "Q", 0, "PCSAFT::PROPANE");
    CHECK(abs((vp_calc / vp) - 1) < 0.01);
}

TEST_CASE("Check PC-SAFT interaction parameter functions", "[pcsaft_binary_interaction]") {
    std::string CAS_water = get_fluid_param_string("WATER", "CAS");
    std::string CAS_aacid = "64-19-7";
    set_mixture_binary_pair_pcsaft(CAS_water, CAS_aacid, "kij", -0.127);
    CHECK(atof(get_mixture_binary_pair_pcsaft(CAS_water, CAS_aacid, "kij").c_str()) == -0.127);
}

TEST_CASE("Check bubble pressures calculated using PC-SAFT", "[pcsaft_bubble_pressure]") {
    double vp = 1816840.45112607;
    double vp_calc = CoolProp::PropsSI("P", "T", 421.05, "Q", 0, "PCSAFT::METHANE[0.0252]&BENZENE[0.9748]");
    CHECK(abs((vp_calc / vp) - 1) < 1e-3);

    vp = 96634.2439079;
    vp_calc = CoolProp::PropsSI("P", "T", 327.48, "Q", 0, "PCSAFT::METHANOL[0.3]&CYCLOHEXANE[0.7]");
    CHECK(abs((vp_calc / vp) - 1) < 1e-3);

    // set binary interaction parameter
    std::string CAS_water = get_fluid_param_string("WATER", "CAS");
    std::string CAS_aacid = "64-19-7";
    set_mixture_binary_pair_pcsaft(CAS_water, CAS_aacid, "kij", -0.127);

    vp = 274890.39985918;
    vp_calc = CoolProp::PropsSI("P", "T", 403.574, "Q", 0, "PCSAFT::WATER[0.9898662364]&ACETIC ACID[0.0101337636]");
    CHECK(abs((vp_calc / vp) - 1) < 1e-2);

    vp = 72915.92217342;
    vp_calc = CoolProp::PropsSI("P", "T", 372.774, "Q", 0, "PCSAFT::WATER[0.2691800943]&ACETIC ACID[0.7308199057]");
    CHECK(abs((vp_calc / vp) - 1) < 2e-2);

    vp = 2387.42669687;
    vp_calc = CoolProp::PropsSI("P", "T", 298.15, "Q", 0, "PCSAFT::Na+[0.0907304774758426]&Cl-[0.0907304774758426]&WATER[0.818539045048315]");
    CHECK(abs((vp_calc / vp) - 1) < 1e-3);
}

TEST_CASE("Check bubble temperatures calculated using PC-SAFT", "[pcsaft_bubble_temperature]") {
    double t = 572.6667;
    double t_calc = CoolProp::PropsSI("T", "P", 3290651.18080112, "Q", 0, "PCSAFT::TOLUENE");
    CHECK(abs((t_calc / t) - 1) < 1e-3);

    t = 362;
    t_calc = CoolProp::PropsSI("T", "P", 66917.67387203, "Q", 0, "PCSAFT::WATER");
    CHECK(abs((t_calc / t) - 1) < 1e-3);

    t = 413.5385;
    t_calc = CoolProp::PropsSI("T", "P", 190061.78088909, "Q", 0, "PCSAFT::ACETIC ACID");
    CHECK(abs((t_calc / t) - 1) < 1e-3);

    t = 300.;
    t_calc = CoolProp::PropsSI("T", "P", 623027.07850612, "Q", 0, "PCSAFT::DIMETHYL ETHER");
    CHECK(abs((t_calc / t) - 1) < 1e-3);

    t = 421.05;
    t_calc = CoolProp::PropsSI("T", "P", 1816840.45112607, "Q", 0, "PCSAFT::METHANE[0.0252]&BENZENE[0.9748]");
    CHECK(abs((t_calc / t) - 1) < 1e-3);

    t = 327.48;
    t_calc = CoolProp::PropsSI("T", "P", 96634.2439079, "Q", 0, "PCSAFT::METHANOL[0.3]&CYCLOHEXANE[0.7]");
    CHECK(abs((t_calc / t) - 1) < 1e-3);

    // set binary interaction parameter
    std::string CAS_water = get_fluid_param_string("WATER", "CAS");
    std::string CAS_aacid = "64-19-7";
    set_mixture_binary_pair_pcsaft(CAS_water, CAS_aacid, "kij", -0.127);

    t = 403.574;
    t_calc = CoolProp::PropsSI("T", "P", 274890.39985918, "Q", 0, "PCSAFT::WATER[0.9898662364]&ACETIC ACID[0.0101337636]");
    CHECK(abs((t_calc / t) - 1) < 1e-3);

    t = 372.774;
    t_calc = CoolProp::PropsSI("T", "P", 72915.92217342, "Q", 0, "PCSAFT::WATER[0.2691800943]&ACETIC ACID[0.7308199057]");
    CHECK(abs((t_calc / t) - 1) < 2e-3);

    t = 298.15;
    t_calc = CoolProp::PropsSI("T", "P", 2387.42669687, "Q", 0, "PCSAFT::Na+[0.0907304774758426]&Cl-[0.0907304774758426]&WATER[0.818539045048315]");
    CHECK(abs((t_calc / t) - 1) < 1e-2);
}

TEST_CASE("Check phase determination for PC-SAFT backend", "[pcsaft_phase]") {
    double den = 9033.114209728405;
    double den_calc = CoolProp::PropsSI("Dmolar", "T", 320., "P", 101325., "PCSAFT::TOLUENE");
    CHECK(abs((den_calc / den) - 1) < 1e-2);
    double phase = CoolProp::PropsSI("Phase", "T", 320., "P", 101325., "PCSAFT::TOLUENE");
    CHECK(phase == get_phase_index("phase_liquid"));

    den = 0.376013;
    den_calc = CoolProp::PropsSI("Dmolar", "T", 320., "P", 1000., "PCSAFT::TOLUENE");
    CHECK(abs((den_calc / den) - 1) < 1e-2);
    phase = CoolProp::PropsSI("Phase", "T", 320., "P", 1000., "PCSAFT::TOLUENE");
    CHECK(phase == get_phase_index("phase_gas"));
}

/*
TEST_CASE("Test that HS solver works for a few fluids", "[HS_solver]")
{
    std::vector<std::string> fluids; fluids.push_back("Propane"); fluids.push_back("D4"); fluids.push_back("Water");
    for (std::size_t i = 0; i < fluids.size(); ++i)
    {
        std::vector<std::string> fl(1,fluids[i]);
        shared_ptr<CoolProp::HelmholtzEOSMixtureBackend> HEOS(new CoolProp::HelmholtzEOSMixtureBackend(fl));
        for (double p = HEOS->p_triple()*10; p < HEOS->pmax(); p *= 10)
        {
            double Tmin = HEOS->Ttriple();
            double Tmax = HEOS->Tmax();
            for (double T = Tmin + 1; T < Tmax-1; T += 10)
            {
                std::ostringstream ss;
                ss << "Check HS for " << fluids[i] << " for T=" << T << ", p=" << p;
                SECTION(ss.str(),"")
                {
                    CHECK_NOTHROW(HEOS->update(PT_INPUTS, p, T));
                    std::ostringstream ss1;
                    ss1 << "h=" << HEOS->hmolar() << ", s=" << HEOS->smolar();
                    SECTION(ss1.str(),"")
                    {
                        CAPTURE(T);
                        CAPTURE(p);
                        CAPTURE(HEOS->hmolar());
                        CAPTURE(HEOS->smolar());
                        CHECK_NOTHROW(HEOS->update(HmolarSmolar_INPUTS, HEOS->hmolar(), HEOS->smolar()));
                        double Terr = HEOS->T()- T;
                        CAPTURE(Terr);
                        CHECK(std::abs(Terr) < 1e-6);
                    }
                }
            }
        }
    }
}
*/
#endif
