/**
 *  This file contains some basic methods to generate
 *  objects that can be used in the test routines.
 *  This makes the tests themselves much more readable
 *  and assures that the objects used for testing are the
 *  same in all places.
 */
#include "TestObjects.h"
#include "DataStructures.h"
#include "IncompressibleFluid.h"
#include "Eigen/Core"

#if defined ENABLE_CATCH

Eigen::MatrixXd CoolPropTesting::makeMatrix(const std::vector<double>& coefficients) {
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
    tmpVector.push_back(coefficients[1] * 100.0);
    tmpVector.push_back(coefficients[7] * 100.0);
    tmpVector.push_back(coefficients[12] * 100.0);
    tmpVector.push_back(coefficients[16] * 100.0);
    matrix.push_back(tmpVector);

    tmpVector.clear();
    tmpVector.push_back(coefficients[2] * 100.0 * 100.0);
    tmpVector.push_back(coefficients[8] * 100.0 * 100.0);
    tmpVector.push_back(coefficients[13] * 100.0 * 100.0);
    tmpVector.push_back(coefficients[17] * 100.0 * 100.0);
    matrix.push_back(tmpVector);

    tmpVector.clear();
    tmpVector.push_back(coefficients[3] * 100.0 * 100.0 * 100.0);
    tmpVector.push_back(coefficients[9] * 100.0 * 100.0 * 100.0);
    tmpVector.push_back(coefficients[14] * 100.0 * 100.0 * 100.0);
    tmpVector.push_back(0.0);
    matrix.push_back(tmpVector);

    tmpVector.clear();
    tmpVector.push_back(coefficients[4] * 100.0 * 100.0 * 100.0 * 100.0);
    tmpVector.push_back(coefficients[10] * 100.0 * 100.0 * 100.0 * 100.0);
    tmpVector.push_back(0.0);
    tmpVector.push_back(0.0);
    matrix.push_back(tmpVector);

    tmpVector.clear();
    tmpVector.push_back(coefficients[5] * 100.0 * 100.0 * 100.0 * 100.0 * 100.0);
    tmpVector.push_back(0.0);
    tmpVector.push_back(0.0);
    tmpVector.push_back(0.0);
    matrix.push_back(tmpVector);

    tmpVector.clear();
    return CoolProp::vec_to_eigen(matrix).transpose();
}

CoolProp::IncompressibleFluid CoolPropTesting::incompressibleFluidObject() {

    std::string tmpStr;
    std::vector<double> tmpVector;
    std::vector<std::vector<double>> tmpMatrix;

    tmpVector.clear();
    tmpVector.push_back(960.24665800);
    tmpVector.push_back(-1.2903839100);
    tmpVector.push_back(-0.0161042520);
    tmpVector.push_back(-0.0001969888);
    tmpVector.push_back(1.131559E-05);
    tmpVector.push_back(9.181999E-08);
    tmpVector.push_back(-0.4020348270);
    tmpVector.push_back(-0.0162463989);
    tmpVector.push_back(0.0001623301);
    tmpVector.push_back(4.367343E-06);
    tmpVector.push_back(1.199000E-08);
    tmpVector.push_back(-0.0025204776);
    tmpVector.push_back(0.0001101514);
    tmpVector.push_back(-2.320217E-07);
    tmpVector.push_back(7.794999E-08);
    tmpVector.push_back(9.937483E-06);
    tmpVector.push_back(-1.346886E-06);
    tmpVector.push_back(4.141999E-08);
    CoolProp::IncompressibleData density;
    density.type = CoolProp::IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL;
    density.coeffs = makeMatrix(tmpVector);

    tmpVector.clear();
    tmpVector.push_back(3822.9712300);
    tmpVector.push_back(-23.122409500);
    tmpVector.push_back(0.0678775826);
    tmpVector.push_back(0.0022413893);
    tmpVector.push_back(-0.0003045332);
    tmpVector.push_back(-4.758000E-06);
    tmpVector.push_back(2.3501449500);
    tmpVector.push_back(0.1788839410);
    tmpVector.push_back(0.0006828000);
    tmpVector.push_back(0.0002101166);
    tmpVector.push_back(-9.812000E-06);
    tmpVector.push_back(-0.0004724176);
    tmpVector.push_back(-0.0003317949);
    tmpVector.push_back(0.0001002032);
    tmpVector.push_back(-5.306000E-06);
    tmpVector.push_back(4.242194E-05);
    tmpVector.push_back(2.347190E-05);
    tmpVector.push_back(-1.894000E-06);
    CoolProp::IncompressibleData specific_heat;
    specific_heat.type = CoolProp::IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL;
    specific_heat.coeffs = makeMatrix(tmpVector);

    tmpVector.clear();
    tmpVector.push_back(0.4082066700);
    tmpVector.push_back(-0.0039816870);
    tmpVector.push_back(1.583368E-05);
    tmpVector.push_back(-3.552049E-07);
    tmpVector.push_back(-9.884176E-10);
    tmpVector.push_back(4.460000E-10);
    tmpVector.push_back(0.0006629321);
    tmpVector.push_back(-2.686475E-05);
    tmpVector.push_back(9.039150E-07);
    tmpVector.push_back(-2.128257E-08);
    tmpVector.push_back(-5.562000E-10);
    tmpVector.push_back(3.685975E-07);
    tmpVector.push_back(7.188416E-08);
    tmpVector.push_back(-1.041773E-08);
    tmpVector.push_back(2.278001E-10);
    tmpVector.push_back(4.703395E-08);
    tmpVector.push_back(7.612361E-11);
    tmpVector.push_back(-2.734000E-10);
    CoolProp::IncompressibleData conductivity;
    conductivity.type = CoolProp::IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL;
    conductivity.coeffs = makeMatrix(tmpVector);

    tmpVector.clear();
    tmpVector.push_back(1.4725525500);
    tmpVector.push_back(0.0022218998);
    tmpVector.push_back(-0.0004406139);
    tmpVector.push_back(6.047984E-06);
    tmpVector.push_back(-1.954730E-07);
    tmpVector.push_back(-2.372000E-09);
    tmpVector.push_back(-0.0411841566);
    tmpVector.push_back(0.0001784479);
    tmpVector.push_back(-3.564413E-06);
    tmpVector.push_back(4.064671E-08);
    tmpVector.push_back(1.915000E-08);
    tmpVector.push_back(0.0002572862);
    tmpVector.push_back(-9.226343E-07);
    tmpVector.push_back(-2.178577E-08);
    tmpVector.push_back(-9.529999E-10);
    tmpVector.push_back(-1.699844E-06);
    tmpVector.push_back(-1.023552E-07);
    tmpVector.push_back(4.482000E-09);
    CoolProp::IncompressibleData viscosity;
    viscosity.type = CoolProp::IncompressibleData::INCOMPRESSIBLE_EXPPOLYNOMIAL;
    viscosity.coeffs = makeMatrix(tmpVector);

    tmpVector.clear();
    tmpVector.push_back(27.755555600 / 100.0);  // reference concentration in per cent
    tmpVector.push_back(-22.973221700 + 273.15);
    tmpVector.push_back(-1.1040507200 * 100.0);
    tmpVector.push_back(-0.0120762281 * 100.0 * 100.0);
    tmpVector.push_back(-9.343458E-05 * 100.0 * 100.0 * 100.0);
    CoolProp::IncompressibleData T_freeze;
    T_freeze.type = CoolProp::IncompressibleData::INCOMPRESSIBLE_POLYOFFSET;
    T_freeze.coeffs = CoolProp::vec_to_eigen(tmpVector);
    T_freeze.coeffs.transposeInPlace();

    // After preparing the coefficients, we have to create the objects
    CoolProp::IncompressibleFluid CH3OH;
    CH3OH.setName("CH3OH-testing");
    CH3OH.setDescription("Methanol solution");
    CH3OH.setReference("SecCool software");
    CH3OH.setTmax(20 + 273.15);
    CH3OH.setTmin(-50 + 273.15);
    CH3OH.setxmax(0.5);
    CH3OH.setxmin(0.0);
    CH3OH.setxid(CoolProp::IFRAC_MASS);
    CH3OH.setTminPsat(20 + 273.15);

    CH3OH.setTbase(-4.48 + 273.15);
    CH3OH.setxbase(31.57 / 100.0);

    /// Setters for the coefficients
    CH3OH.setDensity(density);
    CH3OH.setSpecificHeat(specific_heat);
    CH3OH.setViscosity(viscosity);
    CH3OH.setConductivity(conductivity);
    //CH3OH.setPsat(saturation_pressure);
    CH3OH.setTfreeze(T_freeze);
    //CH3OH.setVolToMass(volume2mass);
    //CH3OH.setMassToMole(mass2mole);

    /// A function to check coefficients and equation types.
    CH3OH.validate();

    return CH3OH;
}
//CoolProp::IncompressibleBackend CoolProp::Testing::incompressibleBackendObject(){
//    return CoolProp::IncompressibleBackend(CoolProp::Testing::incompressibleFluidObject());
//}

#endif  // ENABLE_CATCH
