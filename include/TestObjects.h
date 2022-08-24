/**
 *  This file contains some basic methods to generate
 *  objects that can be used in the test routines.
 *  This makes the tests themselves much more readable
 *  and assures that the objects used for testing are the
 *  same in all places.
 */
#include "IncompressibleFluid.h"
#include "Eigen/Core"
#include "MatrixMath.h"

#if defined ENABLE_CATCH
namespace CoolPropTesting {

Eigen::MatrixXd makeMatrix(const std::vector<double>& coefficients);
//CoolProp::IncompressibleFluid incompressibleFluidObject();
//IncompressibleBackend incompressibleBackendObject();

}  // namespace CoolPropTesting
#endif  // ENABLE_CATCH
