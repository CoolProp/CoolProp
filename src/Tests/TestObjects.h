/**
 *  This file contains some basic methods to generate
 *  objects that can be used in the test routines.
 *  This makes the tests themselves much more readable
 *  and assures that the objects used for testing are the
 *  same in all places.
 */
#include "CoolProp/fluids/IncompressibleFluid.h"
#include "Eigen/Core"
#include "CoolProp/numerics/MatrixMath.h"
#include "CoolProp/Configuration.h"

#if defined ENABLE_CATCH
namespace CoolPropTesting {

Eigen::MatrixXd makeMatrix(const std::vector<double>& coefficients);
//CoolProp::IncompressibleFluid incompressibleFluidObject();
//IncompressibleBackend incompressibleBackendObject();

}  // namespace CoolPropTesting

/// Set INCOMP_PREFER_CHEBYSHEV_CALORIC for one scope and restore the previous
/// value on exit, also when an assertion throws, so no other test inherits it.
struct ChebyshevCaloricSetting
{
    bool saved;
    explicit ChebyshevCaloricSetting(bool value) : saved(CoolProp::get_config_bool(INCOMP_PREFER_CHEBYSHEV_CALORIC)) {
        CoolProp::set_config_bool(INCOMP_PREFER_CHEBYSHEV_CALORIC, value);
    }
    ~ChebyshevCaloricSetting() {
        CoolProp::set_config_bool(INCOMP_PREFER_CHEBYSHEV_CALORIC, saved);
    }
    ChebyshevCaloricSetting(const ChebyshevCaloricSetting&) = delete;
    ChebyshevCaloricSetting& operator=(const ChebyshevCaloricSetting&) = delete;
    ChebyshevCaloricSetting(ChebyshevCaloricSetting&&) = delete;
    ChebyshevCaloricSetting& operator=(ChebyshevCaloricSetting&&) = delete;
};
#endif  // ENABLE_CATCH
