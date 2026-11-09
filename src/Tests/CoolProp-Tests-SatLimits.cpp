// Saturation inputs below the minimum temperature of the equation of state (gh-3394).
//
// With superancillaries enabled (the default), QT_flash and PQ_flash for pure fluids
// returned early without the lower-limit check of the iterative path, so water accepted
// T = 250 K (Tmin = 273.16 K).  Run explicitly:  ./CatchTestRunner "[superanc_limits]"

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include <memory>

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/Configuration.h"
#    include "CoolProp/DataStructures.h"
#    include "CoolProp/Exceptions.h"

TEST_CASE("QT and PQ inputs below Tmin are rejected", "[superanc_limits]") {
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));
    REQUIRE(CoolProp::get_config_bool(ENABLE_SUPERANCILLARIES));
    REQUIRE(CoolProp::get_config_bool(DONT_CHECK_PROPERTY_LIMITS) == false);
    const double Tmin = AS->Tmin();

    SECTION("well below Tmin") {
        CHECK_THROWS_AS(AS->update(CoolProp::QT_INPUTS, 0.0, 250.0), CoolProp::ValueError);
        CHECK_THROWS_AS(AS->update(CoolProp::QT_INPUTS, 1.0, 250.0), CoolProp::ValueError);
        // saturation pressure of the superancillary at 250 K
        CHECK_THROWS_AS(AS->update(CoolProp::PQ_INPUTS, 95.2487, 0.0), CoolProp::ValueError);
    }
    SECTION("at and above Tmin still work") {
        CHECK_NOTHROW(AS->update(CoolProp::QT_INPUTS, 0.0, Tmin));
        const double p_at_Tmin = AS->p();
        CHECK_NOTHROW(AS->update(CoolProp::PQ_INPUTS, p_at_Tmin, 0.0));
        CHECK(AS->T() == Catch::Approx(Tmin).epsilon(1e-9));
        CHECK_NOTHROW(AS->update(CoolProp::QT_INPUTS, 0.5, 300.0));
        // the iterative path already allowed 0.1 K below Tmin; keep that
        CHECK_NOTHROW(AS->update(CoolProp::QT_INPUTS, 0.0, Tmin - 0.05));
    }
    SECTION("the limit check can still be disabled") {
        CoolProp::set_config_bool(DONT_CHECK_PROPERTY_LIMITS, true);
        CHECK_NOTHROW(AS->update(CoolProp::QT_INPUTS, 0.0, 250.0));
        CoolProp::set_config_bool(DONT_CHECK_PROPERTY_LIMITS, false);
    }
}

#endif  // ENABLE_CATCH
