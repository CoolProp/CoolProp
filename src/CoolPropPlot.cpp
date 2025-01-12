/*
 * AbstractState.cpp
 *
 *  Created on: 21 Dec 2013
 *      Author: jowr
 */

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "CoolPropPlot.h"

namespace CoolProp {
namespace Plot {


} /* namespace Plot */
} /* namespace CoolProp */

#ifdef ENABLE_CATCH
#    include <catch2/catch_all.hpp>
#    include <catch2/matchers/catch_matchers_floating_point.hpp>


TEST_CASE("Check CoolPropPlot", "[Plot]") {
    CHECK(1 == 1);
}
TEST_CASE("CoolPropPlot", "[log p-H plot value_at]") {
    CoolProp::Plot::PropertyPlot plot("R134a", CoolProp::iP, CoolProp::iHmass, "ACHP");
    plot.set_dimension_unit(plot.y_index, "bar");
    plot.set_dimension_unit(plot.x_index, "kJ/kg");
    plot.set_dimension_unit(CoolProp::iQ, "%");

    CHECK(*plot.value_at(CoolProp::iP, 300, 2) == 2);
    // EXPECT_NEAR(*plot.value_at(CoolProp::iT, 300, 2), 263.07372753976694, 1e-10);
    // EXPECT_NEAR(*plot.value_at(CoolProp::iQ, 300, 2), 55.044347874344737, 1e-10);
}

#endif
