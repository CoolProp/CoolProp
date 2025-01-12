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

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("Check value_at for log p-h plots", "[ph_plot]") {
    CoolProp::Plot::PropertyPlot plot("R134a", CoolProp::iP, CoolProp::iHmass, "ACHP");

    CHECK_THAT(*plot.value_at(CoolProp::iP, 300000/*Pa*/, 200000/*J/kg*/), WithinAbs(200000, 1e-10));
    CHECK_THAT(*plot.value_at(CoolProp::iHmass, 300000, 200000), WithinAbs(300000, 1e-10));
    CHECK_THAT(*plot.value_at(CoolProp::iT, 300000, 200000), WithinAbs(263.07372753976694, 1e-10));
    CHECK_THAT(*plot.value_at(CoolProp::iQ, 300000, 200000), WithinAbs(0.55044347874344737, 1e-10));
}

TEST_CASE("Check the output is the same as Python", "[ph_plot]") {
    CoolProp::Plot::PropertyPlot plot("HEOS::R134a", CoolProp::iP, CoolProp::iHmass, "ACHP");
    const int isoline_count = 5;
    const int isoline_points = 5;

    CHECK(plot.x_index == CoolProp::iHmass);
    CHECK(plot.y_index == CoolProp::iP);

    CHECK(plot.axis_x_scale() == CoolProp::Plot::Scale::Lin);
    CHECK(plot.axis_y_scale() == CoolProp::Plot::Scale::Log);

    CHECK_THAT(plot.axis_x_range().min, WithinAbs(75373.1, 1));
    CHECK_THAT(plot.axis_x_range().max, WithinAbs(577605, 1));
    CHECK_THAT(plot.axis_y_range().min, WithinAbs(25000, 1));
    CHECK_THAT(plot.axis_y_range().max, WithinAbs(9.133e6, 1));

    auto iso_types = plot.supported_dimensions();
    REQUIRE(iso_types.size() == 4);
    CHECK(iso_types[0] == CoolProp::iQ);
    CHECK(iso_types[1] == CoolProp::iT);
    CHECK(iso_types[2] == CoolProp::iSmass);
    CHECK(iso_types[3] == CoolProp::iDmass);

    {
        // Q isolines
        CoolProp::Plot::Range iso_range = plot.isoline_range(CoolProp::iQ);
        std::vector<double> iso_values = CoolProp::Plot::Detail::generate_ranges(CoolProp::iQ, iso_range.min, iso_range.max, isoline_count);
        auto isoQ = plot.calc_isolines(CoolProp::iQ, iso_values, isoline_points);
        REQUIRE(isoQ.size() == isoline_count);
        CHECK_THAT(isoQ[0].value, WithinAbs(0.0, 1e-10));
        CHECK_THAT(isoQ[1].value, WithinAbs(0.25, 1e-10));
        CHECK_THAT(isoQ[2].value, WithinAbs(0.5, 1e-10));
        CHECK_THAT(isoQ[3].value, WithinAbs(0.75, 1e-10));
        CHECK_THAT(isoQ[4].value, WithinAbs(1.0, 1e-10));
        const double expected_x[isoline_count][isoline_points] = {
            {71455.08256999, 132939.99472497, 198497.0525314, 271576.58908888, 389440.73899019},
            {137326.83116781, 191267.05172013, 248359.90642508, 309538.95484829, 389511.40516519},
            {203198.57976564, 249594.1087153, 298222.76031877, 347501.3206077, 389582.0713402},
            {269070.32836347, 307921.16571046, 348085.61421246, 385463.68636711, 389652.73751521},
            {334942.07696129, 366248.22270562, 397948.46810615, 423426.05212652, 389723.40369022}
        };
        const double expected_y[isoline_count][isoline_points] = {
            {3.89567060e+02, 2.58505756e+04, 2.81105747e+05, 1.31691170e+06, 4.05910826e+06},
            {3.89567060e+02, 2.58505756e+04, 2.81105747e+05, 1.31691170e+06, 4.05910826e+06},
            {3.89567060e+02, 2.58505756e+04, 2.81105747e+05, 1.31691170e+06, 4.05910826e+06},
            {3.89567060e+02, 2.58505756e+04, 2.81105747e+05, 1.31691170e+06, 4.05910826e+06},
            {3.89567060e+02, 2.58505756e+04, 2.81105747e+05, 1.31691170e+06, 4.05910826e+06}
        };
        for (int i = 0; i < isoline_count; ++i) {
            for (int j = 0; j < isoline_points; ++j) {
                CHECK_THAT(isoQ[i].x[j], WithinRel(expected_x[i][j], 1e-8));
                CHECK_THAT(isoQ[i].y[j], WithinRel(expected_y[i][j], 1e-8));
            }
        }
    }
    {
        // T isolines
        CoolProp::Plot::Range iso_range = plot.isoline_range(CoolProp::iT);
        std::vector<double> iso_values = CoolProp::Plot::Detail::generate_ranges(CoolProp::iT, iso_range.min, iso_range.max, isoline_count);
        auto isoT = plot.calc_isolines(CoolProp::iT, iso_values, isoline_points);
        REQUIRE(isoT.size() == isoline_count);
        CHECK_THAT(isoT[0].value, WithinAbs(173.15, 1e-10));
        CHECK_THAT(isoT[1].value, WithinAbs(243.6125, 1e-10));
        CHECK_THAT(isoT[2].value, WithinAbs(314.075, 1e-10));
        CHECK_THAT(isoT[3].value, WithinAbs(384.5375, 1e-10));
        CHECK_THAT(isoT[4].value, WithinAbs(455.0, 1e-10));
        const double expected_x[isoline_count][isoline_points] = {
            {75373.12689908, 75410.99061368, 75576.57734102, 76301.46320034, 79487.71936123},
            {382785.23058756, 161389.44197265, 161516.21527543, 162076.96181603, 164636.92352411},
            {439466.64984328, 438148.19040202, 431912.24266074, 257605.32897193, 257512.80690587},
            {504550.62606561, 503783.53950874, 500331.68636519, 482707.98489058, 366959.25257106},
            {577604.59497521, 577097.07174302, 574850.21206939, 564444.22079795, 507878.75380996}
        };
        const double expected_y[isoline_count][isoline_points] = {
            {25000.0, 109297.03270763, 477833.65434772, 2089032.02192197, 9133000.04909175},
            {25000.0, 109297.03270763, 477833.65434772, 2089032.02192197, 9133000.04909175},
            {25000.0, 109297.03270763, 477833.65434772, 2089032.02192197, 9133000.04909175},
            {25000.0, 109297.03270763, 477833.65434772, 2089032.02192197, 9133000.04909175},
            {25000.0, 109297.03270763, 477833.65434772, 2089032.02192197, 9133000.04909175}
        };
        for (int i = 0; i < isoline_count; ++i) {
            for (int j = 0; j < isoline_points; ++j) {
                CHECK_THAT(isoT[i].x[j], WithinRel(expected_x[i][j], 1e-8));
                CHECK_THAT(isoT[i].y[j], WithinRel(expected_y[i][j], 1e-8));
            }
        }
    }
    {
        // S isolines
        CoolProp::Plot::Range iso_range = plot.isoline_range(CoolProp::iSmass);
        std::vector<double> iso_values = CoolProp::Plot::Detail::generate_ranges(CoolProp::iSmass, iso_range.min, iso_range.max, isoline_count);
        auto isoS = plot.calc_isolines(CoolProp::iSmass, iso_values, isoline_points);
        REQUIRE(isoS.size() == isoline_count);
        CHECK_THAT(isoS[0].value, WithinAbs(426.0098553415565, 1e-10));
        CHECK_THAT(isoS[1].value, WithinAbs(925.2753143522199, 1e-10));
        CHECK_THAT(isoS[2].value, WithinAbs(1424.540773362883, 1e-10));
        CHECK_THAT(isoS[3].value, WithinAbs(1923.8062323735467, 1e-10));
        CHECK_THAT(isoS[4].value, WithinAbs(2423.07169138421, 1e-10));
        const double expected_x[isoline_count][isoline_points] = {
            {73758.20064883, 73811.34928194, 74043.68191583, 75058.90103546, 79487.71936135},
            {176257.41124603, 179794.86552304, 180290.38376303, 181487.99233203, 186690.75978367},
            {286286.21675221, 303984.6485323, 321692.18498473, 335551.56116185, 344087.52152244},
            {399372.58517377, 433400.13264058, 471964.37092222, 513835.0906295, 555823.55477128},
            {577604.59497522, 635257.81956317, 698998.52697158, 768744.13132575, std::nan("")}
        };
        const double expected_y[isoline_count][isoline_points] = {
            {25000.0, 109297.03270763, 477833.65434772, 2089032.02192197, 9133000.04909175},
            {25000.0, 109297.03270763, 477833.65434772, 2089032.02192197, 9133000.04909175},
            {25000.0, 109297.03270763, 477833.65434772, 2089032.02192197, 9133000.04909175},
            {25000.0, 109297.03270763, 477833.65434772, 2089032.02192197, 9133000.04909175},
            {25000.0, 109297.03270763, 477833.65434772, 2089032.02192197, 9133000.04909175}
        };
        for (int i = 0; i < isoline_count; ++i) {
            for (int j = 0; j < isoline_points; ++j) {
                if (std::isnan(isoS[i].x[j]))
                    CHECK(std::isnan(expected_x[i][j]));
                else
                    CHECK_THAT(isoS[i].x[j], WithinRel(expected_x[i][j], 1e-8));
                CHECK_THAT(isoS[i].y[j], WithinRel(expected_y[i][j], 1e-8));
            }
        }
    }
    {
        // D isolines
        CoolProp::Plot::Range iso_range = plot.isoline_range(CoolProp::iDmass);
        std::vector<double> iso_values = CoolProp::Plot::Detail::generate_ranges(CoolProp::iDmass, iso_range.min, iso_range.max, isoline_count);
        auto isoD = plot.calc_isolines(CoolProp::iDmass, iso_values, isoline_points);
        REQUIRE(isoD.size() == isoline_count);
        CHECK_THAT(isoD[0].value, WithinAbs(0.6749779869915704, 1e-10));
        CHECK_THAT(isoD[1].value, WithinAbs(4.704765330619733, 1e-10));
        CHECK_THAT(isoD[2].value, WithinAbs(32.793390662794806, 1e-10));
        CHECK_THAT(isoD[3].value, WithinAbs(228.57813208316188, 1e-10));
        CHECK_THAT(isoD[4].value, WithinAbs(1593.2467308391022, 1e-10));
        const double expected_x[isoline_count][isoline_points] = {
            {5.77604595e+05, 3.65397965e+06, 3.84283606e+07, 4.40954815e+08, 5.22494051e+09},
            {2.02365849e+05, 4.19227802e+05, 1.84512838e+06, 1.78590373e+07, 2.01674048e+08},
            {1.42114492e+05, 2.04387395e+05, 3.51213643e+05, 1.00846111e+06, 8.12669299e+06},
            {1.33470419e+05, 1.72415441e+05, 2.35381904e+05, 3.57488786e+05, 6.69475618e+05},
            {7.05185202e+04, 7.06013993e+04, 7.09637630e+04, 7.25484894e+04, 7.94877194e+04}
        };
        const double expected_y[isoline_count][isoline_points] = {
            {25000.0, 109297.03270763, 477833.65434772, 2089032.02192197, 9133000.04909175},
            {25000.0, 109297.03270763, 477833.65434772, 2089032.02192197, 9133000.04909175},
            {25000.0, 109297.03270763, 477833.65434772, 2089032.02192197, 9133000.04909175},
            {25000.0, 109297.03270763, 477833.65434772, 2089032.02192197, 9133000.04909175},
            {25000.0, 109297.03270763, 477833.65434772, 2089032.02192197, 9133000.04909175}
        };
        for (int i = 0; i < isoline_count; ++i) {
            for (int j = 0; j < isoline_points; ++j) {
                CHECK_THAT(isoD[i].x[j], WithinRel(expected_x[i][j], 1e-8));
                CHECK_THAT(isoD[i].y[j], WithinRel(expected_y[i][j], 1e-8));
            }
        }
    }
}

#endif
