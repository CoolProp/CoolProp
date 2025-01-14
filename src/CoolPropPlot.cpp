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

TEST_CASE("Check value_at for p-h plots", "[ph_plot]") {
    CoolProp::Plot::PropertyPlot plot("R134a", CoolProp::iP, CoolProp::iHmass, "ACHP");

    CHECK_THAT(*plot.value_at(CoolProp::iP, 300000/*Pa*/, 200000/*J/kg*/), WithinAbs(200000, 1e-10));
    CHECK_THAT(*plot.value_at(CoolProp::iHmass, 300000, 200000), WithinAbs(300000, 1e-10));
    CHECK_THAT(*plot.value_at(CoolProp::iT, 300000, 200000), WithinAbs(263.07372753976694, 1e-10));
    CHECK_THAT(*plot.value_at(CoolProp::iQ, 300000, 200000), WithinAbs(0.55044347874344737, 1e-10));
}

TEST_CASE("Check that the isolines are the same as from Python", "[ph_plot]") {
    CoolProp::Plot::PropertyPlot plot("HEOS::R134a", CoolProp::iP, CoolProp::iHmass, "ACHP");
    const int isoline_count = 5;
    const int points_per_isoline = 5;

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
        auto isoQ = plot.calc_isolines(CoolProp::iQ, iso_values, points_per_isoline);
        REQUIRE(isoQ.size() == isoline_count);
        CHECK_THAT(isoQ[0].value, WithinAbs(0.0, 1e-10));
        CHECK_THAT(isoQ[1].value, WithinAbs(0.25, 1e-10));
        CHECK_THAT(isoQ[2].value, WithinAbs(0.5, 1e-10));
        CHECK_THAT(isoQ[3].value, WithinAbs(0.75, 1e-10));
        CHECK_THAT(isoQ[4].value, WithinAbs(1.0, 1e-10));
        const double expected_x[isoline_count][points_per_isoline] = {
            {71455.08256999, 132939.99472497, 198497.0525314, 271576.58908888, 389440.73899019},
            {137326.83116781, 191267.05172013, 248359.90642508, 309538.95484829, 389511.40516519},
            {203198.57976564, 249594.1087153, 298222.76031877, 347501.3206077, 389582.0713402},
            {269070.32836347, 307921.16571046, 348085.61421246, 385463.68636711, 389652.73751521},
            {334942.07696129, 366248.22270562, 397948.46810615, 423426.05212652, 389723.40369022}
        };
        const double expected_y[isoline_count][points_per_isoline] = {
            {3.89567060e+02, 2.58505756e+04, 2.81105747e+05, 1.31691170e+06, 4.05910826e+06},
            {3.89567060e+02, 2.58505756e+04, 2.81105747e+05, 1.31691170e+06, 4.05910826e+06},
            {3.89567060e+02, 2.58505756e+04, 2.81105747e+05, 1.31691170e+06, 4.05910826e+06},
            {3.89567060e+02, 2.58505756e+04, 2.81105747e+05, 1.31691170e+06, 4.05910826e+06},
            {3.89567060e+02, 2.58505756e+04, 2.81105747e+05, 1.31691170e+06, 4.05910826e+06}
        };
        for (int i = 0; i < isoline_count; ++i) {
            for (int j = 0; j < points_per_isoline; ++j) {
                CHECK_THAT(isoQ[i].x[j], WithinRel(expected_x[i][j], 1e-8));
                CHECK_THAT(isoQ[i].y[j], WithinRel(expected_y[i][j], 1e-8));
            }
        }
    }
    {
        // T isolines
        CoolProp::Plot::Range iso_range = plot.isoline_range(CoolProp::iT);
        std::vector<double> iso_values = CoolProp::Plot::Detail::generate_ranges(CoolProp::iT, iso_range.min, iso_range.max, isoline_count);
        auto isoT = plot.calc_isolines(CoolProp::iT, iso_values, points_per_isoline);
        REQUIRE(isoT.size() == isoline_count);
        CHECK_THAT(isoT[0].value, WithinAbs(173.15, 1e-10));
        CHECK_THAT(isoT[1].value, WithinAbs(243.6125, 1e-10));
        CHECK_THAT(isoT[2].value, WithinAbs(314.075, 1e-10));
        CHECK_THAT(isoT[3].value, WithinAbs(384.5375, 1e-10));
        CHECK_THAT(isoT[4].value, WithinAbs(455.0, 1e-10));
        const double expected_x[isoline_count][points_per_isoline] = {
            {75373.12689908, 75410.99061368, 75576.57734102, 76301.46320034, 79487.71936123},
            {382785.23058756, 161389.44197265, 161516.21527543, 162076.96181603, 164636.92352411},
            {439466.64984328, 438148.19040202, 431912.24266074, 257605.32897193, 257512.80690587},
            {504550.62606561, 503783.53950874, 500331.68636519, 482707.98489058, 366959.25257106},
            {577604.59497521, 577097.07174302, 574850.21206939, 564444.22079795, 507878.75380996}
        };
        const double expected_y[isoline_count][points_per_isoline] = {
            {25000.0, 109297.03270763, 477833.65434772, 2089032.02192197, 9133000.04909175},
            {25000.0, 109297.03270763, 477833.65434772, 2089032.02192197, 9133000.04909175},
            {25000.0, 109297.03270763, 477833.65434772, 2089032.02192197, 9133000.04909175},
            {25000.0, 109297.03270763, 477833.65434772, 2089032.02192197, 9133000.04909175},
            {25000.0, 109297.03270763, 477833.65434772, 2089032.02192197, 9133000.04909175}
        };
        for (int i = 0; i < isoline_count; ++i) {
            for (int j = 0; j < points_per_isoline; ++j) {
                CHECK_THAT(isoT[i].x[j], WithinRel(expected_x[i][j], 1e-8));
                CHECK_THAT(isoT[i].y[j], WithinRel(expected_y[i][j], 1e-8));
            }
        }
    }
    {
        // S isolines
        CoolProp::Plot::Range iso_range = plot.isoline_range(CoolProp::iSmass);
        std::vector<double> iso_values = CoolProp::Plot::Detail::generate_ranges(CoolProp::iSmass, iso_range.min, iso_range.max, isoline_count);
        auto isoS = plot.calc_isolines(CoolProp::iSmass, iso_values, points_per_isoline);
        REQUIRE(isoS.size() == isoline_count);
        CHECK_THAT(isoS[0].value, WithinAbs(426.0098553415565, 1e-10));
        CHECK_THAT(isoS[1].value, WithinAbs(925.2753143522199, 1e-10));
        CHECK_THAT(isoS[2].value, WithinAbs(1424.540773362883, 1e-10));
        CHECK_THAT(isoS[3].value, WithinAbs(1923.8062323735467, 1e-10));
        CHECK_THAT(isoS[4].value, WithinAbs(2423.07169138421, 1e-10));
        const double expected_x[isoline_count][points_per_isoline] = {
            {73758.20064883, 73811.34928194, 74043.68191583, 75058.90103546, 79487.71936135},
            {176257.41124603, 179794.86552304, 180290.38376303, 181487.99233203, 186690.75978367},
            {286286.21675221, 303984.6485323, 321692.18498473, 335551.56116185, 344087.52152244},
            {399372.58517377, 433400.13264058, 471964.37092222, 513835.0906295, 555823.55477128},
            {577604.59497522, 635257.81956317, 698998.52697158, 768744.13132575, std::nan("")}
        };
        const double expected_y[isoline_count][points_per_isoline] = {
            {25000.0, 109297.03270763, 477833.65434772, 2089032.02192197, 9133000.04909175},
            {25000.0, 109297.03270763, 477833.65434772, 2089032.02192197, 9133000.04909175},
            {25000.0, 109297.03270763, 477833.65434772, 2089032.02192197, 9133000.04909175},
            {25000.0, 109297.03270763, 477833.65434772, 2089032.02192197, 9133000.04909175},
            {25000.0, 109297.03270763, 477833.65434772, 2089032.02192197, 9133000.04909175}
        };
        for (int i = 0; i < isoline_count; ++i) {
            for (int j = 0; j < points_per_isoline; ++j) {
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
        auto isoD = plot.calc_isolines(CoolProp::iDmass, iso_values, points_per_isoline);
        REQUIRE(isoD.size() == isoline_count);
        CHECK_THAT(isoD[0].value, WithinAbs(0.6749779869915704, 1e-10));
        CHECK_THAT(isoD[1].value, WithinAbs(4.704765330619733, 1e-10));
        CHECK_THAT(isoD[2].value, WithinAbs(32.793390662794806, 1e-10));
        CHECK_THAT(isoD[3].value, WithinAbs(228.57813208316188, 1e-10));
        CHECK_THAT(isoD[4].value, WithinAbs(1593.2467308391022, 1e-10));
        const double expected_x[isoline_count][points_per_isoline] = {
            {5.77604595e+05, 3.65397965e+06, 3.84283606e+07, 4.40954815e+08, 5.22494051e+09},
            {2.02365849e+05, 4.19227802e+05, 1.84512838e+06, 1.78590373e+07, 2.01674048e+08},
            {1.42114492e+05, 2.04387395e+05, 3.51213643e+05, 1.00846111e+06, 8.12669299e+06},
            {1.33470419e+05, 1.72415441e+05, 2.35381904e+05, 3.57488786e+05, 6.69475618e+05},
            {7.05185202e+04, 7.06013993e+04, 7.09637630e+04, 7.25484894e+04, 7.94877194e+04}
        };
        const double expected_y[isoline_count][points_per_isoline] = {
            {25000.0, 109297.03270763, 477833.65434772, 2089032.02192197, 9133000.04909175},
            {25000.0, 109297.03270763, 477833.65434772, 2089032.02192197, 9133000.04909175},
            {25000.0, 109297.03270763, 477833.65434772, 2089032.02192197, 9133000.04909175},
            {25000.0, 109297.03270763, 477833.65434772, 2089032.02192197, 9133000.04909175},
            {25000.0, 109297.03270763, 477833.65434772, 2089032.02192197, 9133000.04909175}
        };
        for (int i = 0; i < isoline_count; ++i) {
            for (int j = 0; j < points_per_isoline; ++j) {
                CHECK_THAT(isoD[i].x[j], WithinRel(expected_x[i][j], 1e-8));
                CHECK_THAT(isoD[i].y[j], WithinRel(expected_y[i][j], 1e-8));
            }
        }
    }
}

TEST_CASE("Basic TS Plot has same output as Python", "[ts_plot]") {
    CoolProp::Plot::PropertyPlot plot("HEOS::R134a", CoolProp::iT, CoolProp::iSmass, "ACHP");
    const int isoline_count = 5;
    const int points_per_isoline = 5;

    CHECK(plot.x_index == CoolProp::iSmass);
    CHECK(plot.y_index == CoolProp::iT);

    CHECK(plot.axis_x_scale() == CoolProp::Plot::Scale::Lin);
    CHECK(plot.axis_y_scale() == CoolProp::Plot::Scale::Lin);

    CHECK_THAT(plot.axis_x_range().min, WithinAbs(426, 1));
    CHECK_THAT(plot.axis_x_range().max, WithinAbs(2423, 1));
    CHECK_THAT(plot.axis_y_range().min, WithinAbs(173, 1));
    CHECK_THAT(plot.axis_y_range().max, WithinAbs(455, 1));

    auto iso_types = plot.supported_dimensions();
    REQUIRE(iso_types.size() == 4);
    CHECK(iso_types[0] == CoolProp::iQ);
    CHECK(iso_types[1] == CoolProp::iP);
    CHECK(iso_types[2] == CoolProp::iHmass);
    CHECK(iso_types[3] == CoolProp::iDmass);

    {
        // Q isolines
        CoolProp::Plot::Range iso_range = plot.isoline_range(CoolProp::iQ);
        std::vector<double> iso_values = CoolProp::Plot::Detail::generate_ranges(CoolProp::iQ, iso_range.min, iso_range.max, isoline_count);
        auto isoQ = plot.calc_isolines(CoolProp::iQ, iso_values, points_per_isoline);
        REQUIRE(isoQ.size() == isoline_count);
        CHECK_THAT(isoQ[0].value, WithinAbs(0.0, 1e-10));
        CHECK_THAT(isoQ[1].value, WithinAbs(0.25, 1e-10));
        CHECK_THAT(isoQ[2].value, WithinAbs(0.5, 1e-10));
        CHECK_THAT(isoQ[3].value, WithinAbs(0.75, 1e-10));
        CHECK_THAT(isoQ[4].value, WithinAbs(1.0, 1e-10));

        const double expected_x[isoline_count][points_per_isoline] = {
            {412.61753823, 728.71208313, 994.51958846, 1237.3122956, 1561.56967158},
            {800.44043827, 992.70703491, 1177.81867482, 1354.79919476, 1561.75851162},
            {1188.26333832, 1256.70198669, 1361.11776119, 1472.28609393, 1561.94735166},
            {1576.08623836, 1520.69693847, 1544.41684755, 1589.77299309, 1562.1361917},
            {1963.9091384, 1784.69189025, 1727.71593392, 1707.25989226, 1562.32503174}
        };
        const double expected_y[isoline_count][points_per_isoline] = {
            {169.85007484, 220.94004678, 272.03001871, 323.11999064, 374.20996258},
            {169.85007484, 220.94004678, 272.03001871, 323.11999064, 374.20996258},
            {169.85007484, 220.94004678, 272.03001871, 323.11999064, 374.20996258},
            {169.85007484, 220.94004678, 272.03001871, 323.11999064, 374.20996258},
            {169.85007484, 220.94004678, 272.03001871, 323.11999064, 374.20996258}
        };

        for(int i = 0; i < isoline_count; i++) {
            for(int j = 0; j < points_per_isoline; j++) {
                CHECK_THAT(isoQ[i].x[j], WithinRel(expected_x[i][j], 1e-8));
                CHECK_THAT(isoQ[i].y[j], WithinRel(expected_y[i][j], 1e-8));
            }
        }
    }

    {
        // P isolines
        CoolProp::Plot::Range iso_range = plot.isoline_range(CoolProp::iP);
        std::vector<double> iso_values = CoolProp::Plot::Detail::generate_ranges(CoolProp::iP, iso_range.min, iso_range.max, isoline_count);
        auto isoP = plot.calc_isolines(CoolProp::iP, iso_values, points_per_isoline);
        REQUIRE(isoP.size() == isoline_count);
        CHECK_THAT(isoP[0].value, WithinAbs(25000.000000000007, 1e-8));
        CHECK_THAT(isoP[1].value, WithinAbs(109297.03270763098, 1e-8));
        CHECK_THAT(isoP[2].value, WithinAbs(477833.65434771683, 1e-8));
        CHECK_THAT(isoP[3].value, WithinAbs(2089032.0219219688, 1e-8));
        CHECK_THAT(isoP[4].value, WithinAbs(9133000.049091753, 1e-8));

        const double expected_x[isoline_count][points_per_isoline] = {
            {426.00985534, 925.27531435, 1424.54077336, 1923.80623237, 2423.07169138},
            {426.00985534, 925.27531435, 1424.54077336, 1923.80623237, 2423.07169138},
            {426.00985534, 925.27531435, 1424.54077336, 1923.80623237, 2423.07169138},
            {426.00985534, 925.27531435, 1424.54077336, 1923.80623237, 2423.07169138},
            {426.00985534, 925.27531435, 1424.54077336, 1923.80623237, 2423.07169138}
        };
        const double expected_y[isoline_count][points_per_isoline] = {
            {171.78612656, 220.38136931, 220.38136931, 265.36250878, 455.0},
            {171.7989644, 248.7449928, 248.7449928, 308.63364605, 506.38739121},
            {171.85504128, 258.19575097, 287.47240852, 355.96420961, 560.29233808},
            {172.09927987, 258.74250558, 342.56046108, 411.32270988, 618.0350615},
            {173.15, 261.02100275, 371.32673696, 484.42563591, std::nan("")}
        };

        for(int i = 0; i < isoline_count; i++) {
            for(int j = 0; j < points_per_isoline; j++) {
                if (std::isnan(expected_y[i][j])) {
                    CHECK(std::isnan(isoP[i].y[j]));
                } else {
                    CHECK_THAT(isoP[i].x[j], WithinRel(expected_x[i][j], 1e-8));
                    CHECK_THAT(isoP[i].y[j], WithinRel(expected_y[i][j], 1e-8));
                }
            }
        }
    }

    {
        // H isolines
        CoolProp::Plot::Range iso_range = plot.isoline_range(CoolProp::iHmass);
        std::vector<double> iso_values = CoolProp::Plot::Detail::generate_ranges(CoolProp::iHmass, iso_range.min, iso_range.max, isoline_count);
        auto isoH = plot.calc_isolines(CoolProp::iHmass, iso_values, points_per_isoline);
        REQUIRE(isoH.size() == isoline_count);
        CHECK_THAT(isoH[0].value, WithinAbs(75373.12689908482, 1e-10));
        CHECK_THAT(isoH[1].value, WithinAbs(200930.99391811725, 1e-10));
        CHECK_THAT(isoH[2].value, WithinAbs(326488.8609371497, 1e-10));
        CHECK_THAT(isoH[3].value, WithinAbs(452046.72795618215, 1e-10));
        CHECK_THAT(isoH[4].value, WithinAbs(577604.5949752146, 1e-10));

        const double expected_x[isoline_count][points_per_isoline] = {
            {426.00985534, 925.27531435, 1424.54077336, 1923.80623237, 2423.07169138},
            {426.00985534, 925.27531435, 1424.54077336, 1923.80623237, 2423.07169138},
            {426.00985534, 925.27531435, 1424.54077336, 1923.80623237, 2423.07169138},
            {426.00985534, 925.27531435, 1424.54077336, 1923.80623237, 2423.07169138},
            {426.00985534, 925.27531435, 1424.54077336, 1923.80623237, 2423.07169138}
        };
        const double expected_y[isoline_count][points_per_isoline] = {
            {172.17461409, std::nan(""), std::nan(""), std::nan(""), std::nan("")},
            {196.07460266, 266.63119128, std::nan(""), std::nan(""), std::nan("")},
            {213.66474036, 299.9847036, 301.72638717, std::nan(""), std::nan("")},
            {228.4112647, 322.84362137, 426.78791645, 331.52116588, 328.04216753},
            {241.56832462, 341.66140042, 458.59389239, std::nan(""), 455.0}
        };

        for(int i = 0; i < isoline_count; i++) {
            for(int j = 0; j < points_per_isoline; j++) {
                if (std::isnan(expected_y[i][j])) {
                    CHECK(std::isnan(isoH[i].y[j]));
                } else {
                    CHECK_THAT(isoH[i].x[j], WithinRel(expected_x[i][j], 1e-8));
                    CHECK_THAT(isoH[i].y[j], WithinRel(expected_y[i][j], 1e-8));
                }
            }
        }
    }

    {
        // D isolines
        CoolProp::Plot::Range iso_range = plot.isoline_range(CoolProp::iDmass);
        std::vector<double> iso_values = CoolProp::Plot::Detail::generate_ranges(CoolProp::iDmass, iso_range.min, iso_range.max, isoline_count);
        auto isoD = plot.calc_isolines(CoolProp::iDmass, iso_values, points_per_isoline);
        REQUIRE(isoD.size() == isoline_count);
        CHECK_THAT(isoD[0].value, WithinAbs(0.6749779869915704, 1e-10));
        CHECK_THAT(isoD[1].value, WithinAbs(4.704765330619733, 1e-10));
        CHECK_THAT(isoD[2].value, WithinAbs(32.793390662794806, 1e-10));
        CHECK_THAT(isoD[3].value, WithinAbs(228.57813208316188, 1e-10));
        CHECK_THAT(isoD[4].value, WithinAbs(1593.2467308391022, 1e-10));

        const double expected_x[isoline_count][points_per_isoline] = {
            {524.17387831, 1911.09303198, 2092.95299736, 2262.71394473, 2423.07169138},
            {448.10309047, 1715.11962047, 1932.46628376, 2103.15612883, 2263.90954344},
            {437.18945214, 972.48976631, 1758.36242394, 1935.75229847, 2099.20644384},
            {435.62370654, 865.94698069, 1292.02342121, 1720.27748872, 1899.3816054},
            {426.00985534, 710.87738683, 946.96900373, 1151.9181105, 1335.56535146}
        };
        const double expected_y[isoline_count][points_per_isoline] = {
            {173.15, 243.6125, 314.075, 384.5375, 455.0},
            {173.15, 243.6125, 314.075, 384.5375, 455.0},
            {173.15, 243.6125, 314.075, 384.5375, 455.0},
            {173.15, 243.6125, 314.075, 384.5375, 455.0},
            {173.15, 243.6125, 314.075, 384.5375, 455.0}
        };

        for(int i = 0; i < isoline_count; i++) {
            for(int j = 0; j < points_per_isoline; j++) {
                CHECK_THAT(isoD[i].x[j], WithinRel(expected_x[i][j], 1e-8));
                CHECK_THAT(isoD[i].y[j], WithinRel(expected_y[i][j], 1e-8));
            }
        }
    }
}

#endif
