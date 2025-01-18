#include "CoolPropPlot.h"

namespace CoolProp {
namespace Plot {

PropertyPlot::PropertyPlot(const std::string& fluid_name, CoolProp::parameters ykey, CoolProp::parameters xkey, const std::string& tp_limits)
    : fluid_name(fluid_name),
      xkey(xkey),
      ykey(ykey),
      xscale(Detail::default_scale(xkey)),
      yscale(Detail::default_scale(ykey))
{
    this->state = Detail::process_fluid_state(fluid_name);
    this->critical_state = Detail::get_critical_point(state);

    // We are just assuming that all inputs and outputs are in SI units. We
    // take care of any conversions before calling the library and after
    // getting the results.
    int out1, out2;
    axis_pair = CoolProp::generate_update_pair(xkey, 0, ykey, 1, out1, out2);
    swap_axis_inputs_for_update = (out1 == 1);

    const double HI_FACTOR = 2.25; // Upper default limits: HI_FACTOR*T_crit and HI_FACTOR*p_crit
    const double LO_FACTOR = 1.01; // Lower default limits: LO_FACTOR*T_triple and LO_FACTOR*p_triple
    if (tp_limits == "NONE")
        this->limits = {Detail::NaN, Detail::NaN, Detail::NaN, Detail::NaN};
    else if (tp_limits == "DEF")
        this->limits = {LO_FACTOR, HI_FACTOR, LO_FACTOR, HI_FACTOR};
    else if (tp_limits == "ACHP")
        this->limits = {173.15, 493.15, 0.25e5, HI_FACTOR};
    else if (tp_limits == "ORC")
        this->limits = {273.15, 673.15, 0.25e5, HI_FACTOR};
    else
        throw CoolProp::ValueError("Invalid tp_limits");

    get_axis_limits();
}

Range PropertyPlot::isoline_range(CoolProp::parameters key)
{
    if (key == CoolProp::iQ)
        return {0., 1.};
    else // TODO: always against iT?
    {
        std::vector<double> range = get_axis_limits(key, CoolProp::iT);
        return {range[0], range[1]};
    }
}

Isolines PropertyPlot::calc_isolines(CoolProp::parameters key, const std::vector<double>& values, int points) const
{
    std::vector<double> xvals = Detail::generate_values_in_range(xscale, xrange.min, xrange.max, points);
    std::vector<double> yvals = Detail::generate_values_in_range(yscale, yrange.min, yrange.max, points);

    Isolines lines;
    for (double val : values)
    {
        Isoline line(key, xkey, ykey, val, state);
        line.calc_range(xvals, yvals);
        // TODO: line.sanitize_data();
        lines.push_back(line);
    }
    return lines;
}

std::vector<CoolProp::parameters> PropertyPlot::supported_dimensions() const
{
    // taken from PropertyPlot::calc_isolines when called with iso_type='all'
    std::vector<CoolProp::parameters> keys;
    for (auto it = Detail::xy_switch.begin(); it != Detail::xy_switch.end(); ++it)
    {
        const std::map<int, Detail::IsolineSupported>& supported = it->second;
        auto supported_xy = supported.find(ykey * 10 + xkey);
        if (supported_xy != supported.end() && supported_xy->second != Detail::IsolineSupported::No)
            keys.push_back(it->first);
    }
    return keys;
}

double PropertyPlot::value_at(CoolProp::parameters key, double axis_x_value, double axis_y_value, CoolProp::phases phase) const
{
    if (key == xkey) return axis_x_value;
    if (key == ykey) return axis_y_value;

    try
    {
        if (swap_axis_inputs_for_update)
            std::swap(axis_x_value, axis_y_value);
        state->specify_phase(phase);
        state->update(axis_pair, axis_x_value, axis_y_value);
        switch (key)
        {
            case CoolProp::iT: return state->T();
            case CoolProp::iP: return state->p();
            case CoolProp::iDmass: return state->rhomass();
            case CoolProp::iHmass: return state->hmass();
            case CoolProp::iSmass: return state->smass();
            case CoolProp::iUmass: return state->umass();
            case CoolProp::iQ: return state->Q();
            default: return Detail::NaN;
        }
    }
    catch (...)
    {
        return Detail::NaN;
    }
}

Range PropertyPlot::get_sat_bounds(CoolProp::parameters key)
{
    // TODO: duplicated code from IsoLine
    double s = 1e-7;
    double t_small = critical_state->keyed_output(CoolProp::iT) * s;
    double p_small = critical_state->keyed_output(CoolProp::iP) * s;

    double t_triple = state->trivial_keyed_output(CoolProp::iT_triple);
    double t_min = state->trivial_keyed_output(CoolProp::iT_min);
    state->update(CoolProp::QT_INPUTS, 0, std::max(t_triple, t_min) + t_small);
    double fluid_min, fluid_max;
    if (key == CoolProp::iP)
    {
        fluid_min = state->keyed_output(CoolProp::iP) + p_small;
        fluid_max = critical_state->keyed_output(CoolProp::iP) - p_small;
    }
    else if (key == CoolProp::iT)
    {
        fluid_min = state->keyed_output(CoolProp::iT) + t_small;
        fluid_max = critical_state->keyed_output(CoolProp::iT) - t_small;
    }
    else
    {
        throw CoolProp::ValueError("Invalid key");
    }
    return {fluid_min, fluid_max};
}

void PropertyPlot::get_Tp_limits(double& T_lo, double& T_hi, double& P_lo, double& P_hi)
{
    T_lo = limits[0];
    T_hi = limits[1];
    P_lo = limits[2];
    P_hi = limits[3];

    double Ts_lo, Ts_hi;
    Range Ts = get_sat_bounds(CoolProp::iT);
    Ts_lo = Ts.min;
    Ts_hi = Ts.max;

    double Ps_lo, Ps_hi;
    Range Ps = get_sat_bounds(CoolProp::iP);
    Ps_lo = Ps.min;
    Ps_hi = Ps.max;

    const double ID_FACTOR = 10.0; // Values below this number are interpreted as factors
    if (std::isnan(T_lo)) T_lo = 0.0;
    else if (T_lo < ID_FACTOR) T_lo *= Ts_lo;
    if (std::isnan(T_hi)) T_hi = 1e6;
    else if (T_hi < ID_FACTOR) T_hi *= Ts_hi;
    if (std::isnan(P_lo)) P_lo = 0.0;
    else if (P_lo < ID_FACTOR) P_lo *= Ps_lo;
    if (std::isnan(P_hi)) P_hi = 1e10;
    else if (P_hi < ID_FACTOR) P_hi *= Ps_hi;

    try { T_lo = std::max(T_lo, state->trivial_keyed_output(CoolProp::iT_min)); } catch (...) {}
    try { T_hi = std::min(T_hi, state->trivial_keyed_output(CoolProp::iT_max)); } catch (...) {}
    try { P_lo = std::max(P_lo, state->trivial_keyed_output(CoolProp::iP_min)); } catch (...) {}
    try { P_hi = std::min(P_hi, state->trivial_keyed_output(CoolProp::iP_max)); } catch (...) {}
}

std::vector<double> PropertyPlot::get_axis_limits(CoolProp::parameters xkey, CoolProp::parameters ykey, bool autoscale)
{
    if (xkey == CoolProp::parameters::iundefined_parameter) xkey = this->xkey;
    if (ykey == CoolProp::parameters::iundefined_parameter) ykey = this->ykey;

    // TODO: double check comparing xkey against ykey is the same as in python
    if (xkey != this->ykey || ykey != this->ykey || autoscale)
    {
        double T_lo, T_hi, P_lo, P_hi;
        get_Tp_limits(T_lo, T_hi, P_lo, P_hi); // TODO
        std::vector<double> limits = {std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest(),
                                      std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()};
        for (double T : {T_lo, T_hi})
        {
            for (double P : {P_lo, P_hi})
            {
                try
                {
                    state->update(CoolProp::PT_INPUTS, P, T);
                    double x = state->keyed_output(xkey);
                    double y = state->keyed_output(ykey);
                    if (x < limits[0]) limits[0] = x;
                    if (x > limits[1]) limits[1] = x;
                    if (y < limits[2]) limits[2] = y;
                    if (y > limits[3]) limits[3] = y;
                }
                catch (...) { }
            }
        }
        if (xkey == this->xkey)
        {
            xrange.min = limits[0];
            xrange.max = limits[1];
        }
        if (ykey == this->ykey)
        {
            yrange.min = limits[2];
            yrange.max = limits[3];
        }
        return limits;
    }
    else
    {
        return {xrange.min, xrange.max, yrange.min, yrange.max};
    }
}

} /* namespace Plot */
} /* namespace CoolProp */

#ifdef ENABLE_CATCH
#    include <catch2/catch_all.hpp>
#    include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("Check value_at for p-h plots", "[ph_plot]") {
    CoolProp::Plot::PropertyPlot plot("R134a", CoolProp::iP, CoolProp::iHmass, "ACHP");

    CHECK_THAT(plot.value_at(CoolProp::iP, 300000/*Pa*/, 200000/*J/kg*/), WithinAbs(200000, 1e-10));
    CHECK_THAT(plot.value_at(CoolProp::iHmass, 300000, 200000), WithinAbs(300000, 1e-10));
    CHECK_THAT(plot.value_at(CoolProp::iT, 300000, 200000), WithinAbs(263.07372753976694, 1e-10));
    CHECK_THAT(plot.value_at(CoolProp::iQ, 300000, 200000), WithinAbs(0.55044347874344737, 1e-10));
}

TEST_CASE("Check that the isolines are the same as from Python", "[ph_plot]") {
    CoolProp::Plot::PropertyPlot plot("HEOS::R134a", CoolProp::iP, CoolProp::iHmass, "ACHP");
    const int isoline_count = 5;
    const int points_per_isoline = 5;

    CHECK(plot.xkey == CoolProp::iHmass);
    CHECK(plot.ykey == CoolProp::iP);

    CHECK(plot.xscale == CoolProp::Plot::Scale::Lin);
    CHECK(plot.yscale == CoolProp::Plot::Scale::Log);

    CHECK_THAT(plot.xrange.min, WithinAbs(75373.1, 1));
    CHECK_THAT(plot.xrange.max, WithinAbs(577605, 1));
    CHECK_THAT(plot.yrange.min, WithinAbs(25000, 1));
    CHECK_THAT(plot.yrange.max, WithinAbs(9.133e6, 1));

    std::vector<CoolProp::parameters> iso_types = plot.supported_dimensions();
    REQUIRE(iso_types.size() == 4);
    CHECK(iso_types[0] == CoolProp::iT);
    CHECK(iso_types[1] == CoolProp::iQ);
    CHECK(iso_types[2] == CoolProp::iDmass);
    CHECK(iso_types[3] == CoolProp::iSmass);

    {
        // Q isolines
        CoolProp::Plot::Range iso_range = plot.isoline_range(CoolProp::iQ);
        std::vector<double> iso_values = CoolProp::Plot::Detail::generate_values_in_range(CoolProp::iQ, iso_range.min, iso_range.max, isoline_count);
        CoolProp::Plot::Isolines q_isolines = plot.calc_isolines(CoolProp::iQ, iso_values, points_per_isoline);
        REQUIRE(q_isolines.size() == isoline_count);
        CHECK_THAT(q_isolines[0].value, WithinAbs(0.0, 1e-10));
        CHECK_THAT(q_isolines[1].value, WithinAbs(0.25, 1e-10));
        CHECK_THAT(q_isolines[2].value, WithinAbs(0.5, 1e-10));
        CHECK_THAT(q_isolines[3].value, WithinAbs(0.75, 1e-10));
        CHECK_THAT(q_isolines[4].value, WithinAbs(1.0, 1e-10));
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
                CHECK_THAT(q_isolines[i].x[j], WithinRel(expected_x[i][j], 1e-8));
                CHECK_THAT(q_isolines[i].y[j], WithinRel(expected_y[i][j], 1e-8));
            }
        }
    }
    {
        // T isolines
        CoolProp::Plot::Range iso_range = plot.isoline_range(CoolProp::iT);
        std::vector<double> iso_values = CoolProp::Plot::Detail::generate_values_in_range(CoolProp::iT, iso_range.min, iso_range.max, isoline_count);
        CoolProp::Plot::Isolines t_isolines = plot.calc_isolines(CoolProp::iT, iso_values, points_per_isoline);
        REQUIRE(t_isolines.size() == isoline_count);
        CHECK_THAT(t_isolines[0].value, WithinAbs(173.15, 1e-10));
        CHECK_THAT(t_isolines[1].value, WithinAbs(243.6125, 1e-10));
        CHECK_THAT(t_isolines[2].value, WithinAbs(314.075, 1e-10));
        CHECK_THAT(t_isolines[3].value, WithinAbs(384.5375, 1e-10));
        CHECK_THAT(t_isolines[4].value, WithinAbs(455.0, 1e-10));
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
                CHECK_THAT(t_isolines[i].x[j], WithinRel(expected_x[i][j], 1e-8));
                CHECK_THAT(t_isolines[i].y[j], WithinRel(expected_y[i][j], 1e-8));
            }
        }
    }
    {
        // S isolines
        CoolProp::Plot::Range iso_range = plot.isoline_range(CoolProp::iSmass);
        std::vector<double> iso_values = CoolProp::Plot::Detail::generate_values_in_range(CoolProp::iSmass, iso_range.min, iso_range.max, isoline_count);
        CoolProp::Plot::Isolines s_isolines = plot.calc_isolines(CoolProp::iSmass, iso_values, points_per_isoline);
        REQUIRE(s_isolines.size() == isoline_count);
        CHECK_THAT(s_isolines[0].value, WithinAbs(426.0098553415565, 1e-10));
        CHECK_THAT(s_isolines[1].value, WithinAbs(925.2753143522199, 1e-10));
        CHECK_THAT(s_isolines[2].value, WithinAbs(1424.540773362883, 1e-10));
        CHECK_THAT(s_isolines[3].value, WithinAbs(1923.8062323735467, 1e-10));
        CHECK_THAT(s_isolines[4].value, WithinAbs(2423.07169138421, 1e-10));
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
                if (std::isnan(s_isolines[i].x[j]))
                    CHECK(std::isnan(expected_x[i][j]));
                else
                    CHECK_THAT(s_isolines[i].x[j], WithinRel(expected_x[i][j], 1e-8));
                CHECK_THAT(s_isolines[i].y[j], WithinRel(expected_y[i][j], 1e-8));
            }
        }
    }
    {
        // D isolines
        CoolProp::Plot::Range iso_range = plot.isoline_range(CoolProp::iDmass);
        std::vector<double> iso_values = CoolProp::Plot::Detail::generate_values_in_range(CoolProp::iDmass, iso_range.min, iso_range.max, isoline_count);
        CoolProp::Plot::Isolines d_isolines = plot.calc_isolines(CoolProp::iDmass, iso_values, points_per_isoline);
        REQUIRE(d_isolines.size() == isoline_count);
        CHECK_THAT(d_isolines[0].value, WithinAbs(0.6749779869915704, 1e-10));
        CHECK_THAT(d_isolines[1].value, WithinAbs(4.704765330619733, 1e-10));
        CHECK_THAT(d_isolines[2].value, WithinAbs(32.793390662794806, 1e-10));
        CHECK_THAT(d_isolines[3].value, WithinAbs(228.57813208316188, 1e-10));
        CHECK_THAT(d_isolines[4].value, WithinAbs(1593.2467308391022, 1e-10));
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
                CHECK_THAT(d_isolines[i].x[j], WithinRel(expected_x[i][j], 1e-8));
                CHECK_THAT(d_isolines[i].y[j], WithinRel(expected_y[i][j], 1e-8));
            }
        }
    }
}

TEST_CASE("Basic TS Plot has same output as Python", "[ts_plot]") {
    CoolProp::Plot::PropertyPlot plot("HEOS::R134a", CoolProp::iT, CoolProp::iSmass, "ACHP");
    const int isoline_count = 5;
    const int points_per_isoline = 5;

    CHECK(plot.xkey == CoolProp::iSmass);
    CHECK(plot.ykey == CoolProp::iT);

    CHECK(plot.xscale == CoolProp::Plot::Scale::Lin);
    CHECK(plot.yscale == CoolProp::Plot::Scale::Lin);

    CHECK_THAT(plot.xrange.min, WithinAbs(426, 1));
    CHECK_THAT(plot.xrange.max, WithinAbs(2423, 1));
    CHECK_THAT(plot.yrange.min, WithinAbs(173, 1));
    CHECK_THAT(plot.yrange.max, WithinAbs(455, 1));

    std::vector<CoolProp::parameters> iso_types = plot.supported_dimensions();
    REQUIRE(iso_types.size() == 4);
    CHECK(iso_types[0] == CoolProp::iP);
    CHECK(iso_types[1] == CoolProp::iQ);
    CHECK(iso_types[2] == CoolProp::iDmass);
    CHECK(iso_types[3] == CoolProp::iHmass);

    {
        // Q isolines
        CoolProp::Plot::Range iso_range = plot.isoline_range(CoolProp::iQ);
        std::vector<double> iso_values = CoolProp::Plot::Detail::generate_values_in_range(CoolProp::iQ, iso_range.min, iso_range.max, isoline_count);
        CoolProp::Plot::Isolines q_isolines = plot.calc_isolines(CoolProp::iQ, iso_values, points_per_isoline);
        REQUIRE(q_isolines.size() == isoline_count);
        CHECK_THAT(q_isolines[0].value, WithinAbs(0.0, 1e-10));
        CHECK_THAT(q_isolines[1].value, WithinAbs(0.25, 1e-10));
        CHECK_THAT(q_isolines[2].value, WithinAbs(0.5, 1e-10));
        CHECK_THAT(q_isolines[3].value, WithinAbs(0.75, 1e-10));
        CHECK_THAT(q_isolines[4].value, WithinAbs(1.0, 1e-10));

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
                CHECK_THAT(q_isolines[i].x[j], WithinRel(expected_x[i][j], 1e-8));
                CHECK_THAT(q_isolines[i].y[j], WithinRel(expected_y[i][j], 1e-8));
            }
        }
    }

    {
        // P isolines
        CoolProp::Plot::Range iso_range = plot.isoline_range(CoolProp::iP);
        std::vector<double> iso_values = CoolProp::Plot::Detail::generate_values_in_range(CoolProp::iP, iso_range.min, iso_range.max, isoline_count);
        CoolProp::Plot::Isolines p_isolines = plot.calc_isolines(CoolProp::iP, iso_values, points_per_isoline);
        REQUIRE(p_isolines.size() == isoline_count);
        CHECK_THAT(p_isolines[0].value, WithinAbs(25000.000000000007, 1e-8));
        CHECK_THAT(p_isolines[1].value, WithinAbs(109297.03270763098, 1e-8));
        CHECK_THAT(p_isolines[2].value, WithinAbs(477833.65434771683, 1e-8));
        CHECK_THAT(p_isolines[3].value, WithinAbs(2089032.0219219688, 1e-8));
        CHECK_THAT(p_isolines[4].value, WithinAbs(9133000.049091753, 1e-8));

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
                    CHECK(std::isnan(p_isolines[i].y[j]));
                } else {
                    CHECK_THAT(p_isolines[i].x[j], WithinRel(expected_x[i][j], 1e-8));
                    CHECK_THAT(p_isolines[i].y[j], WithinRel(expected_y[i][j], 1e-8));
                }
            }
        }
    }

    {
        // H isolines
        CoolProp::Plot::Range iso_range = plot.isoline_range(CoolProp::iHmass);
        std::vector<double> iso_values = CoolProp::Plot::Detail::generate_values_in_range(CoolProp::iHmass, iso_range.min, iso_range.max, isoline_count);
        CoolProp::Plot::Isolines h_isolines = plot.calc_isolines(CoolProp::iHmass, iso_values, points_per_isoline);
        REQUIRE(h_isolines.size() == isoline_count);
        CHECK_THAT(h_isolines[0].value, WithinAbs(75373.12689908482, 1e-10));
        CHECK_THAT(h_isolines[1].value, WithinAbs(200930.99391811725, 1e-10));
        CHECK_THAT(h_isolines[2].value, WithinAbs(326488.8609371497, 1e-10));
        CHECK_THAT(h_isolines[3].value, WithinAbs(452046.72795618215, 1e-10));
        CHECK_THAT(h_isolines[4].value, WithinAbs(577604.5949752146, 1e-10));

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
                    CHECK(std::isnan(h_isolines[i].y[j]));
                } else {
                    CHECK_THAT(h_isolines[i].x[j], WithinRel(expected_x[i][j], 1e-8));
                    CHECK_THAT(h_isolines[i].y[j], WithinRel(expected_y[i][j], 1e-8));
                }
            }
        }
    }

    {
        // D isolines
        CoolProp::Plot::Range iso_range = plot.isoline_range(CoolProp::iDmass);
        std::vector<double> iso_values = CoolProp::Plot::Detail::generate_values_in_range(CoolProp::iDmass, iso_range.min, iso_range.max, isoline_count);
        CoolProp::Plot::Isolines d_isolines = plot.calc_isolines(CoolProp::iDmass, iso_values, points_per_isoline);
        REQUIRE(d_isolines.size() == isoline_count);
        CHECK_THAT(d_isolines[0].value, WithinAbs(0.6749779869915704, 1e-10));
        CHECK_THAT(d_isolines[1].value, WithinAbs(4.704765330619733, 1e-10));
        CHECK_THAT(d_isolines[2].value, WithinAbs(32.793390662794806, 1e-10));
        CHECK_THAT(d_isolines[3].value, WithinAbs(228.57813208316188, 1e-10));
        CHECK_THAT(d_isolines[4].value, WithinAbs(1593.2467308391022, 1e-10));

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
                CHECK_THAT(d_isolines[i].x[j], WithinRel(expected_x[i][j], 1e-8));
                CHECK_THAT(d_isolines[i].y[j], WithinRel(expected_y[i][j], 1e-8));
            }
        }
    }
}

#endif
