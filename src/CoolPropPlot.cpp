#include "CoolPropPlot.h"
#include "CoolProp.h"
#include "CPnumerics.h"
#include <map>

namespace CoolProp {
namespace Plot {

namespace Detail {

const double NaN = std::numeric_limits<double>::quiet_NaN();

const int TS = CoolProp::iT * 10 + CoolProp::iSmass;
const int PH = CoolProp::iP * 10 + CoolProp::iHmass;
const int HS = CoolProp::iHmass * 10 + CoolProp::iSmass;
const int PS = CoolProp::iP * 10 + CoolProp::iSmass;
const int PD = CoolProp::iP * 10 + CoolProp::iDmass;
const int TD = CoolProp::iT * 10 + CoolProp::iDmass;
const int PT = CoolProp::iP * 10 + CoolProp::iT;

enum IsolineSupported
{
    No = 0,
    Yes = 1,
    Flipped = 2
};

const std::map<CoolProp::parameters, std::map<int, IsolineSupported>> xy_switch = {
    {CoolProp::iDmass, {{TS, Flipped}, {PH, Flipped}, {HS, Yes    }, {PS, Flipped}, {PD, No     }, {TD, No     }, {PT, Yes    }}},
    {CoolProp::iHmass, {{TS, Yes    }, {PH, No     }, {HS, No     }, {PS, Flipped}, {PD, Flipped}, {TD, Yes    }, {PT, Yes    }}},
    {CoolProp::iP,     {{TS, Yes    }, {PH, No     }, {HS, Yes    }, {PS, No     }, {PD, No     }, {TD, Yes    }, {PT, No     }}},
    {CoolProp::iSmass, {{TS, No     }, {PH, Flipped}, {HS, No     }, {PS, No     }, {PD, Flipped}, {TD, Yes    }, {PT, Flipped}}},
    {CoolProp::iT,     {{TS, No     }, {PH, Flipped}, {HS, Yes    }, {PS, Yes    }, {PD, Yes    }, {TD, No     }, {PT, No     }}},
    {CoolProp::iQ,     {{TS, Flipped}, {PH, Flipped}, {HS, Flipped}, {PS, Flipped}, {PD, Flipped}, {TD, Flipped}, {PT, Yes    }}}
};

Scale default_scale(CoolProp::parameters key) {
    switch (key) {
        case CoolProp::iDmass: return Scale::Log;
        case CoolProp::iHmass: return Scale::Lin;
        case CoolProp::iP:     return Scale::Log;
        case CoolProp::iSmass: return Scale::Lin;
        case CoolProp::iT:     return Scale::Lin;
        case CoolProp::iUmass: return Scale::Lin;
        case CoolProp::iQ:     return Scale::Lin;
        default:               return Scale::Lin;
    }
}

inline std::shared_ptr<CoolProp::AbstractState> process_fluid_state(const std::string& fluid_ref) {
    std::string backend;
    std::string fluids;
    CoolProp::extract_backend(fluid_ref, backend, fluids);
    std::vector<double> fractions;
    fluids = CoolProp::extract_fractions(fluids, fractions);

    return std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory(backend, fluids));
}

std::shared_ptr<CoolProp::AbstractState> get_critical_point(const std::shared_ptr<CoolProp::AbstractState>& state) {
    CoolProp::CriticalState crit_state;
    crit_state.T = Detail::NaN;
    crit_state.p = Detail::NaN;
    crit_state.rhomolar = Detail::NaN;
    crit_state.rhomolar = Detail::NaN;
    crit_state.stable = false;
    try {
        crit_state.T = state->T_critical();
        crit_state.p = state->p_critical();
        crit_state.rhomolar = state->rhomolar_critical();
        crit_state.stable = true;
    } catch (...) {
        try {
            for (CoolProp::CriticalState crit_state_tmp: state->all_critical_points()) {
                if (crit_state_tmp.stable && (crit_state_tmp.T > crit_state.T || !std::isfinite(crit_state.T))) {
                    crit_state.T = crit_state_tmp.T;
                    crit_state.p = crit_state_tmp.p;
                    crit_state.rhomolar = crit_state_tmp.rhomolar;
                    crit_state.stable = crit_state_tmp.stable;
                }
            }
        } catch (...) {
            throw CoolProp::ValueError("Could not calculate the critical point data.");
        }
    }

    std::shared_ptr<CoolProp::AbstractState> new_state(CoolProp::AbstractState::factory(state->backend_name(), state->fluid_names()));
    std::vector<double> masses = state->get_mass_fractions();
    if (masses.size() > 1)
        new_state->set_mass_fractions(masses);

    if (std::isfinite(crit_state.p) && std::isfinite(crit_state.T)) {
        try {
            new_state->specify_phase(CoolProp::iphase_critical_point);
            new_state->update(CoolProp::PT_INPUTS, crit_state.p, crit_state.T);
            return new_state;
        } catch (...) { }
        try {
            new_state->update(CoolProp::PT_INPUTS, crit_state.p, crit_state.T);
            return new_state;
        } catch (...) { }
    }

    if (std::isfinite(crit_state.rhomolar) && std::isfinite(crit_state.T)) {
        try {
            new_state->specify_phase(CoolProp::iphase_critical_point);
            new_state->update(CoolProp::DmolarT_INPUTS, crit_state.rhomolar, crit_state.T);
            return new_state;
        } catch (...) { }
        try {
            new_state->update(CoolProp::DmolarT_INPUTS, crit_state.rhomolar, crit_state.T);
            return new_state;
        } catch (...) { }
    }
    throw CoolProp::ValueError("Could not calculate the critical point data.");
    return nullptr;
}

} /* namespace Detail */


std::vector<double> generate_values_in_range(Scale scale, const Range& range, int count) {
    if (scale == Scale::Log)
        return logspace(range.min, range.max, count);
    else
        return linspace(range.min, range.max, count);
}

std::vector<double> generate_values_in_range(CoolProp::parameters type, const Range& range, int count) {
    return generate_values_in_range(Detail::default_scale(type), range, count);
}


Isoline::Isoline(CoolProp::parameters key, CoolProp::parameters xkey, CoolProp::parameters ykey, double value, const std::shared_ptr<CoolProp::AbstractState>& state)
    : key_(key),
      xkey_(xkey),
      ykey_(ykey),
      value(value),
      state_(state) {
    this->critical_state_ = Detail::get_critical_point(state);
}

Range Isoline::get_sat_bounds(CoolProp::parameters key) const {
    double s = 1e-7;
    double t_small = critical_state_->keyed_output(CoolProp::iT) * s;
    double p_small = critical_state_->keyed_output(CoolProp::iP) * s;

    double t_triple = state_->trivial_keyed_output(CoolProp::iT_triple);
    double t_min;
    try {
        t_min = state_->trivial_keyed_output(CoolProp::iT_min);
    } catch (...) {
        t_min = t_triple;
    }
    state_->update(CoolProp::QT_INPUTS, 0, std::max(t_triple, t_min) + t_small);
    if (key == CoolProp::iP)
        return {state_->keyed_output(CoolProp::iP) + p_small, critical_state_->keyed_output(CoolProp::iP) - p_small};
    else if (key == CoolProp::iT)
        return {state_->keyed_output(CoolProp::iT) + t_small, critical_state_->keyed_output(CoolProp::iT) - t_small};
    else
        throw CoolProp::ValueError("Invalid key");
}

void Isoline::calc_sat_range(int count) {
    Range t = get_sat_bounds(CoolProp::iT);
    std::vector<double> two = ::linspace(t.min, t.max, count);
    std::vector<double> one(two.size(), value);
    CoolProp::input_pairs input_pair = CoolProp::QT_INPUTS;

    double t_crit = critical_state_->keyed_output(CoolProp::iT);
    double p_crit = critical_state_->keyed_output(CoolProp::iP);
    double x_crit = critical_state_->keyed_output(xkey_);
    double y_crit = critical_state_->keyed_output(ykey_);
    x.resize(one.size());
    y.resize(one.size());
    for (int i = 0; i < one.size(); ++i) {
        try {
            state_->update(input_pair, one[i], two[i]);
            x[i] = state_->keyed_output(xkey_);
            y[i] = state_->keyed_output(ykey_);
        } catch (...) {
            if ((input_pair == CoolProp::QT_INPUTS && abs(two[i] - t_crit) < 1e0)
             || (input_pair == CoolProp::PQ_INPUTS && abs(one[i] - p_crit) < 1e2)) {
                x[i] = x_crit;
                y[i] = y_crit;
                std::cerr << "ERROR near critical inputs" << std::endl;
            } else {
                x[i] = Detail::NaN;
                y[i] = Detail::NaN;
                std::cerr << "ERROR" << std::endl;
            }
        }
    }
}

void Isoline::update_pair(int& ipos, int& xpos, int& ypos, int& pair) {
    Detail::IsolineSupported should_switch = Detail::xy_switch.at(key_).at(ykey_ * 10 + xkey_);
    double out1, out2;
    if (should_switch == Detail::IsolineSupported::No)
        throw CoolProp::ValueError("This isoline cannot be calculated!");
    else if (should_switch == Detail::IsolineSupported::Yes)
        pair = CoolProp::generate_update_pair(key_, 0.0, xkey_, 1.0, out1, out2);
    else if (should_switch == Detail::IsolineSupported::Flipped)
        pair = CoolProp::generate_update_pair(key_, 0.0, ykey_, 1.0, out1, out2);

    bool should_swap = (out1 != 0.0);

    if (should_switch == Detail::IsolineSupported::Yes && !should_swap) {
        ipos = 0;
        xpos = 1;
        ypos = 2;
    } else if (should_switch == Detail::IsolineSupported::Flipped && !should_swap) {
        ipos = 0;
        xpos = 2;
        ypos = 1;
    } else if (should_switch == Detail::IsolineSupported::Yes && should_swap) {
        ipos = 1;
        xpos = 0;
        ypos = 2;
    } else if (should_switch == Detail::IsolineSupported::Flipped && should_swap) {
        ipos = 1;
        xpos = 2;
        ypos = 0;
    } else {
        throw CoolProp::ValueError("Check the code, this should not happen!");
    }
}

void Isoline::calc_range(std::vector<double>& xvals, std::vector<double>& yvals) {
    if (key_ == CoolProp::iQ) {
        calc_sat_range(static_cast<int>(xvals.size()));
    } else {
        int ipos, xpos, ypos, pair;
        update_pair(ipos, xpos, ypos, pair);

        std::vector<double> ivals(xvals.size(), value);
        std::vector<int> order = {ipos, xpos, ypos};
        std::vector<CoolProp::parameters> idxs(3);
        idxs[ipos] = key_;
        idxs[xpos] = xkey_;
        idxs[ypos] = ykey_;
        std::vector<std::vector<double>> vals(3);
        vals[ipos] = ivals;
        vals[xpos] = xvals;
        vals[ypos] = yvals;

        for (int i = 0; i < vals[2].size(); ++i) {
            try {
                state_->update((CoolProp::input_pairs)pair, vals[0][i], vals[1][i]);
                vals[2][i] = state_->keyed_output(idxs[2]);
            } catch (...) {
                vals[2][i] = Detail::NaN;
            }
        }

        for (int i = 0; i < idxs.size(); ++i) {
            if (idxs[i] == xkey_) x = vals[i];
            if (idxs[i] == ykey_) y = vals[i];
        }
    }
}

PropertyPlot::PropertyPlot(const std::string& fluid_name, CoolProp::parameters ykey, CoolProp::parameters xkey, TPLimits tp_limits)
    : xkey_(xkey),
      ykey_(ykey) {
    this->state_ = Detail::process_fluid_state(fluid_name);
    this->critical_state_ = Detail::get_critical_point(state_);

    xaxis.scale = Detail::default_scale(xkey);
    yaxis.scale = Detail::default_scale(ykey);

    // We are just assuming that all inputs and outputs are in SI units. We
    // take care of any conversions before calling the library and after
    // getting the results.
    int out1, out2;
    axis_pair_ = CoolProp::generate_update_pair(xkey, 0, ykey, 1, out1, out2);
    swap_axis_inputs_for_update_ = (out1 == 1);

    const double HI_FACTOR = 2.25; // Upper default limits: HI_FACTOR*T_crit and HI_FACTOR*p_crit
    const double LO_FACTOR = 1.01; // Lower default limits: LO_FACTOR*T_triple and LO_FACTOR*p_triple
    switch (tp_limits) {
        case TPLimits::None: this->Tp_limits_ = {{Detail::NaN, Detail::NaN}, {Detail::NaN, Detail::NaN}}; break;
        case TPLimits::Def:  this->Tp_limits_ = {{LO_FACTOR, HI_FACTOR}, {LO_FACTOR, HI_FACTOR}}; break;
        case TPLimits::Achp: this->Tp_limits_ = {{173.15, 493.15}, {0.25e5, HI_FACTOR}}; break;
        case TPLimits::Orc:  this->Tp_limits_ = {{273.15, 673.15}, {0.25e5, HI_FACTOR}}; break;
    }

    Range2D ranges = get_axis_limits();
    xaxis.range = ranges.x;
    yaxis.range = ranges.y;
}

Range PropertyPlot::isoline_range(CoolProp::parameters key) const {
    if (key == CoolProp::iQ)
        return {0, 1};
    else
        return get_axis_limits(key, CoolProp::iT).x;
}

Isolines PropertyPlot::calc_isolines(CoolProp::parameters key, const std::vector<double>& values, int points) const {
    std::vector<double> xvals = generate_values_in_range(xaxis.scale, xaxis.range, points);
    std::vector<double> yvals = generate_values_in_range(yaxis.scale, yaxis.range, points);

    Isolines lines;
    for (double val : values) {
        Isoline line(key, xkey_, ykey_, val, state_);
        line.calc_range(xvals, yvals);
        lines.push_back(line);
    }
    return lines;
}

std::vector<CoolProp::parameters> PropertyPlot::supported_isoline_keys() const {
    // taken from PropertyPlot::calc_isolines when called with iso_type='all'
    std::vector<CoolProp::parameters> keys;
    for (auto it = Detail::xy_switch.begin(); it != Detail::xy_switch.end(); ++it) {
        const std::map<int, Detail::IsolineSupported>& supported = it->second;
        auto supported_xy = supported.find(ykey_ * 10 + xkey_);
        if (supported_xy != supported.end() && supported_xy->second != Detail::IsolineSupported::No)
            keys.push_back(it->first);
    }
    return keys;
}

double PropertyPlot::value_at(CoolProp::parameters key, double xvalue, double yvalue, CoolProp::phases phase) const {
    if (key == xkey_) return xvalue;
    if (key == ykey_) return yvalue;

    try {
        if (swap_axis_inputs_for_update_)
            std::swap(xvalue, yvalue);
        state_->specify_phase(phase);
        state_->update(axis_pair_, xvalue, yvalue);
        switch (key) {
            case CoolProp::iT: return state_->T();
            case CoolProp::iP: return state_->p();
            case CoolProp::iDmass: return state_->rhomass();
            case CoolProp::iHmass: return state_->hmass();
            case CoolProp::iSmass: return state_->smass();
            case CoolProp::iUmass: return state_->umass();
            case CoolProp::iQ: return state_->Q();
            default: return Detail::NaN;
        }
    } catch (...) {
        return Detail::NaN;
    }
}

Range PropertyPlot::get_sat_bounds(CoolProp::parameters key) const {
    double s = 1e-7;
    double t_small = critical_state_->keyed_output(CoolProp::iT) * s;
    double p_small = critical_state_->keyed_output(CoolProp::iP) * s;

    double t_triple = state_->trivial_keyed_output(CoolProp::iT_triple);
    double t_min;
    try {
        t_min = state_->trivial_keyed_output(CoolProp::iT_min);
    } catch (...) {
        t_min = t_triple;
    }
    state_->update(CoolProp::QT_INPUTS, 0, std::max(t_triple, t_min) + t_small);
    if (key == CoolProp::iP)
        return {state_->keyed_output(CoolProp::iP) + p_small, critical_state_->keyed_output(CoolProp::iP) - p_small};
    else if (key == CoolProp::iT)
        return {state_->keyed_output(CoolProp::iT) + t_small, critical_state_->keyed_output(CoolProp::iT) - t_small};
    else
        throw CoolProp::ValueError("Invalid key");
}

PropertyPlot::Range2D PropertyPlot::get_Tp_limits() const {
    Range t = Tp_limits_.T;
    Range p = Tp_limits_.p;
    Range tsat = get_sat_bounds(CoolProp::iT);
    Range psat = get_sat_bounds(CoolProp::iP);

    const double ID_FACTOR = 10.0; // Values below this number are interpreted as factors
    if (std::isnan(t.min)) t.min = 0.0;
    else if (t.min < ID_FACTOR) t.min *= tsat.min;
    if (std::isnan(t.max)) t.max = 1e6;
    else if (t.max < ID_FACTOR) t.max *= tsat.max;
    if (std::isnan(p.min)) p.min = 0.0;
    else if (p.min < ID_FACTOR) p.min *= psat.min;
    if (std::isnan(p.max)) p.max = 1e10;
    else if (p.max < ID_FACTOR) p.max *= psat.max;

    try { t.min = std::max(t.min, state_->trivial_keyed_output(CoolProp::iT_min)); } catch (...) {}
    try { t.max = std::min(t.max, state_->trivial_keyed_output(CoolProp::iT_max)); } catch (...) {}
    try { p.min = std::max(p.min, state_->trivial_keyed_output(CoolProp::iP_min)); } catch (...) {}
    try { p.max = std::min(p.max, state_->trivial_keyed_output(CoolProp::iP_max)); } catch (...) {}
    return {t, p};
}

PropertyPlot::Range2D PropertyPlot::get_axis_limits(CoolProp::parameters xkey, CoolProp::parameters ykey, bool autoscale) const {
    if (xkey == CoolProp::parameters::iundefined_parameter) xkey = this->xkey_;
    if (ykey == CoolProp::parameters::iundefined_parameter) ykey = this->ykey_;

    if (xkey != this->xkey_ || ykey != this->ykey_ || autoscale) {
        Range2D tp_limits = get_Tp_limits();
        Range xrange = {std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()};
        Range yrange = {std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()};

        for (double T : {tp_limits.T.min, tp_limits.T.max}) {
            for (double p : {tp_limits.p.min, tp_limits.p.max}) {
                try {
                    state_->update(CoolProp::PT_INPUTS, p, T);
                    double x = state_->keyed_output(xkey);
                    double y = state_->keyed_output(ykey);
                    xrange.min = std::min(xrange.min, x);
                    xrange.max = std::max(xrange.max, x);
                    yrange.min = std::min(yrange.min, y);
                    yrange.max = std::max(yrange.max, y);
                } catch (...) { }
            }
        }
        return {xrange, yrange};
    } else {
        return {xaxis.range, yaxis.range};
    }
}

} /* namespace Plot */
} /* namespace CoolProp */


#ifdef ENABLE_CATCH
#    include <catch2/catch_all.hpp>
#    include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("Check value_at for p-h plots", "[Plot]") {
    CoolProp::Plot::PropertyPlot plot("R134a", CoolProp::iP, CoolProp::iHmass, CoolProp::Plot::TPLimits::Achp);

    CHECK_THAT(plot.value_at(CoolProp::iP, 300000/*Pa*/, 200000/*J/kg*/), WithinAbs(200000, 1e-10));
    CHECK_THAT(plot.value_at(CoolProp::iHmass, 300000, 200000), WithinAbs(300000, 1e-10));
    CHECK_THAT(plot.value_at(CoolProp::iT, 300000, 200000), WithinAbs(263.07372753976694, 1e-10));
    CHECK_THAT(plot.value_at(CoolProp::iQ, 300000, 200000), WithinAbs(0.55044347874344737, 1e-10));
}

TEST_CASE("Check that the isolines are the same as from Python", "[Plot]") {
    CoolProp::Plot::PropertyPlot plot("HEOS::R134a", CoolProp::iP, CoolProp::iHmass, CoolProp::Plot::TPLimits::Achp);
    const int isoline_count = 5;
    const int points_per_isoline = 5;

    // CHECK(plot.xkey_ == CoolProp::iHmass);
    // CHECK(plot.ykey_ == CoolProp::iP);

    CHECK(plot.xaxis.scale == CoolProp::Plot::Scale::Lin);
    CHECK(plot.yaxis.scale == CoolProp::Plot::Scale::Log);
    CHECK_THAT(plot.xaxis.min, WithinAbs(75373.1, 1));
    CHECK_THAT(plot.xaxis.max, WithinAbs(577605, 1));
    CHECK_THAT(plot.yaxis.min, WithinAbs(25000, 1));
    CHECK_THAT(plot.yaxis.max, WithinAbs(9.133e6, 1));

    std::vector<CoolProp::parameters> iso_types = plot.supported_isoline_keys();
    REQUIRE(iso_types.size() == 4);
    CHECK(iso_types[0] == CoolProp::iT);
    CHECK(iso_types[1] == CoolProp::iQ);
    CHECK(iso_types[2] == CoolProp::iDmass);
    CHECK(iso_types[3] == CoolProp::iSmass);

    {
        // Q isolines
        CoolProp::Plot::Range q_range = plot.isoline_range(CoolProp::iQ);
        std::vector<double> q_values = CoolProp::Plot::generate_values_in_range(CoolProp::iQ, q_range, isoline_count);
        CoolProp::Plot::Isolines q_isolines = plot.calc_isolines(CoolProp::iQ, q_values, points_per_isoline);
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
        for (int i = 0; i < q_isolines.size(); ++i) {
            REQUIRE(q_isolines[i].size() == points_per_isoline);
            for (int j = 0; j < q_isolines[i].size(); ++j) {
                CHECK_THAT(q_isolines[i].x[j], WithinRel(expected_x[i][j], 1e-8));
                CHECK_THAT(q_isolines[i].y[j], WithinRel(expected_y[i][j], 1e-8));
            }
        }
    }
    {
        // T isolines
        CoolProp::Plot::Range t_range = plot.isoline_range(CoolProp::iT);
        std::vector<double> t_values = CoolProp::Plot::generate_values_in_range(CoolProp::iT, t_range, isoline_count);
        CoolProp::Plot::Isolines t_isolines = plot.calc_isolines(CoolProp::iT, t_values, points_per_isoline);
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
        for (int i = 0; i < t_isolines.size(); ++i) {
            REQUIRE(t_isolines[i].size() == points_per_isoline);
            for (int j = 0; j < t_isolines[i].size(); ++j) {
                CHECK_THAT(t_isolines[i].x[j], WithinRel(expected_x[i][j], 1e-8));
                CHECK_THAT(t_isolines[i].y[j], WithinRel(expected_y[i][j], 1e-8));
            }
        }
    }
    {
        // S isolines
        CoolProp::Plot::Range s_range = plot.isoline_range(CoolProp::iSmass);
        std::vector<double> s_values = CoolProp::Plot::generate_values_in_range(CoolProp::iSmass, s_range, isoline_count);
        CoolProp::Plot::Isolines s_isolines = plot.calc_isolines(CoolProp::iSmass, s_values, points_per_isoline);
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
        for (int i = 0; i < s_isolines.size(); ++i) {
            REQUIRE(s_isolines[i].size() == points_per_isoline);
            for (int j = 0; j < s_isolines[i].size(); ++j) {
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
        CoolProp::Plot::Range d_range = plot.isoline_range(CoolProp::iDmass);
        std::vector<double> d_values = CoolProp::Plot::generate_values_in_range(CoolProp::iDmass, d_range, isoline_count);
        CoolProp::Plot::Isolines d_isolines = plot.calc_isolines(CoolProp::iDmass, d_values, points_per_isoline);
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
        for (int i = 0; i < d_isolines.size(); ++i) {
            REQUIRE(d_isolines[i].size() == points_per_isoline);
            for (int j = 0; j < d_isolines[i].size(); ++j) {
                CHECK_THAT(d_isolines[i].x[j], WithinRel(expected_x[i][j], 1e-8));
                CHECK_THAT(d_isolines[i].y[j], WithinRel(expected_y[i][j], 1e-8));
            }
        }
    }
}

TEST_CASE("Basic TS Plot has same output as Python", "[Plot]") {
    CoolProp::Plot::PropertyPlot plot("HEOS::R134a", CoolProp::iT, CoolProp::iSmass, CoolProp::Plot::TPLimits::Achp);
    const int isoline_count = 5;
    const int points_per_isoline = 5;

    // CHECK(plot.xkey_ == CoolProp::iSmass);
    // CHECK(plot.ykey_ == CoolProp::iT);

    CHECK(plot.xaxis.scale == CoolProp::Plot::Scale::Lin);
    CHECK(plot.yaxis.scale == CoolProp::Plot::Scale::Lin);
    CHECK_THAT(plot.xaxis.min, WithinAbs(426, 1));
    CHECK_THAT(plot.xaxis.max, WithinAbs(2423, 1));
    CHECK_THAT(plot.yaxis.min, WithinAbs(173, 1));
    CHECK_THAT(plot.yaxis.max, WithinAbs(455, 1));

    std::vector<CoolProp::parameters> iso_types = plot.supported_isoline_keys();
    REQUIRE(iso_types.size() == 4);
    CHECK(iso_types[0] == CoolProp::iP);
    CHECK(iso_types[1] == CoolProp::iQ);
    CHECK(iso_types[2] == CoolProp::iDmass);
    CHECK(iso_types[3] == CoolProp::iHmass);

    {
        // Q isolines
        CoolProp::Plot::Range q_range = plot.isoline_range(CoolProp::iQ);
        std::vector<double> q_values = CoolProp::Plot::generate_values_in_range(CoolProp::iQ, q_range, isoline_count);
        CoolProp::Plot::Isolines q_isolines = plot.calc_isolines(CoolProp::iQ, q_values, points_per_isoline);
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

        for(int i = 0; i < q_isolines.size(); i++) {
            REQUIRE(q_isolines[i].size() == points_per_isoline);
            for(int j = 0; j < q_isolines[i].size(); j++) {
                CHECK_THAT(q_isolines[i].x[j], WithinRel(expected_x[i][j], 1e-8));
                CHECK_THAT(q_isolines[i].y[j], WithinRel(expected_y[i][j], 1e-8));
            }
        }
    }
    {
        // P isolines
        CoolProp::Plot::Range p_range = plot.isoline_range(CoolProp::iP);
        std::vector<double> p_values = CoolProp::Plot::generate_values_in_range(CoolProp::iP, p_range, isoline_count);
        CoolProp::Plot::Isolines p_isolines = plot.calc_isolines(CoolProp::iP, p_values, points_per_isoline);
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

        for(int i = 0; i < p_isolines.size(); i++) {
            REQUIRE(p_isolines[i].size() == points_per_isoline);
            for(int j = 0; j < p_isolines[i].size(); j++) {
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
        CoolProp::Plot::Range h_range = plot.isoline_range(CoolProp::iHmass);
        std::vector<double> h_values = CoolProp::Plot::generate_values_in_range(CoolProp::iHmass, h_range, isoline_count);
        CoolProp::Plot::Isolines h_isolines = plot.calc_isolines(CoolProp::iHmass, h_values, points_per_isoline);
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

        for(int i = 0; i < h_isolines.size(); i++) {
            REQUIRE(h_isolines[i].size() == points_per_isoline);
            for(int j = 0; j < h_isolines[i].size(); j++) {
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
        CoolProp::Plot::Range d_range = plot.isoline_range(CoolProp::iDmass);
        std::vector<double> d_values = CoolProp::Plot::generate_values_in_range(CoolProp::iDmass, d_range, isoline_count);
        CoolProp::Plot::Isolines d_isolines = plot.calc_isolines(CoolProp::iDmass, d_values, points_per_isoline);
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

        for(int i = 0; i < d_isolines.size(); i++) {
            REQUIRE(d_isolines[i].size() == points_per_isoline);
            for(int j = 0; j < d_isolines[i].size(); j++) {
                CHECK_THAT(d_isolines[i].x[j], WithinRel(expected_x[i][j], 1e-8));
                CHECK_THAT(d_isolines[i].y[j], WithinRel(expected_y[i][j], 1e-8));
            }
        }
    }
}

#endif /* ENABLE_CATCH */
