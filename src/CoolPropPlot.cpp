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

// The following block of code was auto-generated and inserted by the
// script dev/scripts/generate_Plot_test_data.py.
//
// <autogenerated>
TEST_CASE("Check value_at for p-h plots", "[Plot]") {
    CoolProp::Plot::PropertyPlot plot("R134a", CoolProp::iP, CoolProp::iHmass, CoolProp::Plot::TPLimits::Achp);

    CHECK_THAT(plot.value_at(CoolProp::iP, 300000/*Pa*/, 200000/*J/kg*/), WithinAbs(200000, 1e-10));
    CHECK_THAT(plot.value_at(CoolProp::iHmass, 300000, 200000), WithinAbs(300000, 1e-10));
    CHECK_THAT(plot.value_at(CoolProp::iT, 300000, 200000), WithinAbs(263.0737275397678, 1e-10));
    CHECK_THAT(plot.value_at(CoolProp::iQ, 300000, 200000), WithinAbs(0.5504434787434432, 1e-10));
}
TEST_CASE("Check that the isolines are the same as from Python", "[Plot]") {
    CoolProp::Plot::PropertyPlot plot("HEOS::R134a", CoolProp::iP, CoolProp::iHmass, CoolProp::Plot::TPLimits::Achp);
    const int isoline_count = 5;
    const int points_per_isoline = 5;

    // CHECK(plot.xkey_ == CoolProp::iHmass);
    // CHECK(plot.ykey_ == CoolProp::iP);

    CHECK(plot.xaxis.scale == CoolProp::Plot::Scale::Lin);
    CHECK(plot.yaxis.scale == CoolProp::Plot::Scale::Log);
    CHECK_THAT(plot.xaxis.min, WithinAbs(75373.12689908482, 1));
    CHECK_THAT(plot.xaxis.max, WithinAbs(577604.5949752146, 1));
    CHECK_THAT(plot.yaxis.min, WithinAbs(25000.0, 1));
    CHECK_THAT(plot.yaxis.max, WithinAbs(9133370.875847604, 1));

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
            {71455.0825704527, 132940.602012992, 198498.370551912, 271578.877763124, 389490.979699808},
            {137326.831168219, 191267.585241559, 248361.039003664, 309540.80583791, 389563.709352125},
            {203198.579765986, 249594.568470126, 298223.707455415, 347502.733912697, 389636.439004441},
            {269070.328363753, 307921.551698693, 348086.375907167, 385464.661987484, 389709.168656758},
            {334942.07696152, 366248.53492726, 397949.044358919, 423426.59006227, 389781.898309075},
        };
        const double expected_y[isoline_count][points_per_isoline] = {
            {389.56705952134, 25851.3343934178, 281115.856001781, 1316960.5263817, 4059273.23696491},
            {389.56705952134, 25851.3343934178, 281115.856001781, 1316960.5263817, 4059273.23696491},
            {389.56705952134, 25851.3343934178, 281115.856001781, 1316960.5263817, 4059273.23696491},
            {389.56705952134, 25851.3343934178, 281115.856001781, 1316960.5263817, 4059273.23696491},
            {389.56705952134, 25851.3343934178, 281115.856001781, 1316960.5263817, 4059273.23696491},
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
        CHECK_THAT(t_isolines[2].value, WithinAbs(314.07500000000005, 1e-10));
        CHECK_THAT(t_isolines[3].value, WithinAbs(384.5375, 1e-10));
        CHECK_THAT(t_isolines[4].value, WithinAbs(455.0, 1e-10));
        const double expected_x[isoline_count][points_per_isoline] = {
            {75373.1268990848, 75410.9911120345, 75576.5817006844, 76301.4918515715, 79487.8877883422},
            {382785.230587559, 161389.442353423, 161516.218619848, 162076.984158624, 164637.062377748},
            {439466.649843277, 438148.172824179, 431912.0662387, 257605.319479605, 257512.839247251},
            {504550.626065608, 503783.529360532, 500331.593280543, 482707.178360249, 366958.520785585},
            {577604.594975215, 577097.065048065, 574850.152315662, 564443.789731467, 507875.800635261},
        };
        const double expected_y[isoline_count][points_per_isoline] = {
            {25000, 109298.142136262, 477843.354977538, 2089095.63724813, 9133370.87584761},
            {25000, 109298.142136262, 477843.354977538, 2089095.63724813, 9133370.87584761},
            {25000, 109298.142136262, 477843.354977538, 2089095.63724813, 9133370.87584761},
            {25000, 109298.142136262, 477843.354977538, 2089095.63724813, 9133370.87584761},
            {25000, 109298.142136262, 477843.354977538, 2089095.63724813, 9133370.87584761},
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
        CHECK_THAT(s_isolines[0].value, WithinAbs(426.00948386039755, 1e-10));
        CHECK_THAT(s_isolines[1].value, WithinAbs(925.2750357413507, 1e-10));
        CHECK_THAT(s_isolines[2].value, WithinAbs(1424.540587622304, 1e-10));
        CHECK_THAT(s_isolines[3].value, WithinAbs(1923.806139503257, 1e-10));
        CHECK_THAT(s_isolines[4].value, WithinAbs(2423.07169138421, 1e-10));
        const double expected_x[isoline_count][points_per_isoline] = {
            {73758.1368335347, 73811.2861613466, 74043.6241898207, 75058.8771715961, 79487.8877884637},
            {176257.349845383, 179794.807761573, 180290.319046323, 181487.967471084, 186690.959612256},
            {286286.175818458, 303984.726428782, 321692.362821643, 335551.688987588, 344087.839487745},
            {399372.560529476, 433400.354292387, 471964.89621373, 513835.931064411, 555824.663124966},
            {577604.594975221, 635258.237156301, 698999.445970987, 768745.631252166, std::nan("")},
        };
        const double expected_y[isoline_count][points_per_isoline] = {
            {25000, 109298.142136262, 477843.354977538, 2089095.63724813, 9133370.87584761},
            {25000, 109298.142136262, 477843.354977538, 2089095.63724813, 9133370.87584761},
            {25000, 109298.142136262, 477843.354977538, 2089095.63724813, 9133370.87584761},
            {25000, 109298.142136262, 477843.354977538, 2089095.63724813, 9133370.87584761},
            {25000, 109298.142136262, 477843.354977538, 2089095.63724813, 9133370.87584761},
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
        CHECK_THAT(d_isolines[1].value, WithinAbs(4.704765645219758, 1e-10));
        CHECK_THAT(d_isolines[2].value, WithinAbs(32.79339504847662, 1e-10));
        CHECK_THAT(d_isolines[3].value, WithinAbs(228.57817793711163, 1e-10));
        CHECK_THAT(d_isolines[4].value, WithinAbs(1593.2471569904417, 1e-10));
        const double expected_x[isoline_count][points_per_isoline] = {
            {577604.594975212, std::nan(""), std::nan(""), std::nan(""), std::nan("")},
            {202365.843978511, 419230.112111493, std::nan(""), std::nan(""), std::nan("")},
            {142114.491283644, 204388.004478758, 351216.809707051, std::nan(""), std::nan("")},
            {133470.418481246, 172415.768780675, 235383.044874193, 357492.457483747, 669493.625997729},
            {70518.3287895177, 70601.2088976224, 70963.5807789929, 72548.359197014, 79487.8877879113},
        };
        const double expected_y[isoline_count][points_per_isoline] = {
            {25000, 109298.142136262, 477843.354977538, 2089095.63724813, 9133370.87584761},
            {25000, 109298.142136262, 477843.354977538, 2089095.63724813, 9133370.87584761},
            {25000, 109298.142136262, 477843.354977538, 2089095.63724813, 9133370.87584761},
            {25000, 109298.142136262, 477843.354977538, 2089095.63724813, 9133370.87584761},
            {25000, 109298.142136262, 477843.354977538, 2089095.63724813, 9133370.87584761},
        };
        for (int i = 0; i < d_isolines.size(); ++i) {
            REQUIRE(d_isolines[i].size() == points_per_isoline);
            for (int j = 0; j < d_isolines[i].size(); ++j) {
                if (std::isnan(d_isolines[i].x[j]))
                    CHECK(std::isnan(expected_x[i][j]));
                else
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
    CHECK_THAT(plot.xaxis.min, WithinAbs(426.00948386039755, 1));
    CHECK_THAT(plot.xaxis.max, WithinAbs(2423.07169138421, 1));
    CHECK_THAT(plot.yaxis.min, WithinAbs(173.15, 1));
    CHECK_THAT(plot.yaxis.max, WithinAbs(455.0, 1));

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
            {412.617538232079, 728.71482941326, 994.524404955042, 1237.31924154895, 1561.70306865236},
            {800.440438274308, 992.708859865778, 1177.8221470675, 1354.80424987622, 1561.8974228315},
            {1188.26333831654, 1256.70289031829, 1361.11988917995, 1472.28925820349, 1562.09177701064},
            {1576.08623835876, 1520.69692077081, 1544.4176312924, 1589.77426653076, 1562.28613118978},
            {1963.90913840099, 1784.69095122333, 1727.71537340486, 1707.25927485803, 1562.48048536892},
        };
        const double expected_y[isoline_count][points_per_isoline] = {
            {169.850074842393, 220.940538422734, 272.031002003074, 323.121465583414, 374.211929163755},
            {169.850074842393, 220.940538422734, 272.031002003074, 323.121465583414, 374.211929163755},
            {169.850074842393, 220.940538422734, 272.031002003074, 323.121465583414, 374.211929163755},
            {169.850074842393, 220.940538422734, 272.031002003074, 323.121465583414, 374.211929163755},
            {169.850074842393, 220.940538422734, 272.031002003074, 323.121465583414, 374.211929163755},
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
        // P isolines
        CoolProp::Plot::Range p_range = plot.isoline_range(CoolProp::iP);
        std::vector<double> p_values = CoolProp::Plot::generate_values_in_range(CoolProp::iP, p_range, isoline_count);
        CoolProp::Plot::Isolines p_isolines = plot.calc_isolines(CoolProp::iP, p_values, points_per_isoline);
        REQUIRE(p_isolines.size() == isoline_count);
        CHECK_THAT(p_isolines[0].value, WithinAbs(25000.000000000007, 1e-7));
        CHECK_THAT(p_isolines[1].value, WithinAbs(109298.14213626183, 1e-7));
        CHECK_THAT(p_isolines[2].value, WithinAbs(477843.3549775384, 1e-7));
        CHECK_THAT(p_isolines[3].value, WithinAbs(2089095.6372481277, 1e-7));
        CHECK_THAT(p_isolines[4].value, WithinAbs(9133370.87584761, 1e-7));
        const double expected_x[isoline_count][points_per_isoline] = {
            {426.009483860398, 925.275035741351, 1424.5405876223, 1923.80613950326, 2423.07169138421},
            {426.009483860398, 925.275035741351, 1424.5405876223, 1923.80613950326, 2423.07169138421},
            {426.009483860398, 925.275035741351, 1424.5405876223, 1923.80613950326, 2423.07169138421},
            {426.009483860398, 925.275035741351, 1424.5405876223, 1923.80613950326, 2423.07169138421},
            {426.009483860398, 925.275035741351, 1424.5405876223, 1923.80613950326, 2423.07169138421},
        };
        const double expected_y[isoline_count][points_per_isoline] = {
            {171.786072659192, 220.381369310476, 220.381369310476, 265.362477224881, 455.000000000006},
            {171.798910666292, 248.745218249749, 248.745218249749, 308.633922577123, 506.387752763855},
            {171.854988815762, 258.195699077866, 287.473037337147, 355.964867192619, 560.29310120217},
            {172.099235421196, 258.742471436004, 342.561817261331, 411.323964493198, 618.036314177106},
            {173.15, 261.021061581425, 371.327173900344, 484.427831614361, std::nan("")},
        };
        for (int i = 0; i < p_isolines.size(); ++i) {
            REQUIRE(p_isolines[i].size() == points_per_isoline);
            for (int j = 0; j < p_isolines[i].size(); ++j) {
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
            {426.009483860398, 925.275035741351, 1424.5405876223, 1923.80613950326, 2423.07169138421},
            {426.009483860398, 925.275035741351, 1424.5405876223, 1923.80613950326, 2423.07169138421},
            {426.009483860398, 925.275035741351, 1424.5405876223, 1923.80613950326, 2423.07169138421},
            {426.009483860398, 925.275035741351, 1424.5405876223, 1923.80613950326, 2423.07169138421},
            {426.009483860398, 925.275035741351, 1424.5405876223, 1923.80613950326, 2423.07169138421},
        };
        const double expected_y[isoline_count][points_per_isoline] = {
            {172.174575309065, std::nan(""), std::nan(""), std::nan(""), std::nan("")},
            {196.074550634008, 266.631159312075, std::nan(""), std::nan(""), std::nan("")},
            {213.664681842583, 299.984652703232, 301.726570477946, std::nan(""), std::nan("")},
            {228.411201679534, 322.843563825212, 426.787882130168, 331.521169967777, 328.042167528594},
            {241.568258023047, 341.661338916035, 458.593848045394, std::nan(""), 455.000000000079},
        };
        for (int i = 0; i < h_isolines.size(); ++i) {
            REQUIRE(h_isolines[i].size() == points_per_isoline);
            for (int j = 0; j < h_isolines[i].size(); ++j) {
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
        CHECK_THAT(d_isolines[1].value, WithinAbs(4.704765645219758, 1e-10));
        CHECK_THAT(d_isolines[2].value, WithinAbs(32.79339504847662, 1e-10));
        CHECK_THAT(d_isolines[3].value, WithinAbs(228.57817793711163, 1e-10));
        CHECK_THAT(d_isolines[4].value, WithinAbs(1593.2471569904417, 1e-10));
        const double expected_x[isoline_count][points_per_isoline] = {
            {524.17387831234, 1911.09303197673, 2092.95299735844, 2262.71394473455, 2423.07169138421},
            {448.103089616845, 1715.11956249481, 1932.46627813427, 2103.15612327654, 2263.90953791772},
            {437.189451894057, 972.489749676211, 1758.36241052056, 1935.7522861596, 2099.20643194095},
            {435.623706482622, 865.946977105694, 1292.02339683139, 1720.27746043057, 1899.38158004697},
            {426.009483860398, 710.877062878169, 946.968704707899, 1151.91782375377, 1335.56507098504},
        };
        const double expected_y[isoline_count][points_per_isoline] = {
            {173.15, 243.6125, 314.075, 384.5375, 455},
            {173.15, 243.6125, 314.075, 384.5375, 455},
            {173.15, 243.6125, 314.075, 384.5375, 455},
            {173.15, 243.6125, 314.075, 384.5375, 455},
            {173.15, 243.6125, 314.075, 384.5375, 455},
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

// </autogenerated>

#endif /* ENABLE_CATCH */
