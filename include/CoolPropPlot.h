#ifndef COOLPROPPLOT_H_
#define COOLPROPPLOT_H_

#include "CoolProp.h"
#include "AbstractState.h"
#include "CPnumerics.h"
#include <vector>
#include <memory>
#include <unordered_map>

namespace CoolProp {
namespace Plot {

namespace Detail {

template <typename T>
struct Optional {
    Optional(T value) : value(value), has_value_(true) {}
    Optional() : has_value_(false) {}
    T operator*() const { return value; }
    operator bool() const { return has_value_; }
    bool has_value() const { return has_value_; }

    T value;
    bool has_value_;
};

const double NaN = std::numeric_limits<double>::quiet_NaN();

}

enum class Scale
{
    Lin,
    Log
};

inline Scale default_scale(CoolProp::parameters iso_type)
{
    switch (iso_type)
    {
        case CoolProp::iDmass: return Scale::Log;
        case CoolProp::iHmass: return Scale::Lin;
        case CoolProp::iP:     return Scale::Log;
        case CoolProp::iSmass: return Scale::Lin;
        case CoolProp::iT:     return Scale::Lin;
        case CoolProp::iUmass: return Scale::Lin;
        case CoolProp::iQ:     return Scale::Lin;

        default: return Scale::Lin;
    }
}


namespace Detail
{

inline std::shared_ptr<CoolProp::AbstractState> process_fluid_state(std::string fluid_ref)
{
    std::string backend;
    std::string fluids;
    CoolProp::extract_backend(fluid_ref, backend, fluids);
    std::vector<double> fractions;
    fluids = CoolProp::extract_fractions(fluids, fractions);

    return std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory(backend, fluids));
}

inline std::vector<double> generate_ranges(Scale scale, double start, double end, int num)
{
    if (scale == Scale::Log)
        return logspace(start, end, num);
    else
        return linspace(start, end, num);
}

inline std::vector<double> generate_ranges(CoolProp::parameters type, double start, double end, int num)
{
    return generate_ranges(default_scale(type), start, end, num);
}

inline std::shared_ptr<CoolProp::AbstractState> get_critical_point(std::shared_ptr<CoolProp::AbstractState> state)
{
    CoolProp::CriticalState crit_state;
    crit_state.T = Detail::NaN;
    crit_state.p = Detail::NaN;
    crit_state.rhomolar = Detail::NaN;
    crit_state.rhomolar = Detail::NaN;
    crit_state.stable = false;
    try
    {
        crit_state.T = state->T_critical();
        crit_state.p = state->p_critical();
        crit_state.rhomolar = state->rhomolar_critical();
        crit_state.stable = true;
    }
    catch (...)
    {
        try
        {
            for (CoolProp::CriticalState crit_state_tmp: state->all_critical_points())
            {
                if (crit_state_tmp.stable && (crit_state_tmp.T > crit_state.T || !std::isfinite(crit_state.T)))
                {
                    crit_state.T = crit_state_tmp.T;
                    crit_state.p = crit_state_tmp.p;
                    crit_state.rhomolar = crit_state_tmp.rhomolar;
                    crit_state.stable = crit_state_tmp.stable;
                }
            }
        }
        catch (...)
        {
            throw CoolProp::ValueError("Could not calculate the critical point data.");
        }
    }

    std::shared_ptr<CoolProp::AbstractState> new_state(CoolProp::AbstractState::factory(state->backend_name(), state->fluid_names()));
    std::vector<double> masses = state->get_mass_fractions();
    if (masses.size() > 1)
        new_state->set_mass_fractions(masses);

    if (std::isfinite(crit_state.p) && std::isfinite(crit_state.T))
    {
        try
        {
            new_state->specify_phase(CoolProp::iphase_critical_point);
            new_state->update(CoolProp::PT_INPUTS, crit_state.p, crit_state.T);
            return new_state;
        }
        catch (...) { }
        try
        {
            new_state->update(CoolProp::PT_INPUTS, crit_state.p, crit_state.T);
            return new_state;
        }
        catch (...) { }
    }

    if (std::isfinite(crit_state.rhomolar) && std::isfinite(crit_state.T))
    {
        try
        {
            new_state->specify_phase(CoolProp::iphase_critical_point);
            new_state->update(CoolProp::DmolarT_INPUTS, crit_state.rhomolar, crit_state.T);
            return new_state;
        }
        catch (...) { }
        try
        {
            new_state->update(CoolProp::DmolarT_INPUTS, crit_state.rhomolar, crit_state.T);
            return new_state;
        }
        catch (...) { }
    }
    throw CoolProp::ValueError("Could not calculate the critical point data.");
    return nullptr;
}

} // namespace Detail

struct Range
{
    double min, max;
};


static std::unordered_map<CoolProp::parameters, std::unordered_map<int, Detail::Optional<bool>>> xy_switch;

class IsoLine
{
public:
    std::vector<double> x;
    std::vector<double> y;
    double value;

    IsoLine(CoolProp::parameters i_index, CoolProp::parameters x_index, CoolProp::parameters y_index, double value, std::shared_ptr<CoolProp::AbstractState> state)
    {
        this->state = state;
        this->critical_state = Detail::get_critical_point(state);
        this->i_index = i_index;
        this->x_index = x_index;
        this->y_index = y_index;
        this->value = value;

        fill_xy_switch(); // TODO: remove
    }
    size_t size() const
    {
        return x.size();
    }
private:
    void fill_xy_switch()
    {
        static const int TS = CoolProp::iT * 10 + CoolProp::iSmass;
        static const int PH = CoolProp::iP * 10 + CoolProp::iHmass;
        static const int HS = CoolProp::iHmass * 10 + CoolProp::iSmass;
        static const int PS = CoolProp::iP * 10 + CoolProp::iSmass;
        static const int PD = CoolProp::iP * 10 + CoolProp::iDmass;
        static const int TD = CoolProp::iT * 10 + CoolProp::iDmass;
        static const int PT = CoolProp::iP * 10 + CoolProp::iT;
        static const int PU = CoolProp::iP * 10 + CoolProp::iUmass;

        xy_switch = {
            {CoolProp::iDmass, {{TS, true }, {PH, true}, {HS, false}, {PS, true }, {PD, {}   }, {TD, {}   }, {PT, false}}},
            {CoolProp::iHmass, {{TS, false}, {PH, {}  }, {HS, {}   }, {PS, true }, {PD, true }, {TD, false}, {PT, false}}},
            {CoolProp::iP,     {{TS, false}, {PH, {}  }, {HS, false}, {PS, {}   }, {PD, {}   }, {TD, false}, {PT, {}  }}} ,
            {CoolProp::iSmass, {{TS, {}   }, {PH, true}, {HS, {}   }, {PS, {}   }, {PD, true }, {TD, false}, {PT, true}}} ,
            {CoolProp::iT,     {{TS, {}   }, {PH, true}, {HS, false}, {PS, false}, {PD, false}, {TD, {}   }, {PT, {}  }}} ,
            {CoolProp::iQ,     {{TS, true }, {PH, true}, {HS, true }, {PS, true }, {PD, true }, {TD, true }, {PT, false}}}
        };
    }

    std::shared_ptr<CoolProp::AbstractState> state;
    std::shared_ptr<CoolProp::AbstractState> critical_state;
    CoolProp::parameters x_index;
    CoolProp::parameters y_index;
    CoolProp::parameters i_index;

    Range get_sat_bounds(CoolProp::parameters kind)
    {
        double s = 1e-7;
        double t_small = critical_state->keyed_output(CoolProp::iT) * s;
        double p_small = critical_state->keyed_output(CoolProp::iP) * s;

        double t_triple = state->trivial_keyed_output(CoolProp::iT_triple);
        double t_min = state->trivial_keyed_output(CoolProp::iT_min);
        state->update(CoolProp::QT_INPUTS, 0, std::max(t_triple, t_min) + t_small);
        double fluid_min, fluid_max;
        if (kind == CoolProp::iP)
        {
            fluid_min = state->keyed_output(CoolProp::iP) + p_small;
            fluid_max = critical_state->keyed_output(CoolProp::iP) - p_small;
        }
        else if (kind == CoolProp::iT)
        {
            fluid_min = state->keyed_output(CoolProp::iT) + t_small;
            fluid_max = critical_state->keyed_output(CoolProp::iT) - t_small;
        }
        else
        {
            throw CoolProp::ValueError("Invalid kind");
        }
        double sat_min = fluid_min;
        double sat_max = fluid_max;
        return {sat_min, sat_max};
    }
    void calc_sat_range(int num)
    {
        double t_lo, t_hi;
        auto t = get_sat_bounds(CoolProp::iT);
        t_lo = t.min;
        t_hi = t.max;
        std::vector<double> two = ::linspace(t_lo, t_hi, num);
        std::vector<double> one(two.size(), value);
        CoolProp::input_pairs input_pair = CoolProp::QT_INPUTS;

        double t_crit = critical_state->keyed_output(CoolProp::iT);
        double p_crit = critical_state->keyed_output(CoolProp::iP);
        double x_crit = critical_state->keyed_output(x_index);
        double y_crit = critical_state->keyed_output(y_index);
        x.resize(one.size());
        y.resize(one.size());
        for (int i = 0; i < one.size(); ++i)
        {
            try
            {
                state->update(input_pair, one[i], two[i]);
                x[i] = state->keyed_output(x_index);
                y[i] = state->keyed_output(y_index);
            }
            catch (...)
            {
                if (input_pair == CoolProp::QT_INPUTS && abs(two[i] - t_crit) < 1e0 ||
                    input_pair == CoolProp::PQ_INPUTS && abs(one[i] - p_crit) < 1e2)
                {
                    x[i] = x_crit;
                    y[i] = y_crit;
                    std::cerr << "ERROR near critical inputs" << std::endl;
                }
                else
                {
                    x[i] = Detail::NaN;
                    y[i] = Detail::NaN;
                    std::cerr << "ERROR" << std::endl;
                }
            }
        }
    }

    void update_pair(int& ipos, int& xpos, int& ypos, int& pair)
    {
        Detail::Optional<bool> should_switch = xy_switch.at(i_index).at(y_index * 10 + x_index);
        double out1, out2;
        if (!should_switch)
            throw CoolProp::ValueError("This isoline cannot be calculated!");
        else if (*should_switch == false)
            pair = CoolProp::generate_update_pair(i_index, 0.0, x_index, 1.0, out1, out2);
        else if (*should_switch == true)
            pair = CoolProp::generate_update_pair(i_index, 0.0, y_index, 1.0, out1, out2);

        bool should_swap = (out1 != 0.0);
        if (!(*should_switch) && !should_swap)
        {
            ipos = 0;
            xpos = 1;
            ypos = 2;
        }
        else if (*should_switch && !should_swap)
        {
            ipos = 0;
            xpos = 2;
            ypos = 1;
        }
        else if (!(*should_switch) && should_swap)
        {
            ipos = 1;
            xpos = 0;
            ypos = 2;
        }
        else if (*should_switch && should_swap)
        {
            ipos = 1;
            xpos = 2;
            ypos = 0;
        }
        else
        {
            throw CoolProp::ValueError("Check the code, this should not happen!");
        }
    }

    void calc_range(std::vector<double>& xvals, std::vector<double>& yvals)
    {
        if (i_index == CoolProp::iQ)
        {
            calc_sat_range(xvals.size());
        }
        else
        {
            int ipos, xpos, ypos, pair;
            update_pair(ipos, xpos, ypos, pair);

            std::vector<double> ivals(xvals.size(), value);
            std::vector<int> order = {ipos, xpos, ypos};
            std::vector<CoolProp::parameters> idxs(3);
            idxs[ipos] = i_index;
            idxs[xpos] = x_index;
            idxs[ypos] = y_index;
            std::vector<std::vector<double>> vals(3);
            vals[ipos] = ivals;
            vals[xpos] = xvals;
            vals[ypos] = yvals;

            // TODO: guesses missing

            for (int i = 0; i < vals[2].size(); ++i)
            {
                try
                {
                    state->update((CoolProp::input_pairs)pair, vals[0][i], vals[1][i]);
                    vals[2][i] = state->keyed_output(idxs[2]);
                }
                catch (...)
                {
                    vals[2][i] = Detail::NaN;
                }
            }

            for (int i = 0; i < idxs.size(); ++i)
            {
                if (idxs[i] == x_index) x = vals[i];
                if (idxs[i] == y_index) y = vals[i];
            }
        }
    }

    friend class PropertyPlot;
};

using IsoLines = std::vector<IsoLine>;

class PropertyPlot
{
public:
    CoolProp::parameters x_index;
    CoolProp::parameters y_index;

    PropertyPlot(std::string fluid_name, CoolProp::parameters y_index, CoolProp::parameters x_index, std::string tp_limits)
    {
        this->fluid_name = fluid_name;
        this->state = Detail::process_fluid_state(fluid_name);
        this->critical_state = Detail::get_critical_point(state);
        this->x_index = x_index;
        this->y_index = y_index;
        this->axis_x_scale_ = default_scale(x_index);
        this->axis_y_scale_ = default_scale(y_index);
        // We are just assuming that all inputs and outputs are in SI units. We
        // take care of any conversions before calling the library and after
        // getting the results.
        int out1, out2;
        axis_pair = CoolProp::generate_update_pair(x_index, 0, y_index, 1, out1, out2);
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

    Range isoline_range(CoolProp::parameters iso_index)
    {
        std::vector<double> iso_range;
        if (iso_index == CoolProp::iQ)
            return {0., 1.};
        else // TODO: always against iT?
        {
            std::vector<double> range = get_axis_limits(iso_index, CoolProp::iT);
            return {range[0], range[1]};
        }
    }

    IsoLines calc_isolines(CoolProp::parameters iso_index, const std::vector<double>& iso_values, int points)
    {
        std::vector<double> ixrange = Detail::generate_ranges(axis_x_scale_, axis_x_limits.min, axis_x_limits.max, points);
        std::vector<double> iyrange = Detail::generate_ranges(axis_y_scale_, axis_y_limits.min, axis_y_limits.max, points);

        IsoLines lines;
        for (double iso_value : iso_values)
        {
            IsoLine line(iso_index, x_index, y_index, iso_value, state);
            line.calc_range(ixrange, iyrange);
            // TODO: line.sanitize_data();
            lines.push_back(line);
        }
        return lines;
    }

    std::vector<CoolProp::parameters> supported_dimensions() const
    {
        // static const int TS = CoolProp::iT * 10 + CoolProp::iSmass;
        // static const int PH = CoolProp::iP * 10 + CoolProp::iHmass;
        // static const int HS = CoolProp::iHmass * 10 + CoolProp::iSmass;
        // static const int PS = CoolProp::iP * 10 + CoolProp::iSmass;
        // static const int PD = CoolProp::iP * 10 + CoolProp::iDmass;
        // static const int TD = CoolProp::iT * 10 + CoolProp::iDmass;
        // static const int PT = CoolProp::iP * 10 + CoolProp::iT;
        // static const int PU = CoolProp::iP * 10 + CoolProp::iUmass;
        //
        // xy_switch = {
        //     {CoolProp::iDmass, {{TS, true }, {PH, true}, {HS, false}, {PS, true }, {PD, {}   }, {TD, {}   }, {PT, false}}},
        //     {CoolProp::iHmass, {{TS, false}, {PH, {}  }, {HS, {}   }, {PS, true }, {PD, true }, {TD, false}, {PT, false}}},
        //     {CoolProp::iP,     {{TS, false}, {PH, {}  }, {HS, false}, {PS, {}   }, {PD, {}   }, {TD, false}, {PT, {}  }}} ,
        //     {CoolProp::iSmass, {{TS, {}   }, {PH, true}, {HS, {}   }, {PS, {}   }, {PD, true }, {TD, false}, {PT, true}}} ,
        //     {CoolProp::iT,     {{TS, {}   }, {PH, true}, {HS, false}, {PS, false}, {PD, false}, {TD, {}   }, {PT, {}  }}} ,
        //     {CoolProp::iQ,     {{TS, true }, {PH, true}, {HS, true }, {PS, true }, {PD, true }, {TD, true }, {PT, false}}}
        // };

        if (x_index == CoolProp::iHmass && y_index == CoolProp::iP) return {CoolProp::iQ, CoolProp::iT, CoolProp::iSmass, CoolProp::iDmass};
        if (x_index == CoolProp::iP && y_index == CoolProp::iHmass) return {CoolProp::iQ, CoolProp::iT, CoolProp::iSmass, CoolProp::iDmass};

        if (x_index == CoolProp::iT && y_index == CoolProp::iSmass) return {CoolProp::iQ, CoolProp::iP, CoolProp::iHmass, CoolProp::iDmass};
        if (x_index == CoolProp::iSmass && y_index == CoolProp::iT) return {CoolProp::iQ, CoolProp::iP, CoolProp::iHmass, CoolProp::iDmass};

        return {};
    }

    void set_axis_y_scale(Scale scale)
    {
        axis_y_scale_ = scale;
    }
    Scale axis_y_scale() const
    {
        return axis_y_scale_;
    }

    void set_axis_x_scale(Scale scale)
    {
        axis_x_scale_ = scale;
    }
    Scale axis_x_scale() const
    {
        return axis_x_scale_;
    }

    void set_axis_x_range(Range limits)
    {
        axis_x_limits = limits;
    }
    void set_axis_y_range(Range limits)
    {
        axis_y_limits = limits;
    }

    Range axis_x_range() const
    {
        return axis_x_limits; // TODO: just remove this function
    }
    Range axis_y_range() const
    {
        return axis_y_limits;
    }

    // for value under cursor
    double value_at(CoolProp::parameters iso_type, double axis_x_value, double axis_y_value, CoolProp::phases phase = CoolProp::phases::iphase_not_imposed)
    {
        if (iso_type == x_index) return axis_x_value;
        if (iso_type == y_index) return axis_y_value;

        try
        {
            if (swap_axis_inputs_for_update)
                std::swap(axis_x_value, axis_y_value);
            state->specify_phase(phase);
            state->update(axis_pair, axis_x_value, axis_y_value);
            switch (iso_type)
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

private:
    std::string fluid_name;
    std::string graph_type;
    CoolProp::input_pairs axis_pair;
    bool swap_axis_inputs_for_update;
    std::shared_ptr<CoolProp::AbstractState> state;
    std::shared_ptr<CoolProp::AbstractState> critical_state;
    std::vector<double> limits;
    Range axis_x_limits;
    Range axis_y_limits;
    Scale axis_x_scale_;
    Scale axis_y_scale_;

    Range get_sat_bounds(CoolProp::parameters kind)
    {
        // TODO: duplicated code from IsoLine
        double s = 1e-7;
        double t_small = critical_state->keyed_output(CoolProp::iT) * s;
        double p_small = critical_state->keyed_output(CoolProp::iP) * s;

        double t_triple = state->trivial_keyed_output(CoolProp::iT_triple);
        double t_min = state->trivial_keyed_output(CoolProp::iT_min);
        state->update(CoolProp::QT_INPUTS, 0, std::max(t_triple, t_min) + t_small);
        double fluid_min, fluid_max;
        if (kind == CoolProp::iP)
        {
            fluid_min = state->keyed_output(CoolProp::iP) + p_small;
            fluid_max = critical_state->keyed_output(CoolProp::iP) - p_small;
        }
        else if (kind == CoolProp::iT)
        {
            fluid_min = state->keyed_output(CoolProp::iT) + t_small;
            fluid_max = critical_state->keyed_output(CoolProp::iT) - t_small;
        }
        else
        {
            throw CoolProp::ValueError("Invalid kind");
        }
        double sat_min = fluid_min;
        double sat_max = fluid_max;
        return {sat_min, sat_max};
    }

    void get_Tp_limits(double& T_lo, double& T_hi, double& P_lo, double& P_hi)
    {
        T_lo = limits[0];
        T_hi = limits[1];
        P_lo = limits[2];
        P_hi = limits[3];

        double Ts_lo, Ts_hi;
        auto Ts = get_sat_bounds(CoolProp::iT);
        Ts_lo = Ts.min;
        Ts_hi = Ts.max;

        double Ps_lo, Ps_hi;
        auto Ps = get_sat_bounds(CoolProp::iP);
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

    std::vector<double> get_axis_limits(CoolProp::parameters x_index = CoolProp::parameters::iundefined_parameter, CoolProp::parameters y_index = CoolProp::parameters::iundefined_parameter, bool autoscale = true)
    {
        if (x_index == CoolProp::parameters::iundefined_parameter) x_index = this->x_index;
        if (y_index == CoolProp::parameters::iundefined_parameter) y_index = this->y_index;

        if (x_index != this->y_index || y_index != this->y_index || autoscale)
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
                        double x = state->keyed_output(x_index);
                        double y = state->keyed_output(y_index);
                        if (x < limits[0]) limits[0] = x;
                        if (x > limits[1]) limits[1] = x;
                        if (y < limits[2]) limits[2] = y;
                        if (y > limits[3]) limits[3] = y;
                    }
                    catch (...) { }
                }
            }
            if (x_index == this->x_index)
            {
                axis_x_limits.min = limits[0];
                axis_x_limits.max = limits[1];
            }
            if (y_index == this->y_index)
            {
                axis_y_limits.min = limits[2];
                axis_y_limits.max = limits[3];
            }
            return limits;
        }
        else
        {
            return {axis_x_limits.min, axis_x_limits.max, axis_y_limits.min, axis_y_limits.max};
        }
    }
};

} /* namespace Plot */
} /* namespace CoolProp */

#endif /* COOLPROPPLOT_H_ */





