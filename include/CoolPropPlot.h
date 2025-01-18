#ifndef COOLPROPPLOT_H_
#define COOLPROPPLOT_H_

#include "CoolProp.h"
#include "AbstractState.h"
#include "CPnumerics.h"
#include <vector>
#include <map>

namespace CoolProp {
namespace Plot {

namespace Detail {

const double NaN = std::numeric_limits<double>::quiet_NaN();

enum IsolineSupported
{
    No = 0,
    Yes = 1,
    Flipped = 2
};

static const int TS = CoolProp::iT * 10 + CoolProp::iSmass;
static const int PH = CoolProp::iP * 10 + CoolProp::iHmass;
static const int HS = CoolProp::iHmass * 10 + CoolProp::iSmass;
static const int PS = CoolProp::iP * 10 + CoolProp::iSmass;
static const int PD = CoolProp::iP * 10 + CoolProp::iDmass;
static const int TD = CoolProp::iT * 10 + CoolProp::iDmass;
static const int PT = CoolProp::iP * 10 + CoolProp::iT;
static const int PU = CoolProp::iP * 10 + CoolProp::iUmass;

static std::map<CoolProp::parameters, std::map<int, IsolineSupported>> xy_switch = {
    {CoolProp::iDmass, {{TS, Flipped}, {PH, Flipped}, {HS, Yes    }, {PS, Flipped}, {PD, No     }, {TD, No     }, {PT, Yes    }}},
    {CoolProp::iHmass, {{TS, Yes    }, {PH, No     }, {HS, No     }, {PS, Flipped}, {PD, Flipped}, {TD, Yes    }, {PT, Yes    }}},
    {CoolProp::iP,     {{TS, Yes    }, {PH, No     }, {HS, Yes    }, {PS, No     }, {PD, No     }, {TD, Yes    }, {PT, No     }}},
    {CoolProp::iSmass, {{TS, No     }, {PH, Flipped}, {HS, No     }, {PS, No     }, {PD, Flipped}, {TD, Yes    }, {PT, Flipped}}},
    {CoolProp::iT,     {{TS, No     }, {PH, Flipped}, {HS, Yes    }, {PS, Yes    }, {PD, Yes    }, {TD, No     }, {PT, No     }}},
    {CoolProp::iQ,     {{TS, Flipped}, {PH, Flipped}, {HS, Flipped}, {PS, Flipped}, {PD, Flipped}, {TD, Flipped}, {PT, Yes    }}}
};

}

enum class Scale
{
    Lin,
    Log
};

namespace Detail {

inline Scale default_scale(CoolProp::parameters key)
{
    switch (key)
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

inline std::shared_ptr<CoolProp::AbstractState> process_fluid_state(const std::string& fluid_ref)
{
    std::string backend;
    std::string fluids;
    CoolProp::extract_backend(fluid_ref, backend, fluids);
    std::vector<double> fractions;
    fluids = CoolProp::extract_fractions(fluids, fractions);

    return std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory(backend, fluids));
}

inline std::vector<double> generate_values_in_range(Scale scale, double start, double end, int count)
{
    if (scale == Scale::Log)
        return logspace(start, end, count);
    else
        return linspace(start, end, count);
}

inline std::vector<double> generate_values_in_range(CoolProp::parameters type, double start, double end, int count)
{
    return generate_values_in_range(default_scale(type), start, end, count);
}

inline std::shared_ptr<CoolProp::AbstractState> get_critical_point(const std::shared_ptr<CoolProp::AbstractState>& state)
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


class IsoLine
{
public:
    std::vector<double> x;
    std::vector<double> y;
    double value;

    IsoLine(CoolProp::parameters key, CoolProp::parameters xkey, CoolProp::parameters ykey, double value, const std::shared_ptr<CoolProp::AbstractState>& state)
        : key(key),
          xkey(xkey),
          ykey(ykey),
          value(value),
          state(state)
    {
        this->critical_state = Detail::get_critical_point(state);
    }
private:
    std::shared_ptr<CoolProp::AbstractState> state;
    std::shared_ptr<CoolProp::AbstractState> critical_state;
    CoolProp::parameters xkey;
    CoolProp::parameters ykey;
    CoolProp::parameters key;

    Range get_sat_bounds(CoolProp::parameters key)
    {
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
        double sat_min = fluid_min;
        double sat_max = fluid_max;
        return {sat_min, sat_max};
    }
    void calc_sat_range(int count)
    {
        double t_lo, t_hi;
        Range t = get_sat_bounds(CoolProp::iT);
        t_lo = t.min;
        t_hi = t.max;
        std::vector<double> two = ::linspace(t_lo, t_hi, count);
        std::vector<double> one(two.size(), value);
        CoolProp::input_pairs input_pair = CoolProp::QT_INPUTS;

        double t_crit = critical_state->keyed_output(CoolProp::iT);
        double p_crit = critical_state->keyed_output(CoolProp::iP);
        double x_crit = critical_state->keyed_output(xkey);
        double y_crit = critical_state->keyed_output(ykey);
        x.resize(one.size());
        y.resize(one.size());
        for (int i = 0; i < one.size(); ++i)
        {
            try
            {
                state->update(input_pair, one[i], two[i]);
                x[i] = state->keyed_output(xkey);
                y[i] = state->keyed_output(ykey);
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
        Detail::IsolineSupported should_switch = Detail::xy_switch.at(key).at(ykey * 10 + xkey);
        double out1, out2;
        switch (should_switch)
        {
            case Detail::IsolineSupported::No:
                throw CoolProp::ValueError("This isoline cannot be calculated!");
                break;
            case Detail::IsolineSupported::Yes:
                pair = CoolProp::generate_update_pair(key, 0.0, xkey, 1.0, out1, out2);
                break;
            case Detail::IsolineSupported::Flipped:
                pair = CoolProp::generate_update_pair(key, 0.0, ykey, 1.0, out1, out2);
                break;
        }
        bool should_swap = (out1 != 0.0);

        if (should_switch == Detail::IsolineSupported::Yes && !should_swap)
        {
            ipos = 0;
            xpos = 1;
            ypos = 2;
        }
        else if (should_switch == Detail::IsolineSupported::Flipped && !should_swap)
        {
            ipos = 0;
            xpos = 2;
            ypos = 1;
        }
        else if (should_switch == Detail::IsolineSupported::Yes && should_swap)
        {
            ipos = 1;
            xpos = 0;
            ypos = 2;
        }
        else if (should_switch == Detail::IsolineSupported::Flipped && should_swap)
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
        if (key == CoolProp::iQ)
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
            idxs[ipos] = key;
            idxs[xpos] = xkey;
            idxs[ypos] = ykey;
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
                if (idxs[i] == xkey) x = vals[i];
                if (idxs[i] == ykey) y = vals[i];
            }
        }
    }

    friend class PropertyPlot;
};

using IsoLines = std::vector<IsoLine>;

class PropertyPlot
{
public:
    CoolProp::parameters xkey;
    CoolProp::parameters ykey;
    Scale axis_x_scale;
    Scale axis_y_scale;
    Range axis_x_range;
    Range axis_y_range;

    PropertyPlot(const std::string& fluid_name, CoolProp::parameters ykey, CoolProp::parameters xkey, const std::string& tp_limits);
    Range isoline_range(CoolProp::parameters key);
    IsoLines calc_isolines(CoolProp::parameters key, const std::vector<double>& values, int points) const;
    std::vector<CoolProp::parameters> supported_dimensions() const;
    double value_at(CoolProp::parameters key, double axis_x_value, double axis_y_value, CoolProp::phases phase = CoolProp::phases::iphase_not_imposed) const;

private:
    std::string fluid_name;
    CoolProp::input_pairs axis_pair;
    bool swap_axis_inputs_for_update;
    std::shared_ptr<CoolProp::AbstractState> state;
    std::shared_ptr<CoolProp::AbstractState> critical_state;
    std::vector<double> limits;

    Range get_sat_bounds(CoolProp::parameters key);
    void get_Tp_limits(double& T_lo, double& T_hi, double& P_lo, double& P_hi);
    std::vector<double> get_axis_limits(CoolProp::parameters xkey = CoolProp::parameters::iundefined_parameter, CoolProp::parameters ykey = CoolProp::parameters::iundefined_parameter, bool autoscale = true);
};

} /* namespace Plot */
} /* namespace CoolProp */

#endif /* COOLPROPPLOT_H_ */





