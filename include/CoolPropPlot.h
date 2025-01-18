#ifndef COOLPROPPLOT_H_
#define COOLPROPPLOT_H_

#include "AbstractState.h"
#include "CPnumerics.h"
#include <vector>
#include <map>

namespace CoolProp {
namespace Plot {

enum class Scale
{
    Lin,
    Log
};

struct Range
{
    double min, max;
};

namespace Detail {

const double NaN = std::numeric_limits<double>::quiet_NaN();

enum IsolineSupported
{
    No = 0,
    Yes = 1,
    Flipped = 2
};

 const int TS = CoolProp::iT * 10 + CoolProp::iSmass;
 const int PH = CoolProp::iP * 10 + CoolProp::iHmass;
 const int HS = CoolProp::iHmass * 10 + CoolProp::iSmass;
 const int PS = CoolProp::iP * 10 + CoolProp::iSmass;
 const int PD = CoolProp::iP * 10 + CoolProp::iDmass;
 const int TD = CoolProp::iT * 10 + CoolProp::iDmass;
 const int PT = CoolProp::iP * 10 + CoolProp::iT;
 const int PU = CoolProp::iP * 10 + CoolProp::iUmass;

 const std::map<CoolProp::parameters, std::map<int, IsolineSupported>> xy_switch = {
    {CoolProp::iDmass, {{TS, Flipped}, {PH, Flipped}, {HS, Yes    }, {PS, Flipped}, {PD, No     }, {TD, No     }, {PT, Yes    }}},
    {CoolProp::iHmass, {{TS, Yes    }, {PH, No     }, {HS, No     }, {PS, Flipped}, {PD, Flipped}, {TD, Yes    }, {PT, Yes    }}},
    {CoolProp::iP,     {{TS, Yes    }, {PH, No     }, {HS, Yes    }, {PS, No     }, {PD, No     }, {TD, Yes    }, {PT, No     }}},
    {CoolProp::iSmass, {{TS, No     }, {PH, Flipped}, {HS, No     }, {PS, No     }, {PD, Flipped}, {TD, Yes    }, {PT, Flipped}}},
    {CoolProp::iT,     {{TS, No     }, {PH, Flipped}, {HS, Yes    }, {PS, Yes    }, {PD, Yes    }, {TD, No     }, {PT, No     }}},
    {CoolProp::iQ,     {{TS, Flipped}, {PH, Flipped}, {HS, Flipped}, {PS, Flipped}, {PD, Flipped}, {TD, Flipped}, {PT, Yes    }}}
};

Scale default_scale(CoolProp::parameters key);
std::shared_ptr<CoolProp::AbstractState> process_fluid_state(const std::string& fluid_ref);
inline std::shared_ptr<CoolProp::AbstractState> get_critical_point(const std::shared_ptr<CoolProp::AbstractState>& state);

} /* namespace Detail */

class Isoline
{
public:
    std::vector<double> x;
    std::vector<double> y;
    double value;

private:
    std::shared_ptr<CoolProp::AbstractState> state;
    std::shared_ptr<CoolProp::AbstractState> critical_state;
    CoolProp::parameters xkey;
    CoolProp::parameters ykey;
    CoolProp::parameters key;

    Isoline(CoolProp::parameters key, CoolProp::parameters xkey, CoolProp::parameters ykey, double value, const std::shared_ptr<CoolProp::AbstractState>& state);

    Range get_sat_bounds(CoolProp::parameters key);
    void calc_sat_range(int count);
    void update_pair(int& ipos, int& xpos, int& ypos, int& pair);
    void calc_range(std::vector<double>& xvals, std::vector<double>& yvals);

    friend class PropertyPlot;
};

using Isolines = std::vector<Isoline>;

std::vector<double> generate_values_in_range(Scale scale, const Range& range, int count);
std::vector<double> generate_values_in_range(CoolProp::parameters type, const Range& range, int count);

class PropertyPlot
{
public:
    CoolProp::parameters xkey;
    CoolProp::parameters ykey;
    Scale xscale;
    Scale yscale;
    Range xrange;
    Range yrange;

    PropertyPlot(const std::string& fluid_name, CoolProp::parameters ykey, CoolProp::parameters xkey, const std::string& tp_limits);

    Range isoline_range(CoolProp::parameters key);
    Isolines calc_isolines(CoolProp::parameters key, const std::vector<double>& values, int points) const;
    std::vector<CoolProp::parameters> supported_isoline_keys() const;
    double value_at(CoolProp::parameters key, double xvalue, double yvalue, CoolProp::phases phase = CoolProp::phases::iphase_not_imposed) const;

private:
    struct Range2D
    {
        union
        {
            Range x, T;
        };
        union
        {
            Range y, p;
        };
    };

    std::string fluid_name;
    CoolProp::input_pairs axis_pair;
    bool swap_axis_inputs_for_update;
    std::shared_ptr<CoolProp::AbstractState> state;
    std::shared_ptr<CoolProp::AbstractState> critical_state;
    std::vector<double> limits;

    Range get_sat_bounds(CoolProp::parameters key);
    void get_Tp_limits(double& T_lo, double& T_hi, double& P_lo, double& P_hi);
    Range2D get_axis_limits(CoolProp::parameters xkey = CoolProp::parameters::iundefined_parameter, CoolProp::parameters ykey = CoolProp::parameters::iundefined_parameter, bool autoscale = true);
};

} /* namespace Plot */
} /* namespace CoolProp */

#endif /* COOLPROPPLOT_H_ */





