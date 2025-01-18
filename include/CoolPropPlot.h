#ifndef COOLPROPPLOT_H_
#define COOLPROPPLOT_H_

#include "AbstractState.h"
#include <vector>

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

    size_t size() const { return x.size(); };

   private:
    std::shared_ptr<CoolProp::AbstractState> state_;
    std::shared_ptr<CoolProp::AbstractState> critical_state_;
    CoolProp::parameters xkey_;
    CoolProp::parameters ykey_;
    CoolProp::parameters key_;

    Isoline(CoolProp::parameters key, CoolProp::parameters xkey, CoolProp::parameters ykey, double value, const std::shared_ptr<CoolProp::AbstractState>& state);

    Range get_sat_bounds(CoolProp::parameters key) const;
    void calc_sat_range(int count);
    void update_pair(int& ipos, int& xpos, int& ypos, int& pair);
    void calc_range(std::vector<double>& xvals, std::vector<double>& yvals);

    friend class PropertyPlot;
};

using Isolines = std::vector<Isoline>;

std::vector<double> generate_values_in_range(Scale scale, const Range& range, int count);
std::vector<double> generate_values_in_range(CoolProp::parameters type, const Range& range, int count);

enum class TPLimits
{
    None,
    Def,
    Achp,
    Orc
};

class PropertyPlot
{
   public:
    CoolProp::parameters xkey;
    CoolProp::parameters ykey;
    Scale xscale;
    Scale yscale;
    Range xrange;
    Range yrange;

    PropertyPlot(const std::string& fluid_name, CoolProp::parameters ykey, CoolProp::parameters xkey, CoolProp::Plot::TPLimits tp_limits = CoolProp::Plot::TPLimits::Def);

    Range isoline_range(CoolProp::parameters key) const;
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

    CoolProp::input_pairs axis_pair_;
    bool swap_axis_inputs_for_update_;
    std::shared_ptr<CoolProp::AbstractState> state_;
    std::shared_ptr<CoolProp::AbstractState> critical_state_;
    Range2D Tp_limits_;

    Range get_sat_bounds(CoolProp::parameters key) const;
    Range2D get_Tp_limits() const;
    Range2D get_axis_limits(CoolProp::parameters xkey = CoolProp::parameters::iundefined_parameter, CoolProp::parameters ykey = CoolProp::parameters::iundefined_parameter, bool autoscale = true) const;
};

} /* namespace Plot */
} /* namespace CoolProp */

#endif /* COOLPROPPLOT_H_ */
