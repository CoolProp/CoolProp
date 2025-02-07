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

/**
A class representing a single isoline, returned by ``PropertyPlot::calc_isolines`` method. ``Isoline::x`` and ``Isoline::y`` members contain the x, y coordinates of points on the isoline in SI units of the x and y axis (respectively) of the queried PropertyPlot, and ``Isoline::value`` contains the value of the isoline in SI units.

Either ``x[i]`` or ``y[i]`` may be NaN if the isoline is not defined at that point, and should be treated as a discontinuity in the isoline.
*/
class Isoline
{
   public:
    std::vector<double> x; ///< x positions of the isoline, in SI unit of the x axis of the queried ``PropertyPlot``. May contain NaN values.
    std::vector<double> y; ///< y positions of the isoline, in SI unit of the y axis of the queried ``PropertyPlot``. May contain NaN values.
    double value;          ///< Value of the isoline in SI units.

    /**
         @brief Convenience method to get the number of points in the isoline. ``size() = x.size() = y.size()``.
         */
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

/**
     @brief Generate a list of ``count`` equally (linearly or logarithmically) spaced values in a given ``range``.

     @param scale The spacing of the values
     @param range The minimum and maximum values of generated list of values (inclusive)
     @param count The number of values to generate
     */
std::vector<double> generate_values_in_range(Scale scale, const Range& range, int count);
/**
     @brief Generate a list of ``count`` equally spaced values in a given ``range``, spaced either linearly or logarithmically based on the parameter type.

     @param type The parameter type to generate values for (usually one of the supported isoline keys, e.g. ``CoolProp::iP`` or ``CoolProp::iT``)
     @param range The minimum and maximum values of generated list of values (inclusive)
     @param count The number of values to generate
     */
std::vector<double> generate_values_in_range(CoolProp::parameters type, const Range& range, int count);

enum class TPLimits
{
    None,
    Def,
    Achp,
    Orc
};

/**
A class representing a property plot of a fluid. Used to generate isolines of a given parameter over a given pair of parameters. The code is a C++ reimplementation of the Python CoolProp.Plots.PropertyPlot class, with the actual drawing of the calculated values of the plot left to the user.

Supported plots: "p-h", "T-s", "h-s", "p-s", "p-rho", "T-rho", and "p-T".

Example usage for generating a log p-h plot for R134a with temperature isolines (for more examples, see the tests at the bottom of ``src/CoolPropPlot.cpp``):
\code{.cpp}
    // Create a plot object for R134a with pressure on the y axis and enthalpy on the x axis
    CoolProp::Plot::PropertyPlot plot("HEOS::R134a", CoolProp::iP, CoolProp::iHmass);
    // print the axis properties:
    std::cout << "x axis: h, " << (plot.xaxis.scale == CoolProp::Plot::Scale::Lin ? "lin" : "log") << ", limits [" << plot.xaxis.min << ", " << plot.xaxis.max << "] J/kg" << "\n";
    std::cout << "y axis: p, " << (plot.yaxis.scale == CoolProp::Plot::Scale::Lin ? "lin" : "log") << ", limits [" << plot.yaxis.min << ", " << plot.yaxis.max << "] Pa" << "\n";

    // Generate 5 isolines for temperature in K on this plot, 100 points per isoline
    std::vector<double> t_values = CoolProp::Plot::generate_values_in_range(CoolProp::iT, plot.isoline_range(CoolProp::iT), 5);
    CoolProp::Plot::Isolines t_isolines = plot.calc_isolines(CoolProp::iT, t_values, 100);

    // print the first temperature isoline:
    std::cout << "T: " << t_isolines[0].value << " K\n";
    for (int i = 0; i < t_isolines[0].size(); ++i) {
        std::cout << "h: " << t_isolines[0].x[i] << " J/kg, p: " << t_isolines[0].y[i] << " Pa\n";
    }
\endcode
*/
class PropertyPlot
{
   public:
    struct Axis
    {
        Scale scale;
        union
        {
            Range range;
            struct
            {
                double min, max;
            };
        };
    } xaxis,  ///< The (non-modifiable) properties of the x axis of the plot
      yaxis;  ///< The (non-modifiable) properties of the y axis of the plot

    /**
         @brief Construct a PropertyPlot object for a given fluid and a pair of parameters to plot.

         @param fluid_name The name of the fluid to plot (e.g. "HEOS::R134a", or "Water")
         @param ykey The parameter to plot on the y axis (e.g. for "ph" plots this would be ``CoolProp::iP``, and for "Ts" plots this would be ``CoolProp::iT``)
         @param xkey The parameter to plot on the x axis (e.g. for "ph" plots this would be ``CoolProp::iHmass``, and for "Ts" plots this would be ``CoolProp::iSmass``)
         @param tp_limits The temperature and pressure limits of the plot.
         */
    PropertyPlot(const std::string& fluid_name, CoolProp::parameters ykey, CoolProp::parameters xkey, CoolProp::Plot::TPLimits tp_limits = CoolProp::Plot::TPLimits::Def);

    /**
         @brief Retrieve a valid range of values (inclusive) for the isoline with the given ``key`` in SI units.
         */
    Range isoline_range(CoolProp::parameters key) const;
    /**
         @brief Calculate the isolines for the given ``key`` at given ``values`` in SI unit.

         @param key The parameter to calculate the isolines for (usually one of the supported isoline keys, e.g. ``CoolProp::iP`` or ``CoolProp::iT``. Full list can be obtained with ``PropertyPlot::supported_isoline_keys()``)
         @param values The values of the isolines to calculate in SI units
         @param points The number of points to calculate for each isoline. The larger the number, the smoother the isoline will look
         */
    Isolines calc_isolines(CoolProp::parameters key, const std::vector<double>& values, int points) const;
    /**
         @brief Retrieve a list of supported isoline keys for this plot
         */
    std::vector<CoolProp::parameters> supported_isoline_keys() const;
    /**
         @brief A method for calculating the value of a parameter at a given point in the plot. Useful for "value under cursor" type of queries. Returns NaN if the value is not defined at the given point. Convenience method for ``CoolProp::PropsSI`` based on the configuration of the plot object.

         @param key The parameter to calculate the value for (usually one of the supported isoline keys, e.g. ``CoolProp::iP`` or ``CoolProp::iT``. Full list can be obtained with ``PropertyPlot::supported_isoline_keys()``)
         @param xvalue The x coordinate of the query point in SI units of the x axis of the plot
         @param yvalue The y coordinate of the query point in SI units of the y axis of the plot
         @param phase The phase to impose for the calculation. By default, the phase is not imposed.
         */
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

    CoolProp::parameters xkey_;
    CoolProp::parameters ykey_;
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
