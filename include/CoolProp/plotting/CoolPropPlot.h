#ifndef COOLPROPPLOT_H_
#define COOLPROPPLOT_H_

/**
 * @file CoolPropPlot.h
 * @brief Lines of constant property for property diagrams, in SI units.
 *
 * The numerics of CoolPlot (github.com/CoolProp/CoolPlot, package 0.7 of
 * its road to 1.0), so that every wrapper of CoolProp can draw property
 * diagrams: a line of constant property is stepped through the
 * temperature or pressure range of a chart with steps chosen from the
 * shape of the line in plot coordinates, nothing is computed outside the
 * area the axes show, an isotherm or isobar is followed through the
 * two-phase dome in vapour quality, other lines get the saturation state
 * inserted where their phase changes, and mixtures are traced point by
 * point (a Newton iteration in temperature and density started from the
 * previous point, a secant iteration for the vapour quality inside the
 * dome) instead of flashed.
 *
 * Nothing here knows about units, styles or drawing: the result is
 * polylines in SI units with diagnostics. Drawing is left to the caller.
 */

#include "CoolProp/AbstractState.h"
#include <memory>
#include <string>
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

/// The region of temperature and pressure a chart covers, as in CoolPlot's TpLimits presets.
enum class TPLimits
{
    None,
    Def,
    Achp,
    Orc
};

/**
 * How lines are calculated.
 *
 * ``points`` is the resolution of a line: about the number of points a
 * line gets that runs across the whole chart; bends get more, straight
 * parts fewer. ``solver`` is "auto" (trace mixtures, flash pure fluids),
 * "tracer" or "flash". ``critical_margin`` is the relative distance kept
 * from the triple point and the critical point on the saturation curve.
 * ``saturation_spacing`` of the quality lines of a pure fluid is
 * "critical" (denser towards the critical point) or "linear".
 */
struct CalcSettings
{
    int points = 250;
    int saturation_points = 250;
    std::string solver = "auto";
    double critical_margin = 1e-7;
    std::string saturation_spacing = "critical";
};

/**
 * One line of constant property: ``x`` and ``y`` in the SI units of the
 * axes of the ``PropertyPlot`` that made it, with diagnostics.
 *
 * ``x[i]`` and ``y[i]`` are NaN where a point could not be computed; draw
 * a gap there. ``n_failed`` counts those points, ``n_fallback`` the points
 * the tracer handed to a plain flash, ``n_inserted`` the saturation states
 * inserted where the line crosses the phase boundary, ``n_evaluated`` the
 * states computed for the line and ``n_outside`` how many of them lay
 * outside the area the axes show. ``method`` is "flash", "tracer" or
 * "envelope"; ``reason`` says why a line could not be computed at all or
 * why it has fewer than two points.
 */
class Isoline
{
   public:
    std::vector<double> x;
    std::vector<double> y;
    double value = 0.0;
    int n_failed = 0;
    int n_fallback = 0;
    int n_inserted = 0;
    int n_evaluated = 0;
    int n_outside = 0;
    std::string method = "flash";
    std::string reason;

    /// The number of points; ``x.size() == y.size()``.
    size_t size() const {
        return x.size();
    }
    /// The number of points with finite x and y.
    size_t n_valid() const;
    /// Invalid points between the first and the last valid one.
    size_t n_gaps() const;
};

using Isolines = std::vector<Isoline>;

/// ``count`` values spaced linearly or logarithmically over ``range``, both ends included.
std::vector<double> generate_values_in_range(Scale scale, const Range& range, int count);
/// ``count`` values spaced as the axis of ``type`` is scaled (pressure and density logarithmically).
std::vector<double> generate_values_in_range(CoolProp::parameters type, const Range& range, int count);

namespace Detail {
class Fluid;
Scale default_scale(CoolProp::parameters key);
/// The state of a fluid string such as "HEOS::R32[0.7]&R125[0.3]", with its fractions set.
std::shared_ptr<CoolProp::AbstractState> process_fluid_state(const std::string& fluid_ref);
}  // namespace Detail

/**
 * A property diagram of one fluid: which quantity is on which axis, the
 * region covered, and the lines of constant property over it.
 *
 * Supported diagrams, vertical axis first: p-h, T-s, h-s, p-s, p-rho,
 * T-rho and p-T. The example draws the isotherms of a log(p)-h diagram
 * of R134a:
 * \code{.cpp}
 *     CoolProp::Plot::PropertyPlot plot("HEOS::R134a", CoolProp::iP, CoolProp::iHmass);
 *     std::vector<double> T_values = CoolProp::Plot::generate_values_in_range(CoolProp::iT, plot.isoline_range(CoolProp::iT), 5);
 *     CoolProp::Plot::Isolines isotherms = plot.calc_isolines(CoolProp::iT, T_values, 250);
 *     for (size_t i = 0; i < isotherms[0].size(); ++i) {
 *         std::cout << isotherms[0].x[i] << " J/kg, " << isotherms[0].y[i] << " Pa\n";
 *     }
 * \endcode
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
    } xaxis,  ///< The horizontal axis: scale and range in SI units
      yaxis;  ///< The vertical axis

    /**
     * @param fluid_name A CoolProp fluid string, e.g. "HEOS::R134a", "R32[0.7]&R125[0.3]" or "INCOMP::MEG[0.3]"
     * @param ykey The quantity on the vertical axis (``CoolProp::iP`` for a p-h diagram)
     * @param xkey The quantity on the horizontal axis (``CoolProp::iHmass`` for a p-h diagram)
     * @param tp_limits The temperature and pressure region the lines cover
     */
    PropertyPlot(const std::string& fluid_name, CoolProp::parameters ykey, CoolProp::parameters xkey,
                 CoolProp::Plot::TPLimits tp_limits = CoolProp::Plot::TPLimits::Def);
    /// The same, with absolute limits: ``T_min_K, T_max_K, p_min_Pa, p_max_Pa`` (an incompressible fluid needs them).
    PropertyPlot(const std::string& fluid_name, CoolProp::parameters ykey, CoolProp::parameters xkey, double T_min_K, double T_max_K, double p_min_Pa,
                 double p_max_Pa);
    ~PropertyPlot();
    PropertyPlot(const PropertyPlot&) = delete;
    PropertyPlot& operator=(const PropertyPlot&) = delete;

    /// The range of values a line of constant ``key`` can take over the region, in SI units.
    Range isoline_range(CoolProp::parameters key) const;
    /// The lines of constant ``key`` at the given values in SI units, with the resolution ``points``.
    Isolines calc_isolines(CoolProp::parameters key, const std::vector<double>& values, int points) const;
    /// The same with every calculation setting.
    Isolines calc_isolines(CoolProp::parameters key, const std::vector<double>& values, const CalcSettings& settings) const;
    /// The quantities that can be drawn as lines in this diagram.
    std::vector<CoolProp::parameters> supported_isoline_keys() const;
    /// The quantity that is stepped along a line of constant ``key`` (``iT`` or ``iP``), or ``iundefined_parameter``.
    CoolProp::parameters sweep_key(CoolProp::parameters key) const;
    /// The value of ``key`` at a point of the diagram (the state under the cursor), NaN if there is none.
    double value_at(CoolProp::parameters key, double xvalue, double yvalue, CoolProp::phases phase = CoolProp::phases::iphase_not_imposed) const;
    /// The temperature and pressure region in SI units: T_min, T_max, p_min, p_max.
    std::vector<double> tp_box() const;
    /// Whether the fluid is a mixture.
    bool is_mixture() const;

   private:
    std::unique_ptr<Detail::Fluid> fluid_;
    CoolProp::parameters xkey_;
    CoolProp::parameters ykey_;
    CoolProp::input_pairs axis_pair_;
    bool swap_axis_inputs_for_update_;
    double T_min_K_, T_max_K_, p_min_Pa_, p_max_Pa_;

    void init(const std::string& fluid_name);
    void finish_limits(double T_min, double T_max, double p_min, double p_max, bool relative);
    Range axis_range(CoolProp::parameters key) const;
};

}  // namespace Plot
}  // namespace CoolProp

#endif /* COOLPROPPLOT_H_ */
