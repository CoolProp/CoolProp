#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

#include "CoolProp/region/AxisTransform.h"
#include "CoolProp/region/ConstantCurve.h"
#include "CoolProp/region/CubicSplineCurve.h"
#include "CoolProp/sbtl/SVDSurfaceFactory.h"
#include "CoolProp/sbtl/SatBoundaryFactory.h"
#include "CoolProp/sbtl/SurfaceSpec.h"
#include "DataStructures.h"

namespace CoolProp {
namespace sbtl {
namespace presets {

namespace {

// Shared property list for HmassP_INPUTS: h is the input (secondary
// axis), so the outputs are ρ, T, s, u.  Density is log-transformed
// (EXP) since it varies over orders of magnitude across LIQUID +
// VAPOR; the rest are identity-transformed (additive scale).
std::vector<PropertySpec> ph_properties() {
    return {
      {::CoolProp::iDmass, svd::OutputTransform::EXP},
      {::CoolProp::iT, svd::OutputTransform::IDENTITY},
      {::CoolProp::iSmass, svd::OutputTransform::IDENTITY},
      {::CoolProp::iUmass, svd::OutputTransform::IDENTITY},
    };
}

// Shared property list for PT_INPUTS: T is the input (secondary
// axis), so outputs are ρ, h, s, u.
std::vector<PropertySpec> pt_properties() {
    return {
      {::CoolProp::iDmass, svd::OutputTransform::EXP},
      {::CoolProp::iHmass, svd::OutputTransform::IDENTITY},
      {::CoolProp::iSmass, svd::OutputTransform::IDENTITY},
      {::CoolProp::iUmass, svd::OutputTransform::IDENTITY},
    };
}

}  // namespace

SurfaceSpec ph_subcritical(::CoolProp::AbstractState& heos, std::size_t NT, std::size_t NR, std::int32_t rank) {
    const auto [p_min, p_max] = subcritical_pressure_range(heos);
    const double T_min_eos = std::max(heos.Ttriple(), heos.Tmin());
    const double T_max_eos = heos.Tmax();

    SurfaceSpec spec;
    spec.fluid_name = heos.fluid_names().empty() ? std::string{} : heos.fluid_names().front();
    spec.input_pair = ::CoolProp::HmassP_INPUTS;
    spec.NT = NT;
    spec.NR = NR;
    spec.rank = rank;
    spec.properties = ph_properties();

    // LIQUID region: primary = log p, secondary = h.
    //   b_lo = h(T_min, p) — non-isothermal floor that walks T up
    //                       if T_min < T_melt(p) at high p.
    //   b_hi = h_sat,L(p) — liquid side of the saturation dome.
    {
        auto axis = region::AxisTransform::make(region::AxisScale::LOG, p_min, p_max);
        auto lo = build_h_isotherm_floor(heos, p_min, p_max, T_min_eos);
        auto hi = build_h_sat_L(heos, p_min, p_max);
        spec.regions.emplace_back(axis, std::move(lo), std::move(hi));
    }
    // VAPOR region: primary = log p, secondary = h.
    //   b_lo = h_sat,V(p) — vapor side of the saturation dome.
    //   b_hi = h(T_max - 0.5 K, p) — hot ceiling within the HEOS
    //                                validity envelope.
    {
        auto axis = region::AxisTransform::make(region::AxisScale::LOG, p_min, p_max);
        auto lo = build_h_sat_V(heos, p_min, p_max);
        auto hi = build_h_isotherm_ceiling(heos, p_min, p_max, T_max_eos - 0.5);
        spec.regions.emplace_back(axis, std::move(lo), std::move(hi));
    }

    // Thermodynamics glue: PH update + per-property readout.
    spec.update_state = [](::CoolProp::AbstractState& s, double p, double h) {
        // SurfaceSpec passes (a, b) where a is the primary axis (p) and
        // b is the secondary axis (h).  HmassP_INPUTS takes (h, p) in
        // that order.
        s.update(::CoolProp::HmassP_INPUTS, h, p);
    };
    spec.read_property = [](::CoolProp::AbstractState& s, ::CoolProp::parameters key) -> double {
        switch (key) {
            case ::CoolProp::iDmass:
                return s.rhomass();
            case ::CoolProp::iT:
                return s.T();
            case ::CoolProp::iSmass:
                return s.smass();
            case ::CoolProp::iUmass:
                return s.umass();
            default:
                throw std::invalid_argument("ph_subcritical: unsupported output property");
        }
    };

    return spec;
}

SurfaceSpec pt_subcritical(::CoolProp::AbstractState& heos, std::size_t NT, std::size_t NR, std::int32_t rank) {
    const auto [p_min, p_max] = subcritical_pressure_range(heos);
    const double T_min_eos = std::max(heos.Ttriple(), heos.Tmin());
    const double T_max_eos = heos.Tmax();

    SurfaceSpec spec;
    spec.fluid_name = heos.fluid_names().empty() ? std::string{} : heos.fluid_names().front();
    spec.input_pair = ::CoolProp::PT_INPUTS;
    spec.NT = NT;
    spec.NR = NR;
    spec.rank = rank;
    spec.properties = pt_properties();

    // For PT_INPUTS the secondary axis is T (instead of h for PH).
    // The b_lo / b_hi curves return T values, not h values.  Reuse
    // SatBoundaryFactory's T-side curves: T_sat(p) is the dome,
    // bracketed by [T_min(p), T_sat(p)] for LIQUID and [T_sat(p),
    // T_max(p)] for VAPOR.  T_min(p) walks up over a melting line
    // exactly like the h-floor case.

    // LIQUID region.  T-floor walks T_min(p) up over the melting
    // line, same trick as build_h_isotherm_floor.  Pre-compute the
    // (p, T_floor) knots inline, then hand them to CubicSplineCurve.
    {
        const SatBoundaryBuildOptions sbopts;
        const double log_p_min = std::log(p_min);
        const double log_p_max = std::log(p_max);
        std::vector<double> p_knots(sbopts.n_knots);
        std::vector<double> y(sbopts.n_knots);
        for (std::size_t k = 0; k < sbopts.n_knots; ++k) {
            const double p = std::exp(log_p_min + static_cast<double>(k) * (log_p_max - log_p_min) / static_cast<double>(sbopts.n_knots - 1));
            p_knots[k] = p;
            double T_found = std::nan("");
            for (int s = 0; s < 40; ++s) {
                const double T_try = T_min_eos + 0.5 * static_cast<double>(s);
                try {
                    heos.update(::CoolProp::PT_INPUTS, p, T_try);
                    T_found = T_try;
                    break;
                } catch (...) {  // NOLINT(bugprone-empty-catch)
                }
            }
            if (!std::isfinite(T_found)) {
                throw std::runtime_error("pt_subcritical: T-floor unreachable within 20 K of T_min_eos");
            }
            y[k] = T_found;
        }
        auto axis = region::AxisTransform::make(region::AxisScale::LOG, p_min, p_max);
        auto t_lo = region::CubicSplineCurve::build(std::move(p_knots), std::move(y));
        auto t_hi = build_T_sat(heos, p_min, p_max);
        spec.regions.emplace_back(axis, std::move(t_lo), std::move(t_hi));
    }
    // VAPOR region.
    {
        auto axis = region::AxisTransform::make(region::AxisScale::LOG, p_min, p_max);
        auto t_lo = build_T_sat(heos, p_min, p_max);
        // Constant T_max ceiling — no melting consideration on the
        // hot side, so just a flat isotherm.
        auto t_hi = std::make_unique<region::ConstantCurve>(p_min, p_max, T_max_eos - 0.5);
        spec.regions.emplace_back(axis, std::move(t_lo), std::move(t_hi));
    }

    // Thermodynamics glue: PT update + per-property readout.
    spec.update_state = [](::CoolProp::AbstractState& s, double p, double T) { s.update(::CoolProp::PT_INPUTS, p, T); };
    spec.read_property = [](::CoolProp::AbstractState& s, ::CoolProp::parameters key) -> double {
        switch (key) {
            case ::CoolProp::iDmass:
                return s.rhomass();
            case ::CoolProp::iHmass:
                return s.hmass();
            case ::CoolProp::iSmass:
                return s.smass();
            case ::CoolProp::iUmass:
                return s.umass();
            default:
                throw std::invalid_argument("pt_subcritical: unsupported output property");
        }
    };

    return spec;
}

}  // namespace presets
}  // namespace sbtl
}  // namespace CoolProp
