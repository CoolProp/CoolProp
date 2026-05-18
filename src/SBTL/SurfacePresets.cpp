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
#include "Exceptions.h"

namespace CoolProp {
namespace sbtl {
namespace presets {

namespace {

// Shared property list for HmassP_INPUTS: h is the input (secondary
// axis), so the outputs are ρ, T, s, u, w (speed of sound).  Density
// is log-transformed (EXP) since it varies over orders of magnitude
// across LIQUID + VAPOR; the rest are identity-transformed.
// Returns true if the source backend can compute the given transport
// property at a benign single-phase probe point.  Used to gate
// inclusion of η / λ in the tabulated property list — fluids whose
// HEOS model lacks a transport correlation (e.g. some siloxanes) get
// a slimmer table rather than an all-NaN sampling failure.
bool source_supports_transport_(::CoolProp::AbstractState& src, ::CoolProp::parameters key) {
    const double T = 0.85 * src.T_critical();
    const double p = 0.50 * src.p_critical();
    try {
        src.update(::CoolProp::PT_INPUTS, p, T);
        double v = (key == ::CoolProp::iviscosity) ? src.viscosity() : src.conductivity();
        return std::isfinite(v);
    } catch (...) {  // NOLINT(bugprone-empty-catch)
        return false;
    }
}

std::vector<PropertySpec> ph_properties(::CoolProp::AbstractState& src) {
    std::vector<PropertySpec> ps = {
      {::CoolProp::iDmass, svd::OutputTransform::EXP},
      {::CoolProp::iT, svd::OutputTransform::IDENTITY},
      {::CoolProp::iSmass, svd::OutputTransform::IDENTITY},
      {::CoolProp::iUmass, svd::OutputTransform::IDENTITY},
      {::CoolProp::ispeed_sound, svd::OutputTransform::IDENTITY},
    };
    // Transport properties — IAPWS G13-15 Tables 12 (η) and 13 (λ).
    // Both span ~5 orders of magnitude across the (p, h) envelope so
    // they ride the EXP transform alongside density.  Only added if
    // the source backend actually exposes them (some fluids' HEOS
    // models lack transport correlations entirely).
    if (source_supports_transport_(src, ::CoolProp::iviscosity)) {
        ps.push_back({::CoolProp::iviscosity, svd::OutputTransform::EXP});
    }
    if (source_supports_transport_(src, ::CoolProp::iconductivity)) {
        ps.push_back({::CoolProp::iconductivity, svd::OutputTransform::EXP});
    }
    return ps;
}

// Shared property list for PT_INPUTS: T is the input (secondary
// axis), so outputs are ρ, h, s, u, w (speed of sound).
std::vector<PropertySpec> pt_properties(::CoolProp::AbstractState& src) {
    std::vector<PropertySpec> ps = {
      {::CoolProp::iDmass, svd::OutputTransform::EXP},
      {::CoolProp::iHmass, svd::OutputTransform::IDENTITY},
      {::CoolProp::iSmass, svd::OutputTransform::IDENTITY},
      {::CoolProp::iUmass, svd::OutputTransform::IDENTITY},
      {::CoolProp::ispeed_sound, svd::OutputTransform::IDENTITY},
    };
    if (source_supports_transport_(src, ::CoolProp::iviscosity)) {
        ps.push_back({::CoolProp::iviscosity, svd::OutputTransform::EXP});
    }
    if (source_supports_transport_(src, ::CoolProp::iconductivity)) {
        ps.push_back({::CoolProp::iconductivity, svd::OutputTransform::EXP});
    }
    return ps;
}

}  // namespace

SurfaceSpec ph_subcritical(::CoolProp::AbstractState& heos, std::size_t NT, std::size_t NR, std::int32_t rank) {
    // Name is historical: this preset now also registers a SUPER
    // region above p_crit so the surface covers the full HEOS
    // (p, h) envelope.  Three regions: LIQUID + VAPOR + SUPER.
    const auto [p_min, p_max_sub] = subcritical_pressure_range(heos);
    const auto [p_min_sup, p_max_sup] = supercritical_pressure_range(heos);
    const double T_min_eos = std::max(heos.Ttriple(), heos.Tmin());
    const double T_max_eos = heos.Tmax();

    SurfaceSpec spec;
    spec.fluid_name = heos.fluid_names().empty() ? std::string{} : heos.fluid_names().front();
    spec.input_pair = ::CoolProp::HmassP_INPUTS;
    spec.NT = NT;
    spec.NR = NR;
    spec.rank = rank;
    spec.properties = ph_properties(heos);

    // LIQUID region: primary = log p, secondary = h.
    //   b_lo = h(T_min, p) — non-isothermal floor that walks T up
    //                       if T_min < T_melt(p) at high p.
    //   b_hi = h_sat,L(p) — liquid side of the saturation dome.
    {
        auto axis = region::AxisTransform::make(region::AxisScale::LOG, p_min, p_max_sub);
        auto lo = build_h_isotherm_floor(heos, p_min, p_max_sub, T_min_eos);
        auto hi = build_h_sat_L(heos, p_min, p_max_sub);
        spec.regions.emplace_back(axis, std::move(lo), std::move(hi));
    }
    // VAPOR region: primary = log p, secondary = h.
    //   b_lo = h_sat,V(p) — vapor side of the saturation dome.
    //   b_hi = h(T_max - 0.5 K, p) — hot ceiling within the HEOS
    //                                validity envelope.
    {
        auto axis = region::AxisTransform::make(region::AxisScale::LOG, p_min, p_max_sub);
        auto lo = build_h_sat_V(heos, p_min, p_max_sub);
        auto hi = build_h_isotherm_ceiling(heos, p_min, p_max_sub, T_max_eos - 0.5);
        spec.regions.emplace_back(axis, std::move(lo), std::move(hi));
    }
    // SUPER region: p > p_crit.  Primary = log p over [p_crit*1.001,
    // p_max_eos*0.99], secondary = h.  No sat dome above p_crit, so
    // the secondary bounds are the same low-T / high-T isotherms used
    // by LIQUID / VAPOR — they continue smoothly into supercritical.
    {
        auto axis = region::AxisTransform::make(region::AxisScale::LOG, p_min_sup, p_max_sup);
        auto lo = build_h_isotherm_floor(heos, p_min_sup, p_max_sup, T_min_eos);
        auto hi = build_h_isotherm_ceiling(heos, p_min_sup, p_max_sup, T_max_eos - 0.5);
        spec.regions.emplace_back(axis, std::move(lo), std::move(hi));
    }

    // Thermodynamics glue: PH update + per-property readout.
    //
    // For IF97-source sampling we Newton-refine the backwards-equation
    // T(p, h) against the forward h(T, p) at build time.  The IAPWS
    // R7-97 backwards equations carry a documented ±25 mK residual
    // (see IF97.h's RegionOutputBackward note), which propagates to
    // ~0.25 J/(kg·K) in s via ds ≈ cp·dT/T — well above the IAPWS
    // perm budgets of 1e-3 J/(kg·K).  Newton-refining closes that gap
    // so the SVD grid stores forward-consistent properties.  Costs
    // ~3 extra forward calls per sample (build-time only); query path
    // unchanged.  HEOS / REFPROP sources don't need refinement
    // because their HmassP path converges to forward precision
    // internally.
    //
    // Newton iterates may land within IF97's 3.3e-5 ε-band of the
    // saturation curve.  Pin the phase across iterations so PT_INPUTS
    // doesn't trip the ε-band reject.
    const bool source_is_if97 = (heos.backend_name() == "IF97Backend");
    if (source_is_if97) {
        spec.update_state = [](::CoolProp::AbstractState& s, double p, double h) {
            s.update(::CoolProp::HmassP_INPUTS, h, p);
            const ::CoolProp::phases phase0 = s.phase();
            const bool single_phase =
              (phase0 == ::CoolProp::iphase_liquid || phase0 == ::CoolProp::iphase_gas || phase0 == ::CoolProp::iphase_supercritical_liquid
               || phase0 == ::CoolProp::iphase_supercritical_gas || phase0 == ::CoolProp::iphase_supercritical);
            if (single_phase) s.specify_phase(phase0);
            for (int k = 0; k < 5; ++k) {
                const double h_now = s.hmass();
                const double dh = h_now - h;
                if (std::abs(dh) < 1e-10 * std::abs(h)) break;
                const double cp = s.cpmass();
                if (!std::isfinite(cp) || cp <= 0.0) break;
                const double T_new = s.T() - dh / cp;
                s.update(::CoolProp::PT_INPUTS, p, T_new);
            }
            if (single_phase) s.unspecify_phase();
        };
    } else {
        spec.update_state = [](::CoolProp::AbstractState& s, double p, double h) {
            // SurfaceSpec passes (a, b) where a is the primary axis (p)
            // and b is the secondary axis (h).  HmassP_INPUTS takes
            // (h, p) in that order.
            s.update(::CoolProp::HmassP_INPUTS, h, p);
        };
    }
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
            case ::CoolProp::ispeed_sound:
                return s.speed_sound();
            case ::CoolProp::iviscosity:
                return s.viscosity();
            case ::CoolProp::iconductivity:
                return s.conductivity();
            default:
                throw std::invalid_argument("ph_subcritical: unsupported output property");
        }
    };

    return spec;
}

// Build a T-floor CubicSplineCurve for any pressure range, walking T
// up at each knot to clear the melting line if needed.  Reused by
// both the subcritical LIQUID region and the SUPER region of the PT
// preset.
namespace {
std::unique_ptr<region::CubicSplineCurve> pt_T_floor_curve_(::CoolProp::AbstractState& heos, double p_min, double p_max, double T_min_eos,
                                                            const SatBoundaryBuildOptions& sbopts) {
    const double log_p_min = std::log(p_min);
    const double log_p_max = std::log(p_max);
    std::vector<double> p_knots(sbopts.n_knots);
    std::vector<double> y(sbopts.n_knots);
    const std::size_t walk_steps = sbopts.t_floor_walk_steps;
    for (std::size_t k = 0; k < sbopts.n_knots; ++k) {
        const double p = std::exp(log_p_min + static_cast<double>(k) * (log_p_max - log_p_min) / static_cast<double>(sbopts.n_knots - 1));
        p_knots[k] = p;
        double T_found = std::nan("");
        for (std::size_t s = 0; s < walk_steps; ++s) {
            const double T_try = T_min_eos + 0.5 * static_cast<double>(s);
            try {
                heos.update(::CoolProp::PT_INPUTS, p, T_try);
                T_found = T_try;
                break;
            } catch (const CoolProp::CoolPropBaseError&) {
                // Expected: below the melting line at this p, or
                // otherwise outside the HEOS envelope.  Bump T up and
                // retry; we don't swallow std::exception broadly to
                // avoid masking real bugs (allocator failures,
                // numeric_cast overflows, etc.).
            }
        }
        if (!std::isfinite(T_found)) {
            throw std::runtime_error("pt preset: T-floor unreachable within " + std::to_string(walk_steps / 2) + " K of T_min_eos");
        }
        y[k] = T_found;
    }
    return region::CubicSplineCurve::build(std::move(p_knots), std::move(y));
}
}  // namespace

SurfaceSpec pt_subcritical(::CoolProp::AbstractState& heos, std::size_t NT, std::size_t NR, std::int32_t rank) {
    // Historical name; now covers the full phase diagram via three
    // regions: LIQUID + VAPOR (subcritical) + SUPER (p > p_crit).
    const auto [p_min, p_max] = subcritical_pressure_range(heos);
    const auto [p_min_sup, p_max_sup] = supercritical_pressure_range(heos);
    const double T_min_eos = std::max(heos.Ttriple(), heos.Tmin());
    const double T_max_eos = heos.Tmax();

    SurfaceSpec spec;
    spec.fluid_name = heos.fluid_names().empty() ? std::string{} : heos.fluid_names().front();
    spec.input_pair = ::CoolProp::PT_INPUTS;
    spec.NT = NT;
    spec.NR = NR;
    spec.rank = rank;
    spec.properties = pt_properties(heos);

    // For PT_INPUTS the secondary axis is T (instead of h for PH).
    // The b_lo / b_hi curves return T values, not h values.  Reuse
    // SatBoundaryFactory's T-side curves: T_sat(p) is the dome,
    // bracketed by [T_min(p), T_sat(p)] for LIQUID and [T_sat(p),
    // T_max(p)] for VAPOR.  T_min(p) walks up over a melting line
    // exactly like the h-floor case.

    // LIQUID region.  T-floor walks T_min(p) up over the melting
    // line.  Helper handles the knot loop so the SUPER region below
    // can reuse it.
    {
        auto axis = region::AxisTransform::make(region::AxisScale::LOG, p_min, p_max);
        auto t_lo = pt_T_floor_curve_(heos, p_min, p_max, T_min_eos, SatBoundaryBuildOptions{});
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
    // SUPER region: p > p_crit.  T-floor still walks T_min up over
    // melting line (the melting curve extends well into supercritical
    // for most fluids).  Hot ceiling stays at T_max - 0.5 K.
    {
        auto axis = region::AxisTransform::make(region::AxisScale::LOG, p_min_sup, p_max_sup);
        auto t_lo = pt_T_floor_curve_(heos, p_min_sup, p_max_sup, T_min_eos, SatBoundaryBuildOptions{});
        auto t_hi = std::make_unique<region::ConstantCurve>(p_min_sup, p_max_sup, T_max_eos - 0.5);
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
            case ::CoolProp::ispeed_sound:
                return s.speed_sound();
            case ::CoolProp::iviscosity:
                return s.viscosity();
            case ::CoolProp::iconductivity:
                return s.conductivity();
            default:
                throw std::invalid_argument("pt_subcritical: unsupported output property");
        }
    };

    return spec;
}

}  // namespace presets
}  // namespace sbtl
}  // namespace CoolProp
