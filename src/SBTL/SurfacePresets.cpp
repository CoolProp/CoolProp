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

// IF97 R2/R3 boundary equation, inverted:
//   p_B23(T) = 0.10192970039326e-2 · (T − 572.54459862746)^2 + 13.91883776670  [MPa]
// Closed-form inverse for p ∈ [16.529, 100] MPa (the valid range for
// IF97's R2/R3 boundary) → T ∈ [623.15, 863.15] K.  Outside this p
// range the curve isn't defined; callers must gate.
double if97_T_B23_K(double p_Pa) {
    constexpr double a = 0.10192970039326e-2;
    constexpr double b = 0.57254459862746e3;
    constexpr double c = 0.1391883776670e2;
    const double p_MPa = p_Pa / 1.0e6;
    return b + std::sqrt((p_MPa - c) / a);
}

// Build the h = h_B23(p) boundary curve in (p, h) coords by walking p
// knots and evaluating forward IF97 h(T_B23(p), p) at each.  Used to
// split the SUPER region in ph_subcritical along the IF97 R2/R3 kink
// when the source backend is IF97 (the kink is intrinsic to IF97's
// piecewise EOS; HEOS source has no such internal boundary).  Valid
// only for p ∈ [16.529 MPa, 100 MPa] — caller must clamp.
// Build h = h_IF97(623.15 K, p) along the IF97 R1/R3 isotherm.  Used to
// split SUPER_R3 (above pcrit, below p_B23) into a R1-territory sub-
// region (T < 623.15 K, compressed-liquid supercritical) and a R3-
// proper sub-region (T >= 623.15 K, dense supercritical).  Same
// pattern as build_h_B23_curve but on the T isotherm instead of the
// p_B23(T) curve.  Valid for p > p_sat(623.15) = 16.529 MPa.
std::unique_ptr<region::CubicSplineCurve> build_h_R1R3_curve(::CoolProp::AbstractState& heos, double p_lo, double p_hi) {
    constexpr double kT_R1R3 = 623.15;
    constexpr double kP_R1R3_lo_MPa = 16.529;
    const double p_lo_clamped = std::max(p_lo, kP_R1R3_lo_MPa * 1.0e6);
    if (!(p_lo_clamped < p_hi)) {
        throw std::invalid_argument("build_h_R1R3_curve: p range does not overlap IF97 R1/R3 isotherm's validity (p > 16.529 MPa)");
    }
    constexpr std::size_t n_knots = 64;
    std::vector<double> p_knots(n_knots);
    std::vector<double> h_knots(n_knots);
    const double log_p_lo = std::log(p_lo_clamped);
    const double log_p_hi = std::log(p_hi);
    for (std::size_t k = 0; k < n_knots; ++k) {
        const double p = std::exp(log_p_lo + static_cast<double>(k) * (log_p_hi - log_p_lo) / static_cast<double>(n_knots - 1));
        heos.update(::CoolProp::PT_INPUTS, p, kT_R1R3);
        p_knots[k] = p;
        h_knots[k] = heos.hmass();
    }
    return region::CubicSplineCurve::build(std::move(p_knots), std::move(h_knots));
}

std::unique_ptr<region::CubicSplineCurve> build_h_B23_curve(::CoolProp::AbstractState& heos, double p_lo, double p_hi) {
    constexpr double kP_B23_lo_MPa = 16.529;
    constexpr double kP_B23_hi_MPa = 100.0;
    const double p_lo_clamped = std::max(p_lo, kP_B23_lo_MPa * 1.0e6);
    const double p_hi_clamped = std::min(p_hi, kP_B23_hi_MPa * 1.0e6);
    if (!(p_lo_clamped < p_hi_clamped)) {
        throw std::invalid_argument("build_h_B23_curve: p range does not overlap IF97 B23 curve's [16.529, 100] MPa validity");
    }
    constexpr std::size_t n_knots = 64;
    std::vector<double> p_knots(n_knots);
    std::vector<double> h_knots(n_knots);
    const double log_p_lo = std::log(p_lo_clamped);
    const double log_p_hi = std::log(p_hi_clamped);
    for (std::size_t k = 0; k < n_knots; ++k) {
        const double p = std::exp(log_p_lo + static_cast<double>(k) * (log_p_hi - log_p_lo) / static_cast<double>(n_knots - 1));
        const double T_B23 = if97_T_B23_K(p);
        heos.update(::CoolProp::PT_INPUTS, p, T_B23);
        p_knots[k] = p;
        h_knots[k] = heos.hmass();
    }
    return region::CubicSplineCurve::build(std::move(p_knots), std::move(h_knots));
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
    //
    // For IF97 source backend, split SUPER along the IF97 R2/R3
    // boundary curve h_B23(p).  IF97 uses different basis equations
    // (R3: ρ-T fundamental Helmholtz; R2: g-basis) on the two sides
    // of p = p_B23(T), so the property values match at the boundary
    // but the derivatives have a kink.  A single SVD over a region
    // straddling that kink encodes it into the modes and the
    // Hermite-cubic interpolant overshoots, producing R3 worst-case
    // residuals ~100% in v / ~190 K in T against IAPWS-IF97 — five
    // orders of magnitude beyond IAPWS G13-15 perm.  Splitting the
    // SUPER region at h_B23(p) puts the kink on a region edge instead
    // of inside a cell, so each sub-region's SVD only sees cells on
    // ONE side of the kink and reconstructs smoothly.  (CoolProp-foi.9.5
    // — subagent investigation 2026-05-20 ruled out kernel-level
    // fixes and converged on this atlas-level structural fix.)
    //
    // HEOS source backend has no piecewise structure, so its SUPER
    // stays one region.
    const bool source_is_if97 = (heos.backend_name() == "IF97Backend");
    constexpr double kP_B23_hi = 100.0e6;  // IF97 B23 curve's upper p bound
    if (source_is_if97 && p_min_sup < kP_B23_hi) {
        // p range where the kink applies; clamped at the B23 upper bound.
        const double p_split_hi = std::min(p_max_sup, kP_B23_hi);
        // SUPER_R1_super: IF97 R1 territory at p > pcrit, h ∈ [h_floor(p),
        // h_IF97(623.15 K, p)).  Separating R1 from R3 in this slab
        // (foi.9.9) gives the R1-side compressed-liquid SVD its own η
        // normalization — the previous combined SUPER_R3 region had
        // R1+R3 sharing modes, and R1 conformance suffered (post-foi.9.5
        // R1 worst v 0.164%).  Subcritical-p R1 (the LIQUID region)
        // hits 0.008% max v with the same physics; this split brings
        // SUPER_R3's R1 territory to that level.
        {
            auto axis = region::AxisTransform::make(region::AxisScale::LOG, p_min_sup, p_split_hi);
            auto lo = build_h_isotherm_floor(heos, p_min_sup, p_split_hi, T_min_eos);
            auto hi = build_h_R1R3_curve(heos, p_min_sup, p_split_hi);
            spec.regions.emplace_back(axis, std::move(lo), std::move(hi));
        }
        // SUPER_R3_proper: IF97 R3 territory at p > pcrit, h ∈ [h_R1R3(p),
        // h_B23(p)).  Only the (T, p, h) region where IF97 actually
        // evaluates R3.
        {
            auto axis = region::AxisTransform::make(region::AxisScale::LOG, p_min_sup, p_split_hi);
            auto lo = build_h_R1R3_curve(heos, p_min_sup, p_split_hi);
            auto hi = build_h_B23_curve(heos, p_min_sup, p_split_hi);
            spec.regions.emplace_back(axis, std::move(lo), std::move(hi));
        }
        // SUPER_R2: high-h side of the IF97 R2/R3 kink, same p range.
        {
            auto axis = region::AxisTransform::make(region::AxisScale::LOG, p_min_sup, p_split_hi);
            auto lo = build_h_B23_curve(heos, p_min_sup, p_split_hi);
            auto hi = build_h_isotherm_ceiling(heos, p_min_sup, p_split_hi, T_max_eos - 0.5);
            spec.regions.emplace_back(axis, std::move(lo), std::move(hi));
        }
        // SUPER_HIGH_P (if any): above the IF97 B23 curve's valid p
        // range there's no R2/R3 boundary (IF97 R3 extrapolates), so
        // a single region suffices.
        if (p_max_sup > kP_B23_hi) {
            auto axis = region::AxisTransform::make(region::AxisScale::LOG, kP_B23_hi, p_max_sup);
            auto lo = build_h_isotherm_floor(heos, kP_B23_hi, p_max_sup, T_min_eos);
            auto hi = build_h_isotherm_ceiling(heos, kP_B23_hi, p_max_sup, T_max_eos - 0.5);
            spec.regions.emplace_back(axis, std::move(lo), std::move(hi));
        }
    } else {
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
    if (source_is_if97) {
        spec.update_state = [](::CoolProp::AbstractState& s, double p, double h) {
            // IF97's HmassP_INPUTS only implements closed-form backward
            // T(p, h) for R1, R2, and R5.  R3 (dense supercritical /
            // near-critical) throws "Pressure out of range" because no
            // backward equation is wired up.  This used to be silently
            // masked by sample_grid's NaN-fill / median-fill machinery
            // — every R3 cell in the SUPER region got a garbage stored
            // value — until the foi.9.5 split shrank SUPER_R3's η-axis
            // so much that the bad cells dominated the surface and
            // worst-case T residuals blew out to ~470 K.
            //
            // Fallback: when HmassP_INPUTS throws, find T such that
            // h_IF97(T, p) = h by bracketed Newton iteration starting
            // from T = 700 K (mid-R3 T range).  IF97 forward h(T, p)
            // is monotonic and smooth in T at fixed p for the SUPER
            // envelope so simple Newton with a generous bracket
            // converges in 5-10 steps.
            ::CoolProp::phases phase0 = ::CoolProp::iphase_not_imposed;
            bool single_phase = false;
            bool pinned = false;
            try {
                s.update(::CoolProp::HmassP_INPUTS, h, p);
                phase0 = s.phase();
                single_phase =
                  (phase0 == ::CoolProp::iphase_liquid || phase0 == ::CoolProp::iphase_gas || phase0 == ::CoolProp::iphase_supercritical_liquid
                   || phase0 == ::CoolProp::iphase_supercritical_gas || phase0 == ::CoolProp::iphase_supercritical);
            } catch (...) {  // NOLINT(bugprone-empty-catch)
                // IF97 R3 path: start Newton from T = 700 K.  R3 spans
                // T ∈ [623.15, 863.15] K so the bracket midpoint is a
                // reasonable seed everywhere in R3.
                s.update(::CoolProp::PT_INPUTS, p, 700.0);
            }
            if (single_phase) {
                s.specify_phase(phase0);
                pinned = true;
            }
            // Phase pin must be exception-safe: AbstractState::clear()
            // doesn't reset imposed_phase_index in IF97, so a stale
            // hint can leak to the next sample cell otherwise.
            try {
                for (int k = 0; k < 10; ++k) {
                    const double h_now = s.hmass();
                    const double dh = h_now - h;
                    if (std::abs(dh) < 1e-10 * std::abs(h)) break;
                    const double cp = s.cpmass();
                    if (!std::isfinite(cp) || cp <= 0.0) break;
                    const double T_new = s.T() - dh / cp;
                    s.update(::CoolProp::PT_INPUTS, p, T_new);
                }
            } catch (...) {
                if (pinned) s.unspecify_phase();
                throw;
            }
            if (pinned) s.unspecify_phase();
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
