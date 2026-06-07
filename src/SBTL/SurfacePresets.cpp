#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <utility>

#include <IF97.h>
#include "boost/math/tools/toms748_solve.hpp"

#include "CoolProp/region/AxisTransform.h"
#include "CoolProp/region/ConstantCurve.h"
#include "CoolProp/region/CubicSplineCurve.h"
#include "CoolProp/sbtl/SVDSurfaceFactory.h"
#include "CoolProp/sbtl/SatBoundaryFactory.h"
#include "CoolProp/sbtl/SurfaceSpec.h"
#include "CoolProp/DataStructures.h"
#include "CoolProp/Exceptions.h"

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

// NOTE on the IF97 R2/R3 boundary: T_B23(p) comes straight from the
// if97 library's IF97::Region23_p(p) — no local copy of the B23
// coefficients (henningjp, PR #2938).  CoolProp builds the if97 dep
// in Pa units (p_fact = 1e6), matching the Pa-based IF97::psat97 /
// Tsat97 used in HumidAirProp, so call sites pass p in Pa directly.
// IF97::Region23_p is valid for p ∈ [16.529, 100] MPa → T ∈
// [623.15, 863.15] K; callers must gate to that range.

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
        const double T_B23 = IF97::Region23_p(p);
        heos.update(::CoolProp::PT_INPUTS, p, T_B23);
        p_knots[k] = p;
        h_knots[k] = heos.hmass();
    }
    return region::CubicSplineCurve::build(std::move(p_knots), std::move(h_knots));
}

// Find T such that h(T, p) = h_target via bracketed TOMS748 iteration on
// the source backend's forward h.  Bracketed (vs Newton) so we stay
// robust across IF97 region boundaries where cp jumps cause Newton to
// overshoot into the wrong basin — the foi.9.10 root cause that
// contaminated R3 cells with garbage stored T / s / ρ values up through
// foi.9.9.  Boost's toms748_solve is already used elsewhere in CoolProp
// (superancillary, AbstractState, FlashRoutines), so no new dependency.
//
// Returns the converged T.  If the bracket [T_lo, T_hi] doesn't contain
// a root, returns the nearer endpoint and signals via `ok` so the caller
// can fall through to the NaN-fill machinery downstream.
double solve_T_from_h_toms748(::CoolProp::AbstractState& s, double p, double h_target, double T_lo, double T_hi, bool* ok = nullptr) {
    auto resid = [&s, p, h_target](double T) -> double {
        try {
            s.update(::CoolProp::PT_INPUTS, p, T);
            return s.hmass() - h_target;
        } catch (...) {  // NOLINT(bugprone-empty-catch)
            // Out-of-range / ε-band rejects bubble up as NaN so the
            // sign-bracket check below cleanly fails the bracket test.
            return std::nan("");
        }
    };
    const double r_lo = resid(T_lo);
    const double r_hi = resid(T_hi);
    if (!std::isfinite(r_lo) || !std::isfinite(r_hi) || r_lo * r_hi > 0.0) {
        if (ok != nullptr) {
            *ok = false;
        }
        // Returning the endpoint with smaller residual would re-mutate
        // `s`; callers that want to fall through to the HmassP_INPUTS
        // state must check *ok and skip the subsequent PT_INPUTS update.
        return (std::abs(r_lo) < std::abs(r_hi)) ? T_lo : T_hi;
    }
    std::uintmax_t max_iter = 30;
    const auto bracket = boost::math::tools::toms748_solve(resid, T_lo, T_hi, r_lo, r_hi, boost::math::tools::eps_tolerance<double>(40), max_iter);
    if (ok != nullptr) {
        *ok = true;
    }
    return 0.5 * (bracket.first + bracket.second);
}

}  // namespace

SurfaceSpec ph_subcritical(::CoolProp::AbstractState& heos, std::size_t NT, std::size_t NR, std::int32_t rank) {
    // Name is historical: this preset now also registers a SUPER
    // region above p_crit so the surface covers the full HEOS
    // (p, h) envelope.  Three regions: LIQUID + VAPOR + SUPER, plus
    // an optional NC_LIQUID/NC_VAPOR pair near pc when source is
    // HEOS/REFPROP (CoolProp-4u9).
    const auto [p_min, p_max_sub_raw] = subcritical_pressure_range(heos);
    const auto [p_min_sup_raw, p_max_sup] = supercritical_pressure_range(heos);
    const double T_min_eos = std::max(heos.Ttriple(), heos.Tmin());
    const double T_max_eos = heos.Tmax();
    const double pc = heos.p_critical();

    // Near-critical sub-region geometry (CoolProp-4u9).  HEOS/REFPROP
    // sources gain a dedicated NC_LIQUID/NC_VAPOR pair on p ∈ [0.9·pc,
    // (1 − 1 ppm)·pc] with a POWER(β=1/3) primary axis that crowds
    // grid lines toward pc.  Python diagnostics on R245fa/Water/CO2
    // showed POWER buys 1000–100000× better off-grid max error vs LOG
    // and stays flat in p_hi from 0.999·pc to 0.999999·pc — the
    // critical-scaling singularity gets absorbed by the cube-root
    // transform (Ising β ≈ 0.326 ≈ 1/3).  IF97 source skips this
    // (G13-15 already passes without it).
    const bool source_is_if97 = (heos.backend_name() == "IF97Backend");
    constexpr double kNC_p_lo_mult = 0.9;
    constexpr double kNC_p_hi_mult = 1.0 - 1.0e-10;      // gap of 1e-10·pc below pc — sub-ULP for practical p
    constexpr double kNC_sup_p_lo_mult = 1.0 + 1.0e-10;  // mirror above pc
    constexpr double kNC_sup_p_hi_mult = 1.1;
    const bool nc_enabled = !source_is_if97;
    const double p_nc_lo = nc_enabled ? kNC_p_lo_mult * pc : 0.0;
    const double p_nc_hi = nc_enabled ? kNC_p_hi_mult * pc : 0.0;
    const double p_nc_sup_lo = nc_enabled ? kNC_sup_p_lo_mult * pc : 0.0;
    const double p_nc_sup_hi = nc_enabled ? kNC_sup_p_hi_mult * pc : 0.0;
    // Clip parent LIQUID/VAPOR p-max to the NC band's lower edge when
    // NC is enabled — handoff at p = 0.9·pc.  Clip parent SUPER p-min
    // to NC_SUPER's upper edge.  When NC disabled (IF97 source), parent
    // bounds use legacy values.
    const double p_max_sub = nc_enabled ? std::min(p_max_sub_raw, p_nc_lo) : p_max_sub_raw;
    const double p_min_sup = nc_enabled ? std::max(p_min_sup_raw, p_nc_sup_hi) : p_min_sup_raw;

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
    if (p_min < p_max_sub) {
        auto axis = region::AxisTransform::make(region::AxisScale::LOG, p_min, p_max_sub);
        auto lo = build_h_isotherm_floor(heos, p_min, p_max_sub, T_min_eos);
        auto hi = build_h_sat_L(heos, p_min, p_max_sub);
        spec.regions.emplace_back(axis, std::move(lo), std::move(hi));
    }
    // VAPOR region: primary = log p, secondary = h.
    //   b_lo = h_sat,V(p) — vapor side of the saturation dome.
    //   b_hi = h(T_max - 0.5 K, p) — hot ceiling within the HEOS
    //                                validity envelope.
    if (p_min < p_max_sub) {
        auto axis = region::AxisTransform::make(region::AxisScale::LOG, p_min, p_max_sub);
        auto lo = build_h_sat_V(heos, p_min, p_max_sub);
        auto hi = build_h_isotherm_ceiling(heos, p_min, p_max_sub, T_max_eos - 0.5);
        spec.regions.emplace_back(axis, std::move(lo), std::move(hi));
    }
    // NC_LIQUID + NC_VAPOR (CoolProp-4u9): near-critical sub-regions on
    // p ∈ [0.9·pc, (1 − 1e-10)·pc] with POWER(β=1/3) primary axis.
    // The 1e-10 gap below pc is sub-ULP for practical p — anything
    // closer falls inside the critical patch's auto-cal'd (T, p) bbox.
    // Same boundary curves as the parent regions — the SuperAncillary-
    // backed sat curves are continuous across the p = 0.9·pc handoff.
    if (nc_enabled && p_nc_lo < p_nc_hi) {
        {
            auto axis = region::AxisTransform::make(region::AxisScale::POWER, p_nc_lo, p_nc_hi);
            auto lo = build_h_isotherm_floor(heos, p_nc_lo, p_nc_hi, T_min_eos);
            auto hi = build_h_sat_L(heos, p_nc_lo, p_nc_hi);
            spec.regions.emplace_back(axis, std::move(lo), std::move(hi));
        }
        {
            auto axis = region::AxisTransform::make(region::AxisScale::POWER, p_nc_lo, p_nc_hi);
            auto lo = build_h_sat_V(heos, p_nc_lo, p_nc_hi);
            auto hi = build_h_isotherm_ceiling(heos, p_nc_lo, p_nc_hi, T_max_eos - 0.5);
            spec.regions.emplace_back(axis, std::move(lo), std::move(hi));
        }
    }
    // NC_SUPER (CoolProp-4u9): super-critical near-critical sub-region
    // on p ∈ [(1 + 1e-10)·pc, 1.1·pc] with POWER_LO(β=1/3) primary
    // axis that crowds grid lines toward pc.  No saturation dome above
    // pc, so b_lo / b_hi are the same isotherm-floor / isotherm-ceiling
    // pair used by parent SUPER — they continue smoothly across the
    // p = 1.1·pc handoff.  Single region (not the LIQ/VAP split used
    // sub-critically) because the dome doesn't bisect this band.
    // Closes the (h, p) coverage gap exposed by the h=350 kJ/kg sweep
    // for R245fa: cells at p just above pc with T outside the patch's
    // (T, p) bbox previously fell through every region and returned
    // NaN.  NC_SUPER claims them; SVD residual ≲ 1e-7 in this band.
    if (nc_enabled && p_nc_sup_lo < p_nc_sup_hi) {
        auto axis = region::AxisTransform::make(region::AxisScale::POWER_LO, p_nc_sup_lo, p_nc_sup_hi);
        auto lo = build_h_isotherm_floor(heos, p_nc_sup_lo, p_nc_sup_hi, T_min_eos);
        auto hi = build_h_isotherm_ceiling(heos, p_nc_sup_lo, p_nc_sup_hi, T_max_eos - 0.5);
        spec.regions.emplace_back(axis, std::move(lo), std::move(hi));
    }
    // SUPER region: p > p_crit (clipped up to 1.1·pc when NC enabled,
    // so NC_SUPER owns the strip just above pc).  Primary = log p over
    // [max(p_crit, 1.1·pc), p_max_eos*0.99], secondary = h.  No sat
    // dome above p_crit, so
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
        // SUPER_R5: IAPWS R7-97 Region 5 — extended-T high-temperature
        // slab.  T ∈ [1073.15, 2273.15] K, p ∈ (0, 50] MPa.  Closes
        // CoolProp-pd6.  All other SVDSBTL&IF97 regions are bounded
        // above by IF97::get_Tmax() = 1073.15 K (the R1/R2/R3 ceiling),
        // so without this slab the conformance plots have a hard
        // cutoff at 1073 K and the G13-15 Table 11 claim for R5 isn't
        // actually exercised.  IF97Backend dispatches PT_INPUTS with
        // T > 1073.15 K (and p ≤ 50 MPa) to R5 internally, so the
        // sampling reuses the same update_state lambda as the rest of
        // the IF97 path.  No internal subdivision: R5 is one closed-
        // form expansion in IF97 and smooth enough to fit cleanly in a
        // single SVD region.
        // Nudge above 1073.15 because IF97's RegionDetermination_TP
        // dispatches T == 1073.15 K to R2 (the boundary uses
        // T <= Tmax strictly).  Sampling h at T=1073.15 would build
        // the SVD on R2's reference frame, then runtime queries that
        // bracket at T=1073.15 would converge to R2 instead of R5
        // (R2 and R5 use independent ideal-gas references; no
        // IAPWS R7-97 smoothness condition at the boundary).  One
        // ULP above 1073.15 puts the sample and the runtime bracket
        // both squarely in R5 territory.  Found by code-reviewer
        // adversarial pass.
        const double kT_R5_lo = std::nextafter(1073.15, 2273.15);
        constexpr double kT_R5_hi = 2273.15;
        constexpr double kP_R5_hi = 50.0e6;
        // IF97 Pmin = 611.213 Pa (saturation pressure at T_min =
        // 273.15 K) per IAPWS R7-97; queries below that throw
        // "Pressure out of range".  Use 700 Pa as the SVDSBTL R5
        // lower bound to leave a small margin against numerical
        // boundary effects in the log-uniform grid sampling near
        // p_min.  Practical cost: ~50 Pa of R5 coverage lost on the
        // very low-p side, where nothing physically interesting
        // happens (R5 is supercritical-vapor at extreme T).
        constexpr double kP_R5_lo = 700.0;
        if (kP_R5_lo < kP_R5_hi) {
            const double p_r5_hi = std::min(p_max_sup, kP_R5_hi);
            auto axis = region::AxisTransform::make(region::AxisScale::LOG, kP_R5_lo, p_r5_hi);
            // h_lo / h_hi sampled at the R5 endpoints via the same
            // isotherm-spline path the rest of the preset uses.  R5
            // is monotone in h(T) at fixed p, so the spline through
            // 64 knots reproduces the slab cleanly.
            auto lo = build_h_isotherm_floor(heos, kP_R5_lo, p_r5_hi, kT_R5_lo);
            auto hi = build_h_isotherm_ceiling(heos, kP_R5_lo, p_r5_hi, kT_R5_hi);
            spec.regions.emplace_back(axis, std::move(lo), std::move(hi));
        }
    } else if (p_min_sup < p_max_sup) {
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
            // Two-stage h-inversion via TOMS748 (foi.9.10).  Bracketed
            // iteration is essential — Newton from a fixed T=700 K seed
            // overshoots across IF97 R2/R3 cp-jump boundaries when
            // T_B23(p) < 700 K and fails to converge, leaving stored
            // cell values wrong by tens of K (e.g., 65 K T error at
            // T_target=663.7, p=26.6 MPa, propagating to s, ρ residuals
            // of >20% downstream).
            //
            //   1. HmassP_INPUTS (R7-97 backward) seeds T to ±25 mK,
            //      sets the phase classification.  Throws in R3 where
            //      no backward equation is published.
            //   2. TOMS748 polish on [T_seed − 0.5, T_seed + 0.5], clamped
            //      to IF97's [273.16, 2273.15] T validity, closes the
            //      ±25 mK R7-97 residual to machine precision.  Critical
            //      for R2/R5 s conformance (cp·dT/T propagates 25 mK in
            //      T to ~0.25 J/(kg·K) in s — over the 1e-3 perm budget).
            //      If the clamped bracket doesn't span a sign change
            //      (rare: HmassP already at machine precision, or
            //      near-T_min cell with bracket clipped tight), the
            //      helper signals via *ok and we keep the HmassP state.
            //   3. R3 fallback: TOMS748 on [623.15 K, T_B23(p)].
            ::CoolProp::phases phase0 = ::CoolProp::iphase_not_imposed;
            bool single_phase = false;
            bool seeded_via_hmass = true;
            try {
                s.update(::CoolProp::HmassP_INPUTS, h, p);
                phase0 = s.phase();
                single_phase =
                  (phase0 == ::CoolProp::iphase_liquid || phase0 == ::CoolProp::iphase_gas || phase0 == ::CoolProp::iphase_supercritical_liquid
                   || phase0 == ::CoolProp::iphase_supercritical_gas || phase0 == ::CoolProp::iphase_supercritical);
            } catch (...) {  // NOLINT(bugprone-empty-catch)
                seeded_via_hmass = false;
            }
            // Phase pin (PR E mechanism) so the TOMS748 polish's
            // PT_INPUTS evaluations near the sat dome don't trip
            // IF97's ε-band reject.  Exception-safe — AbstractState::
            // clear() doesn't reset imposed_phase_index in IF97.
            bool pinned = false;
            if (single_phase) {
                s.specify_phase(phase0);
                pinned = true;
            }
            try {
                if (!seeded_via_hmass) {
                    // Two fallback paths for cells where the R7-97
                    // backward equation can't dispatch:
                    //   * R3 — T ∈ [623.15, T_B23(p)], p ∈ [16.529,
                    //     100] MPa (existing path).
                    //   * R5 — T ∈ [1073.15, 2273.15], p ≤ 50 MPa
                    //     (added for CoolProp-pd6 to support the new
                    //     SUPER_R5 region).  R5 has no backward
                    //     equation either, so we forward-solve via
                    //     TOMS748 just like R3.
                    //
                    // Dispatch by h_target magnitude — R5's lowest h
                    // (at T=1073.15) is ~4 MJ/kg, well above R3's
                    // upper h (~2.5 MJ/kg at the R3/R2 boundary).
                    // 3 MJ/kg is a conservative split.
                    // kR5_T_lo nudged one ULP above 1073.15 to skip
                    // IF97's R2/R5 dispatch boundary — see the
                    // matching comment in the SUPER_R5 region build.
                    // Without this, resid(T_lo) evaluates R2 while
                    // resid(T_lo + ε) evaluates R5, and the bracket
                    // straddles a small h discontinuity between the
                    // two formulations.
                    const double kR5_T_lo = std::nextafter(1073.15, 2273.15);
                    constexpr double kR5_T_hi = 2273.15;
                    constexpr double kR5_P_hi = 50.0e6;
                    constexpr double kH_R3_vs_R5_split = 3.0e6;
                    constexpr double kP_B23_lo_Pa = 16.529e6;
                    constexpr double kP_B23_hi_Pa = 100.0e6;

                    bool resolved = false;
                    // R5 fallback (try first when h is high and p is
                    // in R5's slab).
                    if (h >= kH_R3_vs_R5_split && p <= kR5_P_hi) {
                        bool ok = false;
                        const double T_R5 = solve_T_from_h_toms748(s, p, h, kR5_T_lo, kR5_T_hi, &ok);
                        if (ok) {
                            s.update(::CoolProp::PT_INPUTS, p, T_R5);
                            resolved = true;
                        }
                    }
                    if (!resolved) {
                        // R3 fallback: bracket on [623.15 K, T_B23(p)].
                        // IF97::Region23_p is only valid for p ∈ [16.529,
                        // 100] MPa; outside that range it returns garbage
                        // (sub-623.15 K or NaN from the negative-argument
                        // sqrt).  Above 100 MPa is the SUPER_HIGH_P slab
                        // where there's no R2/R3 boundary at all — the R3
                        // fallback shouldn't run there in the first place,
                        // but guard anyway.
                        if (!(p >= kP_B23_lo_Pa && p <= kP_B23_hi_Pa)) {
                            throw ::CoolProp::ValueError(
                              "ph_subcritical IF97 sampling: neither R5 (h<3MJ/kg or p>50MPa) nor R3 (p outside [16.529, 100] MPa) "
                              "fallback applies — HmassP backward seed failed and no bracket is defined");
                        }
                        const double T_R3_lo = 623.15;
                        const double T_R3_hi = IF97::Region23_p(p);
                        if (!(T_R3_lo < T_R3_hi)) {
                            throw ::CoolProp::ValueError("ph_subcritical IF97 sampling: degenerate R3 bracket (T_B23(p) <= 623.15 K)");
                        }
                        bool ok = false;
                        const double T_R3 = solve_T_from_h_toms748(s, p, h, T_R3_lo, T_R3_hi, &ok);
                        if (!ok) {
                            // Bracket doesn't contain h — propagate as
                            // failure so the SurfaceFactory marks this
                            // cell NaN, rather than committing an endpoint
                            // T as a real state.
                            throw ::CoolProp::ValueError("ph_subcritical IF97 sampling: R3 TOMS748 missed bracket; cell will be NaN-filled");
                        }
                        s.update(::CoolProp::PT_INPUTS, p, T_R3);
                    }
                } else {
                    // Polish R7-97 backward seed.
                    const double T_seed = s.T();
                    constexpr double kIF97_T_min = 273.16;
                    constexpr double kIF97_T_max = 2273.15;
                    const double T_lo = std::max(T_seed - 0.5, kIF97_T_min);
                    const double T_hi = std::min(T_seed + 0.5, kIF97_T_max);
                    if (T_lo < T_hi) {
                        bool ok = false;
                        const double T_polished = solve_T_from_h_toms748(s, p, h, T_lo, T_hi, &ok);
                        if (ok) {
                            s.update(::CoolProp::PT_INPUTS, p, T_polished);
                        } else {
                            // solve_T_from_h_toms748 mutates `s` via its
                            // bracket-residual probes; on a miss the
                            // last probe leaves s at PT_INPUTS(p, T_hi
                            // or T_lo) instead of the HmassP seed.
                            // Restore so callers downstream see the
                            // backward-equation seed (±25 mK floor)
                            // rather than a near-endpoint state.
                            try {
                                s.update(::CoolProp::HmassP_INPUTS, h, p);
                            } catch (...) {  // NOLINT(bugprone-empty-catch)
                                // Seed-restore failed (the original
                                // HmassP at line 362 succeeded, so this
                                // is unexpected); leave state alone and
                                // let downstream handle as a bad cell.
                            }
                        }
                    }
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
            } catch (const CoolProp::CoolPropBaseError&) {  // NOLINT(bugprone-empty-catch)
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
    // regions: LIQUID + VAPOR (subcritical) + SUPER (p > p_crit), plus
    // an optional NC_LIQUID/NC_VAPOR pair near pc when source is
    // HEOS/REFPROP (CoolProp-4u9).
    const auto [p_min, p_max_raw] = subcritical_pressure_range(heos);
    const auto [p_min_sup_raw, p_max_sup] = supercritical_pressure_range(heos);
    const double T_min_eos = std::max(heos.Ttriple(), heos.Tmin());
    const double T_max_eos = heos.Tmax();
    const double pc = heos.p_critical();

    // Near-critical sub-region geometry (mirrors ph_subcritical above).
    const bool source_is_if97 = (heos.backend_name() == "IF97Backend");
    constexpr double kNC_p_lo_mult = 0.9;
    constexpr double kNC_p_hi_mult = 1.0 - 1.0e-10;
    constexpr double kNC_sup_p_lo_mult = 1.0 + 1.0e-10;
    constexpr double kNC_sup_p_hi_mult = 1.1;
    const bool nc_enabled = !source_is_if97;
    const double p_nc_lo = nc_enabled ? kNC_p_lo_mult * pc : 0.0;
    const double p_nc_hi = nc_enabled ? kNC_p_hi_mult * pc : 0.0;
    const double p_nc_sup_lo = nc_enabled ? kNC_sup_p_lo_mult * pc : 0.0;
    const double p_nc_sup_hi = nc_enabled ? kNC_sup_p_hi_mult * pc : 0.0;
    const double p_max = nc_enabled ? std::min(p_max_raw, p_nc_lo) : p_max_raw;
    const double p_min_sup = nc_enabled ? std::max(p_min_sup_raw, p_nc_sup_hi) : p_min_sup_raw;

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
    if (p_min < p_max) {
        auto axis = region::AxisTransform::make(region::AxisScale::LOG, p_min, p_max);
        auto t_lo = pt_T_floor_curve_(heos, p_min, p_max, T_min_eos, SatBoundaryBuildOptions{});
        auto t_hi = build_T_sat(heos, p_min, p_max);
        spec.regions.emplace_back(axis, std::move(t_lo), std::move(t_hi));
    }
    // VAPOR region.
    if (p_min < p_max) {
        auto axis = region::AxisTransform::make(region::AxisScale::LOG, p_min, p_max);
        auto t_lo = build_T_sat(heos, p_min, p_max);
        // Constant T_max ceiling — no melting consideration on the
        // hot side, so just a flat isotherm.
        auto t_hi = std::make_unique<region::ConstantCurve>(p_min, p_max, T_max_eos - 0.5);
        spec.regions.emplace_back(axis, std::move(t_lo), std::move(t_hi));
    }
    // NC_LIQUID + NC_VAPOR (CoolProp-4u9): near-critical sub-regions on
    // p ∈ [0.9·pc, (1 − 1e-10)·pc] with POWER(β=1/3) primary axis.
    // The 1e-10 gap below pc is sub-ULP for practical p — anything
    // closer falls inside the critical patch's auto-cal'd bbox.
    if (nc_enabled && p_nc_lo < p_nc_hi) {
        {
            auto axis = region::AxisTransform::make(region::AxisScale::POWER, p_nc_lo, p_nc_hi);
            auto t_lo = pt_T_floor_curve_(heos, p_nc_lo, p_nc_hi, T_min_eos, SatBoundaryBuildOptions{});
            auto t_hi = build_T_sat(heos, p_nc_lo, p_nc_hi);
            spec.regions.emplace_back(axis, std::move(t_lo), std::move(t_hi));
        }
        {
            auto axis = region::AxisTransform::make(region::AxisScale::POWER, p_nc_lo, p_nc_hi);
            auto t_lo = build_T_sat(heos, p_nc_lo, p_nc_hi);
            auto t_hi = std::make_unique<region::ConstantCurve>(p_nc_lo, p_nc_hi, T_max_eos - 0.5);
            spec.regions.emplace_back(axis, std::move(t_lo), std::move(t_hi));
        }
    }
    // NC_SUPER (CoolProp-4u9): super-critical near-critical sub-region
    // on p ∈ [(1 + 1e-10)·pc, 1.1·pc] with POWER_LO(β=1/3) primary
    // axis crowding toward pc.  Single region (no sat dome above pc);
    // b_lo / b_hi are the T-floor / constant T_max-0.5 ceiling.
    if (nc_enabled && p_nc_sup_lo < p_nc_sup_hi) {
        auto axis = region::AxisTransform::make(region::AxisScale::POWER_LO, p_nc_sup_lo, p_nc_sup_hi);
        auto t_lo = pt_T_floor_curve_(heos, p_nc_sup_lo, p_nc_sup_hi, T_min_eos, SatBoundaryBuildOptions{});
        auto t_hi = std::make_unique<region::ConstantCurve>(p_nc_sup_lo, p_nc_sup_hi, T_max_eos - 0.5);
        spec.regions.emplace_back(axis, std::move(t_lo), std::move(t_hi));
    }
    // SUPER region: p > p_crit (clipped up to 1.1·pc when NC enabled
    // so NC_SUPER owns the strip just above pc).  T-floor still walks
    // T_min up over melting line (the melting curve extends well into
    // supercritical for most fluids).  Hot ceiling stays at T_max-0.5 K.
    if (p_min_sup < p_max_sup) {
        auto axis = region::AxisTransform::make(region::AxisScale::LOG, p_min_sup, p_max_sup);
        auto t_lo = pt_T_floor_curve_(heos, p_min_sup, p_max_sup, T_min_eos, SatBoundaryBuildOptions{});
        auto t_hi = std::make_unique<region::ConstantCurve>(p_min_sup, p_max_sup, T_max_eos - 0.5);
        spec.regions.emplace_back(axis, std::move(t_lo), std::move(t_hi));
    }
    // SUPER_R5: IAPWS R7-97 Region 5 — extended-T high-temperature
    // slab.  T ∈ [1073.15, 2273.15] K, p ∈ (0, 50] MPa.  Closes
    // CoolProp-pd6.  Mirror of the SUPER_R5 region in ph_subcritical.
    // IF97-only because R5 has no HEOS counterpart in CoolProp for
    // water (Wagner & Pruss covers 273-1273 K).  Constant-T b_lo /
    // b_hi: R5 is a flat T-slab, no melting curve consideration.
    if (source_is_if97) {
        // Nudge above 1073.15 because IF97's RegionDetermination_TP
        // dispatches T == 1073.15 K to R2 (the boundary uses
        // T <= Tmax strictly).  Sampling h at T=1073.15 would build
        // the SVD on R2's reference frame, then runtime queries that
        // bracket at T=1073.15 would converge to R2 instead of R5
        // (R2 and R5 use independent ideal-gas references; no
        // IAPWS R7-97 smoothness condition at the boundary).  One
        // ULP above 1073.15 puts the sample and the runtime bracket
        // both squarely in R5 territory.  Found by code-reviewer
        // adversarial pass.
        const double kT_R5_lo = std::nextafter(1073.15, 2273.15);
        constexpr double kT_R5_hi = 2273.15;
        constexpr double kP_R5_hi = 50.0e6;
        // Same kP_R5_lo as ph_subcritical — 700 Pa is just above
        // IF97's documented Pmin = 611.213 Pa (sat-p at T_min).
        constexpr double kP_R5_lo = 700.0;
        if (kP_R5_lo < kP_R5_hi) {
            const double p_r5_hi = std::min(p_max_sup, kP_R5_hi);
            auto axis = region::AxisTransform::make(region::AxisScale::LOG, kP_R5_lo, p_r5_hi);
            auto t_lo = std::make_unique<region::ConstantCurve>(kP_R5_lo, p_r5_hi, kT_R5_lo);
            auto t_hi = std::make_unique<region::ConstantCurve>(kP_R5_lo, p_r5_hi, kT_R5_hi);
            spec.regions.emplace_back(axis, std::move(t_lo), std::move(t_hi));
        }
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

namespace {
// Shared property list for DmassT_INPUTS: (D, T) inputs, so outputs
// are p, h, s, u, w (speed of sound).  Pressure rides EXP since it
// can span ~10 orders of magnitude across the (D, T) envelope (from
// rarified vapor at low T to compressed liquid at high p_max_eos).
std::vector<PropertySpec> dt_properties(::CoolProp::AbstractState& src) {
    std::vector<PropertySpec> ps = {
      {::CoolProp::iP, svd::OutputTransform::EXP},
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

// Per-knot D upper-bound: at each T knot, find the max pressure ≤
// p_target that HEOS accepts (binary-search down if p_target itself
// is rejected by the melting line at this T) and return D at that
// (T, p).  Result is a CubicSplineCurve representing the "compressed-
// liquid / supercritical-dense ceiling".
//
// Why per-knot rather than fixed isobar: the melting curve cuts into
// the high-p part of the (T, p) envelope for steep-melting fluids
// (CO2 T_melt(0.8 GPa) ≈ 327 K > Tc = 304 K), so a single
// p_max_eos isobar isn't defined at the cold end of LIQUID or SUPER.
// Per-knot walk-down gives the actual reachable D extent — narrower
// at cold T (where melting limits p), wider at hot T.  The resulting
// curve is monotonically non-decreasing in T-region-of-interest and
// always above ρ_sat,L (since p_target > p_sat for liquid-side
// queries).
//
// p_target is the EOS-validity ceiling (≤ 0.99 * pmax); p_floor is
// the minimum acceptable pressure (typically ≥ pc / 2 or sat pressure
// at T to ensure the curve stays above the dome).
// Per-knot floor source: 'sat_L' anchors the walk-down at saturation
// pressure (Q=0 PQ flash), used for LIQUID-side envelopes that must
// stay above the sat dome.  'critical' uses p_critical, used for
// SUPER-side envelopes that should stay supercritical.
enum class WalkFloor : std::uint8_t
{
    SAT_L,
    CRITICAL,
};

std::unique_ptr<region::CubicSplineCurve> build_rho_max_envelope(::CoolProp::AbstractState& heos, double T_min, double T_max, double p_target,
                                                                 WalkFloor floor_source, std::size_t n_knots = 64) {
    if (!(T_min > 0.0) || !(T_max > T_min) || n_knots < 2) {
        throw std::invalid_argument("build_rho_max_envelope: invalid T range or n_knots");
    }
    std::vector<double> T_knots(n_knots);
    std::vector<double> y(n_knots);
    const double dT = (T_max - T_min) / static_cast<double>(n_knots - 1);
    const double p_crit = heos.p_critical();
    constexpr int kMaxHalvings = 50;
    for (std::size_t k = 0; k < n_knots; ++k) {
        T_knots[k] = T_min + static_cast<double>(k) * dT;
        // Per-knot floor: just above sat pressure for LIQUID-side
        // (don't cross the dome into vapor) or just above pc for
        // SUPER-side (stay supercritical).
        double p_floor = 0.0;
        if (floor_source == WalkFloor::SAT_L) {
            try {
                heos.update(::CoolProp::QT_INPUTS, 0.0, T_knots[k]);
                p_floor = heos.p() * 1.001;
            } catch (const ::CoolProp::CoolPropBaseError&) {  // NOLINT(bugprone-empty-catch)
                p_floor = 1.0;                                // fluid without QT support at this T; defer to walk-down
            }
        } else {
            p_floor = 1.001 * p_crit;
        }
        // Find the HIGHEST admissible pressure ≤ p_target at this T (the
        // melting line caps it for steep-melting fluids like CO2).  Walk
        // down from p_target to bracket a reachable pressure, then bisect
        // upward: stopping at the first successful halving would freeze
        // the envelope at ≤ 0.5·p_target and clip valid dense single-
        // phase DT states near the real HEOS limit.
        double p_try = std::max(p_target, p_floor);
        double rho_found = std::nan("");
        double p_ok = std::nan("");   // highest known-reachable p
        double p_bad = std::nan("");  // lowest known-unreachable p (> p_ok)
        for (int i = 0; i < kMaxHalvings; ++i) {
            bool ok = false;
            try {
                heos.update(::CoolProp::PT_INPUTS, p_try, T_knots[k]);
                const double rho = heos.rhomass();
                if (std::isfinite(rho)) {
                    rho_found = rho;
                    p_ok = p_try;
                    ok = true;
                }
            } catch (const ::CoolProp::CoolPropBaseError&) {  // NOLINT(bugprone-empty-catch)
                // Below the melting line at this (T, p); drop p.
            }
            if (ok) break;
            p_bad = p_try;  // this p is unreachable
            if (p_try <= p_floor) break;
            p_try = std::max(0.5 * p_try, p_floor);
        }
        // If p_target itself was unreachable and we backed off to a lower
        // p_ok, bisect (p_ok, p_bad) to recover the pressure ceiling we
        // skipped past, maximizing rho_found.
        if (std::isfinite(p_ok) && std::isfinite(p_bad) && p_bad > p_ok) {
            for (int i = 0; i < kMaxHalvings; ++i) {
                if ((p_bad - p_ok) <= 1e-6 * p_bad) break;
                const double p_mid = 0.5 * (p_ok + p_bad);
                bool ok = false;
                try {
                    heos.update(::CoolProp::PT_INPUTS, p_mid, T_knots[k]);
                    const double rho = heos.rhomass();
                    if (std::isfinite(rho)) {
                        rho_found = rho;
                        p_ok = p_mid;
                        ok = true;
                    }
                } catch (const ::CoolProp::CoolPropBaseError&) {  // NOLINT(bugprone-empty-catch)
                    // (T, p_mid) below the melting line; shrink the bracket.
                }
                if (!ok) p_bad = p_mid;
            }
        }
        if (!std::isfinite(rho_found)) {
            throw std::runtime_error("build_rho_max_envelope: no reachable p ∈ [" + std::to_string(p_floor) + ", " + std::to_string(p_target)
                                     + "] at T=" + std::to_string(T_knots[k]));
        }
        y[k] = rho_found;
    }
    return region::CubicSplineCurve::build(std::move(T_knots), std::move(y));
}

// Per-knot D lower-bound: D at (T, p_target) where p_target is a
// low-pressure isobar.  No walk-down needed here — low p is always
// reachable for the supported T range (vapor side; ideal-gas limit).
// p_target should be ≤ p_triple for safety.
std::unique_ptr<region::CubicSplineCurve> build_rho_min_envelope(::CoolProp::AbstractState& heos, double T_min, double T_max, double p_target,
                                                                 std::size_t n_knots = 64) {
    if (!(T_min > 0.0) || !(T_max > T_min) || n_knots < 2) {
        throw std::invalid_argument("build_rho_min_envelope: invalid T range or n_knots");
    }
    std::vector<double> T_knots(n_knots);
    std::vector<double> y(n_knots);
    const double dT = (T_max - T_min) / static_cast<double>(n_knots - 1);
    for (std::size_t k = 0; k < n_knots; ++k) {
        T_knots[k] = T_min + static_cast<double>(k) * dT;
        try {
            heos.update(::CoolProp::PT_INPUTS, p_target, T_knots[k]);
            y[k] = heos.rhomass();
        } catch (const ::CoolProp::CoolPropBaseError& e) {
            throw std::runtime_error("build_rho_min_envelope: HEOS rejected (T=" + std::to_string(T_knots[k]) + ", p=" + std::to_string(p_target)
                                     + "): " + e.what());
        }
        if (!std::isfinite(y[k])) {
            throw std::runtime_error("build_rho_min_envelope: non-finite rho at (T=" + std::to_string(T_knots[k]) + ", p=" + std::to_string(p_target)
                                     + ")");
        }
    }
    return region::CubicSplineCurve::build(std::move(T_knots), std::move(y));
}
}  // namespace

SurfaceSpec dt_subcritical(::CoolProp::AbstractState& heos, std::size_t NT, std::size_t NR, std::int32_t rank) {
    // DT-indexed surface — closes CoolProp-i7j and supersedes the
    // BICUBIC invert_single_phase_y patch attempt at #2892 / #1301.
    //
    // (D, T) is the Helmholtz EOS's native coordinate, so every output
    // property (p, h, s, u, w) is a direct evaluation at table-build
    // time.  No inversion, no critical-region stiffness, no "cells
    // with -760 MPa" failure modes.
    //
    // Three primary regions: LIQUID + VAPOR (sub-critical, dome-
    // bounded) and SUPER (T > T_critical, EOS-validity bounded).
    // Anomaly handling: when rho_sat,L(T) has interior extrema (water
    // ~277 K, heavy water ~284 K), LIQUID is split into N+1 sub-
    // regions, each spanning one monotonic piece of rho_sat,L(T).
    const double T_triple = heos.Ttriple();
    const double T_min_eos = std::max(T_triple, heos.Tmin());
    const double T_max_eos = heos.Tmax();
    const double Tc = heos.T_critical();

    // Sub-critical T range with margins away from triple / critical.
    const double T_sub_lo = std::max(T_min_eos, T_triple * 1.001);
    const double T_sub_hi = Tc * 0.999;
    // Super-critical T range with margin above critical / below T_max.
    const double T_sup_lo = Tc * 1.001;
    const double T_sup_hi = T_max_eos - 0.5;
    // KNOWN GAP (CoolProp-i7j follow-up): the [0.999, 1.001]·Tc band is
    // intentionally uncovered.  The LIQUID/VAPOR dome boundaries
    // (ρ_sat,L, ρ_sat,V) pinch to ρ_c there, so the dome-bounded
    // sub-regions degenerate to zero width and the SVD η-normalisation
    // is ill-conditioned.  Unlike PT/HmassP, DmassT is NOT yet wired
    // into build_critical_patch_, so single-phase DmassT queries in
    // this ±0.1%·Tc band currently return OutOfRange rather than
    // routing to the source backend.  Coverage impact is ~0.2% of the
    // (D, T) envelope; closing it (a DmassT critical-patch bbox) is
    // filed as a follow-up and is best done alongside the PT/HP patch
    // machinery.

    // p envelope for the non-dome D extents.  Sub-critical fluids
    // have well-defined p_triple from PQ flash at T_triple; use that
    // as the low-pressure floor for the vapor / super low-D side.
    // p_max_target = 0.99 * heos.pmax() is the EOS-validity ceiling
    // — the build_rho_max_envelope helper handles the per-knot
    // walk-down to clear the melting line where it cuts into high-p.
    heos.update(::CoolProp::QT_INPUTS, 0.0, T_triple * 1.001);
    // p_min_eos must stay >= p_triple regardless of region.  Below
    // p_triple HEOS's PT-input path invokes the melting-curve check
    // (the curve's domain starts at p_triple), which throws on
    // every fluid with a polynomial-in-Θ melting curve (CO2, water,
    // most cryogenic fluids).  The thrown error is the same whether
    // T is sub- or super-critical, so SUPER's low-D bound is
    // similarly constrained.  Low-D ideal-gas territory (ρ < ρ_at_
    // (T, p_triple)) is NOT covered by this preset's atlas; queries
    // there should be added in a follow-up that handles the ideal-
    // gas extension (file when needed).
    const double p_min_eos = heos.p() * 1.01;
    const double p_max_target = 0.99 * heos.pmax();

    // Detect rho_sat,L(T) extrema and build a sorted list of LIQUID
    // sub-region breakpoints: [T_sub_lo, extrema..., T_sub_hi].
    std::vector<double> T_breaks = {T_sub_lo};
    for (double T_anom : find_rho_satL_extrema_T(heos, T_sub_lo, T_sub_hi)) {
        T_breaks.push_back(T_anom);
    }
    T_breaks.push_back(T_sub_hi);

    SurfaceSpec spec;
    spec.fluid_name = heos.fluid_names().empty() ? std::string{} : heos.fluid_names().front();
    spec.input_pair = ::CoolProp::DmassT_INPUTS;
    spec.NT = NT;
    spec.NR = NR;
    spec.rank = rank;
    spec.properties = dt_properties(heos);

    // LIQUID sub-region(s).  Primary = T (LINEAR — T span is at most
    // factor ~3 from triple to critical, no orders-of-magnitude
    // variation that would motivate log scaling).  Secondary = D in
    // [rho_sat,L(T), rho_at(T, p_max_eos)].  Each sub-region spans
    // ONE monotonic piece of rho_sat,L(T), so the CubicSplineCurve
    // fits cleanly without the anomaly's hairpin.
    for (std::size_t i = 0; i + 1 < T_breaks.size(); ++i) {
        const double T_lo = T_breaks[i];
        const double T_hi = T_breaks[i + 1];
        if (!(T_hi > T_lo)) continue;
        auto axis = region::AxisTransform::make(region::AxisScale::LINEAR, T_lo, T_hi);
        auto b_lo = build_rho_sat_L(heos, T_lo, T_hi);                                         // dome side
        auto b_hi = build_rho_max_envelope(heos, T_lo, T_hi, p_max_target, WalkFloor::SAT_L);  // compressed-liquid ceiling
        spec.regions.emplace_back(axis, std::move(b_lo), std::move(b_hi));
    }
    // VAPOR.  Secondary = D in [rho_at(T, p_min_eos), rho_sat,V(T)].
    // No anomaly handling needed — rho_sat,V is monotone in T for
    // every known fluid (increases from triple-vapor toward critical
    // density as T → Tc).
    if (T_sub_hi > T_sub_lo) {
        auto axis = region::AxisTransform::make(region::AxisScale::LINEAR, T_sub_lo, T_sub_hi);
        auto b_lo = build_rho_min_envelope(heos, T_sub_lo, T_sub_hi, p_min_eos);
        auto b_hi = build_rho_sat_V(heos, T_sub_lo, T_sub_hi);
        // LOG secondary: rho spans from the near-ideal-gas floor
        // rho(T, p_min_eos) up to rho_sat,V — orders of magnitude at low
        // T — and p ∝ rho there.  Linear-eta strands that tail below the
        // first grid node (CoolProp-wvtz); log-eta resolves it.
        spec.regions.emplace_back(axis, std::move(b_lo), std::move(b_hi), region::AxisScale::LOG);
    }
    // SUPER (T > T_critical).  Secondary D bounded by the two HEOS-
    // validity isobars.  No dome above Tc — single region across the
    // full D extent.  This is the region that fixes yufang67's #1301
    // CO2 supercritical low-density failure: HEOS at low D / high T
    // evaluates directly (P = ρ²(∂α/∂ρ)RT), no inverter, no
    // rectangular bands, no -760 MPa cells.
    if (T_sup_hi > T_sup_lo) {
        auto axis = region::AxisTransform::make(region::AxisScale::LINEAR, T_sup_lo, T_sup_hi);
        auto b_lo = build_rho_min_envelope(heos, T_sup_lo, T_sup_hi, p_min_eos);
        auto b_hi = build_rho_max_envelope(heos, T_sup_lo, T_sup_hi, p_max_target, WalkFloor::CRITICAL);
        // LOG secondary: rho_hi/rho_lo ~ 1e5 here (ideal-gas floor to
        // compressed-supercritical ceiling) with p ∝ rho in the tail.
        // This is the region whose low-rho edge dominated the DT
        // validation error band (CoolProp-wvtz).
        spec.regions.emplace_back(axis, std::move(b_lo), std::move(b_hi), region::AxisScale::LOG);
    }

    spec.update_state = [](::CoolProp::AbstractState& s, double T, double D) {
        // SurfaceSpec passes (a, b) where a is primary (T) and b is
        // secondary (D).  DmassT_INPUTS takes (D, T) in that order.
        s.update(::CoolProp::DmassT_INPUTS, D, T);
    };
    spec.read_property = [](::CoolProp::AbstractState& s, ::CoolProp::parameters key) -> double {
        switch (key) {
            case ::CoolProp::iP:
                return s.p();
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
                throw std::invalid_argument("dt_subcritical: unsupported output property");
        }
    };

    return spec;
}

}  // namespace presets
}  // namespace sbtl
}  // namespace CoolProp
