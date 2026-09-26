#ifndef BOLLENGIERBACKEND_H_
#define BOLLENGIERBACKEND_H_

#include <cmath>
#include <string>
#include <vector>

#include "Backends/Bollengier/BollengierWaterCoefficients.h"
#include "CoolProp/AbstractState.h"
#include "CoolProp/DataStructures.h"
#include "CoolProp/Exceptions.h"
#include "CoolProp/spline/TensorBSpline2D.h"

namespace CoolProp {

/// Liquid water from Bollengier, Brown & Shaw, J. Chem. Phys. 151, 054501
/// (2019) -- a Gibbs energy surface represented as a tensor B-spline in
/// (P, T), doi:10.1063/1.5097179.
///
/// WHY THIS EXISTS.  CoolProp's IAPWS-95 refuses the entire region below the
/// melting line: 250 K/200 MPa, 300 K/2000 MPa and similar all throw.  This
/// backend covers P in [0, 2300.6] MPa, T in [239, 501] K, which includes
/// that region, and is more accurate than IAPWS-95 above ~100 MPa.
///
/// PT_INPUTS ONLY, by design.  For a Gibbs-explicit model (P, T) is the
/// native pair: one surface evaluation plus five partials yields every
/// property with no iteration, where the Helmholtz backends must solve for
/// density.  Serving any other pair would mean adding a Newton inversion,
/// which is deliberately out of scope (see COO-44).  Everything else throws,
/// exactly as IF97Backend does for the pairs it does not implement.
///
/// LIQUID ONLY.  The representation has no vapour branch, no saturation
/// curve and no critical point, so every saturation-flavoured call throws
/// rather than returning a plausible-looking number.
class BollengierBackend : public AbstractState
{
   private:
    /// The Gibbs surface: G [J/kg] as a function of (P [MPa], T [K]).
    /// Built once and shared; TensorBSpline2D is immutable after
    /// construction and its eval() is const, so this is safe to share.
    static const spline::TensorBSpline2D& surface() {
        static const spline::TensorBSpline2D s(std::vector<double>(std::begin(Bollengier::kKnotsP), std::end(Bollengier::kKnotsP)),
                                               std::vector<double>(std::begin(Bollengier::kKnotsT), std::end(Bollengier::kKnotsT)),
                                               Bollengier::kOrderP, Bollengier::kOrderT,
                                               std::vector<double>(std::begin(Bollengier::kCoefs), std::end(Bollengier::kCoefs)));
        return s;
    }

    /// Reference-state offsets (h0 [J/kg], s0 [J/kg/K]).
    ///
    /// The published surface does not use the IAPWS convention, so without
    /// a shift `Water` would silently disagree with itself depending on the
    /// backend.  These are DERIVED from the surface at run time rather than
    /// hard-coded, so they cannot drift away from the coefficients if those
    /// are ever regenerated.
    ///
    /// Applied as G -> G + (h0 - T*s0), which yields h + h0 and s + s0
    /// together.  Shifting h and s independently would break g = h - Ts.
    struct RefShift
    {
        double h0;
        double s0;
    };
    static const RefShift& ref_shift() {
        static const RefShift r = [] {
            // IAPWS-95 reference: u = s = 0 for saturated liquid at the
            // triple point, giving h = p_t * v there.
            constexpr double Tt = 273.16;          // K
            constexpr double pt_MPa = 611.657e-6;  // 611.657 Pa
            constexpr double h_iapws = 0.611781703;
            constexpr double s_iapws = 0.0;
            const double g_raw = surface().eval(pt_MPa, Tt);
            const double s_raw = -surface().eval(pt_MPa, Tt, 0, 1);
            const double h_raw = g_raw + Tt * s_raw;
            return RefShift{h_iapws - h_raw, s_iapws - s_raw};
        }();
        return r;
    }

    /// Published validity domain, stated in the paper.  Written as literals
    /// rather than read back from the knot vectors so that a truncated or
    /// mis-parsed coefficient set cannot quietly redefine what "in range"
    /// means.  dev/scripts/extract_bollengier_coefficients.py checks the
    /// data against these same values at extraction time.
    static constexpr double kPminMPa = 0.0;
    static constexpr double kPmaxMPa = 2300.6;
    static constexpr double kTminK = 239.0;
    static constexpr double kTmaxK = 501.0;

    /// Molar mass of water, IAPWS-95 value [kg/mol].
    static constexpr double kMolarMass = 0.018015268;

    // Cached state, all mass-based; the surface is in J/kg.
    double _p_Pa = _HUGE, _T_K = _HUGE;
    double _v = _HUGE;   ///< m^3/kg
    double _s = _HUGE;   ///< J/kg/K, reference-shifted
    double _h = _HUGE;   ///< J/kg,   reference-shifted
    double _cp = _HUGE;  ///< J/kg/K
    double _cv = _HUGE;  ///< J/kg/K
    double _w = _HUGE;   ///< m/s
    bool _valid = false;

    void require_valid() const {
        if (!_valid) {
            throw ValueError("BollengierBackend: no state has been set; call update(PT_INPUTS, p, T) first");
        }
    }

   public:
    std::string backend_name() override {
        return get_backend_string(BOLLENGIER_BACKEND);
    }
    std::vector<std::string> calc_fluid_names() override {
        return {"Water"};
    }

    // Pure fluid, mass-based like IF97: the surface is J/kg.
    bool using_mole_fractions() override {
        return false;
    }
    bool using_mass_fractions() override {
        return true;
    }
    bool using_volu_fractions() override {
        return false;
    }
    void set_mole_fractions(const std::vector<CoolPropDbl>& /*mole_fractions*/) override {
        throw NotImplementedError("Mole composition is meaningless for a pure fluid backend");
    }
    void set_mass_fractions(const std::vector<CoolPropDbl>& /*mass_fractions*/) override {}
    void set_volu_fractions(const std::vector<CoolPropDbl>& /*volu_fractions*/) override {}
    const std::vector<CoolPropDbl>& get_mole_fractions() override {
        // Mirrors IF97Backend: composition is meaningless here, and
        // throwing is preferable to handing back a fabricated {1.0} that
        // callers might treat as a real mixture specification.
        throw NotImplementedError("Composition has not been implemented for the Bollengier backend");
    }

    void update(CoolProp::input_pairs input_pair, double value1, double value2) override {
        if (input_pair != PT_INPUTS) {
            // Deliberately exhaustive-by-default.  A liquid-only,
            // Gibbs-explicit surface has no honest answer for a density,
            // enthalpy or quality input without an inversion layer, and
            // returning a plausible number would be worse than refusing.
            throw ValueError(format("BollengierBackend supports PT_INPUTS only; got [%s]. See COO-44: the Newton inversion for other "
                                    "pairs is a deliberate follow-up, not an oversight.",
                                    get_input_pair_short_desc(input_pair).c_str()));
        }
        _valid = false;
        const double p_Pa = value1;
        const double T_K = value2;
        const double p_MPa = p_Pa * 1e-6;

        // Range-check before evaluating.  The spline guards its own domain
        // too, but this reports in the caller's units and terms.
        if (!std::isfinite(p_Pa) || !std::isfinite(T_K)) {
            throw ValueError("BollengierBackend: non-finite pressure or temperature");
        }
        if (p_MPa < kPminMPa || p_MPa > kPmaxMPa) {
            throw ValueError(format("BollengierBackend: pressure %g Pa is outside the published domain [%g, %g] MPa", p_Pa, kPminMPa, kPmaxMPa));
        }
        if (T_K < kTminK || T_K > kTmaxK) {
            throw ValueError(format("BollengierBackend: temperature %g K is outside the published domain [%g, %g] K", T_K, kTminK, kTmaxK));
        }

        const spline::TensorBSpline2D& G = surface();
        // Partials of G [J/kg] with respect to P [MPa] and T [K].
        const double G_P = G.eval(p_MPa, T_K, 1, 0);
        const double G_T = G.eval(p_MPa, T_K, 0, 1);
        const double G_PP = G.eval(p_MPa, T_K, 2, 0);
        const double G_TT = G.eval(p_MPa, T_K, 0, 2);
        const double G_PT = G.eval(p_MPa, T_K, 1, 1);
        const double G_val = G.eval(p_MPa, T_K);

        const RefShift& r = ref_shift();
        // v = (dG/dP)_T.  G is J/kg and P is MPa, so the 1e-6 converts to
        // m^3/kg (1 MPa = 1e6 Pa).
        _v = G_P * 1e-6;
        _s = -G_T + r.s0;
        const double g_shifted = G_val + r.h0 - T_K * r.s0;
        _h = g_shifted + T_K * _s;
        _cp = -T_K * G_TT;
        // cv = cp + T (dv/dT)^2 / (dv/dP); the MPa factors cancel exactly.
        _cv = _cp + T_K * G_PT * G_PT / G_PP;
        // (dv/dP)_s = (dv/dP)_T + T (dv/dT)^2 / cp, all in SI.
        const double dvdP_T = G_PP * 1e-12;
        const double dvdT_P = G_PT * 1e-6;
        const double dvdP_s = dvdP_T + T_K * dvdT_P * dvdT_P / _cp;
        _w = std::sqrt(-_v * _v / dvdP_s);

        _p_Pa = p_Pa;
        _T_K = T_K;
        _phase = iphase_liquid;
        _Q = -1;
        _valid = true;
    }

    // ---- state ----------------------------------------------------------
    CoolPropDbl calc_pressure() override {
        require_valid();
        return _p_Pa;
    }
    CoolPropDbl calc_rhomass() override {
        require_valid();
        return 1.0 / _v;
    }
    CoolPropDbl calc_rhomolar() override {
        return calc_rhomass() / kMolarMass;
    }
    CoolPropDbl calc_hmass() override {
        require_valid();
        return _h;
    }
    CoolPropDbl calc_hmolar() override {
        return calc_hmass() * kMolarMass;
    }
    CoolPropDbl calc_smass() override {
        require_valid();
        return _s;
    }
    CoolPropDbl calc_smolar() override {
        return calc_smass() * kMolarMass;
    }
    CoolPropDbl calc_umass() override {
        require_valid();
        return _h - _p_Pa * _v;
    }
    CoolPropDbl calc_umolar() override {
        return calc_umass() * kMolarMass;
    }
    CoolPropDbl calc_cpmass() override {
        require_valid();
        return _cp;
    }
    CoolPropDbl calc_cpmolar() override {
        return calc_cpmass() * kMolarMass;
    }
    CoolPropDbl calc_cvmass() override {
        require_valid();
        return _cv;
    }
    CoolPropDbl calc_cvmolar() override {
        return calc_cvmass() * kMolarMass;
    }
    CoolPropDbl calc_speed_sound() override {
        require_valid();
        return _w;
    }
    phases calc_phase() override {
        require_valid();
        return iphase_liquid;
    }

    // ---- constants ------------------------------------------------------
    CoolPropDbl calc_molar_mass() override {
        return kMolarMass;
    }
    CoolPropDbl calc_gas_constant() override {
        return 8.31446261815324;
    }
    CoolPropDbl calc_Tmin() override {
        return kTminK;
    }
    CoolPropDbl calc_Tmax() override {
        return kTmaxK;
    }
    CoolPropDbl calc_pmax() override {
        return kPmaxMPa * 1e6;
    }
    CoolPropDbl calc_Ttriple() override {
        return 273.16;
    }

    // ---- deliberately absent -------------------------------------------
    // No critical point, no saturation, no melting line: this is a
    // single-phase liquid surface.  AbstractState's defaults already throw
    // NotImplementedError for those, which is the honest answer, so they
    // are not overridden here to return a fabricated value.
};

}  // namespace CoolProp

#endif  // BOLLENGIERBACKEND_H_
