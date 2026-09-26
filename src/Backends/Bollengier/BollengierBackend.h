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

    /// NOTE ON THE REFERENCE STATE.  This backend reports the published
    /// surface as-is.  Bollengier et al. do not use the IAPWS convention:
    /// at the triple-point saturated-liquid state (273.16 K, 611.657 Pa)
    /// this model gives h = 71.2317874 J/kg and s = 0.2581766650 J/kg/K, where
    /// IAPWS-95 has 0.6117817 and 0.
    ///
    /// No shift is applied, deliberately.  h and s are defined only up to
    /// a constant; silently re-anchoring them would mean CoolProp reporting
    /// numbers that are not the published model's, so anyone checking
    /// against the paper's own tables would see an unexplained ~70 J/kg
    /// offset.  Reference state is a user-level concern in CoolProp, with
    /// set_reference_stateS() / set_reference_stateD() as the supported
    /// mechanism -- though note that as of writing neither reaches this
    /// backend: set_reference_stateS() silently no-ops on an unrecognised
    /// backend prefix and set_reference_stateD() throws, so there is
    /// currently NO supported way to re-anchor it (tracked separately).
    /// That is an argument for fixing that plumbing, not for baking a datum
    /// in here.
    ///
    /// Note also that re-anchoring to IAPWS would NOT buy cross-backend
    /// agreement: it forces equality at the reference point and moves
    /// everything else further away.  Measured at 300 K, 0.1 MPa the
    /// unshifted model differs from IAPWS-95 by -4.25 J/kg in h; shifted it
    /// differs by -74.9.  The residual is a genuine cp difference in the cold
    /// liquid region, not a datum mismatch.

    /// Published validity domain, as stated in the paper.  These are
    /// LITERALS, deliberately -- they are used to validate the coefficient
    /// data, not derived from it, so a truncated or mis-parsed set cannot
    /// quietly redefine what "in range" means.
    static constexpr double kPaperPminMPa = 0.0;
    static constexpr double kPaperPmaxMPa = 2300.6;
    static constexpr double kPaperTminK = 240.0;
    static constexpr double kPaperTmaxK = 500.0;

    /// The advertised domain, and the check that the data supports it.
    ///
    /// TEMPERATURE comes from the PAPER: 240-500 K, as stated in its title
    /// and abstract.  The fitted knots run slightly wider
    /// (238.99999999999991 to 501.00000000000023), but a knot span is a
    /// property of the fit, not a claim of validity, and advertising it
    /// would promise a degree of coverage the authors never did.
    ///
    /// PRESSURE comes from the surface's own support, because the paper's
    /// "2300 MPa" is a rounded statement of a fitted bound that is actually
    /// 2300.5999999999995.  Guarding on the rounded literal would accept
    /// p = 2300.6 MPa and then throw from INSIDE the spline -- a worse
    /// message, from a guard that disagrees with the thing it guards.
    ///
    /// Construction checks that the knots COVER the advertised range
    /// (rather than equal it), so a truncated or wrong-revision coefficient
    /// set still fails loudly, while a legitimate ULP difference does not.
    ///
    /// WHERE THE AUTHORS DECLINE TO REPORT.  Supplementary Material E
    /// tabulates rho, cp and w on a 250-500 K x 0.1-2200 MPa grid, and
    /// OMITS 250 K above 900 MPa in all three properties -- the cold
    /// high-pressure corner, deep in the ice VI/VII field.  The surface is
    /// still evaluable there and this backend serves it, matching the
    /// reference implementation, but the authors publish no values to check
    /// against and the fit is unconstrained; see update().
    struct Domain
    {
        double p_min, p_max, T_min, T_max;
    };
    static const Domain& domain() {
        static const Domain d = [] {
            const spline::TensorBSpline2D& s = surface();
            (void)s;  // built for its validation side effect
            const double kp_lo = Bollengier::kKnotsP[Bollengier::kOrderP - 1];
            const double kp_hi = Bollengier::kKnotsP[Bollengier::kNP];
            const double kt_lo = Bollengier::kKnotsT[Bollengier::kOrderT - 1];
            const double kt_hi = Bollengier::kKnotsT[Bollengier::kNT];
            auto near = [](double a, double b) { return std::fabs(a - b) <= 1e-9 * std::fmax(1.0, std::fabs(b)); };
            // p: the support IS the advertised bound, to within ULP noise.
            if (!near(kp_lo, kPaperPminMPa) || !near(kp_hi, kPaperPmaxMPa)) {
                throw ValueError(format("BollengierBackend: the loaded coefficients span p [%.17g, %.17g] MPa, which does not "
                                        "match the published [%g, %g]",
                                        kp_lo, kp_hi, kPaperPminMPa, kPaperPmaxMPa));
            }
            // T: the support must COVER the paper's range, not equal it.
            if (kt_lo > kPaperTminK || kt_hi < kPaperTmaxK) {
                throw ValueError(format("BollengierBackend: the loaded coefficients span T [%.17g, %.17g] K, which does not "
                                        "cover the published [%g, %g]",
                                        kt_lo, kt_hi, kPaperTminK, kPaperTmaxK));
            }
            return Domain{kp_lo, kp_hi, kPaperTminK, kPaperTmaxK};
        }();
        return d;
    }

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
        // Invalidate FIRST, before anything that can throw.  Two separate
        // hazards, both of which returned stale values from a previous
        // successful update():
        //
        //  - AbstractState caches speed_sound, hmolar, smolar, cpmolar and
        //    friends in its own CacheArray.  Without clear() those freeze at
        //    the first value ever computed, and it reaches PropsSI --
        //    successive queries in one process return the PREVIOUS query's
        //    answer.  IF97Backend::clear() carries a comment describing this
        //    exact symptom ("a constant speed_sound surface").
        //  - _valid must be cleared ahead of the input-pair check below, not
        //    after it: a rejected pair otherwise throws while leaving the
        //    old state readable.
        clear();
        _valid = false;

        if (input_pair != PT_INPUTS) {
            // Deliberately exhaustive-by-default.  A liquid-only,
            // Gibbs-explicit surface has no honest answer for a density,
            // enthalpy or quality input without an inversion layer, and
            // returning a plausible number would be worse than refusing.
            throw ValueError(format("BollengierBackend supports PT_INPUTS only; got [%s]. See COO-44: the Newton inversion for other "
                                    "pairs is a deliberate follow-up, not an oversight.",
                                    get_input_pair_short_desc(input_pair).c_str()));
        }
        const double p_Pa = value1;
        const double T_K = value2;
        const double p_MPa = p_Pa * 1e-6;

        // Range-check before evaluating.  The spline guards its own domain
        // too, but this reports in the caller's units and terms.
        if (!std::isfinite(p_Pa) || !std::isfinite(T_K)) {
            throw ValueError("BollengierBackend: non-finite pressure or temperature");
        }
        // Tested on p_Pa, not p_MPa: a tiny negative pressure underflows to
        // -0.0 under the 1e-6 conversion, and -0.0 < 0.0 is false, so the
        // range check below would admit it.
        if (p_Pa < 0.0) {
            throw ValueError(format("BollengierBackend: negative pressure %g Pa", p_Pa));
        }
        const Domain& dom = domain();
        if (p_MPa < dom.p_min || p_MPa > dom.p_max) {
            throw ValueError(
              format("BollengierBackend: pressure %g Pa is outside the published domain [%g, %g] MPa", p_Pa, kPaperPminMPa, kPaperPmaxMPa));
        }
        if (T_K < dom.T_min || T_K > dom.T_max) {
            throw ValueError(format("BollengierBackend: temperature %g K is outside the published domain [%g, %g] K", T_K, kPaperTminK, kPaperTmaxK));
        }

        const spline::TensorBSpline2D& G = surface();
        // Partials of G [J/kg] with respect to P [MPa] and T [K].
        const double G_P = G.eval(p_MPa, T_K, 1, 0);
        const double G_T = G.eval(p_MPa, T_K, 0, 1);
        const double G_PP = G.eval(p_MPa, T_K, 2, 0);
        const double G_TT = G.eval(p_MPa, T_K, 0, 2);
        const double G_PT = G.eval(p_MPa, T_K, 1, 1);
        const double G_val = G.eval(p_MPa, T_K);

        // v = (dG/dP)_T.  G is J/kg and P is MPa, so the 1e-6 converts to
        // m^3/kg (1 MPa = 1e6 Pa).
        _v = G_P * 1e-6;
        _s = -G_T;
        _h = G_val + T_K * _s;
        _cp = -T_K * G_TT;
        // cv = cp + T (dv/dT)^2 / (dv/dP); the MPa factors cancel exactly.
        _cv = _cp + T_K * G_PT * G_PT / G_PP;
        // (dv/dP)_s = (dv/dP)_T + T (dv/dT)^2 / cp, all in SI.
        const double dvdP_T = G_PP * 1e-12;
        const double dvdT_P = G_PT * 1e-6;
        const double dvdP_s = dvdP_T + T_K * dvdT_P * dvdT_P / _cp;
        _w = std::sqrt(-_v * _v / dvdP_s);

        // Refuse states the representation cannot evaluate, and ONLY
        // those.
        //
        // The reference implementation (SeaFreeze, by a co-author of the
        // paper) does no range-checking and no output sanity-checking at
        // all: getProp() evaluates whatever phase it is asked for, returns
        // NaN outside the parametrisation, and exposes a SEPARATE
        // whichphase() for callers who need to know which phase is actually
        // stable.  Evaluating the liquid surface inside the ice field is a
        // deliberate capability there, not an error -- serving the
        // metastable region is much of the point.
        //
        // So this guard matches that contract on the model and keeps
        // CoolProp's on safety: what SeaFreeze returns as NaN, this throws,
        // because a non-finite property silently propagating is the failure
        // mode TensorBSpline2D exists to prevent.  Nothing more.
        //
        // In particular there is NO sound-speed ceiling.  An earlier version
        // had one, on the theory that absurd stiffness marked the fit
        // breaking down.  It cannot work: measured, an unphysical state at
        // 1661.5 MPa / 239 K has w = 3368 m/s while a perfectly physical one
        // at 2300.6 MPa / 300 K has w = 3587.  The populations overlap in w,
        // so no threshold separates them, and the reference implementation
        // makes no such claim in the first place.
        //
        // WHAT THIS MEANS FOR CALLERS.  In the cold high-pressure corner
        // (p >~ 1660 MPa, T <~ 260 K) the surface is an unconstrained
        // extrapolation -- deep inside the ice VI/VII field, where no liquid
        // data exists to have fitted against.  It stays thermodynamically
        // self-consistent (cv > 0, (dv/dP)_T < 0) but stops describing
        // water: cv falls as low as ~416 J/kg/K against roughly 3800 for
        // real water, with cp/cv ~ 5.4.  Those states are SERVED, as the
        // reference serves them.  A caller who needs to know whether the
        // liquid is the stable phase there must ask separately.
        if (!(G_PP < 0.0) || !std::isfinite(_w) || !std::isfinite(_cv) || _cv <= 0.0) {
            throw ValueError(format("BollengierBackend: the representation is not evaluable at p = %g Pa, T = %g K "
                                    "((dv/dP)_T >= 0, or a non-finite result). The reference implementation returns NaN "
                                    "here; CoolProp throws rather than propagating one.",
                                    p_Pa, T_K));
        }

        _p_Pa = p_Pa;
        _T_K = T_K;
        // B2: the base class's _p/_T back T(), p() and keyed_output(); without
        // these they stay at -_HUGE from clear() and every consumer that does
        // not go through PropsSI's input short-circuit sees -inf.
        _p = p_Pa;
        _T = T_K;
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
        return domain().T_min;
    }
    CoolPropDbl calc_Tmax() override {
        return domain().T_max;
    }
    CoolPropDbl calc_pmax() override {
        return domain().p_max * 1e6;
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
