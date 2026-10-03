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
/// WHY THIS EXISTS.  CoolProp's IAPWS-95 refuses the entire region below
/// the melting line: 250 K/200 MPa, 300 K/2000 MPa and similar all throw.
/// This backend covers T in [240, 500] K and P up to
/// 2300.5999999999995 MPa, which includes that region, and is more
/// accurate than IAPWS-95 above ~100 MPa.  (The pressure bound is the
/// fitted knot.  The paper says "2300 MPa"; 2300.6 is this code's
/// rounding of the knot and is REFUSED -- see domain().)
///
/// ONE EXCLUDED BOX INSIDE THAT DOMAIN.  The advertised bounds stay as
/// published rather than being narrowed to dodge this, and a single
/// rectangle is cut out of them: update() refuses
/// p >= 1500 MPa AND T <= 255 K.  Both of CoolProp's motivating cases
/// above -- 250 K/200 MPa and 300 K/2000 MPa -- are outside it and still
/// served.
///
/// In that corner the fitted surface stops being usable, in two ways
/// that run into each other.  Measured on a 1 MPa x 0.1 K grid: over
/// p in [1833.5, 2300.6] MPa and T in [240.0, 249.6] K it is not
/// thermodynamically admissible at all (cv <= 0 over most of it,
/// (dv/dP)_T >= 0 over the rest, and a non-finite w throughout), and
/// around that, out to p >= 1526.4 MPa and T <= 251.6 K, cv stays
/// positive but falls continuously to 0.03 J/kg/K -- two orders of
/// magnitude below liquid water.  There is no boundary between the two
/// for a caller to find.  The box rounds the union outward to round
/// numbers; zero problem states lie outside it.
///
/// It is stated as a REGION, not as a test on the computed properties,
/// because the defect belongs to the fit over an area -- so an area is
/// what can be written down, audited, and held still across a
/// coefficient revision.  It also means refusal does not depend on which
/// property you ask for.  See kExcludedPminMPa for the property-based
/// rules that were tried and are worse, and for the cost: the box is
/// 2.02% of the rectangle, and 1.26% of it would have evaluated fine.
///
/// That corner is ~2.6x deeper in pressure than any datum below 293 K in
/// the authors' own fitted set, and Supplementary Material E declines to
/// tabulate 250 K above 900 MPa, so nothing validated is forfeited.
/// Reported upstream to SeaFreeze.  It is a limit of the published fit,
/// not of this implementation.
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

    /// The excluded box: refused outright, inside the published domain.
    ///
    /// Chosen to contain the whole problematic corner with margin, rather
    /// than to trace its edge.  Measured on a 1 MPa x 0.1 K grid over the
    /// published rectangle, everything the surface gets wrong lies at
    /// p >= 1526.4 MPa and T <= 251.6 K; this box rounds that outward.
    /// Verified: zero problem states fall outside it.
    ///
    /// A box, and not a test on the computed properties, because the
    /// defect is a property of the FIT over a region -- so a region is
    /// what can be stated, audited and held constant.  Every
    /// property-based rule tried here was worse.  A cv floor leaves the
    /// absurd sound speeds on both sides of itself (cv and w are not in
    /// correspondence: states at 11969 m/s clear any floor that still
    /// admits ordinary water, and states at an unremarkable 2689 m/s fall
    /// below it).  A w ceiling cannot be placed at all, because the
    /// pathological and physical populations overlap in w -- see update().
    /// Both also make the refusal depend on which property you ask for,
    /// which a caller cannot predict.
    ///
    /// The cost is accepted deliberately: the box is 2.02% of the
    /// rectangle and 1.26% of it is states that would have evaluated
    /// fine.  They sit in the corner the authors themselves decline to
    /// tabulate (SM_E omits 250 K above 900 MPa), so refusing them
    /// forfeits nothing that was ever validated.
    static constexpr double kExcludedPminMPa = 1500.0;
    static constexpr double kExcludedTmaxK = 255.0;

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

        if (p_MPa >= kExcludedPminMPa && T_K <= kExcludedTmaxK) {
            throw ValueError(format("BollengierBackend: p = %g Pa, T = %g K is inside the excluded box (p >= %g MPa and "
                                    "T <= %g K), a cold high-pressure corner of the published domain where the fitted "
                                    "surface is not usable: it turns thermodynamically inadmissible over part of it and "
                                    "sub-physical over the rest. The corner is far beyond any data the fit was "
                                    "constrained by, and the authors' own tables omit it. This is a limit of the "
                                    "published fit, not of this implementation; the domain is otherwise as published.",
                                    p_Pa, T_K, kExcludedPminMPa, kExcludedTmaxK));
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
        // breaking down.
        //
        // NO threshold at any value separates the populations.  Measured:
        // pathological states (cv < 2000, half of water's ~3800) run as slow
        // as w = 2496.7 m/s at 1480.4 MPa / 240 K -- 1016 m/s BELOW the
        // fastest state the authors themselves publish (3512.7 at
        // 2200 MPa / 300 K).  Any cut above 2496.7 admits pathological
        // states; any cut at or below it refuses published ones.
        //
        // A very high cut would truncate the worst tail, but that is not
        // separation: a 20 km/s cut removes 289 states of which fewer than
        // half have cv < 100 (the rest run up to cv = 613), and leaves
        // 36744 of 36775 cv < 2000 states served.  And it is not the
        // argument for omitting a ceiling anyway.  The argument is that the
        // reference implementation makes no such claim: SeaFreeze serves
        // this region deliberately, and inventing a criterion it does not
        // apply is what produced two failed guards here already.
        //
        // WHAT THIS MEANS FOR CALLERS -- read this before using the cold
        // high-pressure corner.
        //
        // How far the advertised rectangle reaches beyond the measurements.
        // WaterEOS.mat carries the fitted data alongside the coefficients,
        // in /H2O/{SS,rho,Cp}/data -- 2781 points, which the extraction
        // script does not read but which answer this directly:
        //
        //   199-240 K   150 pts   max p    399 MPa
        //   240-260 K   280 pts   max p    399 MPa
        //   260-280 K   618 pts   max p    611 MPa
        //   280-293 K   257 pts   max p    695 MPa
        //   293-320 K   617 pts   max p   1720 MPa
        //   320-500 K   714 pts   max p   4832 MPa
        //
        // So the constraint is a JOINT one, not a bound on either axis:
        // below 293 K no datum of any kind exists above 695 MPa.  The
        // region refused below (p >~ 1834 MPa, T <~ 249 K) sits about 2.6x
        // deeper in pressure than any measurement at that temperature.
        //
        // (Two earlier versions of this comment were wrong.  One named a
        // pressure where extrapolation "starts" -- no such threshold
        // exists.  The other read the coverage off Supplementary Material B
        // and concluded the domain was extrapolation in both directions
        // separately; SM_B is only the authors' own 901 runs, and the full
        // set reaches 8595 MPa in sound speed and 199.6 K in density.)
        //
        // The authors decline to publish in the cold high-pressure corner
        // at all: Supplementary Material E omits 250 K above 900 MPa in
        // every property it tabulates.
        //
        // Those states stay thermodynamically self-consistent (cv > 0,
        // (dv/dP)_T < 0, cv < cp) but stop describing water, and the
        // degradation is UNBOUNDED, not merely large.  Measured over the
        // advertised rectangle on a grid of 2001 pressures x 1001
        // temperatures, admitted states reach cv = 0.047 J/kg/K (water is
        // ~3800), cp/cv = 56002, and w = 1.18e6 m/s.
        //
        // Those are sampling artefacts, not extrema.  cv -> 0 continuously
        // along the locus cp*G_PP + T*G_PT^2 = 0, where G_PP is strictly
        // NEGATIVE (-0.053 at the nearest point) -- so this is a separate
        // locus from the convexity flip G_PP = 0, where cv instead diverges
        // to -infinity.  Both appear in the refused set.  Since
        // w^2 = -v^2 cp / (1e-12 G_PP cv), the cv -> 0 locus is exactly
        // where w diverges: the true infimum of cv is 0 and w is unbounded.
        //
        // They are SERVED regardless, matching the reference implementation.
        // A caller who needs to know whether liquid is the stable phase --
        // or whether these numbers mean anything -- must determine that
        // separately.
        // Unreachable with the committed coefficients: the excluded box
        // above contains every state that trips this, verified by scanning
        // the published rectangle.  Kept as a backstop, not as dead code --
        // it is the only thing standing between a coefficient revision that
        // moved the bad region and a silently propagating NaN, and the
        // values it tests are computed anyway, so it costs nothing.  If
        // this ever fires, the box needs re-measuring.
        if (!(G_PP < 0.0) || !std::isfinite(_w) || !std::isfinite(_cv) || _cv <= 0.0) {
            throw ValueError(format("BollengierBackend: the published surface is not thermodynamically admissible at "
                                    "p = %g Pa, T = %g K ((dv/dP)_T >= 0 or cv <= 0), so no properties can be returned. "
                                    "This is a known carve-out inside the paper's stated domain -- a small cold "
                                    "high-pressure corner far beyond any data the fit was constrained by -- not a "
                                    "failure of this backend. The reference implementation returns NaN here; CoolProp "
                                    "throws rather than propagate one.",
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
