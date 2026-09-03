/*
Regression guards for the Peneloux volume translation in the cubic backends.

A temperature-independent volume translation v -> v - c is a change of variable, not a new model, so
almost every property either stays *exactly* the same or moves by an amount that is known in closed
form.  Jaubert, Privat, Le Guennec & Coniglio, Fluid Phase Equilibria 419 (2016) 88-95 tabulate
which is which; those identities are exact, which makes them unusually good regression guards.

Before this file there was no coverage of `cm` anywhere in the suite, and four code paths that
reason directly about volume silently ignored it: the Chebyshev superancillary (which made
update(DmolarT) hand back a NEGATIVE pressure for a state inside the dome), the maximum-density
bound, the spinodal, and the critical density.  `cm` also did not survive get_copy().

Run with:  ./CatchTestRunner "[volume_translation]"
*/

#if defined(ENABLE_CATCH)
#    include <catch2/catch_all.hpp>

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/DataStructures.h"
#    include "../Backends/Cubics/CubicBackend.h"
#    include "../Backends/Cubics/GeneralizedCubic.h"

#    include <cmath>
#    include <functional>
#    include <memory>
#    include <string>
#    include <utility>
#    include <vector>

using namespace CoolProp;

namespace {

using ASptr = std::shared_ptr<AbstractState>;

/// The backends that carry a superancillary and are exercised throughout this file.
const std::vector<std::string>& cubic_backends() {
    static const std::vector<std::string> v{"SRK", "PR"};
    return v;
}

/// A small, chemically varied set that all have a reference EOS available.
const std::vector<std::string>& probe_fluids() {
    static const std::vector<std::string> v{"n-Propane", "Water", "CarbonDioxide"};
    return v;
}

AbstractCubicBackend& as_cubic(AbstractState& AS) {
    auto* cb = dynamic_cast<AbstractCubicBackend*>(&AS);
    REQUIRE(cb != nullptr);
    return *cb;
}

/// The pure-fluid covolume, which sets the scale for a safe translation: the repulsive pole sits at
/// v = b - c, so |c| is always chosen as a modest fraction of b rather than as an absolute number
/// that might be fine for propane and past the pole for water.
double covolume(const std::string& backend, const std::string& fluid) {
    ASptr AS(AbstractState::factory(backend, fluid));
    return as_cubic(*AS).get_cubic()->bm_term(std::vector<double>(1, 1.0));
}

ASptr make_translated(const std::string& backend, const std::string& fluid, double c) {
    ASptr AS(AbstractState::factory(backend, fluid));
    if (c != 0.0) {
        AS->set_fluid_parameter_double(0, "cm", c);
    }
    return AS;
}

/// |got - ref| against a relative tolerance with an absolute floor.  Energies and entropies carry an
/// arbitrary reference state and can legitimately sit near zero, where a pure relative test degrades
/// into a division by nothing.
bool close_to(double got, double ref, double rtol, double afloor) {
    return std::abs(got - ref) <= rtol * std::abs(ref) + afloor;
}

struct Probe
{
    const char* name;
    std::function<double(AbstractState&)> get;
    /// Expected value minus untranslated value.  Invariant properties return 0.
    std::function<double(double p, double T, double c)> shift;
    double afloor;  ///< absolute tolerance floor, in the property's own units
};

/// The invariance table of Jaubert et al. (2016), as executable assertions.
const std::vector<Probe>& pure_fluid_probes() {
    static const double R = 8.31446261815324;
    static const std::vector<Probe> v{
      // -- properties that move by a known amount ------------------------------------------------
      {"v", [](AbstractState& S) { return 1.0 / S.rhomolar(); }, [](double, double, double c) { return -c; }, 1e-18},
      {"h", [](AbstractState& S) { return S.hmolar(); }, [](double p, double, double c) { return -p * c; }, 1e-9},
      {"g", [](AbstractState& S) { return S.gibbsmolar(); }, [](double p, double, double c) { return -p * c; }, 1e-9},
      {"ln(phi)", [](AbstractState& S) { return std::log(S.fugacity_coefficient(0)); }, [](double p, double T, double c) { return -p * c / (R * T); },
       1e-12},
      // -- properties that must not move at all --------------------------------------------------
      {"s", [](AbstractState& S) { return S.smolar(); }, [](double, double, double) { return 0.0; }, 1e-10},
      {"u", [](AbstractState& S) { return S.umolar(); }, [](double, double, double) { return 0.0; }, 1e-8},
      {"a", [](AbstractState& S) { return S.helmholtzmolar(); }, [](double, double, double) { return 0.0; }, 1e-8},
      {"cp", [](AbstractState& S) { return S.cpmolar(); }, [](double, double, double) { return 0.0; }, 1e-10},
      {"cv", [](AbstractState& S) { return S.cvmolar(); }, [](double, double, double) { return 0.0; }, 1e-10},
      // Jaubert et al. state these two as the PRODUCTS kappa_T*v and alpha_v*v rather than as the
      // normalised compressibility and expansivity, and that is the sharp form: the derivatives
      // -(dv/dP)_T and (dv/dT)_P are invariant, while dividing either by v does NOT survive,
      // because v is exactly what moved.  Both are asserted -- invariance here, and the
      // normalised pair is checked to genuinely change at the call site.
      {"-(dv/dP)_T", [](AbstractState& S) { return S.isothermal_compressibility() / S.rhomolar(); }, [](double, double, double) { return 0.0; },
       1e-25},
      {"(dv/dT)_P", [](AbstractState& S) { return S.isobaric_expansion_coefficient() / S.rhomolar(); }, [](double, double, double) { return 0.0; },
       1e-20},
    };
    return v;
}

struct StatePoint
{
    std::string label;
    double T, p;
};

/// The quantities that move because v enters them explicitly, asserted as the exact identities
/// they are rather than merely as "different".
///
/// Each shift is known in closed form, so there is no reason to settle for an inequality: asserting
/// only that one of these changed accepts a change of any wrong magnitude, which is what the
/// weaker checks alongside these calls do on their own.
void check_v_scaled_quantities(AbstractState& ref, AbstractState& trn, double c) {
    const double v_ref = 1.0 / ref.rhomolar(), v_trn = 1.0 / trn.rhomolar();

    // w^2 = -v^2*(dp/dv)_s/M with (dp/dv)_s invariant, so w scales with v EXACTLY.  Asserting only
    // that w changed would accept a change of any wrong magnitude.
    CHECK(close_to(trn.speed_sound() / ref.speed_sound(), v_trn / v_ref, 1e-9, 0.0));

    // kappa_T and alpha_v are the invariant derivatives divided by v, so they scale as 1/v.  The
    // documented rows are the bare derivatives; these are the normalised forms the public API
    // actually returns, so both are worth pinning.
    CHECK(close_to(trn.isothermal_compressibility() / ref.isothermal_compressibility(), v_ref / v_trn, 1e-9, 0.0));
    CHECK(close_to(trn.isobaric_expansion_coefficient() / ref.isobaric_expansion_coefficient(), v_ref / v_trn, 1e-9, 0.0));

    // mu_JT = (dT/dp)_h shifts by +c/c_p.
    const double mu_ref = ref.first_partial_deriv(iT, iP, iHmolar);
    const double mu_trn = trn.first_partial_deriv(iT, iP, iHmolar);
    CAPTURE(mu_ref);
    CAPTURE(mu_trn);
    CHECK(close_to(mu_trn, mu_ref + c / ref.cpmolar(), 1e-8, 1e-18));

    // Z = p*v/(RT) inherits the volume shift directly.
    CHECK(close_to(trn.compressibility_factor(), ref.compressibility_factor() - ref.p() * c / (8.31446261815324 * ref.T()), 1e-9, 1e-14));
}

/// Derived rather than hardcoded: a fixed (T, p) that is liquid for propane can be supercritical for
/// carbon dioxide, and a state that lands inside the dome makes the PT flash pick a phase rather
/// than a root.
std::vector<StatePoint> single_phase_points(AbstractState& AS) {
    const double Tc = AS.T_critical(), pc = AS.p_critical();
    const double T_sub = 0.75 * Tc;
    AS.update(QT_INPUTS, 0, T_sub);
    const double psat = AS.p();
    return {
      {"compressed liquid", T_sub, 3.0 * psat},
      {"superheated vapour", T_sub, 0.4 * psat},
      {"supercritical", 1.25 * Tc, 2.0 * pc},
    };
}

}  // namespace

TEST_CASE("Cubic volume translation: pure-fluid property invariances", "[cubic][volume_translation]") {
    for (const auto& backend : cubic_backends()) {
        for (const auto& fluid : probe_fluids()) {
            CAPTURE(backend);
            CAPTURE(fluid);
            const double b = covolume(backend, fluid);

            // Both signs matter: a positive c shrinks the volume and moves the repulsive pole to
            // higher density, a negative c does the reverse.  Both occur in the real parameter sets.
            for (double frac : {0.15, -0.15}) {
                const double c = frac * b;
                CAPTURE(c);

                ASptr ref = make_translated(backend, fluid, 0.0);
                ASptr trn = make_translated(backend, fluid, c);

                for (const auto& pt : single_phase_points(*ref)) {
                    CAPTURE(pt.label);
                    CAPTURE(pt.T);
                    CAPTURE(pt.p);
                    REQUIRE_NOTHROW(ref->update(PT_INPUTS, pt.p, pt.T));
                    REQUIRE_NOTHROW(trn->update(PT_INPUTS, pt.p, pt.T));

                    for (const auto& probe : pure_fluid_probes()) {
                        CAPTURE(probe.name);
                        const double v0 = probe.get(*ref);
                        const double v1 = probe.get(*trn);
                        const double expected = v0 + probe.shift(pt.p, pt.T, c);
                        CAPTURE(v0);
                        CAPTURE(v1);
                        CAPTURE(expected);
                        // 1e-9 relative is loose enough for the two independent density solves and
                        // far tighter than any of the defects this guards, which are O(1) relative.
                        CHECK(close_to(v1, expected, 1e-9, probe.afloor));
                    }

                    // Quantities that must genuinely MOVE.  Without these a suite made only of
                    // invariances would pass with the translation silently disabled.
                    //
                    // w^2 scales with v^2; the normalised compressibility and expansivity are the
                    // invariant derivatives above divided by v, so they inherit the shift even
                    // though their numerators do not.
                    CHECK(std::abs(trn->speed_sound() / ref->speed_sound() - 1.0) > 1e-6);
                    CHECK(std::abs(trn->isothermal_compressibility() / ref->isothermal_compressibility() - 1.0) > 1e-6);
                    CHECK(std::abs(trn->isobaric_expansion_coefficient() / ref->isobaric_expansion_coefficient() - 1.0) > 1e-6);
                    check_v_scaled_quantities(*ref, *trn, c);
                }
            }
        }
    }
}

TEST_CASE("Cubic volume translation: second virial coefficient shifts by -c", "[cubic][volume_translation]") {
    for (const auto& backend : cubic_backends()) {
        for (const auto& fluid : probe_fluids()) {
            const double b = covolume(backend, fluid);
            // Both signs, for the same reason the single-phase suite sweeps both: the two are not
            // mirror images through the code, and both occur in published parameter sets.
            for (double frac : {0.15, -0.15}) {
                CAPTURE(backend);
                CAPTURE(fluid);
                CAPTURE(frac);
                const double c = frac * b;
                ASptr ref = make_translated(backend, fluid, 0.0);
                ASptr trn = make_translated(backend, fluid, c);
                const double T = 1.25 * ref->T_critical();
                ref->update(DmolarT_INPUTS, 1e-3, T);
                trn->update(DmolarT_INPUTS, 1e-3, T);
                CHECK(close_to(trn->Bvirial(), ref->Bvirial() - c, 1e-9, 1e-18));
            }
        }
    }
}

TEST_CASE("Cubic volume translation: saturation invariances", "[cubic][volume_translation]") {
    for (const auto& backend : cubic_backends()) {
        for (const auto& fluid : probe_fluids()) {
            const double b = covolume(backend, fluid);
            // Both signs, and three reduced temperatures rather than one.  The property changes on
            // vaporisation shrink as the two phases converge towards the critical point, so a
            // single mid-dome temperature says nothing about whether the identities still hold
            // where Dvap V and Dvap H are small and the cancellation is correspondingly delicate.
            for (double frac : {0.15, -0.15}) {
                for (double Tr : {0.65, 0.75, 0.85}) {
                    CAPTURE(backend);
                    CAPTURE(fluid);
                    CAPTURE(frac);
                    CAPTURE(Tr);
                    const double c = frac * b;
                    ASptr ref = make_translated(backend, fluid, 0.0);
                    ASptr trn = make_translated(backend, fluid, c);

                    const double T = Tr * ref->T_critical();
                    REQUIRE_NOTHROW(ref->update(QT_INPUTS, 0, T));
                    REQUIRE_NOTHROW(trn->update(QT_INPUTS, 0, T));

                    // The vapour pressure is the headline invariance: it is what lets a translation
                    // be fitted to densities without disturbing anything that was fitted to VLE.
                    CHECK(close_to(trn->p(), ref->p(), 1e-10, 1e-6));

                    // Each saturated phase moves by the same -c, so for a PURE fluid every property
                    // change on vaporisation is invariant -- the volume change included.
                    const auto dvap = [](AbstractState& S, parameters key) {
                        return S.saturated_vapor_keyed_output(key) - S.saturated_liquid_keyed_output(key);
                    };
                    const auto dvap_v = [](AbstractState& S) {
                        return 1.0 / S.saturated_vapor_keyed_output(iDmolar) - 1.0 / S.saturated_liquid_keyed_output(iDmolar);
                    };
                    CHECK(close_to(dvap_v(*trn), dvap_v(*ref), 1e-9, 1e-18));
                    for (parameters key : {iHmolar, iSmolar, iUmolar, iHelmholtzmolar, iCvmolar, iCpmolar}) {
                        CAPTURE(get_parameter_information(key, "short"));
                        CHECK(close_to(dvap(*trn, key), dvap(*ref, key), 1e-9, 1e-8));
                    }

                    // Each saturated phase individually: the volume shifts by -c.
                    //
                    // The corresponding h -> h - p*c identity is deliberately NOT asserted here.
                    // SatL/SatV are set from (rho, T), so each carries its own converged pressure
                    // rather than a shared exact p, and the residual p mismatch feeds straight into
                    // h = u + p*v.  The identity is pinned rigorously in the single-phase test
                    // above, where p is an exact input.
                    for (int Q : {0, 1}) {
                        CAPTURE(Q);
                        const double vr = 1.0 / (Q == 0 ? ref->saturated_liquid_keyed_output(iDmolar) : ref->saturated_vapor_keyed_output(iDmolar));
                        const double vt = 1.0 / (Q == 0 ? trn->saturated_liquid_keyed_output(iDmolar) : trn->saturated_vapor_keyed_output(iDmolar));
                        CHECK(close_to(vt, vr - c, 1e-9, 1e-18));
                    }
                }
            }
        }
    }
}

TEST_CASE("Cubic volume translation: superancillary returns translated densities", "[cubic][volume_translation]") {
    for (const auto& backend : cubic_backends()) {
        for (const auto& fluid : probe_fluids()) {
            CAPTURE(backend);
            CAPTURE(fluid);
            const double c = 0.15 * covolume(backend, fluid);
            ASptr ref = make_translated(backend, fluid, 0.0);
            ASptr trn = make_translated(backend, fluid, c);
            auto& cb_ref = as_cubic(*ref);
            auto& cb_trn = as_cubic(*trn);

            const double T = 0.75 * trn->T_critical();

            const double rhoL_ref = cb_ref.calc_saturation_ancillary(iDmolar, 0, iT, T);
            const double rhoV_ref = cb_ref.calc_saturation_ancillary(iDmolar, 1, iT, T);
            const double rhoL_anc = cb_trn.calc_saturation_ancillary(iDmolar, 0, iT, T);
            const double rhoV_anc = cb_trn.calc_saturation_ancillary(iDmolar, 1, iT, T);
            const double p_anc = cb_trn.calc_saturation_ancillary(iP, 0, iT, T);
            CAPTURE(rhoL_anc);
            CAPTURE(rhoV_anc);

            // The sharp assertion, and the one that fails without the un-translation: the reduced
            // argument Ttilde = R*T*bm/am carries no c, so both calls hit the SAME point of the same
            // Chebyshev expansion.  The two volumes must therefore differ by exactly -c, with no
            // approximation error to hide behind.
            CHECK(close_to(1.0 / rhoL_anc, 1.0 / rhoL_ref - c, 1e-12, 1e-20));
            CHECK(close_to(1.0 / rhoV_anc, 1.0 / rhoV_ref - c, 1e-12, 1e-20));
            // ... and the pressure carries no translation at all.
            CHECK(close_to(p_anc, cb_ref.calc_saturation_ancillary(iP, 0, iT, T), 1e-12, 0.0));

            // Separately, the translated ancillary must still track the translated EOS flash to the
            // superancillary's own accuracy.  1e-3 is the figure the existing [cubic_superanc] test
            // documents; this check is about the un-translation being applied, not about the
            // Chebyshev fit being exact.
            trn->update(QT_INPUTS, 0, T);
            CHECK(close_to(rhoL_anc, trn->saturated_liquid_keyed_output(iDmolar), 1e-3, 0.0));
            CHECK(close_to(rhoV_anc, trn->saturated_vapor_keyed_output(iDmolar), 1e-3, 0.0));
            CHECK(close_to(p_anc, trn->p(), 1e-3, 0.0));
        }
    }
}

TEST_CASE("Cubic volume translation: DmolarT across the dome", "[cubic][volume_translation]") {
    // The sharp case.  With the untranslated superancillary bracketing the dome, a density that is
    // really inside the two-phase region was reported as single-phase liquid with a NEGATIVE
    // pressure -- a point on the unstable branch of the isotherm returned as if it were physical.
    for (const auto& backend : cubic_backends()) {
        for (const auto& fluid : probe_fluids()) {
            CAPTURE(backend);
            CAPTURE(fluid);
            const double c = 0.15 * covolume(backend, fluid);
            ASptr ref = make_translated(backend, fluid, 0.0);
            ASptr trn = make_translated(backend, fluid, c);
            auto& cb = as_cubic(*trn);

            const double T = 0.8 * trn->T_critical();

            // Part 0 -- the narrow band where the defect actually lived.
            //
            // A translated saturated liquid volume is exactly c smaller than the untranslated one,
            // so ANY volume strictly between them is inside the translated dome and outside the
            // untranslated one.  That band is only a fraction of a percent of the dome wide -- for
            // propane at c = 0.15*b it is frac < 0.004 -- so a uniform sweep steps straight over it.
            // Aim at it directly: the midpoint is the state that came back as single-phase liquid
            // at negative pressure.
            ref->update(QT_INPUTS, 0, T);
            const double vL_untranslated = 1.0 / ref->rhomolar();
            const double rho_band = 1.0 / (vL_untranslated - 0.5 * c);
            CAPTURE(vL_untranslated);
            CAPTURE(rho_band);
            REQUIRE_NOTHROW(trn->update(DmolarT_INPUTS, rho_band, T));
            CHECK(trn->phase() == iphase_twophase);
            CHECK(trn->p() > 0);
            CHECK(trn->Q() > 0);
            CHECK(trn->Q() < 1);

            // Part 1 -- the regression guard, and it must be driven from the TRUE dome.
            //
            // The defect was that update_DmolarT bracketed with the UNTRANSLATED superancillary.
            // For a positive c the translated liquid branch sits at higher density than the
            // untranslated one, so a band of genuinely two-phase states fell outside the bracket
            // and came back as single-phase liquid with a NEGATIVE pressure -- a point on the
            // unstable branch of the isotherm, returned as if it were physical.  Sweeping the dome
            // taken from the iterative flash is what exposes that; sweeping the ancillary's own
            // dome would be self-consistent either way and would guard nothing.
            trn->update(QT_INPUTS, 0, T);
            const double vL_true = 1.0 / trn->rhomolar();
            const double psat_true = trn->p();
            trn->update(QT_INPUTS, 1, T);
            const double vV_true = 1.0 / trn->rhomolar();
            REQUIRE(vV_true > vL_true);

            for (double frac : {0.05, 0.25, 0.5, 0.75, 0.95}) {
                const double rho = 1.0 / (vL_true + frac * (vV_true - vL_true));
                CAPTURE(frac);
                CAPTURE(rho);
                REQUIRE_NOTHROW(trn->update(DmolarT_INPUTS, rho, T));
                CHECK(trn->phase() == iphase_twophase);
                CHECK(trn->p() > 0);
                CHECK(trn->Q() > 0);
                CHECK(trn->Q() < 1);
                // The dome is bracketed with the ancillary, so the pressure it reports is the
                // ancillary's psat -- agreement with the flash is bounded by the ~1e-3 the
                // superancillary is documented to achieve, not by solver precision.
                CHECK(close_to(trn->p(), psat_true, 1e-3, 0.0));
            }

            // Part 2 -- the lever rule itself, on the densities update_DmolarT actually uses, so
            // the ancillary's approximation error does not leak into the expected quality.
            const double rhoL_anc = cb.calc_saturation_ancillary(iDmolar, 0, iT, T);
            const double rhoV_anc = cb.calc_saturation_ancillary(iDmolar, 1, iT, T);
            for (double frac : {0.25, 0.5, 0.75}) {
                const double rho = 1.0 / (1.0 / rhoL_anc + frac * (1.0 / rhoV_anc - 1.0 / rhoL_anc));
                CAPTURE(frac);
                REQUIRE_NOTHROW(trn->update(DmolarT_INPUTS, rho, T));
                const double Q_expected = (rhoV_anc * (rhoL_anc - rho)) / (rho * (rhoL_anc - rhoV_anc));
                CHECK(close_to(trn->Q(), Q_expected, 1e-9, 1e-12));
            }
        }
    }
}

TEST_CASE("Cubic volume translation: survives get_copy", "[cubic][volume_translation]") {
    // get_copy() builds a brand-new AbstractCubic, whose constructor zeroes the translation, and
    // copy_internals() replays only the kij matrix and the components vector.  Without copy_cm()
    // the TPD, critical and transient states silently run an untranslated EOS.
    for (const auto& backend : cubic_backends()) {
        CAPTURE(backend);
        const double c = 0.15 * covolume(backend, "n-Propane");
        ASptr trn = make_translated(backend, "n-Propane", c);
        auto& cb = as_cubic(*trn);

        std::shared_ptr<HelmholtzEOSMixtureBackend> clone(cb.get_copy(true));
        auto* clone_cb = dynamic_cast<AbstractCubicBackend*>(clone.get());
        REQUIRE(clone_cb != nullptr);
        CHECK(clone_cb->get_cubic()->get_cm(0) == Catch::Approx(c).epsilon(1e-15));

        // SatL/SatV of the clone must agree too, or a saturation call on the copy silently mixes a
        // translated parent with untranslated daughters.
        //
        // That has to be checked THROUGH a saturation call.  Reading the clone's own cm back routes
        // to the same storage the CHECK above already read, so it asserts nothing new -- deleting
        // copy_cm's linked_states fan-out leaves such a test green.  Drive a QT flash instead:
        // SatL/SatV are where saturation() puts its converged phases, so an untranslated daughter
        // shows up immediately.  The density agrees even when it is broken (DmolarT imposes it), so
        // compare the caloric properties, which are what actually diverge -- by ~11 % here.
        const double T = 0.7 * trn->T_critical();
        REQUIRE_NOTHROW(trn->update(QT_INPUTS, 0, T));
        REQUIRE_NOTHROW(clone->update(QT_INPUTS, 0, T));
        CHECK(close_to(clone->hmolar(), trn->hmolar(), 1e-9, 1e-6));
        CHECK(close_to(clone->smolar(), trn->smolar(), 1e-9, 1e-6));
        CHECK(close_to(clone->umolar(), trn->umolar(), 1e-9, 1e-6));
        CHECK(close_to(clone->rhomolar(), trn->rhomolar(), 1e-9, 1e-6));
        REQUIRE_NOTHROW(trn->update(QT_INPUTS, 1, T));
        REQUIRE_NOTHROW(clone->update(QT_INPUTS, 1, T));
        CHECK(close_to(clone->hmolar(), trn->hmolar(), 1e-9, 1e-6));
        CHECK(close_to(clone->rhomolar(), trn->rhomolar(), 1e-9, 1e-6));
    }
}

TEST_CASE("Cubic volume translation: the critical and reducing volumes move with c", "[cubic][volume_translation]") {
    // Three more paths reason about volume and all had to learn about c, but none of them is
    // reachable through an ordinary property call, so none was covered: reverting the translation in
    // all three at once left the suite green.  They are exercised directly here.
    for (const auto& backend : cubic_backends()) {
        CAPTURE(backend);
        const double b = covolume(backend, "n-Propane");
        const double c = 0.15 * b;
        ASptr ref = make_translated(backend, "n-Propane", 0.0);
        ASptr trn = make_translated(backend, "n-Propane", c);

        // calc_rhomolar_critical(): the Kazakov correlation returns an UNtranslated critical volume,
        // so the translated one is exactly c smaller.  An identity, not an approximation.
        CHECK(close_to(1.0 / trn->rhomolar_critical(), 1.0 / ref->rhomolar_critical() - c, 1e-12, 1e-20));

        // get_linear_reducing_parameters(): the same shift, through the mole-fraction weighted
        // average that the critical-point tracer starts from.  T_r must NOT move.
        double rhor_ref = 0, Tr_ref = 0, rhor_trn = 0, Tr_trn = 0;
        as_cubic(*ref).get_linear_reducing_parameters(rhor_ref, Tr_ref);
        as_cubic(*trn).get_linear_reducing_parameters(rhor_trn, Tr_trn);
        CHECK(close_to(1.0 / rhor_trn, 1.0 / rhor_ref - c, 1e-12, 1e-20));
        CHECK(close_to(Tr_trn, Tr_ref, 1e-14, 0.0));

        // get_critical_is_terminated(): the tracer stops when v < 1.1*(b_m - c_m), so the floor
        // moves with the translation.  Probe just inside and just outside the TRANSLATED floor --
        // the outside probe is the discriminating one, because a version still using the
        // untranslated b_m would see 1.15*(b - c) as below its own (larger) 1.1*b threshold and
        // terminate early.
        auto& cb = as_cubic(*trn);
        const double bmc = b - c;
        for (const auto& probe : {std::make_pair(1.05, true), std::make_pair(1.15, false)}) {
            const double v = probe.first * bmc;
            double delta = 1.0 / (v * cb.rhomolar_reducing()), tau = 1.0;
            CAPTURE(probe.first);
            CHECK(cb.get_critical_is_terminated(delta, tau) == probe.second);
        }
    }
}

TEST_CASE("Cubic volume translation: density bound tracks the pole", "[cubic][volume_translation]") {
    for (const auto& backend : cubic_backends()) {
        CAPTURE(backend);
        const double b = covolume(backend, "n-Propane");
        for (double frac : {0.0, 0.3, -0.3}) {
            const double c = frac * b;
            CAPTURE(c);
            ASptr AS = make_translated(backend, "n-Propane", c);
            auto& cb = as_cubic(*AS);
            const double pole = 1.0 / (b - c);
            // The bound must stay strictly inside the pole for BOTH signs.  With the untranslated
            // 0.9/b_m, a sufficiently negative c put the bound past the pole, where alpha^r takes
            // the log of a negative number but p, h and cp still return finite nonsense.
            CHECK(cb.calc_rhomolar_max_bound() < pole);
            CHECK(cb.calc_rhomolar_max_bound() == Catch::Approx(0.9 * pole).epsilon(1e-12));
        }
    }
}

TEST_CASE("Cubic volume translation: rejects a translation past the covolume", "[cubic][volume_translation]") {
    for (const auto& backend : cubic_backends()) {
        CAPTURE(backend);
        const double b = covolume(backend, "n-Propane");
        ASptr AS(AbstractState::factory(backend, "n-Propane"));
        // b - c <= 0 inverts the repulsive branch; it has to be rejected at set time because the
        // model does not fail loudly afterwards.
        CHECK_THROWS(AS->set_fluid_parameter_double(0, "cm", b));
        CHECK_THROWS(AS->set_fluid_parameter_double(0, "cm", 1.5 * b));
        CHECK_NOTHROW(AS->set_fluid_parameter_double(0, "cm", 0.9 * b));
    }
}

TEST_CASE("Cubic volume translation: moving Tc or pc cannot strand c past the covolume", "[cubic][volume_translation]") {
    // b = Omega_b*R*Tc/pc, so the covolume moves with the critical constants and "set c, then
    // shrink b" reaches the same forbidden state that setting c too large does -- but by a path
    // that used to be unguarded.  It fails silently rather than loudly: once b - c < 0, the
    // bracket 1 - (b_m - c_m)*rho in psi_minus is greater than 1 for EVERY rho, so the log never
    // goes NaN and p, h and c_p keep returning finite nonsense.  update_DmolarT does not consult
    // calc_rhomolar_max_bound either, so nothing downstream catches it.
    for (const auto& backend : cubic_backends()) {
        CAPTURE(backend);
        const double b = covolume(backend, "n-Propane");
        const double c = 0.9 * b;

        for (const std::string key : {std::string("Tcrit"), std::string("pcrit")}) {
            CAPTURE(key);
            ASptr AS = make_translated(backend, "n-Propane", c);
            const double before = AS->get_fluid_parameter_double(0, key);
            // b is proportional to Tc and inversely proportional to pc, so either of these halves
            // the covolume and leaves b - c = -0.4*b.
            const double bad = (key == "Tcrit") ? 0.5 * before : 2.0 * before;
            CHECK_THROWS(AS->set_fluid_parameter_double(0, key, bad));

            // Rejected is not enough -- it has to be rolled back.  set_Tci/set_pci write the
            // member first and refresh the alpha function and the b0 cache from it, so a throw
            // without a rollback would leave a state that is wrong in a different way.
            CHECK(AS->get_fluid_parameter_double(0, key) == Catch::Approx(before).epsilon(1e-15));
            CHECK(AS->get_fluid_parameter_double(0, "cm") == Catch::Approx(c).epsilon(1e-15));

            // And the rolled-back state must still be usable, which it is not if the derived
            // caches were left referring to the value that was rejected.
            REQUIRE_NOTHROW(AS->update(QT_INPUTS, 0, 0.7 * AS->T_critical()));
            CHECK(AS->p() > 0);
            CHECK(1.0 / AS->rhomolar() > b - c);
        }

        // A move in the safe direction still goes through -- the guard must not have frozen Tc/pc.
        ASptr ok = make_translated(backend, "n-Propane", 0.5 * b);
        const double Tc0 = ok->get_fluid_parameter_double(0, "Tcrit");
        REQUIRE_NOTHROW(ok->set_fluid_parameter_double(0, "Tcrit", 1.10 * Tc0));
        CHECK(ok->get_fluid_parameter_double(0, "Tcrit") == Catch::Approx(1.10 * Tc0).epsilon(1e-12));
    }
}

TEST_CASE("Cubic volume translation: the bare AbstractCubic setters enforce the pole condition", "[cubic][volume_translation]") {
    // The backend is not the only way in.  HEOS embeds an AbstractCubic through
    // ResidualHelmholtzGeneralizedCubic -- change_EOS(i, "SRK"/"Peng-Robinson") and the
    // "-SRK"/"-PengRobinson" fluid endings -- and those objects have no AbstractCubicBackend
    // wrapping them, so a check that lives only in set_fluid_parameter_double does not cover them.
    const std::vector<double> Tc{369.89}, pc{4.2512e6}, acentric{0.1521};
    const double R = 8.31446261815324;

    for (const auto& backend : cubic_backends()) {
        CAPTURE(backend);
        std::shared_ptr<AbstractCubic> bare;
        if (backend == "PR") {
            bare = std::make_shared<PengRobinson>(Tc, pc, acentric, R);
        } else {
            bare = std::make_shared<SRK>(Tc, pc, acentric, R);
        }
        const double b = bare->b0_ii(0);
        REQUIRE(b > 0);

        CHECK_THROWS(bare->set_cm(0, b));
        CHECK_THROWS(bare->set_cm(0, 1.5 * b));
        CHECK_THROWS(bare->set_cm(b));  // the broadcast form
        CHECK_THROWS(bare->set_cm_vector(std::vector<double>{b}));
        // A rejected value must leave the translation untouched rather than half-applied.
        CHECK(bare->get_cm(0) == Catch::Approx(0.0).margin(1e-300));

        CHECK_NOTHROW(bare->set_cm_vector(std::vector<double>{0.5 * b}));
        CHECK(bare->get_cm(0) == Catch::Approx(0.5 * b).epsilon(1e-15));
    }
}

TEST_CASE("Cubic volume translation: set/get round-trip", "[cubic][volume_translation]") {
    for (const auto& backend : cubic_backends()) {
        CAPTURE(backend);
        const double c = 0.2 * covolume(backend, "n-Propane");
        ASptr AS(AbstractState::factory(backend, "n-Propane"));
        CHECK(AS->get_fluid_parameter_double(0, "cm") == Catch::Approx(0.0).margin(1e-300));
        for (const char* alias : {"c", "cm", "c_m"}) {
            CAPTURE(alias);
            REQUIRE_NOTHROW(AS->set_fluid_parameter_double(0, alias, c));
            CHECK(AS->get_fluid_parameter_double(0, alias) == Catch::Approx(c).epsilon(1e-15));
        }
    }
}

TEST_CASE("Cubic volume translation: spinodal translates", "[cubic][volume_translation]") {
    for (const auto& backend : cubic_backends()) {
        CAPTURE(backend);
        const double c = 0.15 * covolume(backend, "n-Propane");
        ASptr ref = make_translated(backend, "n-Propane", 0.0);
        ASptr trn = make_translated(backend, "n-Propane", c);
        const double T = 0.8 * ref->T_critical();
        ref->update(QT_INPUTS, 0, T);
        trn->update(QT_INPUTS, 0, T);

        const std::vector<double> s0 = as_cubic(*ref).spinodal_densities();
        const std::vector<double> s1 = as_cubic(*trn).spinodal_densities();
        REQUIRE(s0.size() == s1.size());
        REQUIRE(!s0.empty());
        // -(dv/dP)_T is invariant under the translation, so the spinodal does not need re-deriving:
        // it simply moves with the volume.
        for (std::size_t i = 0; i < s0.size(); ++i) {
            CAPTURE(i);
            CHECK(close_to(1.0 / s1[i], 1.0 / s0[i] - c, 1e-9, 1e-18));
        }
    }
}

// ============================================================================================
// Mixtures
//
// The mixture column of the invariance table.  Two things differ from a pure fluid.  First the
// shifts are in c_m(z) = sum_i x_i c_i, EXCEPT the partial molar pair mu_i and ln(phi_i), which
// shift by the PURE-component c_i -- that is what makes the equifugacity condition cancel, and
// hence what makes VLE invariant for ARBITRARY c_i and not merely for equal ones.  Second, the
// property changes on vaporisation split: the two coexisting phases have different compositions, so
// c_m(y) != c_m(x) unless every c_i is equal, and Dvap V and Dvap H inherit that difference while
// Dvap S, U, A, c_p and c_v do not.
// ============================================================================================

namespace {

const char* kMixture = "Methane&Ethane&n-Propane";

/// Distinct per-component translations with mixed signs -- the case that a scalar c cannot express
/// and that only a correct dc_m/dx_i chain gets right.  Well inside methane's covolume (~2.7e-5).
const std::vector<double>& mixture_c() {
    static const std::vector<double> v{2.5e-6, -1.5e-6, 4.0e-6};
    return v;
}

ASptr make_mixture(const std::string& backend, const std::vector<double>& z, const std::vector<double>& c) {
    ASptr AS(AbstractState::factory(backend, kMixture));
    for (std::size_t i = 0; i < c.size(); ++i) {
        if (c[i] != 0.0) {
            AS->set_fluid_parameter_double(i, "cm", c[i]);
        }
    }
    AS->set_mole_fractions(std::vector<CoolPropDbl>(z.begin(), z.end()));
    return AS;
}

/// The invariant column of the property-change-on-vaporisation row, each with the absolute
/// tolerance floor appropriate to its own magnitude.
const std::vector<std::pair<parameters, double>> kMixtureDvapInvariants{
  {iSmolar, 1e-8}, {iUmolar, 1e-6}, {iHelmholtzmolar, 1e-6}, {iCpmolar, 1e-8}, {iCvmolar, 1e-8}};

/// Single-phase mixture states, derived from the phase envelope rather than hardcoded.
///
/// A fixed 5 MPa at 320 K sat 0.51 % below SRK's dew pressure (5.0256 MPa; PR 5.0931 MPa, 1.86 %).
/// That is closer than a routine tweak to a critical constant in the cubic JSON would move it, and
/// crossing in is not a soft failure: calc_fugacity_coefficient throws outright for 0 < Q < 1, so
/// the ln(phi_i) assertions would abort the whole case rather than fail informatively.
std::vector<StatePoint> mixture_single_phase_points(const std::string& backend, const std::vector<double>& z) {
    ASptr AS(AbstractState::factory(backend, kMixture));
    AS->set_mole_fractions(std::vector<CoolPropDbl>(z.begin(), z.end()));

    std::vector<StatePoint> pts;
    for (double T : {320.0, 450.0}) {
        double p_dew = -1, p_bub = -1;
        try {
            AS->update(QT_INPUTS, 1, T);
            p_dew = AS->p();
            AS->update(QT_INPUTS, 0, T);
            p_bub = AS->p();
        } catch (...) {
            // No envelope at this T: the mixture is supercritical there, so every pressure is
            // single-phase and the absolute fallbacks below cannot land inside a dome.
        }
        if (p_dew > 0 && p_bub > 0) {
            REQUIRE(p_bub > p_dew);
            pts.push_back({"vapour, 40 % below the dew line", T, 0.6 * p_dew});
            pts.push_back({"liquid, 3x the bubble pressure", T, 3.0 * p_bub});
        } else {
            pts.push_back({"supercritical", T, 5e6});
            pts.push_back({"supercritical, dense", T, 2e7});
        }
    }
    return pts;
}

double cm_of(const std::vector<double>& z, const std::vector<double>& c) {
    double s = 0;
    for (std::size_t i = 0; i < z.size(); ++i) {
        s += z[i] * c[i];
    }
    return s;
}

}  // namespace

TEST_CASE("Cubic volume translation: mixture single-phase invariances", "[cubic][volume_translation][mixture]") {
    const std::vector<double> z{0.30, 0.35, 0.35};
    const std::vector<double> c = mixture_c();
    const double cm = cm_of(z, c);
    const double R = 8.31446261815324;

    for (const auto& backend : cubic_backends()) {
        CAPTURE(backend);
        ASptr ref = make_mixture(backend, z, {0.0, 0.0, 0.0});
        ASptr trn = make_mixture(backend, z, c);

        // Single-phase states, so the PT flash returns a root rather than a phase split whose
        // composition would differ between the two runs.
        for (const auto& pt : mixture_single_phase_points(backend, z)) {
            {
                const double T = pt.T, p = pt.p;
                CAPTURE(pt.label);
                CAPTURE(T);
                CAPTURE(p);
                REQUIRE_NOTHROW(ref->update(PT_INPUTS, p, T));
                REQUIRE_NOTHROW(trn->update(PT_INPUTS, p, T));

                CHECK(close_to(1.0 / trn->rhomolar(), 1.0 / ref->rhomolar() - cm, 1e-9, 1e-18));
                CHECK(close_to(trn->hmolar(), ref->hmolar() - p * cm, 1e-9, 1e-8));
                CHECK(close_to(trn->gibbsmolar(), ref->gibbsmolar() - p * cm, 1e-9, 1e-8));
                CHECK(close_to(trn->smolar(), ref->smolar(), 1e-9, 1e-10));
                CHECK(close_to(trn->umolar(), ref->umolar(), 1e-9, 1e-8));
                CHECK(close_to(trn->helmholtzmolar(), ref->helmholtzmolar(), 1e-9, 1e-8));
                CHECK(close_to(trn->cpmolar(), ref->cpmolar(), 1e-9, 1e-10));
                CHECK(close_to(trn->cvmolar(), ref->cvmolar(), 1e-9, 1e-10));
                CHECK(
                  close_to(trn->isothermal_compressibility() / trn->rhomolar(), ref->isothermal_compressibility() / ref->rhomolar(), 1e-9, 1e-25));
                CHECK(close_to(trn->isobaric_expansion_coefficient() / trn->rhomolar(), ref->isobaric_expansion_coefficient() / ref->rhomolar(), 1e-9,
                               1e-20));

                // The one that needs the PURE-component c_i, not c_m -- and the one that the
                // d_A_term_dxi defect got wrong.
                for (std::size_t i = 0; i < z.size(); ++i) {
                    CAPTURE(i);
                    CHECK(close_to(std::log(trn->fugacity_coefficient(i)), std::log(ref->fugacity_coefficient(i)) - p * c[i] / (R * T), 1e-8, 1e-10));
                }

                // mu_i is NOT shifted by -p*c_m, even though g is.  G shifts by -p*sum_i n_i*c_i,
                // so mu_i = d(G)/dn_i picks up the PURE-component c_i -- the same c_i that ln(phi_i)
                // carries, which is the consistency the two have to satisfy.  The distinction is
                // invisible for a pure fluid, where mu = g and c_i = c_m.
                for (std::size_t i = 0; i < z.size(); ++i) {
                    CAPTURE(i);
                    const double dmu = trn->chemical_potential(i) - ref->chemical_potential(i);
                    CAPTURE(dmu);
                    // Sharp enough to separate the two hypotheses: -p*c_i and -p*c_m differ by
                    // 17 J/mol here at p = 2e7 Pa, and the floor is four orders below that.
                    CAPTURE(-p * cm);
                    CHECK(close_to(dmu, -p * c[i], 1e-6, 1e-4));
                }

                // The second virial coefficient carries the mixture translation too.  B is
                // evaluated in the zero-density limit at the current tau, so it depends on T and
                // the composition but not on the state's own density.
                CHECK(close_to(trn->Bvirial(), ref->Bvirial() - cm, 1e-9, 1e-18));

                CHECK(std::abs(trn->speed_sound() / ref->speed_sound() - 1.0) > 1e-6);
                check_v_scaled_quantities(*ref, *trn, cm);
            }
        }
    }
}

TEST_CASE("Cubic volume translation: mixture VLE is invariant for arbitrary c_i", "[cubic][volume_translation][mixture]") {
    // The headline mixture result.  ln(phi_i) shifts by -P*c_i/(RT), which is the SAME in both
    // phases at the same (T, P), so it cancels in the equifugacity condition regardless of whether
    // the c_i are equal.  Bubble and dew pressures and all K-values must therefore be untouched.
    //
    // This is the case that fails on the unfixed A_term derivatives: with c = 1e-5 on a
    // methane/n-decane pair the bubble pressure moved from 10.20 MPa to 6.87 MPa.
    const std::vector<double> z{0.30, 0.35, 0.35};
    const std::vector<double> c = mixture_c();

    for (const auto& backend : cubic_backends()) {
        CAPTURE(backend);
        for (double T : {220.0, 250.0}) {
            CAPTURE(T);
            for (int Q : {0, 1}) {
                CAPTURE(Q);
                ASptr ref = make_mixture(backend, z, {0.0, 0.0, 0.0});
                ASptr trn = make_mixture(backend, z, c);
                REQUIRE_NOTHROW(ref->update(QT_INPUTS, Q, T));
                REQUIRE_NOTHROW(trn->update(QT_INPUTS, Q, T));

                CHECK(close_to(trn->p(), ref->p(), 1e-6, 1.0));

                const std::vector<CoolPropDbl> xr = ref->mole_fractions_liquid(), yr = ref->mole_fractions_vapor();
                const std::vector<CoolPropDbl> xt = trn->mole_fractions_liquid(), yt = trn->mole_fractions_vapor();
                REQUIRE(xr.size() == xt.size());
                for (std::size_t i = 0; i < xr.size(); ++i) {
                    CAPTURE(i);
                    CHECK(close_to(xt[i], xr[i], 1e-6, 1e-9));
                    CHECK(close_to(yt[i], yr[i], 1e-6, 1e-9));
                    // K-values are what a flash actually consumes.
                    CHECK(close_to(yt[i] / xt[i], yr[i] / xr[i], 1e-6, 1e-9));
                }

                // Property changes on vaporisation.  S, U and A survive because the c-dependent
                // parts cancel between h, P*v and T*s, and c_p and c_v survive because each phase's
                // own value does; V and H carry the composition difference in c_m, which is nonzero
                // precisely because the c_i differ.
                const double dcm = cm_of(std::vector<double>(yr.begin(), yr.end()), c) - cm_of(std::vector<double>(xr.begin(), xr.end()), c);
                CAPTURE(dcm);
                const auto dvap = [](AbstractState& S, parameters key) {
                    return S.saturated_vapor_keyed_output(key) - S.saturated_liquid_keyed_output(key);
                };
                // The absolute floors are per property because these are differences of two
                // large numbers: S and the heat capacities are O(10..100) J/mol/K, while U and A are
                // O(1e4) J/mol and so carry ~1e-11 relative round-off of their own.
                for (const auto& probe : kMixtureDvapInvariants) {
                    CAPTURE(get_parameter_information(probe.first, "short"));
                    CHECK(close_to(dvap(*trn, probe.first), dvap(*ref, probe.first), 1e-6, probe.second));
                }
                const double dV_ref = 1.0 / ref->saturated_vapor_keyed_output(iDmolar) - 1.0 / ref->saturated_liquid_keyed_output(iDmolar);
                const double dV_trn = 1.0 / trn->saturated_vapor_keyed_output(iDmolar) - 1.0 / trn->saturated_liquid_keyed_output(iDmolar);
                CHECK(close_to(dV_trn, dV_ref - dcm, 1e-6, 1e-15));
                // Tolerance tied to the SHIFT, not to Dvap H.  The shift is p*dcm ~ 0.27 J/mol at the
                // 220 K dew point while Dvap H is ~1.4e4 J/mol, so a relative test on the latter
                // accepted a 5 % error in the former.  The measured residual across all four
                // (backend, T, Q) points is <= 1.4e-10 J/mol, so this still clears the numerical
                // floor by four orders while pinning the shift itself to 1e-5 relative.
                const double dH_shift = ref->p() * dcm;
                CAPTURE(dH_shift);
                CHECK(close_to(dvap(*trn, iHmolar), dvap(*ref, iHmolar) - dH_shift, 0.0, 1e-5 * std::abs(dH_shift) + 1e-8));
            }
        }
    }
}

TEST_CASE("Cubic volume translation: c_all broadcasts to every component", "[cubic][volume_translation][mixture]") {
    // "c_all" is the spelling that keeps the pre-vectorisation behaviour, where the index was
    // accepted and then ignored.  Nothing else exercises it, so a regression there would silently
    // change what callers relying on those old fluid-wide semantics get back.
    const std::vector<double> z{0.30, 0.35, 0.35};
    const std::size_t N = mixture_c().size();

    for (const auto& backend : cubic_backends()) {
        CAPTURE(backend);
        // Seeded with distinct per-component values, so a broadcast that reaches only one of them
        // is visible rather than being masked by the components already agreeing.
        ASptr AS = make_mixture(backend, z, mixture_c());
        for (std::size_t i = 0; i < N; ++i) {
            CAPTURE(i);
            CHECK(AS->get_fluid_parameter_double(i, "cm") == Catch::Approx(mixture_c()[i]).epsilon(1e-15));
        }

        // The index is ignored for c_all, so naming a different one must give the same result.
        for (std::size_t named = 0; named < N; ++named) {
            const double broadcast = 1.25e-6 * static_cast<double>(named + 1);
            CAPTURE(named);
            CAPTURE(broadcast);
            REQUIRE_NOTHROW(AS->set_fluid_parameter_double(named, "c_all", broadcast));
            for (std::size_t i = 0; i < N; ++i) {
                CAPTURE(i);
                CHECK(AS->get_fluid_parameter_double(i, "cm") == Catch::Approx(broadcast).epsilon(1e-15));
            }
        }
    }
}

// ============================================================================================
// Composition derivatives
//
// These operate on a bare AbstractCubic rather than through a backend, for two reasons: the
// composition-derivative chain is pure algebra that needs no state object around it, and doing it
// this way keeps the check independent of get_copy()/linked_states plumbing.
//
// The reference is a central difference of the *value* function.  For second derivatives the
// reference is a central difference of the *analytic first* derivative -- nesting two finite
// differences would lose far more precision than the defect being hunted.
// ============================================================================================

namespace {

/// A three-component methane/ethane/propane mixture: enough to make every composition index
/// distinct and to give b_m a real spread.
struct CompDerivCase
{
    std::shared_ptr<AbstractCubic> cubic;
    std::vector<double> x{0.3, 0.35, 0.35};
    double delta = 5000.0;  ///< rho_r defaults to 1, so delta IS the molar density [mol/m^3]
    double tau = 1.0 / 250.0;
};

CompDerivCase make_comp_case(const std::string& which, const std::vector<double>& cvec) {
    const std::vector<double> Tc{190.564, 305.322, 369.89};
    const std::vector<double> pc{4.5992e6, 4.8722e6, 4.2512e6};
    const std::vector<double> acentric{0.01142, 0.0995, 0.1521};
    const double R = 8.31446261815324;
    CompDerivCase c;
    if (which == "PR") {
        c.cubic = std::make_shared<PengRobinson>(Tc, pc, acentric, R);
    } else {
        c.cubic = std::make_shared<SRK>(Tc, pc, acentric, R);
    }
    c.cubic->set_cm_vector(cvec);
    return c;
}

/// The translations to sweep.  The uniform vector reproduces what a single scalar c used to mean;
/// the mixed one is the case that actually exercises dc_m/dx_i, and it deliberately mixes signs
/// because both occur in the real parameter sets.  The smallest covolume here is methane's, about
/// 2.7e-5 m^3/mol, so all of these stay well inside the pole.
const std::vector<std::vector<double>>& comp_deriv_translations() {
    static const std::vector<std::vector<double>> v{
      {0.0, 0.0, 0.0},      // control: must pass before and after the fix
      {3e-6, 3e-6, 3e-6},   // uniform
      {2e-6, -1e-6, 5e-6},  // genuinely per-component, mixed signs
    };
    return v;
}

/// Central difference of f with respect to x_i, honouring the same x_N convention the analytic
/// derivatives use: when x_N is dependent, perturbing x_i has to be compensated in x_N.
double fd_dxi(const std::function<double(const std::vector<double>&)>& f, const std::vector<double>& x, std::size_t i, bool xN_independent,
              double dz) {
    std::vector<double> xp = x, xm = x;
    xp[i] += dz;
    xm[i] -= dz;
    if (!xN_independent) {
        xp[xp.size() - 1] -= dz;
        xm[xm.size() - 1] += dz;
    }
    return (f(xp) - f(xm)) / (2 * dz);
}

/// Relative error, falling back to absolute when the analytic value is ~0.
double deriv_err(double numeric, double analytic) {
    return (std::abs(analytic) > 1e-12) ? std::abs(numeric / analytic - 1) : std::abs(numeric - analytic);
}

}  // namespace

TEST_CASE("Cubic composition derivatives carry the volume translation", "[cubic][volume_translation][cubic_compderiv]") {
    // 1e-6 sits comfortably above the floor of a central difference on well-scaled quantities and
    // five orders below the defect it guards: dropping the (1 + c*rho) factor from d_A_term_dxi is
    // a 2.4 % error at c = 5e-6 and 9.1 % at c = 2e-5.
    const double dz = 1e-7, tol = 1e-6;

    for (const auto& which : cubic_backends()) {
        CAPTURE(which);
        // c = 0 must pass both before and after the fix -- it is the control that proves the
        // harness itself is sound.  The nonzero vectors are the reproduction.
        for (const auto& cvec : comp_deriv_translations()) {
            CAPTURE(cvec[0]);
            CAPTURE(cvec[1]);
            CAPTURE(cvec[2]);
            CompDerivCase cs = make_comp_case(which, cvec);
            AbstractCubic& C = *cs.cubic;
            const std::vector<double>& x = cs.x;
            const double d = cs.delta, tau = cs.tau;

            for (bool xNi : {true, false}) {
                CAPTURE(xNi);
                for (std::size_t i = 0; i < x.size(); ++i) {
                    CAPTURE(i);

                    // b_m carries no translation, so this is a second control.
                    CHECK(deriv_err(fd_dxi([&](const std::vector<double>& z) { return C.bm_term(z); }, x, i, xNi, dz), C.d_bm_term_dxi(x, i, xNi))
                          < tol);

                    // psi_minus depends on x only through (b_m - c_m).
                    for (std::size_t id = 0; id <= 4; ++id) {
                        CAPTURE(id);
                        CHECK(deriv_err(fd_dxi([&](const std::vector<double>& z) { return C.psi_minus(d, z, 0, id); }, x, i, xNi, dz),
                                        C.d_psi_minus_dxi(d, x, 0, id, i, xNi))
                              < tol);
                    }

                    // PI_12 carries c_m in both factors.
                    for (std::size_t id = 0; id <= 2; ++id) {
                        CAPTURE(id);
                        CHECK(deriv_err(fd_dxi([&](const std::vector<double>& z) { return C.PI_12(d, z, id); }, x, i, xNi, dz),
                                        C.d_PI_12_dxi(d, x, id, i, xNi))
                              < tol);
                    }

                    // A_term -- the one that drops (1 + c_m*rho).
                    CHECK(deriv_err(fd_dxi([&](const std::vector<double>& z) { return C.A_term(d, z); }, x, i, xNi, dz), C.d_A_term_dxi(d, x, i, xNi))
                          < tol);

                    for (std::size_t id = 0; id <= 4; ++id) {
                        CAPTURE(id);
                        CHECK(deriv_err(fd_dxi([&](const std::vector<double>& z) { return C.psi_plus(d, z, id); }, x, i, xNi, dz),
                                        C.d_psi_plus_dxi(d, x, id, i, xNi))
                              < tol);
                    }

                    // alphar closes the loop: it is what the backend actually consumes.
                    for (std::size_t it = 0; it <= 1; ++it) {
                        for (std::size_t id = 0; id <= 2; ++id) {
                            CAPTURE(it);
                            CAPTURE(id);
                            CHECK(deriv_err(fd_dxi([&](const std::vector<double>& z) { return C.alphar(tau, d, z, it, id); }, x, i, xNi, dz),
                                            C.d_alphar_dxi(tau, d, x, it, id, i, xNi))
                                  < tol);
                        }
                    }

                    // Second derivative of A, referenced against the VALUE function via a
                    // second-order central difference.
                    //
                    // This is deliberately independent of d_A_term_dxi.  Differencing the analytic
                    // first derivative -- which the loop below does, at much better precision --
                    // cannot see an error that is a common multiplicative factor on both orders,
                    // and that is exactly the shape of the missing (1 + c*rho): with only that
                    // check, d2_A_term_dxidxj passed against the unfixed code.
                    {
                        const double h = 1e-4;
                        std::vector<double> xp = x, xm = x;
                        xp[i] += h;
                        xm[i] -= h;
                        if (!xNi) {
                            xp[xp.size() - 1] -= h;
                            xm[xm.size() - 1] += h;
                        }
                        const double d2_num = (C.A_term(d, xp) - 2 * C.A_term(d, x) + C.A_term(d, xm)) / (h * h);
                        // A second difference loses about half the available digits, so 1e-4
                        // relative -- still two orders below the 1.5 % defect at this c and rho.
                        CHECK(deriv_err(d2_num, C.d2_A_term_dxidxj(d, x, i, i, xNi)) < 1e-4);
                    }

                    // Second derivatives, referenced against the analytic first derivative.
                    for (std::size_t j = 0; j < x.size(); ++j) {
                        CAPTURE(j);
                        CHECK(deriv_err(fd_dxi([&](const std::vector<double>& z) { return C.d_A_term_dxi(d, z, i, xNi); }, x, j, xNi, dz),
                                        C.d2_A_term_dxidxj(d, x, i, j, xNi))
                              < tol);
                        CHECK(deriv_err(fd_dxi([&](const std::vector<double>& z) { return C.d_PI_12_dxi(d, z, 0, i, xNi); }, x, j, xNi, dz),
                                        C.d2_PI_12_dxidxj(d, x, 0, i, j, xNi))
                              < tol);
                        CHECK(deriv_err(fd_dxi([&](const std::vector<double>& z) { return C.d_psi_minus_dxi(d, z, 0, 0, i, xNi); }, x, j, xNi, dz),
                                        C.d2_psi_minus_dxidxj(d, x, 0, 0, i, j, xNi))
                              < tol);
                    }
                }
            }
        }
    }
}

#endif
