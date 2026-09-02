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
                }
            }
        }
    }
}

TEST_CASE("Cubic volume translation: second virial coefficient shifts by -c", "[cubic][volume_translation]") {
    for (const auto& backend : cubic_backends()) {
        for (const auto& fluid : probe_fluids()) {
            CAPTURE(backend);
            CAPTURE(fluid);
            const double c = 0.15 * covolume(backend, fluid);
            ASptr ref = make_translated(backend, fluid, 0.0);
            ASptr trn = make_translated(backend, fluid, c);
            const double T = 1.25 * ref->T_critical();
            ref->update(DmolarT_INPUTS, 1e-3, T);
            trn->update(DmolarT_INPUTS, 1e-3, T);
            CHECK(close_to(trn->Bvirial(), ref->Bvirial() - c, 1e-9, 1e-18));
        }
    }
}

TEST_CASE("Cubic volume translation: saturation invariances", "[cubic][volume_translation]") {
    for (const auto& backend : cubic_backends()) {
        for (const auto& fluid : probe_fluids()) {
            CAPTURE(backend);
            CAPTURE(fluid);
            const double c = 0.15 * covolume(backend, fluid);
            ASptr ref = make_translated(backend, fluid, 0.0);
            ASptr trn = make_translated(backend, fluid, c);

            const double T = 0.75 * ref->T_critical();
            REQUIRE_NOTHROW(ref->update(QT_INPUTS, 0, T));
            REQUIRE_NOTHROW(trn->update(QT_INPUTS, 0, T));

            // The vapour pressure is the headline invariance: it is what lets a translation be
            // fitted to densities without disturbing anything that was fitted to VLE.
            CHECK(close_to(trn->p(), ref->p(), 1e-10, 1e-6));

            // Each saturated phase moves by the same -c, so for a PURE fluid every property change
            // on vaporisation is invariant -- the volume change included.
            const auto dvap = [](AbstractState& S, parameters key) {
                return S.saturated_vapor_keyed_output(key) - S.saturated_liquid_keyed_output(key);
            };
            const auto dvap_v = [](AbstractState& S) {
                return 1.0 / S.saturated_vapor_keyed_output(iDmolar) - 1.0 / S.saturated_liquid_keyed_output(iDmolar);
            };
            CHECK(close_to(dvap_v(*trn), dvap_v(*ref), 1e-9, 1e-18));
            for (parameters key : {iHmolar, iSmolar, iUmolar, iCvmolar, iCpmolar}) {
                CAPTURE(get_parameter_information(key, "short"));
                CHECK(close_to(dvap(*trn, key), dvap(*ref, key), 1e-9, 1e-8));
            }

            // Each saturated phase individually: the volume shifts by -c.
            //
            // The corresponding h -> h - p*c identity is deliberately NOT asserted here.  SatL/SatV
            // are set from (rho, T), so each carries its own converged pressure rather than a shared
            // exact p, and the residual p mismatch feeds straight into h = u + p*v.  The identity is
            // pinned rigorously in the single-phase test above, where p is an exact input.
            for (int Q : {0, 1}) {
                CAPTURE(Q);
                const double vr = 1.0 / (Q == 0 ? ref->saturated_liquid_keyed_output(iDmolar) : ref->saturated_vapor_keyed_output(iDmolar));
                const double vt = 1.0 / (Q == 0 ? trn->saturated_liquid_keyed_output(iDmolar) : trn->saturated_vapor_keyed_output(iDmolar));
                CHECK(close_to(vt, vr - c, 1e-9, 1e-18));
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
        CHECK(clone_cb->get_cubic()->get_cm() == Catch::Approx(c).epsilon(1e-15));

        // SatL/SatV of the clone must agree too, or a saturation call on the copy silently mixes a
        // translated parent with untranslated daughters.
        CHECK(clone_cb->get_fluid_parameter_double(0, "cm") == Catch::Approx(c).epsilon(1e-15));
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

TEST_CASE("Cubic volume translation: set/get round-trip", "[cubic][volume_translation]") {
    for (const auto& backend : cubic_backends()) {
        CAPTURE(backend);
        const double c = 0.2 * covolume(backend, "n-Propane");
        ASptr AS(AbstractState::factory(backend, "n-Propane"));
        CHECK(AS->get_fluid_parameter_double(0, "cm") == Catch::Approx(0.0));
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

#endif
