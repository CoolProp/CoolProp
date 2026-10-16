// Pressure-input flashes (P,H / P,S / P,U / P,D) for pure cubic fluids.
//
// GitHub #3441: the subcritical branch of p_phase_determination_pure_or_pseudopure asked the
// EOS-level (Helmholtz) superancillary for Tsat(p) and, finding none on a cubic backend, fell
// through to the classic pL/pV ancillaries -- which cubics never populate -- and threw
// "type not set".  The cubic Chebyshev superancillary now answers that query instead.

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/Configuration.h"
#    include "CoolProp/DataStructures.h"
#    include "../Backends/Cubics/CubicBackend.h"

#    include <memory>
#    include <string>
#    include <vector>

using namespace CoolProp;

namespace {

// Restores ENABLE_SUPERANCILLARIES on scope exit so a failing REQUIRE cannot leak the setting.
struct SuperancFlagGuard
{
    bool saved;
    explicit SuperancFlagGuard(bool value) : saved(get_config_bool(ENABLE_SUPERANCILLARIES)) {
        set_config_bool(ENABLE_SUPERANCILLARIES, value);
    }
    ~SuperancFlagGuard() {
        set_config_bool(ENABLE_SUPERANCILLARIES, saved);
    }
    SuperancFlagGuard(const SuperancFlagGuard&) = delete;
    SuperancFlagGuard& operator=(const SuperancFlagGuard&) = delete;
    SuperancFlagGuard(SuperancFlagGuard&&) = delete;
    SuperancFlagGuard& operator=(SuperancFlagGuard&&) = delete;
};

// Flash the state held by `ref` back from (p, other) and check T, rho and the phase come back.
void check_p_roundtrip(AbstractState& AS, AbstractState& ref, input_pairs pair, double other) {
    CAPTURE(get_input_pair_short_desc(pair), ref.T(), ref.p(), ref.Q(), other);
    // H,P and D,P take the pressure second; P,S and P,U take it first.
    if (pair == HmolarP_INPUTS || pair == DmolarP_INPUTS) {
        AS.update(pair, other, ref.p());
    } else {
        AS.update(pair, ref.p(), other);
    }
    CHECK(AS.T() == Catch::Approx(ref.T()).epsilon(1e-8));
    CHECK(AS.rhomolar() == Catch::Approx(ref.rhomolar()).epsilon(1e-8));
    const bool ref_twophase = ref.phase() == iphase_twophase;
    CHECK((AS.phase() == iphase_twophase) == ref_twophase);
    if (ref_twophase) {
        CHECK(AS.Q() == Catch::Approx(ref.Q()).margin(1e-8));
    }
}

void check_all_p_pairs(AbstractState& AS, AbstractState& ref) {
    check_p_roundtrip(AS, ref, HmolarP_INPUTS, ref.hmolar());
    check_p_roundtrip(AS, ref, PSmolar_INPUTS, ref.smolar());
    check_p_roundtrip(AS, ref, PUmolar_INPUTS, ref.umolar());
    check_p_roundtrip(AS, ref, DmolarP_INPUTS, ref.rhomolar());
}

}  // namespace

TEST_CASE("Cubic P+{H,S,U,D} flashes for pure fluids (#3441)", "[cubic][flash][3441]") {
    const auto backend = GENERATE(as<std::string>{}, "SRK", "PR");
    const auto fluid = GENERATE(as<std::string>{}, "CO2", "Propane", "Water");
    const bool superanc_flag = GENERATE(true, false);
    CAPTURE(backend, fluid, superanc_flag);
    SuperancFlagGuard guard(superanc_flag);

    std::unique_ptr<AbstractState> AS(AbstractState::factory(backend, fluid));
    std::unique_ptr<AbstractState> ref(AbstractState::factory(backend, fluid));
    const double Tc = AS->T_critical();
    const double pc = AS->p_critical();

    SECTION("issue reproducer: superheated SRK CO2 at 50 bar") {
        if (backend == "SRK" && fluid == "CO2") {
            ref->update(PT_INPUTS, 50e5, 320.0);
            check_all_p_pairs(*AS, *ref);
        }
    }
    SECTION("subcritical gas and liquid") {
        for (double Tr : {0.6, 0.8, 0.95}) {
            ref->update(QT_INPUTS, 0, Tr * Tc);
            const double psat = ref->p();
            // Vapour: 5 K above Tsat at psat; liquid: 5 K below.
            ref->update(PT_INPUTS, psat, Tr * Tc + 5);
            check_all_p_pairs(*AS, *ref);
            ref->update(PT_INPUTS, psat, Tr * Tc - 5);
            check_all_p_pairs(*AS, *ref);
        }
    }
    SECTION("inside the dome") {
        // Build the reference from the superancillary densities via D,T rather than from a Q,T
        // update: the cubic's iterative Q,T saturation is only converged to ~1e-7 in p near Tc,
        // whereas the superancillary satisfies equal Gibbs energy to ~1e-15.
        auto* cubic = dynamic_cast<AbstractCubicBackend*>(ref.get());
        REQUIRE(cubic != nullptr);
        for (double Tr : {0.6, 0.8, 0.95}) {
            const double T = Tr * Tc;
            const double rhoL = cubic->calc_saturation_ancillary(iDmolar, 0, iT, T);
            const double rhoV = cubic->calc_saturation_ancillary(iDmolar, 1, iT, T);
            for (double Q : {0.0, 0.3, 1.0}) {
                ref->update(DmolarT_INPUTS, 1 / (Q / rhoV + (1 - Q) / rhoL), T);
                if (Q > 0 && Q < 1) {
                    check_all_p_pairs(*AS, *ref);
                } else {
                    // Saturated endpoints: within rounding of the boundary, so only T and p are
                    // meaningful (phase may legitimately resolve to either side).
                    AS->update(HmolarP_INPUTS, ref->hmolar(), ref->p());
                    CHECK(AS->T() == Catch::Approx(ref->T()).epsilon(1e-7));
                }
            }
        }
    }
    SECTION("near-critical subcritical pressure") {
        // Between the top of the superancillary and the critical pressure.
        ref->update(PT_INPUTS, 0.9999 * pc, 1.01 * Tc);
        check_all_p_pairs(*AS, *ref);
        ref->update(PT_INPUTS, 0.9999 * pc, 0.99 * Tc);
        check_all_p_pairs(*AS, *ref);
    }
    SECTION("hook rejects pressures outside the superancillary") {
        auto* cubic = dynamic_cast<AbstractCubicBackend*>(AS.get());
        REQUIRE(cubic != nullptr);
        CHECK_THROWS_AS(cubic->superanc_saturation_at_p(0), ValueError);
        CHECK_THROWS_AS(cubic->superanc_saturation_at_p(-1), ValueError);
        CHECK_THROWS_AS(cubic->superanc_saturation_at_p(1e-300), ValueError);
        CHECK_THROWS_AS(cubic->superanc_saturation_at_p(1.01 * pc), ValueError);
        const auto sat = cubic->superanc_saturation_at_p(0.5 * pc);
        REQUIRE(sat.has_value());
        if (sat) {  // REQUIRE already aborts; the branch is for clang-tidy's optional-access check
            CHECK(sat->T < Tc);
            CHECK(sat->rhomolarL > sat->rhomolarV);
        }
    }
    SECTION("supercritical pressure") {
        ref->update(PT_INPUTS, 1.5 * pc, 1.2 * Tc);
        check_all_p_pairs(*AS, *ref);
    }
}

#endif
