// Pressure-input (P,H / P,S / P,U / P,D) and temperature-input (H,T / S,T / T,U / D,T) flashes for
// pure cubic fluids.
//
// GitHub #3441: the subcritical branch of p_phase_determination_pure_or_pseudopure asked the
// EOS-level (Helmholtz) superancillary for Tsat(p) and, finding none on a cubic backend, fell
// through to the classic pL/pV ancillaries -- which cubics never populate -- and threw
// "type not set".  T_phase_determination_pure_or_pseudopure did the same on the T side.  The cubic
// Chebyshev superancillary now answers both, through superanc_saturation_at_p / _at_T.

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

TEST_CASE("Cubic P-input flashes for pure fluids (#3441)", "[cubic][flash][3441]") {
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
        // Build the reference directly from the superancillary densities via D,T, so it does not
        // depend on which solver the cubic's Q,T update uses.
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

TEST_CASE("Cubic superancillary hook: cached table bounds follow parameter changes", "[cubic][3441]") {
    // superanc_saturation_at_p caches the temperatures at the ends of the table.  Changing the alpha
    // function on an instance that has already filled that cache must give the same answer as a fresh
    // instance with the new alpha.
    for (const auto& backend : std::vector<std::string>{"SRK", "PR"}) {
        CAPTURE(backend);
        std::unique_ptr<AbstractState> AS(AbstractState::factory(backend, "Propane"));
        auto* cubic = dynamic_cast<AbstractCubicBackend*>(AS.get());
        REQUIRE(cubic != nullptr);
        const double p = 0.3 * AS->p_critical();
        const auto before = cubic->superanc_saturation_at_p(p);  // fills the cache
        REQUIRE(before.has_value());

        AS->set_cubic_alpha_C(0, "MC", 0.9, -0.3, 0.2);
        std::unique_ptr<AbstractState> fresh(AbstractState::factory(backend, "Propane"));
        fresh->set_cubic_alpha_C(0, "MC", 0.9, -0.3, 0.2);
        auto* fresh_cubic = dynamic_cast<AbstractCubicBackend*>(fresh.get());
        REQUIRE(fresh_cubic != nullptr);

        const auto after = cubic->superanc_saturation_at_p(p);
        const auto expected = fresh_cubic->superanc_saturation_at_p(p);
        REQUIRE(after.has_value());
        REQUIRE(expected.has_value());
        if (after && expected && before) {
            CHECK(after->T == Catch::Approx(expected->T).epsilon(1e-13));
            CHECK(after->rhomolarL == Catch::Approx(expected->rhomolarL).epsilon(1e-12));
            CHECK(std::abs(after->T / before->T - 1) > 1e-4);  // the alpha change really moved Tsat
        }
    }
}

// T-input flashes (H,T / S,T / T,U / D,T) for pure cubic fluids.  Subcritical H,T, S,T and T,U threw
// "type not set" for the same reason as #3441, on the T side: T_phase_determination_pure_or_pseudopure
// asked get_superanc() and fell through to the classic ancillaries cubics never populate.
namespace {
void check_T_roundtrip(AbstractState& AS, AbstractState& ref) {
    CAPTURE(ref.T(), ref.p(), ref.Q());
    const double T = ref.T();
    const bool ref_twophase = ref.phase() == iphase_twophase;
    auto check = [&](const char* name, parameters key) {
        CAPTURE(name);
        // Always: the returned state reproduces both inputs.
        CHECK(AS.T() == T);
        CHECK(AS.keyed_output(key) == Catch::Approx(ref.keyed_output(key)).epsilon(1e-9));
        // In a cold liquid (dh/dp)_T = v(1 - T*beta) > 0, so compressing raises h above h_L(T) and a
        // compressed-liquid (T, h) also lies inside the dome: H,T is not unique there, and the flash
        // (like HEOS) returns the two-phase state.  Everywhere else the state must match the reference.
        if (key == iHmolar && !ref_twophase && AS.phase() == iphase_twophase) {
            CHECK(ref.phase() == iphase_liquid);
            // The exemption must be the real non-uniqueness, checked independently of the flash:
            // the reference h lies inside [h_L(T), h_V(T)] from a separate Q,T update.
            std::unique_ptr<AbstractState> sat(AbstractState::factory(AS.backend_name(), AS.fluid_names()));
            sat->update(QT_INPUTS, 0, T);
            const double hL = sat->hmolar();
            sat->update(QT_INPUTS, 1, T);
            const double hV = sat->hmolar();
            CHECK(ref.hmolar() >= hL);
            CHECK(ref.hmolar() <= hV);
            CHECK(AS.Q() >= 0);
            CHECK(AS.Q() <= 1);
            return;
        }
        CHECK(AS.rhomolar() == Catch::Approx(ref.rhomolar()).epsilon(1e-8));
        // 1e-7: a stiff compressed liquid amplifies density round-off ~1e3-1e4-fold in p.
        CHECK(AS.p() == Catch::Approx(ref.p()).epsilon(1e-7));
        CHECK((AS.phase() == iphase_twophase) == ref_twophase);
        if (ref_twophase) {
            CHECK(AS.Q() == Catch::Approx(ref.Q()).margin(1e-8));
        }
    };
    AS.update(HmolarT_INPUTS, ref.hmolar(), T);
    check("HmolarT", iHmolar);
    AS.update(SmolarT_INPUTS, ref.smolar(), T);
    check("SmolarT", iSmolar);
    AS.update(TUmolar_INPUTS, T, ref.umolar());
    check("TUmolar", iUmolar);
    AS.update(DmolarT_INPUTS, ref.rhomolar(), T);
    check("DmolarT", iDmolar);
}
}  // namespace

TEST_CASE("Cubic T-input flashes for pure fluids", "[cubic][flash][cubic_T_flash]") {
    const auto backend = GENERATE(as<std::string>{}, "SRK", "PR");
    const auto fluid = GENERATE(as<std::string>{}, "CO2", "Propane", "Water");
    CAPTURE(backend, fluid);
    std::unique_ptr<AbstractState> AS(AbstractState::factory(backend, fluid));
    std::unique_ptr<AbstractState> ref(AbstractState::factory(backend, fluid));
    auto* cubic = dynamic_cast<AbstractCubicBackend*>(ref.get());
    REQUIRE(cubic != nullptr);
    const double Tc = AS->T_critical(), pc = AS->p_critical();

    SECTION("subcritical gas and liquid") {
        for (double Tr : {0.6, 0.8, 0.95}) {
            const double T = Tr * Tc;
            const double psat = cubic->calc_saturation_ancillary(iP, 0, iT, T);
            ref->update(PT_INPUTS, 0.5 * psat, T);
            check_T_roundtrip(*AS, *ref);
            ref->update(PT_INPUTS, 2 * psat, T);
            check_T_roundtrip(*AS, *ref);
        }
    }
    SECTION("inside the dome") {
        for (double Tr : {0.6, 0.8, 0.95}) {
            const double T = Tr * Tc;
            const double rhoL = cubic->calc_saturation_ancillary(iDmolar, 0, iT, T);
            const double rhoV = cubic->calc_saturation_ancillary(iDmolar, 1, iT, T);
            ref->update(DmolarT_INPUTS, 1 / (0.3 / rhoV + 0.7 / rhoL), T);
            check_T_roundtrip(*AS, *ref);
        }
    }
    SECTION("supercritical temperature") {
        ref->update(PT_INPUTS, 0.5 * pc, 1.2 * Tc);
        check_T_roundtrip(*AS, *ref);
        ref->update(PT_INPUTS, 2 * pc, 1.2 * Tc);
        check_T_roundtrip(*AS, *ref);
    }
}

TEST_CASE("Cubic T-input flashes with an imposed phase", "[cubic][flash][cubic_T_flash]") {
    // With a phase imposed, DHSU_T_flash skips T_phase_determination, so the density seed in
    // solver_for_rho_given_T_oneof_HSU must come from the superancillary hook, not from the
    // classic ancillaries (empty for a cubic: "type not set").
    for (const auto& backend : std::vector<std::string>{"SRK", "PR"}) {
        CAPTURE(backend);
        std::unique_ptr<AbstractState> AS(AbstractState::factory(backend, "Propane"));
        std::unique_ptr<AbstractState> ref(AbstractState::factory(backend, "Propane"));
        const double T = 0.8 * AS->T_critical();
        auto* cubic = dynamic_cast<AbstractCubicBackend*>(ref.get());
        REQUIRE(cubic != nullptr);
        const double psat = cubic->calc_saturation_ancillary(iP, 0, iT, T);
        for (const auto& [p, phase] : std::vector<std::pair<double, phases>>{{2 * psat, iphase_liquid}, {0.5 * psat, iphase_gas}}) {
            CAPTURE(p, phase);
            ref->update(PT_INPUTS, p, T);
            AS->specify_phase(phase);
            REQUIRE_NOTHROW(AS->update(SmolarT_INPUTS, ref->smolar(), T));
            CHECK(AS->rhomolar() == Catch::Approx(ref->rhomolar()).epsilon(1e-8));
            REQUIRE_NOTHROW(AS->update(HmolarT_INPUTS, ref->hmolar(), T));
            CHECK(AS->rhomolar() == Catch::Approx(ref->rhomolar()).epsilon(1e-8));
            AS->unspecify_phase();
        }
    }
}

TEST_CASE("Cubic T,U around T* where the Soave alpha vanishes", "[cubic][flash][cubic_T_flash]") {
    // (du/dv)_T = (a - T da/dT)/((v + D1 b)(v + D2 b)), and for the Soave alpha a - T da/dT vanishes at
    // T* = Tc(1 + 1/m)^2 (COO-145): u does not depend on density there, so T,U must refuse rather than
    // return an arbitrary density (it returned rho ~ 0).  Above T*, u rises with density, and dense
    // states threw "input not in range" because the bracket assumed the opposite.
    std::unique_ptr<AbstractState> AS(AbstractState::factory("SRK", "Water"));
    const double Tc = AS->T_critical(), w = AS->acentric_factor();
    const double m = 0.480 + 1.574 * w - 0.176 * w * w;
    const double Tstar = Tc * (1 + 1 / m) * (1 + 1 / m);
    std::unique_ptr<AbstractState> ref(AbstractState::factory("SRK", "Water"));
    ref->update(DmolarT_INPUTS, 5000, Tstar);
    CHECK_THROWS_WITH(AS->update(TUmolar_INPUTS, Tstar, ref->umolar()), Catch::Matchers::ContainsSubstring("ill-posed"));
    for (double f : {0.99, 1.01, 1.1}) {
        for (double rho : {100.0, 5000.0, 20000.0}) {
            CAPTURE(f, rho);
            ref->update(DmolarT_INPUTS, rho, f * Tstar);
            REQUIRE_NOTHROW(AS->update(TUmolar_INPUTS, f * Tstar, ref->umolar()));
            CHECK(AS->rhomolar() == Catch::Approx(rho).epsilon(1e-8));
        }
    }
}

TEST_CASE("Cubic imposed-phase T-input flash below the table gives the superancillary's error", "[cubic][flash][cubic_T_flash]") {
    // Below the cubic superancillary's low-T end the hook throws a clear range error.  With a phase
    // imposed, solver_for_rho_given_T_oneof_HSU used to swallow it and fall back to the classic
    // ancillaries, which a cubic does not have: "type not set".
    std::unique_ptr<AbstractState> AS(AbstractState::factory("SRK", "Propane"));
    const double T = 0.15 * AS->T_critical();
    AS->specify_phase(iphase_liquid);
    CHECK_THROWS_WITH(AS->update(HmolarT_INPUTS, -20000.0, T), Catch::Matchers::ContainsSubstring("superancillary"));
}

TEST_CASE("Cubic T-input flash rejects densities past the repulsive pole", "[cubic][flash][cubic_T_flash]") {
    // For a cubic, u and h stay finite past rho = 1/(b - c), and above T* u rises with density, so a large
    // enough T,U input used to "solve" beyond the pole with a negative pressure.
    std::unique_ptr<AbstractState> AS(AbstractState::factory("SRK", "Water"));
    auto* cubic = dynamic_cast<AbstractCubicBackend*>(AS.get());
    REQUIRE(cubic != nullptr);
    const double w = AS->acentric_factor(), m = 0.480 + 1.574 * w - 0.176 * w * w;
    const double T = 1.1 * AS->T_critical() * (1 + 1 / m) * (1 + 1 / m);
    const std::vector<double> x = {1.0};
    const double b = cubic->get_cubic()->bm_term(x);
    std::unique_ptr<AbstractState> ref(AbstractState::factory("SRK", "Water"));
    ref->specify_phase(iphase_gas);
    ref->update(DmolarT_INPUTS, 0.999 / b, T);
    const double u_near_pole = ref->umolar();
    ref->update(DmolarT_INPUTS, 0.5 / b, T);
    REQUIRE(u_near_pole > ref->umolar());  // u rises with density here
    // A u that the EOS does reach, but only past the pole (rho = 1.2/b): no physical state has it.
    ref->update(DmolarT_INPUTS, 1.2 / b, T);
    const double u_past_pole = ref->umolar();
    REQUIRE(std::isfinite(u_past_pole));
    REQUIRE(u_past_pole > u_near_pole);
    CHECK_THROWS_WITH(AS->update(TUmolar_INPUTS, T, u_past_pole), Catch::Matchers::ContainsSubstring("No physical state"));
}

TEST_CASE("Cubic P-input flashes at high supercritical pressure", "[cubic][flash][3441]") {
    // At several pc the cubic has three real roots well above Tc, the densest past the repulsive pole,
    // and the P-input bracket runs up to 1.5 Tmax with a phase imposed.  Round trips over that range.
    for (const auto& backend : std::vector<std::string>{"PR", "SRK"}) {
        for (const auto& fluid : std::vector<std::string>{"CO2", "Water"}) {
            std::unique_ptr<AbstractState> AS(AbstractState::factory(backend, fluid)), ref(AbstractState::factory(backend, fluid));
            const double Tc = AS->T_critical(), pc = AS->p_critical();
            for (double pr : {2.0, 5.0}) {
                for (double Tr : {0.9, 1.05, 1.15, 1.5, 3.0}) {
                    const double T = Tr * Tc, p = pr * pc;
                    CAPTURE(backend, fluid, pr, Tr);
                    ref->update(PT_INPUTS, p, T);
                    const double h = ref->hmolar(), s = ref->smolar();
                    REQUIRE(std::isfinite(h));
                    CHECK_NOTHROW(AS->update(HmolarP_INPUTS, h, p));
                    CHECK(AS->T() == Catch::Approx(T).epsilon(1e-8));
                    CHECK_NOTHROW(AS->update(PSmolar_INPUTS, p, s));
                    CHECK(AS->T() == Catch::Approx(T).epsilon(1e-8));
                }
            }
        }
    }
}

#endif
