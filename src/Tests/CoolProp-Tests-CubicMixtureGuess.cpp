// Density solve of cubic mixtures (PR, SRK) with a guess density and no imposed phase (GitHub #3597).
//
// With three density roots AbstractCubicBackend::solver_rho_Tp used to ignore the guess and call p_critical(), which
// for a mixture is a full critical-point search.  The recovery pass of the mixture PT flash solves both phases that
// way, so a PT flash of a superheated gas took from about 0.1 s (2 components) to several seconds (6) wherever the
// feed also has a liquid root, and a split was lost where the critical-point search threw.

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/Configuration.h"
#    include "CoolProp/DataStructures.h"
#    include "../Backends/Cubics/CubicBackend.h"

#    include <cmath>
#    include <limits>
#    include <memory>
#    include <string>
#    include <type_traits>
#    include <utility>
#    include <vector>

using namespace CoolProp;

namespace {

// A cubic backend that counts its critical-point searches, including those of the states derived from it (SatL, SatV,
// the critical-point and transient states), which share the counter through get_copy().
template <class Backend>
class CountingCubic final : public Backend
{
   private:
    std::shared_ptr<int> m_searches;

    // Virtual dispatch does not reach get_copy() below from Backend's constructor, which would create plain Backend
    // saturation states; create them here instead, as AbstractCubicBackend::setup does.
    void add_saturation_states() {
        this->SatL.reset(this->get_copy(false));
        this->SatL->specify_phase(iphase_liquid);
        this->linked_states.push_back(this->SatL);
        this->SatV.reset(this->get_copy(false));
        this->SatV->specify_phase(iphase_gas);
        this->linked_states.push_back(this->SatV);
    }

   public:
    CountingCubic(const std::vector<std::string>& fluids, std::shared_ptr<int> searches)
      : Backend(fluids, get_config_double(R_U_CODATA), false), m_searches(std::move(searches)) {
        add_saturation_states();
    }
    CountingCubic(const std::vector<double>& Tc, const std::vector<double>& pc, const std::vector<double>& acentric, double R_u,
                  std::shared_ptr<int> searches, bool generate_SatL_and_SatV)
      : Backend(Tc, pc, acentric, R_u, false), m_searches(std::move(searches)) {
        if (generate_SatL_and_SatV) add_saturation_states();
    }
    HelmholtzEOSMixtureBackend* get_copy(bool generate_SatL_and_SatV = true) override {
        auto& c = this->get_cubic();
        auto copy = std::make_unique<CountingCubic>(c->get_Tc(), c->get_pc(), c->get_acentric(), c->get_R_u(), m_searches, generate_SatL_and_SatV);
        copy->copy_internals(*this);
        return copy.release();
    }
    std::vector<CriticalState> calc_all_critical_points() override {
        ++*m_searches;
        return Backend::calc_all_critical_points();
    }
};

struct Mixture
{
    std::vector<std::string> fluids;
    std::vector<double> z;
};

struct TwoRootState
{
    Mixture mix;
    double T;
};

std::string joined(const std::vector<std::string>& fluids) {
    std::string s;
    for (const auto& f : fluids) {
        if (!s.empty()) s += "&";
        s += f;
    }
    return s;
}

std::unique_ptr<AbstractState> make_state(const std::string& backend, const Mixture& m) {
    std::unique_ptr<AbstractState> AS(AbstractState::factory(backend, joined(m.fluids)));
    AS->set_mole_fractions(m.z);
    return AS;
}

// Density of the branch selected by an imposed phase
double branch_density(const std::string& backend, const Mixture& m, phases phase, double p, double T) {
    auto AS = make_state(backend, m);
    AS->specify_phase(phase);
    AS->update(PT_INPUTS, p, T);
    return AS->rhomolar();
}

HelmholtzEOSMixtureBackend& as_heos(AbstractState& AS) {
    auto* p = dynamic_cast<HelmholtzEOSMixtureBackend*>(&AS);
    REQUIRE(p != nullptr);
    return *p;
}

Mixture nitrogen_methane() {
    return {{"Nitrogen", "Methane"}, {0.3, 0.7}};
}

// Superheated-gas states at 1 bar whose feed also has a (metastable) liquid root
std::vector<TwoRootState> two_root_gas_states() {
    const Mixture ternary = {{"Nitrogen", "Methane", "Ethane"}, {0.3, 0.4, 0.3}};
    const Mixture quaternary = {{"Nitrogen", "Methane", "Ethane", "Propane"}, {0.3, 0.4, 0.2, 0.1}};
    return {{nitrogen_methane(), 108.5}, {nitrogen_methane(), 127.4}, {nitrogen_methane(), 147.4}, {ternary, 180.0}, {quaternary, 193.0}};
}

}  // namespace

TEST_CASE("Cubic mixture: the density solve follows its guess when no phase is imposed", "[cubic][cubic_mixture_guess]") {
    const std::string backend = GENERATE("PR", "SRK");
    const double p = 1e5;
    for (const auto& state : two_root_gas_states()) {
        const double T = state.T;
        const double rho_gas = branch_density(backend, state.mix, iphase_gas, p, T);
        const double rho_liq = branch_density(backend, state.mix, iphase_liquid, p, T);
        CAPTURE(backend, joined(state.mix.fluids), T, rho_gas, rho_liq);
        REQUIRE(rho_liq > 5 * rho_gas);  // two distinct roots, or the checks below are vacuous
        auto AS = make_state(backend, state.mix);
        auto& eos = as_heos(*AS);
        for (const double factor : {0.5, 0.95, 1.0, 1.05, 2.0}) {
            CAPTURE(factor);
            eos.update_TP_guessrho(T, p, factor * rho_gas);
            CHECK(eos.rhomolar() == Catch::Approx(rho_gas).epsilon(1e-12));
            CHECK(eos.p() == Catch::Approx(p).epsilon(1e-6));
            eos.update_TP_guessrho(T, p, factor * rho_liq);
            CHECK(eos.rhomolar() == Catch::Approx(rho_liq).epsilon(1e-12));
            CHECK(eos.p() == Catch::Approx(p).epsilon(1e-6));
        }
        // A guess between the two roots never returns the mechanically unstable middle root
        for (const double w : {0.1, 0.3, 0.5, 0.7, 0.9}) {
            CAPTURE(w);
            eos.update_TP_guessrho(T, p, std::exp((1 - w) * std::log(rho_gas) + w * std::log(rho_liq)));
            const bool on_gas = eos.rhomolar() == Catch::Approx(rho_gas).epsilon(1e-12);
            const bool on_liq = eos.rhomolar() == Catch::Approx(rho_liq).epsilon(1e-12);
            CHECK((on_gas || on_liq));
        }
        // Extreme but finite guesses: squaring them underflows or overflows
        eos.update_TP_guessrho(T, p, 1e-200);
        CHECK(eos.rhomolar() == Catch::Approx(rho_gas).epsilon(1e-12));
        eos.update_TP_guessrho(T, p, 1e200);
        CHECK(eos.rhomolar() == Catch::Approx(rho_liq).epsilon(1e-12));
        // An imposed phase still wins over the guess
        AS->specify_phase(iphase_gas);
        eos.update_TP_guessrho(T, p, rho_liq);
        CHECK(eos.rhomolar() == Catch::Approx(rho_gas).epsilon(1e-12));
        AS->specify_phase(iphase_liquid);
        eos.update_TP_guessrho(T, p, rho_gas);
        CHECK(eos.rhomolar() == Catch::Approx(rho_liq).epsilon(1e-12));
    }
}

TEST_CASE("Cubic mixture: an unusable guess is treated as no guess", "[cubic][cubic_mixture_guess]") {
    // Without a guess the root is the gas one above the dew point and the liquid one below the bubble point (86.5 K).
    // Both states are needed: a followed -1 or 0 always selects the gas root and a followed NaN or infinity the liquid
    // one, so each is wrong in only one of them.
    const std::string backend = GENERATE("PR", "SRK");
    const double p = 1e5;
    const Mixture mix = nitrogen_methane();
    auto AS = make_state(backend, mix);
    for (const bool liquid : {false, true}) {
        const double T = liquid ? 75.0 : 127.4;
        const double rho_gas = branch_density(backend, mix, iphase_gas, p, T);
        const double rho_liq = branch_density(backend, mix, iphase_liquid, p, T);
        REQUIRE(rho_liq > 5 * rho_gas);
        for (const double guess : {-1.0, 0.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()}) {
            CAPTURE(backend, T, guess);
            as_heos(*AS).update_TP_guessrho(T, p, guess);
            CHECK(AS->rhomolar() == Catch::Approx(liquid ? rho_liq : rho_gas).epsilon(1e-12));
        }
    }
}

TEST_CASE("Cubic mixture: the physical domain of a translated cubic ends at 1/(b - c)", "[cubic][cubic_mixture_guess]") {
    // The translation is large enough that the liquid is denser than 1/b, so a root filter that forgot c would drop it
    const std::string backend = GENERATE("PR", "SRK");
    const double p = 1e5, T = 75;
    const Mixture mix = nitrogen_methane();
    auto translated = [&]() {
        auto state = make_state(backend, mix);
        auto* ACB = dynamic_cast<AbstractCubicBackend*>(state.get());
        REQUIRE(ACB != nullptr);
        const double b = ACB->get_cubic()->bm_term(mix.z);
        for (std::size_t i = 0; i < mix.z.size(); ++i) {
            state->set_fluid_parameter_double(i, "cm", 0.2 * b);
        }
        return state;
    };
    std::vector<double> rho;
    for (const phases ph : {iphase_gas, iphase_liquid}) {
        auto state = translated();
        state->specify_phase(ph);
        state->update(PT_INPUTS, p, T);
        rho.push_back(state->rhomolar());
    }
    auto AS = translated();
    auto* ACB = dynamic_cast<AbstractCubicBackend*>(AS.get());
    REQUIRE(ACB != nullptr);
    const double b = ACB->get_cubic()->bm_term(mix.z), c = ACB->get_cubic()->cm_term(mix.z);
    CAPTURE(backend, rho[0], rho[1], b, c);
    REQUIRE(rho[1] > 5 * rho[0]);
    REQUIRE(rho[1] * b > 1);
    REQUIRE(rho[1] * (b - c) < 1);
    auto& eos = as_heos(*AS);
    for (const double factor : {0.9, 1.0, 1.1}) {
        CAPTURE(factor);
        eos.update_TP_guessrho(T, p, factor * rho[0]);
        CHECK(eos.rhomolar() == Catch::Approx(rho[0]).epsilon(1e-12));
        eos.update_TP_guessrho(T, p, factor * rho[1]);
        CHECK(eos.rhomolar() == Catch::Approx(rho[1]).epsilon(1e-12));
    }
}

TEST_CASE("Cubic mixture: a guess does not select a root at negative pressure", "[cubic][cubic_mixture_guess]") {
    // Of the two physical roots there the lighter is mechanically unstable; the solve refuses rather than choose
    const std::string backend = GENERATE("PR", "SRK");
    const double p = -1e5, T = 75;
    auto AS = make_state(backend, nitrogen_methane());
    auto* ACB = dynamic_cast<AbstractCubicBackend*>(AS.get());
    REQUIRE(ACB != nullptr);
    int Nsoln = 0;
    double r0 = 0, r1 = 0, r2 = 0;
    ACB->rho_Tp_cubic(T, p, Nsoln, r0, r1, r2);
    CAPTURE(backend, Nsoln, r0, r1, r2);
    REQUIRE(Nsoln == 3);
    REQUIRE(r0 < 0);
    REQUIRE(r1 > 0);
    for (const double guess : {r1, r2}) {
        CHECK_THROWS_WITH(as_heos(*AS).update_TP_guessrho(T, p, guess), Catch::Matchers::ContainsSubstring("2 physical density roots"));
    }
}

TEST_CASE("Pure cubic: a guess does not override the choice of the stable root", "[cubic][cubic_mixture_guess]") {
    // A pure SRK/PR picks its root from the superancillary's saturation pressure; only mixtures follow the guess
    const std::string backend = GENERATE("PR", "SRK");
    const Mixture propane = {{"n-Propane"}, {1.0}};
    const double T = 300;
    auto sat = make_state(backend, propane);
    sat->update(QT_INPUTS, 0, T);
    const double psat = sat->p();
    for (const double ratio : {0.9, 1.1}) {
        const double p = ratio * psat;
        const double rho_gas = branch_density(backend, propane, iphase_gas, p, T);
        const double rho_liq = branch_density(backend, propane, iphase_liquid, p, T);
        CAPTURE(backend, ratio, rho_gas, rho_liq);
        REQUIRE(rho_liq > 5 * rho_gas);
        auto AS = make_state(backend, propane);
        for (const double guess : {rho_gas, rho_liq}) {
            CAPTURE(guess);
            as_heos(*AS).update_TP_guessrho(T, p, guess);
            CHECK(AS->rhomolar() == Catch::Approx(ratio < 1 ? rho_gas : rho_liq).epsilon(1e-12));
        }
    }
}

TEST_CASE("Cubic mixture: a guess does not change a state with a single root", "[cubic][cubic_mixture_guess]") {
    const std::string backend = GENERATE("PR", "SRK");
    const double p = 1e5, T = 200;
    const Mixture mix = nitrogen_methane();
    const double rho = branch_density(backend, mix, iphase_gas, p, T);
    CAPTURE(backend, rho);
    REQUIRE(branch_density(backend, mix, iphase_liquid, p, T) == Catch::Approx(rho).epsilon(1e-12));
    auto AS = make_state(backend, mix);
    for (const double guess : {0.5 * rho, rho, 2 * rho, 30000.0}) {
        CAPTURE(guess);
        as_heos(*AS).update_TP_guessrho(T, p, guess);
        CHECK(AS->rhomolar() == Catch::Approx(rho).epsilon(1e-12));
    }
}

TEST_CASE("Cubic mixture: with three real roots but one physical, every guess returns the physical one", "[cubic][cubic_mixture_guess]") {
    // Hot and compressed, the cubic in density can have three real roots of which only one lies in the physical
    // domain; the others are negative or past the 1/(b - c) pole and must not be taken for a liquid.
    const std::string backend = GENERATE("PR", "SRK");
    const Mixture mix = nitrogen_methane();
    int n_three_root_states = 0;
    for (const double T : {500.0, 800.0, 1200.0}) {
        for (const double p : {1e7, 3e7}) {
            auto AS = make_state(backend, mix);
            auto* ACB = dynamic_cast<AbstractCubicBackend*>(AS.get());
            REQUIRE(ACB != nullptr);
            int Nsoln = 0;
            double r0 = 0, r1 = 0, r2 = 0;
            ACB->rho_Tp_cubic(T, p, Nsoln, r0, r1, r2);
            n_three_root_states += (Nsoln == 3) ? 1 : 0;
            const double b_minus_c = ACB->get_cubic()->bm_term(mix.z) - ACB->get_cubic()->cm_term(mix.z);
            const double rho_ideal = p / (8.314462618 * T);
            for (const double guess : {0.1 * rho_ideal, rho_ideal, 10 * rho_ideal, 1e6}) {
                CAPTURE(backend, T, p, Nsoln, guess);
                as_heos(*AS).update_TP_guessrho(T, p, guess);
                CHECK(AS->rhomolar() > 0);
                CHECK(AS->rhomolar() * b_minus_c < 1);
                CHECK(AS->p() == Catch::Approx(p).epsilon(1e-8));
            }
        }
    }
    CHECK(n_three_root_states >= 1);  // the grid must reach the regime this test is about
}

TEMPLATE_TEST_CASE("Cubic mixture: PT and HP flashes of a superheated gas do not search for critical points", "[cubic][cubic_mixture_guess]",
                   SRKBackend, PengRobinsonBackend) {
    const std::string backend = std::is_same_v<TestType, SRKBackend> ? "SRK" : "PR";
    const double p = 1e5;
    for (const auto& state : two_root_gas_states()) {
        const double T = state.T;
        const double rho_gas = branch_density(backend, state.mix, iphase_gas, p, T);
        const double rho_liq = branch_density(backend, state.mix, iphase_liquid, p, T);
        CAPTURE(backend, joined(state.mix.fluids), T, rho_gas, rho_liq);
        REQUIRE(rho_liq > 5 * rho_gas);

        auto searches = std::make_shared<int>(0);
        CountingCubic<TestType> AS(state.mix.fluids, searches);
        AS.set_mole_fractions(state.mix.z);

        // The density solve itself, on a sub-state without an imposed phase, as the recovery pass of the PT flash does
        // it.  Does not depend on the flash taking that pass.
        std::unique_ptr<HelmholtzEOSMixtureBackend> sub(AS.get_copy(false));
        sub->set_mole_fractions(state.mix.z);
        sub->update_TP_guessrho(T, p, rho_liq);
        sub->update_TP_guessrho(T, p, rho_gas);
        CHECK(*searches == 0);

        *searches = 0;
        AS.update(PT_INPUTS, p, T);
        CHECK(*searches == 0);
        CHECK(AS.rhomolar() == Catch::Approx(rho_gas).epsilon(1e-10));
        CHECK(AS.phase() != iphase_twophase);

        // The temperature search of the HP flash runs through the same region
        const double h = AS.hmolar();
        *searches = 0;
        AS.update(HmolarP_INPUTS, h, p);
        CHECK(*searches == 0);
        CHECK(AS.T() == Catch::Approx(T).epsilon(1e-8));

        // The counter sees a search, also one made by a sub-state: otherwise the checks above could pass because
        // nothing is counted
        *searches = 0;
        (void)AS.p_critical();
        CHECK(*searches == 1);
        (void)sub->p_critical();
        CHECK(*searches == 2);
    }
}

#endif  // ENABLE_CATCH
