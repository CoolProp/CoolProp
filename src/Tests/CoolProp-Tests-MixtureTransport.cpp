#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>
#    include "../Backends/Helmholtz/HelmholtzEOSBackend.h"

#    include <cmath>
#    include <memory>
#    include <string>
#    include <vector>

using namespace CoolProp;

namespace {

class InspectableMixtureTransport : public HelmholtzEOSMixtureBackend
{
   public:
    using HelmholtzEOSMixtureBackend::HelmholtzEOSMixtureBackend;

    const std::vector<shared_ptr<HelmholtzEOSMixtureBackend>>& transport_states() const {
        return component_transport_states;
    }

    const std::vector<shared_ptr<HelmholtzEOSMixtureBackend>>& linked_helpers() const {
        return linked_states;
    }

    void create_linked_helpers() {
        add_TPD_state();
        add_critical_state();
        add_transient_pure_state();
    }
};

// Preserve the previous implementation as an equivalence reference, not as
// experimental validation of the approximate mixture transport rules.
double fresh_component_transport(const HelmholtzEOSMixtureBackend& mixture, const std::vector<CoolPropDbl>& fractions, double rho, double T,
                                 parameters property) {
    double sum = 0.0;
    for (std::size_t i = 0; i < fractions.size(); ++i) {
        auto pure = std::make_shared<HelmholtzEOSBackend>(mixture.get_components()[i]);
        pure->update(DmolarT_INPUTS, rho, T);
        const double value = pure->keyed_output(property);
        sum += fractions[i] * (property == iviscosity ? std::log(value) : value);
    }
    return property == iviscosity ? std::exp(sum) : sum;
}

void check_transport(HelmholtzEOSMixtureBackend& mixture, const std::vector<CoolPropDbl>& fractions, double rho, double T,
                     bool conductivity_first = false) {
    mixture.set_mole_fractions(fractions);
    mixture.update(DmolarT_INPUTS, rho, T);
    const auto first = conductivity_first ? iconductivity : iviscosity;
    const auto second = conductivity_first ? iviscosity : iconductivity;
    for (const auto property : {first, second}) {
        CAPTURE(property, rho, T, fractions);
        const double expected = fresh_component_transport(mixture, fractions, rho, T, property);
        CHECK(mixture.keyed_output(property) == Catch::Approx(expected).epsilon(1e-12));
    }
}

}  // namespace

TEST_CASE("Mixture transport agrees with fresh component backends across updates", "[transport][mixture][mixture-transport]") {
    struct Mixture
    {
        std::vector<std::string> names;
        std::vector<CoolPropDbl> fractions;
        double temperature;
    };
    for (const auto& entry : {Mixture{{"Nitrogen", "Oxygen", "Argon", "CarbonDioxide", "Water"}, {0.70, 0.08, 0.01, 0.10, 0.11}, 700.0},
                              Mixture{{"Methane", "Ethane", "Propane"}, {0.7, 0.2, 0.1}, 400.0}, Mixture{{"R125", "R143a"}, {0.5, 0.5}, 400.0},
                              Mixture{{"Nitrogen", "Oxygen"}, {0.8, 0.2}, 110.0}}) {
        CAPTURE(entry.names);
        HelmholtzEOSMixtureBackend mixture(entry.names);
        for (const double rho : {10.0, 200.0, 1000.0, 10.0}) {
            check_transport(mixture, entry.fractions, rho, entry.temperature);
            check_transport(mixture, entry.fractions, rho, entry.temperature + 50.0, true);
        }
        auto changed = entry.fractions;
        changed[0] -= 0.1;
        changed[1] += 0.1;
        check_transport(mixture, changed, 20.0, entry.temperature);
        changed.assign(changed.size(), 0.0);
        changed[0] = 1.0;
        check_transport(mixture, changed, 20.0, entry.temperature, true);
        check_transport(mixture, entry.fractions, 10.0, entry.temperature);
    }
}

TEST_CASE("Excess matrices remain square when the component count changes", "[transport][mixture][mixture-transport]") {
    ExcessTerm excess;
    for (const std::size_t count : {2, 3, 2, 0, 3}) {
        CAPTURE(count);
        excess.resize(count);
        REQUIRE(excess.N == count);
        REQUIRE(excess.F.size() == count);
        REQUIRE(excess.DepartureFunctionMatrix.size() == count);
        for (std::size_t i = 0; i < count; ++i) {
            REQUIRE(excess.F[i].size() == count);
            REQUIRE(excess.DepartureFunctionMatrix[i].size() == count);
        }
    }
}

TEST_CASE("Mixture transport rejects missing or mismatched mole fractions", "[transport][mixture][mixture-transport]") {
    InspectableMixtureTransport mixture(std::vector<std::string>{"Methane", "Ethane"});
    SECTION("mole fractions have not been set") {}
    SECTION("too few mole fractions after replacing components") {
        check_transport(mixture, {0.6, 0.4}, 100.0, 400.0);
        const HelmholtzEOSMixtureBackend replacement(std::vector<std::string>{"Methane", "Ethane", "Propane"});
        mixture.set_components(replacement.get_components());
    }
    SECTION("too many mole fractions after replacing components") {
        const HelmholtzEOSMixtureBackend replacement(std::vector<std::string>{"Methane", "Ethane", "Propane"});
        mixture.set_components(replacement.get_components());
        check_transport(mixture, {0.6, 0.3, 0.1}, 100.0, 400.0);
        const HelmholtzEOSMixtureBackend original(std::vector<std::string>{"Methane", "Ethane"});
        mixture.set_components(original.get_components());
    }
    // Invalidate AbstractState's cached outputs without a flash, which would
    // itself require valid mole fractions before reaching the transport path.
    mixture.clear();
    for (const auto property : {iviscosity, iconductivity}) {
        CHECK_THROWS_WITH(mixture.keyed_output(property),
                          "Mole fractions must be set and match the component count before evaluating mixture transport");
        CHECK(mixture.transport_states().empty());
    }
    const auto count = static_cast<const HelmholtzEOSMixtureBackend&>(mixture).get_components().size();
    check_transport(mixture, std::vector<CoolPropDbl>(count, 1.0 / count), 100.0, 400.0);
}

TEST_CASE("Mixture transport follows replaced components and EOS", "[transport][mixture][mixture-transport]") {
    InspectableMixtureTransport mixture(std::vector<std::string>{"Methane", "Ethane"});
    check_transport(mixture, {0.6, 0.4}, 100.0, 400.0);

    SECTION("component replacement with the same count") {
        const HelmholtzEOSMixtureBackend replacement(std::vector<std::string>{"Nitrogen", "Oxygen"});
        mixture.set_components(replacement.get_components());
        check_transport(mixture, {0.6, 0.4}, 100.0, 400.0, true);
    }
    SECTION("component replacement with a different count") {
        const HelmholtzEOSMixtureBackend replacement(std::vector<std::string>{"Methane", "Ethane", "Propane"});
        mixture.set_components(replacement.get_components());
        check_transport(mixture, {0.6, 0.3, 0.1}, 100.0, 400.0);
    }
    SECTION("component EOS replacement") {
        mixture.change_EOS(0, "SRK");
        CHECK(mixture.transport_states().empty());
        check_transport(mixture, {0.6, 0.4}, 1000.0, 400.0, true);
    }
    SECTION("named reference state change") {
        mixture.set_reference_stateS("NBP");
        CHECK(mixture.transport_states().empty());
        check_transport(mixture, {0.6, 0.4}, 100.0, 400.0);
    }
    SECTION("custom reference state change") {
        mixture.set_reference_stateD(400.0, 100.0, 1000.0, 10.0);
        CHECK(mixture.transport_states().empty());
        check_transport(mixture, {0.6, 0.4}, 100.0, 400.0);
    }
    SECTION("independent backend and get_copy") {
        HelmholtzEOSMixtureBackend other(std::vector<std::string>{"Methane", "Ethane"});
        std::unique_ptr<HelmholtzEOSMixtureBackend> copy;
        {
            HelmholtzEOSMixtureBackend source(std::vector<std::string>{"Methane", "Ethane"});
            check_transport(source, {0.6, 0.4}, 100.0, 400.0);
            copy.reset(source.get_copy());
        }
        check_transport(other, {0.2, 0.8}, 200.0, 450.0);
        check_transport(*copy, {0.4, 0.6}, 400.0, 500.0, true);
        check_transport(mixture, {0.6, 0.4}, 100.0, 400.0);
    }
}

TEST_CASE("Component replacement discards old linked helper models", "[transport][mixture][mixture-transport]") {
    InspectableMixtureTransport mixture(std::vector<std::string>{"Methane", "Ethane"});
    check_transport(mixture, {0.6, 0.4}, 100.0, 400.0);
    mixture.create_linked_helpers();
    const auto old_helpers = mixture.linked_helpers();
    REQUIRE(old_helpers.size() == 5);
    const HelmholtzEOSMixtureBackend replacement(std::vector<std::string>{"Methane", "Ethane", "Propane"});

    SECTION("replacement with saturation helpers") {
        mixture.set_components(replacement.get_components());
        REQUIRE(mixture.linked_helpers().size() == 2);
        check_transport(mixture, {0.6, 0.3, 0.1}, 100.0, 400.0);
        // With an old two-component helper still linked, this write passes
        // its updated N check but indexes past its unchanged 2x2 matrix.
        REQUIRE_NOTHROW(mixture.set_binary_interaction_double(0, 2, "Fij", 1.0));
        CHECK(mixture.SatL->get_binary_interaction_double(0, 2, "Fij") == 1.0);
        CHECK(mixture.SatV->get_binary_interaction_double(0, 2, "Fij") == 1.0);
        mixture.create_linked_helpers();
        REQUIRE(mixture.linked_helpers().size() == 5);
        for (const auto& helper : mixture.linked_helpers()) {
            CHECK(static_cast<const HelmholtzEOSMixtureBackend&>(*helper).get_components().size() == 3);
        }
    }
    SECTION("replacement without saturation helpers") {
        mixture.set_components(replacement.get_components(), false);
        CHECK(mixture.linked_helpers().empty());
        CHECK_FALSE(mixture.SatL);
        CHECK_FALSE(mixture.SatV);
    }
    for (const auto& helper : old_helpers) {
        CHECK(static_cast<const HelmholtzEOSMixtureBackend&>(*helper).get_components().size() == 2);
        CHECK_THROWS(helper->set_binary_interaction_double(0, 2, "Fij", 1.0));
    }
}

TEST_CASE("Invalid mole fractions discard already warmed transport helpers", "[transport][mixture][mixture-transport]") {
    for (const auto property : {iviscosity, iconductivity}) {
        for (const std::size_t count : {0, 1, 3}) {
            CAPTURE(property, count);
            InspectableMixtureTransport mixture(std::vector<std::string>{"Methane", "Ethane"});
            check_transport(mixture, {0.6, 0.4}, 100.0, 400.0);
            const auto old_helpers = mixture.transport_states();
            REQUIRE(old_helpers.size() == 2);
            mixture.get_mole_fractions_ref().resize(count);
            mixture.clear();
            REQUIRE(mixture.transport_states() == old_helpers);
            CHECK_THROWS_WITH(mixture.keyed_output(property),
                              "Mole fractions must be set and match the component count before evaluating mixture transport");
            CHECK(mixture.transport_states().empty());
            check_transport(mixture, {0.6, 0.4}, 100.0, 400.0);
            CHECK(mixture.transport_states() != old_helpers);
        }
    }
}

TEST_CASE("Mixture transport reuses helpers without linking their component counts", "[transport][mixture][mixture-transport]") {
    InspectableMixtureTransport mixture(std::vector<std::string>{"Methane", "Ethane"});
    REQUIRE(mixture.transport_states().empty());
    check_transport(mixture, {0.6, 0.4}, 100.0, 400.0);
    const auto original = mixture.transport_states();
    REQUIRE(original.size() == 2);
    REQUIRE(original[0]);
    REQUIRE(original[1]);

    CHECK(mixture.fluid_param_string("CAS") == "74-82-8");
    CHECK(mixture.transport_states() == original);
    mixture.clear();
    check_transport(mixture, {0.2, 0.8}, 200.0, 450.0, true);
    CHECK(mixture.transport_states() == original);
    mixture.set_binary_interaction_double(0, 1, "betaT", 1.01);
    check_transport(mixture, {0.6, 0.4}, 100.0, 400.0);
    CHECK(mixture.transport_states() == original);
    for (const auto& pure : original) {
        CHECK(pure->get_mole_fractions().size() == 1);
    }
    InspectableMixtureTransport other(std::vector<std::string>{"Methane", "Ethane"});
    check_transport(other, {0.6, 0.4}, 100.0, 400.0);
    CHECK(other.transport_states()[0] != original[0]);
    CHECK(other.transport_states()[1] != original[1]);
}

TEST_CASE("Mixture transport handles persistent mutable component aliases", "[transport][mixture][mixture-transport]") {
    InspectableMixtureTransport mixture(std::vector<std::string>{"Methane", "Ethane"});
    check_transport(mixture, {0.6, 0.4}, 100.0, 400.0);
    auto& alias = mixture.get_components();
    CHECK(mixture.transport_states().empty());

    SECTION("alias retained across repeated evaluations") {
        check_transport(mixture, {0.6, 0.4}, 100.0, 400.0);
        const auto previous = mixture.transport_states();
        check_transport(mixture, {0.6, 0.4}, 100.0, 400.0);
        CHECK(mixture.transport_states() != previous);
    }
    SECTION("vector alias retained across component replacement") {
        const HelmholtzEOSMixtureBackend replacement(std::vector<std::string>{"Nitrogen", "Oxygen"});
        mixture.set_components(replacement.get_components());
        check_transport(mixture, {0.6, 0.4}, 100.0, 400.0);
    }

    alias[1].transport.viscosity_model_provided = false;
    mixture.update(DmolarT_INPUTS, 100.0, 400.0);
    CHECK_THROWS(mixture.viscosity());
    CHECK(mixture.transport_states().empty());
    alias[1].transport.viscosity_model_provided = true;
    check_transport(mixture, {0.6, 0.4}, 100.0, 400.0);

    alias[1].transport.conductivity_model_provided = false;
    mixture.update(DmolarT_INPUTS, 100.0, 400.0);
    CHECK_THROWS(mixture.conductivity());
    alias[1].transport.conductivity_model_provided = true;
    check_transport(mixture, {0.6, 0.4}, 100.0, 400.0, true);
}

TEST_CASE("Mixture transport discards helpers after failed evaluation or model mutation", "[transport][mixture][mixture-transport]") {
    InspectableMixtureTransport mixture(std::vector<std::string>{"Methane", "Ethane"});
    check_transport(mixture, {0.6, 0.4}, 100.0, 400.0);
    const auto original = mixture.transport_states();
    REQUIRE(original.size() == 2);

    SECTION("transport failure after an earlier component succeeds") {
        original[1]->get_components()[0].transport.viscosity_model_provided = false;
        mixture.update(DmolarT_INPUTS, 100.0, 400.0);
        CHECK_THROWS(mixture.viscosity());
    }
    SECTION("invalid EOS name") {
        CHECK_THROWS(mixture.change_EOS(0, "invalid"));
    }
    SECTION("invalid reference state") {
        CHECK_THROWS(mixture.set_reference_stateS("invalid"));
    }
    CHECK(mixture.transport_states().empty());
    check_transport(mixture, {0.6, 0.4}, 100.0, 400.0, true);
    CHECK(mixture.transport_states() != original);
}

TEST_CASE("Mixture transport timing", "[mixture-transport][!benchmark]") {
    const std::vector<CoolPropDbl> fractions{0.70, 0.08, 0.01, 0.10, 0.11};
    HelmholtzEOSMixtureBackend mixture(std::vector<std::string>{"Nitrogen", "Oxygen", "Argon", "CarbonDioxide", "Water"});
    mixture.set_mole_fractions(fractions);
    // Each iteration updates the state, invalidating AbstractState's property
    // cache. Both paths use identical update calls and the same mixing rules.
    mixture.update(DmolarT_INPUTS, 20.0, 700.0);
    const double expected = fresh_component_transport(mixture, fractions, 20.0, 700.0, iviscosity)
                            + fresh_component_transport(mixture, fractions, 20.0, 700.0, iconductivity);
    REQUIRE(mixture.viscosity() + mixture.conductivity() == Catch::Approx(expected).epsilon(1e-12));

    BENCHMARK("fresh components: rho,T update + viscosity + conductivity") {
        mixture.update(DmolarT_INPUTS, 20.0, 700.0);
        return fresh_component_transport(mixture, fractions, 20.0, 700.0, iviscosity)
               + fresh_component_transport(mixture, fractions, 20.0, 700.0, iconductivity);
    };
    BENCHMARK("reused components: rho,T update + viscosity + conductivity") {
        mixture.update(DmolarT_INPUTS, 20.0, 700.0);
        return mixture.viscosity() + mixture.conductivity();
    };
    // Also include phase determination/density solving, as in a PT-based caller.
    mixture.update(PT_INPUTS, 101325.0, 700.0);
    BENCHMARK("fresh components: p,T update + viscosity + conductivity") {
        mixture.update(PT_INPUTS, 101325.0, 700.0);
        return fresh_component_transport(mixture, fractions, mixture.rhomolar(), mixture.T(), iviscosity)
               + fresh_component_transport(mixture, fractions, mixture.rhomolar(), mixture.T(), iconductivity);
    };
    BENCHMARK("reused components: p,T update + viscosity + conductivity") {
        mixture.update(PT_INPUTS, 101325.0, 700.0);
        return mixture.viscosity() + mixture.conductivity();
    };
}

#endif
