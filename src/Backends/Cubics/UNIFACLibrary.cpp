#include "UNIFACLibrary.h"
#include "Backends/Helmholtz/Fluids/FluidLibrary.h"
#include "CoolProp/Configuration.h"

namespace UNIFACLibrary {

void UNIFACParameterLibrary::jsonize(std::string& s, nlohmann::json& d) {
    d = cpjson::parse(s);
}
void UNIFACParameterLibrary::populate(const nlohmann::json& group_data, const nlohmann::json& interaction_data, const nlohmann::json& comp_data) {
    if (CoolProp::get_config_bool(VTPR_ALWAYS_RELOAD_LIBRARY)) {
        groups.clear();
        interaction_parameters.clear();
        components.clear();
    }
    // Schema should have been used to validate the data already, so by this point we are can safely consume the data without checking ...
    for (const auto& el : group_data) {
        Group g;
        g.sgi = el.at("sgi").get<int>();
        g.mgi = el.at("mgi").get<int>();
        g.R_k = el.at("R_k").get<double>();
        g.Q_k = el.at("Q_k").get<double>();
        groups.push_back(g);
    }
    for (const auto& el : interaction_data) {
        InteractionParameters ip;
        ip.mgi1 = el.at("mgi1").get<int>();
        ip.mgi2 = el.at("mgi2").get<int>();
        ip.a_ij = el.at("a_ij").get<double>();
        ip.a_ji = el.at("a_ji").get<double>();
        ip.b_ij = el.at("b_ij").get<double>();
        ip.b_ji = el.at("b_ji").get<double>();
        ip.c_ij = el.at("c_ij").get<double>();
        ip.c_ji = el.at("c_ji").get<double>();
        interaction_parameters.push_back(ip);
    }
    for (const auto& el : comp_data) {
        Component c;
        c.inchikey = el.at("inchikey").get<std::string>();
        c.registry_number = el.at("registry_number").get<std::string>();
        c.name = el.at("name").get<std::string>();
        c.Tc = el.at("Tc").get<double>();
        c.pc = el.at("pc").get<double>();
        c.acentric = el.at("acentric").get<double>();
        c.molemass = el.at("molemass").get<double>();
        // userid is an optional user identifier
        if (el.contains("userid")) {
            c.userid = el.at("userid").get<std::string>();
        }
        // If provided, store information about the alpha function in use
        if (el.contains("alpha") && el.at("alpha").is_object()) {
            const nlohmann::json& alpha = el.at("alpha");
            c.alpha_type = cpjson::get_string(alpha, "type");
            c.alpha_coeffs = cpjson::get_double_array(alpha, "c");
        } else {
            c.alpha_type = "default";
        }
        if (el.contains("alpha0") && el.at("alpha0").is_array()) {
            c.alpha0 = CoolProp::JSONFluidLibrary::parse_alpha0(el.at("alpha0"));
        }
        const nlohmann::json& comp_groups = el.at("groups");
        for (const auto& g : comp_groups) {
            int count = g.at("count").get<int>();
            int sgi = g.at("sgi").get<int>();
            if (has_group(sgi)) {
                ComponentGroup cg(count, get_group(sgi));
                c.groups.push_back(cg);
            }
        }
        components.push_back(c);
    }
}
void UNIFACParameterLibrary::populate(std::string& group_data, std::string& interaction_data, std::string& decomp_data) {
    nlohmann::json group_JSON, interaction_JSON, decomp_JSON;
    jsonize(group_data, group_JSON);
    jsonize(interaction_data, interaction_JSON);
    jsonize(decomp_data, decomp_JSON);
    populate(group_JSON, interaction_JSON, decomp_JSON);
    m_populated = true;
}
Group UNIFACParameterLibrary::get_group(int sgi) const {
    for (const auto& group : groups) {
        if (group.sgi == sgi) {
            return group;
        }
    }
    throw CoolProp::ValueError("Could not find group");
}
bool UNIFACParameterLibrary::has_group(int sgi) const {
    for (const auto& group : groups) {
        if (group.sgi == sgi) {
            return true;
        }
    }
    return false;
}

InteractionParameters UNIFACParameterLibrary::get_interaction_parameters(int mgi1, int mgi2) const {

    // If both mgi are the same, yield all zeros for the interaction parameters
    if (mgi1 == mgi2) {
        InteractionParameters ip;
        ip.mgi1 = mgi1;
        ip.mgi2 = mgi2;
        ip.zero_out();
        return ip;
    }
    for (const auto& interaction_parameter : interaction_parameters) {
        if (interaction_parameter.mgi1 == mgi1 && interaction_parameter.mgi2 == mgi2) {
            // Correct order, return it
            return interaction_parameter;
        }
        if (interaction_parameter.mgi2 == mgi1 && interaction_parameter.mgi1 == mgi2) {
            // Backwards, swap the parameters
            InteractionParameters ip = interaction_parameter;
            ip.swap();
            return ip;
        }
    }
    throw CoolProp::ValueError(format("Could not find interaction between pair mgi[%d]-mgi[%d]", static_cast<int>(mgi1), static_cast<int>(mgi2)));
}

Component UNIFACParameterLibrary::get_component(const std::string& identifier, const std::string& value) const {
    if (identifier == "name") {
        for (const auto& component : components) {
            if (component.name == value) {
                return component;
            }
        }
    }
    throw CoolProp::ValueError(format("Could not find component: %s with identifier: %s", value.c_str(), identifier.c_str()));
}

}; /* namespace UNIFACLibrary */

#if defined(ENABLE_CATCH)
#    include <catch2/catch_all.hpp>

#    include "UNIFAC.h"

TEST_CASE("Check Poling example for UNIFAC", "[UNIFAC]") {
    std::string acetone_pentane_groups =
      "[{ \"Tc\": 508.1, \"acentric\": 0.3071, \"groups\": [ { \"count\": 1,  \"sgi\": 1 },  {\"count\": 1, \"sgi\": 18 } ],  \"molemass\": 0.44, "
      "\"inchikey\": \"?????????????\",  \"name\": \"Acetone\", \"pc\": 4700000.0, \"registry_number\": \"67-64-1\", \"userid\": \"\" },  { \"Tc\": "
      "469.7000000000001,  \"acentric\": 0.251,  \"molemass\": 0.44,    \"groups\": [ { \"count\": 2, \"sgi\": 1 }, { \"count\": 3, \"sgi\": 2 } ],  "
      "\"inchikey\": \"?????????????\", \"name\": \"n-Pentane\", \"pc\": 3370000.0,  \"registry_number\": \"109-66-0\",  \"userid\": \"\" } ]";
    std::string groups = "[{\"Q_k\": 0.848, \"R_k\": 0.9011, \"maingroup_name\": \"CH2\", \"mgi\": 1, \"sgi\": 1, \"subgroup_name\": \"CH3\"},"
                         "{\"Q_k\": 0.540, \"R_k\": 0.6744, \"maingroup_name\": \"CH2\", \"mgi\": 1, \"sgi\": 2, \"subgroup_name\": \"CH2\"},"
                         "{\"Q_k\": 1.488, \"R_k\": 1.6724, \"maingroup_name\": \"CH2CO\", \"mgi\": 9, \"sgi\": 18, \"subgroup_name\": \"CH3CO\"}]";
    std::string interactions = R"([{"a_ij": 476.4, "a_ji": 26.76, "b_ij": 0.0, "b_ji": 0.0,  "c_ij": 0.0, "c_ji": 0.0, "mgi1": 1, "mgi2": 9}])";

    SECTION("Validate AC for acetone + n-pentane") {
        UNIFACLibrary::UNIFACParameterLibrary lib;
        CHECK_NOTHROW(lib.populate(groups, interactions, acetone_pentane_groups));
        UNIFAC::UNIFACMixture mix(lib, 1.0);
        std::vector<std::string> names;
        names.emplace_back("Acetone");
        names.emplace_back("n-Pentane");
        mix.set_components("name", names);
        mix.set_interaction_parameters();

        std::vector<double> z(2, 0.047);
        z[1] = 1 - z[0];
        mix.set_mole_fractions(z);
        CHECK_NOTHROW(mix.set_temperature(307));

        double lngammaR0 = mix.ln_gamma_R(1.0 / 307, 0, 0);
        double lngammaR1 = mix.ln_gamma_R(1.0 / 307, 1, 0);
        CAPTURE(lngammaR0);
        CAPTURE(lngammaR1);
        CHECK(std::abs(lngammaR0 - 1.66) < 1e-2);
        CHECK(std::abs(lngammaR1 - 5.68e-3) < 1e-3);

        std::vector<double> gamma(2);
        mix.activity_coefficients(1.0 / 307, z, gamma);
        CAPTURE(gamma[0]);
        CAPTURE(gamma[1]);
        CHECK(std::abs(gamma[0] - 4.99) < 1e-2);
        CHECK(std::abs(gamma[1] - 1.005) < 1e-3);
    };
};

#endif
