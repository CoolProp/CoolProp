#include "UNIFAQLibrary.h"

namespace UNIFAQLibrary{

    void UNIFAQParameterLibrary::jsonize(std::string &s, rapidjson::Document &d)
    {
        d.Parse<0>(s.c_str());
        if (d.HasParseError()) {
            throw -1;
        }
        else {
            return;
        }
    }
    void UNIFAQParameterLibrary::populate(rapidjson::Value &group_data, rapidjson::Value &interaction_data, rapidjson::Value &comp_data)
    {
        // Schema should have been used to validate the data already, so by this point we are can safely consume the data without checking ...
        for (rapidjson::Value::ValueIterator itr = group_data.Begin(); itr != group_data.End(); ++itr)
        {
            Group g;
            g.sgi = (*itr)["sgi"].GetInt();
            g.mgi = (*itr)["mgi"].GetInt();
            g.R_k = (*itr)["R_k"].GetDouble();
            g.Q_k = (*itr)["Q_k"].GetDouble();
            groups.push_back(g);
        }
        for (rapidjson::Value::ValueIterator itr = interaction_data.Begin(); itr != interaction_data.End(); ++itr)
        {
            InteractionParameters ip;
            ip.mgi1 = (*itr)["mgi1"].GetInt();
            ip.mgi2 = (*itr)["mgi2"].GetInt();
            ip.a_ij = (*itr)["a_ij"].GetDouble();
            ip.a_ji = (*itr)["a_ji"].GetDouble();
            ip.b_ij = (*itr)["b_ij"].GetDouble();
            ip.b_ji = (*itr)["b_ji"].GetDouble();
            ip.c_ij = (*itr)["c_ij"].GetDouble();
            ip.c_ji = (*itr)["c_ji"].GetDouble();
            interaction_parameters.push_back(ip);
        }
        for (rapidjson::Value::ValueIterator itr = comp_data.Begin(); itr != comp_data.End(); ++itr)
        {
            Component c;
            c.inchikey = (*itr)["inchikey"].GetString();
            c.registry_number = (*itr)["registry_number"].GetString();
            c.name = (*itr)["name"].GetString();
            c.Tc = (*itr)["Tc"].GetDouble();
            c.pc = (*itr)["pc"].GetDouble();
            c.acentric = (*itr)["acentric"].GetDouble();
            c.molemass = (*itr)["molemass"].GetDouble();
            // userid is an optional user identifier
            if ((*itr).HasMember("userid")){
                c.userid = (*itr)["userid"].GetString();
            }
            // If provided, store information about the alpha function in use
            if ((*itr).HasMember("alpha") && (*itr)["alpha"].IsObject()){
                rapidjson::Value &alpha = (*itr)["alpha"];
                c.alpha_type = cpjson::get_string(alpha, "type");
                c.alpha_coeffs = cpjson::get_double_array(alpha, "c");
            }
            else{
                c.alpha_type = "default";
            }
            
            rapidjson::Value &groups = (*itr)["groups"];
            for (rapidjson::Value::ValueIterator itrg = groups.Begin(); itrg != groups.End(); ++itrg)
            {
                int count = (*itrg)["count"].GetInt();
                int sgi = (*itrg)["sgi"].GetInt();
                if  (has_group(sgi)){
                    ComponentGroup cg(count, get_group(sgi));
                    c.groups.push_back(cg);
                }
            }
            components.push_back(c);
        }
    }
    void UNIFAQParameterLibrary::populate(std::string &group_data, std::string &interaction_data, std::string &decomp_data)
    {
        rapidjson::Document group_JSON; jsonize(group_data, group_JSON);
        rapidjson::Document interaction_JSON; jsonize(interaction_data, interaction_JSON);
        rapidjson::Document decomp_JSON; jsonize(decomp_data, decomp_JSON);
        populate(group_JSON, interaction_JSON, decomp_JSON);
    }
    Group UNIFAQParameterLibrary::get_group(int sgi) const {
        for (std::vector<Group>::const_iterator it = groups.begin(); it != groups.end(); ++it) {
            if (it->sgi == sgi) { return *it; }
        }
        throw CoolProp::ValueError("Could not find group");
    }
    bool UNIFAQParameterLibrary::has_group(int sgi) const {
        for (std::vector<Group>::const_iterator it = groups.begin(); it != groups.end(); ++it) {
            if (it->sgi == sgi) { return true; }
        }
        return false;
    }

    InteractionParameters UNIFAQParameterLibrary::get_interaction_parameters(int mgi1, int mgi2) const {

        // If both mgi are the same, yield all zeros for the interaction parameters
        if (mgi1 == mgi2){
            InteractionParameters ip; ip.mgi1 = mgi1; ip.mgi2 = mgi2;
            ip.zero_out();
            return ip;
        }
        for (std::vector<InteractionParameters>::const_iterator it = interaction_parameters.begin(); it != interaction_parameters.end(); ++it) {
            if (it->mgi1 == mgi1 && it->mgi2 == mgi2) {
                // Correct order, return it
                return *it;
            }
            if (it->mgi2 == mgi1 && it->mgi1 == mgi2) {
                // Backwards, swap the parameters
                InteractionParameters ip = *it;
                ip.swap();
                return ip;
            }
        }
        throw CoolProp::ValueError(format("Could not find interaction between pair mgi[%d]-mgi[%d]", static_cast<int>(mgi1), static_cast<int>(mgi2)));
    }

    Component UNIFAQParameterLibrary::get_component(const std::string &identifier, const std::string &value) const {
        if (identifier == "name"){
            for (std::vector<Component>::const_iterator it = components.begin(); it != components.end(); ++it ){
                if (it->name == value ){ return *it; }
            }
        }
        throw CoolProp::ValueError(format("Could not find component: %s with identifier: %s", value.c_str(), identifier.c_str()));
    }

}; /* namespace UNIFAQLibrary */

#if defined(ENABLE_CATCH)
#include "catch.hpp"

#include "UNIFAQ.h"

TEST_CASE("Check Poling example for UNIFAQ", "[UNIFAQ]")
{
    std::string acetone_pentane_groups = "[{ \"Tc\": 508.1, \"acentric\": 0.3071, \"groups\": [ { \"count\": 1,  \"sgi\": 1 },  {\"count\": 1, \"sgi\": 18 } ],  \"inchikey\": \"?????????????\",  \"name\": \"Acetone\", \"pc\": 4700000.0, \"registry_number\": \"67-64-1\", \"userid\": \"\" },  { \"Tc\": 469.7000000000001,  \"acentric\": 0.251,    \"groups\": [ { \"count\": 2, \"sgi\": 1 }, { \"count\": 3, \"sgi\": 2 } ],  \"inchikey\": \"?????????????\", \"name\": \"n-Pentane\", \"pc\": 3370000.0,  \"registry_number\": \"109-66-0\",  \"userid\": \"\" } ]";
    std::string groups = "[{\"Q_k\": 0.848, \"R_k\": 0.9011, \"maingroup_name\": \"CH2\", \"mgi\": 1, \"sgi\": 1, \"subgroup_name\": \"CH3\"},"
        "{\"Q_k\": 0.540, \"R_k\": 0.6744, \"maingroup_name\": \"CH2\", \"mgi\": 1, \"sgi\": 2, \"subgroup_name\": \"CH2\"},"
        "{\"Q_k\": 1.488, \"R_k\": 1.6724, \"maingroup_name\": \"CH2CO\", \"mgi\": 9, \"sgi\": 18, \"subgroup_name\": \"CH3CO\"}]";
    std::string interactions = "[{\"a_ij\": 476.4, \"a_ji\": 26.76, \"b_ij\": 0.0, \"b_ji\": 0.0,  \"c_ij\": 0.0, \"c_ji\": 0.0, \"mgi1\": 1, \"mgi2\": 9}]";

    SECTION("Validate AC for acetone + n-pentane") {
        UNIFAQLibrary::UNIFAQParameterLibrary lib;
        CHECK_NOTHROW(lib.populate(groups, interactions, acetone_pentane_groups););
        UNIFAQ::UNIFAQMixture mix(lib);
        std::vector<std::string> names; names.push_back("Acetone"); names.push_back("n-Pentane");
        mix.set_components("name",names);
        mix.set_interaction_parameters();
       
        std::vector<double> z(2,0.047); z[1] = 1-z[0];
        mix.set_mole_fractions(z);
        CHECK_NOTHROW(mix.set_temperature(307, z););
        
        double lngammaR0 = mix.ln_gamma_R(0);
        double lngammaR1 = mix.ln_gamma_R(1);
        CAPTURE(lngammaR0);
        CAPTURE(lngammaR1);
        CHECK(std::abs(lngammaR0 - 1.66) < 1e-2);
        CHECK(std::abs(lngammaR1 - 5.68e-3) < 1e-3);

        double gamma0 = mix.activity_coefficient(0);
        double gamma1 = mix.activity_coefficient(1);
        CAPTURE(gamma0);
        CAPTURE(gamma1);
        CHECK(std::abs(gamma0 - 4.99) < 1e-2);
        CHECK(std::abs(gamma1 - 1.005) < 1e-3);
    };
};

#endif
