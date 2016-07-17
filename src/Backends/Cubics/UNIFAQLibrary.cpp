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
            // userid is an optional user identifier
            if ((*itr).HasMember("userid")){
                c.userid = (*itr)["userid"].GetString();
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
        for (std::vector<Group>::const_iterator it = groups.cbegin(); it != groups.cend(); ++it) {
            if (it->sgi == sgi) { return *it; }
        }
        throw CoolProp::ValueError("Could not find group");
    }
    bool UNIFAQParameterLibrary::has_group(int sgi) const {
        for (std::vector<Group>::const_iterator it = groups.cbegin(); it != groups.cend(); ++it) {
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
        for (std::vector<InteractionParameters>::const_iterator it = interaction_parameters.cbegin(); it != interaction_parameters.cend(); ++it) {
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
        throw CoolProp::ValueError("Could not find interaction pair");
    }

    Component UNIFAQParameterLibrary::get_component(const std::string &identifier, const std::string &value) const {
        if (identifier == "name"){
            for (std::vector<Component>::const_iterator it = components.cbegin(); it != components.cend(); ++it ){
                if (it->name == value ){ return *it; }
            }
        }
        throw CoolProp::ValueError("Could not find component");
    }

}; /* namespace UNIFAQLibrary */