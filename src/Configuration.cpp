#include "Configuration.h"

namespace CoolProp
{
    
std::string config_key_to_string(configuration_keys keys)
{
    // Using this pragma, any missing keys will cause compilation error rather than warning
    #pragma GCC diagnostic error "-Wswitch"
    switch (keys)
    {
        case NORMALIZE_GAS_CONSTANTS:
            return "NORMALIZE_GAS_CONSTANTS";
        case CRITICAL_SPLINES_ENABLED:
            return "CRITICAL_SPLINES_ENABLED";
    }
    // Back to a warning again
    #pragma GCC diagnostic warning "-Wswitch"
    throw ValueError();
}; 

/// Go from string to enum key
configuration_keys config_string_to_key(std::string &s)
{
    if (s == "NORMALIZE_GAS_CONSTANTS"){
        return NORMALIZE_GAS_CONSTANTS;
    }
    else if (s == "CRITICAL_SPLINES_ENABLED"){
        return CRITICAL_SPLINES_ENABLED;
    }
    else{
        throw ValueError();
    }
};
    
static Configuration config;

bool get_config_bool(configuration_keys key){ 
    return static_cast<bool>(config.get_item(key)); 
}
double get_config_double(configuration_keys key){ 
    return static_cast<double>(config.get_item(key)); 
}
std::string get_config_string(configuration_keys key){ 
    return static_cast<std::string>(config.get_item(key)); 
}
rapidjson::Document get_config_as_json(){
    // Init the rapidjson doc
    rapidjson::Document doc;
    doc.SetObject();
    
    // Get the items
    std::map<configuration_keys, ConfigurationItem> items = config.get_items();
    for (std::map<configuration_keys, ConfigurationItem>::iterator it = items.begin(); it != items.end(); ++it){
        it->second.add_to_json(doc, doc);
    }
    // Convert to string
    return doc;
}
std::string get_config_as_json_string(){
    rapidjson::Document doc = get_config_as_json();
    return cpjson::to_string(doc);
}
void set_config_as_json(rapidjson::Value &val){
    
    // First check that all keys are valid
    for (rapidjson::Value::MemberIterator it = val.MemberBegin(); it != val.MemberEnd(); ++it){
        try{
            // Try to get the key for the string
            std::string s = std::string(it->name.GetString());
            configuration_keys key = config_string_to_key(s);
            // Try to retrieve the item from the config for this key
            ConfigurationItem & item = config.get_item(key);
        }
        catch(std::exception &e)
        {
            throw ValueError(format("Unable to parse json file with error: %s", e.what()));
        }
    }
    
    // Now we actually set the values
    for (rapidjson::Value::MemberIterator it = val.MemberBegin(); it != val.MemberEnd(); ++it){
        // Try to get the key for the string
        std::string s = std::string(it->name.GetString());
        configuration_keys key = config_string_to_key(s);
        // Try to retrieve the item from the config for this key
        ConfigurationItem & item = config.get_item(key);
        try{
            // Set the value from what is stored in the json value
            item.set_from_json(it->value);
        }
        catch(std::exception &e)
        {
            throw ValueError(format("Unable to parse json file with error: %s", e.what()));
        }
    }
}
void set_config_as_json_string(std::string &s){
    // Init the rapidjson doc
    rapidjson::Document doc;
    doc.Parse<0>(s.c_str());
    set_config_as_json(doc);
}

}

