#include "CoolProp/Configuration.h"
#include "CoolProp/detail/json.h"
#include "src/Backends/REFPROP/REFPROPMixtureBackend.h"

namespace {

void item_to_json(const CoolProp::ConfigurationItem& item, nlohmann::json& obj) {
    const std::string name = CoolProp::config_key_to_string(item.get_key());
    switch (item.get_type()) {
        case CONFIGURATION_BOOL_TYPE:
            obj[name] = static_cast<bool>(item);
            break;
        case CONFIGURATION_INTEGER_TYPE:
            obj[name] = static_cast<int>(item);
            break;
        case CONFIGURATION_DOUBLE_TYPE:
            obj[name] = static_cast<double>(item);
            break;
        case CONFIGURATION_STRING_TYPE:
            obj[name] = static_cast<std::string>(item);
            break;
        case CONFIGURATION_ENDOFLIST_TYPE:
        case CONFIGURATION_NOT_DEFINED_TYPE:
            throw CoolProp::ValueError();
    }
}

void item_from_json(CoolProp::ConfigurationItem& item, const nlohmann::json& val) {
    switch (item.get_type()) {
        case CONFIGURATION_BOOL_TYPE:
            if (!val.is_boolean()) {
                throw CoolProp::ValueError(format("Input is not boolean"));
            }
            item.set_bool(val.get<bool>());
            break;
        case CONFIGURATION_INTEGER_TYPE:
            if (!val.is_number_integer()) {
                throw CoolProp::ValueError(format("Input is not integer"));
            }
            item.set_integer(val.get<int>());
            break;
        case CONFIGURATION_DOUBLE_TYPE:
            if (!val.is_number()) {
                throw CoolProp::ValueError(
                  format("Input [%s] is not double (or something that can be cast to double)", val.dump().c_str()));
            }
            item.set_double(val.get<double>());
            break;
        case CONFIGURATION_STRING_TYPE:
            if (!val.is_string()) {
                throw CoolProp::ValueError(format("Input is not string"));
            }
            item.set_string(val.get<std::string>());
            break;
        case CONFIGURATION_ENDOFLIST_TYPE:
        case CONFIGURATION_NOT_DEFINED_TYPE:
            throw CoolProp::ValueError();
    }
}

}  // namespace

namespace CoolProp {

std::string config_key_to_string(configuration_keys keys) {
    switch (keys) {
        /* ***MAGIC WARNING**!!
         * See http://stackoverflow.com/a/148610
         * See http://stackoverflow.com/questions/147267/easy-way-to-use-variables-of-enum-types-as-string-in-c#202511
         */
#define X(Enum, String, Default, Desc) \
    case Enum:                         \
        return String;                 \
        break;
        CONFIGURATION_KEYS_ENUM
#undef X
    }
    return "";  // will never get here, just to make compiler happy
};

std::string config_key_description(configuration_keys keys) {
    switch (keys) {
/* ***MAGIC WARNING**!!
        * See http://stackoverflow.com/a/148610
        * See http://stackoverflow.com/questions/147267/easy-way-to-use-variables-of-enum-types-as-string-in-c#202511
        */
#define X(Enum, String, Default, Desc) \
    case Enum:                         \
        return Desc;                   \
        break;
        CONFIGURATION_KEYS_ENUM
#undef X
    }
    return "";  // will never get here, just to make compiler happy
};

std::string config_key_description(const std::string& key) {
/* ***MAGIC WARNING**!!
    * See http://stackoverflow.com/a/148610
    * See http://stackoverflow.com/questions/147267/easy-way-to-use-variables-of-enum-types-as-string-in-c#202511
    */
#define X(Enum, String, Default, Desc) \
    if (key == (String)) {             \
        return (Desc);                 \
    }
    CONFIGURATION_KEYS_ENUM
#undef X
    return "INVALID KEY";
};

/// Go from string to enum key
configuration_keys config_string_to_key(const std::string& s) {
/* See http://stackoverflow.com/a/148610
     * See http://stackoverflow.com/questions/147267/easy-way-to-use-variables-of-enum-types-as-string-in-c#202511
     */
#define X(Enum, String, Default, Desc) \
    if (s == (String)) {               \
        return Enum;                   \
    }
    CONFIGURATION_KEYS_ENUM
#undef X

    // Nothing else has fired
    throw ValueError();
};

std::unique_ptr<Configuration> pconfig;
/// A helper function to ensure that configuration is not accessed before it is initialized (was formerly static)
Configuration* _get_config() {
    if (!pconfig) {
        pconfig = std::make_unique<Configuration>();
    }
    return pconfig.get();
}

void set_config_bool(configuration_keys key, bool val) {
    _get_config()->get_item(key).set_bool(val);
}
void set_config_int(configuration_keys key, int val) {
    _get_config()->get_item(key).set_integer(val);
}
void set_config_double(configuration_keys key, double val) {
    _get_config()->get_item(key).set_double(val);
}
void set_config_string(configuration_keys key, const std::string& val) {
    _get_config()->get_item(key).set_string(val);
    if (key == ALTERNATIVE_REFPROP_PATH || key == ALTERNATIVE_REFPROP_HMX_BNC_PATH || key == ALTERNATIVE_REFPROP_LIBRARY_PATH) {
        CoolProp::force_unload_REFPROP();
    }
}

bool get_config_bool(configuration_keys key) {
    return static_cast<bool>(_get_config()->get_item(key));
}
int get_config_int(configuration_keys key) {
    return static_cast<int>(_get_config()->get_item(key));
}
double get_config_double(configuration_keys key) {
    return static_cast<double>(_get_config()->get_item(key));
}
std::string get_config_string(configuration_keys key) {
    return static_cast<std::string>(_get_config()->get_item(key));
}
std::string get_config_as_json_string() {
    nlohmann::json doc = nlohmann::json::object();
    for (auto& kv : _get_config()->get_items()) {
        item_to_json(kv.second, doc);
    }
    return doc.dump();
}
void set_config_as_json_string(const std::string& s) {
    nlohmann::json doc = cpjson::parse(s);

    // First pass: validate all keys
    for (auto& [name, value] : doc.items()) {
        try {
            _get_config()->get_item(config_string_to_key(name));
        } catch (std::exception& e) {
            throw ValueError(format("Unable to parse json file with error: %s", e.what()));
        }
    }

    // Second pass: set the values
    for (auto& [name, value] : doc.items()) {
        ConfigurationItem& item = _get_config()->get_item(config_string_to_key(name));
        try {
            item_from_json(item, value);
        } catch (std::exception& e) {
            throw ValueError(format("Unable to parse json file with error: %s", e.what()));
        }
    }
}

}  // namespace CoolProp
