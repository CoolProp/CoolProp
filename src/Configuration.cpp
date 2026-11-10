#include "CoolProp/Configuration.h"
#include "CoolProp/detail/json.h"
#include "src/Backends/REFPROP/REFPROPMixtureBackend.h"
#include <mutex>

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
                throw CoolProp::ValueError(format("Input [%s] is not double (or something that can be cast to double)", val.dump().c_str()));
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
// Thread-safe lazy initialization of the global configuration, matching the
// FluidLibrary singleton.  The previous `if (!pconfig) pconfig = ...` pattern
// raced when threads first-called get_/set_config_* concurrently: each could
// construct its own instance, the last assignment won, and settings made
// through a losing instance were silently lost.
static std::once_flag pconfig_init_flag;
/// A helper function to ensure that configuration is not accessed before it is initialized (was formerly static)
Configuration* _get_config() {
    std::call_once(pconfig_init_flag, [] { pconfig = std::make_unique<Configuration>(); });
    return pconfig.get();
}

/// Act on a configuration value that has just been changed.
/** Every setter calls this after storing the new value, including
 *  set_config_as_json_string.  Keep all such side effects here: a side effect
 *  wired into only one setter is silently skipped when the same key is set
 *  through another one (the REFPROP unload used to live in set_config_string
 *  alone, so setting ALTERNATIVE_REFPROP_PATH through JSON left the old
 *  REFPROP library loaded).
 */
static void on_config_changed(configuration_keys key) {
    if (key == ALTERNATIVE_REFPROP_PATH || key == ALTERNATIVE_REFPROP_HMX_BNC_PATH || key == ALTERNATIVE_REFPROP_LIBRARY_PATH) {
        // One lock across both steps: another thread could otherwise make a
        // load attempt between them that is cached as failed and then
        // cleared again, printing the load failure twice.
        const REFPROPLock refprop_lock(REFPROP_mutex());
        CoolProp::force_unload_REFPROP();
        // The new path deserves a new load attempt, even if an earlier
        // attempt (with the old path) failed and was cached.
        CoolProp::REFPROPMixtureBackend::forget_REFPROP_load_failure();
    }
}

void set_config_bool(configuration_keys key, bool val) {
    _get_config()->get_item(key).set_bool(val);
    on_config_changed(key);
}
void set_config_int(configuration_keys key, int val) {
    _get_config()->get_item(key).set_integer(val);
    on_config_changed(key);
}
void set_config_double(configuration_keys key, double val) {
    _get_config()->get_item(key).set_double(val);
    on_config_changed(key);
}
void set_config_string(configuration_keys key, const std::string& val) {
    _get_config()->get_item(key).set_string(val);
    on_config_changed(key);
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
        const configuration_keys key = config_string_to_key(name);
        ConfigurationItem& item = _get_config()->get_item(key);
        // Snapshot the old value: a JSON document is commonly a full dump
        // from get_config_as_json_string being restored, and firing the
        // side effects for every unchanged key would, for example, unload
        // REFPROP on every save/restore round-trip.
        nlohmann::json old_value = nlohmann::json::object();
        item_to_json(item, old_value);
        try {
            item_from_json(item, value);
        } catch (std::exception& e) {
            throw ValueError(format("Unable to parse json file with error: %s", e.what()));
        }
        nlohmann::json new_value = nlohmann::json::object();
        item_to_json(item, new_value);
        if (new_value != old_value) {
            on_config_changed(key);
        }
    }
}

}  // namespace CoolProp
