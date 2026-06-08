#ifndef COOLPROP_CONFIGURATION
#define COOLPROP_CONFIGURATION

#include "CoolProp/Exceptions.h"
#include "CoolProp/detail/tools.h"
#include <cstdlib>
#include <unordered_map>

/* See http://stackoverflow.com/a/148610
 * See http://stackoverflow.com/questions/147267/easy-way-to-use-variables-of-enum-types-as-string-in-c#202511
 * This will be used to generate an enum like:
 * enum configuration_keys {NORMALIZE_GAS_CONSTANTS, CRITICAL_SPLINES_ENABLED};
 *
 * The values in this list are given by:
 * enum, string representation of enum, default value, description
 *
 * The type of the default value specifies the only type that will be accepted for this parameter
 */
#include "CoolProp/detail/configuration_keys.h"

// Evidently SWIG+MATLAB cannot properly wrap enums within classes
enum ConfigurationDataTypes
{
    CONFIGURATION_NOT_DEFINED_TYPE = 0,
    CONFIGURATION_BOOL_TYPE,
    CONFIGURATION_DOUBLE_TYPE,
    CONFIGURATION_INTEGER_TYPE,
    CONFIGURATION_STRING_TYPE,
    CONFIGURATION_ENDOFLIST_TYPE
};

namespace CoolProp {

/// Convert the configuration key to a string in a 1-1 representation.
std::string config_key_to_string(configuration_keys keys);

/// Convert a string description to a configuration key
configuration_keys config_string_to_key(const std::string& s);

/// Return a string description of the configuration key
std::string config_key_description(configuration_keys keys);

/// Return a string description of the configuration key (with the key passed as a string)
std::string config_key_description(const std::string& key);

/// A class that contains one entry in configuration
/// Can be cast to yield the output value
class ConfigurationItem
{
   public:
    ConfigurationDataTypes get_type() const {
        return type;
    }

    /// Cast to boolean
    operator bool() const {
        check_data_type(CONFIGURATION_BOOL_TYPE);
        return v_bool;
    };
    /// Cast to double
    operator double() const {
        check_data_type(CONFIGURATION_DOUBLE_TYPE);
        return v_double;
    };
    /// Cast to string
    operator std::string() const {
        check_data_type(CONFIGURATION_STRING_TYPE);
        return v_string;
    };
    /// Cast to integer
    operator int() const {
        check_data_type(CONFIGURATION_INTEGER_TYPE);
        return v_integer;
    };
    // Initializer for bool
    ConfigurationItem(configuration_keys key, bool val) : key(key), type(CONFIGURATION_BOOL_TYPE) {

        v_bool = val;
    };
    // Initializer for integer
    ConfigurationItem(configuration_keys key, int val) : key(key), type(CONFIGURATION_INTEGER_TYPE) {

        v_integer = val;
    };
    // Initializer for double
    ConfigurationItem(configuration_keys key, double val) : key(key), type(CONFIGURATION_DOUBLE_TYPE) {

        v_double = val;
    };
    // Initializer for const char *
    ConfigurationItem(configuration_keys key, const char* val)
      : key(key), type(CONFIGURATION_STRING_TYPE), v_string(val) {

        };
    // Initializer for string
    ConfigurationItem(configuration_keys key, const std::string& val)
      : key(key), type(CONFIGURATION_STRING_TYPE), v_string(val) {

        };
    void set_bool(bool val) {
        check_data_type(CONFIGURATION_BOOL_TYPE);
        v_bool = val;
    }
    void set_integer(int val) {
        check_data_type(CONFIGURATION_INTEGER_TYPE);
        v_integer = val;
    }
    void set_double(double val) {
        check_data_type(CONFIGURATION_DOUBLE_TYPE);
        v_double = val;
    }
    void set_string(const std::string& val) {
        check_data_type(CONFIGURATION_STRING_TYPE);
        v_string = val;
    }

    configuration_keys get_key() const {
        return this->key;
    }

   private:
    void check_data_type(ConfigurationDataTypes type) const {
        if (type != this->type) {
            throw ValueError(format("type does not match"));
        }
    };
    ConfigurationDataTypes type;
    union
    {
        double v_double;
        bool v_bool;
        int v_integer;
    };
    std::string v_string;
    configuration_keys key;
};

class Configuration
{
   protected:
    std::unordered_map<configuration_keys, ConfigurationItem> items;

   public:
    Configuration() {
        set_defaults();
    };
    ~Configuration() = default;

    /// Get an item from the configuration
    ConfigurationItem& get_item(configuration_keys key) {
        // Try to find it
        auto it = items.find(key);
        // If equal to end, not found
        if (it != items.end()) {
            // Found, return it
            return it->second;
        } else {
            throw ValueError(format("invalid item"));
        }
    }
    /// Add an item to the configuration
    void add_item(const ConfigurationItem& item) {
        std::pair<configuration_keys, ConfigurationItem> pair(item.get_key(), item);
        items.insert(pair);
    };

    /// Return a reference to all of the items
    std::unordered_map<configuration_keys, ConfigurationItem>& get_items() {
        return items;
    };

    bool possibly_set_from_env(configuration_keys key) {
        /// Try to get from environment variable with the key name, prefixed by "COOLPROP_"
        std::string envkey = "COOLPROP_" + config_key_to_string(key);
        const char* envval = std::getenv(envkey.c_str());
        if (envval) {
            auto tobool = [](const std::string& x) {
                if (x == "True" || x == "true") {
                    return true;
                }
                if (x == "False" || x == "false") {
                    return false;
                }
                throw ValueError(x);
            };
            switch (get_item(key).get_type()) {
                case ConfigurationDataTypes::CONFIGURATION_STRING_TYPE:
                    items.erase(key);
                    items.emplace(key, ConfigurationItem(key, std::string(envval)));
                    break;
                case ConfigurationDataTypes::CONFIGURATION_INTEGER_TYPE: {
                    int i;
                    try {
                        i = std::stoi(envval);
                    } catch (...) {
                        auto skey = config_key_to_string(key);
                        std::string msg = "Unable to convert \"" + std::string(envval) + "\" to int for key [" + skey + "]";
                        std::cerr << msg << '\n';
                        throw ValueError(msg);
                    }
                    items.erase(key);
                    items.emplace(key, ConfigurationItem(key, i));
                    break;
                }
                case ConfigurationDataTypes::CONFIGURATION_DOUBLE_TYPE: {
                    double d;
                    try {
                        d = std::stod(envval);
                    } catch (...) {
                        auto skey = config_key_to_string(key);
                        std::string msg = "Unable to convert \"" + std::string(envval) + "\" to double for key [" + skey + "]";
                        std::cerr << msg << '\n';
                        throw ValueError(msg);
                    }
                    items.erase(key);
                    items.emplace(key, ConfigurationItem(key, d));
                    break;
                }
                case ConfigurationDataTypes::CONFIGURATION_BOOL_TYPE: {
                    bool b;
                    try {
                        b = tobool(envval);
                    } catch (...) {
                        auto skey = config_key_to_string(key);
                        std::string msg = "Unable to convert \"" + std::string(envval) + "\" to bool for key [" + skey + "]";
                        std::cerr << msg << '\n';
                        throw ValueError(msg);
                    }
                    items.erase(key);
                    items.emplace(key, ConfigurationItem(key, b));
                    break;
                }
                default: {
                    auto skey = config_key_to_string(key);
                    throw ValueError("This key [" + skey + "] has the wrong type; value was " + std::string(envval) + " ");
                }
            }
            return true;
        }
        return false;
    }

    /// Set the default values in the configuration
    void set_defaults() {
/* ***MAGIC WARNING**!!
             * See http://stackoverflow.com/a/148610
             * See http://stackoverflow.com/questions/147267/easy-way-to-use-variables-of-enum-types-as-string-in-c#202511
             */
#define X(Enum, String, Default, Desc) add_item(ConfigurationItem(Enum, Default));
        CONFIGURATION_KEYS_ENUM
#undef X

    // See if the variable is already present as environment variable
#define X(Enum, String, Default, Desc) possibly_set_from_env(Enum);
          CONFIGURATION_KEYS_ENUM
#undef X

    };
};

/// *********************************************************
///                      GETTERS
/// *********************************************************

/// Return the value of a boolean key from the configuration
bool get_config_bool(configuration_keys key);
/// Return the value of an integer key from the configuration
int get_config_int(configuration_keys key);
/// Return the value of a double configuration key
double get_config_double(configuration_keys key);
/// Return the value of a string configuration key
std::string get_config_string(configuration_keys key);
/// Get all the values in the configuration as a json-formatted string
std::string get_config_as_json_string();

/// *********************************************************
///                      SETTERS
/// *********************************************************

/// Set the value of a boolean configuration value
void set_config_bool(configuration_keys key, bool val);
/// Set the value of an integer configuration value
void set_config_int(configuration_keys key, int val);
/// Set the value of a double configuration value
void set_config_double(configuration_keys key, double val);
/// Set the value of a string configuration value
void set_config_string(configuration_keys key, const std::string& val);
/// Set the entire configuration based on a json-formatted string
void set_config_as_json_string(const std::string& s);
}  // namespace CoolProp

#endif  // COOLPROP_CONFIGURATION
