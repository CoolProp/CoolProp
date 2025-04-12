#ifndef COOLPROP_CONFIGURATION
#define COOLPROP_CONFIGURATION

#include "Exceptions.h"
#include "CoolPropTools.h"
#include <cstdlib>

#if !defined(SWIG)  // Hide this for swig - Swig gets confused
#    include "rapidjson_include.h"
#endif

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
#define CONFIGURATION_KEYS_ENUM                                                                                                                      \
    X(NORMALIZE_GAS_CONSTANTS, "NORMALIZE_GAS_CONSTANTS", true, "If true, for mixtures, the molar gas constant (R) will be set to the CODATA value") \
    X(CRITICAL_WITHIN_1UK, "CRITICAL_WITHIN_1UK", true,                                                                                              \
      "If true, any temperature within 1 uK of the critical temperature will be considered to be AT the critical point")                             \
    X(CRITICAL_SPLINES_ENABLED, "CRITICAL_SPLINES_ENABLED", true,                                                                                    \
      "If true, the critical splines will be used in the near-vicinity of the critical point")                                                       \
    X(SAVE_RAW_TABLES, "SAVE_RAW_TABLES", false, "If true, the raw, uncompressed tables will also be written to file")                               \
    X(ALTERNATIVE_TABLES_DIRECTORY, "ALTERNATIVE_TABLES_DIRECTORY", "",                                                                              \
      "If provided, this path will be the root directory for the tabular data.  Otherwise, ${HOME}/.CoolProp/Tables is used")                        \
    X(ALTERNATIVE_REFPROP_PATH, "ALTERNATIVE_REFPROP_PATH", "",                                                                                      \
      "An alternative path to be provided to the directory that contains REFPROP's fluids and mixtures directories.  If provided, the SETPATH "      \
      "function will be called with this directory prior to calling any REFPROP functions.")                                                         \
    X(ALTERNATIVE_REFPROP_HMX_BNC_PATH, "ALTERNATIVE_REFPROP_HMX_BNC_PATH", "",                                                                      \
      "An alternative path to the HMX.BNC file.  If provided, it will be passed into REFPROP's SETUP or SETMIX routines")                            \
    X(ALTERNATIVE_REFPROP_LIBRARY_PATH, "ALTERNATIVE_REFPROP_LIBRARY_PATH", "",                                                                      \
      "An alternative path to the shared library file.  If provided, it will be used to load REFPROP")                                               \
    X(REFPROP_DONT_ESTIMATE_INTERACTION_PARAMETERS, "REFPROP_DONT_ESTIMATE_INTERACTION_PARAMETERS", false,                                           \
      "If true, if the binary interaction parameters in REFPROP are estimated, throw an error rather than silently continuing")                      \
    X(REFPROP_IGNORE_ERROR_ESTIMATED_INTERACTION_PARAMETERS, "REFPROP_IGNORE_ERROR_ESTIMATED_INTERACTION_PARAMETERS", false,                         \
      "If true, if the binary interaction parameters in REFPROP are unable to be estimated, silently continue rather than failing")                  \
    X(REFPROP_USE_GERG, "REFPROP_USE_GERG", false,                                                                                                   \
      "If true, rather than using the highly-accurate pure fluid equations of state, use the pure-fluid EOS from GERG-2008")                         \
    X(REFPROP_ERROR_THRESHOLD, "REFPROP_ERROR_THRESHOLD", static_cast<int>(0), "The highest acceptable error code without throwing an exception")    \
    X(REFPROP_USE_PENGROBINSON, "REFPROP_USE_PENGROBINSON", false,                                                                                   \
      "If true, rather than using the highly-accurate pure fluid equations of state, use the Peng-Robinson EOS")                                     \
    X(MAXIMUM_TABLE_DIRECTORY_SIZE_IN_GB, "MAXIMUM_TABLE_DIRECTORY_SIZE_IN_GB", 1.0,                                                                 \
      "The maximum allowed size of the directory that is used to store tabular data")                                                                \
    X(DONT_CHECK_PROPERTY_LIMITS, "DONT_CHECK_PROPERTY_LIMITS", false,                                                                               \
      "If true, when possible, CoolProp will skip checking whether values are inside the property limits")                                           \
    X(HENRYS_LAW_TO_GENERATE_VLE_GUESSES, "HENRYS_LAW_TO_GENERATE_VLE_GUESSES", false,                                                               \
      "If true, when doing water-based mixture dewpoint calculations, use Henry's Law to generate guesses for liquid-phase composition")             \
    X(PHASE_ENVELOPE_STARTING_PRESSURE_PA, "PHASE_ENVELOPE_STARTING_PRESSURE_PA", 100.0, "Starting pressure [Pa] for phase envelope construction")   \
    X(R_U_CODATA, "R_U_CODATA", 8.3144598,                                                                                                           \
      "The value for the ideal gas constant in J/mol/K according to CODATA 2014.  This value is used to harmonize all the ideal gas constants. "     \
      "This is especially important in the critical region.")                                                                                        \
    X(VTPR_UNIFAC_PATH, "VTPR_UNIFAC_PATH", "", "The path to the directory containing the UNIFAC JSON files.  Should be slash terminated")           \
    X(SPINODAL_MINIMUM_DELTA, "SPINODAL_MINIMUM_DELTA", 0.5,                                                                                         \
      "The minimal delta to be used in tracing out the spinodal; make sure that the EOS has a spinodal at this value of delta=rho/rho_r")            \
    X(OVERWRITE_FLUIDS, "OVERWRITE_FLUIDS", false,                                                                                                   \
      "If true, and a fluid is added to the fluids library that is already there, rather than not adding the fluid (and probably throwing an "       \
      "exception), overwrite it")                                                                                                                    \
    X(OVERWRITE_DEPARTURE_FUNCTION, "OVERWRITE_DEPARTURE_FUNCTION", false,                                                                           \
      "If true, and a departure function to be added is already there, rather than not adding the departure function (and probably throwing an "     \
      "exception), overwrite it")                                                                                                                    \
    X(OVERWRITE_BINARY_INTERACTION, "OVERWRITE_BINARY_INTERACTION", false,                                                                           \
      "If true, and a pair of binary interaction pairs to be added is already there, rather than not adding the binary interaction pair (and "       \
      "probably throwing an exception), overwrite it")                                                                                               \
    X(USE_GUESSES_IN_PROPSSI, "USE_GUESSES_IN_PROPSSI", false,                                                                                       \
      "If true, calls to the vectorized versions of PropsSI use the previous state as guess value while looping over the input vectors, only makes " \
      "sense when working with a single fluid and with points that are not too far from each other.")                                                \
    X(ASSUME_CRITICAL_POINT_STABLE, "ASSUME_CRITICAL_POINT_STABLE", false,                                                                           \
      "If true, evaluation of the stability of critical point will be skipped and point will be assumed to be stable")                               \
    X(VTPR_ALWAYS_RELOAD_LIBRARY, "VTPR_ALWAYS_RELOAD_LIBRARY", false,                                                                               \
      "If true, the library will always be reloaded, no matter what is currently loaded")                                                            \
    X(FLOAT_PUNCTUATION, "FLOAT_PUNCTUATION", ".", "The first character of this string will be used as the separator between the number fraction.")  \
    X(LIST_STRING_DELIMITER, "LIST_STRING_DELIMITER", ",", "The delimiter to be used when converting a list of strings to a string")

// Use preprocessor to create the Enum
enum configuration_keys
{
#define X(Enum, String, Default, Desc) Enum,
    CONFIGURATION_KEYS_ENUM
#undef X
};

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
    ConfigurationDataTypes get_type() const { return type; }
    
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
    ConfigurationItem(configuration_keys key, bool val) {
        this->key = key;
        type = CONFIGURATION_BOOL_TYPE;
        v_bool = val;
    };
    // Initializer for integer
    ConfigurationItem(configuration_keys key, int val) {
        this->key = key;
        type = CONFIGURATION_INTEGER_TYPE;
        v_integer = val;
    };
    // Initializer for double
    ConfigurationItem(configuration_keys key, double val) {
        this->key = key;
        type = CONFIGURATION_DOUBLE_TYPE;
        v_double = val;
    };
    // Initializer for const char *
    ConfigurationItem(configuration_keys key, const char* val) {
        this->key = key;
        type = CONFIGURATION_STRING_TYPE;
        v_string = val;
    };
    // Initializer for string
    ConfigurationItem(configuration_keys key, const std::string& val) {
        this->key = key;
        type = CONFIGURATION_STRING_TYPE;
        v_string = val;
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

    configuration_keys get_key(void) const {
        return this->key;
    }
#if !defined(SWIG)
    /// Cast to rapidjson::Value
    void add_to_json(rapidjson::Value& val, rapidjson::Document& d) const {
        std::string name_string = config_key_to_string(key);
        rapidjson::Value name(name_string.c_str(), d.GetAllocator());
        switch (type) {
            case CONFIGURATION_BOOL_TYPE: {
                rapidjson::Value v(v_bool);
                val.AddMember(name, v, d.GetAllocator());
                break;
            }
            case CONFIGURATION_INTEGER_TYPE: {
                rapidjson::Value v(v_integer);
                val.AddMember(name, v, d.GetAllocator());
                break;
            }
            case CONFIGURATION_DOUBLE_TYPE: {
                rapidjson::Value v(v_double);  // Try to upcast
                val.AddMember(name, v, d.GetAllocator());
                break;
            }
            case CONFIGURATION_STRING_TYPE: {
                rapidjson::Value v(v_string.c_str(), d.GetAllocator());
                val.AddMember(name, v, d.GetAllocator());
                break;
            }
            case CONFIGURATION_ENDOFLIST_TYPE:
            case CONFIGURATION_NOT_DEFINED_TYPE:
                throw ValueError();
        }
    }
    void set_from_json(rapidjson::Value& val) {
        switch (type) {
            case CONFIGURATION_BOOL_TYPE:
                if (!val.IsBool()) {
                    throw ValueError(format("Input is not boolean"));
                };
                v_bool = val.GetBool();
                break;
            case CONFIGURATION_INTEGER_TYPE:
                if (!val.IsInt()) {
                    throw ValueError(format("Input is not integer"));
                };
                v_integer = val.GetInt();
                break;
            case CONFIGURATION_DOUBLE_TYPE: {
                if (!val.IsDouble() && !val.IsInt()) {
                    throw ValueError(format("Input [%s] is not double (or something that can be cast to double)", cpjson::to_string(val).c_str()));
                };
                if (val.IsDouble()) {
                    v_double = val.GetDouble();
                } else {
                    v_double = static_cast<double>(val.GetInt());
                }
                break;
            }
            case CONFIGURATION_STRING_TYPE:
                if (!val.IsString()) {
                    throw ValueError(format("Input is not string"));
                };
                v_string = val.GetString();
                break;
            case CONFIGURATION_ENDOFLIST_TYPE:
            case CONFIGURATION_NOT_DEFINED_TYPE:
                throw ValueError();
        }
    }
#endif  // !defined(SWIG)

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
    std::map<configuration_keys, ConfigurationItem> items;

   public:
    Configuration() {
        set_defaults();
    };
    ~Configuration(){};

    /// Get an item from the configuration
    ConfigurationItem& get_item(configuration_keys key) {
        // Try to find it
        std::map<configuration_keys, ConfigurationItem>::iterator it = items.find(key);
        // If equal to end, not found
        if (it != items.end()) {
            // Found, return it
            return it->second;
        } else {
            throw ValueError(format("invalid item"));
        }
    }
    /// Add an item to the configuration
    void add_item(ConfigurationItem item) {
        std::pair<configuration_keys, ConfigurationItem> pair(item.get_key(), item);
        items.insert(pair);
    };

    /// Return a reference to all of the items
    std::map<configuration_keys, ConfigurationItem>& get_items(void) {
        return items;
    };
    
    bool possibly_set_from_env(configuration_keys key){
        /// Try to get from environment variable with the key name, prefixed by "COOLPROP_"
        std::string envkey = "COOLPROP_" + config_key_to_string(key);
        const char *envval = std::getenv(envkey.c_str());
        if (envval){
            auto tobool = [](const std::string x){
                if (x == "True" || x == "true"){ return true;}
                if (x == "False" || x == "false"){ return false;}
                throw ValueError(x);
            };
            switch (get_item(key).get_type()){
                case ConfigurationDataTypes::CONFIGURATION_STRING_TYPE:
                    items.erase(key); items.emplace(key, ConfigurationItem(key, std::string(envval)));
                    break;
                case ConfigurationDataTypes::CONFIGURATION_INTEGER_TYPE:
                    int i;
                    try{
                        i = std::stoi(envval);
                    }
                    catch(...){
                        auto skey = config_key_to_string(key);
                        std::string msg = "Unable to convert \""+std::string(envval)+"\" to int for key ["+skey+"]";
                        std::cerr << msg << std::endl;
                        throw ValueError(msg);
                    }
                    items.erase(key); items.emplace(key, ConfigurationItem(key, i));
                    break;
                case ConfigurationDataTypes::CONFIGURATION_DOUBLE_TYPE:
                    int d;
                    try{
                        d = std::stod(envval);
                    }
                    catch(...){
                        auto skey = config_key_to_string(key);
                        std::string msg = "Unable to convert \""+std::string(envval)+"\" to double for key ["+skey+"]";
                        std::cerr << msg << std::endl;
                        throw ValueError(msg);
                    }
                    items.erase(key); items.emplace(key, ConfigurationItem(key, d));
                    break;
                case ConfigurationDataTypes::CONFIGURATION_BOOL_TYPE:
                    int b;
                    try{
                        b = tobool(envval);
                    }
                    catch(...){
                        auto skey = config_key_to_string(key);
                        std::string msg = "Unable to convert \""+std::string(envval)+"\" to bool for key ["+skey+"]";
                        std::cerr << msg << std::endl;
                        throw ValueError(msg);
                    }
                    items.erase(key); items.emplace(key, ConfigurationItem(key, b));
                    break;
                default:
                    auto skey = config_key_to_string(key);
                    throw ValueError("This key ["+skey+"] has the wrong type; value was "+std::string(envval)+" ");
            }
            return true;
        }
        return false;
    }

    /// Set the default values in the configuration
    void set_defaults(void) {
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
#if !defined(SWIG)  // Hide this for swig - Swig gets confused
void get_config_as_json(rapidjson::Document& doc);
#endif
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
/// Set values in the configuration based on a json file
#if !defined(SWIG)  // Hide this for swig - Swig gets confused
void set_config_json(rapidjson::Document& doc);
#endif
/// Set the entire configuration based on a json-formatted string
void set_config_as_json_string(const std::string& s);
}  // namespace CoolProp

#endif  // COOLPROP_CONFIGURATION
