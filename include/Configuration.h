#ifndef COOLPROP_CONFIGURATION
#define COOLPROP_CONFIGURATION

#include "Exceptions.h"
#include "CoolPropTools.h"

#if !defined(SWIG) // Hide this for swig - Swig gets confused
#include "rapidjson/rapidjson_include.h"
#endif

/* See http://stackoverflow.com/a/148610
 * See http://stackoverflow.com/questions/147267/easy-way-to-use-variables-of-enum-types-as-string-in-c#202511
 * This will be used to generate an enum like:
 * enum configuration_keys {NORMALIZE_GAS_CONSTANTS, CRITICAL_SPLINES_ENABLED};
 * 
 * The values in this list are given by:
 * enum, string representation of enum, default value
 * 
 * The type of the default value specifies the only type that will be accepted for this parameter
 */
#define CONFIGURATION_KEYS_ENUM \
    X(NORMALIZE_GAS_CONSTANTS, "NORMALIZE_GAS_CONSTANTS", true) \
    X(CRITICAL_WITHIN_1UK, "CRITICAL_WITHIN_1UK", true) \
    X(CRITICAL_SPLINES_ENABLED, "CRITICAL_SPLINES_ENABLED", true) \
	X(ALTERNATIVE_REFPROP_PATH, "ALTERNATIVE_REFPROP_PATH", "") \
    X(ALTERNATIVE_REFPROP_HMX_BNC_PATH, "ALTERNATIVE_REFPROP_HMX_BNC_PATH", "") \
    X(MAXIMUM_TABLE_DIRECTORY_SIZE_IN_GB, "MAXIMUM_TABLE_DIRECTORY_SIZE_IN_GB", 1.0) \

 // Use preprocessor to create the Enum
 enum configuration_keys{
  #define X(Enum, String, Default)       Enum,
   CONFIGURATION_KEYS_ENUM
  #undef X
 };

namespace CoolProp
{

/// Convert the configuration key to a string in a 1-1 representation.
std::string config_key_to_string(configuration_keys keys);
    
/// A class that contains one entry in configuration
/// Can be cast to yield the output value
class ConfigurationItem
{
    public:
        enum ConfigurationDataTypes {CONFIGURATION_NOT_DEFINED_TYPE = 0, 
                                     CONFIGURATION_BOOL_TYPE, 
                                     CONFIGURATION_DOUBLE_TYPE, 
                                     CONFIGURATION_INTEGER_TYPE, 
                                     CONFIGURATION_STRING_TYPE,
                                     CONFIGURATION_ENDOFLIST_TYPE};
        ConfigurationDataTypes type;
        
        /// Cast to boolean
        operator bool() const { check_data_type(CONFIGURATION_BOOL_TYPE);  return v_bool; };
        /// Cast to double
        operator double() const { check_data_type(CONFIGURATION_DOUBLE_TYPE);  return v_double; };
        /// Cast to string
        operator std::string() const { check_data_type(CONFIGURATION_STRING_TYPE);  return v_string; };
        // Initializer for bool
        ConfigurationItem(configuration_keys key, bool val){
            this->key = key; type = CONFIGURATION_BOOL_TYPE; v_bool = val;
        };
        // Initializer for integer
        ConfigurationItem(configuration_keys key, int val){
            this->key = key; type = CONFIGURATION_INTEGER_TYPE; v_integer = val;
        };
        // Initializer for double
        ConfigurationItem(configuration_keys key, double val){ 
            this->key = key; type = CONFIGURATION_DOUBLE_TYPE; v_double = val;
        };
		// Initializer for const char *
        ConfigurationItem(configuration_keys key, const char *val){
            this->key = key; type = CONFIGURATION_STRING_TYPE; v_string = val;
        };
        // Initializer for string
        ConfigurationItem(configuration_keys key, const std::string &val){
            this->key = key; type = CONFIGURATION_STRING_TYPE; v_string = val;
        };
		void set_bool(bool val){
			check_data_type(CONFIGURATION_BOOL_TYPE);
			v_bool = val;
		}
		void set_integer(int val){
			check_data_type(CONFIGURATION_INTEGER_TYPE);
			v_integer = val;
		}
		void set_double(double val){
			check_data_type(CONFIGURATION_DOUBLE_TYPE);
			v_double = val;
		}
		void set_string(const std::string &val){
			check_data_type(CONFIGURATION_STRING_TYPE);
			v_string = val;
		}
		
        configuration_keys get_key(void) const {
            return this->key;
        }
		#if !defined(SWIG)
        /// Cast to rapidjson::Value
        void add_to_json(rapidjson::Value &val, rapidjson::Document &d) const {
            std::string name_string = config_key_to_string(key);
            rapidjson::Value name(name_string.c_str(), d.GetAllocator());
            switch (type){
                case CONFIGURATION_BOOL_TYPE:
                {
                    rapidjson::Value v(v_bool);
                    val.AddMember(name, v, d.GetAllocator()); break;
                }
                case CONFIGURATION_INTEGER_TYPE:
                {
                    rapidjson::Value v(v_integer);
                    val.AddMember(name, v, d.GetAllocator()); break;
                }
                case CONFIGURATION_DOUBLE_TYPE:
                {
                    rapidjson::Value v(v_double);
                    val.AddMember(name, v, d.GetAllocator()); break;
                }
                case CONFIGURATION_STRING_TYPE:
                {
                    rapidjson::Value v(v_string.c_str(), d.GetAllocator());
                    val.AddMember(name, v, d.GetAllocator()); break;
                }
                case CONFIGURATION_ENDOFLIST_TYPE:
                case CONFIGURATION_NOT_DEFINED_TYPE:
                    throw ValueError();
            }
        }
        void set_from_json(rapidjson::Value &val){
            switch (type){
                case CONFIGURATION_BOOL_TYPE: if (!val.IsBool()){throw ValueError(format("Input is not boolean"));}; v_bool = val.GetBool(); break;
                case CONFIGURATION_INTEGER_TYPE: if (!val.IsInt()){throw ValueError(format("Input is not integer"));}; v_integer = val.GetInt(); break;
                case CONFIGURATION_DOUBLE_TYPE: if (!val.IsDouble()){throw ValueError(format("Input is not double"));}; v_double = val.GetDouble(); break;
                case CONFIGURATION_STRING_TYPE: if (!val.IsString()){throw ValueError(format("Input is not string"));}; v_string = val.GetString(); break; 
                case CONFIGURATION_ENDOFLIST_TYPE:
                case CONFIGURATION_NOT_DEFINED_TYPE:
                    throw ValueError();
            }
        }
		#endif // !defined(SWIG)
         
    protected:
        void check_data_type(ConfigurationDataTypes type) const {
            if (type != this->type){
                throw ValueError(format("type does not match"));
            }
        };
        double v_double;
        bool v_bool;
        int v_integer;
        std::string v_string;
        configuration_keys key;
};

class Configuration
{
    protected:
        std::map<configuration_keys,ConfigurationItem> items;
        
    public:
        Configuration(){set_defaults();};
        ~Configuration(){};
        
        /// Get an item from the configuration
        ConfigurationItem &get_item(configuration_keys key){
            // Try to find it
            std::map<configuration_keys,ConfigurationItem>::iterator it = items.find(key);
            // If equal to end, not found
            if (it != items.end()){
                // Found, return it
                return it->second;
            }
            else{
                throw ValueError(format("invalid item"));
            }
        }
        /// Add an item to the configuration
        void add_item(ConfigurationItem item)
        {
            std::pair<configuration_keys, ConfigurationItem> pair(item.get_key(), item);
            items.insert(pair);
        };
        
        /// Return a reference to all of the items
        std::map<configuration_keys, ConfigurationItem> & get_items(void){return items;};
        
        /// Set the default values in the configuration
        void set_defaults(void)
        {
            /* ***MAGIC WARNING**!!
             * See http://stackoverflow.com/a/148610
             * See http://stackoverflow.com/questions/147267/easy-way-to-use-variables-of-enum-types-as-string-in-c#202511
             */
            #define X(Enum, String, Default) \
                add_item(ConfigurationItem(Enum, Default));
                CONFIGURATION_KEYS_ENUM
            #undef X
        };
};

/// *********************************************************
///                      GETTERS
/// *********************************************************

/// Return the value of a boolean key from the configuration
bool get_config_bool(configuration_keys key);
double get_config_double(configuration_keys key);
std::string get_config_string(configuration_keys key);
#if !defined(SWIG) // Hide this for swig - Swig gets confused
void get_config_as_json(rapidjson::Document &doc);
#endif
/// Get values in the configuration based as json data in string format
std::string get_config_as_json_string();

/// *********************************************************
///                      SETTERS
/// *********************************************************

void set_config_bool(configuration_keys key, bool val);
void set_config_double(configuration_keys key, double val);
void set_config_string(configuration_keys key, const std::string &val);
/// Set values in the configuration based on a json file
#if !defined(SWIG) // Hide this for swig - Swig gets confused
void set_config_json(rapidjson::Document &doc);
#endif
void set_config_as_json_string(const std::string &s);
}

#endif // COOLPROP_CONFIGURATION
