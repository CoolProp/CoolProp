#ifndef RAPIDJSON_COOLPROP_H
#define RAPIDJSON_COOLPROP_H

// On PowerPC, we are going to use the stdint.h integer types and not let rapidjson use its own
#if defined(__powerpc__)
typedef unsigned int UINT32;
#include "stdint.h"
#define RAPIDJSON_NO_INT64DEFINE
#endif

#include "Exceptions.h"
#include "CoolPropTools.h"

#include "rapidjson/rapidjson.h"
#include "rapidjson/document.h"
#include "rapidjson/filestream.h"	// wrapper of C stream for prettywriter as output
#include "rapidjson/prettywriter.h"	// for stringify JSON
#include "rapidjson/stringbuffer.h" // for string buffer

#include <assert.h>

namespace cpjson
{
	/// A convenience function to get a double from a JSON value, including error checking
	inline double get_double(rapidjson::Value &v, std::string m)
	{
		if (!v.HasMember(m.c_str())){ throw CoolProp::ValueError(format("Does not have member [%s]",m.c_str())); }
		rapidjson::Value &el = v[m.c_str()];
		if (!el.IsNumber()){  throw CoolProp::ValueError(format("Member [%s] is not a number",m.c_str())); }
		else
		{
			return el.GetDouble();
		}
	};
    /// A convenience function to get a bool from a JSON value, including error checking
	inline bool get_bool(rapidjson::Value &v, std::string m)
	{
		if (!v.HasMember(m.c_str())){ throw CoolProp::ValueError(format("Does not have member [%s]",m.c_str())); }
		rapidjson::Value &el = v[m.c_str()];
        if (!el.IsBool()){  throw CoolProp::ValueError(format("Member [%s] is not a boolean",m.c_str())); }
		else
		{
            return el.GetBool();
		}
	};
     /// A convenience function to get a string from a JSON value, including error checking
	inline std::string get_string(rapidjson::Value &v, std::string m)
	{
		if (!v.HasMember(m.c_str())){ throw CoolProp::ValueError(format("Does not have member [%s]",m.c_str())); }
		rapidjson::Value &el = v[m.c_str()];
        if (!el.IsString()){  throw CoolProp::ValueError(format("Member [%s] is not a string",m.c_str())); }
		else
		{
            return el.GetString();
		}
	};

	/// A convenience function to get a double array compactly
	inline std::vector<double> get_double_array(rapidjson::Value &v)
	{
		std::vector<double> out;
		if (!v.IsArray()) { throw CoolProp::ValueError("input is not an array"); }
		for (rapidjson::Value::ValueIterator itr = v.Begin(); itr != v.End(); ++itr)
		{
            if (!itr->IsNumber()){throw CoolProp::ValueError("input is not a number");}
			out.push_back(itr->GetDouble());
		}
		return out;
	};

    /// A convenience function to get a long double array compactly
	inline std::vector<long double> get_long_double_array(rapidjson::Value &v)
	{
		std::vector<long double> out;
		if (!v.IsArray()) { throw CoolProp::ValueError("input is not an array"); }
		for (rapidjson::Value::ValueIterator itr = v.Begin(); itr != v.End(); ++itr)
		{
            if (!itr->IsNumber()){throw CoolProp::ValueError("input is not a number");}
			out.push_back(itr->GetDouble());
		}
		return out;
	};

	/// A convenience function to get a string array compactly
	inline std::vector<std::string> get_string_array(rapidjson::Value &v)
	{
		std::vector<std::string> out;
		if (!v.IsArray()) { throw CoolProp::ValueError("input is not an array"); }
		for (rapidjson::Value::ValueIterator itr = v.Begin(); itr != v.End(); ++itr)
		{
			out.push_back(itr->GetString());
		}
		return out;
	};

	/// A convenience function to get a std::string from a JSON value
	inline std::string to_string(rapidjson::Value &v)
	{
		rapidjson::StringBuffer buffer;
		rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(buffer);
		v.Accept(writer);
		return buffer.GetString();
	};

    /// A convenience function to set a double array compactly
    inline void set_double_array(const char *key, const std::vector<double> &vec, rapidjson::Value &value, rapidjson::Document &doc)
    {
        rapidjson::Value _v(rapidjson::kArrayType);
        for (unsigned int i = 0; i < vec.size(); ++i)
        {
            _v.PushBack(vec[i],doc.GetAllocator());
        }
        value.AddMember(key, _v, doc.GetAllocator());
    };

    /// A convenience function to set a double array compactly
    inline void set_long_double_array(const char *key, const std::vector<long double> &vec, rapidjson::Value &value, rapidjson::Document &doc)
    {
        rapidjson::Value _v(rapidjson::kArrayType);
        for (unsigned int i = 0; i < vec.size(); ++i)
        {
            _v.PushBack(static_cast<double>(vec[i]), doc.GetAllocator());
        }
        value.AddMember(key, _v, doc.GetAllocator());
    };

}
#endif
