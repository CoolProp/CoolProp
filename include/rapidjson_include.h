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

#include "externals/rapidjson/include/rapidjson/rapidjson.h"
#include "externals/rapidjson/include/rapidjson/document.h"
#include "externals/rapidjson/include/rapidjson/filewritestream.h"    // wrapper of C stream for prettywriter as output
#include "externals/rapidjson/include/rapidjson/prettywriter.h"    // for stringify JSON
#include "externals/rapidjson/include/rapidjson/stringbuffer.h" // for string buffer
#include "externals/rapidjson/include/rapidjson/schema.h"

#include <cassert>

namespace cpjson
{
    struct value_information{
        bool isnull, isfalse, istrue, isbool, isobject, isarray, isnumber, isint, isint64, isuint, isuint64, isdouble, isstring;
    };
    inline value_information get_information(rapidjson::Value &v)
    {
        value_information i;
        i.isnull = v.IsNull();
        i.isfalse = v.IsFalse();
        i.istrue = v.IsTrue();
        i.isbool = v.IsBool();
        i.isobject = v.IsObject();
        i.isarray = v.IsArray();
        i.isnumber = v.IsNumber();
        i.isint = v.IsInt();
        i.isuint = v.IsUint();
        i.isint64 = v.IsInt64();
        i.isuint64 = v.IsUint64();
        i.isdouble = v.IsDouble();
        i.isstring = v.IsString();
        return i;
    };
    
    inline std::string json2string(rapidjson::Value &v)
    {
        rapidjson::StringBuffer buffer;
        rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(buffer);

        v.Accept(writer);
        return buffer.GetString();
    }
    /// A convenience function to get a double from a JSON value, including error checking
    inline int get_integer(rapidjson::Value &v, std::string m)
    {
        if (!v.HasMember(m.c_str())){ throw CoolProp::ValueError(format("Does not have member [%s]",m.c_str())); }
        rapidjson::Value &el = v[m.c_str()];
        if (!el.IsInt()){  throw CoolProp::ValueError(format("Member [%s] is not an integer",m.c_str())); }
        else
        {
            return el.GetInt();
        }
    };
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
    
    /// A convenience function to get a double array compactly
    inline std::vector<double> get_double_array(rapidjson::Value &v, std::string m)
    {
        if (!v.HasMember(m.c_str())){ throw CoolProp::ValueError(format("Does not have member [%s]",m.c_str())); }
        else{
            return get_double_array(v[m.c_str()]);
        }
    };    

    /// A convenience function to get a long double array compactly
    inline std::vector<CoolPropDbl> get_long_double_array(rapidjson::Value &v)
    {
        std::vector<CoolPropDbl> out;
        if (!v.IsArray()) { throw CoolProp::ValueError("input is not an array"); }
        for (rapidjson::Value::ValueIterator itr = v.Begin(); itr != v.End(); ++itr)
        {
            if (!itr->IsNumber()){throw CoolProp::ValueError("input is not a number");}
            out.push_back(itr->GetDouble());
        }
        return out;
    };

    /// A convenience function to get a 2D double array compactly
    inline std::vector< std::vector<double> > get_double_array2D(rapidjson::Value &v)
    {
        std::vector< std::vector<double> > out;
        std::vector<double> tmp;
        if (!v.IsArray()) {throw CoolProp::ValueError("input is not an array"); }
        for (rapidjson::Value::ValueIterator itr = v.Begin(); itr != v.End(); ++itr)
        {
            // This is here for debugging purposes
            // cpjson::value_information vi = cpjson::get_information((*itr));
            if (!(itr->IsArray())) {
                throw CoolProp::ValueError(format("input \"%s\" is not a 2D array",cpjson::json2string(v).c_str())); 
            }
            tmp.clear();
            for (rapidjson::Value::ValueIterator i = itr->Begin(); i != itr->End(); ++i)
            {
                if (!i->IsNumber()){throw CoolProp::ValueError("input is not a number");}
                tmp.push_back(i->GetDouble());
            }
            out.push_back(tmp);
        }
        return out;
    };

    /// A convenience function to get a 2D long double array compactly
    inline std::vector< std::vector<CoolPropDbl> >  get_long_double_array2D(rapidjson::Value &v)
    {
        std::vector< std::vector<CoolPropDbl> > out;
        std::vector<CoolPropDbl> tmp;
        if (!v.IsArray()) { throw CoolProp::ValueError("input is not an array"); }
        for (rapidjson::Value::ValueIterator itr = v.Begin(); itr != v.End(); ++itr)
        {
            if (!itr->IsArray()) { throw CoolProp::ValueError("input is not a 2D array"); }
            tmp.clear();
            for (rapidjson::Value::ValueIterator i = itr->Begin(); i != itr->End(); ++i)
            {
                if (!i->IsNumber()){throw CoolProp::ValueError("input is not a number");}
                tmp.push_back(i->GetDouble());
            }
            out.push_back(tmp);
        }
        return out;
    };

    /// A convenience function to get a long double array compactly
    inline std::vector<CoolPropDbl> get_long_double_array(rapidjson::Value &v, std::string name)
    {
        std::vector<CoolPropDbl> out;
        if (!v.HasMember(name.c_str())){ throw CoolProp::ValueError(format("Does not have member [%s]",name.c_str())); }
        if (!v[name.c_str()].IsArray()) { throw CoolProp::ValueError("input is not an array"); }
        for (rapidjson::Value::ValueIterator itr = v[name.c_str()].Begin(); itr != v[name.c_str()].End(); ++itr)
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
    
    /// A convenience function to get a double array compactly
    inline std::vector<std::string> get_string_array(rapidjson::Value &v, std::string m)
    {
        if (!v.HasMember(m.c_str())){ throw CoolProp::ValueError(format("Does not have member [%s]",m.c_str())); }
        else{
            return get_string_array(v[m.c_str()]);
        }
    };

    /// A convenience function to get a std::string from a JSON value
    inline std::string to_string(rapidjson::Value &v)
    {
        rapidjson::StringBuffer buffer;
        rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(buffer);
        v.Accept(writer);
        return buffer.GetString();
    };

//    /// A convenience function to set an array compactly
//    template<typename T>
//    inline void set_array(const char *key, const std::vector<T> &vec, rapidjson::Value &value, rapidjson::Document &doc)
//    {
//        rapidjson::Value _v(rapidjson::kArrayType);
//        for (unsigned int i = 0; i < vec.size(); ++i)
//        {
//            _v.PushBack(vec[i],doc.GetAllocator());
//        }
//        value.AddMember(key, _v, doc.GetAllocator());
//    };
//
//    /// A convenience function to set an array compactly
//    template<typename T>
//    inline void set_array2D(const char *key, const std::vector< std::vector<T> > &vec, rapidjson::Value &value, rapidjson::Document &doc)
//    {
//        rapidjson::Value _i(rapidjson::kArrayType);
//        rapidjson::Value _j;
//        for (unsigned int i = 0; i < vec.size(); ++i) {
//            _j = rapidjson::Value(rapidjson::kArrayType);
//            for (unsigned int j = 0; j < vec[i].size(); ++j) {
//                _j.PushBack(vec[j],doc.GetAllocator());
//            }
//            _i.PushBack(_j,doc.GetAllocator());
//        }
//        value.AddMember(key, _i, doc.GetAllocator());
//    };

    /// A convenience function to set a double array compactly
    inline void set_double_array(const char *key, const std::vector<double> &vec, rapidjson::Value &value, rapidjson::Document &doc)
    {
        rapidjson::Value _v(rapidjson::kArrayType);
        for (unsigned int i = 0; i < vec.size(); ++i)
        {
            _v.PushBack(vec[i],doc.GetAllocator());
        }
        value.AddMember(rapidjson::Value(key, doc.GetAllocator()).Move(), 
                        _v,
                        doc.GetAllocator());
    };

    /// A convenience function to set a double array compactly
    inline void set_long_double_array(const char * const key, const std::vector<CoolPropDbl> &vec, rapidjson::Value &value, rapidjson::Document &doc)
    {
        rapidjson::Value _v(rapidjson::kArrayType);
        for (unsigned int i = 0; i < vec.size(); ++i)
        {
            _v.PushBack(static_cast<double>(vec[i]), doc.GetAllocator());
        }
        value.AddMember(rapidjson::Value(key, doc.GetAllocator()).Move(), 
                        _v,
                        doc.GetAllocator());
    };
    
    enum schema_validation_code{
        SCHEMA_VALIDATION_OK = 0,
        SCHEMA_INVALID_JSON,
        INPUT_INVALID_JSON,
        SCHEMA_NOT_VALIDATED
    };
    /**
     * Validate a JSON-formatted string against a JSON-formatted schema string
     */
    inline schema_validation_code validate_schema(const std::string &schemaJson, const std::string &inputJson, std::string &errstr){
        rapidjson::Document sd;
        sd.Parse(schemaJson.c_str());
        if (sd.HasParseError()) { errstr = fmt::sprintf("Invalid schema: %s\n", schemaJson.c_str()); return SCHEMA_INVALID_JSON; }
        rapidjson::SchemaDocument schema(sd); // Compile a Document to SchemaDocument
        
        rapidjson::Document d;
        d.Parse(inputJson.c_str());
        if (d.HasParseError()) { errstr = fmt::sprintf("Invalid input json: %s\n", inputJson.c_str()); return INPUT_INVALID_JSON; }
        
        rapidjson::SchemaValidator validator(schema);
        if (!d.Accept(validator)) {
            // Input JSON is invalid according to the schema
            // Output diagnostic information
            rapidjson::StringBuffer sb;
            validator.GetInvalidSchemaPointer().StringifyUriFragment(sb);
            errstr = fmt::sprintf("Invalid schema: %s\n", sb.GetString());
            errstr += fmt::sprintf("Invalid keyword: %s\n", validator.GetInvalidSchemaKeyword());
            sb.Clear();
            return SCHEMA_NOT_VALIDATED;
        }
        return SCHEMA_VALIDATION_OK;
    }

}
#endif
