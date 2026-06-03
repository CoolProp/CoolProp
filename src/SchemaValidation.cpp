#include "CoolProp/SchemaValidation.h"

#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>

#include "CoolProp/Exceptions.h"
#include "CoolProp/detail/rapidjson.h"

namespace CoolProp {

namespace {

// Recursively copy `src` into `dst`, sorting object members by key
// (lexicographic byte-wise) at every level.  Arrays preserve order.
// `allocator` is the target document's allocator.
void copy_sorted(const rapidjson::Value& src, rapidjson::Value& dst, rapidjson::Document::AllocatorType& allocator) {
    if (src.IsObject()) {
        dst.SetObject();
        // Collect member references, sort by key, then copy in order.
        std::vector<std::pair<std::string, const rapidjson::Value*>> kv;
        kv.reserve(src.MemberCount());
        for (auto it = src.MemberBegin(); it != src.MemberEnd(); ++it) {
            kv.emplace_back(std::string(it->name.GetString(), it->name.GetStringLength()), &it->value);
        }
        std::sort(kv.begin(), kv.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
        for (const auto& [key, val_ptr] : kv) {
            rapidjson::Value name_copy;
            name_copy.SetString(key.c_str(), static_cast<rapidjson::SizeType>(key.size()), allocator);
            rapidjson::Value sub;
            copy_sorted(*val_ptr, sub, allocator);
            dst.AddMember(name_copy, sub, allocator);
        }
        return;
    }
    if (src.IsArray()) {
        dst.SetArray();
        for (auto it = src.Begin(); it != src.End(); ++it) {
            rapidjson::Value sub;
            copy_sorted(*it, sub, allocator);
            dst.PushBack(sub, allocator);
        }
        return;
    }
    // Scalar: copy via the document's allocator.
    dst.CopyFrom(src, allocator);
}

// RapidJSON validator helper — convert a SchemaPointer / DocumentPointer
// to a string for inclusion in the error message.
std::string pointer_to_string(const rapidjson::GenericPointer<rapidjson::Value>& p) {
    rapidjson::StringBuffer sb;
    p.StringifyUriFragment(sb);
    return {sb.GetString(), sb.GetSize()};
}

rapidjson::Document parse_json(const std::string& s, const char* what) {
    rapidjson::Document d;
    if (d.Parse(s.c_str(), s.size()).HasParseError()) {
        throw ValueError(std::string("parse_json: failed to parse ") + what + " (offset " + std::to_string(d.GetErrorOffset()) + ")");
    }
    return d;
}

}  // namespace

// ---------------------------------------------------------------------------
// Schema validation
// ---------------------------------------------------------------------------

void validate_against_schema(const rapidjson::Value& instance, const rapidjson::Document& schema) {
    rapidjson::SchemaDocument schema_doc(schema);
    rapidjson::SchemaValidator validator(schema_doc);
    if (instance.Accept(validator)) {
        return;
    }
    // Build a useful error: schema-pointer + instance-pointer + keyword.
    const std::string sp = pointer_to_string(validator.GetInvalidSchemaPointer());
    const std::string ip = pointer_to_string(validator.GetInvalidDocumentPointer());
    const char* kw = validator.GetInvalidSchemaKeyword();
    throw ValueError(std::string("schema validation failed: instance ") + (ip.empty() ? "/" : ip) + " violates schema " + (sp.empty() ? "/" : sp)
                     + " (keyword: " + (kw ? kw : "<unknown>") + ")");
}

void validate_against_schema(const rapidjson::Value& instance, const std::string& schema_json) {
    rapidjson::Document schema = parse_json(schema_json, "schema");
    validate_against_schema(instance, schema);
}

// ---------------------------------------------------------------------------
// Canonical JSON
// ---------------------------------------------------------------------------

std::string to_canonical_json(const rapidjson::Value& value) {
    rapidjson::Document sorted;
    rapidjson::Document::AllocatorType& alloc = sorted.GetAllocator();
    rapidjson::Value root;
    copy_sorted(value, root, alloc);
    sorted.Swap(root);

    rapidjson::StringBuffer sb;
    rapidjson::Writer<rapidjson::StringBuffer> writer(sb);
    sorted.Accept(writer);
    return {sb.GetString(), sb.GetSize()};
}

std::string to_canonical_json(const std::string& json_str) {
    rapidjson::Document d = parse_json(json_str, "input");
    return to_canonical_json(d);
}

// ---------------------------------------------------------------------------
// SWIG-friendly string overloads
// ---------------------------------------------------------------------------

void validate_json_against_schema(const std::string& instance_json, const std::string& schema_json) {
    rapidjson::Document instance = parse_json(instance_json, "instance");
    validate_against_schema(instance, schema_json);
}

std::string to_canonical_json_str(const std::string& json_str) {
    return to_canonical_json(json_str);
}

}  // namespace CoolProp
