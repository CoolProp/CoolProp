#include "CoolProp/expression/ExpressionBlock.h"

#include <map>

#include "CoolProp/Exceptions.h"
#include "CoolProp/detail/json.h"
#include "CoolProp/detail/strings.h"
#include "CoolProp/expression/ExpressionCorrelation.h"  // evaluate_at

namespace CoolProp {
namespace expression {

Program compile_block(const std::string& json_text, const std::string& context) {
    // Everything (JSON extraction + compile) inside the try so that a malformed
    // constants/arrays/formula entry is reported with the caller's context too.
    try {
        nlohmann::json j = cpjson::parse(json_text);
        std::map<std::string, double> constants;
        std::map<std::string, std::vector<double>> arrays;
        if (j.contains("constants")) {
            for (auto it = j["constants"].begin(); it != j["constants"].end(); ++it)
                constants[it.key()] = it.value().get<double>();
        }
        if (j.contains("arrays")) {
            for (auto it = j["arrays"].begin(); it != j["arrays"].end(); ++it)
                arrays[it.key()] = it.value().get<std::vector<double>>();
        }
        std::vector<std::string> state_variables;
        if (j.contains("state_variables")) {
            // Must be an array: nlohmann happily iterates a bare string as a single
            // element, so "state_variables": "T" would be silently accepted.
            if (!j["state_variables"].is_array()) throw ValueError("\"state_variables\" must be an array of names");
            for (const auto& v : j["state_variables"]) {
                if (!v.is_string()) throw ValueError("\"state_variables\" must contain only names");
                state_variables.push_back(v.get<std::string>());
            }
        }
        return compile(cpjson::get_string(j, "formula"), constants, arrays, state_variables);
    } catch (std::exception& e) {
        const std::string where = context.empty() ? std::string() : " for " + context;
        throw ValueError(format("expression block failed%s: %s", where.c_str(), e.what()));
    }
}

std::vector<std::string> ExpressionBlock::required_inputs() const {
    std::vector<std::string> out;
    const std::vector<parameters>& keys = m_program.requiredInputs();
    out.reserve(keys.size());
    for (parameters k : keys)
        out.push_back(inputName(k));
    return out;
}

double ExpressionBlock::evaluate(AbstractState& AS) const {
    // Deliberately the same call the fluid-library correlation makes, so a block
    // evaluated here and the same block loaded into a fluid cannot disagree.  The
    // caller owns the state: whatever input pair, backend, or composition `AS` was
    // set up with is what the formula sees.
    return evaluate_at(m_program, AS);
}

}  // namespace expression
}  // namespace CoolProp
