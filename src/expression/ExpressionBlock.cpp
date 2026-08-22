#include "CoolProp/expression/ExpressionBlock.h"

#include <map>
#include <memory>

#include "CoolProp/AbstractState.h"
#include "CoolProp/CoolProp.h"  // extract_backend
#include "CoolProp/DataStructures.h"
#include "CoolProp/Exceptions.h"
#include "CoolProp/detail/json.h"
#include "CoolProp/detail/strings.h"

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
        return compile(cpjson::get_string(j, "formula"), constants, arrays);
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

double ExpressionBlock::evaluate(double T, double rhomolar, const std::string& fluid) const {
    const std::vector<parameters>& keys = m_program.requiredInputs();
    std::vector<double> vals(keys.size(), 0.0);
    if (fluid.empty()) {
        // No fluid, no EOS: only the two independent variables the caller passed
        // are knowable.  Anything else is a hard error rather than a silent zero.
        for (std::size_t i = 0; i < keys.size(); ++i) {
            if (keys[i] == iT) {
                vals[i] = T;
            } else if (keys[i] == iDmolar) {
                vals[i] = rhomolar;
            } else {
                throw ValueError(format("expression input '%s' needs a fluid; pass one to evaluate()", inputName(keys[i]).c_str()));
            }
        }
    } else {
        std::string backend, fluids;
        extract_backend(fluid, backend, fluids);
        std::shared_ptr<AbstractState> AS(AbstractState::factory(backend, fluids));
        AS->update(DmolarT_INPUTS, rhomolar, T);
        for (std::size_t i = 0; i < keys.size(); ++i)
            vals[i] = AS->keyed_output(keys[i]);
    }
    return m_program.evaluate(vals.empty() ? nullptr : vals.data());
}

}  // namespace expression
}  // namespace CoolProp
