#ifndef COOLPROP_EXPRESSION_H
#define COOLPROP_EXPRESSION_H

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "CoolProp/DataStructures.h"  // CoolProp::parameters

namespace CoolProp {
namespace expression {

/// Thermodynamic quantities a formula may reference are keyed by the existing
/// CoolProp::parameters enum -- one bucket, no DSL-private enum.  Whatever a
/// program asks for, the host fills by calling AbstractState::keyed_output() with
/// that key.  Some keys are free (T, rhomolar) and some cost an EOS call (p); the
/// evaluator does not care, and Program itself stays EOS-free -- it only reports
/// which keys it needs and reads back the values the host supplies.
///
/// The DSL spellings that resolve to a thermodynamic input, paired with the
/// CoolProp::parameters key each binds to, in name-resolution order.
///
/// This is deliberately a curated allowlist rather than an open door onto
/// get_parameter_index():
///  * the DSL spells state variables the way the C++ members do (`rhomolar`,
///    `rhomass`), which get_parameter_index() does not recognise -- it spells
///    them `Dmolar`/`Dmass` -- and those spellings are baked into fluid JSON;
///  * an open door would let a viscosity formula reference `V`, whose
///    keyed_output() re-enters the very correlation being defined.
const std::vector<std::pair<std::string, parameters>>& inputTable();

/// DSL spelling for `key`.  Throws CoolProp::ValueError if `key` is not in inputTable().
std::string inputName(parameters key);

namespace detail {
struct ProgramData;
}

/// An immutable, compiled expression. Cheap to copy (shared, refcounted body).
class Program
{
   public:
    /// Evaluate.  `inputVals` holds one value per entry of requiredInputs(), in
    /// that order; pass an empty vector when none are required.  A size mismatch
    /// throws CoolProp::ValueError rather than reading past the end.
    [[nodiscard]] double evaluate(const std::vector<double>& inputVals) const;
    /// Thermodynamic inputs this program references, in the order evaluate() expects them.
    [[nodiscard]] const std::vector<parameters>& requiredInputs() const;

   private:
    /// Throws CoolProp::ValueError if this Program was default-constructed (no
    /// compiled body); every public accessor goes through it before dereferencing.
    const detail::ProgramData& data() const;

    friend Program compile(const std::string&, const std::map<std::string, double>&, const std::map<std::string, std::vector<double>>&);
    std::shared_ptr<const detail::ProgramData> m_data;
};

/// Numeric-domain policy: division and `^`/`pow`/`log` follow IEEE 754 / `std::`
/// semantics and may yield `inf` or `nan` at domain edges (e.g., log(0), 0/0).
/// This is BY DESIGN so that compiled DSL expressions reproduce the existing
/// hardcoded C++ correlations bit-for-bit, including at domain edges.  No runtime
/// domain guards are added; adding them would diverge from the hardcoded routines.
///
/// Compile a formula string. `constants` are scalar names -> SI values; `arrays`
/// are vector names -> values. Throws CoolProp::ValueError on any lex/parse/bind error.
Program compile(const std::string& source, const std::map<std::string, double>& constants, const std::map<std::string, std::vector<double>>& arrays);

}  // namespace expression
}  // namespace CoolProp

#endif
