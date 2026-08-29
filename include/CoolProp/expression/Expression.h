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
/// Resolve `name` as a state variable a correlation may declare.
///
/// The DSL does NOT maintain its own list of thermodynamic quantities.  A name is
/// resolved by CoolProp::is_valid_parameter(), so the DSL's vocabulary IS CoolProp's
/// vocabulary -- every quantity keyed_output() can produce is reachable, and adding
/// a new one never requires touching this library.  The spellings are CoolProp's
/// canonical ones (`P`, `Dmolar`, `Dmass`); the DSL invents no aliases of its own,
/// because an invented lowercase `p` is precisely what once collided with the
/// exponent array every viscosity paper calls p_i.
///
/// Returns true and sets `key` when `name` is resolvable and permitted.  Otherwise
/// returns false and sets `reason` to a message fit for a compile error.  Two
/// classes are refused even though CoolProp resolves them:
///
///  * transport outputs (`V`, `L`, `Prandtl`, ...) -- keyed_output() for these
///    re-enters the very correlation being defined;
///  * the critical point and the EOS reducing state -- calc_T_critical() and
///    calc_rhomolar_critical() return the NUMERICAL critical point when
///    ENABLE_SUPERANCILLARIES is set, so a correlation reducing on them would
///    silently change answer with configuration.  Freeze those as `constants`
///    instead, at the value the paper's authors regressed against.
bool resolveStateVariable(const std::string& name, parameters& key, std::string& reason);

/// CoolProp's canonical spelling for `key`, as used in error messages and in
/// ExpressionBlock::required_inputs().  Since only the canonical spelling is
/// accepted as a declaration, this round-trips whatever the author wrote.  Throws
/// CoolProp::ValueError for a key outside the parameter-information table.
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
    /// Thermodynamic inputs this program references, in the order evaluate() expects
    /// them: the block's DECLARED `state_variables`, in order of FIRST REFERENCE in
    /// the formula -- which is not necessarily the order they were declared in.  The
    /// declaration gates what may be read; the formula fixes the order.
    [[nodiscard]] const std::vector<parameters>& requiredInputs() const;

   private:
    /// Throws CoolProp::ValueError if this Program was default-constructed (no
    /// compiled body); every public accessor goes through it before dereferencing.
    const detail::ProgramData& data() const;

    friend Program compile(const std::string&, const std::map<std::string, double>&, const std::map<std::string, std::vector<double>>&,
                           const std::vector<std::string>&);
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
/// `state_variables` names the thermodynamic quantities this formula is allowed to
/// read, in CoolProp's spelling.  It is opt-in and per-formula: a name not declared
/// here is never state, so a block that does not ask for pressure keeps `P` -- and,
/// more to the point, keeps every lowercase coefficient name -- for its own use.
/// Declaring a name the formula never reads is an error, as is declaring one that
/// also appears in `constants` or `arrays`.
Program compile(const std::string& source, const std::map<std::string, double>& constants, const std::map<std::string, std::vector<double>>& arrays,
                const std::vector<std::string>& state_variables = {});

}  // namespace expression
}  // namespace CoolProp

#endif
