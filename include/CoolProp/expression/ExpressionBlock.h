#ifndef COOLPROP_EXPRESSION_BLOCK_H
#define COOLPROP_EXPRESSION_BLOCK_H

#include <string>
#include <vector>

#include "CoolProp/expression/Expression.h"

namespace CoolProp {

class AbstractState;  // forward decl keeps this header light

namespace expression {

/// Compile a transport `"type": "expression"` block from its JSON text:
///
///     {"formula": "...", "constants": {...}, "arrays": {...}}
///
/// Any other key -- `"type"` included -- is ignored, so the text can be lifted
/// verbatim out of (or pasted verbatim into) a fluid JSON file.  `context`, when
/// non-empty, is folded into the error message (the fluid library passes the
/// fluid name).  Throws CoolProp::ValueError on any JSON/lex/parse/bind error.
///
/// Text in, not nlohmann::json in, on purpose: this header is reached from the
/// scripting wrappers, and keeping it JSON-library-free keeps nlohmann out of
/// those translation units.
Program compile_block(const std::string& json_text, const std::string& context = "");

/// A standalone, compiled expression block, evaluable against any AbstractState
/// without having to graft it into a fluid file first.  This is the authoring
/// entry point the scripting wrappers expose: paste the same JSON block you would
/// put in the fluid JSON, evaluate it at a state, and see the number.
class ExpressionBlock
{
   public:
    /// `json_text` is the block as it appears in fluid JSON (see compile_block()).
    explicit ExpressionBlock(const std::string& json_text) : m_program(compile_block(json_text)) {}
    /// DSL names of the thermodynamic inputs the formula references, in the order
    /// they are fetched.  Empty for a formula built only from constants and arrays.
    [[nodiscard]] std::vector<std::string> required_inputs() const;
    /// Evaluate at the state `AS` is currently sitting at; the result is in whatever
    /// base-SI unit the formula produces.  The caller sets the state through the
    /// ordinary AbstractState API, so any input pair, backend, or mixture
    /// composition works, and nothing here duplicates state-setting logic:
    ///
    ///     std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", "R123"));
    ///     AS->update(DmolarT_INPUTS, rhomolar, T);
    ///     block.evaluate(*AS);
    ///
    /// Throws CoolProp::ValueError if an input reads back non-finite -- in practice,
    /// if `AS` was never update()d.
    [[nodiscard]] double evaluate(AbstractState& AS) const;

   private:
    Program m_program;
};

}  // namespace expression
}  // namespace CoolProp

#endif
