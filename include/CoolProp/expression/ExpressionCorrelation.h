#ifndef COOLPROP_EXPRESSION_CORRELATION_H
#define COOLPROP_EXPRESSION_CORRELATION_H

#include <memory>

#include "CoolProp/expression/Expression.h"

namespace CoolProp {

class AbstractState;  // forward decls keep this header light
class HelmholtzEOSMixtureBackend;

namespace expression {

/// Fill `prog`'s required inputs from the state `AS` is currently sitting at --
/// one AbstractState::keyed_output() per key, in requiredInputs() order -- and
/// evaluate.  This is THE host-side evaluation path: the fluid-library
/// correlation below and the standalone expression::ExpressionBlock both go
/// through it, so neither can drift from the other.
///
/// Throws CoolProp::ValueError if an input reads back non-finite.  That guards the
/// INPUTS only; inside the formula, pow/log domain edges still yield inf/nan
/// exactly as the hardcoded routines do.
double evaluate_at(const Program& prog, AbstractState& AS);

/// Host-side wrapper: owns a compiled Program, fetches every thermodynamic input
/// the program asked for from an EOS backend (one keyed_output() per key), then
/// evaluates.
class ExpressionCorrelation
{
   public:
    ExpressionCorrelation() = default;
    explicit ExpressionCorrelation(Program prog) : m_program(std::move(prog)) {}
    /// Evaluate the formula at the backend's current state; returns base-SI result.
    /// A default-constructed correlation throws CoolProp::ValueError (its Program
    /// has no compiled body); "is there a correlation here?" is answered by the
    /// owning ExpressionData::correlation shared_ptr, which the dispatch checks.
    [[nodiscard]] double eval(HelmholtzEOSMixtureBackend& HEOS) const;

   private:
    Program m_program{};
};

}  // namespace expression

/// Stored in each transport sub-correlation container. Empty unless type==expression.
struct ExpressionData
{
    std::shared_ptr<expression::ExpressionCorrelation> correlation{};
};

}  // namespace CoolProp

#endif
