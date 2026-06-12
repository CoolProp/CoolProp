#ifndef COOLPROP_EXPRESSION_CORRELATION_H
#define COOLPROP_EXPRESSION_CORRELATION_H

#include <memory>

#include "CoolProp/expression/Expression.h"

namespace CoolProp {

class HelmholtzEOSMixtureBackend;  // forward decl keeps this header light

namespace expression {

/// Host-side wrapper: owns a compiled Program and fetches the intrinsic state
/// and registered derived quantities from an EOS backend, then evaluates.
class ExpressionCorrelation
{
   public:
    ExpressionCorrelation() = default;
    explicit ExpressionCorrelation(Program prog) : m_program(std::move(prog)), m_set(true) {}
    [[nodiscard]] bool is_set() const { return m_set; }
    /// Evaluate the formula at the backend's current state; returns base-SI result.
    [[nodiscard]] double eval(HelmholtzEOSMixtureBackend& HEOS) const;

   private:
    Program m_program{};
    bool m_set = false;
};

}  // namespace expression

/// Stored in each transport sub-correlation container. Empty unless type==expression.
struct ExpressionData
{
    std::shared_ptr<expression::ExpressionCorrelation> correlation{};
};

}  // namespace CoolProp

#endif
