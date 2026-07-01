#ifndef COOLPROP_EXPRESSION_H
#define COOLPROP_EXPRESSION_H

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace CoolProp {
namespace expression {

/// EOS-free independent state and pure fluid metadata, filled by the host per eval.
enum class Intrinsic : std::uint8_t { T, rhomolar, rhomass, molar_mass };
/// State-dependent quantities the EOS must compute; v1 registry holds only p.
enum class Derived : std::uint8_t { p };

namespace detail {
struct ProgramData;
}

/// An immutable, compiled expression. Cheap to copy (shared, refcounted body).
class Program
{
   public:
    /// Evaluate. intrinsicVals/derivedVals are arrays in the order given by
    /// requiredIntrinsics()/requiredDerived(); pass nullptr when none are required.
    [[nodiscard]] double evaluate(const double* intrinsicVals, const double* derivedVals) const;
    /// Intrinsics this program references, in the order evaluate() expects them.
    [[nodiscard]] const std::vector<Intrinsic>& requiredIntrinsics() const;
    /// Derived quantities this program references, in the order evaluate() expects them.
    [[nodiscard]] const std::vector<Derived>& requiredDerived() const;

   private:
    friend Program compile(const std::string&, const std::map<std::string, double>&,
                           const std::map<std::string, std::vector<double>>&);
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
Program compile(const std::string& source, const std::map<std::string, double>& constants,
                const std::map<std::string, std::vector<double>>& arrays);

}  // namespace expression
}  // namespace CoolProp

#endif
