#include "CoolProp/expression/Expression.h"
#include "CoolProp/Exceptions.h"
#include "CoolProp/detail/strings.h"

namespace CoolProp {
namespace expression {
namespace detail {
struct ProgramData
{
    std::vector<Intrinsic> intrinsics;
    std::vector<Derived> deriveds;
};
}  // namespace detail

double Program::evaluate(const double*, const double*) const {
    // Stub: real implementation lands in Task 5.
    return 7.0;
}
const std::vector<Intrinsic>& Program::requiredIntrinsics() const {
    return m_data->intrinsics;
}
const std::vector<Derived>& Program::requiredDerived() const {
    return m_data->deriveds;
}

Program compile(const std::string&, const std::map<std::string, double>&,
                const std::map<std::string, std::vector<double>>&) {
    Program p;
    p.m_data = std::make_shared<detail::ProgramData>();
    return p;
}

}  // namespace expression
}  // namespace CoolProp
