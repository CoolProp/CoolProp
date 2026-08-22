#include "CoolProp/expression/ExpressionCorrelation.h"

#include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"

namespace CoolProp {
namespace expression {

double ExpressionCorrelation::eval(HelmholtzEOSMixtureBackend& HEOS) const {
    // One bucket: every input the program asked for is a CoolProp::parameters key,
    // so the whole host side is a keyed_output() per key -- no per-quantity switch
    // to extend when a new input name is added to the table in Expression.cpp.
    const std::vector<parameters>& keys = m_program.requiredInputs();
    std::vector<double> vals(keys.size());
    for (std::size_t i = 0; i < keys.size(); ++i)
        vals[i] = HEOS.keyed_output(keys[i]);
    return m_program.evaluate(vals);
}

}  // namespace expression
}  // namespace CoolProp
