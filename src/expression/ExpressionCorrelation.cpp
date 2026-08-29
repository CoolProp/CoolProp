#include "CoolProp/expression/ExpressionCorrelation.h"

#include <vector>

#include "CoolProp/AbstractState.h"
#include "CoolProp/Exceptions.h"
#include "CoolProp/detail/strings.h"
#include "CoolProp/numerics/numerics.h"  // ValidNumber
#include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"

namespace CoolProp {
namespace expression {

double evaluate_at(const Program& prog, AbstractState& AS) {
    // One bucket: every input the program asked for is a CoolProp::parameters key,
    // so the whole host side is a keyed_output() per key -- no per-quantity switch
    // to extend when a block declares a new state variable.
    const std::vector<parameters>& keys = prog.requiredInputs();
    std::vector<double> vals(keys.size());
    for (std::size_t i = 0; i < keys.size(); ++i) {
        vals[i] = AS.keyed_output(keys[i]);
        // An AbstractState that was never update()d hands back +/-_HUGE rather than
        // throwing, and the formula would happily propagate it into a plausible
        // -inf result.  Fail here, naming the input, instead of downstream.  (This
        // guards the INPUTS only; the DSL's own numeric policy deliberately lets
        // pow/log produce inf/nan mid-formula to match the hardcoded routines.)
        if (!ValidNumber(vals[i])) {
            throw ValueError(
              format("expression input '%s' is %g, not a finite number -- has the state been updated?", inputName(keys[i]).c_str(), vals[i]));
        }
    }
    return prog.evaluate(vals);
}

double ExpressionCorrelation::eval(HelmholtzEOSMixtureBackend& HEOS) const {
    return evaluate_at(m_program, HEOS);
}

}  // namespace expression
}  // namespace CoolProp
