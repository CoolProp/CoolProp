#include "CoolProp/expression/ExpressionCorrelation.h"

#include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"

namespace CoolProp {
namespace expression {

static double fetchIntrinsic(HelmholtzEOSMixtureBackend& H, Intrinsic k) {
    switch (k) {
        case Intrinsic::T: return H.T();
        case Intrinsic::rhomolar: return H.rhomolar();
        case Intrinsic::rhomass: return H.rhomass();
        case Intrinsic::molar_mass: return H.molar_mass();
        default: throw ValueError("internal: unknown intrinsic");
    }
}
static double fetchDerived(HelmholtzEOSMixtureBackend& H, Derived k) {
    switch (k) {
        case Derived::p: return H.p();
        default: throw ValueError("internal: unknown derived quantity");
    }
}

double ExpressionCorrelation::eval(HelmholtzEOSMixtureBackend& HEOS) const {
    const std::vector<Intrinsic>& ik = m_program.requiredIntrinsics();
    const std::vector<Derived>& dk = m_program.requiredDerived();
    std::vector<double> iv(ik.size()), dv(dk.size());
    for (std::size_t i = 0; i < ik.size(); ++i) iv[i] = fetchIntrinsic(HEOS, ik[i]);
    for (std::size_t i = 0; i < dk.size(); ++i) dv[i] = fetchDerived(HEOS, dk[i]);
    return m_program.evaluate(iv.empty() ? nullptr : iv.data(), dv.empty() ? nullptr : dv.data());
}

}  // namespace expression
}  // namespace CoolProp
