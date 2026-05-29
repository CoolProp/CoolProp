#include "MeltingCaloric.h"
// NOTE: match the include paths/style that FlashRoutines.cpp (same directory) uses
// for these two headers:
#include "HelmholtzEOSMixtureBackend.h"
#include "DataStructures.h"
#include <cmath>
#include <stdexcept>

namespace CoolProp {

void MeltingCaloric::sample(HelmholtzEOSMixtureBackend& H, std::size_t n) {
    if (!H.has_melting_line()) {
        throw ValueError("MeltingCaloric::sample: fluid has no melting line");
    }
    if (n < 8) n = 8;
    const double p_lo = H.calc_melting_line(iP_min, iT, 0.0);
    const double p_hi = H.calc_melting_line(iP_max, iT, 0.0);
    const double lo = std::log(p_lo), hi = std::log(p_hi);

    m_lnp.clear(); m_T.clear(); m_rho.clear(); m_h.clear(); m_s.clear();
    m_lnp.reserve(n); m_T.reserve(n); m_rho.reserve(n); m_h.reserve(n); m_s.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double lnp = lo + (hi - lo) * static_cast<double>(i) / static_cast<double>(n - 1);
        const double p = std::exp(lnp);
        double Tm;
        try {
            Tm = H.calc_melting_line(iT, iP, p);
        } catch (...) {
            continue;
        }
        try {
            H.update(PT_INPUTS, p, Tm);
        } catch (...) {
            continue;
        }
        const double rho = H.rhomolar(), h = H.hmolar(), s = H.smolar();
        if (!std::isfinite(rho) || !std::isfinite(h) || !std::isfinite(s)) continue;
        m_lnp.push_back(lnp); m_T.push_back(Tm); m_rho.push_back(rho);
        m_h.push_back(h); m_s.push_back(s);
    }
}

}  // namespace CoolProp
