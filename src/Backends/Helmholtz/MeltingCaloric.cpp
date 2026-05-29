#include "MeltingCaloric.h"
// NOTE: match the include paths/style that FlashRoutines.cpp (same directory) uses
// for these two headers:
#include "HelmholtzEOSMixtureBackend.h"
#include "DataStructures.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <Eigen/Dense>

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

namespace {
using CE_t = std::vector<CoolProp::superancillary::ChebyshevExpansion<Eigen::ArrayXd>>;

void fit_part(CoolProp::HelmholtzEOSMixtureBackend& H, double ln_lo, double ln_hi,
              const std::function<double(CoolProp::HelmholtzEOSMixtureBackend&)>& pick, CE_t& out) {
    using namespace CoolProp;
    // Shrink interval by a tiny factor to avoid landing exactly on segment
    // boundaries where floating-point round-trip (double -> long double p_min)
    // can cause is_in_closed_range() to fail.
    const double span = ln_hi - ln_lo;
    const double eps_rel = 1e-9;
    const double adj_lo = ln_lo + eps_rel * span;
    const double adj_hi = ln_hi - eps_rel * span;
    auto func = [&](long /*i*/, long /*Npts*/, double lnp) -> double {
        const double p = std::exp(lnp);
        const double Tm = H.calc_melting_line(iT, iP, p);
        H.update(PT_INPUTS, p, Tm);
        return pick(H);
    };
    auto exps = superancillary::detail::dyadic_splitting<decltype(func), CE_t>(8, func, adj_lo, adj_hi, 3, 1e-10, 2);
    for (auto& e : exps) out.push_back(std::move(e));
}

/// Probe whether H can evaluate all caloric properties at lnp. Returns true on success.
bool can_eval_all(CoolProp::HelmholtzEOSMixtureBackend& H, double lnp) {
    try {
        const double p = std::exp(lnp);
        const double Tm = H.calc_melting_line(CoolProp::iT, CoolProp::iP, p);
        H.update(CoolProp::PT_INPUTS, p, Tm);
        return std::isfinite(H.T()) && std::isfinite(H.rhomolar()) && std::isfinite(H.hmolar())
               && std::isfinite(H.smolar());
    } catch (...) {
        return false;
    }
}

/// Filter pranges to only those parts where the EOS can be evaluated at both endpoints.
std::vector<std::pair<double, double>>
filter_valid_parts(CoolProp::HelmholtzEOSMixtureBackend& H,
                   const std::vector<std::pair<double, double>>& pranges) {
    std::vector<std::pair<double, double>> valid;
    for (const auto& pr : pranges) {
        const double ln_lo = std::log(pr.first), ln_hi = std::log(pr.second);
        if (!(ln_hi > ln_lo)) continue;
        const double eps = 1e-9 * (ln_hi - ln_lo);
        if (can_eval_all(H, ln_lo + eps) && can_eval_all(H, ln_hi - eps)) {
            valid.push_back(pr);
        }
    }
    return valid;
}

CoolProp::superancillary::ChebyshevApproximation1D<Eigen::ArrayXd>
fit_all_parts(CoolProp::HelmholtzEOSMixtureBackend& H,
              const std::vector<std::pair<double, double>>& pranges,
              const std::function<double(CoolProp::HelmholtzEOSMixtureBackend&)>& pick) {
    using namespace CoolProp;
    CE_t exps;
    for (const auto& pr : pranges) {
        const double ln_lo = std::log(pr.first), ln_hi = std::log(pr.second);
        if (!(ln_hi > ln_lo)) continue;
        fit_part(H, ln_lo, ln_hi, pick, exps);
    }
    std::sort(exps.begin(), exps.end(), [](const auto& a, const auto& b) { return a.xmin() < b.xmin(); });
    return superancillary::ChebyshevApproximation1D<Eigen::ArrayXd>(std::move(exps));
}
}  // namespace

namespace CoolProp {
void MeltingCaloric::build(HelmholtzEOSMixtureBackend& H) {
    if (m_built) return;
    if (!H.has_melting_line()) return;
    auto pranges = H.get_melting_line_part_pranges();
    const double p_lo_eos = H.calc_melting_line(iP_min, iT, 0.0);
    const double p_hi_eos = H.calc_melting_line(iP_max, iT, 0.0);
    for (auto& pr : pranges) {
        pr.first = std::max(pr.first, p_lo_eos);
        pr.second = std::min(pr.second, p_hi_eos);
    }
    pranges.erase(std::remove_if(pranges.begin(), pranges.end(),
                                 [](const std::pair<double, double>& pr) { return !(pr.second > pr.first); }),
                  pranges.end());
    if (pranges.empty()) return;

    // Filter to only parts where the EOS can evaluate all caloric properties at both endpoints.
    // Some melting-curve parts extend to pressures beyond the EOS validity domain.
    pranges = filter_valid_parts(H, pranges);
    if (pranges.empty()) return;

    m_T_approx.emplace(fit_all_parts(H, pranges, [](HelmholtzEOSMixtureBackend& S) { return S.T(); }));
    m_rho_approx.emplace(fit_all_parts(H, pranges, [](HelmholtzEOSMixtureBackend& S) { return S.rhomolar(); }));
    m_h_approx.emplace(fit_all_parts(H, pranges, [](HelmholtzEOSMixtureBackend& S) { return S.hmolar(); }));
    m_s_approx.emplace(fit_all_parts(H, pranges, [](HelmholtzEOSMixtureBackend& S) { return S.smolar(); }));
    m_built = true;
}
}  // namespace CoolProp
