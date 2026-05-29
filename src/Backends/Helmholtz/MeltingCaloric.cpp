#include "MeltingCaloric.h"
// NOTE: match the include paths/style that FlashRoutines.cpp (same directory) uses
// for these two headers:
#include "HelmholtzEOSMixtureBackend.h"
#include "DataStructures.h"
#include "FlashRoutines.h"
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
/// Fallback bracket-scan + bisection on h(lnp) - h_cache.  Used only when the
/// Chebyshev monotone-interval inverter returns nothing (should not happen for
/// water since h is monotonic, but kept for robustness with other fluids).
double scan_bisect_h(const CoolProp::superancillary::ChebyshevApproximation1D<Eigen::ArrayXd>& approx,
                     double h_cache) {
    const double lo = approx.xmin(), hi = approx.xmax();
    constexpr int Nscan = 200;
    double best_lnp = lo, best_gap = 1e300;
    double prev_lnp = lo, prev_h = approx.eval(lo) - h_cache;
    for (int k = 1; k <= Nscan; ++k) {
        const double cur_lnp = lo + (hi - lo) * k / Nscan;
        const double cur_h = approx.eval(cur_lnp) - h_cache;
        if (prev_h * cur_h <= 0.0) {
            double a = prev_lnp, b = cur_lnp, fa = prev_h, fb = cur_h;
            for (int it = 0; it < 60; ++it) {
                const double m = 0.5 * (a + b);
                const double fm = approx.eval(m) - h_cache;
                if (fa * fm <= 0.0) { b = m; fb = fm; } else { a = m; fa = fm; }
                if ((b - a) < 1e-14 * (std::abs(a) + std::abs(b) + 1e-30)) break;
            }
            const double root_lnp = 0.5 * (a + b);
            const double gap = std::abs(approx.eval(root_lnp) - h_cache);
            if (std::isfinite(gap) && gap < best_gap) { best_gap = gap; best_lnp = root_lnp; }
        }
        prev_lnp = cur_lnp; prev_h = cur_h;
    }
    return best_lnp;
}
}  // namespace

namespace CoolProp {
bool MeltingCaloric::seed_for_hs(double s_cache, double h_cache, double& T0, double& rho0) const {
    if (!m_built || !m_s_approx || !m_h_approx || !m_T_approx || !m_rho_approx) return false;
    // Enthalpy is monotonic in ln(p) along the melting curve (entropy and T are
    // NOT -- they fold back at the ice-Ih/III junction), so intersect on h: this
    // yields a single, well-conditioned root via the Chebyshev monotone-interval
    // inverter. Entropy is used only to disambiguate in the (not expected) event
    // of multiple roots.
    std::vector<double> cand_lnp;
    for (const auto& pr : m_h_approx->get_x_for_y(h_cache, 48, 100, 1e-12)) {
        const double lnp = pr.first;
        if (lnp >= m_h_approx->xmin() && lnp <= m_h_approx->xmax()) cand_lnp.push_back(lnp);
    }
    if (cand_lnp.empty()) {
        // Fallback: bracket-scan h(lnp) - h_cache for a sign change (robust if the
        // assembled-expansion metadata ever misbehaves or h is out of range).
        const double fb_lnp = scan_bisect_h(*m_h_approx, h_cache);
        const double fb_gap = std::abs(m_h_approx->eval(fb_lnp) - h_cache);
        // Only accept the fallback if it actually found a near-zero residual
        if (std::isfinite(fb_gap) && fb_gap < std::abs(h_cache) * 0.01 + 1.0) {
            cand_lnp.push_back(fb_lnp);
        }
    }
    if (cand_lnp.empty()) return false;
    // Disambiguate by entropy closeness (usually exactly one candidate).
    double best_lnp = cand_lnp.front(), best_gap = 1e300;
    bool found = false;
    for (double lnp : cand_lnp) {
        const double sgap = std::abs(m_s_approx->eval(lnp) - s_cache);
        if (std::isfinite(sgap) && sgap < best_gap) { best_gap = sgap; best_lnp = lnp; found = true; }
    }
    if (!found) return false;
    T0 = m_T_approx->eval(best_lnp);
    rho0 = m_rho_approx->eval(best_lnp);
    return std::isfinite(T0) && std::isfinite(rho0) && rho0 > 0;
}

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
    m_stamp = FlashRoutines::alpha0_offset_total(H);
    m_built = true;
}
}  // namespace CoolProp
