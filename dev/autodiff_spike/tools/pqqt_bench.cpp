// PQ / QT mixture flash: time per call, failures, solver fallbacks (warnstring), and the published state's
// mass-balance and equal-fugacity errors.  One summary line per (mixture, input pair).
//   ./pqbench
#include "CoolProp/AbstractState.h"
#include "CoolProp/CoolProp.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <cstdlib>
#include <vector>
using namespace CoolProp;
static double now_us() { return std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
int main() {
    struct Mix { std::string name, be, f; std::vector<double> z; double plo, phi, Tlo, Thi; };
    const std::vector<Mix> mixes = {
      {"C1/C2 50/50", "HEOS", "Methane&Ethane", {0.5, 0.5}, 2e4, 4e6, 150, 250},
      {"C1/C2/C3 50/30/20", "HEOS", "Methane&Ethane&Propane", {0.5, 0.3, 0.2}, 2e4, 5e6, 150, 290},
      {"NG5 N2/C1/C2/C3/nC4", "HEOS", "Nitrogen&Methane&Ethane&Propane&n-Butane", {0.03, 0.85, 0.07, 0.035, 0.015}, 2e4, 4e6, 130, 230},
      {"R454B", "HEOS", "R32&R1234yf", {0.8292, 0.1708}, 2e4, 4e6, 220, 330}};
    const std::vector<double> Qs = {0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95};
    const int NX = 12;
    std::printf("%-22s %-3s | %5s %6s %6s | %9s %9s | %10s %10s\n", "mixture", "in", "n", "fail", "fallbk", "med_us", "mean_us", "max|mb|", "max|dlnf|");
    for (const auto& m : mixes) {
        for (int pair = 0; pair < 2; ++pair) {
            std::unique_ptr<AbstractState> AS(AbstractState::factory(m.be, m.f));
            AS->set_mole_fractions(m.z);
            std::vector<double> t;
            int nfail = 0, nfb = 0;
            double mbmax = 0, dfmax = 0;
            for (int k = 0; k < NX; ++k) {
                const double s = (k + 0.5) / NX;
                const double X = pair == 0 ? std::exp(std::log(m.plo) + s * (std::log(m.phi) - std::log(m.plo))) : m.Tlo + s * (m.Thi - m.Tlo);
                for (double Q : Qs) {
                    get_global_param_string("warnstring");  // clear
                    double best = 1e300;
                    bool ok = true;
                    for (int rep = 0; rep < 3; ++rep) {
                        const double t0 = now_us();
                        try { AS->update(pair == 0 ? PQ_INPUTS : QT_INPUTS, pair == 0 ? X : Q, pair == 0 ? Q : X); } catch (...) { ok = false; }
                        best = std::min(best, now_us() - t0);
                    }
                    if (!ok) { ++nfail; continue; }
                    t.push_back(best);
                    const std::string warn = get_global_param_string("warnstring");
                    if (!warn.empty()) ++nfb;
                    double mb_i = 0, df_i = 0;
                    try {
                        const auto x = AS->mole_fractions_liquid_double(), y = AS->mole_fractions_vapor_double();
                        const double q = AS->Q();
                        for (std::size_t i = 0; i < m.z.size(); ++i) mb_i = std::max(mb_i, std::abs(m.z[i] - ((1 - q) * x[i] + q * y[i])));
                        mbmax = std::max(mbmax, mb_i);
                        std::unique_ptr<AbstractState> L(AbstractState::factory(m.be, m.f)), V(AbstractState::factory(m.be, m.f));
                        L->set_mole_fractions(x); V->set_mole_fractions(y);
                        L->update(DmolarT_INPUTS, AS->saturated_liquid_keyed_output(iDmolar), AS->T());
                        V->update(DmolarT_INPUTS, AS->saturated_vapor_keyed_output(iDmolar), AS->T());
                        for (std::size_t i = 0; i < m.z.size(); ++i) df_i = std::max(df_i, std::abs(std::log(V->fugacity(i) / L->fugacity(i))));
                        dfmax = std::max(dfmax, df_i);
                        if (std::abs(q - Q) > 1e-9) mbmax = std::max(mbmax, 1.0);  // published Q is not the requested Q
                    } catch (std::exception& e) { dfmax = std::max(dfmax, 1.0); df_i = 1.0; std::fprintf(stderr, "   check threw: %s\n", e.what()); }
                    if ((mb_i > 1e-6 || df_i > 1e-6 || !warn.empty()) && std::getenv("PQ_DETAIL"))
                        std::fprintf(stderr, "   %s %s X=%.6g Q=%.2f -> T=%.4f p=%.6g Q=%.6g rhoL=%.6g rhoV=%.6g mb=%.2e dlnf=%.2e %s\n", m.name.c_str(), pair == 0 ? "PQ" : "QT", X, Q, AS->T(), AS->p(), AS->Q(),
                                     AS->saturated_liquid_keyed_output(iDmolar), AS->saturated_vapor_keyed_output(iDmolar), mb_i, df_i, warn.empty() ? "" : ("WARN: " + warn.substr(0, 120)).c_str());
                }
            }
            std::vector<double> ts = t;
            std::nth_element(ts.begin(), ts.begin() + ts.size() / 2, ts.end());
            double mean = 0; for (double v : t) mean += v; mean /= std::max<std::size_t>(t.size(), 1);
            std::printf("%-22s %-3s | %5d %6d %6d | %9.1f %9.1f | %10.2e %10.2e\n", m.name.c_str(), pair == 0 ? "PQ" : "QT", NX * (int)Qs.size(), nfail, nfb,
                        ts.empty() ? NAN : ts[ts.size() / 2], mean, mbmax, dfmax);
            std::fflush(stdout);
        }
    }
}
