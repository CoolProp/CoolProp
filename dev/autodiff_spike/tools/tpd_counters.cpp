// TPD-speed spike driver: per-mixture median PT-flash time + mean work counters per flash.
//   ./tpdspike <tag>   -> prints table to stdout, per-state results to tpd_<tag>.txt
#include "CoolProp/AbstractState.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <memory>
#include <random>
#include <string>
#include <vector>
namespace CoolProp { extern long spike_counts[16]; }
using CoolProp::spike_counts;
static double now_us() { return std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
int main(int argc, char** argv) {
    const std::string tag = argc > 1 ? argv[1] : "x";
    const int NS = 2000;
    struct Mix { std::string name, f; std::vector<double> x; double Tlo, Thi; };
    std::vector<Mix> mixes = {{"C1C2","Methane&Ethane",{0.5,0.5},150,400},
      {"C1C2C3","Methane&Ethane&Propane",{0.5,0.3,0.2},150,450},
      {"Amarillo","Methane&Nitrogen&CarbonDioxide&Ethane&Propane&IsoButane&n-Butane&Isopentane&n-Pentane&n-Hexane",
        {0.906724,0.031284,0.004676,0.045279,0.00828,0.001037,0.001563,0.000321,0.000443,0.000393},150,400},
      {"C1H2S","Methane&HydrogenSulfide",{0.5,0.5},180,450},
      {"humidair","Nitrogen&Oxygen&Argon&CarbonDioxide&Water",{0.7654,0.2053,0.0090,0.0003,0.02},250,450}};
    FILE* fo = std::fopen(("tpd_" + tag + ".txt").c_str(), "w");
    std::printf("%-9s %8s %8s | per flash: %6s %6s %6s %6s | %6s %6s | %6s %6s | %7s\n", "mixture", "med_us", "mean_us", "trial", "warmOK", "global", "kdir", "SSstab", "TPDit", "SSspl", "Gibbs", "alldrv");
    for (auto& mx : mixes) {
        if (argc > 2 && mx.name != argv[2]) continue;
        std::unique_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("GERG2008", mx.f));
        AS->set_mole_fractions(mx.x);
        std::mt19937_64 g(42);
        std::vector<double> Ts(NS), ps(NS), t(NS, 1e300);
        for (int i = 0; i < NS; ++i) {
            Ts[i] = mx.Tlo + (mx.Thi - mx.Tlo) * std::uniform_real_distribution<double>(0, 1)(g);
            ps[i] = std::exp(std::log(1e4) + (std::log(3e7) - std::log(1e4)) * std::uniform_real_distribution<double>(0, 1)(g));
        }
        std::vector<double> sum(16, 0.0);
        for (int pass = 0; pass < 2; ++pass)
            for (int i = 0; i < NS; ++i) {
                std::fill(spike_counts, spike_counts + 16, 0);
                double rho = NAN, Q = NAN;
                const double t0 = now_us();
                try { AS->update(CoolProp::PT_INPUTS, ps[i], Ts[i]); rho = AS->rhomolar(); Q = AS->Q(); } catch (...) {}
                t[i] = std::min(t[i], now_us() - t0);
                if (pass == 1) {
                    for (int k = 0; k < 16; ++k) sum[k] += spike_counts[k];
                    std::fprintf(fo, "%s %d %.10g %.10g %.12g %g\n", mx.name.c_str(), i, Ts[i], ps[i], rho, Q);
                }
            }
        std::vector<double> ts = t;
        std::nth_element(ts.begin(), ts.begin() + NS / 2, ts.end());
        double mean = 0; for (double v : t) mean += v; mean /= NS;
        std::printf("%-9s %8.1f %8.1f | per flash: %6.1f %6.1f %6.1f %6.1f | %6.1f %6.1f | %6.1f %6.1f | %7.1f\n", mx.name.c_str(), ts[NS / 2], mean,
                    sum[0] / NS, sum[1] / NS, sum[2] / NS, sum[3] / NS, sum[4] / NS, sum[5] / NS, sum[6] / NS, sum[7] / NS, sum[8] / NS);
        std::printf("%-9s  TPD calls/flash %.2f: converged %.2f, quick-unstable %.2f, step-fail %.2f, max-iter %.2f, density-fail %.2f; SS-verdict skips %.2f\n", "", sum[9]/NS, sum[10]/NS, sum[11]/NS, sum[12]/NS, sum[13]/NS, sum[14]/NS, sum[15]/NS);
        std::fflush(stdout);
    }
    std::fclose(fo);
}
