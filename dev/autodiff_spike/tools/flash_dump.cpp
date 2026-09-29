// Sweep the bench_tpflash states through the PT flash: per-state result on stdout, median time per mixture on stderr.
//   ./rhodump2 <backend>
#include "CoolProp/AbstractState.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <memory>
#include <random>
#include <string>
#include <vector>
static double now_us() {
    return std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
int main(int argc, char** argv) {
    const std::string be = argc > 1 ? argv[1] : "GERG2008";
    const int NS = 2000;
    std::vector<std::string> NG10 = {"Methane",   "Nitrogen", "CarbonDioxide", "Ethane",    "Propane",
                                     "IsoButane", "n-Butane", "Isopentane",    "n-Pentane", "n-Hexane"};
    std::vector<double> AMA = {0.906724, 0.031284, 0.004676, 0.045279, 0.00828, 0.001037, 0.001563, 0.000321, 0.000443, 0.000393};
    struct Mix
    {
        std::string name, f;
        std::vector<double> x;
        double Tlo, Thi;
    };
    std::vector<Mix> mixes = {
      {"C1C2", "Methane&Ethane", {0.5, 0.5}, 150, 400},
      {"C1C2C3", "Methane&Ethane&Propane", {0.5, 0.3, 0.2}, 150, 450},
      {"Amarillo", "Methane&Nitrogen&CarbonDioxide&Ethane&Propane&IsoButane&n-Butane&Isopentane&n-Pentane&n-Hexane", AMA, 150, 400},
      {"C1H2S", "Methane&HydrogenSulfide", {0.5, 0.5}, 180, 450},
      {"humidair", "Nitrogen&Oxygen&Argon&CarbonDioxide&Water", {0.7654, 0.2053, 0.0090, 0.0003, 0.02}, 250, 450}};
    for (auto& mx : mixes) {
        std::unique_ptr<CoolProp::AbstractState> AS;
        try {
            AS.reset(CoolProp::AbstractState::factory(be, mx.f));
        } catch (std::exception& e) {
            std::fprintf(stderr, "%s %s: factory failed: %s\n", be.c_str(), mx.name.c_str(), e.what());
            continue;
        }
        AS->set_mole_fractions(mx.x);
        std::mt19937_64 g(42);
        std::vector<double> Ts(NS), ps(NS), t(NS);
        for (int i = 0; i < NS; ++i) {
            Ts[i] = mx.Tlo + (mx.Thi - mx.Tlo) * std::uniform_real_distribution<double>(0, 1)(g);
            ps[i] = std::exp(std::log(1e4) + (std::log(3e7) - std::log(1e4)) * std::uniform_real_distribution<double>(0, 1)(g));
        }
        int nfail = 0;
        for (int i = 0; i < NS; ++i) {
            double rho = NAN, Q = NAN;
            const double t0 = now_us();
            try {
                AS->update(CoolProp::PT_INPUTS, ps[i], Ts[i]);
                rho = AS->rhomolar();
                Q = AS->Q();
            } catch (...) {
                ++nfail;
            }
            t[i] = now_us() - t0;
            std::printf("%s %s %d %.10g %.10g %.12g %g\n", be.c_str(), mx.name.c_str(), i, Ts[i], ps[i], rho, Q);
        }
        std::nth_element(t.begin(), t.begin() + NS / 2, t.end());
        std::fprintf(stderr, "%-9s %-9s median %8.1f us  fails %d\n", be.c_str(), mx.name.c_str(), t[NS / 2], nfail);
    }
}
