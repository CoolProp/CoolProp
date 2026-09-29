// Sweep the bench_tpflash states (one pass, CoolProp only) and tag each state on stderr.
#include "CoolProp/AbstractState.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <memory>
#include <random>
#include <string>
#include <vector>
int main(int argc, char** argv) {
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
        std::unique_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("GERG2008", mx.f));
        AS->set_mole_fractions(mx.x);
        std::mt19937_64 g(42);
        std::vector<double> Ts(NS), ps(NS);
        for (int i = 0; i < NS; ++i) {
            Ts[i] = mx.Tlo + (mx.Thi - mx.Tlo) * std::uniform_real_distribution<double>(0, 1)(g);
            ps[i] = std::exp(std::log(1e4) + (std::log(3e7) - std::log(1e4)) * std::uniform_real_distribution<double>(0, 1)(g));
        }
        for (int i = 0; i < NS; ++i) {
            if (argc > 2 && (mx.name != argv[1] || i != std::atoi(argv[2]))) {
                if (argc > 3) continue;
            }
            std::fprintf(stderr, "STATE %s %d T=%.10g p=%.10g\n", mx.name.c_str(), i, Ts[i], ps[i]);
            try {
                AS->update(CoolProp::PT_INPUTS, ps[i], Ts[i]);
            } catch (...) {
            }
        }
    }
}
