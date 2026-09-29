// Judge REFPROP TPFLSH two-phase answers with REFPROP's own routines only (no CoolProp).
#include <dlfcn.h>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
using SETPATH_t = void (*)(char*, long);
using SETUP_t = void (*)(int*, char*, char*, char*, int*, char*, long, long, long, long);
using FLAGS_t = void (*)(char*, int*, int*, int*, char*, long, long);
using TPFLSH_t = void (*)(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*,
                          double*, double*, int*, char*, long);
using FG_t = void (*)(double*, double*, double*, double*, int*, char*, long);
using PR_t = void (*)(double*, double*, double*, double*);
SETPATH_t SP;
SETUP_t SU;
FLAGS_t FL;
TPFLSH_t TP;
FG_t FG;
PR_t PR;
std::map<std::string, std::string> RN = {{"Methane", "METHANE"},   {"Ethane", "ETHANE"},      {"Propane", "PROPANE"}, {"HydrogenSulfide", "H2S"},
                                         {"Nitrogen", "NITROGEN"}, {"Oxygen", "OXYGEN"},      {"Argon", "ARGON"},     {"CarbonDioxide", "CO2"},
                                         {"Water", "WATER"},       {"IsoButane", "ISOBUTAN"}, {"n-Butane", "BUTANE"}, {"Isopentane", "IPENTANE"},
                                         {"n-Pentane", "PENTANE"}, {"n-Hexane", "HEXANE"},    {"Helium", "HELIUM"},   {"Hydrogen", "HYDROGEN"}};
void setup(const std::vector<std::string>& f) {
    char h0[256] = "GERG", e0[256] = {};
    int j = 1, k = 0, ie = 0;
    FL(h0, &j, &k, &ie, e0, 255, 255);
    std::string s;
    for (size_t i = 0; i < f.size(); ++i)
        s += (i ? "|" : "") + RN.at(f[i]) + ".FLD";
    std::vector<char> hf(10000, 0), he(255, 0);
    memcpy(hf.data(), s.data(), s.size());
    char hmx[255] = "HMX.BNC", hrf[3] = {'D', 'E', 'F'};
    int n = f.size(), ierr = 0;
    SU(&n, hf.data(), hmx, hrf, &ierr, he.data(), 10000, 255, 3, 255);
}
int main() {
    void* h = dlopen("/Users/ianbell/REFPROP10/librefprop.dylib", RTLD_NOW);
    SP = (SETPATH_t)dlsym(h, "SETPATHdll");
    SU = (SETUP_t)dlsym(h, "SETUPdll");
    FL = (FLAGS_t)dlsym(h, "FLAGSdll");
    TP = (TPFLSH_t)dlsym(h, "TPFLSHdll");
    FG = (FG_t)dlsym(h, "FGCTY2dll");
    PR = (PR_t)dlsym(h, "PRESSdll");
    char path[256] = "/Users/ianbell/REFPROP10/";
    SP(path, 255);
    std::vector<std::string> NG10 = {"Methane",   "Nitrogen", "CarbonDioxide", "Ethane",    "Propane",
                                     "IsoButane", "n-Butane", "Isopentane",    "n-Pentane", "n-Hexane"};
    std::map<std::string, std::pair<std::vector<std::string>, std::vector<double>>> M = {
      {"C1/C2 50/50", {{"Methane", "Ethane"}, {0.5, 0.5}}},
      {"C1/C2/C3 50/30/20", {{"Methane", "Ethane", "Propane"}, {0.5, 0.3, 0.2}}},
      {"Amarillo (10)", {NG10, {0.906724, 0.031284, 0.004676, 0.045279, 0.00828, 0.001037, 0.001563, 0.000321, 0.000443, 0.000393}}},
      {"C1/H2S 50/50", {{"Methane", "HydrogenSulfide"}, {0.5, 0.5}}},
      {"humid air x_w=0.02", {{"Nitrogen", "Oxygen", "Argon", "CarbonDioxide", "Water"}, {0.7654, 0.2053, 0.0090, 0.0003, 0.02}}},
      {"N2/C1/C2/nC4/nC5", {{"Nitrogen", "Methane", "Ethane", "n-Butane", "n-Pentane"}, {0.3797, 0.3225, 0.278, 0.0014, 0.0184}}},
      {"CO2/H2 96.77/3.23", {{"CarbonDioxide", "Hydrogen"}, {0.9677, 0.0323}}},
      {"CO2/N2/O2/He",
       {{"CarbonDioxide", "Nitrogen", "Oxygen", "Helium"}, {0.9403 / 1.00002, 0.0582 / 1.00002, 0.00127 / 1.00002, 0.00025 / 1.00002}}},
      {"N2/O2/Ar (O2-enriched air)", {{"Nitrogen", "Oxygen", "Argon"}, {0.609067, 0.370414, 0.0205193}}}};
    std::ifstream in("rpwrong_states2.tsv");
    std::string line, cur;
    std::map<std::string, std::vector<int>> bins;  // per mixture: [dlnf<1e-8, <1e-6, <1e-5, <1e-3, >=1e-3, not reproduced single-phase]
    while (std::getline(in, line)) {
        std::stringstream ss(line);
        std::string name, Ts, ps, qs;
        std::getline(ss, name, '\t');
        std::getline(ss, Ts, '\t');
        std::getline(ss, ps, '\t');
        std::string is;
        std::getline(ss, is, '\t');
        if (name != cur) {
            setup(M[name].first);
            cur = name;
            bins[name].assign(7, 0);
        }
        auto z = M[name].second;
        size_t N = z.size();
        z.resize(20, 0);
        double T = std::stod(Ts), p = std::stod(ps) / 1000, D, Dl, Dv, xl[20] = {}, yv[20] = {}, q, e, hh, s, cv, cp, w;
        int ie = 0;
        char he[256];
        TP(&T, &p, z.data(), &D, &Dl, &Dv, xl, yv, &q, &e, &hh, &s, &cv, &cp, &w, &ie, he, 255);
        if (!(q >= 0 && q <= 1)) {
            bins[name][5]++;
            fprintf(stderr, "%s\t%s\tnot2ph\n", name.c_str(), is.c_str());
            continue;
        }
        double fl[20], fv[20], pl, pv;
        int e1 = 0, e2 = 0;
        char h1[256], h2[256];
        FG(&T, &Dl, xl, fl, &e1, h1, 255);
        FG(&T, &Dv, yv, fv, &e2, h2, 255);
        PR(&T, &Dl, xl, &pl);
        PR(&T, &Dv, yv, &pv);
        double dl = 0;
        for (size_t i = 0; i < N; ++i)
            if (xl[i] > 0 || yv[i] > 0) dl = std::max(dl, std::fabs(std::log(fv[i] / fl[i])));
        double dp = std::max(std::fabs(pl / p - 1), std::fabs(pv / p - 1));
        double worst = std::max(dl, dp);
        int b = worst < 1e-8 ? 0 : worst < 1e-6 ? 1 : worst < 1e-5 ? 2 : worst < 1e-3 ? 3 : 4;
        bins[name][b]++;
        fprintf(stderr, "%s\t%s\t%.3e\n", name.c_str(), is.c_str(), worst);
        if (worst < 1e-5) bins[name][6]++;
    }
    printf("%-28s %8s %8s %8s %8s %8s %10s\n", "mixture (REFPROP judging itself)", "<1e-8", "<1e-6", "<1e-5", "<1e-3", ">=1e-3", "now-1phase");
    for (auto& kv : bins)
        printf("%-28s %8d %8d %8d %8d %8d %10d\n", kv.first.c_str(), kv.second[0], kv.second[1], kv.second[2], kv.second[3], kv.second[4],
               kv.second[5]);
}
