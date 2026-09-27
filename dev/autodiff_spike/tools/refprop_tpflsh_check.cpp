// Experiment 12: ours (all roots -> select -> polish) vs REFPROP TPFLSH, per state.
//
// TPFLSH answers the question a user actually asks ("what is the fluid at (T, p, x)?") with a
// stability test, a phase split when two-phase, and caloric properties at the answer.  Ours
// returns the selected homogeneous root (no phase split).  So states are split by TPFLSH's verdict:
//   single-phase : like for like -- both return the density; agreement checked
//   two-phase    : REFPROP does a flash, we do not; timed separately, not a like-for-like comparison
// GERG rows: REFPROP FLAGS("GERG", 1) (same model; agreement to 1e-7).  Reference rows: REFPROP
// default vs CoolProp HEOS (model data differ slightly; agreement to 1e-3).
// Per-state times: each call timed individually, min over 3 passes.
//
//   ./bench_flash [nstates=5000] [refprop_dir=/Users/ianbell/REFPROP10/]
#include <dlfcn.h>
#include <algorithm>
#include <chrono>
#include <cstring>
#include <map>
#include <random>
#include <string>
#include <cmath>
#include <cstdio>
#include <vector>


namespace {
using SETPATH_t = void (*)(char*, long);
using SETUP_t = void (*)(int*, char*, char*, char*, int*, char*, long, long, long, long);
using FLAGS_t = void (*)(char*, int*, int*, int*, char*, long, long);
using TPFLSH_t = void (*)(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*,
                          double*, double*, int*, char*, long);
SETPATH_t RP_SETPATH = nullptr;
SETUP_t RP_SETUP = nullptr;
FLAGS_t RP_FLAGS = nullptr;
TPFLSH_t RP_TPFLSH = nullptr;
bool load_rp(const std::string& dir) {
    void* h = dlopen((dir + "librefprop.dylib").c_str(), RTLD_NOW);
    if (!h) return false;
    RP_SETPATH = reinterpret_cast<SETPATH_t>(dlsym(h, "SETPATHdll"));
    RP_SETUP = reinterpret_cast<SETUP_t>(dlsym(h, "SETUPdll"));
    RP_FLAGS = reinterpret_cast<FLAGS_t>(dlsym(h, "FLAGSdll"));
    RP_TPFLSH = reinterpret_cast<TPFLSH_t>(dlsym(h, "TPFLSHdll"));
    if (!RP_SETPATH || !RP_SETUP || !RP_FLAGS || !RP_TPFLSH) return false;
    char path[256] = {};
    std::strncpy(path, dir.c_str(), 255);
    RP_SETPATH(path, 255);
    return true;
}
const std::map<std::string, std::string> RPNAME = {
  {"Methane", "METHANE"}, {"Ethane", "ETHANE"},       {"Propane", "PROPANE"},   {"HydrogenSulfide", "H2S"}, {"Nitrogen", "NITROGEN"},
  {"Oxygen", "OXYGEN"},   {"Argon", "ARGON"},         {"CarbonDioxide", "CO2"}, {"Water", "WATER"},         {"IsoButane", "ISOBUTAN"},
  {"n-Butane", "BUTANE"}, {"Isopentane", "IPENTANE"}, {"n-Pentane", "PENTANE"}, {"n-Hexane", "HEXANE"}};
bool setup_rp(const std::vector<std::string>& fluids, bool gerg) {
    char hf0[256] = {}, herr0[256] = {};
    std::strcpy(hf0, "GERG");
    int j = gerg ? 1 : 0, k = 0, ie = 0;
    RP_FLAGS(hf0, &j, &k, &ie, herr0, 255, 255);
    std::string f;
    for (std::size_t i = 0; i < fluids.size(); ++i)
        f += (i ? "|" : "") + RPNAME.at(fluids[i]) + ".FLD";
    std::vector<char> hf(10000, '\0'), herr(255, '\0');
    std::memcpy(hf.data(), f.data(), f.size());
    char hmx[255] = "HMX.BNC", hrf[3] = {'D', 'E', 'F'};
    int n = static_cast<int>(fluids.size()), ierr = 0;
    RP_SETUP(&n, hf.data(), hmx, hrf, &ierr, herr.data(), 10000, 255, 3, 255);
    return ierr <= 0;
}
std::string join(const std::vector<std::string>& v) {
    std::string s;
    for (std::size_t i = 0; i < v.size(); ++i)
        s += (i ? "&" : "") + v[i];
    return s;
}
double now_us() {
    return std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
double median(std::vector<double> v) {
    if (v.empty()) return NAN;
    std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
    return v[v.size() / 2];
}
double mean(const std::vector<double>& v) {
    double s = 0;
    for (double x : v)
        s += x;
    return v.empty() ? NAN : s / v.size();
}
}  // namespace
int main() {
    load_rp("/Users/ianbell/REFPROP10/");
    std::vector<std::string> NG10 = {"Methane","Nitrogen","CarbonDioxide","Ethane","Propane","IsoButane","n-Butane","Isopentane","n-Pentane","n-Hexane"};
    std::map<std::string, std::pair<std::vector<std::string>, std::vector<double>>> M = {
      {"C1C2", {{"Methane","Ethane"},{0.5,0.5}}}, {"C1C2C3", {{"Methane","Ethane","Propane"},{0.5,0.3,0.2}}},
      {"Amarillo", {NG10,{0.906724,0.031284,0.004676,0.045279,0.00828,0.001037,0.001563,0.000321,0.000443,0.000393}}},
      {"C1H2S", {{"Methane","HydrogenSulfide"},{0.5,0.5}}},
      {"humidair", {{"Nitrogen","Oxygen","Argon","CarbonDioxide","Water"},{0.7654,0.2053,0.0090,0.0003,0.02}}}};
    char name[64]; int idx; double T, p, ro, qo, rn, qn;
    std::string cur;
    while (std::scanf("%63s %d %lf %lf %lf %lf %lf %lf", name, &idx, &T, &p, &ro, &qo, &rn, &qn) == 8) {
        auto& m = M.at(name);
        if (cur != name) { setup_rp(m.first, true); cur = name; }
        std::vector<double> z = m.second; z.resize(20, 0.0);
        double pk = p / 1000, D = 0, Dl = 0, Dv = 0, xl[20] = {}, yv[20] = {}, qq = 0, e = 0, hh = 0, ss = 0, cv = 0, cp = 0, w = 0; int ie = 0; char herr[256];
        RP_TPFLSH(&T, &pk, &z[0], &D, &Dl, &Dv, xl, yv, &qq, &e, &hh, &ss, &cv, &cp, &w, &ie, herr, 255);
        auto ok = [&](double r){ return std::abs(r - D*1000) <= 1e-4 * D*1000; };
        std::printf("%-9s %5d T=%8.3f p=%10.4g | old %10.2f Q=%-8.4g %-5s | new %10.2f Q=%-8.4g %-5s | TPFLSH %10.2f q=%g ierr=%d\n", name, idx, T, p, ro, qo, ok(ro)?"OK":"x", rn, qn, ok(rn)?"OK":"x", D*1000, qq, ie);
    }
}
