// Per-state verdicts for CoolProp/TPFLSH disagreements (see main).  Build like the other tools, with -Isrc -Idev/autodiff_spike.
#include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#include "Backends/Helmholtz/MixtureDerivatives.h"
#include "gerg_cheb_solver.hpp"
#include <cmath>
#include <fstream>
#include <dlfcn.h>
#include <algorithm>
#include <chrono>
#include <cstring>
#include <map>
#include <random>
#include <string>
#include "AbstractState.h"

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
  {"n-Butane", "BUTANE"}, {"Isopentane", "IPENTANE"}, {"n-Pentane", "PENTANE"}, {"n-Hexane", "HEXANE"},
  {"R32", "R32"},         {"R1234yf", "R1234YF"}, {"Helium", "HELIUM"},       {"Hydrogen", "HYDROGEN"}};
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

double gRT(CoolProp::HelmholtzEOSMixtureBackend& H) {
    double g = 0; auto x = H.get_mole_fractions();
    for (std::size_t i = 0; i < x.size(); ++i) if (x[i] > 0) g += x[i] * (std::log(x[i]) + CoolProp::MixtureDerivatives::ln_fugacity_coefficient(H, i, CoolProp::XN_INDEPENDENT));
    return g;
}
std::vector<std::string> split_csv(const std::string& l) {
    std::vector<std::string> out; std::string cur; bool q = false;
    for (char c : l) { if (c == '"') q = !q; else if (c == ',' && !q) { out.push_back(cur); cur.clear(); } else cur += c; }
    out.push_back(cur); return out;
}
}  // namespace

// Per-state verdict for every CoolProp/TPFLSH disagreement in a bench_ratio CSV.
//   ./verdict in.csv out_verdicts.csv
// Verdicts: cp_wrong, rp_wrong (non-equilibrium / non-reproducible TPFLSH split, or wrong single-phase root),
// rp_scope (TPFLSH misses a validated CoolProp split, e.g. LLE), unresolved.  Agreeing states are not listed.
int main(int argc, char** argv) {
    if (argc < 3) { std::printf("usage: verdict in.csv out.csv\n"); return 1; }
    load_rp("/Users/ianbell/REFPROP10/");
    struct Mix { std::string backend; std::vector<std::string> f; std::vector<double> z; };
    const std::vector<std::string> NG10 = {"Methane","Nitrogen","CarbonDioxide","Ethane","Propane","IsoButane","n-Butane","Isopentane","n-Pentane","n-Hexane"};
    std::map<std::string, Mix> M = {
      {"C1/C2 50/50", {"GERG2008", {"Methane","Ethane"}, {0.5,0.5}}},
      {"C1/C2/C3 50/30/20", {"GERG2008", {"Methane","Ethane","Propane"}, {0.5,0.3,0.2}}},
      {"Amarillo (10)", {"GERG2008", NG10, {0.906724,0.031284,0.004676,0.045279,0.00828,0.001037,0.001563,0.000321,0.000443,0.000393}}},
      {"C1/H2S 50/50", {"GERG2008", {"Methane","HydrogenSulfide"}, {0.5,0.5}}},
      {"humid air x_w=0.02", {"GERG2008", {"Nitrogen","Oxygen","Argon","CarbonDioxide","Water"}, {0.7654,0.2053,0.0090,0.0003,0.02}}},
      {"R454B", {"HEOS", {"R32","R1234yf"}, {0.8292,0.1708}}},
      {"N2/C1/C2/nC4/nC5", {"GERG2008", {"Nitrogen","Methane","Ethane","n-Butane","n-Pentane"}, {0.3797,0.3225,0.278,0.0014,0.0184}}},
      {"CO2/H2 96.77/3.23", {"GERG2008", {"CarbonDioxide","Hydrogen"}, {0.9677,0.0323}}},
      {"CO2/N2/O2/He", {"GERG2008", {"CarbonDioxide","Nitrogen","Oxygen","Helium"}, {0.9403/1.00002,0.0582/1.00002,0.00127/1.00002,0.00025/1.00002}}},
      {"N2/O2/Ar (O2-enriched air)", {"GERG2008", {"Nitrogen","Oxygen","Argon"}, {0.609067,0.370414,0.0205193}}}};
    std::ifstream in(argv[1]); std::ofstream out(argv[2]);
    out << "mixture,i,verdict,detail\n";
    std::string line; std::getline(in, line);
    std::string cur; std::unique_ptr<CoolProp::AbstractState> AS, K, L, V; gergcheb::Solver sv; bool gerg = true; Mix mx;
    std::map<std::string, int> cnt;
    while (std::getline(in, line)) {
        auto c = split_csv(line);
        const std::string name = c[0]; const int idx = std::stoi(c[1]);
        const double T = std::stod(c[2]), p = std::stod(c[3]), rcp = std::stod(c[6]), Qcp = std::stod(c[7]), rrp = std::stod(c[9]), qrp = std::stod(c[10]);
        const int cpf = std::stoi(c[8]), ie = std::stoi(c[11]);
        if (cpf || ie > 0) continue;  // hard failures are reported separately
        if (name != cur) {
            mx = M.at(name); cur = name; gerg = mx.backend == "GERG2008";
            AS.reset(CoolProp::AbstractState::factory(mx.backend, join(mx.f))); K.reset(CoolProp::AbstractState::factory(mx.backend, join(mx.f)));
            L.reset(CoolProp::AbstractState::factory(mx.backend, join(mx.f))); V.reset(CoolProp::AbstractState::factory(mx.backend, join(mx.f)));
            K->set_mole_fractions(mx.z);
            auto* HK = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(K.get());
            sv = gergcheb::Solver(); sv.build(HK, 4.0, 1.3 * HK->Reducing->Tr(mx.z) / 60.0, 1e-6);
            setup_rp(mx.f, gerg);
        }
        const bool cp2 = Qcp > 0 && Qcp < 1, rp2 = qrp >= 0 && qrp <= 1;
        const double tol = gerg ? 1e-4 : 1e-3;
        if (cp2 == rp2 && std::abs(rcp - rrp) <= tol * rrp) continue;  // agree
        try {
        auto& HK = *dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(K.get());
        // reference single phase: kernel spinodal-branch selection
        gergcheb::Solver::State S; sv.assemble(T, mx.z, S); gergcheb::Solver::Root r[gergcheb::MAXROOTS];
        const int n = sv.roots(S, p, r, gergcheb::Solver::Polish::All); const int ks = sv.select(S, p, r, n);
        if (ks < 0) { ++cnt["unresolved"]; out << "\"" << name << "\"," << idx << ",unresolved,no kernel root\n"; continue; }
        HK.set_mole_fractions(mx.z); HK.update_DmolarT_direct(r[ks].rho, T); const double gsel = gRT(HK); const double rsel = r[ks].rho;
        std::string v, d;
        // CoolProp side (valid for every backend)
        double gcp = NAN; bool cp_ok = true;
        std::unique_ptr<CoolProp::AbstractState> FR(CoolProp::AbstractState::factory(mx.backend, join(mx.f)));
        FR->set_mole_fractions(mx.z);
        try { FR->update(CoolProp::PT_INPUTS, p, T); } catch (...) {}
        const bool fr2 = FR->phase() == CoolProp::iphase_twophase;
        if (fr2 != cp2 || (!cp2 && std::abs(FR->rhomolar() - rcp) > 1e-7 * rcp)) {
            ++cnt["cp_history"];
            out << "\"" << name << "\"," << idx << ",cp_history,\"CoolProp answer depends on call history (fresh object: " << (fr2 ? "two-phase" : "single phase") << ")\"\n";
            continue;
        }
        if (cp2) {
            auto& AS = FR;
            const auto x = AS->mole_fractions_liquid_double(), y = AS->mole_fractions_vapor_double(); const double Q = AS->Q();
            L->set_mole_fractions(x); V->set_mole_fractions(y); L->specify_phase(CoolProp::iphase_liquid); V->specify_phase(CoolProp::iphase_gas);
            L->update(CoolProp::DmolarT_INPUTS, AS->saturated_liquid_keyed_output(CoolProp::iDmolar), T);
            V->update(CoolProp::DmolarT_INPUTS, AS->saturated_vapor_keyed_output(CoolProp::iDmolar), T);
            double dl = 0; for (std::size_t i = 0; i < x.size(); ++i) dl = std::max(dl, std::abs(std::log(V->fugacity(i) / L->fugacity(i))));
            gcp = (1 - Q) * gRT(*dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(L.get())) + Q * gRT(*dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(V.get()));
            cp_ok = dl < 1e-5 && gcp < gsel;
            if (!cp_ok) { v = "cp_wrong"; d = "CoolProp split not an equilibrium or not below the single phase"; }
        } else {
            cp_ok = std::abs(rcp - rsel) <= 1e-6 * rsel;
            gcp = gsel;
            if (!cp_ok) { v = "cp_wrong"; d = "CoolProp single-phase root is not the spinodal-branch root"; }
        }
        if (v.empty() && !gerg) { v = "unresolved"; d = "different models (HEOS vs REFPROP default)"; }
        if (v.empty()) {
            // REFPROP side (same model)
            std::vector<double> z20 = mx.z; z20.resize(20, 0.0);
            double TT = T, pk = p / 1000, D = 0, Dl = 0, Dv = 0, xl[20] = {}, yv[20] = {}, q = 0, e = 0, hh = 0, ss = 0, cv = 0, cpp = 0, w = 0; int ie2 = 0; char herr[256];
            RP_TPFLSH(&TT, &pk, &z20[0], &D, &Dl, &Dv, xl, yv, &q, &e, &hh, &ss, &cv, &cpp, &w, &ie2, herr, 255);
            const bool rp2b = q >= 0 && q <= 1;
            if (rp2 && !rp2b) { v = "rp_wrong"; d = "TPFLSH two-phase answer not reproducible (second call single phase)"; }
            else if (rp2) {
                const std::size_t N = mx.z.size();
                std::vector<double> x(xl, xl + N), y(yv, yv + N);
                auto& HL = *dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(L.get()); auto& HV = *dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(V.get());
                L->set_mole_fractions(x); V->set_mole_fractions(y); L->specify_phase(CoolProp::iphase_liquid); V->specify_phase(CoolProp::iphase_gas);
                L->update(CoolProp::DmolarT_INPUTS, Dl * 1000, T); V->update(CoolProp::DmolarT_INPUTS, Dv * 1000, T);
                double dl = 0; for (std::size_t i = 0; i < N; ++i) dl = std::max(dl, std::abs(std::log(HV.fugacity(i) / HL.fugacity(i))));
                const bool peq = std::abs(HL.p() / p - 1) < 1e-5 && std::abs(HV.p() / p - 1) < 1e-5;
                const double grp = (1 - q) * gRT(HL) + q * gRT(HV);
                if (dl > 1e-5 || !peq) { v = "rp_wrong"; d = "TPFLSH split is not an equilibrium in the model"; }
                else if (!cp2 && grp < gcp - 1e-10) { v = "cp_wrong"; d = "CoolProp missed a genuine lower-G split"; }
                else if (cp2 && grp < gcp - 1e-10) { v = "cp_wrong"; d = "TPFLSH split lower in G than CoolProp's split"; }
                else if (cp2) { v = "rp_wrong"; d = "CoolProp split lower in G than TPFLSH's split"; }
                else { v = "rp_wrong"; d = "TPFLSH split not below the single phase"; }
            } else {
                if (cp2) { v = "rp_scope"; d = "TPFLSH single phase; CoolProp split validated (equal f, lower G) - LLE or missed VLE"; }
                else {
                    const bool rp_is_sel = std::abs(rrp - rsel) <= 1e-6 * rsel;
                    v = rp_is_sel ? "unresolved" : "rp_wrong"; d = rp_is_sel ? "both roots? CoolProp equals select, REFPROP too" : "TPFLSH single-phase root is not the spinodal-branch root";
                }
            }
        }
        ++cnt[v];
        out << "\"" << name << "\"," << idx << "," << v << ",\"" << d << "\"\n";
        } catch (std::exception& e) {
            ++cnt["unresolved"];
            out << "\"" << name << "\"," << idx << ",unresolved,\"check threw: " << e.what() << "\"\n";
        }
    }
    for (auto& kv : cnt) std::printf("%6d %s\n", kv.second, kv.first.c_str());
}
