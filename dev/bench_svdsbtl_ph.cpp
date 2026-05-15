// State-point benchmark for the SVDSBTL backend in (h, p) coordinates.
//
// Walks a (T, p) grid via HEOS to get reference values, then queries
// SVDSBTL at (h_truth, p) and records: rel-err of rho / T / s / u,
// and the average per-call wall time for SVDSBTL update+rhomass and
// HEOS update+rhomass at that point.
//
// Output: a CSV at /tmp/bench_svdsbtl_ph_<fluid>.csv with one row per
// single-phase grid cell.
//
// Build:
//   cmake -B build -DCOOLPROP_BUILD_SVDSBTL_BENCH=ON
//   cmake --build build --target bench_svdsbtl_ph -j
//
// Run:
//   ./build/bench_svdsbtl_ph                  # defaults to Water
//   ./build/bench_svdsbtl_ph Water Propane    # subset
//
// Env:
//   BENCH_NT (default 80) — temperature grid points
//   BENCH_NP (default 60) — pressure grid points
//   BENCH_REPEATS (default 200) — calls per point for timing
//
// NOLINTBEGIN(cppcoreguidelines-pro-type-vararg,cert-err33-c,hicpp-vararg,bugprone-empty-catch)

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "AbstractState.h"
#include "Configuration.h"
#include "DataStructures.h"

#include <IF97.h>

// Direct REFPROP call path.  REFPROPMixtureBackend.cpp already
// dlopen()s the REFPROP library and SETUPdll's the chosen fluid when
// AbstractState::factory("REFPROP", ...) succeeds — the Fortran
// COMMON-block state that holds is then visible to any other call
// site against the same library handle.  All we need here is the
// PHFLSHdll symbol itself, fetched via dlsym so we can call it
// without going through CoolProp's update() dispatch and
// CachedElement plumbing.  Including REFPROP_lib.h with
// REFPROP_IMPLEMENTATION here would collide with the already-defined
// symbols in REFPROPMixtureBackend's TU, and wrapping in an
// anonymous namespace breaks the standard-library template lookups
// the header transitively needs.  So we do this surgically:
#if defined(__unix__) || defined(__APPLE__)
#    include <dlfcn.h>
#endif

extern "C"
{
    // Fortran-mangled C signature of PHFLSHdll, matching REFPROP_lib.h's
    // PHFLSHdll_ARGS expansion with REFPROP_CSTYLE_REFERENCES enabled.
    using PHFLSHdll_t = void (*)(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*,
                                 double*, double*, int*, char*, std::size_t);

    // Same convention, PHFL1dll signature: single-phase fast path that
    // expects the caller to have already classified the phase (kph =
    // 1 for liquid-like, 2 for vapor-like).  No phase iteration, no
    // sat-property output.
    //   args: p, h, x, kph, T, D, ierr, herr, herr_len
    using PHFL1dll_t = void (*)(double*, double*, double*, int*, double*, double*, int*, char*, std::size_t);
}

namespace {

std::size_t env_size_t(const char* name, std::size_t fallback) {
    const char* v = std::getenv(name);
    if (v == nullptr || v[0] == 0) {
        return fallback;
    }
    try {
        return static_cast<std::size_t>(std::stoul(v));
    } catch (...) {
        return fallback;
    }
}

struct Row
{
    double h;
    double p;
    double T;
    double rho_truth;
    double rho_pred;
    double T_pred;
    double s_truth;
    double s_pred;
    double u_truth;
    double u_pred;
    double w_truth;
    double w_pred;
    double rel_err_rho;
    double rel_err_T;
    double rel_err_s;
    double rel_err_u;
    double rel_err_w;
    double ns_per_call_svd;
    double ns_per_call_heos;
    double ns_per_call_refprop;         // via AbstractState wrapper
    double ns_per_call_refprop_direct;  // direct PHFLSHdll call (phase-determining flash)
    double ns_per_call_refprop_phfl1;   // direct PHFL1dll call (single-phase fast path, kph pre-supplied)
    double ns_per_call_if97;            // via AbstractState wrapper (only when fluid is water)
    double ns_per_call_if97_direct;     // direct IF97::rhomass_phmass (no wrapper)
    double rho_refprop;
    double rel_err_rho_refprop;
};

bool has_refprop_at_(const std::string& path) {
    return !path.empty();
}

// Run `repeats` (update + rhomass) cycles and return the average
// wall-clock per call in nanoseconds.  Uses volatile aggregation so
// the compiler can't optimise the work away.
double time_update_rhomass(::CoolProp::AbstractState& AS, ::CoolProp::input_pairs pair, double v1, double v2, std::size_t repeats) {
    volatile double sink = 0.0;
    const auto t0 = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < repeats; ++i) {
        AS.update(pair, v1, v2);
        sink += AS.rhomass();
    }
    const auto t1 = std::chrono::steady_clock::now();
    (void)sink;
    return std::chrono::duration<double, std::nano>(t1 - t0).count() / static_cast<double>(repeats);
}

// dlopen the REFPROP library (refcounted; returns existing handle if
// already loaded) and dlsym a named symbol, trying the canonical
// Fortran-mangling variants in order.  Returns nullptr if either step
// fails on this platform.
void* refprop_sym_(const std::string& refprop_path, const std::vector<const char*>& names) {
#if defined(__unix__) || defined(__APPLE__)
    const std::string lib_path = refprop_path + "librefprop.dylib";
    void* handle = dlopen(lib_path.c_str(), RTLD_NOW);
    if (handle == nullptr) {
        handle = dlopen((refprop_path + "librefprop.so").c_str(), RTLD_NOW);
    }
    if (handle == nullptr) {
        return nullptr;
    }
    for (const char* sym : names) {
        void* p = dlsym(handle, sym);
        if (p != nullptr) {
            return p;
        }
    }
    return nullptr;
#else
    (void)refprop_path;
    (void)names;
    return nullptr;
#endif
}

PHFLSHdll_t resolve_phflshdll_(const std::string& refprop_path) {
    return reinterpret_cast<PHFLSHdll_t>(refprop_sym_(refprop_path, {"PHFLSHdll", "phflshdll", "phflshdll_"}));
}

PHFL1dll_t resolve_phfl1dll_(const std::string& refprop_path) {
    return reinterpret_cast<PHFL1dll_t>(refprop_sym_(refprop_path, {"PHFL1dll", "phfl1dll", "phfl1dll_"}));
}

// Direct IAPWS-IF97 timing: pure header-only inline call, no
// AbstractState dispatch.  This is the closest we can get to the
// claim in the IF97 release CTR -- the only overhead beyond the math
// is the call-site / loop / volatile sink.
double time_if97_rhomass_phmass_direct(double p_Pa, double h_J_per_kg, std::size_t repeats) {
    volatile double sink = 0.0;
    const auto t0 = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < repeats; ++i) {
        sink += IF97::rhomass_phmass(p_Pa, h_J_per_kg);
    }
    const auto t1 = std::chrono::steady_clock::now();
    (void)sink;
    return std::chrono::duration<double, std::nano>(t1 - t0).count() / static_cast<double>(repeats);
}

// Direct REFPROP PHFLSHdll timing.  Calls the Fortran flash routine
// straight through the dlsym'd function pointer — no AbstractState
// dispatch, no CachedElement.clear(), no kg <-> mol round-trips that
// the wrapper does on every call.  Returns ns per (p, h_mol) flash.
//
// The caller must have already SETUPdll'd the right fluid (which our
// CoolProp factory("REFPROP", fluid) call has done by the time we
// reach this).
// PHFL1dll single-phase fast path.  Caller supplies kph (1=liquid,
// 2=vapor) -- the routine skips REFPROP's internal phase determination.
double time_phfl1_direct(PHFL1dll_t phfl1, double p_kPa, double h_mol, double* mole_fractions, int kph, std::size_t repeats) {
    double T = 0;
    double d = 0;
    int ierr = 0;
    std::array<char, 256> herr{};
    volatile double sink = 0.0;
    const auto t0 = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < repeats; ++i) {
        phfl1(&p_kPa, &h_mol, mole_fractions, &kph, &T, &d, &ierr, herr.data(), 255);
        sink += d;
    }
    const auto t1 = std::chrono::steady_clock::now();
    (void)sink;
    return std::chrono::duration<double, std::nano>(t1 - t0).count() / static_cast<double>(repeats);
}

double time_phflsh_direct(PHFLSHdll_t phflsh, double p_kPa, double h_mol, double* mole_fractions, std::size_t repeats) {
    double T = 0;
    double d = 0;
    double dl = 0;
    double dv = 0;
    std::array<double, 20> xliq{};
    std::array<double, 20> xvap{};
    double q = 0;
    double e = 0;
    double s = 0;
    double cv = 0;
    double cp = 0;
    double w = 0;
    int ierr = 0;
    std::array<char, 256> herr{};
    volatile double sink = 0.0;
    const auto t0 = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < repeats; ++i) {
        phflsh(&p_kPa, &h_mol, mole_fractions, &T, &d, &dl, &dv, xliq.data(), xvap.data(), &q, &e, &s, &cv, &cp, &w, &ierr, herr.data(), 255);
        sink += d;
    }
    const auto t1 = std::chrono::steady_clock::now();
    (void)sink;
    return std::chrono::duration<double, std::nano>(t1 - t0).count() / static_cast<double>(repeats);
}

void bench_one(const std::string& fluid, std::size_t NT, std::size_t NP, std::size_t repeats, const std::string& refprop_path) {
    auto svd = std::shared_ptr<::CoolProp::AbstractState>(::CoolProp::AbstractState::factory("SVDSBTL", fluid));
    auto heos = std::shared_ptr<::CoolProp::AbstractState>(::CoolProp::AbstractState::factory("HEOS", fluid));

    // IF97 only applies to water; lazily allocate so other fluids skip it.
    std::shared_ptr<::CoolProp::AbstractState> if97;
    try {
        if97.reset(::CoolProp::AbstractState::factory("IF97", fluid));
    } catch (const std::exception&) {  // NOLINT(bugprone-empty-catch)
        if97.reset();
    }

    std::shared_ptr<::CoolProp::AbstractState> refprop;
    double refprop_molar_mass = 0.0;
    double refprop_rho_crit_mol_L = 0.0;
    std::array<double, 20> refprop_mole_fractions{};
    PHFLSHdll_t phflsh_direct = nullptr;
    PHFL1dll_t phfl1_direct = nullptr;
    if (has_refprop_at_(refprop_path)) {
        ::CoolProp::set_config_string(ALTERNATIVE_REFPROP_PATH, refprop_path);
        try {
            refprop.reset(::CoolProp::AbstractState::factory("REFPROP", fluid));
            refprop_molar_mass = refprop->molar_mass();                                        // kg/mol
            refprop_rho_crit_mol_L = refprop->rhomass_critical() / refprop_molar_mass * 1e-3;  // mol/L
            refprop_mole_fractions[0] = 1.0;
            phflsh_direct = resolve_phflshdll_(refprop_path);
            phfl1_direct = resolve_phfl1dll_(refprop_path);
            if (phflsh_direct == nullptr) {
                std::printf("  REFPROP direct PHFLSHdll resolve failed\n");
            }
            if (phfl1_direct == nullptr) {
                std::printf("  REFPROP direct PHFL1dll resolve failed\n");
            }
        } catch (const std::exception& e) {
            std::printf("  REFPROP unavailable: %s\n", e.what());
            refprop.reset();
        }
    }

    // Full phase-diagram coverage: triple-point pressure -> backend's
    // upper limit, triple-point temperature -> Tmax.  Goes well past
    // p_crit into supercritical territory; SVDSBTL only covers
    // subcritical so above-p_crit cells will record rho_pred = NaN.
    // We keep them in the output anyway so the plot can show where
    // coverage stops.
    const double T_min = std::max(heos->Ttriple(), heos->Tmin()) + 0.5;
    const double T_max = heos->Tmax() - 0.5;
    const double p_min = heos->p_triple() * 1.1;
    const double p_max = heos->pmax() * 0.99;

    std::vector<Row> rows;
    rows.reserve(NT * NP);

    for (std::size_t i = 0; i < NT; ++i) {
        const double T = T_min + (T_max - T_min) * static_cast<double>(i) / static_cast<double>(NT - 1);
        for (std::size_t j = 0; j < NP; ++j) {
            const double log_p = std::log(p_min) + (std::log(p_max) - std::log(p_min)) * static_cast<double>(j) / static_cast<double>(NP - 1);
            const double p = std::exp(log_p);
            try {
                heos->update(::CoolProp::PT_INPUTS, p, T);
                if (heos->Q() > 0.0 && heos->Q() < 1.0) {
                    continue;  // two-phase, not covered by SVDSBTL subcritical regions
                }
                Row r{};
                r.h = heos->hmass();
                r.p = p;
                r.T = T;
                r.rho_truth = heos->rhomass();
                r.s_truth = heos->smass();
                r.u_truth = heos->umass();
                r.w_truth = heos->speed_sound();

                svd->update(::CoolProp::HmassP_INPUTS, r.h, r.p);
                r.rho_pred = svd->rhomass();
                r.T_pred = svd->T();
                r.s_pred = svd->smass();
                r.u_pred = svd->umass();
                // speed_sound throws in two-phase; the bench skips
                // two-phase cells via the Q check above, so the
                // single-phase eval is always safe here.
                r.w_pred = svd->speed_sound();

                const bool svd_in_domain = !std::isnan(r.rho_pred);
                if (svd_in_domain) {
                    r.rel_err_rho = std::abs(r.rho_pred - r.rho_truth) / r.rho_truth;
                    r.rel_err_T = std::abs(r.T_pred - r.T) / r.T;
                    r.rel_err_s = (r.s_truth != 0.0) ? std::abs(r.s_pred - r.s_truth) / std::abs(r.s_truth) : 0.0;
                    r.rel_err_u = (r.u_truth != 0.0) ? std::abs(r.u_pred - r.u_truth) / std::abs(r.u_truth) : 0.0;
                    r.rel_err_w = (r.w_truth != 0.0) ? std::abs(r.w_pred - r.w_truth) / std::abs(r.w_truth) : 0.0;
                    r.ns_per_call_svd = time_update_rhomass(*svd, ::CoolProp::HmassP_INPUTS, r.h, r.p, repeats);
                } else {
                    r.rel_err_rho = std::nan("");
                    r.rel_err_T = std::nan("");
                    r.rel_err_s = std::nan("");
                    r.rel_err_u = std::nan("");
                    r.rel_err_w = std::nan("");
                    r.ns_per_call_svd = std::nan("");
                }
                r.ns_per_call_heos = time_update_rhomass(*heos, ::CoolProp::PT_INPUTS, p, T, repeats);

                if (if97) {
                    try {
                        r.ns_per_call_if97 = time_update_rhomass(*if97, ::CoolProp::HmassP_INPUTS, r.h, r.p, repeats);
                        r.ns_per_call_if97_direct = time_if97_rhomass_phmass_direct(r.p, r.h, repeats);
                    } catch (...) {  // NOLINT(bugprone-empty-catch)
                        r.ns_per_call_if97 = std::nan("");
                        r.ns_per_call_if97_direct = std::nan("");
                    }
                } else {
                    r.ns_per_call_if97 = std::nan("");
                    r.ns_per_call_if97_direct = std::nan("");
                }

                if (refprop) {
                    try {
                        // Use HmassP for a fair shape comparison with SVDSBTL.
                        refprop->update(::CoolProp::HmassP_INPUTS, r.h, p);
                        r.rho_refprop = refprop->rhomass();
                        r.rel_err_rho_refprop = std::abs(r.rho_pred - r.rho_refprop) / r.rho_refprop;
                        r.ns_per_call_refprop = time_update_rhomass(*refprop, ::CoolProp::HmassP_INPUTS, r.h, r.p, repeats);
                        const double p_kPa = r.p * 1e-3;
                        const double h_mol = r.h * refprop_molar_mass;  // J/kg * kg/mol -> J/mol
                        if (phflsh_direct != nullptr) {
                            r.ns_per_call_refprop_direct = time_phflsh_direct(phflsh_direct, p_kPa, h_mol, refprop_mole_fractions.data(), repeats);
                        } else {
                            r.ns_per_call_refprop_direct = std::nan("");
                        }
                        if (phfl1_direct != nullptr) {
                            // kph = 1 (liquid-like) when denser than rho_crit, 2 (vapor-like) otherwise.
                            // Use REFPROP's own rho output -- safe single-phase signal we just got.
                            const double rho_mol_L = r.rho_refprop / refprop_molar_mass * 1e-3;
                            const int kph = (rho_mol_L > refprop_rho_crit_mol_L) ? 1 : 2;
                            r.ns_per_call_refprop_phfl1 = time_phfl1_direct(phfl1_direct, p_kPa, h_mol, refprop_mole_fractions.data(), kph, repeats);
                        } else {
                            r.ns_per_call_refprop_phfl1 = std::nan("");
                        }
                    } catch (...) {  // NOLINT(bugprone-empty-catch)
                        r.rho_refprop = std::nan("");
                        r.rel_err_rho_refprop = std::nan("");
                        r.ns_per_call_refprop = std::nan("");
                        r.ns_per_call_refprop_direct = std::nan("");
                        r.ns_per_call_refprop_phfl1 = std::nan("");
                    }
                } else {
                    r.rho_refprop = std::nan("");
                    r.rel_err_rho_refprop = std::nan("");
                    r.ns_per_call_refprop = std::nan("");
                    r.ns_per_call_refprop_direct = std::nan("");
                    r.ns_per_call_refprop_phfl1 = std::nan("");
                }

                rows.push_back(r);
            } catch (...) {  // NOLINT(bugprone-empty-catch)
            }
        }
    }

    const std::string out_path = std::string("/tmp/bench_svdsbtl_ph_") + fluid + ".csv";
    std::ofstream ofs(out_path);
    if (!ofs) {
        throw std::runtime_error("bench_one: failed to open " + out_path);
    }
    ofs << "h,p,T,rho_truth,rho_pred,T_pred,s_truth,s_pred,u_truth,u_pred,w_truth,w_pred,"
        << "rel_err_rho,rel_err_T,rel_err_s,rel_err_u,rel_err_w,"
        << "ns_per_call_svd,ns_per_call_heos,ns_per_call_refprop,ns_per_call_refprop_direct,ns_per_call_refprop_phfl1,"
        << "ns_per_call_if97,ns_per_call_if97_direct,"
        << "rho_refprop,rel_err_rho_refprop\n";
    ofs.precision(17);
    for (const auto& r : rows) {
        ofs << r.h << "," << r.p << "," << r.T << "," << r.rho_truth << "," << r.rho_pred << "," << r.T_pred << "," << r.s_truth << "," << r.s_pred
            << "," << r.u_truth << "," << r.u_pred << "," << r.w_truth << "," << r.w_pred << "," << r.rel_err_rho << "," << r.rel_err_T << ","
            << r.rel_err_s << "," << r.rel_err_u << "," << r.rel_err_w << "," << r.ns_per_call_svd << "," << r.ns_per_call_heos << ","
            << r.ns_per_call_refprop << "," << r.ns_per_call_refprop_direct << "," << r.ns_per_call_refprop_phfl1 << "," << r.ns_per_call_if97 << ","
            << r.ns_per_call_if97_direct << "," << r.rho_refprop << "," << r.rel_err_rho_refprop << "\n";
    }

    // Console summary.
    double max_rho = 0;
    double max_w = 0;
    double max_T = 0;
    double max_s = 0;
    double max_u = 0;
    double sum_ns_svd = 0;
    double sum_ns_heos = 0;
    double sum_ns_refprop = 0;
    double sum_ns_refprop_direct = 0;
    double sum_ns_refprop_phfl1 = 0;
    double sum_ns_if97 = 0;
    double sum_ns_if97_direct = 0;
    std::size_t n_svd = 0;  // cells where SVDSBTL is in-domain
    std::size_t n_refprop = 0;
    std::size_t n_refprop_direct = 0;
    std::size_t n_refprop_phfl1 = 0;
    std::size_t n_if97 = 0;
    std::size_t n_if97_direct = 0;
    for (const auto& r : rows) {
        if (!std::isnan(r.rel_err_rho)) {
            max_rho = std::max(max_rho, r.rel_err_rho);
            if (!std::isnan(r.rel_err_w)) {
                max_w = std::max(max_w, r.rel_err_w);
            }
            max_T = std::max(max_T, r.rel_err_T);
            max_s = std::max(max_s, r.rel_err_s);
            max_u = std::max(max_u, r.rel_err_u);
            sum_ns_svd += r.ns_per_call_svd;
            ++n_svd;
        }
        sum_ns_heos += r.ns_per_call_heos;
        if (!std::isnan(r.ns_per_call_refprop)) {
            sum_ns_refprop += r.ns_per_call_refprop;
            ++n_refprop;
        }
        if (!std::isnan(r.ns_per_call_refprop_direct)) {
            sum_ns_refprop_direct += r.ns_per_call_refprop_direct;
            ++n_refprop_direct;
        }
        if (!std::isnan(r.ns_per_call_refprop_phfl1)) {
            sum_ns_refprop_phfl1 += r.ns_per_call_refprop_phfl1;
            ++n_refprop_phfl1;
        }
        if (!std::isnan(r.ns_per_call_if97)) {
            sum_ns_if97 += r.ns_per_call_if97;
            ++n_if97;
        }
        if (!std::isnan(r.ns_per_call_if97_direct)) {
            sum_ns_if97_direct += r.ns_per_call_if97_direct;
            ++n_if97_direct;
        }
    }
    const auto n = static_cast<double>(rows.size());
    const auto n_in_domain = static_cast<double>(n_svd);
    const double mean_svd = (n_svd > 0) ? sum_ns_svd / n_in_domain : std::nan("");
    const double mean_heos = sum_ns_heos / n;
    std::printf("%-12s N=%5zu (SVDSBTL in-domain=%zu, out-of-domain=%zu)\n", fluid.c_str(), rows.size(), n_svd, rows.size() - n_svd);
    std::printf("            max_rel: rho=%.2e T=%.2e s=%.2e u=%.2e w=%.2e\n", max_rho, max_T, max_s, max_u, max_w);
    std::printf("            mean ns/call: svd=%.0f  heos=%.0f", mean_svd, mean_heos);
    if (n_refprop > 0) {
        std::printf("  refprop(via AS)=%.0f", sum_ns_refprop / static_cast<double>(n_refprop));
    }
    if (n_refprop_direct > 0) {
        std::printf("  refprop(direct PHFLSHdll)=%.0f", sum_ns_refprop_direct / static_cast<double>(n_refprop_direct));
    }
    if (n_refprop_phfl1 > 0) {
        std::printf("  refprop(direct PHFL1dll)=%.0f", sum_ns_refprop_phfl1 / static_cast<double>(n_refprop_phfl1));
    }
    if (n_if97 > 0) {
        std::printf("  if97(via AS)=%.0f", sum_ns_if97 / static_cast<double>(n_if97));
    }
    if (n_if97_direct > 0) {
        std::printf("  if97(direct rhomass_phmass)=%.0f", sum_ns_if97_direct / static_cast<double>(n_if97_direct));
    }
    std::printf("\n            speedup vs heos=%.1fx", mean_heos / mean_svd);
    if (n_refprop > 0) {
        std::printf("  vs refprop(AS)=%.1fx", (sum_ns_refprop / static_cast<double>(n_refprop)) / mean_svd);
    }
    if (n_refprop_direct > 0) {
        std::printf("  vs refprop(direct)=%.1fx", (sum_ns_refprop_direct / static_cast<double>(n_refprop_direct)) / mean_svd);
    }
    if (n_refprop_phfl1 > 0) {
        std::printf("  vs refprop(PHFL1)=%.1fx", (sum_ns_refprop_phfl1 / static_cast<double>(n_refprop_phfl1)) / mean_svd);
    }
    if (n_if97 > 0) {
        std::printf("  vs if97(AS)=%.2fx", (sum_ns_if97 / static_cast<double>(n_if97)) / mean_svd);
    }
    if (n_if97_direct > 0) {
        std::printf("  vs if97(direct)=%.2fx", (sum_ns_if97_direct / static_cast<double>(n_if97_direct)) / mean_svd);
    }
    std::printf("\n            -> %s\n", out_path.c_str());
}

}  // namespace

int main(int argc, char** argv) {
    const std::size_t NT = env_size_t("BENCH_NT", 200);
    const std::size_t NP = env_size_t("BENCH_NP", 150);
    const std::size_t repeats = env_size_t("BENCH_REPEATS", 100);
    const char* refprop_env = std::getenv("BENCH_REFPROP_PATH");
    const std::string refprop_path = (refprop_env != nullptr) ? std::string(refprop_env) : std::string();

    std::vector<std::string> fluids;
    if (argc > 1) {
        for (int i = 1; i < argc; ++i) {
            fluids.emplace_back(argv[i]);
        }
    } else {
        fluids = {"Water"};
    }

    std::printf("SVDSBTL PH benchmark (NT=%zu  NP=%zu  repeats=%zu  refprop=%s)\n\n", NT, NP, repeats,
                refprop_path.empty() ? "off" : refprop_path.c_str());
    int failed = 0;
    for (const auto& fluid : fluids) {
        try {
            bench_one(fluid, NT, NP, repeats, refprop_path);
        } catch (const std::exception& e) {
            std::printf("%-12s  ERROR: %s\n", fluid.c_str(), e.what());
            ++failed;
        }
    }
    return failed == 0 ? 0 : 1;
}

// NOLINTEND(cppcoreguidelines-pro-type-vararg,cert-err33-c,hicpp-vararg,bugprone-empty-catch)
