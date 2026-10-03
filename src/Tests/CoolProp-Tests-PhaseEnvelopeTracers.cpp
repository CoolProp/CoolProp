// Fast tests for the experimental isopleth tracers behind PHASE_ENVELOPE_ALGORITHM.
// The torture corpus over all predefined mixtures lives in CoolProp-Tests-PhaseEnvelopeTorture.cpp.

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/Configuration.h"
#    include "CoolProp/DataStructures.h"
#    include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#    include "Backends/Helmholtz/PhaseEnvelopeTracers.h"
#    include "Backends/Helmholtz/MixtureDerivatives.h"
#    include "CoolProp/CoolProp.h"
#    include "CoolProp/detail/strings.h"
#    include "CoolProp/detail/tools.h"
#    include <cstdlib>
#    include <iostream>

#    include <cmath>
#    include <memory>
#    include <string>
#    include <vector>

using namespace CoolProp;

namespace {

/// Restores PHASE_ENVELOPE_ALGORITHM to its default when a test section ends, pass or fail.
struct AlgorithmGuard
{
    explicit AlgorithmGuard(const std::string& algorithm) {
        set_config_string(PHASE_ENVELOPE_ALGORITHM, algorithm);
    }
    ~AlgorithmGuard() {
        set_config_string(PHASE_ENVELOPE_ALGORITHM, "legacy");
    }
    AlgorithmGuard(const AlgorithmGuard&) = delete;
    AlgorithmGuard& operator=(const AlgorithmGuard&) = delete;
    AlgorithmGuard(AlgorithmGuard&&) = delete;
    AlgorithmGuard& operator=(AlgorithmGuard&&) = delete;
};

std::shared_ptr<HelmholtzEOSMixtureBackend> make_heos(const std::string& fluids, const std::vector<double>& z) {
    std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", fluids));
    AS->set_mole_fractions(z);
    auto HEOS = std::dynamic_pointer_cast<HelmholtzEOSMixtureBackend>(AS);
    REQUIRE(HEOS);
    return HEOS;
}

/// Central finite-difference check of a system Jacobian at the starting point of the trace.
void check_jacobian(const std::string& fluids, const std::vector<double>& z, const std::string& algorithm) {
    auto HEOS = make_heos(fluids, z);
    auto sys = PhaseEnvelopeTracers::make_system(*HEOS, algorithm);
    PhaseEnvelopeTracers::Options opts;
    auto s0 = PhaseEnvelopeTracers::starting_point(*HEOS, opts);
    Eigen::VectorXd X = sys->pack(s0);
    sys->set_p_ref(s0.p);
    const auto n = static_cast<Eigen::Index>(sys->size());
    const std::size_t ns = sys->index_marching();
    Eigen::VectorXd F(n);
    Eigen::MatrixXd J(n, n);
    sys->residual_jacobian(X, ns, X[static_cast<Eigen::Index>(ns)], F, J);
    REQUIRE(F.allFinite());
    REQUIRE(J.allFinite());
    for (Eigen::Index j = 0; j < n; ++j) {
        // Central differences at three step sizes; a real Jacobian error persists across all of them
        Eigen::VectorXd best_err = Eigen::VectorXd::Constant(n, 1e300), best_fd = Eigen::VectorXd::Zero(n);
        for (double h : {1e-4, 1e-5, 1e-6}) {
            Eigen::VectorXd Xp = X, Xm = X, Fp(n), Fm(n);
            Eigen::MatrixXd Jd(n, n);
            Xp[j] += h;
            Xm[j] -= h;
            sys->residual_jacobian(Xp, ns, X[static_cast<Eigen::Index>(ns)], Fp, Jd);
            sys->residual_jacobian(Xm, ns, X[static_cast<Eigen::Index>(ns)], Fm, Jd);
            for (Eigen::Index i = 0; i < n; ++i) {
                const double fd = (Fp[i] - Fm[i]) / (2 * h);
                const double err = std::abs(fd - J(i, j));
                if (err < best_err[i]) {
                    best_err[i] = err;
                    best_fd[i] = fd;
                }
            }
        }
        for (Eigen::Index i = 0; i < n; ++i) {
            CAPTURE(algorithm, i, j, best_fd[i], J(i, j));
            CHECK(best_err[i] <= 1e-5 * std::abs(J(i, j)) + 1e-7);
        }
    }
}

/// Linear estimate of the critical point from the stored ln K crossing (icrit-1 -> icrit).
void critical_estimate(const PhaseEnvelopeData& env, double& Tc, double& pc) {
    REQUIRE(env.icrit >= 1);
    REQUIRE(env.icrit < env.T.size());
    const std::size_t i1 = env.icrit, i0 = i1 - 1;
    // Use the component with the largest |delta lnK| for the interpolation fraction
    double best = -1, a = 0;
    for (const auto& lnK : env.lnK) {
        const double d = std::abs(lnK[i1] - lnK[i0]);
        if (d > best) {
            best = d;
            a = -lnK[i0] / (lnK[i1] - lnK[i0]);
        }
    }
    Tc = env.T[i0] + a * (env.T[i1] - env.T[i0]);
    pc = env.p[i0] + a * (env.p[i1] - env.p[i0]);
}

}  // namespace

TEST_CASE("Phase envelope tracers: raw composition derivatives vs finite differences", "[phase_envelope][tracers][diag]") {
    // Isolates dpdxj__constT_V_xi and dln_fugacity_dxj__constT_rho_xi (XN_DEPENDENT) at the
    // incipient liquid of the quaternary start point.
    auto HEOS = make_heos("Nitrogen&Methane&Ethane&Propane", {0.10, 0.34, 0.41, 0.15});
    PhaseEnvelopeTracers::Options opts;
    auto s0 = PhaseEnvelopeTracers::starting_point(*HEOS, opts);
    HelmholtzEOSMixtureBackend& inc = *HEOS->SatL;
    const std::size_t N = s0.x.size();
    auto eval = [&](const std::vector<CoolPropDbl>& x, std::size_t i, double& p, double& lnf) {
        inc.set_mole_fractions(x);
        inc.update(DmolarT_INPUTS, s0.rhomolar_liq, s0.T);
        p = inc.p();
        lnf = std::log(MixtureDerivatives::fugacity_i(inc, i, XN_DEPENDENT));
    };
    for (std::size_t j = 0; j + 1 < N; ++j) {
        if (s0.x[j] < 1e-4) {
            continue;  // a symmetric step would push a trace component negative
        }
        double p0 = 0, lnf0 = 0;
        eval(s0.x, 0, p0, lnf0);
        const double dp_an = MixtureDerivatives::dpdxj__constT_V_xi(inc, j, XN_DEPENDENT);
        const double dlnf_an = MixtureDerivatives::dln_fugacity_dxj__constT_rho_xi(inc, 0, j, XN_DEPENDENT);
        const double h = 1e-6;
        std::vector<CoolPropDbl> xp = s0.x, xm = s0.x;
        xp[j] += h;
        xp[N - 1] -= h;
        xm[j] -= h;
        xm[N - 1] += h;
        double pp = 0, lp = 0, pm = 0, lm = 0;
        eval(xp, 0, pp, lp);
        eval(xm, 0, pm, lm);
        const double dp_fd = (pp - pm) / (2 * h), dlnf_fd = (lp - lm) / (2 * h);
        CAPTURE(j, dp_an, dp_fd, dlnf_an, dlnf_fd, s0.x[j]);
        CHECK(std::abs(dp_fd - dp_an) <= 1e-5 * std::abs(dp_an) + 1e-3);
        CHECK(std::abs(dlnf_fd - dlnf_an) <= 1e-5 * std::abs(dlnf_an) + 1e-7);
    }
}

TEST_CASE("Phase envelope tracers: diagnostic trace from environment", "[phase_envelope][tracers][diag][.]") {
    // Hidden test: COOLPROP_TRACER_DIAG="Methane&Ethane|0.85,0.15|lnK_density[|debuglevel]" prints the full trace.
    // [.] hides this from a default run, but an explicit tag filter (preflight passes
    // [phase_envelope]) still selects it, so an absent variable must skip rather than fail.
    const char* spec = std::getenv("COOLPROP_TRACER_DIAG");
    if (spec == nullptr) {
        SKIP("set COOLPROP_TRACER_DIAG=\"<fluids>|<z or ->|<algorithm>[|debuglevel]\" to run this");
    }
    std::vector<std::string> parts = strsplit(spec, '|');
    REQUIRE(parts.size() >= 3);
    std::vector<double> z;
    if (parts[1] != "-" && !parts[1].empty()) {  // "-" means a predefined mixture, which carries its own composition
        for (const auto& s : strsplit(parts[1], ',')) {
            z.push_back(std::stod(s));
        }
    }
    AlgorithmGuard guard(parts[2]);
    std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", parts[0]));
    if (!z.empty()) {
        AS->set_mole_fractions(z);
    }
    set_debug_level(parts.size() > 3 ? std::stoi(parts[3]) : 1);
    try {
        AS->build_phase_envelope("");
    } catch (std::exception& e) {
        std::cout << "EXCEPTION: " << e.what() << '\n';
    }
    set_debug_level(0);
    const PhaseEnvelopeData& env = AS->get_phase_envelope_data();
    std::cout << "built=" << env.built << " closed=" << env.closed << " n=" << env.T.size() << " icrit=" << static_cast<long long>(env.icrit) << '\n';
    for (std::size_t k = 0; k < env.T.size(); ++k) {
        std::cout << format("%4d T=%10.4f p=%14.6g rho_inc=%12.5g rho_feed=%12.5g x=", static_cast<int>(k), env.T[k], env.p[k], env.rhomolar_liq[k],
                            env.rhomolar_vap[k]);
        for (const auto& xj : env.x) {
            std::cout << format("%9.4g ", xj[k]);
        }
        std::cout << '\n';
    }
}

TEST_CASE("Phase envelope tracers: configuration key", "[phase_envelope][tracers]") {
    CHECK(get_config_string(PHASE_ENVELOPE_ALGORITHM) == "legacy");
    auto HEOS = make_heos("Methane&Ethane", {0.85, 0.15});
    AlgorithmGuard guard("no_such_tracer");
    CHECK_THROWS_AS(HEOS->build_phase_envelope(""), ValueError);
}

TEST_CASE("Phase envelope tracers: Jacobians match finite differences", "[phase_envelope][tracers]") {
    for (const char* alg : {"lnK_density", "lnK_pressure"}) {
        SECTION(std::string(alg) + " CH4/C2") {
            check_jacobian("Methane&Ethane", {0.85, 0.15}, alg);
        }
        SECTION(std::string(alg) + " N2/CH4/C2/C3") {
            check_jacobian("Nitrogen&Methane&Ethane&Propane", {0.10, 0.34, 0.41, 0.15}, alg);
        }
    }
}

TEST_CASE("Phase envelope tracers: methane/ethane closes and matches blind dew points", "[phase_envelope][tracers]") {
    const std::vector<double> z = {0.85, 0.15};
    for (const char* alg : {"lnK_density", "lnK_pressure"}) {
        SECTION(alg) {
            AlgorithmGuard guard(alg);
            auto HEOS = make_heos("Methane&Ethane", z);
            REQUIRE_NOTHROW(HEOS->build_phase_envelope(""));
            const PhaseEnvelopeData& env = HEOS->get_phase_envelope_data();
            CHECK(env.T.size() > 50);
            if (std::string(alg) == "lnK_density") {
                CHECK(env.closed);
                CHECK(env.built);
            } else {
                WARN(alg << ": closed=" << env.closed << " (the (T,p) form loses the density root near the critical point)");
            }
            REQUIRE(env.icrit < env.T.size());
            REQUIRE(env.icrit > 5);
            CHECK(env.T.size() > 30);

            // Every stored dew-side point must agree with a blind PQ flash on a separate instance
            std::shared_ptr<AbstractState> blind(AbstractState::factory("HEOS", "Methane&Ethane"));
            blind->set_mole_fractions(z);
            const std::size_t icrit = env.icrit;
            for (std::size_t k = 2; k + 4 < icrit; k += std::max<std::size_t>(1, icrit / 6)) {
                CAPTURE(k, env.p[k], env.T[k]);
                REQUIRE_NOTHROW(blind->update(PQ_INPUTS, env.p[k], 1.0));
                CHECK(blind->T() == Catch::Approx(env.T[k]).epsilon(1e-6));
            }

            // Critical point estimate from the ln K crossing against the exact locator
            double Tc = 0, pc = 0;
            critical_estimate(env, Tc, pc);
            auto crit = HEOS->all_critical_points();
            REQUIRE(!crit.empty());
            CHECK(std::abs(Tc - crit[0].T) < 0.5);
            CHECK(std::abs(pc / crit[0].p - 1) < 0.01);
        }
    }
}

TEST_CASE("Phase envelope tracers: natural-gas-like quaternary terminates at Tmin, not silently", "[phase_envelope][tracers]") {
    const std::vector<double> z = {0.10, 0.34, 0.41, 0.15};
    for (const char* alg : {"lnK_density", "lnK_pressure"}) {
        SECTION(alg) {
            AlgorithmGuard guard(alg);
            auto HEOS = make_heos("Nitrogen&Methane&Ethane&Propane", z);
            REQUIRE_NOTHROW(HEOS->build_phase_envelope(""));
            const PhaseEnvelopeData& env = HEOS->get_phase_envelope_data();
            // This one need not close; what matters is that it traced a usable stretch and said
            // why it stopped, where the legacy tracer used to return nothing and raise nothing.
            CHECK(env.T.size() > 50);
            CHECK(env.icrit < env.T.size());
            double pmax = 0;
            for (double p : env.p) {
                REQUIRE(std::isfinite(p));
                pmax = std::max(pmax, p);
            }
            CHECK(pmax < 20e6);
            // Whatever ends the trace, it is reported rather than swallowed.  The legacy
            // tracer returns here with built = false and no error at all.
            const std::string& stop = PhaseEnvelopeTracers::last_stop_reason();
            CAPTURE(stop, env.T.back(), env.T.size());
            CHECK(stop != "none");
            CHECK(stop != "max_points");
            CHECK(env.T.size() > 50);
            CHECK(env.T.size() < 1000);
        }
    }
}

TEST_CASE("Phase envelope tracers: trivial solution and runaway pressure are rejected", "[phase_envelope][tracers]") {
    struct Case
    {
        const char* fluids;
        std::vector<double> z;
        double
          pmax_allowed;  ///< both have a genuine open high-pressure branch (legacy reaches 1516 and 8656 MPa); the trace may take one step past the 1 GPa ceiling before stopping, but must never reach a trivial solution
    };
    const std::vector<Case> cases = {{"Methane&n-Decane", {0.7, 0.3}, 1.1e9},
                                     {"CarbonDioxide&Nitrogen&Oxygen&Argon", {0.9, 0.05, 0.03, 0.02}, 1.1e9}};
    for (const char* alg : {"lnK_density", "lnK_pressure"}) {
        for (const auto& c : cases) {
            SECTION(std::string(alg) + " " + c.fluids) {
                AlgorithmGuard guard(alg);
                auto HEOS = make_heos(c.fluids, c.z);
                REQUIRE_NOTHROW(HEOS->build_phase_envelope(""));
                const PhaseEnvelopeData& env = HEOS->get_phase_envelope_data();
                CHECK(env.T.size() > 20);
                // A bounded stop, not a runaway: the ceiling is what ends an open branch
                const std::string& stop = PhaseEnvelopeTracers::last_stop_reason();
                CAPTURE(stop);
                CHECK(stop != "none");
                CHECK(stop != "max_points");
                for (std::size_t k = 0; k < env.T.size(); ++k) {
                    CAPTURE(k, env.T[k], env.p[k]);
                    REQUIRE(std::isfinite(env.p[k]));
                    CHECK(env.p[k] < c.pmax_allowed);
                    double maxlnK = 0;
                    for (const auto& lnK : env.lnK) {
                        maxlnK = std::max(maxlnK, std::abs(lnK[k]));
                    }
                    CHECK(maxlnK > 1e-8);
                }
            }
        }
    }
}

TEST_CASE("Phase envelope tracers: heavy trace components stay finite", "[phase_envelope][tracers]") {
    const std::string fluids =
      "Methane&Nitrogen&CarbonDioxide&Ethane&Propane&IsoButane&n-Butane&Isopentane&n-Pentane&n-Hexane&n-Heptane&n-Octane&Hydrogen&Helium&Oxygen";
    const std::vector<double> z = {0.85, 0.02, 0.01, 0.04, 0.015, 0.004, 0.004, 0.003, 0.003, 0.001, 0.001, 0.001, 0.02, 0.01, 0.018};
    for (const char* alg : {"lnK_density", "lnK_pressure"}) {
        SECTION(alg) {
            AlgorithmGuard guard(alg);
            auto HEOS = make_heos(fluids, z);
            REQUIRE_NOTHROW(HEOS->build_phase_envelope(""));
            const PhaseEnvelopeData& env = HEOS->get_phase_envelope_data();
            CHECK(env.T.size() > 20);
            for (const auto& xj : env.x) {
                for (double v : xj) {
                    REQUIRE(std::isfinite(v));
                    CHECK(v >= 0);
                    CHECK(v <= 1);
                }
            }
        }
    }
}

TEST_CASE("Phase envelope tracers: start pressure fallback", "[phase_envelope][tracers]") {
    // R508A and Amarillo cannot solve the 100 Pa dew point; the tracers retry a decade higher.
    for (const char* mix : {"R508A.mix", "Amarillo.mix"}) {
        for (const char* alg : {"lnK_density", "lnK_pressure"}) {
            SECTION(std::string(alg) + " " + mix) {
                AlgorithmGuard guard(alg);
                std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", mix));
                REQUIRE_NOTHROW(AS->build_phase_envelope(""));
                const PhaseEnvelopeData& env = AS->get_phase_envelope_data();
                CHECK(env.T.size() > 20);
                CHECK(env.p.front() > 100.0);
            }
        }
    }
}

#endif  // ENABLE_CATCH
