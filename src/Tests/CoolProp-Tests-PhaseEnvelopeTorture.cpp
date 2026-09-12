// Torture corpus for the isopleth tracers: every unique predefined mixture plus a set of
// hard binaries and multicomponent cases (azeotropes, open envelopes, wide-boiling pairs,
// heavy trace components, near-pure limits), traced with every algorithm accepted by
// PHASE_ENVELOPE_ALGORITHM.
//
// Prints one line per mixture per algorithm and a per-algorithm tally.  When the
// environment variable COOLPROP_PHASE_ENVELOPE_TORTURE_CSV names a file, the full table is
// also written there.  The legacy tally on the predefined subset is pinned to the baseline
// measured on 2026-09-11 so that it cannot regress; candidate algorithms are informational
// until promoted, but every one must finish each mixture with finite stored values and only
// CoolProp exceptions.

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/Configuration.h"
#    include "CoolProp/CoolProp.h"
#    include "CoolProp/DataStructures.h"
#    include "CoolProp/detail/strings.h"
#    include "CoolProp/detail/tools.h"
#    include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#    include "Backends/Helmholtz/MixtureDerivatives.h"
#    include "Backends/Helmholtz/PhaseEnvelopeTracers.h"

#    include <algorithm>
#    include <chrono>
#    include <cmath>
#    include <limits>
#    include <cstdlib>
#    include <fstream>
#    include <iostream>
#    include <map>
#    include <memory>
#    include <numeric>
#    include <string>
#    include <vector>

using namespace CoolProp;

namespace {

struct CorpusCase
{
    std::string label;
    std::string fluids;     ///< "&"-joined component names, or a predefined mixture name
    std::vector<double> z;  ///< empty for a predefined mixture
    bool predefined = false;
};

std::vector<CorpusCase> build_corpus() {
    std::vector<CorpusCase> corpus;
    // Predefined mixtures, deduplicated case-insensitively (the list carries FOO.MIX and Foo.mix twins)
    std::vector<std::string> names = strsplit(get_global_param_string("predefined_mixtures"), ',');
    std::map<std::string, std::string> seen;
    for (const auto& n : names) {
        std::string key = n;
        std::transform(key.begin(), key.end(), key.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        if (seen.insert({key, n}).second) {
            corpus.push_back({n, n, {}, true});
        }
    }
    // Hard cases
    const std::vector<std::pair<std::string, std::vector<double>>> hard = {
      {"Methane&Ethane", {0.02, 0.98}},
      {"Methane&Ethane", {0.2, 0.8}},
      {"Methane&Ethane", {0.5, 0.5}},
      {"Methane&Ethane", {0.8, 0.2}},
      {"Methane&Ethane", {0.98, 0.02}},
      {"Methane&Ethane", {0.99, 0.01}},
      {"Methane&Ethane", {0.999, 0.001}},
      {"Ethane&CarbonDioxide", {0.3, 0.7}},  // azeotrope
      {"Ethane&CarbonDioxide", {0.5, 0.5}},
      {"Ethane&CarbonDioxide", {0.7, 0.3}},
      {"R125&R143a", {0.5, 0.5}},  // azeotropic (R507A)
      {"R23&R116", {0.5, 0.5}},    // azeotropic (R508)
      {"R32&R125", {0.5, 0.5}},
      {"R1234yf&R134a", {0.5, 0.5}},
      {"R1234ze(E)&R32", {0.5, 0.5}},
      {"R32&Propane", {0.5, 0.5}},
      {"Propane&R134a", {0.5, 0.5}},
      {"Propane&IsoButane", {0.5, 0.5}},
      {"Propylene&Propane", {0.5, 0.5}},
      {"Ethane&Propane", {0.5, 0.5}},
      {"n-Butane&n-Pentane", {0.5, 0.5}},
      {"Nitrogen&Argon", {0.5, 0.5}},
      {"Neon&Argon", {0.5, 0.5}},
      {"Nitrogen&Oxygen&Argon", {0.78, 0.21, 0.01}},
      {"CarbonDioxide&Nitrogen", {0.5, 0.5}},  // runaway pressure in legacy
      {"CarbonDioxide&Nitrogen", {0.9, 0.1}},
      {"CarbonDioxide&Argon", {0.5, 0.5}},
      {"CarbonDioxide&Ethane&Propane", {0.5, 0.3, 0.2}},
      {"CarbonDioxide&Nitrogen&Oxygen&Argon", {0.9, 0.05, 0.03, 0.02}},
      {"CarbonDioxide&Water", {0.98, 0.02}},  // start fails in legacy
      {"Methane&Nitrogen", {0.5, 0.5}},
      {"Methane&Hydrogen", {0.9, 0.1}},  // open envelope (Type III)
      {"Hydrogen&Nitrogen", {0.5, 0.5}},
      {"Nitrogen&Helium", {0.9, 0.1}},
      {"Methane&n-Hexane", {0.9, 0.1}},  // wide boiling
      {"Methane&n-Decane", {0.7, 0.3}},
      {"Methane&n-Decane", {0.3, 0.7}},
      {"Ethanol&Water", {0.5, 0.5}},
      {"Nitrogen&Methane&Ethane&Propane", {0.10, 0.34, 0.41, 0.15}},
      {"Methane&Ethane&Propane&n-Butane&n-Pentane", {0.9, 0.05, 0.03, 0.015, 0.005}},
      {"Methane&Nitrogen&CarbonDioxide&Ethane&Propane&IsoButane&n-Butane&Isopentane&n-Pentane&n-Hexane&n-Heptane&n-Octane&Hydrogen&Helium&Oxygen",
       {0.85, 0.02, 0.01, 0.04, 0.015, 0.004, 0.004, 0.003, 0.003, 0.001, 0.001, 0.001, 0.02, 0.01, 0.018}},
    };
    for (const auto& h : hard) {
        std::string label = h.first;
        if (label.size() > 40) {
            label = label.substr(0, 37) + "...";
        }
        label += " [";
        for (std::size_t i = 0; i < h.second.size(); ++i) {
            label += (i > 0 ? "," : "") + format("%g", h.second[i]);
        }
        label += "]";
        corpus.push_back({label, h.first, h.second, false});
    }
    return corpus;
}

struct Row
{
    std::string label, algorithm, error, stop;
    bool predefined = false, constructed = false, built = false, closed = false, finite = true, coolprop_error = true;
    std::size_t n = 0;
    long long icrit = -1;
    std::size_t dev_samples = 0, dev_failures = 0, fug_samples = 0;
    double fug = -1;  ///< max |ln f_i(liq) - ln f_i(vap)| over sampled stored points; -1 when unavailable
    bool dev_nonfinite = false;
    double seconds = 0, pmax = 0, Tmin = 0, Tmax = 0,
           dev = -1;  ///< max relative pressure deviation vs a blind QT flash on both branches, -1 when unavailable
};

Row run_case(const CorpusCase& c, const std::string& algorithm) {
    Row r;
    r.label = c.label;
    r.algorithm = algorithm;
    r.predefined = c.predefined;
    std::shared_ptr<AbstractState> AS, blind;
    try {
        AS.reset(AbstractState::factory("HEOS", c.fluids));
        blind.reset(AbstractState::factory("HEOS", c.fluids));
        if (!c.z.empty()) {
            AS->set_mole_fractions(c.z);
            blind->set_mole_fractions(c.z);
        }
        r.constructed = true;
    } catch (std::exception& e) {
        r.error = std::string("construct: ") + e.what();
        return r;
    }
    set_config_string(PHASE_ENVELOPE_ALGORITHM, algorithm);
    const auto t0 = std::chrono::steady_clock::now();
    try {
        AS->build_phase_envelope("");
    } catch (CoolPropBaseError& e) {
        r.error = e.what();
    } catch (std::exception& e) {
        r.error = std::string("NON-COOLPROP: ") + e.what();
        r.coolprop_error = false;
    }
    r.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    set_config_string(PHASE_ENVELOPE_ALGORITHM, "legacy");

    const PhaseEnvelopeData& env = AS->get_phase_envelope_data();
    r.n = env.T.size();
    r.built = env.built;
    r.closed = env.closed;
    r.stop = algorithm == "legacy" ? (env.closed ? "closed" : "-") : PhaseEnvelopeTracers::last_stop_reason();
    r.icrit = env.icrit < env.T.size() ? static_cast<long long>(env.icrit) : -1;
    if (r.n > 0) {
        r.pmax = *std::max_element(env.p.begin(), env.p.end());
        r.Tmin = *std::min_element(env.T.begin(), env.T.end());
        r.Tmax = *std::max_element(env.T.begin(), env.T.end());
        for (std::size_t k = 0; k < r.n; ++k) {
            // The derived arrays are where corruption actually lands: store_variables computes
            // K = y/x, ln K, ln T and ln p, so a merely-finite x of exactly zero yields an
            // infinite K and ln K while T, p and the densities all still look clean.  Checking
            // only the state variables would let that through and report the run as finite.
            r.finite = r.finite && std::isfinite(env.T[k]) && std::isfinite(env.p[k]) && std::isfinite(env.rhomolar_liq[k])
                       && std::isfinite(env.rhomolar_vap[k]) && std::isfinite(env.lnT[k]) && std::isfinite(env.lnp[k]);
            for (const auto& xj : env.x) {
                r.finite = r.finite && std::isfinite(xj[k]);
            }
            for (const auto& yj : env.y) {
                r.finite = r.finite && std::isfinite(yj[k]);
            }
            for (const auto& Kj : env.K) {
                r.finite = r.finite && std::isfinite(Kj[k]);
            }
            for (const auto& lnKj : env.lnK) {
                r.finite = r.finite && std::isfinite(lnKj[k]);
            }
        }
        // Consistency: stored points against a blind QT flash on a separate instance.  BOTH
        // branches are sampled, and within each branch the LAST matching points as well as the
        // first: a false closure manifests where the trace ends, so sampling only the opening
        // stretch of each branch is exactly the wrong place to look.  The branch label comes
        // from env.Q, which store_variables sets from the density ordering (1 where the
        // incipient phase is the denser one, i.e. a dew point), so it is correct for every
        // algorithm and assumes nothing about point ordering.
        auto sample = [&](double Q) {
            std::vector<std::size_t> idx;
            for (std::size_t k = 1; k + 1 < r.n; ++k) {
                if (env.Q[k] == Q) {
                    idx.push_back(k);
                }
            }
            if (idx.empty()) {
                return;
            }
            const std::size_t picks[] = {idx.front(), idx[idx.size() / 2], idx.back()};
            for (std::size_t k : picks) {
                try {
                    blind->update(QT_INPUTS, Q, env.T[k]);
                    const double d = std::abs(blind->p() / env.p[k] - 1);
                    // std::max(-1.0, NaN) returns -1.0, so a NaN deviation would leave dev at its
                    // "unavailable" sentinel and take the silent exit below.  Reject it explicitly.
                    if (std::isfinite(d)) {
                        r.dev = std::max(r.dev, d);
                    } else {
                        r.dev_nonfinite = true;
                    }
                    ++r.dev_samples;
                } catch (...) {  // NOLINT(bugprone-empty-catch)
                    ++r.dev_failures;
                }
            }
        };
        sample(1.0);  // dew branch: incipient phase is the liquid
        sample(0.0);  // bubble branch: incipient phase is the vapor

        // Self-contained correctness of the stored points.  Equality of component fugacities at
        // (T, rho', x) and (T, rho'', y) IS the definition of a phase-boundary point, so this
        // says whether each stored point is really on the boundary without asking any other
        // solver.  The blind-flash comparison above is a different question -- whether the
        // envelope traced the right branch -- and where the two disagree it is worth knowing
        // which of the two routines is at fault.
        if (auto* heos = dynamic_cast<HelmholtzEOSMixtureBackend*>(AS.get())) {
            if (heos->SatL && heos->SatV) {
                const std::size_t stride = std::max<std::size_t>(1, r.n / 8);
                for (std::size_t k = 1; k + 1 < r.n; k += stride) {
                    std::vector<CoolPropDbl> xk, yk;
                    for (const auto& xj : env.x) {
                        xk.push_back(xj[k]);
                    }
                    for (const auto& yj : env.y) {
                        yk.push_back(yj[k]);
                    }
                    try {
                        heos->SatL->set_mole_fractions(xk);
                        heos->SatL->update(DmolarT_INPUTS, env.rhomolar_liq[k], env.T[k]);
                        heos->SatV->set_mole_fractions(yk);
                        heos->SatV->update(DmolarT_INPUTS, env.rhomolar_vap[k], env.T[k]);
                        double worst = 0;
                        for (std::size_t i = 0; i < xk.size(); ++i) {
                            const double lnfl = std::log(MixtureDerivatives::fugacity_i(*heos->SatL, i, XN_DEPENDENT));
                            const double lnfv = std::log(MixtureDerivatives::fugacity_i(*heos->SatV, i, XN_DEPENDENT));
                            const double d = std::abs(lnfl - lnfv);
                            if (!std::isfinite(d)) {
                                worst = std::numeric_limits<double>::infinity();
                                break;
                            }
                            worst = std::max(worst, d);
                        }
                        r.fug = std::max(r.fug, worst);
                        ++r.fug_samples;
                    } catch (...) {  // NOLINT(bugprone-empty-catch)
                        // EOS could not be evaluated at this stored point; leave fug as is
                    }
                }
            }
        }
    }
    return r;
}

std::string row_line(const Row& r) {
    return format("%-58s %-13s %s %s %s %-10s n=%4d icrit=%4lld t=%7.3fs pmax=%10.3f MPa T=[%6.1f,%6.1f] dev=%s fug=%s %s",
                  r.label.substr(0, 58).c_str(), r.algorithm.c_str(), r.constructed ? "C" : "-", r.built ? "B" : "-", r.closed ? "K" : "-",
                  r.stop.c_str(), static_cast<int>(r.n), r.icrit, r.seconds, r.pmax / 1e6, r.Tmin, r.Tmax,
                  r.dev < 0 ? "  n/a  " : format("%7.1e", r.dev).c_str(), r.fug < 0 ? "  n/a  " : format("%7.1e", r.fug).c_str(),
                  r.error.substr(0, 60).c_str());
}

}  // namespace

TEST_CASE("Phase envelope torture corpus: all predefined mixtures and hard cases, every algorithm", "[phase_envelope][torture][slow]") {
    const std::vector<CorpusCase> corpus = build_corpus();
    const std::vector<std::string> algorithms = PhaseEnvelopeTracers::algorithms();
    std::vector<Row> rows;
    rows.reserve(corpus.size() * algorithms.size());

    std::cout << "\n=== Phase envelope torture corpus: " << corpus.size() << " mixtures x " << algorithms.size() << " algorithms ===\n";
    for (const auto& c : corpus) {
        for (const auto& alg : algorithms) {
            // Unbuffered progress on stderr: stdout is block-buffered when redirected, so a case
            // that runs away leaves no trace of which mixture it was stuck on.
            // std::cerr is unit-buffered, so this reaches the terminal immediately
            std::cerr << "[torture] " << c.label.substr(0, 48) << " / " << alg << '\n';
            rows.push_back(run_case(c, alg));
            std::cout << row_line(rows.back()) << '\n';
        }
    }

    // Tallies
    struct Tally
    {
        std::size_t constructed = 0, built = 0, closed = 0, complete = 0, consistent = 0, checked = 0, closed_consistent = 0, false_closure = 0,
                    closed_unverified = 0, fug_checked = 0, traced = 0, closed_unmeasured = 0, predefined_constructed = 0, predefined_closed = 0;
        std::vector<double> seconds;
        double fug_max = 0;  ///< worst fugacity-equality residual over every sampled stored point
    };
    std::map<std::string, Tally> tally;
    for (const auto& r : rows) {
        Tally& t = tally[r.algorithm];
        if (!r.constructed) continue;
        ++t.constructed;
        t.seconds.push_back(r.seconds);
        // "traced" is the honest reach metric: a run that produced enough points to plot or
        // inspect, whether or not it closed and whether or not `built` was set.  Before this
        // work 24 of these returned nothing at all AND raised nothing.
        if (r.n >= 20) {
            ++t.traced;
        }
        if (r.built) ++t.built;
        if (r.closed) ++t.closed;
        if (r.closed || r.stop == "floor" || r.stop == "degenerate" || r.stop == "pure") ++t.complete;
        // A closed envelope the blind flash could never evaluate is NOT evidence of correctness.
        // Counting it in neither bucket would let the worst case -- a closure so wrong that the
        // flash fails everywhere on it -- slip past a ceiling-only false-closure pin.
        if (r.closed && (r.dev < 0 || r.dev_nonfinite)) {
            ++t.closed_unverified;
        }
        if (r.fug >= 0) {
            ++t.fug_checked;
            t.fug_max = std::max(t.fug_max, r.fug);
        } else if (r.closed) {
            // A closed envelope that could not be fugacity-checked at all is not evidence of
            // anything; count it so an unmeasurable run cannot masquerade as a clean one.
            ++t.closed_unmeasured;
        }
        if (r.dev >= 0) {
            ++t.checked;
            if (r.dev < 1e-3) ++t.consistent;
            // A closed envelope whose own dew points disagree with a blind flash is a FALSE
            // closure: the trace closed some small loop rather than the real boundary.  Closure
            // alone is therefore not the quality metric; closed AND consistent is.
            if (r.closed) {
                if (r.dev < 1e-3) {
                    ++t.closed_consistent;
                } else {
                    ++t.false_closure;
                }
            }
        }
        if (r.predefined) {
            ++t.predefined_constructed;
            if (r.closed) ++t.predefined_closed;
        }
    }
    std::cout << "\n--- tally (C constructed, B built, K closed) ---\n";
    for (auto& kv : tally) {
        Tally& t = kv.second;
        std::sort(t.seconds.begin(), t.seconds.end());
        const double median = t.seconds.empty() ? 0 : t.seconds[t.seconds.size() / 2];
        const double total = std::accumulate(t.seconds.begin(), t.seconds.end(), 0.0);
        std::cout << format("%-13s constructed=%3d traced=%3d built=%3d closed=%3d closed+consistent=%3d FALSE-closures=%2d unverified=%2d "
                            "complete=%3d consistent=%3d/%3d "
                            "predefined closed=%3d/%3d worst-fugacity=%8.2e median=%.4fs total=%.1fs\n",
                            kv.first.c_str(), static_cast<int>(t.constructed), static_cast<int>(t.traced), static_cast<int>(t.built),
                            static_cast<int>(t.closed), static_cast<int>(t.closed_consistent), static_cast<int>(t.false_closure),
                            static_cast<int>(t.closed_unverified), static_cast<int>(t.complete), static_cast<int>(t.consistent),
                            static_cast<int>(t.checked), static_cast<int>(t.predefined_closed), static_cast<int>(t.predefined_constructed), t.fug_max,
                            median, total);
    }

    std::cout << "\n--- closed envelopes that DISAGREE with a blind QT flash (diagnostic; see the fugacity residual before blaming the tracer) ---\n";
    for (const auto& r : rows) {
        if (r.closed && r.dev >= 1e-3) {
            std::cout << format("%-58s %-13s dev=%8.4f n=%4d pmax=%9.3f MPa\n", r.label.substr(0, 58).c_str(), r.algorithm.c_str(), r.dev,
                                static_cast<int>(r.n), r.pmax / 1e6);
        }
    }

    if (const char* csv = std::getenv("COOLPROP_PHASE_ENVELOPE_TORTURE_CSV")) {
        std::ofstream f(csv);
        if (!f) {
            std::cout << "could not open " << csv << " for writing\n";
            return;
        }
        f << "label,algorithm,predefined,constructed,built,closed,stop,n,icrit,seconds,pmax_Pa,Tmin_K,Tmax_K,dev,dev_samples,dev_failures,fug,fug_"
             "samples,error\n";
        for (const auto& r : rows) {
            std::string err = r.error;
            std::replace(err.begin(), err.end(), ',', ';');
            std::replace(err.begin(), err.end(), '\n', ' ');
            f << '"' << r.label << "\"," << r.algorithm << ',' << r.predefined << ',' << r.constructed << ',' << r.built << ',' << r.closed << ','
              << r.stop << ',' << r.n << ',' << r.icrit << ',' << r.seconds << ',' << r.pmax << ',' << r.Tmin << ',' << r.Tmax << ',' << r.dev << ','
              << r.dev_samples << ',' << r.dev_failures << ',' << r.fug << ',' << r.fug_samples << ",\"" << err << "\"\n";
        }
        f.flush();
        std::cout << (f ? "wrote " : "FAILED to write ") << csv << '\n';
    }

    // No algorithm may store a non-finite value.  This started as a pinned known defect of the
    // DEFAULT tracer (NaN mole fractions on a wide-boiling 15-component gas); the legacy insert
    // sites now validate everything they store, so the count is zero and stays zero.
    for (const auto& r : rows) {
        if (r.constructed && !r.finite) {
            std::cout << "NON-FINITE: " << r.label << " / " << r.algorithm << '\n';
        }
        CAPTURE(r.label, r.algorithm);
        CHECK((!r.constructed || r.finite));
    }

    // Pins.  Legacy on the predefined subset: baseline 2026-09-11 was 116 constructed, 106 closed.
    const Tally& legacy = tally["legacy"];
    CHECK(legacy.predefined_constructed >= 116);
    CHECK(legacy.predefined_closed >= 107);
    // Quality pins measured 2026-09-11.  Closure alone is not enough: a false closure is a
    // silently wrong envelope, which is worse than an honest failure, so it is bounded too.
    // `built` keeps its old meaning (a closed, interpolatable envelope), so it barely moves;
    // the reach improvement shows up in `traced`, which is 154 here against 131 on master.
    CHECK(legacy.built >= 134);
    CHECK(legacy.closed >= 131);
    CHECK(legacy.traced >= 154);
    CHECK(tally["lnK_density"].traced >= 150);
    CHECK(tally["lnK_pressure"].traced >= 145);
    for (auto& kv : tally) {
        CAPTURE(kv.first);
        CHECK(kv.second.closed_unverified == 0);
        CHECK(kv.second.closed_unmeasured == 0);
    }

    // Correctness gate.  Equality of component fugacities is the definition of a phase-boundary
    // point, so it is the one measure here that needs no second solver.  Measured 2026-09-12:
    // legacy 8.4e-4 worst, both candidates 1e-9.  Note what each bound actually protects: for
    // the DEFAULT this largely re-checks its own store-time gate (which rejects above 1e-3), so
    // it mainly guards against that gate being removed or bypassed; for the two candidates,
    // which have no store-time gate, it is the only thing standing between a bad point and the
    // stored envelope, so they are held an order of magnitude tighter.  The count of checked
    // rows is pinned too, so an algorithm that produced nothing cannot pass by being
    // unmeasurable.
    for (auto& kv : tally) {
        CAPTURE(kv.first, kv.second.fug_max);
        CHECK(kv.second.fug_max < (kv.first == "legacy" ? 5e-3 : 1e-6));
        CHECK(kv.second.fug_checked >= 150);
    }

    // Agreement with a blind QT flash is a DIAGNOSTIC, not a correctness gate.  Every point
    // behind these "disagreements" satisfies fugacity equality to ~1e-13, so a disagreement
    // means the tracer and the flash landed on different states at the same (T, Q) -- an
    // incomplete envelope, or a flash that found another root -- and does not by itself show
    // the envelope is wrong.  Bounded generously so a large regression is still visible.
    CHECK(legacy.closed_consistent >= 105);
    CHECK(legacy.false_closure <= 26);
    CHECK(tally["lnK_density"].closed_consistent >= 110);
    CHECK(tally["lnK_density"].false_closure <= 21);
    // Floors for every algorithm, so a candidate that regressed to producing nothing at all
    // cannot pass: the per-row checks below are all vacuously true for an empty envelope.
    // `built` now means "closed" for every algorithm, so reach is floored through `traced`.
    CHECK(tally["lnK_density"].closed >= 125);
    CHECK(tally["lnK_pressure"].constructed >= 157);
    // Every algorithm: finite stored values, CoolProp exceptions only, bounded point count.
    for (const auto& r : rows) {
        CAPTURE(r.label, r.algorithm, r.error);
        CHECK(r.coolprop_error);
        // The tracers cap themselves at Options::max_points (1000, +2 from finalize's maxima).
        // Legacy has no absolute cap, only refine's 4x growth limit, so this is an empirical
        // bound: the observed maximum over the whole corpus is 224, and 512 leaves room for
        // fluid-data churn while still firing on a runaway.
        CHECK(r.n <= (r.algorithm == "legacy" ? 512u : 1002u));
    }
}

#endif  // ENABLE_CATCH
