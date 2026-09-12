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
#    include "Backends/Helmholtz/PhaseEnvelopeTracers.h"

#    include <algorithm>
#    include <chrono>
#    include <cmath>
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
            label += (i ? "," : "") + format("%g", h.second[i]);
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
            r.finite = r.finite && std::isfinite(env.T[k]) && std::isfinite(env.p[k]) && std::isfinite(env.rhomolar_liq[k])
                       && std::isfinite(env.rhomolar_vap[k]);
            for (const auto& xj : env.x) {
                r.finite = r.finite && std::isfinite(xj[k]);
            }
        }
        // Consistency: stored points against a blind QT flash on a separate instance.  BOTH
        // branches are sampled; checking only the dew side misses a bubble branch that has
        // wandered off the boundary, which is exactly how a false closure arises.  The branch
        // of each stored point comes from env.Q, which store_variables sets from the density
        // ordering (1 where the incipient phase is the denser one, i.e. a dew point), so it is
        // correct for every algorithm and needs no assumption about point ordering.
        auto sample = [&](double Q) {
            std::size_t checked = 0;
            const std::size_t stride = std::max<std::size_t>(1, r.n / 12);
            for (std::size_t k = 1; k + 1 < r.n && checked < 3; k += stride) {
                if (env.Q[k] != Q) {
                    continue;
                }
                try {
                    blind->update(QT_INPUTS, Q, env.T[k]);
                    r.dev = std::max(r.dev, std::abs(blind->p() / env.p[k] - 1));
                    ++checked;
                } catch (...) {  // NOLINT(bugprone-empty-catch)
                    // blind flash unavailable at this point; leave dev as is
                }
            }
        };
        sample(1.0);  // dew branch: incipient phase is the liquid
        sample(0.0);  // bubble branch: incipient phase is the vapor
    }
    return r;
}

std::string row_line(const Row& r) {
    return format("%-58s %-13s %s %s %s %-10s n=%4d icrit=%4lld t=%7.3fs pmax=%10.3f MPa T=[%6.1f,%6.1f] dev=%s %s", r.label.substr(0, 58).c_str(),
                  r.algorithm.c_str(), r.constructed ? "C" : "-", r.built ? "B" : "-", r.closed ? "K" : "-", r.stop.c_str(), static_cast<int>(r.n),
                  r.icrit, r.seconds, r.pmax / 1e6, r.Tmin, r.Tmax, r.dev < 0 ? "  n/a  " : format("%7.1e", r.dev).c_str(),
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
            rows.push_back(run_case(c, alg));
            std::cout << row_line(rows.back()) << '\n';
        }
    }

    // Tallies
    struct Tally
    {
        std::size_t constructed = 0, built = 0, closed = 0, complete = 0, consistent = 0, checked = 0, closed_consistent = 0, false_closure = 0,
                    predefined_constructed = 0, predefined_closed = 0;
        std::vector<double> seconds;
    };
    std::map<std::string, Tally> tally;
    for (const auto& r : rows) {
        Tally& t = tally[r.algorithm];
        if (!r.constructed) continue;
        ++t.constructed;
        t.seconds.push_back(r.seconds);
        if (r.built) ++t.built;
        if (r.closed) ++t.closed;
        if (r.closed || r.stop == "floor" || r.stop == "degenerate" || r.stop == "pure") ++t.complete;
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
        std::cout << format("%-13s constructed=%3d built=%3d closed=%3d closed+consistent=%3d FALSE-closures=%2d complete=%3d consistent=%3d/%3d "
                            "predefined closed=%3d/%3d median=%.4fs total=%.1fs\n",
                            kv.first.c_str(), static_cast<int>(t.constructed), static_cast<int>(t.built), static_cast<int>(t.closed),
                            static_cast<int>(t.closed_consistent), static_cast<int>(t.false_closure), static_cast<int>(t.complete),
                            static_cast<int>(t.consistent), static_cast<int>(t.checked), static_cast<int>(t.predefined_closed),
                            static_cast<int>(t.predefined_constructed), median, total);
    }

    std::cout << "\n--- FALSE closures (reported closed, but the stored dew points disagree with a blind flash) ---\n";
    for (const auto& r : rows) {
        if (r.closed && r.dev >= 1e-3) {
            std::cout << format("%-58s %-13s dev=%8.4f n=%4d pmax=%9.3f MPa\n", r.label.substr(0, 58).c_str(), r.algorithm.c_str(), r.dev,
                                static_cast<int>(r.n), r.pmax / 1e6);
        }
    }

    if (const char* csv = std::getenv("COOLPROP_PHASE_ENVELOPE_TORTURE_CSV")) {
        std::ofstream f(csv);
        f << "label,algorithm,predefined,constructed,built,closed,stop,n,icrit,seconds,pmax_Pa,Tmin_K,Tmax_K,dev,error\n";
        for (const auto& r : rows) {
            std::string err = r.error;
            std::replace(err.begin(), err.end(), ',', ';');
            std::replace(err.begin(), err.end(), '\n', ' ');
            f << '"' << r.label << "\"," << r.algorithm << ',' << r.predefined << ',' << r.constructed << ',' << r.built << ',' << r.closed << ','
              << r.stop << ',' << r.n << ',' << r.icrit << ',' << r.seconds << ',' << r.pmax << ',' << r.Tmin << ',' << r.Tmax << ',' << r.dev
              << ",\"" << err << "\"\n";
        }
        std::cout << "wrote " << csv << '\n';
    }

    // Known pre-existing defect of the DEFAULT algorithm, found by this corpus: the legacy
    // tracer stores NaN mole fractions for a 15-component gas with heavy traces, because its
    // 100 Pa dew start is a nearly pure n-octane liquid below n-octane's triple point.  Pinned
    // rather than ignored: the count may not grow, and no candidate algorithm may join it.
    // Tracked as a bd bug; delete this pin when the legacy tracer is fixed.
    std::size_t legacy_nonfinite = 0;
    for (const auto& r : rows) {
        if (r.constructed && !r.finite) {
            CAPTURE(r.label, r.algorithm);
            CHECK(r.algorithm == "legacy");
            if (r.algorithm == "legacy") {
                ++legacy_nonfinite;
            }
        }
    }
    CHECK(legacy_nonfinite <= 1);

    // Pins.  Legacy on the predefined subset: baseline 2026-09-11 was 116 constructed, 106 closed.
    const Tally& legacy = tally["legacy"];
    CHECK(legacy.predefined_constructed >= 116);
    CHECK(legacy.predefined_closed >= 106);
    // Quality pins measured 2026-09-11.  Closure alone is not enough: a false closure is a
    // silently wrong envelope, which is worse than an honest failure, so it is bounded too.
    CHECK(legacy.closed_consistent >= 123);
    CHECK(legacy.false_closure <= 7);
    CHECK(tally["lnK_density"].closed_consistent >= 116);
    CHECK(tally["lnK_density"].false_closure <= 15);
    // Every algorithm: finite stored values, CoolProp exceptions only, bounded point count.
    for (const auto& r : rows) {
        CAPTURE(r.label, r.algorithm, r.error);
        CHECK(r.coolprop_error);
        CHECK(r.n <= 1002);  // finalize may insert the two maxima points
    }
}

#endif  // ENABLE_CATCH
