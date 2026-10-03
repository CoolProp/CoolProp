// Tests of the Chebyshev tables of the pressure equation (src/Backends/Helmholtz/ChebDensitySolver.h).
//
//  - the regrouped model (no tables) reproduces CoolProp's own p and alphar;
//  - the per-piece margin really bounds |G_tables - G_true| on dense grids;
//  - all-roots parity: every sign change of the true G found by a dense scan lies in an interval reported by
//    ChebyshevBernstein::real_roots on the tables (with the margin as tolerance), and the true G changes sign
//    across every certified interval;
//  - declines: non-analytic and other non-factoring terms, cubic backends, invalid options, (T, x) outside the rectangle.
#if defined(ENABLE_CATCH)

#    include <algorithm>
#    include <cfloat>
#    include <cmath>
#    include <memory>
#    include <random>
#    include <string>
#    include <vector>

#    include <catch2/catch_all.hpp>

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/DataStructures.h"
#    include "../Backends/Helmholtz/ChebDensitySolver.h"
#    include "../Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#    include "CoolProp/numerics/ChebyshevBernstein.h"

using namespace CoolProp;
namespace CB = CoolProp::ChebyshevBernstein;
namespace CD = CoolProp::ChebDensity;

namespace {

// Lower bounds on the per-piece error/margin statistics of the margin test, ~10x below the measured ones (median
// 0.006-0.019, 10th percentile 0.003-0.008 over the cases below): an inflated margin (e.g. a roundoff term 1e6 times
// too large) fails them, while the worst-case check alone would still pass on the fit-error-dominated pieces.
constexpr double MEDIAN_MIN = 1e-3, P10_MIN = 3e-4;

struct Case
{
    std::string backend, fluids;
    std::vector<double> z;
};

const std::vector<Case>& cases() {
    static const std::vector<Case> c = {
      {"GERG2008", "Methane&Ethane", {0.6, 0.4}},
      {"GERG2008",
       "Methane&Nitrogen&CarbonDioxide&Ethane&Propane&IsoButane&n-Butane&Isopentane&n-Pentane&n-Hexane",
       {0.906724, 0.031284, 0.004676, 0.045279, 0.00828, 0.001037, 0.001563, 0.000321, 0.000443, 0.000393}},  // Amarillo
      {"GERG2008", "Methane&CarbonDioxide", {0.5, 0.5}},
      {"GERG2008", "Nitrogen&n-Butane", {0.3, 0.7}},
      {"HEOS", "Methane&Ethane&Propane&n-Butane", {0.7, 0.15, 0.1, 0.05}},
      {"HEOS", "Propane&n-Hexane", {0.4, 0.6}},
      {"HEOS", "Methane", {1.0}},
      {"GERG2008", "Nitrogen&Argon&Oxygen", {0.7812, 0.0092, 0.2096}},
    };
    return c;
}

struct Built
{
    std::shared_ptr<AbstractState> AS;
    HelmholtzEOSMixtureBackend* heos = nullptr;
    std::shared_ptr<const CD::Tables> tab;
    double Tmin = 0, Tcmax = 0;
};

// The library parameters: T >= max(60 K, 0.3 Tc,max), tau up to 1.02 Tc,max / Tmin, delta up to 4
Built build(const Case& c) {
    Built b;
    b.AS.reset(AbstractState::factory(c.backend, c.fluids));
    b.heos = dynamic_cast<HelmholtzEOSMixtureBackend*>(b.AS.get());
    REQUIRE(b.heos != nullptr);
    b.AS->set_mole_fractions(c.z);
    for (std::size_t i = 0; i < c.z.size(); ++i)
        b.Tcmax = std::max(b.Tcmax, b.AS->get_fluid_constant(i, iT_critical));
    b.Tmin = std::max(60.0, 0.3 * b.Tcmax);
    CD::BuildOptions opt;
    opt.tau_max = 1.02 * b.Tcmax / b.Tmin;
    opt.delta_max = 4.0;
    std::string why;
    b.tab = CD::Tables::build(*b.heos, opt, &why);
    INFO(c.backend << "::" << c.fluids << " declined: " << why);
    REQUIRE(b.tab != nullptr);
    return b;
}

std::vector<double> temperatures(const Built& b) {
    std::vector<double> T;
    for (double f : {1.0, 1.15, 1.4, 1.8, 2.3, 2.9})
        T.push_back(b.Tmin * f);
    for (double f : {0.97, 1.0, 1.03, 1.5, 2.5})
        T.push_back(b.Tcmax * f);
    std::sort(T.begin(), T.end());
    return T;
}

// a finer grid from Tmin up, where the terms cancel most and the margins are tightest
std::vector<double> temperatures_fine(const Built& b) {
    std::vector<double> T = temperatures(b);
    for (int k = 0; k < 16; ++k)
        T.push_back(b.Tmin * (1.03 + 0.11 * k));
    std::sort(T.begin(), T.end());
    return T;
}

// compositions: the case's own, plus n random ones (fixed seed; shape 0.5 puts weight near the edges)
std::vector<std::vector<double>> compositions(const Case& c, int n) {
    std::vector<std::vector<double>> xs = {c.z};
    if (c.z.size() > 1) {
        std::mt19937_64 g(42);
        std::gamma_distribution<double> G(0.5);
        for (int k = 0; k < n; ++k) {
            std::vector<double> x(c.z.size());
            double s = 0;
            for (double& v : x)
                s += (v = G(g));
            for (double& v : x)
                v /= s;
            xs.push_back(x);
        }
    }
    return xs;
}

}  // namespace

TEST_CASE("ChebDensity: regrouped model reproduces CoolProp p and alphar", "[cheb_density]") {
    for (const auto& c : cases()) {
        Built b = build(c);
        CD::Tables::State S;
        double worst_p = 0, worst_a = 0;
        int n = 0;
        // the EOS itself at every (T, rho), also inside the dome (where an unimposed update returns the two-phase state)
        b.AS->specify_phase(iphase_gas);
        for (double T : temperatures(b)) {
            REQUIRE(b.tab->assemble(T, c.z, S));
            for (double D : {1e-4, 0.01, 0.2, 0.7, 1.0, 1.6, 2.4, 3.2}) {
                const double rho = D * S.rhor;
                try {
                    b.AS->update(DmolarT_INPUTS, rho, T);
                } catch (...) {
                    continue;
                }
                const double p = b.AS->p();
                double G, dG, sc;
                b.tab->true_G(S, D, p * S.t_scale, G, dG, sc);
                worst_p = std::max(worst_p, std::abs(G) / sc);
                worst_a = std::max(worst_a, std::abs(b.tab->alphar(S, D) - b.AS->alphar()) / (1 + std::abs(b.AS->alphar())));
                ++n;
            }
        }
        INFO(c.backend << "::" << c.fluids << ": " << n << " states, worst |G|/scale " << worst_p << ", worst alphar " << worst_a);
        CHECK(n > 40);
        CHECK(worst_p < 1e-14);
        CHECK(worst_a < 1e-12);
    }
}

TEST_CASE("ChebDensity: margin bounds the table error on dense grids", "[cheb_density]") {
    for (const auto& c : cases()) {
        Built b = build(c);
        CD::Tables::State S;
        double worst = 0, worst_margin = 0, worst_G = 0;
        std::vector<double> piece_ratio;  // max error/margin on each (state, piece)
        long n = 0;
        for (const auto& x : compositions(c, 4))
            for (double T : temperatures_fine(b)) {
                if (!b.tab->assemble(T, x, S)) continue;  // a random x can push tau = Tr(x)/T outside the rectangle
                for (int p = 0; p < b.tab->n_pieces(); ++p) {
                    double pr = 0;
                    // dense: the error peaks between the interpolation nodes, in spots a few hundred points per piece miss
                    for (int k = 0; k <= 1000; ++k) {
                        const double u = -1 + 2.0 * k / 1000, D = b.tab->delta_of(p, u);
                        double G, dG, sc;
                        b.tab->true_G(S, D, 0, G, dG, sc);
                        const double e = std::abs(CB::clenshaw<CD::NG>(S.G[p], u) - G);
                        pr = std::max(pr, e / S.margin[p]);
                        worst_G = std::max(worst_G, std::abs(G));
                        ++n;
                    }
                    worst = std::max(worst, pr);
                    piece_ratio.push_back(pr);
                    worst_margin = std::max(worst_margin, S.margin[p]);
                }
            }
        std::sort(piece_ratio.begin(), piece_ratio.end());
        const double median = piece_ratio.empty() ? 0 : piece_ratio[piece_ratio.size() / 2];
        const double p10 = piece_ratio.empty() ? 0 : piece_ratio[piece_ratio.size() / 10];
        INFO(c.backend << "::" << c.fluids << ": per-piece error/margin median " << median << ", 10th percentile " << p10);
        INFO(c.backend << "::" << c.fluids << ": " << b.tab->n_pieces() << " pieces, " << b.tab->n_groups() << " groups of " << b.tab->n_terms()
                       << " terms; " << n << " points, worst error/margin " << worst << ", largest margin " << worst_margin << ", largest |G| "
                       << worst_G);
        CHECK(n > 1000);
        CHECK(worst <= 1.0);
        // the margin must stay useful, not merely valid: an inflated margin would pass the check above trivially.
        // Per piece, not only the worst one: a margin inflated where roundoff dominates must show up too.
        CHECK(worst > 0.01);
        CHECK(median > MEDIAN_MIN);
        CHECK(p10 > P10_MIN);
    }
}

TEST_CASE("ChebDensity: all-roots parity with a dense scan of the true equation", "[cheb_density]") {
    for (const auto& c : cases()) {
        Built b = build(c);
        CD::Tables::State S;
        const double dmax = b.tab->options().delta_max;
        const int NS = 20000;
        long n_true = 0, n_cert = 0, n_unres = 0, n_states = 0;
        std::vector<CB::Root> rr;
        for (const auto& x : compositions(c, 1))
            for (double T : temperatures(b)) {
                if (!b.tab->assemble(T, x, S)) {
                    CHECK(x != c.z);  // only a random x may fall outside the tau range
                    continue;
                }
                // F = delta Z on the scan grid (t = 0)
                std::vector<double> Dg(NS + 1), Fg(NS + 1);
                for (int k = 0; k <= NS; ++k) {
                    double G, dG, sc;
                    Dg[k] = dmax * k / NS;
                    b.tab->true_G(S, Dg[k], 0, G, dG, sc);
                    Fg[k] = G;
                }
                // pressures: F at a spread of densities (so the dome's three-root states are included), and fixed ones
                std::vector<double> ts;
                for (int k = 1; k <= 24; ++k)
                    ts.push_back(Fg[k * NS / 25]);
                for (double p : {1e3, 1e5, 1e6, 5e6, 2e7, 1e8})
                    ts.push_back(p * S.t_scale);
                for (double t : ts) {
                    if (!(t > 0)) continue;
                    ++n_states;
                    // reported intervals in delta, over all pieces
                    std::vector<std::pair<double, double>> rep;
                    for (int p = 0; p < b.tab->n_pieces(); ++p) {
                        CD::CoeffsG g = S.G[p];
                        g[0] -= t;
                        CB::real_roots<CD::NG>(g, S.margin[p], rr);
                        for (const auto& r : rr) {
                            const double Da = b.tab->delta_of(p, r.ua), Db = b.tab->delta_of(p, r.ub);
                            rep.emplace_back(Da, Db);
                            if (r.certified) {
                                ++n_cert;
                                // certified: the TRUE G changes sign across the interval
                                double Ga, Gb, dG, sc;
                                b.tab->true_G(S, Da, t, Ga, dG, sc);
                                b.tab->true_G(S, Db, t, Gb, dG, sc);
                                INFO(c.fluids << " T=" << T << " t=" << t << " certified [" << Da << ", " << Db << "]: G " << Ga << " " << Gb);
                                CHECK((Ga < 0) != (Gb < 0));
                            } else {
                                ++n_unres;
                            }
                        }
                    }
                    // every true root bracketed by the scan lies in a reported interval: bisect the true G down to the
                    // last bit, then require containment (a few ulp of slack for the rounding of the delta mapping)
                    for (int k = 0; k < NS; ++k) {
                        const double a = Fg[k] - t, bb = Fg[k + 1] - t;
                        if ((a < 0) == (bb < 0)) continue;
                        ++n_true;
                        double lo = Dg[k], hi = Dg[k + 1];
                        const bool neg_lo = a < 0;
                        for (int it = 0; it < 200; ++it) {
                            const double mid = 0.5 * (lo + hi);
                            if (!(lo < mid && mid < hi)) break;
                            double G, dG, sc;
                            b.tab->true_G(S, mid, t, G, dG, sc);
                            ((G < 0) == neg_lo ? lo : hi) = mid;
                        }
                        const double root = 0.5 * (lo + hi), slack = 8 * DBL_EPSILON * root;
                        const bool covered = std::any_of(rep.begin(), rep.end(), [&](const std::pair<double, double>& iv) {
                            return iv.first - slack <= root && root <= iv.second + slack;
                        });
                        INFO(c.fluids << " T=" << T << " t=" << t << ": true root " << root << " not in any reported interval");
                        CHECK(covered);
                    }
                }
            }
        INFO(c.backend << "::" << c.fluids << ": " << n_states << " (T, p), " << n_true << " true roots, " << n_cert << " certified, " << n_unres
                       << " unresolved");
        CHECK(n_true > n_states);  // some states have several roots
        // every root isolated and certified: no wide unresolved interval can make the containment check above vacuous
        CHECK(n_cert == n_true);
        CHECK(n_unres == 0);
    }
}

TEST_CASE("ChebDensity: declines", "[cheb_density]") {
    CD::BuildOptions opt;
    opt.tau_max = 3;
    SECTION("non-analytic terms") {
        for (const std::string f : {"Water", "CarbonDioxide", "CarbonDioxide&Water"}) {
            std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", f));
            std::string why;
            CHECK(CD::Tables::build(*dynamic_cast<HelmholtzEOSMixtureBackend*>(AS.get()), opt, &why) == nullptr);
            CHECK(why.find("non-analytic") != std::string::npos);
        }
    }
    SECTION("cubic backends") {
        for (const std::string f : {"Methane", "Methane&Ethane"}) {
            std::shared_ptr<AbstractState> AS(AbstractState::factory("PR", f));
            auto* heos = dynamic_cast<HelmholtzEOSMixtureBackend*>(AS.get());
            REQUIRE(heos != nullptr);
            AS->set_mole_fractions(f == "Methane" ? std::vector<double>{1.0} : std::vector<double>{0.5, 0.5});
            std::string why;
            CHECK(CD::Tables::build(*heos, opt, &why) == nullptr);
            INFO(f << ": " << why);
            CHECK_FALSE(why.empty());
            CHECK(why != "mole fractions not set");
        }
    }
    SECTION("other residual term types") {
        std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", "Ammonia"));  // Gao et al. B terms
        AS->set_mole_fractions({1.0});
        std::string why;
        CHECK(CD::Tables::build(*dynamic_cast<HelmholtzEOSMixtureBackend*>(AS.get()), opt, &why) == nullptr);
        CHECK(why.find("other than generalized exponential") != std::string::npos);
    }
    SECTION("invalid options throw") {
        std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", "Methane"));
        auto& heos = *dynamic_cast<HelmholtzEOSMixtureBackend*>(AS.get());
        AS->set_mole_fractions({1.0});
        CD::BuildOptions bad;
        CHECK_THROWS(CD::Tables::build(heos, bad));  // tau_max not set
        bad.tau_max = 2;
        bad.tau_min = 3;
        CHECK_THROWS(CD::Tables::build(heos, bad));
        bad.tau_min = 0;
        bad.min_width = 1e-15;  // would bisect the delta range into ulp-wide pieces
        CHECK_THROWS(CD::Tables::build(heos, bad));
        bad.min_width = 1e-3;
        bad.tol = NAN;
        CHECK_THROWS(CD::Tables::build(heos, bad));
        bad.tol = 1e-6;
        CHECK(CD::Tables::build(heos, bad) != nullptr);
    }
    SECTION("outside the rectangle") {
        std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", "Methane&Ethane"));
        auto& heos = *dynamic_cast<HelmholtzEOSMixtureBackend*>(AS.get());
        CD::BuildOptions o;
        o.tau_min = 0.5;
        o.tau_max = 2.0;
        std::string why;
        CHECK(CD::Tables::build(heos, o, &why) == nullptr);
        CHECK(why == "mole fractions not set");
        AS->set_mole_fractions({0.5, 0.5});
        auto tab = CD::Tables::build(heos, o);
        REQUIRE(tab != nullptr);
        CD::Tables::State S;
        const std::vector<double> x = {0.5, 0.5};
        const double Tr = heos.Reducing->Tr(x);
        CHECK(tab->assemble(Tr / 1.0, x, S));
        CHECK_FALSE(tab->assemble(Tr / 2.1, x, S));  // tau > tau_max
        CHECK_FALSE(tab->assemble(Tr / 0.4, x, S));  // tau < tau_min
        CHECK_FALSE(tab->assemble(-1, x, S));
        CHECK_FALSE(tab->assemble(Tr, {1.0}, S));  // wrong length
        CHECK_FALSE(tab->assemble(Tr, {NAN, 0.5}, S));
    }
}

#endif
