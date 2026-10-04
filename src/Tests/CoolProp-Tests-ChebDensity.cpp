// Tests of the Chebyshev tables of the pressure equation (src/Backends/Helmholtz/ChebDensitySolver.h).
//
//  - the regrouped model (no tables) reproduces CoolProp's own p and alphar;
//  - the per-piece margin really bounds |G_tables - G_true| on dense grids;
//  - all-roots parity: every sign change of the true G found by a dense scan lies in an interval reported by
//    ChebyshevBernstein::real_roots on the tables (with the margin as tolerance), and the true G changes sign
//    across every certified interval;
//  - the non-analytic add-in: closed-form evaluation vs CoolProp, the rigorous bound used to skip it, the fit path used;
//  - declines: non-factoring terms, cubic backends, invalid options, (T, x) outside the rectangle.
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
#    include "CoolProp/Configuration.h"
#    include "CoolProp/DataStructures.h"
#    include "../Backends/Helmholtz/ChebDensitySolver.h"
#    include "../Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#    include "CoolProp/numerics/ChebyshevBernstein.h"

using namespace CoolProp;
namespace CB = CoolProp::ChebyshevBernstein;
namespace CD = CoolProp::ChebDensity;

namespace {

// Lower bounds on the per-piece error/margin statistics of the margin test.  Measured (macOS arm64 clang) over the
// cases below: median 0.0058-0.019, 10th percentile 0.0027-0.0075, so ~6x and ~9x headroom for other toolchains
// (FMA contraction changes the roundoff).  An inflated margin (e.g. a roundoff term 1e6 times too large) fails them,
// while the worst-case check alone would still pass on the fit-error-dominated pieces.
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
      // non-analytic critical-region terms (IAPWS-95 water, Span-Wagner CO2)
      {"HEOS", "Water", {1.0}},
      {"HEOS", "CarbonDioxide", {1.0}},
      {"HEOS", "CarbonDioxide&Water", {0.9, 0.1}},
      {"HEOS", "Nitrogen&Oxygen&Argon&Water", {0.77, 0.2, 0.01, 0.02}},  // humid air
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
// temperatures packed around tau = Tr(x)/T = 1, where the non-analytic terms are largest and least smooth; the band
// 1e-7 <= |tau - 1| <= 1e-5 is where an ungraded fit next to delta = 1 lost roots (pre-PR review of COO-70)
std::vector<double> temperatures_critical(const Built& b, const std::vector<double>& x) {
    std::vector<double> T;
    const double Tr = b.heos->Reducing->Tr(x);
    for (double dt :
         {-0.1, -0.03, -0.01, -1e-3, -1e-5, -1e-6, -1.58e-7, -1.26e-7, -1e-7, 0.0, 1e-7, 1.26e-7, 1.58e-7, 1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.03, 0.1})
        T.push_back(Tr / (1 + dt));
    return T;
}

// plus temperatures around tau = 1
std::vector<double> temperatures_fine(const Built& b, const std::vector<double>& x) {
    std::vector<double> T = temperatures(b);
    for (int k = 0; k < 16; ++k)
        T.push_back(b.Tmin * (1.03 + 0.11 * k));
    for (double T_c : temperatures_critical(b, x))
        T.push_back(T_c);
    std::sort(T.begin(), T.end());
    return T;
}

// compositions: the case's own, plus n random ones.  Built from the raw mt19937_64 output (portable: the standard
// distributions are implementation-defined, so std::gamma_distribution gives different x on libstdc++, libc++ and
// MSVC).  Squared exponential variates put weight near the edges of the simplex.
std::vector<std::vector<double>> compositions(const Case& c, int n) {
    std::vector<std::vector<double>> xs = {c.z};
    if (c.z.size() > 1) {
        std::mt19937_64 g(42);
        for (int k = 0; k < n; ++k) {
            std::vector<double> x(c.z.size());
            double s = 0;
            for (double& v : x) {
                const double u = (static_cast<double>(g() >> 11U) + 0.5) * 0x1.0p-53, e = -std::log(u);
                s += (v = e * e);
            }
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
        // fewer random compositions where the non-analytic add-in makes every point expensive
        for (const auto& x : compositions(c, b.tab->has_nonanalytic() ? 2 : 4))
            for (double T : temperatures_fine(b, x)) {
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
        double max_unres_width = 0;  // widest uncertified interval, in delta
        const auto xs = compositions(c, 1);
        std::vector<long> n_assembled(xs.size(), 0);
        std::vector<CB::Root> rr;
        for (std::size_t ix = 0; ix < xs.size(); ++ix)
            for (double T : b.tab->has_nonanalytic() ? temperatures_fine(b, xs[ix]) : temperatures(b)) {
                const auto& x = xs[ix];
                if (!b.tab->assemble(T, x, S)) {
                    CHECK(ix != 0);  // only a random x may fall outside the tau range
                    continue;
                }
                ++n_assembled[ix];
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
                                max_unres_width = std::max(max_unres_width, Db - Da);
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
        // (>=: two roots closer than the scan step count once in n_true.)  Uncertified intervals are allowed only
        // rarely (e.g. a tangency, the critical point) and only narrow.
        CHECK(n_cert + n_unres >= n_true);
        CHECK(n_unres <= n_true / 100);
        INFO("widest uncertified interval " << max_unres_width << " in delta");
        // At a critical point G ~ c (delta - delta_c)^3, so a margin m localizes the root only to ~(m / c)^(1/3): ~1e-3
        // for pure CO2.  This cap does not stop a whole graded piece near delta = 1 (~1e-3 wide) from being reported
        // uncertified; the count cap above (<= 1 % of the roots) and the containment check are what bind there.
        CHECK(max_unres_width <= 3e-3);
        for (std::size_t ix = 0; ix < xs.size(); ++ix) {
            INFO("composition " << ix);
            CHECK(n_assembled[ix] > 0);  // each composition actually tested
        }
    }
}

TEST_CASE("ChebDensity: declines", "[cheb_density]") {
    CD::BuildOptions opt;
    opt.tau_max = 3;
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

TEST_CASE("ChebDensity: non-analytic terms", "[cheb_density]") {
    for (const std::string f : {"Water", "CarbonDioxide"}) {
        std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", f));
        auto* heos = dynamic_cast<HelmholtzEOSMixtureBackend*>(AS.get());
        REQUIRE(heos != nullptr);
        auto& na = heos->get_components()[0].EOS().alphar.NonAnalytic;
        REQUIRE(na.N > 0);
        std::vector<CD::NonAnalyticTerm> terms;
        terms.reserve(na.elements.size());
        for (const auto& el : na.elements)
            terms.push_back({el.n, el.a, el.b, el.beta, el.A, el.B, el.C, el.D});

        SECTION(f + ": closed form agrees with ResidualHelmholtzNonAnalytic") {
            double worst = 0;
            for (double tau : {0.5, 0.9, 0.99, 0.9999, 1.0, 1.0001, 1.01, 1.2, 2.5})
                for (double D : {0.05, 0.5, 0.9, 0.99, 0.9999, 1.0, 1.0001, 1.01, 1.3, 2.5, 3.9}) {
                    HelmholtzDerivatives d;
                    na.all_deltaonly(tau, D, d);
                    const auto v = CD::eval_nonanalytic(terms, tau, D);
                    const double chi = D * d.dalphar_ddelta, dchi = d.dalphar_ddelta + D * d.d2alphar_ddelta2;
                    INFO(f << " tau=" << tau << " delta=" << D << ": alphar " << v.alphar << " vs " << d.alphar << ", chi " << v.chi << " vs " << chi
                           << ", dchi " << v.dchi << " vs " << dchi);
                    // the same expressions; differences are FMA contraction only, measured against the roundoff scale
                    CHECK(std::abs(v.alphar - d.alphar) <= 1e-13 * (std::abs(d.alphar) + 1e-300) + 1e-14 * v.parts);
                    CHECK(std::abs(v.chi - chi) <= 1e-14 * v.parts + 1e-300);
                    CHECK(std::abs(v.dchi - dchi) <= 1e-12 * std::abs(dchi) + 1e-14 * v.parts / std::max(std::abs(D - 1), 1e-4));
                    worst = std::max(worst, std::abs(v.chi - chi) / (v.parts + 1e-300));
                }
            INFO(f << ": worst |chi diff| / parts " << worst);
            CHECK(worst < 1e-14);
        }

        SECTION(f + ": the bound used to skip the fit holds") {
            // dense in delta on the table pieces, over tau in [0.3, 3.4], per term
            CD::BuildOptions o;
            o.tau_min = 0.3;
            o.tau_max = 3.4;
            AS->set_mole_fractions({1.0});
            auto tab = CD::Tables::build(*heos, o);
            REQUIRE(tab != nullptr);
            const double dtau_max = 2.4;
            // both over the whole tau range and at the actual |1 - tau| (assemble() takes the first rung of a ladder at or
            // above it; the bound is monotone in |1 - tau|, so the actual value is the strictest test), with tau packed
            // toward 1 where the bound is tightest
            std::vector<double> taus;
            for (int it = 0; it <= 40; ++it)
                taus.push_back(0.3 + 3.1 * it / 40.0);
            for (int k = 1; k <= 9; ++k) {
                taus.push_back(1 + std::pow(10.0, -k));
                taus.push_back(1 - std::pow(10.0, -k));
            }
            double worst = 0, worst_actual = 0;
            long n = 0;
            for (int p = 0; p < tab->n_pieces(); ++p)
                for (const auto& term : terms) {
                    const double K = term.bound_factor(tab->edges()[p], tab->edges()[p + 1], dtau_max);
                    REQUIRE(std::isfinite(K));
                    for (double tau : taus) {
                        const double Ka = term.bound_factor(tab->edges()[p], tab->edges()[p + 1], std::abs(tau - 1));
                        for (int j = 0; j <= 400; ++j) {
                            const double D = tab->delta_of(p, -1 + 2.0 * j / 400);
                            const auto v = CD::eval_nonanalytic({term}, tau, D);
                            const double ef = std::exp(-term.D * (tau - 1) * (tau - 1));
                            worst = std::max(worst, std::abs(v.chi) / (K * ef));
                            worst_actual = std::max(worst_actual, std::abs(v.chi) / (Ka * ef));
                            ++n;
                        }
                    }
                }
            INFO(f << ": " << n << " points, worst |chi| / bound " << worst << ", at the actual |1 - tau| " << worst_actual);
            CHECK(worst <= 1.0);
            CHECK(worst_actual <= 1.0);
            CHECK(worst_actual > 0.1);  // not vacuous: the bound is tight near tau = 1
        }

        SECTION(f + ": the add-in paths: 2-D table, fit at tau, skipped") {
            CD::BuildOptions o;
            o.tau_min = 0.3;
            o.tau_max = 3.4;
            AS->set_mole_fractions({1.0});
            auto tab = CD::Tables::build(*heos, o);
            REQUIRE(tab != nullptr);
            CD::Tables::State S;
            const double Tc = AS->get_fluid_constant(0, iT_reducing);
            REQUIRE(tab->assemble(Tc, {1.0}, S));  // tau = 1: every non-negligible piece from the 2-D table
            CHECK(S.na_table > 0);
            CHECK(S.na_fits == 0);
            // tau = 1.01: the valley Delta ~ 0 at tau - 1 = A |delta - 1|^(1/beta) crosses a few pieces, whose cells keep
            // no table: fit at this tau
            REQUIRE(tab->assemble(Tc / 1.01, {1.0}, S));
            CHECK(S.na_table > 0);
            CHECK(S.na_fits > 0);
            CHECK(S.na_fits < S.na_table);
            CHECK(S.na_fits <= 8);                       // a few (2 today); the exact count depends on rounding near the acceptance
            REQUIRE(tab->assemble(Tc / 3.0, {1.0}, S));  // tau = 3: exp(-D (tau-1)^2) is negligible
            CHECK(S.na_fits == 0);
            CHECK(S.na_table == 0);
            // delta = 1 is a piece edge
            CHECK(std::find(tab->edges().begin(), tab->edges().end(), 1.0) != tab->edges().end());
        }
    }
}

// COO-125: the spike's cached solver kept a raw pointer to the non-analytic terms of the backend that built it -- a
// heap use-after-free once that backend was destroyed, whose garbage made phase verdicts vary with ASLR.  The tables
// own copies of everything; this use after the backend's destruction is what the ASan CI job checks.
TEST_CASE("ChebDensity: tables outlive the backend that built them", "[cheb_density]") {
    std::shared_ptr<const CD::Tables> tab;
    std::vector<double> ref_G;
    const std::vector<double> x = {0.9, 0.1};
    const double T = 310;
    {
        std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", "CarbonDioxide&Water"));
        AS->set_mole_fractions(x);
        CD::BuildOptions o;
        o.tau_min = 0.3;
        o.tau_max = 3.4;
        tab = CD::Tables::build(*dynamic_cast<HelmholtzEOSMixtureBackend*>(AS.get()), o);
        REQUIRE(tab != nullptr);
        CD::Tables::State S;
        REQUIRE(tab->assemble(T, x, S));
        for (double D : {0.1, 0.9, 1.1, 2.0}) {
            double G, dG, sc;
            tab->true_G(S, D, 0, G, dG, sc);
            ref_G.push_back(G);
        }
    }  // backend destroyed here
    std::shared_ptr<AbstractState> other(AbstractState::factory("HEOS", "Water&CarbonDioxide"));  // reuse the freed heap
    CD::Tables::State S;
    REQUIRE(tab->assemble(T, x, S));
    CHECK(S.na_fits + S.na_table > 0);
    std::size_t k = 0;
    for (double D : {0.1, 0.9, 1.1, 2.0}) {
        double G, dG, sc;
        tab->true_G(S, D, 0, G, dG, sc);
        CHECK(G == ref_G[k++]);  // same object, same inputs: bit-identical
    }
}

TEST_CASE("ChebDensity: non-analytic tables are shared across mixtures", "[cheb_density]") {
    // the 2-D tables depend only on a fluid's non-analytic terms and the tolerance: built once per process, shared
    auto make = [](const std::string& f, const std::vector<double>& x, double tol) {
        std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", f));
        AS->set_mole_fractions(x);
        CD::BuildOptions o;
        o.tau_min = 0.3;
        o.tau_max = 3.4;
        o.tol = tol;
        auto t = CD::Tables::build(*dynamic_cast<HelmholtzEOSMixtureBackend*>(AS.get()), o);
        REQUIRE(t != nullptr);
        return t;
    };
    const auto co2w = make("CarbonDioxide&Water", {0.9, 0.1}, 1e-6);                      // non-analytic: CO2 (0), water (1)
    const auto air = make("Nitrogen&Oxygen&Argon&Water", {0.77, 0.2, 0.01, 0.02}, 1e-6);  // water (0)
    const auto co2 = make("CarbonDioxide", {1.0}, 1e-6);
    REQUIRE(co2w->nonanalytic_table(0) != nullptr);
    REQUIRE(co2w->nonanalytic_table(1) != nullptr);
    CHECK(co2w->nonanalytic_table(0) != co2w->nonanalytic_table(1));
    CHECK(air->nonanalytic_table(0) == co2w->nonanalytic_table(1));  // water's
    CHECK(co2->nonanalytic_table(0) == co2w->nonanalytic_table(0));  // CO2's
    CHECK(air->nonanalytic_table(1) == nullptr);
    // a different tolerance is a different table
    const auto co2_loose = make("CarbonDioxide", {1.0}, 1e-4);
    CHECK(co2_loose->nonanalytic_table(0) != co2->nonanalytic_table(0));
    // and the shared table gives the same tabulated G in each: CO2 at the same (T, x) through both
    CD::Tables::State S1, S2;
    REQUIRE(co2->assemble(300.0, {1.0}, S1));
    REQUIRE(co2w->assemble(300.0, {1.0, 0.0}, S2));
    CHECK(S1.na_table > 0);
    CHECK(S2.na_table > 0);
    auto table_G = [](const CD::Tables& t, const CD::Tables::State& S, double D, double& margin) {
        const auto& e = t.edges();
        const int p = static_cast<int>(std::upper_bound(e.begin(), e.end(), D) - e.begin()) - 1;
        margin = S.margin[p];
        return CB::clenshaw<CD::NG>(S.G[p], 2 * (D - e[p]) / (e[p + 1] - e[p]) - 1);
    };
    for (double D : {0.6, 0.97, 1.03, 1.2, 1.7}) {
        double m1 = 0, m2 = 0;
        const double G1 = table_G(*co2, S1, D, m1), G2 = table_G(*co2w, S2, D, m2);
        INFO("delta " << D << ": " << G1 << " vs " << G2 << ", margins " << m1 << " " << m2);
        CHECK(std::abs(G1 - G2) <= m1 + m2);
        CHECK(m1 < 1e-6);  // not vacuous: small, finite margins (the GE fits put ~2e-7 here)
        CHECK(m2 < 1e-6);
    }
}

namespace {
// sets a configuration bool for the scope, restoring the previous value
struct ConfigBoolScope
{
    configuration_keys key;
    bool old;
    ConfigBoolScope(configuration_keys k, bool v) : key(k), old(get_config_bool(k)) {
        set_config_bool(k, v);
    }
    ~ConfigBoolScope() {
        set_config_bool(key, old);
    }
    ConfigBoolScope(const ConfigBoolScope&) = delete;
    ConfigBoolScope& operator=(const ConfigBoolScope&) = delete;
};
}  // namespace

TEST_CASE("ChebDensity: PT flash with CHEBYSHEV_DENSITY_SOLVER", "[cheb_density][flash]") {
    struct C
    {
        std::string backend, fluids;
        std::vector<double> z;
        double Tlo, Thi;
    };
    const std::vector<C> cs = {
      {"GERG2008", "Methane&Ethane", {0.5, 0.5}, 150, 400},
      {"GERG2008", "Methane&Nitrogen&CarbonDioxide&Ethane&Propane", {0.85, 0.05, 0.03, 0.05, 0.02}, 150, 400},
      {"HEOS", "Methane&Ethane&Propane", {0.7, 0.2, 0.1}, 150, 400},
      {"HEOS", "CarbonDioxide&Water", {0.95, 0.05}, 280, 450},
    };
    for (const auto& c : cs) {
        std::shared_ptr<AbstractState> off(AbstractState::factory(c.backend, c.fluids)), on(AbstractState::factory(c.backend, c.fluids));
        off->set_mole_fractions(c.z);
        on->set_mole_fractions(c.z);
        std::mt19937_64 g(7);
        int n = 0, same = 0, single = 0, answered = 0;
        for (int k = 0; k < 150; ++k) {
            const double u1 = (static_cast<double>(g() >> 11U) + 0.5) * 0x1.0p-53, u2 = (static_cast<double>(g() >> 11U) + 0.5) * 0x1.0p-53;
            const double T = c.Tlo + (c.Thi - c.Tlo) * u1, p = std::exp(std::log(1e4) + (std::log(3e7) - std::log(1e4)) * u2);
            double rho_off = NAN, Q_off = NAN, rho_on = NAN, Q_on = NAN;
            try {
                off->update(PT_INPUTS, p, T);
                rho_off = off->rhomolar();
                Q_off = off->Q();
            } catch (...) {
            }
            {
                ConfigBoolScope flag(CHEBYSHEV_DENSITY_SOLVER, true);
                try {
                    on->update(PT_INPUTS, p, T);
                    rho_on = on->rhomolar();
                    Q_on = on->Q();
                } catch (...) {
                }
                auto* heos = dynamic_cast<HelmholtzEOSMixtureBackend*>(on.get());
                REQUIRE(heos != nullptr);
                if (heos->solver_rho_Tp_cheb(T, p) > 0) ++answered;
            }
            if (!std::isfinite(rho_off) || !std::isfinite(rho_on)) continue;
            ++n;
            const bool one_phase = !(Q_off > 0 && Q_off < 1) && !(Q_on > 0 && Q_on < 1);
            if (one_phase) {
                ++single;
                // both single-phase: the same stable root (the legacy solver's own convergence is ~1e-8)
                if (std::abs(rho_on / rho_off - 1) < 1e-6) ++same;
            }
        }
        INFO(c.backend << "::" << c.fluids << ": " << n << " states, " << single << " single-phase, " << same << " same density, kernel answered "
                       << answered);
        CHECK(n > 120);
        CHECK(answered > n / 2);    // the solver really serves the flash
        CHECK(same >= single - 2);  // at most a couple of states where the legacy solver took another root
    }
}

TEST_CASE("ChebDensity: tables follow changed interaction parameters", "[cheb_density][flash]") {
    ConfigBoolScope flag(CHEBYSHEV_DENSITY_SOLVER, true);
    const std::vector<double> z = {0.5, 0.5};
    std::shared_ptr<AbstractState> A(AbstractState::factory("HEOS", "Methane&Ethane")), B(AbstractState::factory("HEOS", "Methane&Ethane"));
    A->set_mole_fractions(z);
    B->set_mole_fractions(z);
    const double T = 250, p = 5e6;
    auto* ha = dynamic_cast<HelmholtzEOSMixtureBackend*>(A.get());
    auto* hb = dynamic_cast<HelmholtzEOSMixtureBackend*>(B.get());
    const double ra = ha->solver_rho_Tp_cheb(T, p);
    REQUIRE(ra > 0);
    // change B's reducing function; A keeps the original model
    B->set_binary_interaction_double(0, 1, "betaT", 1.05 * B->get_binary_interaction_double(0, 1, "betaT"));
    const double rb = hb->solver_rho_Tp_cheb(T, p);
    REQUIRE(rb > 0);
    CHECK(std::abs(rb / ra - 1) > 1e-4);  // a different model, different tables
    // and B's root is a root of B's model, and A's of A's: the pressure at the returned density
    auto p_at = [&](AbstractState& S, double rho) {
        S.specify_phase(iphase_gas);
        S.update(DmolarT_INPUTS, rho, T);
        const double pp = S.p();
        S.unspecify_phase();
        return pp;
    };
    CHECK(std::abs(p_at(*B, rb) / p - 1) < 1e-9);
    CHECK(std::abs(p_at(*A, ra) / p - 1) < 1e-9);
    CHECK(std::abs(p_at(*B, ra) / p - 1) > 1e-4);                    // A's root is not a root of B's model
    CHECK(std::abs(ha->solver_rho_Tp_cheb(T, p) / ra - 1) < 1e-14);  // A unchanged
}

TEST_CASE("ChebDensity: cubic backends are not served", "[cheb_density][flash]") {
    ConfigBoolScope flag(CHEBYSHEV_DENSITY_SOLVER, true);
    std::shared_ptr<AbstractState> A(AbstractState::factory("PR", "Methane&Ethane"));
    A->set_mole_fractions({0.5, 0.5});
    auto* h = dynamic_cast<HelmholtzEOSMixtureBackend*>(A.get());
    REQUIRE(h != nullptr);
    CHECK(h->solver_rho_Tp_cheb(250, 1e6) < 0);
    CHECK_NOTHROW(A->update(PT_INPUTS, 1e6, 250));
}

#endif
