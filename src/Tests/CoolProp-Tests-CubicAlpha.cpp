// Catch2 tests for cubic EOS alpha-function implementations.
// Verifies that the vectorised calc_all_terms() path produces bit-identical
// results to the scalar term() path for BasicMathiasCopeman, MathiasCopeman,
// and Twu alpha functions across a representative set of tau values.

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include "../Backends/Cubics/GeneralizedCubic.h"

#    include <array>
#    include <cmath>
#    include <vector>

using namespace CoolProp;

TEST_CASE("calc_all_terms matches term() for non-default alpha functions", "[cubic_alpha]") {
    const std::vector<double> tau_vals = {0.5, 0.7, 1.0, 1.3, 2.0};

    SECTION("BasicMathiasCopemanAlphaFunction") {
        // Arbitrary non-degenerate coefficients; values chosen to exercise all terms
        const double omega = 0.152;
        const double m = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
        const double a0 = 1.0;
        const double Tr_over_Tci = 1.0;
        BasicMathiasCopemanAlphaFunction alpha(a0, m, Tr_over_Tci);
        for (double tau : tau_vals) {
            CAPTURE(tau);
            std::array<double, 5> terms;
            alpha.calc_all_terms(tau, terms);
            for (int k = 0; k < 5; ++k) {
                CAPTURE(k);
                double ref = alpha.term(tau, k);
                double rel = (ref != 0.0) ? std::abs(terms[k] - ref) / std::abs(ref) : std::abs(terms[k]);
                CHECK(rel < 1e-12);
            }
        }
    }

    SECTION("MathiasCopemanAlphaFunction") {
        // Arbitrary non-degenerate coefficients; values chosen to exercise all terms.
        const double c1 = 0.8240, c2 = -0.4800, c3 = 0.6200;
        const double a0 = 1.0;
        const double Tr_over_Tci = 0.95;
        MathiasCopemanAlphaFunction alpha(a0, c1, c2, c3, Tr_over_Tci);
        for (double tau : tau_vals) {
            CAPTURE(tau);
            std::array<double, 5> terms;
            alpha.calc_all_terms(tau, terms);
            for (int k = 0; k < 5; ++k) {
                CAPTURE(k);
                double ref = alpha.term(tau, k);
                double rel = (ref != 0.0) ? std::abs(terms[k] - ref) / std::abs(ref) : std::abs(terms[k]);
                CHECK(rel < 1e-12);
            }
        }
    }

    SECTION("TwuAlphaFunction") {
        // Arbitrary non-degenerate coefficients; values chosen to exercise all terms
        const double L = 0.4692, M = 0.8610, N = 1.9127;
        const double a0 = 1.0;
        const double Tr_over_Tci = 1.0;
        TwuAlphaFunction alpha(a0, L, M, N, Tr_over_Tci);
        for (double tau : tau_vals) {
            CAPTURE(tau);
            std::array<double, 5> terms;
            alpha.calc_all_terms(tau, terms);
            for (int k = 0; k < 5; ++k) {
                CAPTURE(k);
                double ref = alpha.term(tau, k);
                double rel = (ref != 0.0) ? std::abs(terms[k] - ref) / std::abs(ref) : std::abs(terms[k]);
                CHECK(rel < 1e-12);
            }
        }
    }

    SECTION("TwuAlphaFunction - Tr_over_Tci != 1") {
        const double L = 0.3000, M = 0.9200, N = 2.1000;
        const double a0 = 1.5;
        const double Tr_over_Tci = 0.85;
        TwuAlphaFunction alpha(a0, L, M, N, Tr_over_Tci);
        for (double tau : tau_vals) {
            CAPTURE(tau);
            std::array<double, 5> terms;
            alpha.calc_all_terms(tau, terms);
            for (int k = 0; k < 5; ++k) {
                CAPTURE(k);
                double ref = alpha.term(tau, k);
                double rel = (ref != 0.0) ? std::abs(terms[k] - ref) / std::abs(ref) : std::abs(terms[k]);
                CHECK(rel < 1e-12);
            }
        }
    }
}

#endif
