// Catch2 tests for the runtime transport-property expression DSL (tag [expression]).
//
// Split out of CoolProp-Tests.cpp so the DSL's lexer/parser/binder/evaluator
// coverage, its golden regression against the hardcoded Tier-A routines, and the
// end-to-end JSON-load round-trip all live in one translation unit.

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include <chrono>
#    include <cmath>
#    include <cstddef>
#    include <map>
#    include <memory>
#    include <string>
#    include <vector>

#    include "CoolProp/CoolProp.h"
#    include "CoolProp/CoolPropFluid.h"
#    include "CoolProp/DataStructures.h"
#    include "CoolProp/Exceptions.h"
#    include "CoolProp/detail/json.h"
#    include "CoolProp/expression/Expression.h"
#    include "CoolProp/expression/ExpressionBlock.h"
#    include "CoolProp/expression/ExpressionCorrelation.h"
#    include "CoolProp/expression/detail/Lexer.h"
#    include "CoolProp/numerics/numerics.h"  // ValidNumber
#    include "CoolProp/AbstractState.h"
#    include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#    include "Backends/Helmholtz/TransportRoutines.h"

TEST_CASE("expression module links and evaluates a constant", "[expression]") {
    using namespace CoolProp::expression;
    Program prog = compile("3 + 4", {}, {});
    CHECK(prog.evaluate({}) == Catch::Approx(7.0));
}

TEST_CASE("lexer tokenizes numbers, idents, operators", "[expression]") {
    using namespace CoolProp::expression::detail;
    std::vector<Token> t = lex("let x = a[i]*T^2.66958e-08 + sum(i: b)");
    REQUIRE(t.front().type == TokenType::Keyword_let);
    bool sawScientific = false, sawCaret = false, sawColon = false, sawLBracket = false;
    for (const auto& tok : t) {
        if (tok.type == TokenType::Number && tok.number == Catch::Approx(2.66958e-08)) sawScientific = true;
        if (tok.type == TokenType::Caret) sawCaret = true;
        if (tok.type == TokenType::Colon) sawColon = true;
        if (tok.type == TokenType::LBracket) sawLBracket = true;
    }
    CHECK(sawScientific);
    CHECK(sawCaret);
    CHECK(sawColon);
    CHECK(sawLBracket);
    CHECK(t.back().type == TokenType::End);
}

TEST_CASE("lexer rejects an illegal character", "[expression]") {
    using namespace CoolProp::expression::detail;
    CHECK_THROWS_AS(lex("a @ b"), CoolProp::ValueError);
}

TEST_CASE("compile evaluates arithmetic, precedence, right-assoc power", "[expression]") {
    using namespace CoolProp::expression;
    CHECK(compile("2 + 3*4", {}, {}).evaluate({}) == Catch::Approx(14.0));
    CHECK(compile("2^3^2", {}, {}).evaluate({}) == Catch::Approx(512.0));
    CHECK(compile("-2^2", {}, {}).evaluate({}) == Catch::Approx(-4.0));
}

TEST_CASE("compile resolves constants and let bindings", "[expression]") {
    using namespace CoolProp::expression;
    Program p = compile("let y = a*2\ny + 1", {{"a", 5.0}}, {});
    CHECK(p.evaluate({}) == Catch::Approx(11.0));
}

TEST_CASE("compile sum over co-indexed arrays", "[expression]") {
    using namespace CoolProp::expression;
    Program p = compile("sum(i: a[i]*b[i])", {}, {{"a", {1, 2, 3}}, {"b", {4, 5, 6}}});
    CHECK(p.evaluate({}) == Catch::Approx(32.0));
}

TEST_CASE("compile binds inputs in required order", "[expression]") {
    using namespace CoolProp::expression;
    Program p = compile("T + rhomolar", {}, {});
    REQUIRE(p.requiredInputs().size() == 2);
    std::vector<double> iv(2);
    for (std::size_t k = 0; k < 2; ++k)
        iv[k] = (p.requiredInputs()[k] == CoolProp::iT) ? 300.0 : 50.0;
    CHECK(p.evaluate(iv) == Catch::Approx(350.0));
}

TEST_CASE("compile errors: unknown var, sum length mismatch, bad arity", "[expression]") {
    using namespace CoolProp::expression;
    CHECK_THROWS_AS(compile("nope + 1", {}, {}), CoolProp::ValueError);
    CHECK_THROWS_AS(compile("sum(i: a[i]*b[i])", {}, {{"a", {1, 2}}, {"b", {1, 2, 3}}}), CoolProp::ValueError);
    CHECK_THROWS_AS(compile("exp(1, 2)", {}, {}), CoolProp::ValueError);
    CHECK_THROWS_AS(compile("a[i]", {}, {{"a", {1.0}}}), CoolProp::ValueError);
}

TEST_CASE("DSL scientific-notation evaluates; pow two-arg; trig", "[expression]") {
    using namespace CoolProp::expression;
    CHECK(compile("1e3 + 1", {}, {}).evaluate({}) == Catch::Approx(1001.0));
    CHECK(compile("pow(2, 10)", {}, {}).evaluate({}) == Catch::Approx(1024.0));
    CHECK(compile("sqrt(2)^2", {}, {}).evaluate({}) == Catch::Approx(2.0));
}

TEST_CASE("DSL let chaining: a let sees earlier lets, not itself", "[expression]") {
    using namespace CoolProp::expression;
    CHECK(compile("let x = a*2\nlet y = x+1\ny", {{"a", 5.0}}, {}).evaluate({}) == Catch::Approx(11.0));
    CHECK_THROWS_AS(compile("let x = x+1\nx", {}, {}), CoolProp::ValueError);  // self-ref unknown
}

TEST_CASE("DSL nested summation is rejected in v1", "[expression]") {
    using namespace CoolProp::expression;
    CHECK_THROWS_AS(compile("sum(i: sum(j: a[i]*b[j]))", {}, {{"a", {1, 2}}, {"b", {3, 4}}}), CoolProp::ValueError);
}

TEST_CASE("DSL pressure and temperature share one input bucket", "[expression]") {
    using namespace CoolProp::expression;
    Program p = compile("p*2 + T", {}, {});
    // p costs an EOS call and T does not, but both are plain CoolProp::parameters
    // keys to the program: one vector, reported in first-reference order.
    const std::vector<CoolProp::parameters>& req = p.requiredInputs();
    REQUIRE(req.size() == 2);
    CHECK(req[0] == CoolProp::iP);
    CHECK(req[1] == CoolProp::iT);
    std::vector<double> iv(2);
    for (std::size_t k = 0; k < 2; ++k)
        iv[k] = (req[k] == CoolProp::iP) ? 1e5 : 300.0;
    CHECK(p.evaluate(iv) == Catch::Approx(2e5 + 300.0));
}

TEST_CASE("expression input names map onto CoolProp::parameters", "[expression]") {
    using namespace CoolProp::expression;
    // Every DSL spelling must round-trip through the shared parameters enum, and
    // the state variables must keep their CoolProp member spellings (NOT the
    // PropsSI shorthands) so existing fluid JSON keeps compiling.
    for (const auto& e : inputTable())
        CHECK(inputName(e.second) == e.first);
    CHECK(compile("T", {}, {}).requiredInputs()[0] == CoolProp::iT);
    CHECK(compile("rhomolar", {}, {}).requiredInputs()[0] == CoolProp::iDmolar);
    CHECK(compile("rhomass", {}, {}).requiredInputs()[0] == CoolProp::iDmass);
    CHECK(compile("molar_mass", {}, {}).requiredInputs()[0] == CoolProp::imolar_mass);
    CHECK(compile("p", {}, {}).requiredInputs()[0] == CoolProp::iP);
    CHECK_THROWS_AS(inputName(CoolProp::iHmolar), CoolProp::ValueError);
}

TEST_CASE("only allowlisted parameter names resolve", "[expression]") {
    using namespace CoolProp::expression;
    // `V` and `L` are valid get_parameter_index() names but MUST NOT resolve: a
    // transport formula naming its own output would recurse through keyed_output().
    CHECK_THROWS_AS(compile("V", {}, {}), CoolProp::ValueError);
    CHECK_THROWS_AS(compile("L", {}, {}), CoolProp::ValueError);
    CHECK_THROWS_AS(compile("Hmolar", {}, {}), CoolProp::ValueError);
    // ...unless the fluid JSON declares it as a constant, which is the normal
    // way an author introduces an arbitrary name.
    CHECK(compile("V", {{"V", 4.0}}, {}).evaluate({}) == Catch::Approx(4.0));
}

TEST_CASE("a JSON constant that collides with an input is rejected", "[expression]") {
    using namespace CoolProp::expression;
    // Inputs resolve before constants, so a constant named `T` would be dead data
    // and the formula would quietly mean the state temperature.  Silently
    // discarding what the author wrote is the wrong answer: reject it at compile
    // time, for every name in the table.
    for (const auto& e : inputTable()) {
        CAPTURE(e.first);
        CHECK_THROWS_AS(compile("1", {{e.first, 1.0}}, {}), CoolProp::ValueError);
    }
    // A constant whose name is merely close to an input is fine.
    CHECK(compile("T_reduce", {{"T_reduce", 132.0}}, {}).evaluate({}) == Catch::Approx(132.0));
}

TEST_CASE("evaluate rejects a wrong-sized input vector", "[expression]") {
    using namespace CoolProp::expression;
    Program p = compile("T + p", {}, {});
    REQUIRE(p.requiredInputs().size() == 2);
    CHECK_THROWS_AS(p.evaluate({}), CoolProp::ValueError);
    CHECK_THROWS_AS(p.evaluate({300.0}), CoolProp::ValueError);
    CHECK_THROWS_AS(p.evaluate({300.0, 1e5, 7.0}), CoolProp::ValueError);
    CHECK(p.evaluate({1e5, 300.0}) == Catch::Approx(1e5 + 300.0));
}

TEST_CASE("a default-constructed Program throws instead of dereferencing null", "[expression]") {
    using namespace CoolProp::expression;
    Program p;  // no compile() ever ran
    CHECK_THROWS_AS(p.evaluate({}), CoolProp::ValueError);
    CHECK_THROWS_AS(p.requiredInputs(), CoolProp::ValueError);
}

TEST_CASE("DSL inputs rhomass and molar_mass resolve", "[expression]") {
    using namespace CoolProp::expression;
    Program p = compile("rhomass * molar_mass", {}, {});
    REQUIRE(p.requiredInputs().size() == 2);
    std::vector<double> iv(2);
    for (std::size_t k = 0; k < 2; ++k)
        iv[k] = (p.requiredInputs()[k] == CoolProp::iDmass) ? 2.0 : 3.0;
    CHECK(p.evaluate(iv) == Catch::Approx(6.0));
}

TEST_CASE("ExpressionData default-constructs unset", "[expression]") {
    CoolProp::ExpressionData d;
    CHECK(!d.correlation);
}

TEST_CASE("expression block compile path (constants+arrays) yields expected value", "[expression]") {
    using namespace CoolProp::expression;
    std::map<std::string, double> consts{{"T_reduce", 132.0}, {"rhomolar_reduce", 10000.0}};
    std::map<std::string, std::vector<double>> arrays{{"a", {1.0e-5}}, {"d1", {1.0}}, {"t1", {0.2}}};
    Program p = compile("let delta = rhomolar/rhomolar_reduce\nlet tau = T_reduce/T\nsum(i: a[i]*delta^d1[i]*tau^t1[i])", consts, arrays);
    std::vector<double> iv(p.requiredInputs().size());
    for (std::size_t k = 0; k < iv.size(); ++k) {
        const CoolProp::parameters key = p.requiredInputs()[k];
        iv[k] = (key == CoolProp::iT) ? 300.0 : (key == CoolProp::iDmolar) ? 5000.0 : 0.0;
    }
    double expected = 1.0e-5 * std::pow(0.5, 1.0) * std::pow(132.0 / 300.0, 0.2);
    CHECK(p.evaluate(iv) == Catch::Approx(expected));
}

// ---------------------------------------------------------------------------
// Task 9: GOLDEN REGRESSION — DSL reproduces the hardcoded Tier-A transport
// correlations.  For each Tier-A form we translate the exact C++ algebra in
// TransportRoutines.cpp into a DSL formula built from the SAME coefficients
// read from the fluid's transport struct, evaluate both at a grid of (T,rho),
// and CHECK they agree.  This is the feature-completeness proof.
//
// Tolerances:
//   * EVERY form matches to <1e-14 relative (the real completeness proof) at
//     every grid point — asserted on all 8 forms.
//   * Three forms (powers_of_Tr, collision_integral, polynomial_and_exponential)
//     additionally match BIT-EXACTLY and are asserted with `got == expected`.
//   * The other five (powers_of_T, modified_Batschinski_Hildebrand,
//     ratio_of_polynomials, eta0_and_poly, residual polynomial) match to ~1-3
//     ULP (observed max relative error ~2.4e-15) but NOT bit-for-bit.  The DSL
//     and the routine perform the identical op sequence; the last-ULP
//     divergence is FMA/`-ffp-contract`-class rounding (the optimized routine
//     may fuse `a[i]*pow(...)`+accumulate; the DSL evaluator rounds each op
//     separately).  This is a compiler codegen artifact, NOT an algebra
//     difference, so we keep the 1e-14 assertion (which they pass with margin)
//     rather than loosening it — and do NOT claim bit-exactness where it is
//     not achieved.
//
// TEST-ONLY: no production code is touched.
// ---------------------------------------------------------------------------

namespace {
std::shared_ptr<CoolProp::HelmholtzEOSMixtureBackend> make_HEOS_for(const std::string& fluid) {
    std::vector<std::string> names(1, fluid);
    return std::make_shared<CoolProp::HelmholtzEOSMixtureBackend>(names);
}
// Fill the input array for a Program in its requiredInputs() order -- the same
// one-liner the production host (ExpressionCorrelation::eval) uses.
void fill_inputs(const CoolProp::expression::Program& p, CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::vector<double>& iv) {
    const std::vector<CoolProp::parameters>& req = p.requiredInputs();
    iv.assign(req.size(), 0.0);
    for (std::size_t k = 0; k < req.size(); ++k)
        iv[k] = HEOS.keyed_output(req[k]);
}
}  // namespace

TEST_CASE("golden: viscosity dilute powers_of_T", "[expression][golden]") {
    using namespace CoolProp::expression;
    auto HEOS = make_HEOS_for("R123");  // viscosity/dilute type == powers_of_T
    auto& data = HEOS->get_components()[0].transport.viscosity_dilute.powers_of_T;
    REQUIRE(!data.a.empty());
    // C++: summer += a[i]*pow(T, t[i]);  return summer;
    Program p = compile("sum(i: a[i]*T^t[i])", {},
                        {{"a", std::vector<double>(data.a.begin(), data.a.end())}, {"t", std::vector<double>(data.t.begin(), data.t.end())}});
    int checks = 0;
    for (double T : {250.0, 300.0, 400.0, 500.0}) {
        for (double rho : {0.1, 100.0, 5000.0}) {
            HEOS->update(CoolProp::DmolarT_INPUTS, rho, T);
            double expected = CoolProp::TransportRoutines::viscosity_dilute_powers_of_T(*HEOS);
            std::vector<double> iv;
            fill_inputs(p, *HEOS, iv);
            double got = p.evaluate(iv);
            CAPTURE(T, rho, expected, got);
            CHECK(got == Catch::Approx(expected).epsilon(1e-14));
            // matches to ~1-3 ULP; not bit-exact (FMA/-ffp-contract-class rounding)
            ++checks;
        }
    }
    CHECK(checks > 0);
}

TEST_CASE("golden: viscosity dilute powers_of_Tr", "[expression][golden]") {
    using namespace CoolProp::expression;
    auto HEOS = make_HEOS_for("Methane");  // viscosity/dilute type == powers_of_Tr
    auto& data = HEOS->get_components()[0].transport.viscosity_dilute.powers_of_Tr;
    REQUIRE(!data.a.empty());
    // C++: Tr = T/T_reducing; summer += a[i]*pow(Tr, t[i]);
    Program p = compile("let Tr = T/T_reducing\nsum(i: a[i]*Tr^t[i])", {{"T_reducing", static_cast<double>(data.T_reducing)}},
                        {{"a", std::vector<double>(data.a.begin(), data.a.end())}, {"t", std::vector<double>(data.t.begin(), data.t.end())}});
    int checks = 0;
    for (double T : {120.0, 200.0, 300.0, 450.0}) {
        for (double rho : {0.1, 100.0, 10000.0}) {
            HEOS->update(CoolProp::DmolarT_INPUTS, rho, T);
            double expected = CoolProp::TransportRoutines::viscosity_dilute_powers_of_Tr(*HEOS);
            std::vector<double> iv;
            fill_inputs(p, *HEOS, iv);
            double got = p.evaluate(iv);
            CAPTURE(T, rho, expected, got);
            CHECK(got == Catch::Approx(expected).epsilon(1e-14));
            CHECK(got == expected);  // bit-exact: same op sequence reproduces hardcoded value exactly
            ++checks;
        }
    }
    CHECK(checks > 0);
}

TEST_CASE("golden: viscosity dilute collision_integral", "[expression][golden]") {
    using namespace CoolProp::expression;
    auto HEOS = make_HEOS_for("IsoButane");  // viscosity/dilute type == collision_integral
    auto& data = HEOS->get_components()[0].transport.viscosity_dilute.collision_integral;
    REQUIRE(!data.a.empty());
    const auto eps_over_k = static_cast<double>(HEOS->get_components()[0].transport.epsilon_over_k);
    const auto sigma_eta = static_cast<double>(HEOS->get_components()[0].transport.sigma_eta);
    // C++:
    //   Tstar = T/epsilon_over_k;  sigma_nm = sigma_eta*1e9;  mm_kgkmol = molar_mass*1000;
    //   summer += a[i]*pow(log(Tstar), t[i]);  S = exp(summer);
    //   return C*sqrt(mm_kgkmol*T)/(pow(sigma_nm,2)*S);
    Program p = compile("let Tstar = T/epsilon_over_k\n"
                        "let sigma_nm = sigma_eta*1e9\n"
                        "let mm_kgkmol = molar_mass_data*1000\n"
                        "let lnTstar = ln(Tstar)\n"
                        "let S = exp(sum(i: a[i]*lnTstar^t[i]))\n"
                        "C*sqrt(mm_kgkmol*T)/(sigma_nm^2*S)",
                        {{"epsilon_over_k", eps_over_k},
                         {"sigma_eta", sigma_eta},
                         {"molar_mass_data", static_cast<double>(data.molar_mass)},
                         {"C", static_cast<double>(data.C)}},
                        {{"a", std::vector<double>(data.a.begin(), data.a.end())}, {"t", std::vector<double>(data.t.begin(), data.t.end())}});
    int checks = 0;
    for (double T : {200.0, 300.0, 400.0, 500.0}) {
        for (double rho : {0.1, 100.0, 5000.0}) {
            HEOS->update(CoolProp::DmolarT_INPUTS, rho, T);
            double expected = CoolProp::TransportRoutines::viscosity_dilute_collision_integral(*HEOS);
            std::vector<double> iv;
            fill_inputs(p, *HEOS, iv);
            double got = p.evaluate(iv);
            CAPTURE(T, rho, expected, got);
            CHECK(got == Catch::Approx(expected).epsilon(1e-14));
            CHECK(got == expected);  // bit-exact: ln/exp/sqrt/pow same order reproduces value exactly
            ++checks;
        }
    }
    CHECK(checks > 0);
}

TEST_CASE("golden: viscosity higher_order modified_Batschinski_Hildebrand", "[expression][golden]") {
    using namespace CoolProp::expression;
    auto HEOS = make_HEOS_for("IsoButane");  // viscosity/higher_order type == modified_Batschinski_Hildebrand
    auto& HO = HEOS->get_components()[0].transport.viscosity_higher_order.modified_Batschinski_Hildebrand;
    REQUIRE(!HO.a.empty());
    // C++:
    //   delta = rhomolar/rhomolar_reduce;  tau = T_reduce/T;
    //   S = sum a[i]*delta^d1[i]*tau^t1[i]*exp(gamma[i]*delta^l[i]);
    //   F = sum f[i]*delta^d2[i]*tau^t2[i];
    //   delta0 = (sum g[i]*tau^h[i]) / (sum p[i]*tau^q[i]);
    //   return S + F*(1/(delta0-delta) - 1/delta0);
    Program p = compile("let delta = rhomolar/rhomolar_reduce\n"
                        "let tau = T_reduce/T\n"
                        "let S = sum(i: a[i]*delta^d1[i]*tau^t1[i]*exp(gamma[i]*delta^l[i]))\n"
                        "let F = sum(i: f[i]*delta^d2[i]*tau^t2[i])\n"
                        "let numer = sum(i: g[i]*tau^h[i])\n"
                        "let denom = sum(i: pp[i]*tau^q[i])\n"
                        "let delta0 = numer/denom\n"
                        "S + F*(1/(delta0-delta) - 1/delta0)",
                        {{"rhomolar_reduce", static_cast<double>(HO.rhomolar_reduce)}, {"T_reduce", static_cast<double>(HO.T_reduce)}},
                        {{"a", std::vector<double>(HO.a.begin(), HO.a.end())},
                         {"d1", std::vector<double>(HO.d1.begin(), HO.d1.end())},
                         {"t1", std::vector<double>(HO.t1.begin(), HO.t1.end())},
                         {"gamma", std::vector<double>(HO.gamma.begin(), HO.gamma.end())},
                         {"l", std::vector<double>(HO.l.begin(), HO.l.end())},
                         {"f", std::vector<double>(HO.f.begin(), HO.f.end())},
                         {"d2", std::vector<double>(HO.d2.begin(), HO.d2.end())},
                         {"t2", std::vector<double>(HO.t2.begin(), HO.t2.end())},
                         {"g", std::vector<double>(HO.g.begin(), HO.g.end())},
                         {"h", std::vector<double>(HO.h.begin(), HO.h.end())},
                         {"pp", std::vector<double>(HO.p.begin(), HO.p.end())},
                         {"q", std::vector<double>(HO.q.begin(), HO.q.end())}});
    int checks = 0;
    for (double T : {200.0, 300.0, 400.0, 500.0}) {
        for (double rho : {100.0, 2000.0, 8000.0, 12000.0}) {
            HEOS->update(CoolProp::DmolarT_INPUTS, rho, T);
            double expected = CoolProp::TransportRoutines::viscosity_higher_order_modified_Batschinski_Hildebrand(*HEOS);
            std::vector<double> iv;
            fill_inputs(p, *HEOS, iv);
            double got = p.evaluate(iv);
            CAPTURE(T, rho, expected, got);
            CHECK(got == Catch::Approx(expected).epsilon(1e-14));
            // matches to ~1-3 ULP; not bit-exact (FMA/-ffp-contract-class rounding)
            ++checks;
        }
    }
    CHECK(checks > 0);
}

TEST_CASE("DSL throughput is comparable to the hardcoded routine", "[expression][!benchmark]") {
    using namespace CoolProp;
    using namespace CoolProp::expression;
    auto HEOS = make_HEOS_for("IsoButane");
    HEOS->update(DmolarT_INPUTS, 5000.0, 300.0);

    auto& HO = HEOS->get_components()[0].transport.viscosity_higher_order.modified_Batschinski_Hildebrand;
    REQUIRE(!HO.a.empty());
    Program prog = compile("let delta = rhomolar/rhomolar_reduce\n"
                           "let tau = T_reduce/T\n"
                           "let S = sum(i: a[i]*delta^d1[i]*tau^t1[i]*exp(gamma[i]*delta^l[i]))\n"
                           "let F = sum(i: f[i]*delta^d2[i]*tau^t2[i])\n"
                           "let numer = sum(i: g[i]*tau^h[i])\n"
                           "let denom = sum(i: pp[i]*tau^q[i])\n"
                           "let delta0 = numer/denom\n"
                           "S + F*(1/(delta0-delta) - 1/delta0)",
                           {{"rhomolar_reduce", static_cast<double>(HO.rhomolar_reduce)}, {"T_reduce", static_cast<double>(HO.T_reduce)}},
                           {{"a", std::vector<double>(HO.a.begin(), HO.a.end())},
                            {"d1", std::vector<double>(HO.d1.begin(), HO.d1.end())},
                            {"t1", std::vector<double>(HO.t1.begin(), HO.t1.end())},
                            {"gamma", std::vector<double>(HO.gamma.begin(), HO.gamma.end())},
                            {"l", std::vector<double>(HO.l.begin(), HO.l.end())},
                            {"f", std::vector<double>(HO.f.begin(), HO.f.end())},
                            {"d2", std::vector<double>(HO.d2.begin(), HO.d2.end())},
                            {"t2", std::vector<double>(HO.t2.begin(), HO.t2.end())},
                            {"g", std::vector<double>(HO.g.begin(), HO.g.end())},
                            {"h", std::vector<double>(HO.h.begin(), HO.h.end())},
                            {"pp", std::vector<double>(HO.p.begin(), HO.p.end())},
                            {"q", std::vector<double>(HO.q.begin(), HO.q.end())}});
    auto corr = std::make_shared<ExpressionCorrelation>(std::move(prog));

    const int N = 2'000'000;
    volatile double sink = 0.0;
    auto t0 = std::chrono::steady_clock::now();
    for (int k = 0; k < N; ++k)
        sink += TransportRoutines::viscosity_higher_order_modified_Batschinski_Hildebrand(*HEOS);
    auto t1 = std::chrono::steady_clock::now();
    for (int k = 0; k < N; ++k)
        sink += corr->eval(*HEOS);
    auto t2 = std::chrono::steady_clock::now();

    double ns_hard = std::chrono::duration<double, std::nano>(t1 - t0).count() / N;
    double ns_dsl = std::chrono::duration<double, std::nano>(t2 - t1).count() / N;
    double ratio = ns_dsl / ns_hard;
    WARN("hardcoded = " << ns_hard << " ns/eval; DSL = " << ns_dsl << " ns/eval; ratio = " << ratio << "x");
    CHECK(ratio < 5.0);
}

TEST_CASE("golden: conductivity dilute ratio_of_polynomials", "[expression][golden]") {
    using namespace CoolProp::expression;
    auto HEOS = make_HEOS_for("IsoButane");  // conductivity/dilute type == ratio_of_polynomials
    auto& data = HEOS->get_components()[0].transport.conductivity_dilute.ratio_polynomials;
    REQUIRE(!data.A.empty());
    // C++: Tr = T/T_reducing; sum1 = A[i]*Tr^n[i]; sum2 = B[i]*Tr^m[i]; return sum1/sum2;
    Program p = compile("let Tr = T/T_reducing\n"
                        "sum(i: A[i]*Tr^n[i]) / sum(i: B[i]*Tr^m[i])",
                        {{"T_reducing", static_cast<double>(data.T_reducing)}},
                        {{"A", std::vector<double>(data.A.begin(), data.A.end())},
                         {"n", std::vector<double>(data.n.begin(), data.n.end())},
                         {"B", std::vector<double>(data.B.begin(), data.B.end())},
                         {"m", std::vector<double>(data.m.begin(), data.m.end())}});
    int checks = 0;
    for (double T : {150.0, 250.0, 350.0, 500.0}) {
        for (double rho : {0.1, 100.0, 5000.0}) {
            HEOS->update(CoolProp::DmolarT_INPUTS, rho, T);
            double expected = CoolProp::TransportRoutines::conductivity_dilute_ratio_polynomials(*HEOS);
            std::vector<double> iv;
            fill_inputs(p, *HEOS, iv);
            double got = p.evaluate(iv);
            CAPTURE(T, rho, expected, got);
            CHECK(got == Catch::Approx(expected).epsilon(1e-14));
            // matches to ~1-3 ULP; not bit-exact (FMA/-ffp-contract-class rounding)
            ++checks;
        }
    }
    CHECK(checks > 0);
}

TEST_CASE("golden: conductivity dilute eta0_and_poly", "[expression][golden]") {
    using namespace CoolProp::expression;
    auto HEOS = make_HEOS_for("Oxygen");  // conductivity/dilute type == eta0_and_poly
    auto& E = HEOS->get_components()[0].transport.conductivity_dilute.eta0_and_poly;
    REQUIRE(E.A.size() >= 1);
    // C++:
    //   eta0_uPas = calc_viscosity_dilute()*1e6;
    //   summer = A[0]*eta0_uPas; for i>=1: summer += A[i]*pow(tau, t[i]);  return summer;
    // eta0_uPas calls another routine (dilute viscosity) -> pass as a per-state constant.
    // tau == HEOS.tau() == T_reduce/T (EOS reducing T); pass T_reduce as a constant.
    // Split index 0 out (it multiplies eta0, not tau^t[0]) to keep the exact op order.
    std::vector<double> A_rest(E.A.begin() + 1, E.A.end());
    std::vector<double> t_rest(E.t.begin() + 1, E.t.end());
    int checks = 0;
    for (double T : {80.0, 120.0, 200.0, 300.0}) {
        for (double rho : {0.1, 100.0, 5000.0}) {
            HEOS->update(CoolProp::DmolarT_INPUTS, rho, T);
            const double eta0_uPas = static_cast<double>(HEOS->calc_viscosity_dilute()) * 1e6;
            const double T_reduce = HEOS->T_reducing();
            Program p = compile("let tau = T_reduce/T\n"
                                "A0*eta0_uPas + sum(i: A_rest[i]*tau^t_rest[i])",
                                {{"A0", static_cast<double>(E.A[0])}, {"eta0_uPas", eta0_uPas}, {"T_reduce", T_reduce}},
                                {{"A_rest", A_rest}, {"t_rest", t_rest}});
            double expected = CoolProp::TransportRoutines::conductivity_dilute_eta0_and_poly(*HEOS);
            std::vector<double> iv;
            fill_inputs(p, *HEOS, iv);
            double got = p.evaluate(iv);
            CAPTURE(T, rho, expected, got);
            CHECK(got == Catch::Approx(expected).epsilon(1e-14));
            // matches to ~1-3 ULP; not bit-exact (FMA/-ffp-contract-class rounding)
            ++checks;
        }
    }
    CHECK(checks > 0);
}

TEST_CASE("golden: conductivity residual polynomial", "[expression][golden]") {
    using namespace CoolProp::expression;
    auto HEOS = make_HEOS_for("IsoButane");  // conductivity/residual type == polynomial
    auto& data = HEOS->get_components()[0].transport.conductivity_residual.polynomials;
    REQUIRE(!data.B.empty());
    // C++:
    //   tau = T_reducing/T;  delta = rhomass/rhomass_reducing;
    //   summer += B[i]*pow(tau, t[i])*pow(delta, d[i]);
    Program p = compile("let tau = T_reducing/T\n"
                        "let delta = rhomass/rhomass_reducing\n"
                        "sum(i: B[i]*tau^t[i]*delta^d[i])",
                        {{"T_reducing", static_cast<double>(data.T_reducing)}, {"rhomass_reducing", static_cast<double>(data.rhomass_reducing)}},
                        {{"B", std::vector<double>(data.B.begin(), data.B.end())},
                         {"t", std::vector<double>(data.t.begin(), data.t.end())},
                         {"d", std::vector<double>(data.d.begin(), data.d.end())}});
    int checks = 0;
    for (double T : {200.0, 300.0, 400.0, 500.0}) {
        for (double rho : {100.0, 2000.0, 8000.0}) {
            HEOS->update(CoolProp::DmolarT_INPUTS, rho, T);
            double expected = CoolProp::TransportRoutines::conductivity_residual_polynomial(*HEOS);
            std::vector<double> iv;
            fill_inputs(p, *HEOS, iv);
            double got = p.evaluate(iv);
            CAPTURE(T, rho, expected, got);
            CHECK(got == Catch::Approx(expected).epsilon(1e-14));
            // matches to ~1-3 ULP; not bit-exact (FMA/-ffp-contract-class rounding)
            ++checks;
        }
    }
    CHECK(checks > 0);
}

TEST_CASE("golden: conductivity residual polynomial_and_exponential", "[expression][golden]") {
    using namespace CoolProp::expression;
    auto HEOS = make_HEOS_for("Oxygen");  // conductivity/residual type == polynomial_and_exponential
    auto& data = HEOS->get_components()[0].transport.conductivity_residual.polynomial_and_exponential;
    REQUIRE(!data.A.empty());
    // C++:
    //   tau = HEOS.tau();  delta = HEOS.delta();  (EOS reducing T and rhomolar)
    //   summer += A[i]*pow(tau, t[i])*pow(delta, d[i])*exp(-gamma[i]*pow(delta, l[i]));
    // HEOS.tau() == T_reduce/T, HEOS.delta() == rhomolar/rhomolar_reduce; pass reducing values.
    const double T_reduce = HEOS->T_reducing();
    const double rhomolar_reduce = HEOS->rhomolar_reducing();
    Program p = compile("let tau = T_reduce/T\n"
                        "let delta = rhomolar/rhomolar_reduce\n"
                        "sum(i: A[i]*tau^t[i]*delta^d[i]*exp(-gamma[i]*delta^l[i]))",
                        {{"T_reduce", T_reduce}, {"rhomolar_reduce", rhomolar_reduce}},
                        {{"A", std::vector<double>(data.A.begin(), data.A.end())},
                         {"t", std::vector<double>(data.t.begin(), data.t.end())},
                         {"d", std::vector<double>(data.d.begin(), data.d.end())},
                         {"gamma", std::vector<double>(data.gamma.begin(), data.gamma.end())},
                         {"l", std::vector<double>(data.l.begin(), data.l.end())}});
    int checks = 0;
    for (double T : {80.0, 120.0, 200.0, 300.0}) {
        for (double rho : {100.0, 5000.0, 20000.0}) {
            HEOS->update(CoolProp::DmolarT_INPUTS, rho, T);
            double expected = CoolProp::TransportRoutines::conductivity_residual_polynomial_and_exponential(*HEOS);
            std::vector<double> iv;
            fill_inputs(p, *HEOS, iv);
            double got = p.evaluate(iv);
            CAPTURE(T, rho, expected, got);
            CHECK(got == Catch::Approx(expected).epsilon(1e-14));
            CHECK(got == expected);  // bit-exact: same op sequence reproduces hardcoded value exactly
            ++checks;
        }
    }
    CHECK(checks > 0);
}

// ---------------------------------------------------------------------------
// ExpressionBlock: the standalone/scripting entry point.  compile_block() is the
// ONE definition of a `"type": "expression"` block -- the fluid library and this
// class both go through it -- so these tests pin the block contract itself.
// ---------------------------------------------------------------------------
TEST_CASE("ExpressionBlock compiles a fluid-JSON block verbatim", "[expression]") {
    using namespace CoolProp::expression;
    // The `"type"` key is carried along and ignored, so the text can be copied
    // straight out of a fluid file.
    ExpressionBlock b(R"({"type": "expression",
                          "formula": "sum(i: a[i]*T^t[i]) + c",
                          "constants": {"c": 1.0},
                          "arrays": {"a": [2.0, 3.0], "t": [0.0, 1.0]}})");
    REQUIRE(b.required_inputs() == std::vector<std::string>{"T"});
    auto HEOS = make_HEOS_for("R123");
    HEOS->update(CoolProp::DmolarT_INPUTS, 100.0, 400.0);
    CHECK(b.evaluate(*HEOS) == Catch::Approx(2.0 + 3.0 * 400.0 + 1.0));
}

TEST_CASE("ExpressionBlock reports bad blocks as ValueError", "[expression]") {
    using namespace CoolProp::expression;
    CHECK_THROWS_AS(ExpressionBlock("{not json}"), CoolProp::ValueError);
    CHECK_THROWS_AS(ExpressionBlock(R"({"type": "expression"})"), CoolProp::ValueError);  // no formula
    CHECK_THROWS_AS(ExpressionBlock(R"({"formula": "2 +"})"), CoolProp::ValueError);
    // The fluid library passes a context string; it must reach the message.
    try {
        compile_block(R"({"formula": "2 +"})", "fluid Bogus");
        FAIL("expected a ValueError");
    } catch (const CoolProp::ValueError& e) {
        CHECK(std::string(e.what()).find("fluid Bogus") != std::string::npos);
    }
}

TEST_CASE("ExpressionBlock reads whatever state the AbstractState is sitting at", "[expression]") {
    using namespace CoolProp::expression;
    ExpressionBlock b(R"({"formula": "p"})");
    REQUIRE(b.required_inputs() == std::vector<std::string>{"p"});
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "R123"));
    AS->update(CoolProp::DmolarT_INPUTS, 2000.0, 400.0);
    CHECK(b.evaluate(*AS) == Catch::Approx(CoolProp::PropsSI("P", "T", 400.0, "Dmolar", 2000.0, "R123")).epsilon(1e-14));
    // Re-set the state through the ordinary API and the same block follows it --
    // including through an input pair the block knows nothing about.
    AS->update(CoolProp::PT_INPUTS, 1e5, 350.0);
    CHECK(b.evaluate(*AS) == Catch::Approx(1e5).epsilon(1e-12));
}

TEST_CASE("ExpressionBlock refuses a state that was never updated", "[expression]") {
    using namespace CoolProp::expression;
    // AbstractState::p() on a fresh state returns -_HUGE rather than throwing, and
    // the formula would propagate it into a plausible-looking -inf.
    ExpressionBlock b(R"({"formula": "p*2"})");
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "R123"));
    // Pin the message, not just the type: any unrelated ValueError would otherwise
    // satisfy a bare CHECK_THROWS_AS and the guard could quietly stop working.
    CHECK_THROWS_WITH(b.evaluate(*AS), Catch::Matchers::ContainsSubstring("finite"));
    AS->update(CoolProp::DmolarT_INPUTS, 2000.0, 400.0);
    CHECK(ValidNumber(b.evaluate(*AS)));
}

TEST_CASE("ExpressionBlock works on a mixture state the caller set up", "[expression]") {
    using namespace CoolProp::expression;
    // Nothing in the block is pure-fluid-specific: the caller owns the state, so a
    // mixture composition is just another state the formula reads.
    ExpressionBlock b(R"({"formula": "molar_mass*rhomolar"})");
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "R32&R125"));
    AS->set_mole_fractions({0.4, 0.6});
    AS->update(CoolProp::PT_INPUTS, 1e5, 350.0);  // HEOS mixtures do not take DmolarT
    // NOT compared against AS->rhomass(): calc_rhomass() IS rhomolar()*molar_mass(),
    // so that would be true by construction even if the mixture molar mass were
    // wrong.  Build the expected molar mass from the pure-component values instead.
    const double M_mix = 0.4 * CoolProp::Props1SI("R32", "molar_mass") + 0.6 * CoolProp::Props1SI("R125", "molar_mass");
    CHECK(AS->rhomolar() > 0.0);
    CHECK(b.evaluate(*AS) == Catch::Approx(M_mix * AS->rhomolar()).epsilon(1e-14));
}

// ---------------------------------------------------------------------------
// Task 11: END-TO-END round-trip — a real fluid carrying a type:expression
// dilute-viscosity block flows through the ACTUAL load + dispatch path.
//
// Unlike the golden tests (which call compile()+the static routine directly),
// this test exercises:
//   JSON ("type":"expression") -> FluidLibrary::parse_expression_block
//     -> ExpressionData on the loaded CoolPropFluid
//     -> add_fluids_as_JSON("HEOS", ...) registers a NEW fluid
//     -> PropsSI("V", ...) -> calc_viscosity_dilute()'s
//        `case VISCOSITY_DILUTE_EXPRESSION:` dispatch arm
//        -> ExpressionCorrelation::eval.
//
// Construction: take R123's own full library JSON (dilute viscosity is the
// `powers_of_T` form, summer += a[i]*pow(T,t[i])), copy its a/t coefficients
// verbatim into a {"type":"expression"} block whose formula is the exact DSL
// translation `sum(i: a[i]*T^t[i])`, give the fluid a fresh NAME/CAS/ALIASES so
// it registers as a distinct entry, then compare PropsSI("V") of the new fluid
// against the original over a (T,rho) grid.  Every other EOS/transport stage is
// byte-identical to R123, so any disagreement is isolated to the dispatch arm.
// Agreement must be < 1e-14 relative -- tight enough that a lossy JSON
// round-trip of the coefficients (the dump()/parse hop the fluid library now
// takes) would show up rather than hide under the tolerance.
// ---------------------------------------------------------------------------
TEST_CASE("expression end-to-end: JSON load + dispatch round-trip (dilute viscosity)", "[expression]") {
    using nlohmann::json;

    // Pull R123's full library JSON (returned as a one-element array).
    json arr = json::parse(CoolProp::get_fluid_param_string("R123", "JSON"));
    REQUIRE(arr.is_array());
    REQUIRE(arr.size() == 1);
    json fluid = arr[0];

    // Sanity: the dilute viscosity block is the powers_of_T form we translate.
    REQUIRE(fluid.contains("TRANSPORT"));
    REQUIRE(fluid["TRANSPORT"].contains("viscosity"));
    json& visc = fluid["TRANSPORT"]["viscosity"];
    REQUIRE(visc.contains("dilute"));
    REQUIRE(visc["dilute"].at("type").get<std::string>() == "powers_of_T");

    // Copy the a/t coefficients verbatim from the original block.
    std::vector<double> a = visc["dilute"].at("a").get<std::vector<double>>();
    std::vector<double> t = visc["dilute"].at("t").get<std::vector<double>>();
    REQUIRE(!a.empty());
    REQUIRE(a.size() == t.size());

    // Replace the dilute block with an equivalent type:expression block.
    // C++ routine: summer += a[i]*pow(T, t[i]); => DSL: sum(i: a[i]*T^t[i]).
    json expr_block;
    expr_block["type"] = "expression";
    expr_block["formula"] = "sum(i: a[i]*T^t[i])";
    expr_block["arrays"]["a"] = a;
    expr_block["arrays"]["t"] = t;
    visc["dilute"] = expr_block;

    // Give the fluid a fresh identity so it registers as a NEW entry (no
    // collision with R123) -- name, CAS, and all aliases must be unique.
    fluid["INFO"]["NAME"] = "R123_EXPR_E2E";
    fluid["INFO"]["CAS"] = "999-99-90";  // synthetic, collision-free
    fluid["INFO"]["ALIASES"] = json::array({"R123_EXPR_E2E_ALIAS"});

    // Register through the real public load path.
    json doc = json::array();
    doc.push_back(fluid);
    REQUIRE(CoolProp::add_fluids_as_JSON("HEOS", doc.dump()));

    // Compare viscosity of the expression-backed fluid vs the original over a
    // (T,rho) grid -- exercising the VISCOSITY_DILUTE_EXPRESSION dispatch arm.
    int checks = 0;
    for (double T : {250.0, 300.0, 400.0, 500.0}) {
        for (double rho : {0.1, 100.0, 5000.0}) {
            double v_expr = CoolProp::PropsSI("V", "T", T, "Dmolar", rho, "R123_EXPR_E2E");
            double v_orig = CoolProp::PropsSI("V", "T", T, "Dmolar", rho, "R123");
            CAPTURE(T, rho, v_expr, v_orig);
            REQUIRE(ValidNumber(v_expr));
            REQUIRE(ValidNumber(v_orig));
            REQUIRE(v_orig != 0.0);
            CHECK(std::abs(v_expr - v_orig) / std::abs(v_orig) < 1e-14);
            ++checks;
        }
    }
    CHECK(checks > 0);
}

// ---------------------------------------------------------------------------
// Task 10/11 follow-up: derived-`p` HOST-path check.  Build an
// ExpressionCorrelation from compile("p", {}, {}) and evaluate it against a
// backend at a known state; the registry `p` getter is wired to HEOS.p(), so
// the result must equal HEOS->p() bit-for-bit (same getter, no algebra).  This
// proves the `p` input threads through ExpressionCorrelation::eval ->
// HEOS.keyed_output(iP) -> HEOS.p().
// ---------------------------------------------------------------------------
TEST_CASE("expression host path: p equals HEOS.p()", "[expression]") {
    using namespace CoolProp::expression;
    auto HEOS = make_HEOS_for("R123");
    Program prog = compile("p", {}, {});
    REQUIRE(prog.requiredInputs().size() == 1);
    REQUIRE(prog.requiredInputs()[0] == CoolProp::iP);
    ExpressionCorrelation corr(std::move(prog));
    int checks = 0;
    for (double T : {300.0, 400.0, 500.0}) {
        for (double rho : {100.0, 5000.0}) {
            HEOS->update(CoolProp::DmolarT_INPUTS, rho, T);
            double got = corr.eval(*HEOS);
            auto expected = static_cast<double>(HEOS->p());
            // Independently computed, so the check cannot pass just because both
            // sides read the same cached member: HEOS->p() and keyed_output(iP)
            // do, and comparing them alone would be a tautology.
            double independent = CoolProp::PropsSI("P", "T", T, "Dmolar", rho, "R123");
            CAPTURE(T, rho, got, expected, independent);
            REQUIRE(ValidNumber(got));
            REQUIRE(independent > 0.0);
            CHECK(got == expected);  // identical getter, no rounding divergence
            CHECK(got == Catch::Approx(independent).epsilon(1e-14));
            ++checks;
        }
    }
    CHECK(checks > 0);
}

#endif  // ENABLE_CATCH
