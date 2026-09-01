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

#    include "CoolProp/Configuration.h"
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
    Program p = compile("T + Dmolar", {}, {}, {"T", "Dmolar"});
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
    Program p = compile("P*2 + T", {}, {}, {"P", "T"});
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

TEST_CASE("declared state variables resolve through CoolProp's own vocabulary", "[expression]") {
    using namespace CoolProp::expression;
    // Names are CoolProp's own, resolved to CoolProp's own keys.  The set is opt-in
    // and deliberately small (see "the DSL exposes an explicit set" below); a state
    // function that is legitimate but unused -- Hmolar, say -- is simply not exposed
    // until a correlation needs it, which is one line rather than a design decision.
    CHECK(compile("T", {}, {}, {"T"}).requiredInputs()[0] == CoolProp::iT);
    CHECK(compile("Dmolar", {}, {}, {"Dmolar"}).requiredInputs()[0] == CoolProp::iDmolar);
    CHECK(compile("Dmass", {}, {}, {"Dmass"}).requiredInputs()[0] == CoolProp::iDmass);
    CHECK(compile("molar_mass", {}, {}, {"molar_mass"}).requiredInputs()[0] == CoolProp::imolar_mass);
    CHECK(compile("P", {}, {}, {"P"}).requiredInputs()[0] == CoolProp::iP);
    CHECK_THROWS_AS(compile("Hmolar", {}, {}, {"Hmolar"}), CoolProp::ValueError);
    CHECK(inputName(CoolProp::iDmolar) == "Dmolar");

    // The spellings are CoolProp's, with nothing invented.  The DSL used to accept
    // `rhomolar`/`rhomass`/`p`; the lowercase `p` in particular is what collided with
    // the exponent array every viscosity paper writes as p_i, and inventing it bought
    // nothing.  They are gone, and the error says so.
    for (const char* legacy : {"rhomolar", "rhomass", "p"}) {
        CAPTURE(legacy);
        CHECK_THROWS_AS(compile(legacy, {}, {}, {legacy}), CoolProp::ValueError);
    }
}

TEST_CASE("the DSL exposes an explicit set of state variables", "[expression]") {
    using namespace CoolProp::expression;
    // Opt-in, and small on purpose.  Each name must be one of CoolProp's own AND its
    // canonical spelling -- an alias here would silently re-admit `A` for the speed of
    // sound and the single-letter coefficient names papers use.
    CoolProp::parameters key;
    std::string reason;
    const std::vector<std::string> expected = {"T", "P", "Dmolar", "Dmass", "molar_mass", "Smolar_residual", "Bvirial", "dBvirial_dT"};
    for (const auto& nm : expected) {
        CAPTURE(nm);
        REQUIRE(resolveStateVariable(nm, key, reason));
        CHECK(inputName(key) == nm);                                    // canonical, not an alias
        CHECK(compile(nm, {}, {}, {nm}).requiredInputs().size() == 1);  // and actually declarable
    }
}

TEST_CASE("anything outside that set is refused, and the message names the set", "[expression]") {
    using namespace CoolProp::expression;
    // A denylist was tried first and fails OPEN: 59 of CoolProp's 92 canonical names
    // slipped past one that had already grown to four categories.  Each name below is
    // a distinct way that went wrong, and all are refused by construction now.
    const std::pair<const char*, const char*> refused[] = {
      {"viscosity", "a transport output; keyed_output() re-enters the correlation being defined"},
      {"V", "the same, spelled as an alias"},
      {"T_critical", "configuration-dependent under ENABLE_SUPERANCILLARIES"},
      {"T_reducing", "the EOS reducing state is not the correlation's own fitted parameter"},
      {"Tau", "the reducing state wearing another name: _reducing.T/_T"},
      {"Delta", "likewise: _rhomolar/_reducing.rhomolar"},
      {"alphar", "an EOS internal in REDUCED variables, so it carries the reducing state too"},
      {"dalphar_ddelta_consttau", "the same, and a denylist missed it"},
      {"Phase", "an enum ordinal widened to double"},
      {"Q", "-1 outside the dome: a sentinel entering arithmetic as data"},
      {"HFORMATION", "fluid metadata; evaluates to a number that means nothing here"},
      {"P_min", "a range limit, not state -- and it reads back p_triple()"},
      {"T_triple", "likewise a fluid constant"},
      {"gas_constant", "a fluid constant; freeze it if the correlation needs it"},
      {"A", "CoolProp's alias for speed_of_sound -- and a coefficient name in every paper"},
      {"M", "CoolProp's alias for molar_mass"},
      {"DMOLAR", "an upper-cased spelling; one name per quantity"},
    };
    CoolProp::parameters key;
    std::string reason;
    for (const auto& r : refused) {
        CAPTURE(r.first, r.second);
        REQUIRE_FALSE(resolveStateVariable(r.first, key, reason));
        // The message answers "then what do I write?" by naming the whole set.
        CHECK(reason.find("Dmolar") != std::string::npos);
        CHECK_THROWS_AS(compile(r.first, {}, {}, {r.first}), CoolProp::ValueError);
        // ...and every one of them stays usable as the author's own frozen constant,
        // which is where a correlation's reducing parameters belong anyway.
        CHECK(compile(r.first, {{r.first, 7.0}}, {}).evaluate({}) == Catch::Approx(7.0));
    }
    // The DSL's own retired spellings land here too, and the same message answers
    // them: `rhomolar` is met with `Dmolar` in the available set.
    for (const char* old : {"rhomolar", "rhomass", "p"}) {
        CAPTURE(old);
        REQUIRE_FALSE(resolveStateVariable(old, key, reason));
        CHECK_THROWS_AS(compile(old, {}, {}, {old}), CoolProp::ValueError);
    }
    // `p` in particular stays the author's, which is the collision this all fixed.
    CHECK(compile("sum(i: n[i]*p[i])", {}, {{"n", {2.0, 3.0}}, {"p", {5.0, 7.0}}}).evaluate({}) == Catch::Approx(2.0 * 5.0 + 3.0 * 7.0));
}

TEST_CASE("two names for one quantity are rejected rather than fetched twice", "[expression]") {
    using namespace CoolProp::expression;
    // Dedupe is on the resolved key: two slots for one quantity would cost two
    // keyed_output() calls per evaluation.
    CHECK_THROWS_AS(compile("T + T", {}, {}, {"T", "T"}), CoolProp::ValueError);
}

TEST_CASE("a name the block did not declare is the author's to use", "[expression]") {
    using namespace CoolProp::expression;
    // THE point of opt-in.  Propylene glycol's dilute term is eta_0 = sum n_i Tr^p_i;
    // under one shared namespace its exponent array had to be renamed to np_i for a
    // collision with pressure, a quantity that block never reads.  Now `p` is simply
    // the author's, because the block never asked for pressure.
    Program q = compile("sum(i: n[i]*p[i])", {}, {{"n", {2.0, 3.0}}, {"p", {5.0, 7.0}}});
    CHECK(q.requiredInputs().empty());
    CHECK(q.evaluate({}) == Catch::Approx(2.0 * 5.0 + 3.0 * 7.0));
    // Same name, same formula, but declared: now it is pressure, and the array of
    // that name is a collision rather than a silent shadowing.
    CHECK_THROWS_AS(compile("sum(i: n[i]*p[i])", {}, {{"n", {2.0}}, {"p", {5.0}}}, {"P"}), CoolProp::ValueError);
    CHECK_THROWS_AS(compile("T", {{"T", 1.0}}, {}, {"T"}), CoolProp::ValueError);
    CHECK_THROWS_AS(compile("T", {}, {{"T", {1.0}}}, {"T"}), CoolProp::ValueError);

    // Reading a thermodynamic quantity you did not declare is an error, and the
    // message says which fix is wanted rather than "unknown variable".
    try {
        compile("Bvirial", {}, {});
        FAIL("expected a throw");
    } catch (const CoolProp::ValueError& e) {
        const std::string msg = e.what();
        CHECK(msg.find("state_variables") != std::string::npos);
    }
    // Declaring something the formula never reads is also an error: it would cost a
    // keyed_output() per evaluation and misstate the block's dependencies.
    CHECK_THROWS_AS(compile("1", {}, {}, {"T"}), CoolProp::ValueError);
    CHECK_THROWS_AS(compile("T", {}, {}, {"T", "T"}), CoolProp::ValueError);
}

TEST_CASE("`let` renames a state variable to whatever reads best", "[expression]") {
    using namespace CoolProp::expression;
    // Dropping the invented aliases (`rhomolar`, `rhomass`, `p`) costs nothing
    // ergonomically: renaming is already in the language, and a `let` can carry the
    // source paper's own symbol rather than either CoolProp's spelling or a second
    // vocabulary maintained here.  This is the answer to "but I preferred rhomolar".
    Program p = compile("let rho = Dmolar\nlet tau = Tc/T\nrho*tau", {{"Tc", 456.83}}, {}, {"Dmolar", "T"});
    REQUIRE(p.requiredInputs().size() == 2);
    CHECK(p.requiredInputs()[0] == CoolProp::iDmolar);
    CHECK(p.requiredInputs()[1] == CoolProp::iT);
    CHECK(p.evaluate({1e4, 300.0}) == Catch::Approx(1e4 * 456.83 / 300.0));

    // A `let` may not shadow a DECLARED state variable at all.  Leaning on the
    // "declared but never used" rule to catch this was not enough: that only fires
    // when the declared name is read NOWHERE, so the formula below compiled, with
    // `T` meaning the state on one line and 500 on the next.
    CHECK_THROWS_AS(compile("let T = 5\nT", {}, {}, {"T"}), CoolProp::ValueError);
    CHECK_THROWS_AS(compile("let a = T*2\nlet T = 500\na + T", {}, {}, {"T"}), CoolProp::ValueError);
    // Undeclared, that same name is simply a local, with no state involved at all.
    Program q = compile("let T = 5\nT", {}, {});
    CHECK(q.requiredInputs().empty());
    CHECK(q.evaluate({}) == Catch::Approx(5.0));
}

TEST_CASE("evaluate rejects a wrong-sized input vector", "[expression]") {
    using namespace CoolProp::expression;
    Program p = compile("T + P", {}, {}, {"T", "P"});
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
    Program p = compile("Dmass * molar_mass", {}, {}, {"Dmass", "molar_mass"});
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
    Program p =
      compile("let delta = Dmolar/rhomolar_reduce\nlet tau = T_reduce/T\nsum(i: a[i]*delta^d1[i]*tau^t1[i])", consts, arrays, {"Dmolar", "T"});
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
                        {{"a", std::vector<double>(data.a.begin(), data.a.end())}, {"t", std::vector<double>(data.t.begin(), data.t.end())}}, {"T"});
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
    // n-Pentane, not methane: methane's viscosity is now the Sotiriadou (2025)
    // expression model, so it no longer carries a powers_of_Tr dilute block for this
    // golden test to read.  n-Pentane is the only remaining fluid with that form.
    auto HEOS = make_HEOS_for("n-Pentane");  // viscosity/dilute type == powers_of_Tr
    auto& data = HEOS->get_components()[0].transport.viscosity_dilute.powers_of_Tr;
    REQUIRE(!data.a.empty());
    // C++: Tr = T/T_reducing; summer += a[i]*pow(Tr, t[i]);
    Program p = compile("let Tr = T/T_reducing\nsum(i: a[i]*Tr^t[i])", {{"T_reducing", static_cast<double>(data.T_reducing)}},
                        {{"a", std::vector<double>(data.a.begin(), data.a.end())}, {"t", std::vector<double>(data.t.begin(), data.t.end())}}, {"T"});
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
                        {{"a", std::vector<double>(data.a.begin(), data.a.end())}, {"t", std::vector<double>(data.t.begin(), data.t.end())}}, {"T"});
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
    Program p = compile("let delta = Dmolar/rhomolar_reduce\n"
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
                         {"q", std::vector<double>(HO.q.begin(), HO.q.end())}},
                        {"Dmolar", "T"});
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
    Program prog = compile("let delta = Dmolar/rhomolar_reduce\n"
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
                            {"q", std::vector<double>(HO.q.begin(), HO.q.end())}},
                           {"Dmolar", "T"});
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
                         {"m", std::vector<double>(data.m.begin(), data.m.end())}},
                        {"T"});
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
                                {{"A_rest", A_rest}, {"t_rest", t_rest}}, {"T"});
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
                        "let delta = Dmass/rhomass_reducing\n"
                        "sum(i: B[i]*tau^t[i]*delta^d[i])",
                        {{"T_reducing", static_cast<double>(data.T_reducing)}, {"rhomass_reducing", static_cast<double>(data.rhomass_reducing)}},
                        {{"B", std::vector<double>(data.B.begin(), data.B.end())},
                         {"t", std::vector<double>(data.t.begin(), data.t.end())},
                         {"d", std::vector<double>(data.d.begin(), data.d.end())}},
                        {"T", "Dmass"});
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
                        "let delta = Dmolar/rhomolar_reduce\n"
                        "sum(i: A[i]*tau^t[i]*delta^d[i]*exp(-gamma[i]*delta^l[i]))",
                        {{"T_reduce", T_reduce}, {"rhomolar_reduce", rhomolar_reduce}},
                        {{"A", std::vector<double>(data.A.begin(), data.A.end())},
                         {"t", std::vector<double>(data.t.begin(), data.t.end())},
                         {"d", std::vector<double>(data.d.begin(), data.d.end())},
                         {"gamma", std::vector<double>(data.gamma.begin(), data.gamma.end())},
                         {"l", std::vector<double>(data.l.begin(), data.l.end())}},
                        {"T", "Dmolar"});
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
  "state_variables": ["T"], "constants": {"c": 1.0},
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
    ExpressionBlock b(R"({"formula": "P",
  "state_variables": ["P"] })");
    REQUIRE(b.required_inputs() == std::vector<std::string>{"P"});
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
    ExpressionBlock b(R"({"formula": "P*2",
  "state_variables": ["P"] })");
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
    ExpressionBlock b(R"({"formula": "molar_mass*Dmolar",
  "state_variables": ["molar_mass", "Dmolar"] })");
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
    // The block declares what it reads; `a` and `t` stay the author's names.
    expr_block["state_variables"] = json::array({"T"});
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
// ExpressionCorrelation from compile("P", {}, {}, {"P"}) and evaluate it against a
// backend at a known state; the registry `p` getter is wired to HEOS.p(), so
// the result must equal HEOS->p() bit-for-bit (same getter, no algebra).  This
// proves the `p` input threads through ExpressionCorrelation::eval ->
// HEOS.keyed_output(iP) -> HEOS.p().
// ---------------------------------------------------------------------------
TEST_CASE("expression host path: p equals HEOS.p()", "[expression]") {
    using namespace CoolProp::expression;
    auto HEOS = make_HEOS_for("R123");
    Program prog = compile("P", {}, {}, {"P"});
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

// ---------------------------------------------------------------------------
// GOLDEN: the 5th transport stage, viscosity initial_density (empirical form).
//
// CycloHexane is the only shipped fluid using the empirical form; the other 19
// initial-density fluids use Rainwater-Friend, whose host arm scales the routine's
// B_eta by eta_dilute*rho.  eta_dilute is a within-correlation intermediate that
// the DSL cannot see, so an expression block yields this stage's contribution
// DIRECTLY (the empirical convention).  See the dispatch comment.
// ---------------------------------------------------------------------------
TEST_CASE("golden: viscosity initial_density empirical", "[expression][golden]") {
    using namespace CoolProp::expression;
    auto HEOS = make_HEOS_for("CycloHexane");
    auto& data = HEOS->get_components()[0].transport.viscosity_initial.empirical;
    REQUIRE(!data.n.empty());
    // C++: tau = T_reducing/T; delta = rhomolar/rhomolar_reducing;
    //      summer += n[i]*pow(delta,d[i])*pow(tau,t[i]);
    // The constants keep the routine's own names: the EOS reducing state is
    // deliberately NOT in the input table, precisely so `T_reducing` stays available
    // to correlations that mean their own fitted value (see the test below).
    Program p = compile("let tau = T_reducing/T\nlet delta = Dmolar/rhomolar_reducing\nsum(i: n[i]*delta^d[i]*tau^t[i])",
                        {{"T_reducing", static_cast<double>(data.T_reducing)}, {"rhomolar_reducing", static_cast<double>(data.rhomolar_reducing)}},
                        {{"n", std::vector<double>(data.n.begin(), data.n.end())},
                         {"d", std::vector<double>(data.d.begin(), data.d.end())},
                         {"t", std::vector<double>(data.t.begin(), data.t.end())}},
                        {"T", "Dmolar"});
    int checks = 0;
    for (double T : {300.0, 400.0, 500.0}) {
        for (double rho : {100.0, 2000.0, 8000.0}) {
            HEOS->update(CoolProp::DmolarT_INPUTS, rho, T);
            double expected = CoolProp::TransportRoutines::viscosity_initial_density_dependence_empirical(*HEOS);
            std::vector<double> iv;
            fill_inputs(p, *HEOS, iv);
            double got = p.evaluate(iv);
            CAPTURE(T, rho, expected, got);
            CHECK(got == Catch::Approx(expected).epsilon(1e-14));
            ++checks;
        }
    }
    CHECK(checks > 0);
}

// ---------------------------------------------------------------------------
// Huber, Perkins & Lemmon, Int. J. Thermophys. 45:146 (2024), CC BY 4.0 --
// "Reference Correlation for the Viscosity of Nitrogen from the Triple Point to
// 1000 K and Pressures up to 2200 MPa".  A correlation published AFTER the DSL
// was designed, implemented as pure JSON data with no C++ added.
//
//   eta = (eta_0(T) + Delta_eta_res(T,rho)) * Delta_eta_c            Eq. 1
//   eta_0(T) = eta_0(298.15 K) * exp{ sum_i a_i [ln(T/298.15)]^i }   Eq. 3
//   Delta_eta_res = sum_i N_i Tr^t_i rhor^d_i                        Eq. 5
//
// Delta_eta_c (the crossover critical enhancement) is NOT expressible today and is
// out of scope here; the paper's Table 7 sets it to 1, which is what we check.
// Reducing parameters are frozen constants.  Both Span et al. and this correlation
// use rho_c = 11.1839 mol/L; CoolProp's Nitrogen.json currently carries
// 11183.901464580624 mol/m^3 from a 2020 unit-conversion defect (under audit), and
// rho_r^8.4 turns that difference into a ~1e-6 shift -- enough to lose the paper's
// own check values.
// ---------------------------------------------------------------------------
namespace {
// clang-format off
const char* const N2_DILUTE_2024 = R"JSON({
  "type": "expression",
  "formula": "let L = ln(T/T_ref)\n1e-6*eta_ref*exp(sum(i: a[i]*L^p[i]))",
  "state_variables": ["T"],
  "constants": {"T_ref": 298.15, "eta_ref": 17.7494},
  "arrays": {
    "a": [7.734578e-1, -9.310761e-2, 2.716958e-2, 6.175553e-3, -7.201594e-3,
          2.094372e-3, 1.922676e-4, -3.454323e-4, 1.051771e-4, -1.126739e-5],
    "p": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]}
})JSON";

const char* const N2_RESIDUAL_2024 = R"JSON({
  "type": "expression",
  "formula": "let Tr = T/Tc\nlet rhor = Dmolar/1000/rhoc\n1e-6*sum(i: N[i]*Tr^t[i]*rhor^d[i])",
  "state_variables": ["T", "Dmolar"],
  "constants": {"Tc": 126.192, "rhoc": 11.1839},
  "arrays": {
    "N": [9.955235691668, -6.165266404871, 0.213120936996, -8.473713006806, 10.013103356639,
          0.638966874603, 0.311620258213, 9.241856768911, -5.252828814854, -0.667072279228],
    "t": [-0.77512218631, -2.00109608805, -5.814455, -0.67602596996, 0,
          -0.71613185456, -1.14069399597, -2.3652583598, -2.5463699255, -1.00794034515],
    "d": [1, 1, 2, 2, 2, 6.96423363249, 8.38651257501, 2.99224857471, 3.42744617975, 7.97371767411]}
})JSON";
// clang-format on
constexpr double MUPAS = 1e-6;  // the paper tabulates muPa.s; the DSL yields Pa.s
}  // namespace

TEST_CASE("Nitrogen 2024: dilute and residual blocks match the paper's Table 8", "[expression][golden]") {
    using namespace CoolProp::expression;
    // Table 8 lists eta_0 and eta_res SEPARATELY, so each block is pinned on its own.
    ExpressionBlock dilute(N2_DILUTE_2024), residual(N2_RESIDUAL_2024);
    CHECK(dilute.required_inputs() == std::vector<std::string>{"T"});
    CHECK(residual.required_inputs() == std::vector<std::string>{"T", "Dmolar"});

    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Nitrogen"));
    // T (K), rho (kg/m^3), eta_0 (muPa.s), eta_res (muPa.s)
    const double tab8[3][4] = {
      {126.192, 265, 8.43716205, 7.20032032}, {126.212, 333, 8.43843818, 11.03805357}, {126.952, 300, 8.48562461, 9.04720938}};
    for (const auto& row : tab8) {
        AS->update(CoolProp::DmassT_INPUTS, row[1], row[0]);
        CAPTURE(row[0], row[1]);
        CHECK(dilute.evaluate(*AS) == Catch::Approx(row[2] * MUPAS).epsilon(1e-9));
        CHECK(residual.evaluate(*AS) == Catch::Approx(row[3] * MUPAS).epsilon(1e-9));
    }
}

TEST_CASE("Nitrogen 2024: end-to-end through the fluid library matches Table 7", "[expression]") {
    using nlohmann::json;
    // Graft the two blocks onto a copy of Nitrogen and register it as a new fluid,
    // so PropsSI("V") exercises the real load + dispatch path, not compile() direct.
    json fluid = json::parse(CoolProp::get_fluid_param_string("Nitrogen", "JSON"))[0];
    json visc = json::object();
    // The loader requires a BibTeX member on the viscosity block; cite the new paper.
    visc["BibTeX"] = "Huber-IJT-2024";
    visc["dilute"] = json::parse(N2_DILUTE_2024);
    visc["higher_order"] = json::parse(N2_RESIDUAL_2024);
    fluid["TRANSPORT"]["viscosity"] = visc;
    fluid["INFO"]["NAME"] = "N2_HUBER_2024";
    fluid["INFO"]["CAS"] = "999-99-91";
    fluid["INFO"]["ALIASES"] = json::array({"N2_HUBER_2024_ALIAS"});
    REQUIRE(CoolProp::add_fluids_as_JSON("HEOS", json::array({fluid}).dump()));

    // Table 7, the paper's computer-verification values (Delta_eta_c = 1).  The
    // rho = 0 rows are covered by the dilute-block test above; PropsSI needs a
    // physical state, so only the finite-density rows go through here.
    const double tab7[3][3] = {{90, 756, 108.42550781}, {300, 28, 18.23478803}, {300, 560, 50.59605975}};
    for (const auto& row : tab7) {
        double got = CoolProp::PropsSI("V", "T", row[0], "Dmass", row[1], "N2_HUBER_2024");
        CAPTURE(row[0], row[1], row[2]);
        REQUIRE(ValidNumber(got));
        CHECK(got == Catch::Approx(row[2] * MUPAS).epsilon(1e-9));
    }
}

// ---------------------------------------------------------------------------
// Sotiriadou, Antoniadis, Assael & Huber, Int. J. Thermophys. 46:133 (2025) --
// "Reference Correlation of the Viscosity of Argon".  Three stages, all data:
//
//   eta   = eta_0(T) + eta_1(T)*rho + Delta_eta(rho,T)
//   eta_0 = eta_0(298.15 K) exp( sum_i a_i [ln(T/298.15)]^i )          Eq. 2
//   eta_1 = eta_0(T) B_eta(T),  B_eta = B*_eta N_A sigma^3             Eqs. 3,4
//   B*_eta(T*) = sum_{i=0..6} c_i (T*)^-i,  T* = T/(eps/kB)            Eq. 5
//   Delta_eta = rhor^(2/3) Tr^(1/2) { f1 rhor + f2 rhor^2/Tr
//               + (f1 rhor - rhor^2)/Tr^5
//               + (rhor - f3 rhor^5)/(rhor - f4 - Tr) - f5 }           Eq. 6
//
// Three things this exercises that nitrogen did not:
//
//  * A REAL Rainwater-Friend initial-density stage.  eta_1 needs eta_0, which is
//    a within-correlation intermediate the DSL cannot see -- so the block simply
//    RECOMPUTES eta_0.  Duplicated data, zero new C++.  Collapsing all three
//    stages into one block would be tidier but wrong: calc_viscosity_dilute() is
//    consumed independently by conductivity models (TransportRoutines.cpp:845,855)
//    and viscosity_contributions() is public API, so each stage must report only
//    its own contribution.
//  * TWO sums in one formula (12-term dilute, 7-term virial) in separate lets.
//    Only NESTED sums are forbidden; sequential ones each get their own index.
//  * Eq. 6 came out of SYMBOLIC REGRESSION -- not a sum of power terms at all,
//    but a rational expression with fractional exponents and a genuine pole at
//    rhor - f4 - Tr = 0, which the paper warns about.  Our IEEE-semantics policy
//    reproduces that rather than papering over it.
//
// WHY THE TABLE 9 GRID IS THE TEST, not the three verification points:
// Section 3.2 gives three points for checking an implementation, and all three
// are at T = 300 K.  There ln(T/298.15) ~ 6.2e-3, so the i=6 term of Eq. 2
// contributes ~1e-14 relative and the high-order coefficients are unexercised.
// A transcription error in a_6 (10^-5 misread as 10^-3) still matched all three
// points to 6e-8 while being 36 % wrong at 2000 K; only Table 9, which spans
// 100-2000 K, caught it.  Hence the 41-point grid below.
// ---------------------------------------------------------------------------
namespace {
// clang-format off
const char* const AR_A_COEFFS =
  R"([8.395115e-1, -1.062564e-1, 1.065796e-2, 1.879809e-2, -8.881774e-3, -9.613779e-5,
      1.404406e-3, -4.321739e-4, -2.544782e-5, 4.398471e-5, -9.997908e-6, 7.753453e-7])";
const char* const AR_P_POWERS = "[1,2,3,4,5,6,7,8,9,10,11,12]";

std::string ar_dilute() {
    return std::string(R"JSON({"type": "expression",
      "formula": "let L = ln(T/T_ref)\n1e-6*eta_ref*exp(sum(i: a[i]*L^p[i]))",
  "state_variables": ["T"],
      "constants": {"T_ref": 298.15, "eta_ref": 22.5666},
      "arrays": {"a": )JSON") + AR_A_COEFFS + R"JSON(, "p": )JSON" + AR_P_POWERS + R"JSON(}})JSON";
}
std::string ar_initial_density() {
    return std::string(R"JSON({"type": "expression",
      "formula": "let L = ln(T/T_ref)\nlet eta0 = 1e-6*eta_ref*exp(sum(i: a[i]*L^p[i]))\nlet Tstar = T/epsilon_over_k\nlet Bstar = sum(i: c[i]*Tstar^q[i])\neta0*Bstar*N_A*sigma^3*Dmolar",
  "state_variables": ["T", "Dmolar"],
      "constants": {"T_ref": 298.15, "eta_ref": 22.5666, "epsilon_over_k": 143.235,
                    "sigma": 0.33501e-9, "N_A": 6.02214076e23},
      "arrays": {"a": )JSON") + AR_A_COEFFS + R"JSON(, "p": )JSON" + AR_P_POWERS + R"JSON(,
                 "c": [-0.2571, 3.033, 1.144, -5.586, 3.089, -0.8824, -0.03856],
                 "q": [0, -1, -2, -3, -4, -5, -6]}})JSON";
}
const char* const AR_RESIDUAL = R"JSON({"type": "expression",
  "formula": "let Tr = T/Tc\nlet rhor = Dmolar/rhoc\n1e-6*rhor^(2/3)*Tr^0.5*(f1*rhor + f2*rhor^2/Tr + (f1*rhor - rhor^2)/Tr^5 + (rhor - f3*rhor^5)/(rhor - f4 - Tr) - f5)",
  "state_variables": ["T", "Dmolar"],
  "constants": {"Tc": 150.687, "rhoc": 13407.42965855612,
                "f1": 3.62648753859904, "f2": 6.655428299399591, "f3": 0.39751160825739,
                "f4": 2.6697983930209,  "f5": 0.0472018570860789}})JSON";
// clang-format on
}  // namespace

TEST_CASE("Argon 2025: two sums in one formula, in separate lets", "[expression]") {
    using namespace CoolProp::expression;
    // The initial-density block needs a 12-term and a 7-term sum in one formula.
    // Nested sums are rejected; sequential ones must each get their own index scope.
    ExpressionBlock b(ar_initial_density());
    CHECK(b.required_inputs() == std::vector<std::string>{"T", "Dmolar"});
    Program p = compile("let u = sum(i: a[i])\nlet v = sum(i: b[i])\nu*v", {}, {{"a", {1, 2, 3}}, {"b", {10, 20}}});
    CHECK(p.evaluate({}) == Catch::Approx(6.0 * 30.0));
}

TEST_CASE("Argon 2025: stages match the paper's computer-verification points", "[expression][golden]") {
    using namespace CoolProp::expression;
    ExpressionBlock dilute(ar_dilute()), initial(ar_initial_density()), residual(AR_RESIDUAL);
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Argon"));
    // Section 3.2: T (K), rho (kg/m^3), eta (muPa.s).  rho = 0 is the dilute limit.
    const double verif[3][3] = {{300, 0.0, 22.6840}, {300, 4.0, 22.7334}, {300, 700.0, 49.3360}};
    for (const auto& row : verif) {
        // rho = 0 is not a state the EOS will take; evaluate the stages at a valid
        // state and zero the density-dependent ones by hand for that row.
        const double rho = (row[1] > 0) ? row[1] : 1e-6;
        AS->update(CoolProp::DmassT_INPUTS, rho, row[0]);
        double got = dilute.evaluate(*AS);
        if (row[1] > 0) got += initial.evaluate(*AS) + residual.evaluate(*AS);
        CAPTURE(row[0], row[1]);
        CHECK(got == Catch::Approx(row[2] * 1e-6).epsilon(1e-5));
    }
}

TEST_CASE("Argon 2025: end-to-end over the paper's Table 9 grid", "[expression]") {
    using nlohmann::json;
    json fluid = json::parse(CoolProp::get_fluid_param_string("Argon", "JSON"))[0];
    json visc = json::object();
    visc["BibTeX"] = "Sotiriadou-IJT-2025";
    visc["dilute"] = json::parse(ar_dilute());
    visc["initial_density"] = json::parse(ar_initial_density());
    visc["higher_order"] = json::parse(AR_RESIDUAL);
    fluid["TRANSPORT"]["viscosity"] = visc;
    fluid["INFO"]["NAME"] = "AR_SOTIRIADOU_2025";
    fluid["INFO"]["CAS"] = "999-99-92";
    fluid["INFO"]["ALIASES"] = json::array({"AR_SOTIRIADOU_2025_ALIAS"});
    REQUIRE(CoolProp::add_fluids_as_JSON("HEOS", json::array({fluid}).dump()));

    // Table 9: T (K), rho (kg/m^3), eta (muPa.s), printed to 5 significant figures.
    const double tab9[][3] = {
      {100, 4.9152, 8.0810},  {150, 3.2255, 12.095},   {200, 2.4093, 15.889},   {400, 1.2012, 28.642},   {600, 0.80058, 38.804},
      {800, 0.60042, 47.571}, {1000, 0.48034, 55.485}, {1500, 0.32025, 73.039}, {2000, 0.24019, 88.656}, {100, 1349.4, 204.28},
      {150, 964.88, 67.573},  {200, 337.74, 23.007},   {400, 119.43, 30.618},   {600, 78.026, 39.869},   {800, 58.472, 48.214},
      {1000, 46.880, 55.893}, {1500, 31.435, 73.167},  {2000, 23.672, 88.666},  {100, 1448.5, 290.45},   {150, 1234.3, 129.45},
      {200, 1023.7, 79.043},  {400, 511.19, 43.770},   {600, 342.22, 46.526},   {800, 261.22, 52.428},   {1000, 212.57, 58.820},
      {1500, 146.25, 74.560}, {2000, 111.88, 89.378},  {150, 1363.4, 187.48},   {200, 1213.1, 121.19},   {400, 787.37, 61.681},
      {600, 574.15, 56.530},  {800, 454.97, 59.113},   {1000, 378.65, 63.713},  {1500, 269.06, 77.260},  {2000, 209.61, 91.054},
      {150, 1510.1, 309.68},  {200, 1398.8, 198.36},   {400, 1065.5, 93.995},   {600, 856.29, 76.295},   {800, 717.53, 73.178},
      {1000, 619.24, 74.551}};
    // Two ways of hitting the same 41 points, with different error floors.
    //
    // Table 9 prints rho to five significant figures, and in the dense liquid
    // dln(eta)/dln(rho) ~ 6, so feeding the TABLE's rho amplifies its rounding
    // (drho/rho ~ 2.6e-5 at 50 MPa / 100 K) into ~1.6e-4 in viscosity.  That is
    // the table's precision, not ours.  Letting the EOS supply rho from (p,T) --
    // CoolProp's argon EOS is Tegeler-JPCRD-1999 with R = 8.31451, exactly the
    // equation and gas constant of the paper's Table 1 -- drops the worst
    // deviation to ~4e-5.  The implementation's true fidelity is better still:
    // the Section 3.2 verification points, whose densities are exact by
    // construction, agree to 6e-8 / 4.6e-7 / 8.7e-7 (see the test above).
    const double pres[41] = {0.1e6, 0.1e6, 0.1e6, 0.1e6, 0.1e6, 0.1e6, 0.1e6, 0.1e6, 0.1e6, 10e6,  10e6,  10e6,  10e6, 10e6,
                             10e6,  10e6,  10e6,  10e6,  50e6,  50e6,  50e6,  50e6,  50e6,  50e6,  50e6,  50e6,  50e6, 100e6,
                             100e6, 100e6, 100e6, 100e6, 100e6, 100e6, 100e6, 200e6, 200e6, 200e6, 200e6, 200e6, 200e6};
    static_assert(std::size(tab9) == 41, "tab9 and pres must stay in lockstep");
    static_assert(std::size(pres) == std::size(tab9), "tab9 and pres must stay in lockstep");
    double worst_d = 0, worst_pt = 0;
    int checks = 0;
    for (std::size_t k = 0; k < std::size(tab9); ++k) {
        const double T = tab9[k][0], rho = tab9[k][1], ref = tab9[k][2] * 1e-6;
        const double gd = CoolProp::PropsSI("V", "T", T, "Dmass", rho, "AR_SOTIRIADOU_2025");
        const double gp = CoolProp::PropsSI("V", "T", T, "P", pres[k], "AR_SOTIRIADOU_2025");
        CAPTURE(T, rho, pres[k], ref);
        REQUIRE(ValidNumber(gd));
        REQUIRE(ValidNumber(gp));
        worst_d = std::max(worst_d, std::abs(gd - ref) / ref);
        worst_pt = std::max(worst_pt, std::abs(gp - ref) / ref);
        CHECK(std::abs(gd - ref) / ref < 2e-4);  // limited by the table's rho
        CHECK(std::abs(gp - ref) / ref < 1e-4);  // EOS-supplied rho
        ++checks;
    }
    WARN("Argon 2025 vs Table 9 (41 pts): worst " << worst_d << " with table rho, " << worst_pt << " with EOS rho from (p,T)");
    CHECK(checks == 41);
}

// ---------------------------------------------------------------------------
// SHIPPED CORRELATIONS.  Unlike the nitrogen/argon fixtures above, these
// live in dev/fluids/*.json and are loaded through the ordinary fluid-data path
// (JSON -> all_fluids CBOR -> FluidLibrary), so PropsSI("V", ..., "<fluid>") uses
// them directly.  Purely additive: CoolProp shipped no viscosity for these fluids.
// ---------------------------------------------------------------------------

// Sotiriadou, Ntonti, Assael, Perkins & Huber, Int. J. Thermophys. 45(6):87 (2024),
// "Reference Correlation of the Viscosity of Ethene".
//
// Table 8 of that paper is the best verification data in this file: it resolves
// the STAGES at 283 K, listing eta_0, (eta_1*rho + Delta_eta), their sum, the
// critical-enhancement factor and the total, at seven significant figures for
// rho = 0 to 550 kg/m^3.  We implement the background only (Delta_eta_c is the
// Bhattacharjee crossover, which REFPROP omits for all but water and D2O), so the
// comparison is against the background column.
//
// Note on the molar mass: the initial-density term needs molar density, and the
// block freezes M = 0.02805316 kg/mol -- the paper's value.  Ethylene.json's EOS
// carries the same number today, so the two agree; freezing it is not correcting a
// mismatch but pinning the correlation to the value its authors regressed against,
// so a later revision of the fluid file's molar mass cannot silently move it.
// Velliadou, Tasidou, Antoniadis, Assael, Perkins & Huber, Int. J. Thermophys.
// 42(5):74 (2021), "Reference Correlation for the Viscosity of Xenon from the Triple
// Point to 750 K and up to 86 MPa".  CoolProp shipped no viscosity for xenon.
//
// Eq. 6 is implemented AS CORRECTED: the third residual term carries rho_r^12, not
// the rho_r^2 first printed.
//
// Tc = 289.733 K and rho_c = 8400 mol/m^3 are FROZEN constants, not read from the
// fluid.  Xenon.json carries exactly those values, so reading them looked safe --
// but with superancillaries on (the default) the critical-point accessors return the
// NUMERICAL critical point and Eq. 6 moves by 7e-5.  Reducing parameters belong to
// the correlation, not to whatever the EOS currently reports.
//
// The paper writes B_eta in m^3/kg and scales by mass density; on a molar basis the
// molar mass cancels (rho_mass/M == rho_molar), which is why none appears.
TEST_CASE("Xenon: shipped viscosity matches the paper's verification points", "[expression][golden]") {
    // Section 4, background viscosity: T (K), rho (kg/m^3), eta (muPa.s).  These are
    // the values the authors publish precisely so an implementation can be checked.
    const double verif[3][3] = {{300.0, 0.0, 23.1561}, {300.0, 6.0, 23.3186}, {300.0, 2500.0, 206.449}};
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Xenon"));
    double worst = 0;
    for (const auto& row : verif) {
        AS->update(CoolProp::DmassT_INPUTS, (row[1] > 0) ? row[1] : 1e-9, row[0]);
        const double got = static_cast<double>(AS->viscosity()) * 1e6, ref = row[2];
        CAPTURE(row[0], row[1], ref);
        REQUIRE(ValidNumber(got));
        const double rel = std::abs(got - ref) / ref;
        worst = std::max(worst, rel);
        CHECK(rel < 1e-5);
    }
    WARN("Xenon vs Section 4 verification points: worst relative deviation " << worst << " over 3 points");
}

TEST_CASE("Xenon: PropsSI viscosity now works", "[expression]") {
    const double eta = CoolProp::PropsSI("V", "T", 300.0, "Dmass", 2500.0, "Xenon");
    REQUIRE(ValidNumber(eta));
    CHECK(eta == Catch::Approx(206.449e-6).epsilon(1e-5));
    // ...and the stages sum to it, so the dispatch really is running all three.
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Xenon"));
    AS->update(CoolProp::DmassT_INPUTS, 2500.0, 300.0);
    CoolPropDbl d = 0, i = 0, r = 0, c = 0;
    AS->viscosity_contributions(d, i, r, c);
    CHECK(static_cast<double>(d + i + r + c) == Catch::Approx(eta).epsilon(1e-12));
    CHECK(static_cast<double>(d) > 0);
    CHECK(static_cast<double>(r) > 0);
}

// Velliadou, Assael & Huber, Int. J. Thermophys. 43(7):105 (2022) -- R-134a,
// superseding Huber, Laesecke & Perkins (2003).
//
// The one correlation in this batch taken from the PAPER rather than the FLD: its
// REFPROP encoding is a raw RPN stack building a custom G(T), where every operand is
// an unlabelled `cnst` and mis-assigning one is easy and silent.  The paper states
// Eqs. 2, 3 and 7 plainly.  That is the same two-source rule as everywhere else, just
// with the roles swapped -- take the form from whichever source states it
// unambiguously, then make the paper's own check values decide.
TEST_CASE("R134a: shipped viscosity matches the paper's verification points", "[expression][golden]") {
    // T (K), rho (kg/m^3), eta (muPa.s), given inline in Section 3.
    const double verif[3][3] = {{350.0, 0.0, 13.77874}, {350.0, 100.0, 14.70183}, {350.0, 1000.0, 107.98464}};
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "R134a"));
    double worst = 0;
    for (const auto& row : verif) {
        AS->update(CoolProp::DmassT_INPUTS, (row[1] > 0) ? row[1] : 1e-9, row[0]);
        const double got = static_cast<double>(AS->viscosity()) * 1e6, ref = row[2];
        CAPTURE(row[0], row[1], ref);
        REQUIRE(ValidNumber(got));
        const double rel = std::abs(got - ref) / ref;
        worst = std::max(worst, rel);
        CHECK(rel < 1e-6);
    }
    WARN("R134a vs verification points: worst relative deviation " << worst << " over 3 points");

    // Table 7 is a much stronger check than the three inline points: the whole
    // saturation line, 170-370 K, both branches.  It is what settles that the -8.8 %
    // this correlation moves the 185 K saturated liquid away from Huber, Laesecke &
    // Perkins (2003) is a real difference between the two correlations and not a
    // transcription error here -- the 2022 scheme's own saturation values are
    // reproduced across the entire line.
    // T (K), rho_liq, rho_vap (kg/m^3), eta_liq, eta_vap (muPa.s)
    const double tab7[11][5] = {
      {170.0, 1590.7, 0.028625, 1627.8, 6.76}, {190.0, 1537.5, 0.18259, 1040.6, 7.56}, {210.0, 1483.1, 0.76222, 702.2, 8.35},
      {230.0, 1426.8, 2.3660, 498.0, 9.12},    {250.0, 1367.9, 5.9546, 368.6, 9.86},   {270.0, 1305.1, 12.908, 281.4, 10.61},
      {290.0, 1236.8, 25.187, 218.6, 11.39},   {310.0, 1159.9, 45.786, 170.0, 12.29},  {330.0, 1069.1, 80.094, 129.8, 13.48},
      {350.0, 951.32, 140.99, 93.93, 15.44},   {370.0, 740.32, 293.90, 55.12, 21.04}};
    double worst7 = 0;
    for (const auto& row : tab7) {
        for (int br = 0; br < 2; ++br) {
            AS->update(CoolProp::DmassT_INPUTS, row[1 + br], row[0]);
            const double got = static_cast<double>(AS->viscosity()) * 1e6, ref = row[3 + br];
            CAPTURE(row[0], br, ref);
            REQUIRE(ValidNumber(got));
            const double rel = std::abs(got - ref) / ref;
            worst7 = std::max(worst7, rel);
            CHECK(rel < 1e-3);  // the table prints 3-4 significant figures
        }
    }
    WARN("R134a vs Table 7 saturation line: worst relative deviation " << worst7 << " over 22 points");
}

// Sotiriadou, Ntonti, Velliadou, Antoniadis, Assael & Huber, Int. J. Thermophys.
// 44(3):40 (2023) -- ethanol, superseding Kiselev et al. (2005).
// Sotiriadou, Antoniadis, Assael, Martinek & Huber, Int. J. Thermophys. 47(1):18
// (2025) -- methane, superseding Quinones-Cisneros & Deiters (2006) f-theory.
//
// Neither fluid's dilute conductivity is eta0_and_poly, so replacing their viscosity
// moves nothing else.
TEST_CASE("Ethanol: shipped viscosity matches the paper's verification points", "[expression][golden]") {
    // T (K), rho (kg/m^3), eta (muPa.s), given inline in Section 3.
    const double verif[3][3] = {{300.0, 0.0, 8.9893}, {300.0, 10.0, 8.9382}, {300.0, 850.0, 1682.72}};
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Ethanol"));
    double worst = 0;
    for (const auto& row : verif) {
        AS->update(CoolProp::DmassT_INPUTS, (row[1] > 0) ? row[1] : 1e-9, row[0]);
        const double got = static_cast<double>(AS->viscosity()) * 1e6, ref = row[2];
        CAPTURE(row[0], row[1], ref);
        REQUIRE(ValidNumber(got));
        const double rel = std::abs(got - ref) / ref;
        worst = std::max(worst, rel);
        CHECK(rel < 1e-5);
    }
    WARN("Ethanol vs verification points: worst relative deviation " << worst << " over 3 points");
}

TEST_CASE("Methane: shipped viscosity matches the paper's verification points", "[expression][golden]") {
    // T (K), rho (kg/m^3), eta (muPa.s), from the Computer-Program Verification section.
    const double verif[3][3] = {{300.0, 0.0, 11.1230}, {300.0, 3.2, 11.1891}, {300.0, 75.0, 13.7130}};
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Methane"));
    double worst = 0;
    for (const auto& row : verif) {
        AS->update(CoolProp::DmassT_INPUTS, (row[1] > 0) ? row[1] : 1e-9, row[0]);
        const double got = static_cast<double>(AS->viscosity()) * 1e6, ref = row[2];
        CAPTURE(row[0], row[1], ref);
        REQUIRE(ValidNumber(got));
        const double rel = std::abs(got - ref) / ref;
        worst = std::max(worst, rel);
        CHECK(rel < 1e-5);
    }
    WARN("Methane vs verification points: worst relative deviation " << worst << " over 3 points");
}

// Huber & Assael, Int. J. Refrigeration 71:39-45 (2016), one paper covering BOTH
// R-1234yf and R-1234ze(E).  Supersessions: each fluid previously carried a LIST of
// two models, Bell's rho_s r-CS (2016) plus an ECS fallback.
//
// Their dilute terms are ratios of polynomials in ABSOLUTE temperature, not reduced
// -- the FLD writes `$DG SUM:4 SUM:3 /` with no RED, and the coefficients
// (-836950, 6336.28, ...) only make sense that way.
//
// The two fluids agree to within 0.04 % at 10.522 mol/L, which looks like a
// transcription error and is not: they are isomers, both C3H2F4, and the paper uses
// one set of test densities for both.  Checked rather than assumed.
TEST_CASE("R1234yf and R1234ze(E): shipped viscosity matches the paper's verification points", "[expression][golden]") {
    // Section 2.4, Computer-Program Verification, inline rather than tabulated.
    // T = 300 K throughout; rho in mol/L; eta in muPa.s.
    struct Row
    {
        const char* fluid;
        double rho_molL, eta;
    };
    const Row rows[6] = {{"R1234yf", 0.0, 11.579},    {"R1234yf", 0.044, 11.549},    {"R1234yf", 10.522, 217.97},
                         {"R1234ze(E)", 0.0, 11.777}, {"R1234ze(E)", 0.044, 12.041}, {"R1234ze(E)", 10.522, 217.89}};
    double worst = 0;
    for (const auto& r : rows) {
        std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", r.fluid));
        AS->update(CoolProp::DmolarT_INPUTS, (r.rho_molL > 0) ? r.rho_molL * 1000.0 : 1e-9, 300.0);
        const double got = static_cast<double>(AS->viscosity()) * 1e6;
        CAPTURE(r.fluid, r.rho_molL, r.eta);
        REQUIRE(ValidNumber(got));
        const double rel = std::abs(got - r.eta) / r.eta;
        worst = std::max(worst, rel);
        CHECK(rel < 1e-4);
    }
    WARN("R1234yf/ze(E) vs Section 2.4: worst relative deviation " << worst << " over 6 points");
}

// Velliadou, Assael, Huber et al., Int. J. Thermophys. 43:129 (2022), "Reference
// Correlation for the Viscosity of R-32".  A SUPERSESSION: R32.json previously
// carried a LIST of two models, Bell's rho_s r-CS (2016) and an ECS fallback, both
// of which this replaces.
//
// Safe: R-32's dilute conductivity is not eta0_and_poly, so nothing else moves.  The
// four fluids where it IS -- Air, Argon, Nitrogen, Oxygen -- are all either out of
// this audit's scope or deliberately deferred (CoolProp-f1ez).
TEST_CASE("R32: shipped viscosity matches the paper's verification points", "[expression][golden]") {
    // Section 3 gives three points inline: T (K), rho (kg/m^3), eta (muPa.s).
    const double verif[3][3] = {{300.0, 0.0, 12.6170}, {300.0, 10.0, 12.6333}, {300.0, 1100.0, 173.431}};
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "R32"));
    double worst = 0;
    for (const auto& row : verif) {
        AS->update(CoolProp::DmassT_INPUTS, (row[1] > 0) ? row[1] : 1e-9, row[0]);
        const double got = static_cast<double>(AS->viscosity()) * 1e6, ref = row[2];
        CAPTURE(row[0], row[1], ref);
        REQUIRE(ValidNumber(got));
        const double rel = std::abs(got - ref) / ref;
        worst = std::max(worst, rel);
        CHECK(rel < 1e-5);
    }
    WARN("R32 vs Section 3 verification points: worst relative deviation " << worst << " over 3 points");
}

// Assael, Papalas & Huber, J. Phys. Chem. Ref. Data 46(3):033103 (2017), "Reference
// Correlations for the Viscosity and Thermal Conductivity of n-Undecane".  CoolProp
// shipped no viscosity for n-undecane.
//
// TWO stages, not three: the paper states that n-undecane is too large and
// non-spherical for the Rainwater-Friend treatment, so there is no initial-density
// term at all and its effect is absorbed into the residual.  The DSL expresses that
// by simply not declaring the stage.
//
// The residual's temperature factor is MULTIPLIED, not divided.  Eq. 8 renders as
// "rho_r^(2/3) / T_r^(1/2)", and implementing it that way misses the paper's own
// 550 K / 600 kg/m^3 point by 15 %.  REFPROP's C11.FLD has the numerator term at
// tau^0.5 with tau = T/T_red, i.e. multiplication, which reproduces the table to
// 1e-5.  Third correlation in this file where the rendered equation is wrong and the
// machine-readable form is right.
TEST_CASE("n-Undecane: shipped viscosity matches the paper's Table 11", "[expression][golden]") {
    // Table 11, sample points for computer verification: T (K), rho (kg/m^3), eta (muPa.s).
    const double tab11[5][3] = {{550.0, 0.0, 8.935}, {550.0, 10.0, 10.702}, {550.0, 600.0, 188.68}, {635.0, 0.0, 10.252}, {635.0, 325.0, 49.077}};
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "n-Undecane"));
    double worst = 0;
    for (const auto& row : tab11) {
        AS->update(CoolProp::DmassT_INPUTS, (row[1] > 0) ? row[1] : 1e-9, row[0]);
        const double got = static_cast<double>(AS->viscosity()) * 1e6, ref = row[2];
        CAPTURE(row[0], row[1], ref);
        REQUIRE(ValidNumber(got));
        const double rel = std::abs(got - ref) / ref;
        worst = std::max(worst, rel);
        CHECK(rel < 1e-4);
    }
    WARN("n-Undecane vs Table 11: worst relative deviation " << worst << " over 5 points");
}

// Monogenidou, Assael & Huber, J. Phys. Chem. Ref. Data 47(2):023102 (2018),
// "Reference Correlation for the Viscosity of Ammonia from the Triple Point to 725 K
// and up to 50 MPa".  A SUPERSESSION: this replaces Fenghour et al. (1995).
//
// Safe to replace because ammonia's dilute thermal conductivity is
// `ratio_of_polynomials`, which does not consume the dilute viscosity.  Nitrogen and
// argon are the two fluids in this audit where it does (`eta0_and_poly`), and they
// are deliberately NOT part of this batch -- see CoolProp-f1ez.
TEST_CASE("Ammonia: shipped viscosity matches the paper's verification points", "[expression][golden]") {
    // Section 3 gives three points inline rather than as a numbered table:
    // T (K), rho (kg/m^3), eta (muPa.s).
    const double verif[3][3] = {{300.0, 0.0, 10.1812}, {300.0, 8.0, 9.9219}, {300.0, 609.0, 133.3937}};
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Ammonia"));
    double worst = 0;
    for (const auto& row : verif) {
        AS->update(CoolProp::DmassT_INPUTS, (row[1] > 0) ? row[1] : 1e-9, row[0]);
        const double got = static_cast<double>(AS->viscosity()) * 1e6, ref = row[2];
        CAPTURE(row[0], row[1], ref);
        REQUIRE(ValidNumber(got));
        const double rel = std::abs(got - ref) / ref;
        worst = std::max(worst, rel);
        CHECK(rel < 1e-5);
    }
    WARN("Ammonia vs Section 3 verification points: worst relative deviation " << worst << " over 3 points");
}

// Perkins, Huber & Assael, J. Chem. Eng. Data 61(9):3286-3294 (2016).  A supersession
// of the rho_s r-CS + ECS predictive model; R-245fa's dilute conductivity does not
// consume the dilute viscosity, so nothing else moves.
TEST_CASE("R245fa: shipped viscosity matches the paper's Table 8", "[expression][golden]") {
    // Table 8, sample points for computer verification: T (K), rho (kg/m^3), eta (muPa.s).
    const double tab8[4][3] = {{250.0, 0.0, 8.6291}, {250.0, 1500.0, 1085.562}, {430.0, 0.0, 14.630}, {430.0, 530.0, 30.632}};
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "R245fa"));
    double worst = 0;
    for (const auto& row : tab8) {
        AS->update(CoolProp::DmassT_INPUTS, (row[1] > 0) ? row[1] : 1e-9, row[0]);
        const double got = static_cast<double>(AS->viscosity()) * 1e6, ref = row[2];
        CAPTURE(row[0], row[1], ref);
        REQUIRE(ValidNumber(got));
        const double rel = std::abs(got - ref) / ref;
        worst = std::max(worst, rel);
        CHECK(rel < 1e-4);
    }
    WARN("R245fa vs Table 8: worst relative deviation " << worst << " over 4 points");
}

// Wen, Meng, Huber & Wu, J. Chem. Eng. Data 62(10):3603-3609 (2017), "Measurement and
// Correlation of the Viscosity of 1,1,1,2,2,4,5,5,5-Nonafluoro-4-(trifluoromethyl)-3-
// pentanone" (Novec 649).  CoolProp shipped no viscosity for it.
//
// The dilute term is the Chung method over the Neufeld collision integral -- the one
// place in this file where the DSL's `sin` is load-bearing, since Neufeld's Omega*
// carries a sine correction term.  The Chapman-Enskog constant is taken from
// REFPROP's NOVEC649.FLD as a single frozen number (0.412899 = 0.02669*sqrt(MW)*Fc)
// rather than rebuilt from the paper's omega and reduced dipole moment: those give
// Fc = 0.87221 against the FLD's 0.87020, and the FLD's value is the one that
// reproduces the paper's own table.
//
// The residual is again a case where the printed equation cannot be implemented as
// rendered -- Eq. 12 reads "c0 + c1 c2 + c3 rho_r + ..." because the fraction bar is
// lost.  The FLD's descriptor `0 1 1 5 0 0` says c1 sits over a FIVE-term
// denominator, which is what reproduces Table 4.
TEST_CASE("Novec649: shipped viscosity matches the paper's Table 4", "[expression][golden]") {
    // Table 4, sample points for computer verification: T (K), rho (kg/m^3), eta (muPa.s).
    const double tab4[9][3] = {{250.0, 0.0, 8.09},       {250.0, 0.41, 8.33}, {250.0, 1809.77, 2377.5}, {300.0, 0.0, 9.77},      {300.0, 3.89, 10.85},
                               {300.0, 1701.48, 1059.7}, {350.0, 0.0, 11.43}, {350.0, 4.42, 12.65},     {350.0, 1595.99, 587.87}};
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Novec649"));
    double worst = 0;
    for (const auto& row : tab4) {
        AS->update(CoolProp::DmassT_INPUTS, (row[1] > 0) ? row[1] : 1e-9, row[0]);
        const double got = static_cast<double>(AS->viscosity()) * 1e6, ref = row[2];
        CAPTURE(row[0], row[1], ref);
        REQUIRE(ValidNumber(got));
        const double rel = std::abs(got - ref) / ref;
        worst = std::max(worst, rel);
        // The table prints three significant figures for the dilute rows, so half a
        // unit in the last place is 5e-4 there; the dense rows are far tighter.
        CHECK(rel < 6e-4);
    }
    WARN("Novec649 vs Table 4: worst relative deviation " << worst << " over 9 points");
}

TEST_CASE("Novec649: PropsSI viscosity now works", "[expression]") {
    const double eta = CoolProp::PropsSI("V", "T", 300.0, "Dmass", 1701.48, "Novec649");
    REQUIRE(ValidNumber(eta));
    CHECK(eta == Catch::Approx(1059.7e-6).epsilon(1e-4));
}

// Tsolakidou, Assael, Huber & Perkins, J. Phys. Chem. Ref. Data 46(2):023103 (2017),
// "Correlations for the Viscosity and Thermal Conductivity of Ethyl Fluoride (R161)".
// CoolProp shipped no viscosity for R-161.
//
// The residual form was taken from REFPROP 10.1's R161.FLD, not from the paper's
// printed Eq. 8, and the difference is not cosmetic.  Two independent readings of
// Eq. 8 as rendered both fail the paper's OWN verification table by 73 % and 91 % at
// 250 K / 850 kg/m^3.  The FLD's term descriptor `0 4 2 2 0 0` says the residual is
// FOUR simple-polynomial terms PLUS a separate two-over-two rational part -- not one
// quotient as the printed equation reads -- and its last denominator coefficient is
// -1.0 where the paper prints a plus.  With that structure the same coefficients
// reproduce the table to table precision.  The lesson is the one from krypton: the
// machine-readable source is authoritative for the FORM, and the paper's own check
// values are the independent test of it.
TEST_CASE("R161: shipped viscosity matches the paper's Table 11", "[expression][golden]") {
    // Table 11, sample points for computer verification: T (K), rho (kg/m^3),
    // eta (muPa.s).  Viscosity carries no critical enhancement here -- the table's
    // two 375 K / 229 kg/m^3 rows differ only in the conductivity column.
    const double tab11[5][3] = {{250.0, 0.0, 8.280}, {250.0, 1.0, 8.255}, {250.0, 850.0, 308.22}, {375.0, 0.0, 12.171}, {375.0, 229.0, 20.859}};
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "R161"));
    double worst = 0;
    for (const auto& row : tab11) {
        AS->update(CoolProp::DmassT_INPUTS, (row[1] > 0) ? row[1] : 1e-9, row[0]);
        const double got = static_cast<double>(AS->viscosity()) * 1e6, ref = row[2];
        CAPTURE(row[0], row[1], ref);
        REQUIRE(ValidNumber(got));
        const double rel = std::abs(got - ref) / ref;
        worst = std::max(worst, rel);
        CHECK(rel < 5e-5);  // the table is printed to 3-5 significant figures
    }
    WARN("R161 vs Table 11: worst relative deviation " << worst << " over 5 points");
}

TEST_CASE("R161: PropsSI viscosity now works", "[expression]") {
    const double eta = CoolProp::PropsSI("V", "T", 250.0, "Dmass", 850.0, "R161");
    REQUIRE(ValidNumber(eta));
    CHECK(eta == Catch::Approx(308.22e-6).epsilon(5e-5));
}

TEST_CASE("Ethylene: shipped viscosity matches the paper's Table 8", "[expression][golden]") {
    // T (K) = 283 throughout; rho (kg/m^3), eta_0, (eta_1 rho + Delta_eta), background (muPa.s)
    const double tab8[12][4] = {
      {0., 9.753447, 0.000000, 9.753447},     {50., 9.753447, 0.883644, 10.637091},   {100., 9.753447, 2.850625, 12.604071},
      {150., 9.753447, 5.758114, 15.511560},  {200., 9.753447, 9.627286, 19.380732},  {250., 9.753447, 14.637185, 24.390632},
      {300., 9.753447, 21.181254, 30.934701}, {350., 9.753447, 29.948520, 39.701967}, {400., 9.753447, 42.021090, 51.774536},
      {450., 9.753447, 58.985091, 68.738538}, {500., 9.753447, 83.053687, 92.807133}, {550., 9.753447, 117.201315, 126.954762}};
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Ethylene"));
    double worst = 0;
    for (const auto& row : tab8) {
        const double rho = (row[0] > 0) ? row[0] : 1e-8;  // the EOS will not take rho = 0
        AS->update(CoolProp::DmassT_INPUTS, rho, 283.0);
        CoolPropDbl dilute = 0, initial = 0, residual = 0, critical = 0;
        AS->viscosity_contributions(dilute, initial, residual, critical);
        const double total = AS->viscosity();
        CAPTURE(row[0], row[3]);
        REQUIRE(ValidNumber(total));
        // The dilute stage is temperature-only, so it must be the same at every row.
        CHECK(static_cast<double>(dilute) == Catch::Approx(row[1] * 1e-6).epsilon(1e-7));
        // initial_density + higher_order together are the paper's fourth column.
        // margin, not epsilon, for the rho = 0 row: we substitute rho = 1e-8 kg/m^3
        // there, and the residual's rhor^(2/3) leaves ~1e-13 Pa.s rather than exactly
        // zero.  1e-12 is still six orders below the smallest non-zero entry.
        CHECK(static_cast<double>(initial + residual) == Catch::Approx(row[2] * 1e-6).epsilon(1e-6).margin(1e-12));
        const double rel = std::abs(total - row[3] * 1e-6) / (row[3] * 1e-6);
        worst = std::max(worst, rel);
        CHECK(rel < 1e-6);
        CHECK(critical == 0.0);  // Delta_eta_c deliberately not implemented
    }
    WARN("Ethylene vs Table 8 background: worst relative deviation " << worst << " over 12 points");
}

TEST_CASE("Ethylene: PropsSI viscosity now works and is stage-consistent", "[expression]") {
    // Before this correlation shipped, PropsSI("V", ..., "Ethylene") threw.
    const double eta = CoolProp::PropsSI("V", "T", 283.0, "Dmass", 300.0, "Ethylene");
    REQUIRE(ValidNumber(eta));
    CHECK(eta == Catch::Approx(30.934701e-6).epsilon(1e-6));
    // And the three stages sum to the whole.
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Ethylene"));
    AS->update(CoolProp::DmassT_INPUTS, 300.0, 283.0);
    CoolPropDbl d = 0, i = 0, r = 0, c = 0;
    AS->viscosity_contributions(d, i, r, c);
    CHECK(static_cast<double>(d + i + r + c) == Catch::Approx(eta).epsilon(1e-14));
}

// Velliadou, Antoniadis, Assael & Huber, Int. J. Thermophys. 43(3):42 (2022), "Reference
// Correlation for the Viscosity of Propane-1,2-diol (Propylene Glycol)".
//
// New form for the DSL: the residual (Eq. 9) is an EXPONENTIAL,
//   Delta_eta = eta_ref (rhor^(2/3) Tr^(1/2)) exp{c0 + c1 rhor + c2 rhor^2/Tr
//                                                 + c3 rhor^3/Tr^2 + c4 Tr + c5 Tr^2},
// which is how a glycol spans four orders of magnitude in viscosity between 450 K
// and the triple point.  The dilute term (Eq. 5) is a plain polynomial in T/Tc --
// a Chapman-Enskog scheme refit for convenience -- and the initial-density term
// reuses the same universal Vogel/Bich B*(T*) as ethene and xenon.
TEST_CASE("PropyleneGlycol: shipped viscosity matches the paper's values", "[expression][golden]") {
    // Section 3 computer-verification points: T (K), rho (kg/m^3), eta (muPa.s).
    const double verif[3][3] = {{350, 0.0, 9.051368}, {350, 0.02, 9.058162}, {350, 1000.0, 5135.986461}};
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "PropyleneGlycol"));
    for (const auto& row : verif) {
        AS->update(CoolProp::DmassT_INPUTS, (row[1] > 0) ? row[1] : 1e-8, row[0]);
        CAPTURE(row[0], row[1]);
        CHECK(static_cast<double>(AS->viscosity()) == Catch::Approx(row[2] * 1e-6).epsilon(1e-6));
    }
    // Table 7, spanning 245-450 K and four orders of magnitude in viscosity.
    const double tab7[8][3] = {{245, 1072.3, 4683847.}, {300, 1031.2, 39323.}, {450, 903.51, 674.67}, {245, 1091.6, 9060255.},
                               {450, 945.09, 885.40},   {300, 1072.7, 88176.}, {450, 998.75, 1289.5}, {450, 1033.4, 1672.3}};
    double worst = 0;
    for (const auto& row : tab7) {
        AS->update(CoolProp::DmassT_INPUTS, row[1], row[0]);
        const double got = AS->viscosity(), ref = row[2] * 1e-6;
        CAPTURE(row[0], row[1], row[2]);
        REQUIRE(ValidNumber(got));
        const double rel = std::abs(got - ref) / ref;
        worst = std::max(worst, rel);
        CHECK(rel < 1e-4);  // the table is printed to five significant figures
    }
    WARN("PropyleneGlycol vs Table 7: worst relative deviation " << worst << " over 8 points");
}

TEST_CASE("PropyleneGlycol: PropsSI viscosity now works", "[expression]") {
    // Threw before this correlation shipped.
    const double eta = CoolProp::PropsSI("V", "T", 350.0, "Dmass", 1000.0, "PropyleneGlycol");
    REQUIRE(ValidNumber(eta));
    CHECK(eta == Catch::Approx(5135.986461e-6).epsilon(1e-6));
}

// Sotiriadou, Ntonti, Assael, Antoniadis & Huber, Int. J. Thermophys. 45(9):123 (2024),
// "Correlations for the Viscosity and Thermal Conductivity of Tetrahydrofuran".
// Viscosity only; the paper's thermal conductivity is not implemented here.
//
// Third rational-polynomial dilute term, and a residual (Eq. 8) whose bracket is a
// polynomial PLUS a ratio of polynomials:
//   Delta_eta = (rhor^(2/3) Tr^(1/2)) { f0 + f1 Tr + f2 rhor
//                                       + (f3 + f4 Tr)/(f5 + f6 rhor + rhor^2) }
// The initial-density term reuses the same universal Vogel/Bich B*(T*) as ethene,
// xenon and propylene glycol -- that nine-term quarter-power form now serves four
// fluids and has never needed a code change.
//
// TOLERANCE: the checks are against the paper's TABULATED values, and the bound for
// each point is derived from that point's own printed precision rather than a
// blanket epsilon.  Both columns matter and neither dominates everywhere:
//   * eta is printed to 3-5 significant figures, so eta_vap = 10.1 at 350 K carries
//     5e-3 of rounding on its own, while eta_liq = 2021.1 carries 2.5e-5;
//   * rho is printed to 5-6 figures, but dln(eta)/dln(rho) reaches 26 in the
//     saturated liquid at 200 K, which turns 5e-6 of density rounding into 1.3e-4.
// So the tolerance is halfulp(eta) + the model's own response to halfulp(rho),
// evaluated point by point.  A single global epsilon would be far too loose for the
// dense liquid and too tight for the dilute vapour.  All 35 points land inside it,
// the worst at 0.94 of its bound.
namespace {
struct ThfRow
{
    double T, rho, rho_halfulp, eta, eta_halfulp;  // K, kg/m^3, muPa.s
};
// Table 10 (0.1 / 10 / 25 MPa; the paper omits eta above 25 MPa, its validated
// limit being 30 MPa) followed by Table 9's saturation boundary, liquid and vapour.
const ThfRow THF_TABLE[] = {{200, 986.60, 0.005, 2021.1, 0.05},
                            {250, 933.91, 0.005, 826.0, 0.05},
                            {300, 879.97, 0.005, 452.9, 0.05},
                            {350, 2.5418, 0.00005, 10.04, 0.005},
                            {400, 2.2048, 0.00005, 11.56, 0.005},
                            {450, 1.9497, 0.00005, 13.04, 0.005},
                            {500, 1.7489, 0.00005, 14.49, 0.005},
                            {200, 991.26, 0.005, 2290.9, 0.05},
                            {250, 940.08, 0.005, 920.6, 0.05},
                            {300, 888.43, 0.005, 504.8, 0.05},
                            {350, 834.98, 0.005, 322.3, 0.05},
                            {400, 778.12, 0.005, 224.1, 0.05},
                            {450, 715.04, 0.005, 163.1, 0.05},
                            {500, 639.54, 0.005, 119.8, 0.05},
                            {200, 997.97, 0.005, 2768.9, 0.05},
                            {250, 948.79, 0.005, 1080.1, 0.05},
                            {300, 899.99, 0.005, 589.3, 0.05},
                            {350, 850.70, 0.005, 379.5, 0.05},
                            {400, 800.14, 0.005, 269.1, 0.05},
                            {450, 747.42, 0.005, 203.0, 0.05},
                            {500, 691.46, 0.005, 159.4, 0.05},
                            {200, 986.55, 0.005, 2018.5, 0.05},
                            {200, 0.00091066, 0.000000005, 5.57, 0.005},
                            {250, 933.85, 0.005, 825.0, 0.05},
                            {250, 0.055411, 0.0000005, 6.96, 0.005},
                            {300, 879.90, 0.005, 452.5, 0.05},
                            {300, 0.68415, 0.000005, 8.44, 0.005},
                            {350, 822.94, 0.005, 286.1, 0.05},
                            {350, 3.6651, 0.00005, 10.1, 0.05},
                            {400, 760.49, 0.005, 195.1, 0.05},
                            {400, 12.424, 0.0005, 12.1, 0.05},
                            {450, 687.67, 0.005, 137.4, 0.05},
                            {450, 32.984, 0.0005, 14.9, 0.05},
                            {500, 589.08, 0.005, 93.8, 0.05},
                            {500, 81.621, 0.0005, 19.7, 0.05}};
}  // namespace

TEST_CASE("Tetrahydrofuran: shipped viscosity matches the paper's tables", "[expression][golden]") {
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Tetrahydrofuran"));
    auto eta_at = [&](double T, double rho) {
        AS->update(CoolProp::DmassT_INPUTS, rho, T);
        return static_cast<double>(AS->viscosity()) * 1e6;  // muPa.s
    };
    // Section 4.2 computer-verification points, whose densities are exact.
    CHECK(eta_at(300.0, 1e-10) == Catch::Approx(8.3705).epsilon(1e-5));
    CHECK(eta_at(300.0, 900.0) == Catch::Approx(589.3956).epsilon(1e-6));

    double worst_frac = 0;
    int checks = 0;
    for (const auto& r : THF_TABLE) {
        const double got = eta_at(r.T, r.rho);
        // Propagate BOTH printed columns: the entry's own half-ulp, plus the model's
        // response to a half-ulp of density.  See the note above for why a single
        // epsilon cannot serve both the dense liquid and the dilute vapour.
        const double drho = 0.5 * std::abs(eta_at(r.T, r.rho + r.rho_halfulp) - eta_at(r.T, std::max(r.rho - r.rho_halfulp, 1e-12)));
        const double tol = r.eta_halfulp + drho;
        CAPTURE(r.T, r.rho, r.eta, tol);
        REQUIRE(ValidNumber(got));
        CHECK(std::abs(got - r.eta) <= tol);
        worst_frac = std::max(worst_frac, std::abs(got - r.eta) / tol);
        ++checks;
    }
    CHECK(checks == 35);
    WARN("THF vs Tables 9+10: worst deviation is " << worst_frac << " of the tabulated precision, over " << checks << " points");
}

TEST_CASE("Tetrahydrofuran: PropsSI viscosity now works", "[expression]") {
    const double eta = CoolProp::PropsSI("V", "T", 300.0, "Dmass", 900.0, "Tetrahydrofuran");
    REQUIRE(ValidNumber(eta));
    CHECK(eta == Catch::Approx(589.3956e-6).epsilon(1e-6));
}

// Polychroniadou, Antoniadis, Assael & Bell, Int. J. Thermophys. 43(1):6 (2022),
// "A Reference Correlation for the Viscosity of Krypton From Entropy Scaling".
//
// A different FAMILY of correlation, not another polynomial.  There is no
// dilute/initial-density/residual decomposition: the residual reduces on the
// residual entropy s+ = -Smolar_residual/R and on Theta2 = B2 + T dB2/dT, with a
// Lennard-Jones-derived scaled viscosity,
//
//   eta = rho_N^(2/3) sqrt(m kB T) / (s+)^(2/3) * (1.05 eta+_LJres + eta+_(rho->0))
//   eta+_LJres = exp( sum_i d_i (s+)^i ) - 1                          (Eqs. 12-13)
//
// This is what the three new DSL inputs -- Smolar_residual, Bvirial and
// dBvirial_dT -- were added for.  They are genuine state functions, so unlike the
// critical-point inputs that were tried and removed they carry no configuration
// dependence.  Bvirial and dBvirial_dT are temperature-only, verified bit-identical
// across 14 orders of magnitude in density, so reading them at the current state
// reproduces the reference implementation's rho = 1e-10 evaluation exactly.
//
// Stage split follows REFPROP's: eta_0 in `dilute`, and `higher_order` carries the
// whole entropy-scaling expression MINUS eta_0, recomputing eta_0 internally, so
// each stage still reports only its own contribution.
//
// THREE THINGS that each silently change the answer, all pinned below:
//  * R is the EOS's 8.314472, NOT CODATA.  N_A = R/kB is how the correlation
//    defines it; CODATA's 8.31446261815324 shifts the answer.
//  * The paper's printed Eq. 13 reads "-1 + sum d_i (s+)^i", omitting the exp that
//    its own reference implementation (Fig. 8) and REFPROP both apply.  Fig. 6
//    plots the quantity on a log axis rising to ~10, which settles it.
//  * Table 2 and Fig. 8 give the twelve a_i to DIFFERENT precision: Table 2 rounds
//    to 7 significant figures, Fig. 8's script carries 8 (9.1297123e-1, not
//    9.129712e-1).  Table 3 was generated by that script, so the 8-digit values are
//    the ones to ship -- they reproduce Table 3 to 6e-16 where Table 2's reproduce
//    it only to 1.6e-8.  This was checked, not assumed: sweeping the eighth digit
//    0-9 across all twelve coefficients gives a sharp minimum at 3 (6e-16) with
//    ~1e-8 either side, so the digit is real and not a rendering artifact.
//    REFPROP's KRYPTON.FLD carries Table 2's 7-digit values, and so reproduces
//    Table 3 only to ~1e-8.
TEST_CASE("Krypton: shipped entropy-scaling viscosity matches the paper's Table 3", "[expression][golden]") {
    // Table 3 is quoted to 17 digits precisely so an implementation can be checked.
    const double tab3[5][3] = {{200.0, 1e-6, 17.33865170451214},
                               {200.0, 13020.0, 56.4476422453026},
                               {298.15, 1e-6, 25.306200000810886},
                               {400.0, 1e-6, 32.795558620965195},
                               {400.0, 13020.0, 64.8014771396677}};
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Krypton"));
    double worst = 0;
    for (const auto& row : tab3) {
        AS->update(CoolProp::DmolarT_INPUTS, row[1], row[0]);
        const double got = static_cast<double>(AS->viscosity()) * 1e6, ref = row[2];
        CAPTURE(row[0], row[1], ref);
        REQUIRE(ValidNumber(got));
        const double rel = std::abs(got - ref) / ref;
        worst = std::max(worst, rel);
        CHECK(rel < 1e-12);
    }
    WARN("Krypton vs Table 3: worst relative deviation " << worst << " over 5 points");
}

TEST_CASE("Krypton: finite and positive everywhere the state is stable", "[expression][golden]") {
    // An entropy-scaling residual is exp(sum_i d_i (s+)^i), so it CAN overflow where
    // s+ is large -- and s+ is large on a mechanically unstable EOS root.  Sweeping
    // (T, rho) through the two-phase dome finds +inf, but that asks the correlation a
    // meaningless question: inside the dome a single (T, rho) root is not a state the
    // fluid can occupy, and no transport correlation is defined there.
    //
    // So the domain that must hold is the STABLE one, and it is pinned here the way a
    // caller actually addresses it: by (p, T), which always lands on the stable root,
    // plus both saturation branches, which bound it.
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Krypton"));
    const double Tt = AS->trivial_keyed_output(CoolProp::iT_triple);
    const double Tc = AS->trivial_keyed_output(CoolProp::iT_critical);
    int checked = 0;
    for (int i = 0; i < 40; ++i) {
        const double T = Tt + (750.0 - Tt) * i / 39.0;
        for (int j = 0; j < 40; ++j) {
            const double p = 1.0e3 * std::pow(2.0e5 / 1.0e3, j / 39.0);  // 1 kPa .. 200 MPa
            try {
                AS->update(CoolProp::PT_INPUTS, p, T);
            } catch (...) {
                continue;  // outside the EOS range; not this test's business
            }
            const double eta = AS->viscosity();
            CAPTURE(T, p, eta);
            REQUIRE(ValidNumber(eta));
            REQUIRE(eta > 0);
            ++checked;
        }
    }
    for (int i = 1; i < 100; ++i) {
        const double T = Tt + (Tc - Tt) * i / 100.0;
        for (double Q : {0.0, 1.0}) {
            try {
                AS->update(CoolProp::QT_INPUTS, Q, T);
            } catch (...) {
                continue;
            }
            const double eta = AS->viscosity();
            CAPTURE(T, Q, eta);
            REQUIRE(ValidNumber(eta));
            REQUIRE(eta > 0);
            ++checked;
        }
    }
    CHECK(checked > 1000);
}

TEST_CASE("Krypton: the entropy-scaling inputs are what the block asks for", "[expression]") {
    using namespace CoolProp::expression;
    // The residual block is the first thing in the tree to need EOS-derived state
    // functions beyond p.  Pin that it really does read them.
    ExpressionBlock probe(R"JSON({"formula": "Smolar_residual + Bvirial + dBvirial_dT",
  "state_variables": ["Smolar_residual", "Bvirial", "dBvirial_dT"]
  })JSON");
    CHECK(probe.required_inputs() == std::vector<std::string>{"Smolar_residual", "Bvirial", "dBvirial_dT"});
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Krypton"));
    AS->update(CoolProp::DmolarT_INPUTS, 13020.0, 200.0);
    CHECK(
      probe.evaluate(*AS)
      == Catch::Approx(AS->keyed_output(CoolProp::iSmolar_residual) + AS->keyed_output(CoolProp::iBvirial) + AS->keyed_output(CoolProp::idBvirial_dT))
           .epsilon(1e-14));
    // Bvirial and dBvirial_dT are temperature-only; the reference implementation
    // evaluates them at rho = 1e-10, and reading them at the state must agree.
    ExpressionBlock bv(R"JSON({"formula": "Bvirial",
  "state_variables": ["Bvirial"]
  })JSON");
    const double at_dense = bv.evaluate(*AS);
    AS->update(CoolProp::DmolarT_INPUTS, 1e-10, 200.0);
    CHECK(bv.evaluate(*AS) == at_dense);
}

TEST_CASE("Krypton: PropsSI viscosity now works", "[expression]") {
    const double eta = CoolProp::PropsSI("V", "T", 400.0, "Dmolar", 13020.0, "Krypton");
    REQUIRE(ValidNumber(eta));
    CHECK(eta == Catch::Approx(64.8014771396677e-6).epsilon(1e-7));
}

#endif  // ENABLE_CATCH
