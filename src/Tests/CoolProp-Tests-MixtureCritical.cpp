// Critical point of predefined HEOS mixtures through T_critical()/p_critical()/rhomolar_critical().
//
// The L1*=0 tracer in calc_all_critical_points() also returns unstable critical points
// (typically at liquid-like density and negative pressure); the scalar getters must
// skip those and return the single stable one (GitHub #3528, R454B).
//
// The [slow] sweep runs every predefined refrigerant blend and requires exactly one
// stable critical point.  Blends that cannot meet that today are listed explicitly,
// each with the Linear issue that tracks it, and the sweep checks they still fail in
// the listed way -- so a fix shows up as a test failure asking for the list to shrink,
// and a newly broken blend can never hide in a skip.  The exception is the list of
// blends with a marginal stability verdict, which only checks the getters' contract.
//
// Built into CatchTestRunner when COOLPROP_CATCH_MODULE=ON.  Run via
//   ./build_catch/CatchTestRunner "[mixture_critical]"

#include "CoolProp/AbstractState.h"
#include "CoolProp/CoolProp.h"
#include "CoolProp/DataStructures.h"
#include "CPstrings.h"

#if defined(ENABLE_CATCH)
#    include <catch2/catch_all.hpp>

#    include <algorithm>
#    include <memory>
#    include <set>
#    include <string>
#    include <vector>

namespace {

std::size_t count_stable(const std::vector<CoolProp::CriticalState>& pts) {
    return static_cast<std::size_t>(std::count_if(pts.begin(), pts.end(), [](const CoolProp::CriticalState& c) { return c.stable; }));
}

}  // namespace

TEST_CASE("R454B critical point skips the unstable negative-pressure point (#3528)", "[mixture_critical]") {
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "R454B.mix"));
    const std::vector<CoolProp::CriticalState> pts = AS->all_critical_points();
    // The tracer finds the real point and a spurious one; only the real one is stable
    REQUIRE(pts.size() >= 2);
    REQUIRE(count_stable(pts) == 1);

    // Regression values from the tracer itself (not reference data): the stable point
    // reported in #3528, 351.43 K / 5.304 MPa / 7054.7 mol/m^3
    CHECK(AS->T_critical() == Catch::Approx(351.43).epsilon(1e-4));
    CHECK(AS->p_critical() == Catch::Approx(5.304e6).epsilon(1e-3));
    CHECK(AS->rhomolar_critical() == Catch::Approx(7054.7).epsilon(1e-3));

    // The high-level route reported as failing in #3528
    CHECK(CoolProp::PropsSI("Pcrit", "", 0, "", 0, "R454B.mix") == Catch::Approx(AS->p_critical()).epsilon(1e-12));
}

TEST_CASE("Critical point of a mixture throws when there is more than one stable point", "[mixture_critical]") {
    // R452A gives three stable points within 1 K of each other (COO-154); the getters
    // must refuse to pick one rather than silently return the first
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "R452A.mix"));
    REQUIRE(count_stable(AS->all_critical_points()) > 1);
    CHECK_THROWS_WITH(AS->p_critical(), Catch::Matchers::ContainsSubstring("stable critical points"));
}

TEST_CASE("Every predefined refrigerant blend has exactly one stable critical point", "[mixture_critical][slow]") {
    // Blends whose constructor throws: a binary pair has no interaction parameters (COO-157)
    const std::set<std::string> unbuildable = {"R401A.mix", "R401B.mix", "R401C.mix", "R402A.mix", "R402B.mix", "R403A.mix", "R403B.mix", "R405A.mix",
                                               "R406A.mix", "R408A.mix", "R409A.mix", "R409B.mix", "R412A.mix", "R413A.mix", "R414A.mix", "R414B.mix",
                                               "R416A.mix", "R424A.mix", "R426A.mix", "R429A.mix", "R435A.mix", "R437A.mix", "R438A.mix", "R446A.mix",
                                               "R453A.mix", "R458A.mix", "R461A.mix", "R468A.mix", "R468B.mix", "R468C.mix", "R473A.mix"};
    // Blends with several stable critical points within 1 K of each other (COO-154)
    const std::set<std::string> several_stable = {"R452A.mix"};
    // Blends whose tracer also finds a positive-pressure point at ~70-105 K and liquid-like
    // density, whose TPD verdict differs between builds (COO-155: R470A's is stable in the
    // 8.0.1dev wheel, unstable on master).  The stable count is not portable, so only the
    // getters' contract is checked: throw on several stable points, else return the one.
    // Listed whether the point is stable on master (first five) or not (the rest).
    const std::set<std::string> marginal_stability = {"R407H.mix", "R448B.mix", "R449C.mix", "R460A.mix", "R470A.mix", "R407E.mix",
                                                      "R407F.mix", "R407G.mix", "R407I.mix", "R460B.mix", "R501.mix",  "R503.mix"};
    // Blends whose critical-point search throws (COO-156)
    const std::set<std::string> search_throws = {"R452C.mix"};

    std::vector<std::string> names;
    for (const std::string& name : strsplit(CoolProp::get_global_param_string("predefined_mixtures"), ',')) {
        // Each blend is listed as both .MIX and .mix; refrigerants are the R-numbered ones
        if (name.size() > 4 && name[0] == 'R' && name.compare(name.size() - 4, 4, ".mix") == 0) {
            names.push_back(name);
        }
    }
    // Guard against the filter (or the list) silently going empty
    REQUIRE(names.size() > 100);
    // Every listed exception must be a blend the sweep visits, so no entry goes stale unseen
    for (const auto* list : {&unbuildable, &several_stable, &marginal_stability, &search_throws}) {
        for (const std::string& listed : *list) {
            CAPTURE(listed);
            CHECK(std::find(names.begin(), names.end(), listed) != names.end());
        }
    }

    for (const std::string& name : names) {
        CAPTURE(name);
        std::shared_ptr<CoolProp::AbstractState> AS;
        if (unbuildable.count(name) > 0) {
            CHECK_THROWS_WITH(AS.reset(CoolProp::AbstractState::factory("HEOS", name)),
                              Catch::Matchers::ContainsSubstring("binary pair"));  // now buildable? remove it from the list
            continue;
        }
        CHECK_NOTHROW(AS.reset(CoolProp::AbstractState::factory("HEOS", name)));
        if (!AS) {
            continue;
        }
        if (search_throws.count(name) > 0) {
            CHECK_THROWS_WITH(AS->all_critical_points(), Catch::Matchers::ContainsSubstring("bracket"));  // now succeeds? remove it from the list
            continue;
        }
        std::vector<CoolProp::CriticalState> pts;
        try {
            pts = AS->all_critical_points();
        } catch (const std::exception& e) {
            FAIL_CHECK("all_critical_points() threw: " << e.what());
            continue;
        }
        const std::size_t n_stable = count_stable(pts);
        CAPTURE(pts.size(), n_stable);
        if (several_stable.count(name) > 0) {
            CHECK(n_stable > 1);  // now unique? remove it from the list
            CHECK_THROWS_WITH(AS->p_critical(), Catch::Matchers::ContainsSubstring("stable critical points"));
            continue;
        }
        if (marginal_stability.count(name) > 0 && n_stable > 1) {
            CHECK_THROWS_WITH(AS->p_critical(), Catch::Matchers::ContainsSubstring("stable critical points"));
            continue;
        }
        CHECK(n_stable == 1);
        if (n_stable != 1) {
            continue;
        }
        const auto it = std::find_if(pts.begin(), pts.end(), [](const CoolProp::CriticalState& c) { return c.stable; });
        CHECK(AS->T_critical() == Catch::Approx(it->T).epsilon(1e-12));
        CHECK(AS->p_critical() == Catch::Approx(it->p).epsilon(1e-12));
        CHECK(AS->rhomolar_critical() == Catch::Approx(it->rhomolar).epsilon(1e-12));
    }
}

#endif  // ENABLE_CATCH
