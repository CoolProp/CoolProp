// PT flash inside the two-phase dome of NEAR-AZEOTROPIC refrigerant blends (COO-193).
//
// For a near-azeotrope the incipient phase is close to the feed and its tangent-plane distance is shallow.
// The stability test's successive-substitution proximity heuristic ("trial converging to the trivial
// solution") also fires on such a trial.  When that heuristic alone was taken as proof of stability
// (#3427), the PT flash reported single-phase for most dome states of R410A and R507A.  The heuristic is now
// trusted only for a trial on the feed's own density branch, and a trial stuck on that branch is re-run on the
// other one (with the Chebyshev density solver on; see the second test case).  Each state here is
// put on the dome with a QT flash, and the PT flash at that same (T, p) must find the split with the same
// quality.  The PT flash only reaches the stability test because no phase envelope is built here (a built,
// closed envelope bypasses it in PT_flash_mixtures).

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/Configuration.h"
#    include "CoolProp/DataStructures.h"
#    include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"

#    include <cmath>
#    include <memory>
#    include <string>
#    include <vector>

using namespace CoolProp;

namespace {

struct Blend
{
    const char* name;
    const char* fluids;
    std::vector<double> mass_fractions;
    double T_min, T_max;  // [K]
};

// Switches the Chebyshev density solver on for one test and restores the previous setting
struct ChebyshevSolverOn
{
    const bool previous = get_config_bool(CHEBYSHEV_DENSITY_SOLVER);
    ChebyshevSolverOn() {
        set_config_bool(CHEBYSHEV_DENSITY_SOLVER, true);
    }
    ~ChebyshevSolverOn() {
        set_config_bool(CHEBYSHEV_DENSITY_SOLVER, previous);
    }
    ChebyshevSolverOn(const ChebyshevSolverOn&) = delete;
    ChebyshevSolverOn& operator=(const ChebyshevSolverOn&) = delete;
    ChebyshevSolverOn(ChebyshevSolverOn&&) = delete;
    ChebyshevSolverOn& operator=(ChebyshevSolverOn&&) = delete;
};

struct DomeIsotherms
{
    const char* name;
    const char* fluids;
    std::vector<double> mole_fractions;  // empty: the predefined mixture's own
    std::vector<double> T;               // [K]
};

}  // namespace

TEST_CASE("Near-azeotrope: PT flash finds the split at every QT dome state", "[flash][mixture]") {
    // T ranges stay below the blends' critical regions.  R507A below 235 K is left out: those states need the
    // opposite-branch trial, i.e. the Chebyshev density solver (COO-195; covered by the next test case).
    const std::vector<Blend> blends = {
      {"R410A", "R32&R125", {0.5, 0.5}, 215.0, 305.0},
      {"R454B", "R32&R1234yf", {0.689, 0.311}, 215.0, 305.0},
      {"R404A", "R125&R143a&R134a", {0.44, 0.52, 0.04}, 215.0, 305.0},
      {"R507A", "R125&R143a", {0.5, 0.5}, 235.0, 305.0},
    };
    for (const auto& b : blends) {
        DYNAMIC_SECTION(b.name) {
            std::unique_ptr<AbstractState> qt(AbstractState::factory("HEOS", b.fluids));
            std::unique_ptr<AbstractState> pt(AbstractState::factory("HEOS", b.fluids));
            qt->set_mass_fractions(b.mass_fractions);
            pt->set_mass_fractions(b.mass_fractions);
            int n = 0;
            for (int k = 0; b.T_min + 10.0 * k <= b.T_max + 1e-9; ++k) {
                const double T = b.T_min + 10.0 * k;
                for (double Q : {0.05, 0.25, 0.5, 0.75, 0.95}) {
                    CAPTURE(T, Q);
                    qt->update(QT_INPUTS, Q, T);
                    REQUIRE(qt->phase() == iphase_twophase);
                    const double p = qt->p();
                    CAPTURE(p);
                    pt->update(PT_INPUTS, p, T);
                    CHECK(pt->phase() == iphase_twophase);
                    if (pt->phase() == iphase_twophase) {
                        CHECK(std::abs(pt->Q() - Q) < 1e-3);
                    }
                    ++n;
                }
            }
            CHECK(n >= 40);  // the grid was not empty
        }
    }
}

TEST_CASE("PT flash: incipient phase on the branch opposite to the feed's (COO-193)", "[flash][mixture]") {
    // The Wilson trial for the other phase starts so close to the feed that its lowest-Gibbs root is on the
    // FEED's density branch, and successive substitution there collapses to the trivial solution: a liquid feed's
    // vapor trial (CO2/ethane, R440A, R451B at low Q) or a vapor feed's liquid trial (R451B, R507A at high Q).  The
    // stability test then re-runs that trial on the other branch, with roots from the Chebyshev solver, so these
    // states are recovered only with that solver on.  Isotherms from the 203-mixture sweep of COO-193; R507A at
    // 189-225 K and R513A at 315 K are the states of COO-195.
    const ChebyshevSolverOn cheb;
    const std::vector<DomeIsotherms> cases = {
      {"CO2/ethane 20/80", "CarbonDioxide&Ethane", {0.2, 0.8}, {177.69, 207.30, 222.11}},
      {"CO2/ethane 80/20", "CarbonDioxide&Ethane", {0.8, 0.2}, {205.79}},
      {"R440A", "R440A.mix", {}, {231.43, 289.28}},
      {"R451B", "R451B.mix", {}, {312.45, 330.83}},
      {"R507A", "R507A.mix", {}, {189.45, 206.67, 223.89, 225.0}},
      {"R513A", "R513A.mix", {}, {315.0}},
    };
    for (const auto& c : cases) {
        DYNAMIC_SECTION(c.name) {
            std::unique_ptr<AbstractState> qt(AbstractState::factory("HEOS", c.fluids));
            std::unique_ptr<AbstractState> pt(AbstractState::factory("HEOS", c.fluids));
            if (!c.mole_fractions.empty()) {
                qt->set_mole_fractions(c.mole_fractions);
                pt->set_mole_fractions(c.mole_fractions);
            }
            int n = 0;
            for (double T : c.T) {
                for (double Q : {0.05, 0.25, 0.5, 0.75, 0.95}) {
                    CAPTURE(T, Q);
                    qt->update(QT_INPUTS, Q, T);
                    REQUIRE(qt->phase() == iphase_twophase);
                    const double p = qt->p();
                    CAPTURE(p);
                    pt->update(PT_INPUTS, p, T);
                    CHECK(pt->phase() == iphase_twophase);
                    if (pt->phase() == iphase_twophase) {
                        CHECK(std::abs(pt->Q() - Q) < 1e-3);
                    }
                    ++n;
                }
            }
            CHECK(n == 5 * static_cast<int>(c.T.size()));
        }
    }
}

TEST_CASE("PT flash: no opposite-branch split from a spurious density root (COO-193)", "[flash][mixture]") {
    // A phase-imposed Newton solve for the vapor-branch root of this compressed liquid lands in an alpha^r well
    // (a mechanically stable root at the right pressure, tm = -1e12) and once published a split into pure R32 and
    // pure R1234yf.  The Chebyshev solver certifies there is no vapor-branch root, so the state stays liquid.
    const ChebyshevSolverOn cheb;
    std::unique_ptr<AbstractState> pt(AbstractState::factory("HEOS", "R454B.mix"));
    shared_ptr<AbstractState> probe(AbstractState::factory("HEOS", "R454B.mix"));
    auto& heos = dynamic_cast<HelmholtzEOSMixtureBackend&>(*probe);
    for (const auto& Tp :
         std::vector<std::pair<double, double>>{{246.374596, 20.4031194e6}, {219.699301, 16.5283263e6}, {277.840023, 26.6381821e6}}) {
        CAPTURE(Tp.first, Tp.second);
        // the certified path is the one exercised: the feed has no vapor-branch root here
        REQUIRE(heos.solver_rho_Tp_cheb_branch(Tp.first, Tp.second, true) == HelmholtzEOSMixtureBackend::CHEB_NO_ROOT_ON_BRANCH);
        pt->update(PT_INPUTS, Tp.second, Tp.first);
        CHECK(pt->phase() != iphase_twophase);
        CHECK(pt->rhomolar() > 15000);
    }
}

TEST_CASE("Chebyshev solver: root on a requested spinodal branch (COO-193)", "[flash][mixture]") {
    const ChebyshevSolverOn cheb;
    shared_ptr<AbstractState> as(AbstractState::factory("HEOS", "CarbonDioxide&Ethane"));
    auto& heos = dynamic_cast<HelmholtzEOSMixtureBackend&>(*as);
    heos.set_mole_fractions(std::vector<CoolPropDbl>{0.2, 0.8});
    // Inside the dome both branches have a root; regression values from the phase-imposed solves
    const double T = 207.30, p = 395233.33;
    const CoolPropDbl rV = heos.solver_rho_Tp_cheb_branch(T, p, true), rL = heos.solver_rho_Tp_cheb_branch(T, p, false);
    CHECK(rV == Catch::Approx(251.45).epsilon(1e-4));
    CHECK(rL == Catch::Approx(18130.37).epsilon(1e-5));
    // The stable root is the liquid one (the feed's lower Gibbs energy)
    CHECK(heos.solver_rho_Tp_cheb(T, p) == Catch::Approx(rL).epsilon(1e-12));
    // A compressed liquid well below the critical temperature has no vapor-branch root at 20 MPa
    CHECK(heos.solver_rho_Tp_cheb_branch(220.0, 20e6, true) == HelmholtzEOSMixtureBackend::CHEB_NO_ROOT_ON_BRANCH);
}

#endif
