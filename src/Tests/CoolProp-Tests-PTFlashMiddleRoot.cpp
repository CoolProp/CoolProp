// Focused regression test for the blind (no-envelope) mixture PT flash returning a
// spurious middle-branch density root (GitHub #3283).
//
// Runs in the default suite (tag [PTflash], NOT [.]-hidden).  Run explicitly:
//   ./CatchTestRunner "[PTflash]"
//
// Background (#3283): for Methane&Ethane&Propane at z=[0.05,0.90,0.05], p=2.0 MPa the
// state is subcooled single-phase LIQUID (bubble p ~ 8.98 bar << 20 bar).  At scattered
// "knife-edge" temperatures the flash returned a spurious rho ~ 6751 mol/m^3 classified
// as gas (h ~ -666 kJ/mol) instead of the liquid root (rho ~ 16400, h ~ +3 kJ/mol).
//
// Root cause: in FlashRoutines::PT_flash_mixtures the SRK pre-screen has a single real
// root at these conditions (solver_rho_Tp_SRK returns the same value for both phase
// requests), but solve_single_phase treated gas_ok && liq_ok as "two roots", solved a
// spurious gas HEOS branch seeded from the liquid-like SRK density (converging onto a
// secondary loop of the GERG isotherm), and the Gibbs-minimising selection picked that
// artifact root because its unphysical Gibbs energy was lower than the true liquid's.
// The fix only Gibbs-compares when the SRK roots are genuinely distinct.

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include <memory>
#    include <string>
#    include <vector>

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/DataStructures.h"

TEST_CASE("Mixture PT flash selects the liquid root for subcooled liquid (was #3283)", "[PTflash]") {
    // p = 2.0 MPa >> bubble pressure (~0.9 MPa) at all these T -> single-phase liquid.
    // Includes the reported knife-edge temperatures (223.15, 212.5, 235 K) that failed,
    // plus their passing neighbours as controls.
    const double T = GENERATE(212.5, 222.5, 223.15, 225.0, 235.0, 248.15, 193.15);
    CAPTURE(T);
    const double p = 2.0e6;

    // Reference: the correct liquid root, obtained by explicitly imposing the liquid phase.
    auto ref = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Methane&Ethane&Propane"));
    ref->set_mole_fractions({0.05, 0.90, 0.05});
    ref->specify_phase(CoolProp::iphase_liquid);
    ref->update(CoolProp::PT_INPUTS, p, T);
    const double rho_liquid = ref->rhomolar();

    // Sanity on the reference: a real ethane-rich liquid near these T is ~15000-18000 mol/m^3.
    REQUIRE(rho_liquid > 14000.0);
    REQUIRE(rho_liquid < 19000.0);

    // The blind flash (no phase imposed, fresh state) must land on the same liquid root,
    // not the spurious middle-branch root (~6751 mol/m^3).
    auto wrk = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Methane&Ethane&Propane"));
    wrk->set_mole_fractions({0.05, 0.90, 0.05});
    wrk->update(CoolProp::PT_INPUTS, p, T);

    // Loose epsilon: this compares two different code paths (blind SRK-seeded solve vs the
    // imposed-liquid reference); the spurious root is ~4x away, so 1e-4 guards it with headroom.
    CHECK(wrk->rhomolar() == Catch::Approx(rho_liquid).epsilon(1e-4));
    CHECK(wrk->phase() == CoolProp::iphase_liquid);
    // Guard the specific failure signature: the spurious root sat well below the liquid branch.
    CHECK(wrk->rhomolar() > 13000.0);
}

#endif
