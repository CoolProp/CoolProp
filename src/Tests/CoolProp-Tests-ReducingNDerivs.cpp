// The per-backend cache of the reducing-function composition derivatives n*d(rho_r)/dn_i and n*d(T_r)/dn_i
// (HelmholtzEOSMixtureBackend::get_reducing_nderivs, COO-124).  MixtureDerivatives reads these at ~20 call sites per
// fugacity vector; the cache computes them once per (reducing-function version, x_N flag, composition).  What has to
// hold: the cached entries are the reducing function's own values bit for bit, and no change to the composition or to
// the reducing function -- a binary-parameter setter, on the state or through a linked state, or a replaced reducing
// function -- is ever answered from a stale entry.

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/DataStructures.h"
#    include "../Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#    include "../Backends/Helmholtz/ReducingFunctions.h"

#    include <array>
#    include <cmath>
#    include <cstdint>
#    include <cstring>
#    include <memory>
#    include <string>
#    include <vector>

using namespace CoolProp;

namespace {

std::array<unsigned char, sizeof(CoolPropDbl)> bits_of(CoolPropDbl d) {
    std::array<unsigned char, sizeof(CoolPropDbl)> u{};
    std::memcpy(u.data(), &d, sizeof(CoolPropDbl));
    return u;
}

// The cached vectors of `heos` against the reducing function evaluated directly, entry by entry and bit for bit.
// Returns the cached n*d(rho_r)/dn_i so a caller can check that a change actually moved them.
std::vector<CoolPropDbl> check_against_direct(HelmholtzEOSMixtureBackend& heos, x_N_dependency_flag flag) {
    const std::vector<CoolPropDbl> x = heos.get_mole_fractions();
    const auto& c = heos.get_reducing_nderivs(flag);
    REQUIRE(c.ndrhor.size() == x.size());
    REQUIRE(c.ndTr.size() == x.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
        CAPTURE(i, static_cast<int>(flag));
        const CoolPropDbl rho_direct = heos.Reducing->ndrhorbardni__constnj(x, i, flag);
        const CoolPropDbl T_direct = heos.Reducing->ndTrdni__constnj(x, i, flag);
        CHECK(std::isfinite(rho_direct));
        CHECK(std::isfinite(T_direct));
        CHECK(bits_of(c.ndrhor[i]) == bits_of(rho_direct));
        CHECK(bits_of(c.ndTr[i]) == bits_of(T_direct));
    }
    return c.ndrhor;
}

}  // namespace

TEST_CASE("Cached reducing-function composition derivatives equal the reducing function's, bit for bit", "[mixture][reducing_nderivs]") {
    struct Case
    {
        std::string backend, fluids;
        std::vector<std::vector<double>> zs;
        double T, p;  // a single-phase state for every composition, so the fugacity coefficients are defined
    };
    const std::vector<Case> cases = {
      {"GERG2008",
       "Methane&Ethane&Propane&n-Butane&Nitrogen",
       {{0.80, 0.10, 0.05, 0.03, 0.02}, {0.2, 0.2, 0.2, 0.2, 0.2}, {0.05, 0.05, 0.1, 0.3, 0.5}},
       350,
       1e6},
      {"HEOS", "R32&R1234yf", {{0.8292, 0.1708}, {0.3, 0.7}}, 350, 0.5e6},
      {"HEOS", "CarbonDioxide&Water", {{0.99, 0.01}, {0.5, 0.5}}, 450, 0.2e6},
    };
    for (const auto& c : cases) {
        std::shared_ptr<AbstractState> AS(AbstractState::factory(c.backend, c.fluids));
        auto* heos = dynamic_cast<HelmholtzEOSMixtureBackend*>(AS.get());
        REQUIRE(heos != nullptr);
        // The compositions run on ONE backend, so a cache that ignored the composition would serve the previous one's
        // entries and fail the bit comparison.
        for (const auto& z : c.zs) {
            CAPTURE(c.fluids, z);
            AS->set_mole_fractions(z);
            check_against_direct(*heos, XN_DEPENDENT);
            check_against_direct(*heos, XN_INDEPENDENT);
            // and the same through the flash, which is where MixtureDerivatives reads them
            AS->update(PT_INPUTS, c.p, c.T);
            REQUIRE(AS->phase() != iphase_twophase);
            for (std::size_t i = 0; i < z.size(); ++i) {
                CHECK(std::isfinite(AS->fugacity_coefficient(i)));
            }
            check_against_direct(*heos, XN_DEPENDENT);
        }
    }
}

TEST_CASE("Cached reducing-function composition derivatives follow binary-parameter changes, on linked states too", "[mixture][reducing_nderivs]") {
    std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", "Methane&Ethane"));
    auto* heos = dynamic_cast<HelmholtzEOSMixtureBackend*>(AS.get());
    REQUIRE(heos != nullptr);
    AS->set_mole_fractions(std::vector<double>{0.3, 0.7});
    auto* satL = heos->SatL.get();
    REQUIRE(satL != nullptr);
    satL->set_mole_fractions(std::vector<CoolPropDbl>{0.3, 0.7});

    // Fill the caches of the state and of a linked state, then change a parameter through the state.
    const std::vector<CoolPropDbl> before = check_against_direct(*heos, XN_DEPENDENT);
    const std::vector<CoolPropDbl> before_L = check_against_direct(*satL, XN_DEPENDENT);
    const double betaV = AS->get_binary_interaction_double(0, 1, "betaV");
    AS->set_binary_interaction_double(0, 1, "betaV", betaV * 1.05);

    const std::vector<CoolPropDbl> after = check_against_direct(*heos, XN_DEPENDENT);
    const std::vector<CoolPropDbl> after_L = check_against_direct(*satL, XN_DEPENDENT);
    // The change must have reached the linked state (it holds its own copy of the reducing function) and moved the
    // derivatives on both, or the checks above prove nothing about invalidation.
    CHECK(satL->get_binary_interaction_double(0, 1, "betaV") == betaV * 1.05);
    CHECK(bits_of(after[0]) != bits_of(before[0]));
    CHECK(bits_of(after_L[0]) != bits_of(before_L[0]));
}

TEST_CASE("A replaced reducing function never matches a cached entry", "[mixture][reducing_nderivs]") {
    // Caches key on ReducingFunction::get_version() alone, so it must be unique per object and per parameter set: a backend
    // whose Reducing is replaced (set_mixture_parameters, sync_linked_states via copy()) can get a new object at the
    // freed address, which a pointer-plus-per-object-counter key would mistake for the old one.
    std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", "Methane&Ethane"));
    auto* heos = dynamic_cast<HelmholtzEOSMixtureBackend*>(AS.get());
    REQUIRE(heos != nullptr);
    AS->set_mole_fractions(std::vector<double>{0.3, 0.7});
    check_against_direct(*heos, XN_DEPENDENT);

    const std::uint64_t v0 = heos->Reducing->get_version();
    std::shared_ptr<ReducingFunction> replacement(heos->Reducing->copy());
    CHECK(replacement->get_version() != v0);
    replacement->set_binary_interaction_double(0, 1, "betaT", 1.07);
    const std::uint64_t v1 = replacement->get_version();
    CHECK(v1 != v0);
    heos->Reducing = replacement;
    check_against_direct(*heos, XN_DEPENDENT);
    CHECK(heos->get_reducing_nderivs(XN_DEPENDENT).version == v1);
}

#endif
