// The per-backend cache of the reducing-function composition derivatives (HelmholtzEOSMixtureBackend::get_reducing_derivs,
// COO-124): dT_r/dx_j, drho_r/dx_j, n*dT_r/dn_i, n*drho_r/dn_i and the two Hessians d(n*dY_r/dn_i)/dx_j.  MixtureDerivatives
// reads these at ~20 call sites per fugacity vector; the cache fills them in one pass (ReducingFunction::
// composition_derivatives) once per (reducing-function version, x_N flag, composition).  What has to hold: every cached
// entry is the reducing function's own per-element value, and no change to the composition or to the reducing function
// -- a binary-parameter setter, on the state or through a linked state, or a replaced reducing function -- is ever
// answered from a stale entry.
//
// "The same value" is to 1e-11 of the largest entry of that field, not bit for bit.  The one-pass fill evaluates the
// same expressions on the same operands in the same order, and at most compositions it IS bit-identical; but the
// GERG terms it calls (dYrdxi__constxj, d2vrmolardxidxj) get inlined into its loops, where the compiler makes different
// FMA-contraction choices than in the out-of-line per-element functions -- the effect documented in
// CoolProp-Tests-DeltaOnly.cpp.  Measured over the random trace compositions below: at most 7.7e-13 relative to the
// entry itself (a 0.59 entry left by cancellation), 3.6e-12 absolute, and below 4e-14 of the field's largest entry.  A
// formula or indexing error moves entries by O(1) of that scale, so 1e-11 still catches it.

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/DataStructures.h"
#    include "../Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#    include "../Backends/Helmholtz/ReducingFunctions.h"

#    include <algorithm>
#    include <array>
#    include <cmath>
#    include <cstdint>
#    include <cstring>
#    include <memory>
#    include <random>
#    include <string>
#    include <vector>

using namespace CoolProp;

namespace {

std::array<unsigned char, sizeof(CoolPropDbl)> bits_of(CoolPropDbl d) {
    std::array<unsigned char, sizeof(CoolPropDbl)> u{};
    std::memcpy(u.data(), &d, sizeof(CoolPropDbl));
    return u;
}

// Every cached entry of `heos` against the reducing function's per-element functions (to 1e-11 of the field's largest
// entry, see the header).  The first call asks
// for the first derivatives only, the second adds the Hessians to the same slot, so both fill paths are covered.
// Returns the cached n*d(rho_r)/dn_i so a caller can check that a change actually moved them.
std::vector<CoolPropDbl> check_against_direct(HelmholtzEOSMixtureBackend& heos, x_N_dependency_flag flag) {
    const std::vector<CoolPropDbl> x = heos.get_mole_fractions();
    const std::size_t N = x.size();
    const ReducingFunction& R = *heos.Reducing;
    // `scale`: the largest |entry| of the field being compared, so cancellation-small entries are judged on the field's scale
    auto same = [](CoolPropDbl cached, CoolPropDbl direct, CoolPropDbl scale) {
        CHECK(std::isfinite(direct));
        CHECK(std::abs(cached - direct) <= 1e-11 * scale);
    };
    auto vec_scale = [](const std::vector<CoolPropDbl>& v) {
        CoolPropDbl m = 0;
        for (CoolPropDbl e : v) {
            m = std::max(m, std::abs(e));
        }
        return m;
    };
    auto mat_scale = [&](const STLMatrix& M) {
        CoolPropDbl m = 0;
        for (const auto& row : M) {
            m = std::max(m, vec_scale(row));
        }
        return m;
    };
    const auto& c1 = heos.get_reducing_derivs(flag);
    REQUIRE(c1.dTr.size() == N);
    REQUIRE(c1.drhor.size() == N);
    REQUIRE(c1.ndTr.size() == N);
    REQUIRE(c1.ndrhor.size() == N);
    // NaN in a cached entry would make the comparisons below false, but a scale of NaN would not -- check it first
    const CoolPropDbl s_dTr = vec_scale(c1.dTr), s_drhor = vec_scale(c1.drhor), s_ndTr = vec_scale(c1.ndTr), s_ndrhor = vec_scale(c1.ndrhor);
    REQUIRE((std::isfinite(s_dTr) && std::isfinite(s_drhor) && std::isfinite(s_ndTr) && std::isfinite(s_ndrhor)));
    for (std::size_t i = 0; i < N; ++i) {
        CAPTURE(i, static_cast<int>(flag));
        same(c1.dTr[i], R.dTrdxi__constxj(x, i, flag), s_dTr);
        same(c1.drhor[i], R.drhormolardxi__constxj(x, i, flag), s_drhor);
        same(c1.ndTr[i], R.ndTrdni__constnj(x, i, flag), s_ndTr);
        same(c1.ndrhor[i], R.ndrhorbardni__constnj(x, i, flag), s_ndrhor);
    }
    const std::vector<CoolPropDbl> ndrhor = c1.ndrhor;
    const auto& c2 = heos.get_reducing_derivs(flag, /*with_hessian=*/true);
    REQUIRE(c2.d_ndTr_dxj.size() == N);
    REQUIRE(c2.d_ndrhor_dxj.size() == N);
    for (std::size_t i = 0; i < N; ++i) {
        REQUIRE(c2.d_ndTr_dxj[i].size() == N);
        REQUIRE(c2.d_ndrhor_dxj[i].size() == N);
    }
    const CoolPropDbl s_HT = mat_scale(c2.d_ndTr_dxj), s_Hrho = mat_scale(c2.d_ndrhor_dxj);
    REQUIRE((std::isfinite(s_HT) && std::isfinite(s_Hrho)));
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            CAPTURE(i, j, static_cast<int>(flag));
            same(c2.d_ndTr_dxj[i][j], R.d_ndTrdni_dxj__constxi(x, i, j, flag), s_HT);
            same(c2.d_ndrhor_dxj[i][j], R.d_ndrhorbardni_dxj__constxi(x, i, j, flag), s_Hrho);
        }
        CHECK(bits_of(c2.ndrhor[i]) == bits_of(ndrhor[i]));  // adding the Hessians left the first derivatives alone
    }
    return ndrhor;
}

}  // namespace

TEST_CASE("Cached reducing-function composition derivatives equal the reducing function's per-element values", "[mixture][reducing_nderivs]") {
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

TEST_CASE("Cached reducing-function composition derivatives, random compositions with trace components", "[mixture][reducing_nderivs]") {
    // The flash visits compositions far from the feed (stability-test trial phases, near-pure phases), so cover many,
    // with some components at trace levels where the GERG terms cancel most
    const std::vector<std::string> mixes = {"Nitrogen&Methane&Ethane&n-Butane&n-Pentane",
                                            "Methane&Nitrogen&CarbonDioxide&Ethane&Propane&IsoButane&n-Butane&Isopentane&n-Pentane&n-Hexane"};
    std::mt19937 g(7);
    std::uniform_real_distribution<double> U(0, 1);
    for (const auto& fluids : mixes) {
        std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", fluids));
        auto* heos = dynamic_cast<HelmholtzEOSMixtureBackend*>(AS.get());
        REQUIRE(heos != nullptr);
        const std::size_t N = heos->get_mole_fractions().size() ? heos->get_mole_fractions().size() : heos->get_components().size();
        for (int trial = 0; trial < 20; ++trial) {
            std::vector<double> z(N);
            double sum = 0;
            for (auto& zi : z) {
                zi = U(g) < 0.3 ? std::pow(10.0, -8 * U(g)) * 1e-4 : U(g);  // ~30 % trace components
                sum += zi;
            }
            for (auto& zi : z) {
                zi /= sum;
            }
            CAPTURE(fluids, z);
            AS->set_mole_fractions(z);
            check_against_direct(*heos, XN_DEPENDENT);
            check_against_direct(*heos, XN_INDEPENDENT);
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
    CHECK(heos->reducing_derivs[1].version == v1);  // the XN_DEPENDENT slot
}

#endif
