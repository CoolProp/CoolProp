// VTPR (Volume-Translated Peng-Robinson, group-contribution) MIXTURE tests.
//
// CoolProp ships the VTPR machinery but not the UNIFAC parameters, so these tests embed a minimal
// UNIFAC dataset (methanol + water: their two subgroups, R_k/Q_k, the CH3OH<->H2O interaction, and
// the group decompositions), write it to a temporary directory, point VTPR_UNIFAC_PATH at it, and
// run a mixture flash.  This is the first mixture-level coverage for the VTPR backend and guards the
// UNIFAC evaluation path (the caching + dense-array work in UNIFAC.cpp): the asserted bubble
// pressure is a regression value that changes if the UNIFAC group-activity math changes.

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/Configuration.h"
#    include "CoolProp/DataStructures.h"

#    include <filesystem>
#    include <fstream>
#    include <memory>
#    include <string>
#    include <vector>

using namespace CoolProp;

namespace {

const char* const GROUP_DATA = R"([
 {"sgi": 15, "mgi": 6, "R_k": 0.0, "Q_k": 0.8779, "maingroup_name": "CH3OH", "subgroup_name": "CH3OH"},
 {"sgi": 16, "mgi": 7, "R_k": 0.0, "Q_k": 1.5576, "maingroup_name": "H2O",   "subgroup_name": "H2O"}
])";

const char* const INTERACTION_DATA = R"([
 {"mgi1": 6, "mgi2": 7, "a_ij": -387.4, "b_ij": 1.9621, "c_ij": -0.0034336,
  "a_ji": -168.82, "b_ji": 0.6674, "c_ji": 0.0041881}
])";

const char* const DECOMP_DATA = R"([
 {"name": "Methanol", "registry_number": "67-56-1", "inchikey": "OKKJLVBELUTLKV-UHFFFAOYSA-N", "userid": "",
  "Tc": 512.5, "pc": 8084000.0, "acentric": 0.5658, "molemass": 0.03204216,
  "groups": [{"count": 1, "sgi": 15}], "alpha": {"type": "Twu", "c": [0.6755, 0.9141, 1.7586]}},
 {"name": "Water", "registry_number": "7732-18-5", "inchikey": "XLYOFNOQVPJJNP-UHFFFAOYSA-N", "userid": "",
  "Tc": 647.1, "pc": 22064000.0, "acentric": 0.3449, "molemass": 0.018015268,
  "groups": [{"count": 1, "sgi": 16}], "alpha": {"type": "Twu", "c": [0.3865, 0.872, 1.9693]}}
])";

// Writes the minimal UNIFAC dataset to a fresh temp directory and returns its (slash-terminated) path.
std::string write_unifac_dataset() {
    std::filesystem::path dir = std::filesystem::temp_directory_path() / "coolprop_vtpr_test_unifac";
    std::filesystem::create_directories(dir);
    auto put = [&](const char* fname, const char* body) {
        std::ofstream f((dir / fname).string(), std::ios::binary | std::ios::trunc);
        f << body;
    };
    put("group_data.json", GROUP_DATA);
    put("interaction_parameters.json", INTERACTION_DATA);
    put("decompositions.json", DECOMP_DATA);
    std::string p = dir.string();
    if (!p.empty() && p.back() != '/' && p.back() != '\\') p += '/';
    return p;
}

// Restores VTPR config on scope exit so the test does not leak state.  The VTPR UNIFAC library is a
// process-global that this test overwrites with its minimal dataset; restoring only the config keys
// would leave that dataset resident, so a later VTPR test could pick it up.  We therefore restore the
// path but FORCE reload-on-next-use (rather than restoring the old flag), so any subsequent VTPR
// construction reloads from its own path instead of reusing this test's library.
struct VTPRConfigGuard
{
    std::string old_path = get_config_string(VTPR_UNIFAC_PATH);
    VTPRConfigGuard() = default;
    VTPRConfigGuard(const VTPRConfigGuard&) = delete;
    VTPRConfigGuard& operator=(const VTPRConfigGuard&) = delete;
    VTPRConfigGuard(VTPRConfigGuard&&) = delete;
    VTPRConfigGuard& operator=(VTPRConfigGuard&&) = delete;
    ~VTPRConfigGuard() {
        set_config_string(VTPR_UNIFAC_PATH, old_path);
        set_config_bool(VTPR_ALWAYS_RELOAD_LIBRARY, true);
    }
};

}  // namespace

TEST_CASE("VTPR methanol-water mixture flash", "[VTPR][cubic][mixture]") {
    VTPRConfigGuard guard;
    set_config_string(VTPR_UNIFAC_PATH, write_unifac_dataset());
    set_config_bool(VTPR_ALWAYS_RELOAD_LIBRARY, true);  // force our test dataset to load

    SECTION("bubble pressure regression at 313.15 K, x_MeOH = 0.5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("VTPR", "Methanol&Water"));
        AS->set_mole_fractions({0.5, 0.5});
        AS->update(QT_INPUTS, 0.0, 313.15);  // bubble point
        // Regression value from the UNIFAC group-activity evaluation (VTPR, this dataset).
        CHECK(AS->p() == Catch::Approx(24669.872).epsilon(1e-4));
        // The incipient vapour is methanol-rich (positive deviation from Raoult's law).
        CHECK(AS->mole_fractions_vapor()[0] == Catch::Approx(0.822922).epsilon(1e-4));
    }

    SECTION("blind PT flash inside the two-phase region returns two phases") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("VTPR", "Methanol&Water"));
        AS->set_mole_fractions({0.5, 0.5});
        AS->update(PT_INPUTS, 20000.0, 313.15);  // ~20 kPa is between dew and bubble here
        CHECK(AS->phase() == iphase_twophase);
        CHECK(AS->Q() > 0.0);
        CHECK(AS->Q() < 1.0);
    }

    SECTION("pure-endpoint bubble pressures bracket the mixture") {
        // Water is the heavy end, methanol the light end; the mixture bubble P sits between them.
        auto W = std::shared_ptr<AbstractState>(AbstractState::factory("VTPR", "Water"));
        W->update(QT_INPUTS, 0.0, 313.15);
        auto M = std::shared_ptr<AbstractState>(AbstractState::factory("VTPR", "Methanol"));
        M->update(QT_INPUTS, 0.0, 313.15);
        CHECK(W->p() < M->p());
        CHECK(W->p() == Catch::Approx(7453.9).epsilon(1e-3));
        CHECK(M->p() == Catch::Approx(36061.0).epsilon(1e-3));
    }
}

#endif
