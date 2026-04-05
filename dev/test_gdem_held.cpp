#include <iostream>
#include <vector>
#include "CoolProp.h"
#include "AbstractState.h"
#include "Backends/Helmholtz/VLERoutines.h"
#include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"

using namespace CoolProp;

static void check_pt(const std::string& fluids, const std::vector<double>& z,
                     double T, double p, const std::string& label) {
    try {
        auto AS = shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
        AS->set_mole_fractions(z);
        AS->update(PT_INPUTS, p, T);
        std::string ph = (AS->phase() == iphase_twophase) ? "2ph" : "1ph";
        std::cout << "[PT]   " << label << ": " << ph;
        if (AS->phase() == iphase_twophase)
            std::cout << "  Q=" << AS->Q() << "  rho=" << AS->rhomolar() << " mol/m3";
        std::cout << "\n";
    } catch (const std::exception& e) {
        std::cout << "[PT]   " << label << ": ERROR " << e.what() << "\n";
    }
}

// Initialise HEOS to (T,p,z) via PT_INPUTS so T and p are properly set internally
static HelmholtzEOSMixtureBackend* make_heos(const std::string& fluids,
                                              const std::vector<double>& z,
                                              double T, double p,
                                              shared_ptr<AbstractState>& AS) {
    AS = shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
    AS->set_mole_fractions(z);
    AS->update(PT_INPUTS, p, T);  // sets _T, _p inside HEOS
    return static_cast<HelmholtzEOSMixtureBackend*>(AS.get());
}

static void test_gdem(const std::string& fluids, const std::vector<double>& z,
                      double T, double p, const std::string& label) {
    try {
        shared_ptr<AbstractState> AS;
        HelmholtzEOSMixtureBackend& HEOS = *make_heos(fluids, z, T, p, AS);

        SaturationSolvers::GDEMSuccessiveSubstitution gdem(HEOS);
        gdem.run();
        std::cout << "[GDEM] " << label << ": "
                  << (gdem.unstable ? "UNSTABLE" : "stable")
                  << "  iters=" << gdem.iterations << "\n";
        if (gdem.unstable) {
            std::cout << "         x[0]=" << gdem.x[0] << "  y[0]=" << gdem.y[0]
                      << "  rhoL=" << gdem.rhomolar_liq << "  rhoV=" << gdem.rhomolar_vap << "\n";
        }
    } catch (const std::exception& e) {
        std::cout << "[GDEM] " << label << ": ERROR " << e.what() << "\n";
    }
}

static void test_held(const std::string& fluids, const std::vector<double>& z,
                      double T, double p, const std::string& label) {
    try {
        shared_ptr<AbstractState> AS;
        HelmholtzEOSMixtureBackend& HEOS = *make_heos(fluids, z, T, p, AS);

        StabilityRoutines::HELDStabilityClass held(HEOS);
        bool stable = held.run();
        std::cout << "[HELD] " << label << ": "
                  << (stable ? "STABLE" : "UNSTABLE")
                  << "  n_stationary=" << held.stationary_points.size() << "\n";
        for (std::size_t i = 0; i < held.stationary_points.size(); ++i) {
            const auto& sp = held.stationary_points[i];
            std::cout << "         [" << i << "] tpd=" << sp.tpd
                      << "  trivial=" << sp.is_trivial
                      << "  rho=" << sp.rhomolar << " mol/m3\n";
        }
        if (!stable) {
            std::cout << "         seed_x[0]=" << held.x_seed[0]
                      << "  seed_y[0]=" << held.y_seed[0] << "\n";
        }
    } catch (const std::exception& e) {
        std::cout << "[HELD] " << label << ": ERROR " << e.what() << "\n";
    }
}

// Use QT or PQ flash to find a confirmed two-phase (T,p) point for a mixture
static std::pair<double,double> bubble_point(const std::string& fluids,
                                              const std::vector<double>& z, double T) {
    auto AS = shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
    AS->set_mole_fractions(z);
    AS->update(QT_INPUTS, 0.0, T);  // Q=0: bubble point at T
    return {AS->T(), AS->p()};
}

int main() {
    // Helper: find Q=0.5 state (midpoint of two-phase dome) at a given T
    auto midpoint_state = [](const std::string& fluids, const std::vector<double>& z,
                              double T) -> std::pair<double,double> {
        auto AS = shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
        AS->set_mole_fractions(z);
        AS->update(QT_INPUTS, 0.5, T);  // midpoint of two-phase dome
        return {AS->T(), AS->p()};
    };

    // --- Propane/Ethane 50/50 ---
    {
        auto [T_bub, p_bub] = bubble_point("Propane&Ethane", {0.5, 0.5}, 260.0);
        auto [T_mid, p_mid] = midpoint_state("Propane&Ethane", {0.5, 0.5}, 260.0);
        std::cout << "\n=== Propane/Ethane 50/50 @ 260 K ===\n";
        std::cout << "  Bubble pressure: " << p_bub/1e5 << " bar\n";
        std::cout << "  Midpoint (Q=0.5) pressure: " << p_mid/1e5 << " bar\n";

        // Deep inside two-phase dome (Q=0.5)
        std::cout << "  [expect 2ph]\n";
        check_pt ("Propane&Ethane", {0.5,0.5}, T_mid, p_mid, "Q=0.5 midpoint");
        test_gdem("Propane&Ethane", {0.5,0.5}, T_mid, p_mid, "Q=0.5 midpoint");
        test_held("Propane&Ethane", {0.5,0.5}, T_mid, p_mid, "Q=0.5 midpoint");

        // Clearly above bubble point (compressed liquid, 1ph)
        double p_above = p_bub * 2.0;
        std::cout << "  [expect 1ph, above bubble]\n";
        check_pt ("Propane&Ethane", {0.5,0.5}, 260.0, p_above, "260K,2*p_bub");
        test_gdem("Propane&Ethane", {0.5,0.5}, 260.0, p_above, "260K,2*p_bub");
        test_held("Propane&Ethane", {0.5,0.5}, 260.0, p_above, "260K,2*p_bub");
    }

    // --- Propane/Ethane supercritical ---
    std::cout << "\n=== Propane/Ethane 50/50 supercritical @ 400 K, 50 bar (expect 1ph) ===\n";
    check_pt ("Propane&Ethane", {0.5,0.5}, 400.0, 5e6, "400K,50bar");
    test_gdem("Propane&Ethane", {0.5,0.5}, 400.0, 5e6, "400K,50bar");
    test_held("Propane&Ethane", {0.5,0.5}, 400.0, 5e6, "400K,50bar");

    // --- CO2/CH4/N2 3-component ---
    // Use CH4-rich composition to keep bubble pressure low and avoid density solver issues
    try {
        auto [T_mid, p_mid] = midpoint_state("CarbonDioxide&Methane&Nitrogen",
                                              {0.1, 0.8, 0.1}, 200.0);
        std::cout << "\n=== CO2/CH4/N2 (10/80/10) @ 200 K ===\n";
        std::cout << "  Midpoint (Q=0.5) pressure: " << p_mid/1e5 << " bar\n";
        std::cout << "  [expect 2ph]\n";
        check_pt ("CarbonDioxide&Methane&Nitrogen", {0.1,0.8,0.1}, T_mid, p_mid, "Q=0.5");
        test_gdem("CarbonDioxide&Methane&Nitrogen", {0.1,0.8,0.1}, T_mid, p_mid, "Q=0.5");
        test_held("CarbonDioxide&Methane&Nitrogen", {0.1,0.8,0.1}, T_mid, p_mid, "Q=0.5");
    } catch (const std::exception& e) {
        std::cout << "CO2/CH4/N2 midpoint: ERROR " << e.what() << "\n";
    }

    // --- CO2/CH4/N2 single-phase ---
    std::cout << "\n=== CO2/CH4/N2 equimolar @ 300 K, 200 bar (expect 1ph) ===\n";
    check_pt ("CarbonDioxide&Methane&Nitrogen", {1./3,1./3,1./3}, 300.0, 20e6, "300K,200bar");
    test_gdem("CarbonDioxide&Methane&Nitrogen", {1./3,1./3,1./3}, 300.0, 20e6, "300K,200bar");
    test_held("CarbonDioxide&Methane&Nitrogen", {1./3,1./3,1./3}, 300.0, 20e6, "300K,200bar");

    std::cout << "\nAll tests done.\n";
    return 0;
}
