// Regression guard for GH #3287: the cubic (PR/SRK) entropy *value* must be
// self-consistent with CoolProp's own analytic derivatives.
//
// Root cause of #3287: HelmholtzEOSMixtureBackend::calc_all_alpha0_derivs_nocache
// (added in PR #2625) returned the ideal-gas alpha0 derivatives from
// E.alpha0.all(taustar, ...) WITHOUT the chain-rule scaling (Tc/Tr for tau
// derivatives, rhor/rhomolarc for delta derivatives) that the old
// calc_alpha0_deriv_nocache and the mixture branch both apply.  calc_smolar()
// consumed the unscaled da0/dtau, so cubic Smolar/Smass VALUES were too steep in
// T (~1.5x methane .. ~3x water), while Hmolar and d(Smolar)/d(T)|P stayed
// correct.  HEOS is unaffected because Tc/Tr == 1 there.
//
// Two internal self-consistency checks (no external reference data, so these are
// independent of the entropy reference state):
//   1. Smass(T2) - Smass(T1) == integral of the analytic d(Smass)/d(T)|P
//   2. T * (dS/dT|P of the value) == Cpmass  (from T*(dS/dT)_P = (dH/dT)_P = Cp)

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include "CoolProp.h"

#    include <string>
#    include <vector>

namespace {

// Trapezoid integral of the ANALYTIC derivative d(Smass)/d(T)|P over [T1, T2].
double integ_dSdT(const std::string& fluid, double P, double T1, double T2, int n = 2000) {
    const double h = (T2 - T1) / n;
    double s = 0;
    for (int i = 0; i <= n; ++i) {
        const double w = (i == 0 || i == n) ? 0.5 : 1.0;
        s += w * CoolProp::PropsSI("d(Smass)/d(T)|P", "T", T1 + i * h, "P", P, fluid);
    }
    return h * s;
}

}  // namespace

TEST_CASE("Cubic entropy value is self-consistent with its own derivative (GH #3287)", "[cubic][entropy][3287]") {
    // Single-phase, superheated, well away from saturation for each fluid.
    struct Case
    {
        std::string fluid;
        double P, T1, T2;
    };
    const std::vector<Case> cases = {
      {"Water", 1e5, 450.0, 550.0},
      {"Methane", 1e6, 200.0, 300.0},
      {"Ammonia", 1e5, 300.0, 400.0},
    };

    for (const std::string& be : {std::string("PR"), std::string("SRK")}) {
        for (const Case& c : cases) {
            const std::string f = be + "::" + c.fluid;
            CAPTURE(f, c.P, c.T1, c.T2);

            // Check 1: value increment integrates to the analytic derivative.
            const double dS_value = CoolProp::PropsSI("Smass", "T", c.T2, "P", c.P, f) - CoolProp::PropsSI("Smass", "T", c.T1, "P", c.P, f);
            const double dS_integral = integ_dSdT(f, c.P, c.T1, c.T2);
            // 0.2 % of the analytic increment covers trapezoid truncation with room to spare;
            // the #3287 defect made these differ by 50 %..200 %.
            CHECK(dS_value == Catch::Approx(dS_integral).epsilon(2e-3));

            // Check 2: T*(dS/dT)_P from a central difference of the VALUE equals Cp.
            const double Tm = 0.5 * (c.T1 + c.T2);
            const double dT = 0.5;
            const double dSdT =
              (CoolProp::PropsSI("Smass", "T", Tm + dT, "P", c.P, f) - CoolProp::PropsSI("Smass", "T", Tm - dT, "P", c.P, f)) / (2 * dT);
            const double cp = CoolProp::PropsSI("Cpmass", "T", Tm, "P", c.P, f);
            CHECK(Tm * dSdT == Catch::Approx(cp).epsilon(2e-3));
        }
    }
}

#endif  // ENABLE_CATCH
