// Regression test for the binary hydrogen mixture models of Beckmueller, Thol, Bell, Lemmon & Span,
// "New Equations of State for Binary Hydrogen-Mixtures Containing Methane, Nitrogen, Carbon Monoxide,
// and Carbon Dioxide", J. Phys. Chem. Ref. Data 50, 013102 (2021); GitHub issue #2263.
//
// The models add binary-specific departure functions (GERG-2008 "Gaussian+Exponential" form) and
// updated reducing parameters for H2 + {CH4, N2, CO, CO2}.  This test reproduces the single-phase
// computer-implementation test values from Table S6 of the Supplementary Material (x_H2 = 0.4).
//
// Tolerance: the paper's reference values use the GERG-2008 gas constant R = 8.314472 J/mol/K,
// whereas CoolProp's mixture gas constant defaults to the CODATA value R = 8.314462618 J/mol/K
// (config NORMALIZE_GAS_CONSTANTS; the alternative is a mole-fraction-weighted average of the pure
// fluids' R).  That ~1.1 ppm difference in R, plus the limited digits of the published parameters,
// gives a residual agreement of ~1-4 ppm, which this test checks at 2e-5.

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/Configuration.h"
#    include "CoolProp/DataStructures.h"
#    include "CoolProp/Exceptions.h"

#    include <memory>
#    include <string>
#    include <vector>

using namespace CoolProp;

namespace {
struct S6Point
{
    std::string mix;
    double T, rho;  // T [K], rho [mol/m3]; mixture = X(0.6) + H2(0.4)
    double p_MPa;   // pressure [MPa]
    double cp;      // isobaric heat capacity [J/mol/K]
    double w;       // speed of sound [m/s]
    double h;       // enthalpy [J/mol]
    double s;       // entropy [J/mol/K]
    double a;       // Helmholtz energy [J/mol]
};

// Table S6 (Supplementary Material), x_H2 = 0.4 -- full set of tabulated properties.
const std::vector<S6Point> S6 = {
  {"Methane&Hydrogen", 150, 500, 0.600540651, 31.76549083, 398.9153393, 7137.631101, 75.97586506, -5459.82996},
  {"Methane&Hydrogen", 150, 30000, 68.07281093, 40.6797545, 1249.339824, 4383.493622, 24.11332328, -1502.598568},
  {"Methane&Hydrogen", 250, 20000, 55.76246684, 43.32736532, 971.765658, 8529.65422, 47.72321056, -6189.271762},
  {"Methane&Hydrogen", 400, 2000, 6.727960274, 37.47427531, 662.5067436, 15234.47136, 87.50804965, -23132.72864},
  {"Nitrogen&Hydrogen", 90, 500, 0.35662153, 28.53557196, 243.6864251, 2425.964703, 119.0220232, -8999.260448},
  {"Nitrogen&Hydrogen", 90, 30000, 28.47064604, 41.35984046, 833.4991366, -55.46034725, 66.66850703, -7004.647515},
  {"Nitrogen&Hydrogen", 250, 20000, 65.7920363, 35.75318056, 759.9801809, 6660.746484, 100.2488922, -21691.07839},
  {"Nitrogen&Hydrogen", 400, 2000, 6.864178328, 30.05258955, 535.9935405, 11305.91635, 136.7045991, -46808.01243},
  {"CarbonMonoxide&Hydrogen", 100, 500, 0.398173739, 28.35657729, 257.1835191, 4928.015238, 73.30014651, -3198.34689},
  {"CarbonMonoxide&Hydrogen", 100, 30000, 26.39999072, 44.10138546, 743.0344713, 2251.190618, 22.26652008, -855.4610815},
  {"CarbonMonoxide&Hydrogen", 250, 20000, 62.38126662, 37.09831231, 727.5895874, 8628.010993, 52.51656521, -7620.19364},
  {"CarbonMonoxide&Hydrogen", 400, 2000, 6.854666627, 30.18519282, 535.5716536, 13528.69233, 89.01922705, -25506.33181},
  {"CarbonDioxide&Hydrogen", 260, 500, 1.04702001, 34.87609934, 322.0620575, 15033.56108, 96.29236811, -12096.49465},
  {"CarbonDioxide&Hydrogen", 260, 28000, 84.9744746, 52.24665504, 800.5687864, 9911.106009, 45.40690374, -4929.491628},
  {"CarbonDioxide&Hydrogen", 350, 20000, 69.11206089, 52.36675624, 638.5496277, 14555.74799, 62.91398482, -10919.74974},
  {"CarbonDioxide&Hydrogen", 400, 2000, 6.488944258, 39.89263375, 400.6053202, 19576.10784, 95.52311039, -21877.60844},
};
}  // namespace

TEST_CASE("Beckmueller-2021 H2 binary mixtures reproduce Table S6", "[mixtures][hydrogen][2263]") {
    // p, cp and w are reference-state-independent.  h, s and a carry the reference-state offset of
    // the pure-fluid equations (Table S1), which CoolProp reproduces because it uses those same
    // reference states -- so they match too.  The tolerance (2e-5 relative) accommodates the
    // gas-constant convention (see file header); h additionally gets a small absolute margin for the
    // one point where the tabulated enthalpy is near zero (N2+H2, 90 K, h = -55.46 J/mol).
    for (const auto& pt : S6) {
        CAPTURE(pt.mix, pt.T, pt.rho);
        std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", pt.mix));
        AS->set_mole_fractions({0.6, 0.4});  // X = 0.6, H2 = 0.4
        AS->update(DmolarT_INPUTS, pt.rho, pt.T);
        CHECK(AS->p() / 1e6 == Catch::Approx(pt.p_MPa).epsilon(2e-5));
        CHECK(AS->cpmolar() == Catch::Approx(pt.cp).epsilon(2e-5));
        CHECK(AS->speed_sound() == Catch::Approx(pt.w).epsilon(2e-5));
        CHECK(AS->hmolar() == Catch::Approx(pt.h).epsilon(2e-5).margin(0.1));
        CHECK(AS->smolar() == Catch::Approx(pt.s).epsilon(2e-5));
        CHECK(AS->helmholtzmolar() == Catch::Approx(pt.a).epsilon(2e-5));
    }
}

TEST_CASE("Beckmueller-2021 is the default model for its H2 pairs", "[mixtures][hydrogen][2263]") {
    // The superseded Kunz-JCED-2012 records for these four pairs are still present in
    // dev/mixtures/mixture_binary_pairs.json -- they are not deleted, they are kept as alternates
    // behind the Beckmueller records.  Which of the two is in force is decided purely by document
    // order: the loader keeps index 0, the FIRST record for a pair, as the model in force.  These
    // checks fail if the Beckmueller records are ever moved below the Kunz ones, which would
    // otherwise silently revert every number in the Table S6 case above.
    struct Pair
    {
        std::string CAS1, CAS2;  // sorted CAS order, as the library keys the map
        std::string function;    // departure function the default record must point at
        double gammaT;           // Beckmueller reducing parameter (unaffected by CAS swapping)
    };
    const std::vector<Pair> pairs = {
      {"1333-74-0", "74-82-8", "Methane-Hydrogen-Beckmueller2021", 1.44},
      {"1333-74-0", "7727-37-9", "Nitrogen-Hydrogen-Beckmueller2021", 1.24},
      {"124-38-9", "1333-74-0", "CarbonDioxide-Hydrogen-Beckmueller2021", 1.961},
      {"1333-74-0", "630-08-0", "CarbonMonoxide-Hydrogen-Beckmueller2021", 1.105},
    };
    for (const auto& pr : pairs) {
        CAPTURE(pr.CAS1, pr.CAS2);
        CHECK(get_mixture_binary_pair_data(pr.CAS1, pr.CAS2, "BibTeX") == "Beckmueller-JPCRD-2021");
        CHECK(get_mixture_binary_pair_data(pr.CAS1, pr.CAS2, "function") == pr.function);
        CHECK(std::stod(get_mixture_binary_pair_data(pr.CAS1, pr.CAS2, "F")) == Catch::Approx(1.0));
        CHECK(std::stod(get_mixture_binary_pair_data(pr.CAS1, pr.CAS2, "gammaT")) == Catch::Approx(pr.gammaT));
    }
}

TEST_CASE("Run-time binary parameters for a known pair are refused, not filed away", "[mixtures][hydrogen][2263]") {
    // Guards the fail-open half of allowing several models per pair.  The library shipped with
    // CoolProp may carry more than one record for a pair, with the first one in force; parameters
    // handed in at run time must NOT get that treatment, because becoming an unread alternate looks
    // exactly like success to the caller while changing nothing.  With OVERWRITE_BINARY_INTERACTION
    // off -- its default -- supplying a pair that is already known has to throw.
    const std::string CAS1 = "74-82-8", CAS2 = "7727-37-9";  // Methane & Nitrogen, sorted CAS order
    const std::string before = get_mixture_binary_pair_data(CAS1, CAS2, "gammaT");

    // A record for that very pair, carrying a gammaT that differs from the shipped one so that a
    // silent acceptance would be visible.
    const std::string doc = R"([{"CAS1": "74-82-8", "CAS2": "7727-37-9",)"
                            R"( "Name1": "Methane", "Name2": "Nitrogen",)"
                            R"( "BibTeX": "test-not-a-real-model", "F": 0.0,)"
                            R"( "betaT": 1.0, "betaV": 1.0, "gammaT": 1.234, "gammaV": 1.234}])";

    const bool overwrite_was = get_config_bool(OVERWRITE_BINARY_INTERACTION);
    set_config_bool(OVERWRITE_BINARY_INTERACTION, false);
    CHECK_THROWS_AS(set_interaction_parameters(doc), CoolProp::ValueError);
    set_config_bool(OVERWRITE_BINARY_INTERACTION, overwrite_was);

    // ...and the refusal left the pair exactly as it was.
    CHECK(get_mixture_binary_pair_data(CAS1, CAS2, "gammaT") == before);
}

#endif
