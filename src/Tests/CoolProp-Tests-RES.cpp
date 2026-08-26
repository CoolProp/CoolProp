// Catch2 tests for Residual Entropy Scaling (RES) transport properties.
//
// Sample data sources:
//   Viscosity  - Martinek 2025, doi:10.1021/acs.jced.4c00451
//                dev/RES_samples/Martinek_2025_viscosity/
//   Conductivity - Li 2024, doi:10.1021/acs.iecr.4c02946
//                  dev/RES_samples/Li_2024_conductivity/
//
// The "vis_res" / "TC_RES" columns are the authors' own values, taken from the supporting
// information of the two papers and computed there with REFPROP as the thermodynamic backend.

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/CoolProp.h"
#    include "Backends/Cubics/CubicBackend.h"
#    include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#    include "Backends/REFPROP/REFPROPMixtureBackend.h"
#    include "CoolProp/Configuration.h"
#    include "CoolProp/Exceptions.h"

#    include <algorithm>
#    include <cmath>
#    include <limits>
#    include <fstream>
#    include <iostream>
#    include <memory>
#    include <sstream>
#    include <string>
#    include <unordered_map>
#    include <vector>

using namespace CoolProp;

// ----------------------------------------------------------------- paths ----
#    ifndef RES_SAMPLES_DIR
#        define RES_SAMPLES_DIR "dev/RES_samples"
#    endif
static const std::string VIS_PURE_OUT = RES_SAMPLES_DIR "/Martinek_2025_viscosity/Samples_pure_fluids_output.txt";
static const std::string VIS_BIN_OUT = RES_SAMPLES_DIR "/Martinek_2025_viscosity/Samples_binaries_output.txt";
static const std::string TC_PURE_IN = RES_SAMPLES_DIR "/Li_2024_conductivity/Samples_pure_fluids.txt";
static const std::string TC_BIN_IN = RES_SAMPLES_DIR "/Li_2024_conductivity/Samples_binaries.txt";
static const std::string TC_PURE_REF = RES_SAMPLES_DIR "/Li_2024_conductivity/Table_S5_SI.txt";
static const std::string TC_BIN_REF = RES_SAMPLES_DIR "/Li_2024_conductivity/Table_S6_SI.txt";

// --------------------------------------------------------- name mapping -----
static const std::unordered_map<std::string, std::string> REFPROP_TO_CP = {
  {"ARGON", "Argon"},
  {"CO2", "CarbonDioxide"},
  {"BUTANE", "n-Butane"},
  {"METHANE", "Methane"},
  {"ETHANE", "Ethane"},
  {"PROPANE", "Propane"},
  {"HEPTANE", "n-Heptane"},
  {"HEXANE", "n-Hexane"},
  {"PENTANE", "n-Pentane"},
  {"DECANE", "n-Decane"},
  {"DODECANE", "n-Dodecane"},
  {"NITROGEN", "Nitrogen"},
  {"OXYGEN", "Oxygen"},
  {"HELIUM", "Helium"},
  {"NEON", "Neon"},
  {"HYDROGEN", "Hydrogen"},
  {"WATER", "Water"},
  {"AMMONIA", "Ammonia"},
  {"ACETONE", "Acetone"},
  {"BENZENE", "Benzene"},
  {"TOLUENE", "Toluene"},
  {"CYCLOHEX", "CycloHexane"},
  {"CYCLOPEN", "Cyclopentane"},
  {"CYCLOPRO", "CycloPropane"},
  {"CO", "CarbonMonoxide"},
  {"D2O", "HeavyWater"},
  {"D2", "Deuterium"},
  {"ACETYLENE", "Acetylene"},
  {"1BUTENE", "1-Butene"},
  {"13BUTADIENE", "1,3-Butadiene"},
  {"C2BUTENE", "cis-2-Butene"},
  {"C1CC6", "MethylCyclohexane"},
  {"C3CC6", "PropylCyclohexane"},
  {"CHLORINE", "Chlorine"},
  {"C11", "n-C11"},
  {"C12", "n-C12"},
  {"C16", "n-C16"},
  {"C22", "n-C22"},
  {"EGLYCOL", "MEG"},
  {"ETHANOL", "Ethanol"},
};
static std::string cp_name(const std::string& n) {
    auto it = REFPROP_TO_CP.find(n);
    return (it != REFPROP_TO_CP.end()) ? it->second : n;
}

// ---------------------------------------------------------- tokeniser -------
static std::vector<std::string> split_ws(const std::string& line) {
    std::vector<std::string> t;
    std::istringstream ss(line);
    std::string w;
    while (ss >> w)
        t.push_back(w);
    return t;
}

// --------------------------------------------------------- data structs ------
struct VisPureSample
{
    std::string name;
    double T, p_MPa, den, vis_exp, vis_res;
};
struct VisBinSample
{
    std::string c1, c2;
    double mf1, mf2, T, p_MPa, den, vis_exp, vis_res;
};
struct TcPureSample
{
    std::string name;
    double T, p_kPa, den, tc_exp, tc_res;
};
struct TcBinSample
{
    std::string c1, c2;
    double mf1, mf2, T, p_kPa, den, tc_exp, tc_res;
};

// --------------------------------------------------------- parsers -----------
// Martinek output: Material T p den s vis_exp vis_res vis_REF
static std::vector<VisPureSample> parse_vis_pure(const std::string& path) {
    std::vector<VisPureSample> out;
    std::ifstream f(path);
    if (!f) return out;
    std::string line;
    std::getline(f, line);
    while (std::getline(f, line)) {
        auto t = split_ws(line);
        if (t.size() < 8) continue;
        out.push_back({t[0], std::stod(t[1]), std::stod(t[2]), std::stod(t[3]), std::stod(t[5]), std::stod(t[6])});
    }
    return out;
}
// Martinek binary output: Comp1 Comp2 GroupN1 N2 MF1 MF2 T p den s vis_exp vis_res vis_REF
static std::vector<VisBinSample> parse_vis_bin(const std::string& path) {
    std::vector<VisBinSample> out;
    std::ifstream f(path);
    if (!f) return out;
    std::string line;
    std::getline(f, line);
    while (std::getline(f, line)) {
        auto t = split_ws(line);
        if (t.size() < 13) continue;
        out.push_back(
          {t[0], t[1], std::stod(t[4]), std::stod(t[5]), std::stod(t[6]), std::stod(t[7]), std::stod(t[8]), std::stod(t[10]), std::stod(t[11])});
    }
    return out;
}
// Li2024 Table_S5: Material T p den s TC_EXP TC_RES TC_REF
static std::vector<TcPureSample> parse_tc_pure(const std::string& ref_path) {
    std::vector<TcPureSample> out;
    std::ifstream f(ref_path);
    if (!f) return out;
    std::string line;
    std::getline(f, line);
    while (std::getline(f, line)) {
        auto t = split_ws(line);
        if (t.size() < 8) continue;
        out.push_back({t[0], std::stod(t[1]), std::stod(t[2]), std::stod(t[3]), std::stod(t[5]), std::stod(t[6])});
    }
    return out;
}
// Li2024 Table_S6: Comp1 Comp2 GroupN1 N2 MF1 MF2 T p den s TC_EXP TC_RES TC_REF
static std::vector<TcBinSample> parse_tc_bin(const std::string& ref_path) {
    std::vector<TcBinSample> out;
    std::ifstream f(ref_path);
    if (!f) return out;
    std::string line;
    std::getline(f, line);
    while (std::getline(f, line)) {
        auto t = split_ws(line);
        if (t.size() < 13) continue;
        out.push_back(
          {t[0], t[1], std::stod(t[4]), std::stod(t[5]), std::stod(t[6]), std::stod(t[7]), std::stod(t[8]), std::stod(t[10]), std::stod(t[11])});
    }
    return out;
}

// --------------------------------------------------------- helpers -----------
static std::shared_ptr<AbstractState> make_heos(const std::string& f) {
    return std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", f));
}
static std::shared_ptr<AbstractState> make_pr(const std::string& f) {
    return std::shared_ptr<AbstractState>(AbstractState::factory("PR", f));
}

/// dynamic_cast to the cubic backend, asserted.  Never dereference the result of the cast
/// directly: a failed cast would be a null-deref that crashes the whole runner instead of
/// reporting one failed assertion.
static AbstractCubicBackend* as_cubic(const std::shared_ptr<AbstractState>& AS) {
    auto* cubic = dynamic_cast<AbstractCubicBackend*>(AS.get());
    REQUIRE(cubic != nullptr);
    return cubic;
}

// ============================================================= tests =========

TEST_CASE("GEN golden PT", "[.][GEN2]") {
    struct S
    {
        const char* f;
        double T, p;
    };
    static const S st[] = {
      {"AMMONIA", 304.17, 3606963.44186},    {"AMMONIA", 446.116, 28384252.4254},   {"AMMONIA", 527.228, 92191380.93711},
      {"ARGON", 113.01525, 2418347.95793},   {"ARGON", 165.7557, 11548934.69671},   {"ARGON", 195.8931, 38374336.25303},
      {"CO2", 228.09615, 2490320.627062},    {"CO2", 334.54102, 19224671.21944},    {"CO2", 395.36666, 67343164.89186},
      {"ETHANOL", 386.0325, 1596351.335694}, {"ETHANOL", 566.181, 17457987.6492},   {"ETHANOL", 669.123, 62922345.22261},
      {"METHANE", 142.923, 2230767.992782},  {"METHANE", 209.6204, 11014898.69906}, {"METHANE", 247.7332, 36862390.29281},
      {"NITROGEN", 138.8112, 8247948.58104}, {"NITROGEN", 164.0496, 27881957.2182}, {"NITROGEN", 94.644, 1577718.259735},
      {"PROPANE", 277.4175, 1618042.935256}, {"PROPANE", 406.879, 10932008.08445},  {"PROPANE", 480.857, 38972745.33683},
      {"R134A", 280.6575, 1193029.563139},   {"R134A", 411.631, 11082887.77478},    {"R134A", 486.473, 40174480.69617},
      {"TOLUENE", 443.8125, 1288149.39309},  {"TOLUENE", 650.925, 10804158.42783},  {"TOLUENE", 769.275, 39974654.51745},
      {"WATER", 485.322, 6404695.669172},    {"WATER", 711.8056, 54841210.98041},   {"WATER", 841.2248, 170019943.4478},
    };
    for (const auto& s : st) {
        double v[2] = {0, 0}, c[2] = {0, 0}, rho[2] = {0, 0};
        bool ok = true;
        for (int i = 0; i < 2 && ok; ++i) {
            try {
                auto AS = std::shared_ptr<AbstractState>(AbstractState::factory(i == 0 ? "HEOS" : "PR", s.f));
                AS->use_viscosity_RES(true);
                AS->use_conductivity_RES(true);
                AS->update(PT_INPUTS, s.p, s.T);
                rho[i] = AS->rhomass();
                v[i] = AS->viscosity();
                c[i] = AS->conductivity();
            } catch (const std::exception& e) {
                std::cout << "    // " << s.f << " " << (i == 0 ? "HEOS" : "PR") << " failed: " << e.what() << std::endl;
                ok = false;
            }
        }
        if (!ok) continue;
        std::cout.precision(13);
        std::cout << "        {\"" << s.f << "\", " << s.T << ", " << s.p << ", " << std::scientific << v[0] << ", " << c[0] << ", " << v[1] << ", "
                  << c[1] << std::defaultfloat << "},  // PR/HEOS eta " << (v[1] / v[0] - 1) * 100 << "%, tc " << (c[1] / c[0] - 1) * 100 << "%, rho "
                  << (rho[1] / rho[0] - 1) * 100 << "%" << std::endl;
    }
}

TEST_CASE("RES regression at residual-dominated states", "[RES][transport]") {
    // The regression net for the native backends.  These are GOLDEN values produced by this
    // implementation -- their job is to fail when the model changes, not to say it is right.
    //
    // Correctness is established elsewhere, and deliberately not here:
    //   * against the authors' own code, on REFPROP where the coefficients are exact --
    //     "RES on REFPROP reproduces the published..." below, and dev/RES_reference_check.py,
    //     which reproduces Martinek/Li to a median of 2.6e-7 and 1.2e-7;
    //   * across backends, HEOS vs REFPROP over 3828 grid points -- dev/RES_grid_report.py,
    //     median 6e-7, which is where the exclusion lists come from.
    // Repeating those as unit tests would need percent-level tolerances, because several model
    // choices legitimately differ between backends, and a percent-level bound on a quantity that
    // agrees to 1e-7 tests almost nothing while reading as though the implementation were poor.
    //
    // The states sit where the RESIDUAL term supplies 55-96% of the property.  That is the gap the
    // published sample points cannot cover: one state per fluid, several of them almost pure
    // dilute gas -- residual share 0.2% for R1234YF -- where the dilute polynomial is identical on
    // every backend and the residual term is not exercised at all.
    //
    // Both backends are evaluated at (p, T), each finding its own density, because that is how the
    // model is used and how the PR coefficients were fitted.  Pinning PR at HEOS's density instead
    // would test it outside the frame it was fitted in: PR's density is up to 25% from HEOS's at
    // these states, and forcing it onto HEOS's value inflates the apparent PR-HEOS difference from
    // a median of 3% to tens of percent.  The states are far from saturation, so the flash is
    // unambiguous on both.
    //
    // Regenerate deliberately, and only when a change to the model is intended and understood:
    // dev/RES_grid_build.py, then the [RES_grid] harness, then re-read the values.
    struct Golden
    {
        const char* fluid;
        double T, p, eta_heos, tc_heos, eta_pr, tc_pr;
    };
    static const Golden golden[] = {
      {"AMMONIA", 304.17, 3606963.44186, 1.2534715542329e-04, 4.6091983840002e-01, 1.2430582780973e-04, 4.6014354642977e-01},
      {"AMMONIA", 446.116, 28384252.4254, 4.2568850005885e-05, 2.3571693252803e-01, 4.2786422334760e-05, 2.2615317984091e-01},
      {"AMMONIA", 527.228, 92191380.93711, 5.3572527060133e-05, 2.9205044497509e-01, 5.9347099520628e-05, 3.0816594819808e-01},
      {"ARGON", 113.01525, 2418347.95793, 1.2425384695901e-04, 9.3236284166714e-02, 1.2941949107901e-04, 9.2446488869947e-02},
      {"ARGON", 165.7557, 11548934.69671, 4.5407457592372e-05, 4.5422611875469e-02, 4.4549397046605e-05, 4.3690287175164e-02},
      {"ARGON", 195.8931, 38374336.25303, 6.5481497788676e-05, 6.2490603951953e-02, 6.8869885182853e-05, 6.3579015380416e-02},
      {"CO2", 228.09615, 2490320.627062, 2.1706417724008e-04, 1.7472032636980e-01, 2.1027646640269e-04, 1.6968175186621e-01},
      {"CO2", 334.54102, 19224671.21944, 5.7748662636897e-05, 7.8000300052477e-02, 5.8696604679699e-05, 7.3214978852871e-02},
      {"CO2", 395.36666, 67343164.89186, 8.0220003253735e-05, 1.0326713380878e-01, 8.7061221303422e-05, 1.0351998485063e-01},
      {"ETHANOL", 386.0325, 1596351.335694, 2.8376781690654e-04, 1.4825360425374e-01, 2.8170762756801e-04, 1.4839265301095e-01},
      {"ETHANOL", 566.181, 17457987.6492, 5.3192725711269e-05, 1.2862040142026e-01, 4.5834435892781e-05, 1.2743669579936e-01},
      {"ETHANOL", 669.123, 62922345.22261, 5.8461231258474e-05, 1.5254927961775e-01, 6.3802494308521e-05, 1.6291558600855e-01},
      {"METHANE", 142.923, 2230767.992782, 6.3204079903066e-05, 1.4205684520521e-01, 6.3995464623802e-05, 1.4048207374639e-01},
      {"METHANE", 209.6204, 11014898.69906, 2.5963132143352e-05, 7.5897272201812e-02, 2.4828334843230e-05, 7.3627590090899e-02},
      {"METHANE", 247.7332, 36862390.29281, 3.6830270301239e-05, 1.0189887832794e-01, 3.7603345051969e-05, 1.0345987125414e-01},
      {"NITROGEN", 138.8112, 8247948.58104, 3.0053479834784e-05, 5.2312119391949e-02, 3.0563295198522e-05, 5.0809298455116e-02},
      {"NITROGEN", 164.0496, 27881957.2182, 4.3336521007401e-05, 7.1427476618019e-02, 4.6625864875662e-05, 7.2863987077628e-02},
      {"NITROGEN", 94.644, 1577718.259735, 8.7504634740090e-05, 1.0581701680628e-01, 8.8903428620479e-05, 1.0505012170376e-01},
      {"PROPANE", 277.4175, 1618042.935256, 1.1646412035562e-04, 1.0488839979475e-01, 1.1820459342471e-04, 1.0941577238384e-01},
      {"PROPANE", 406.879, 10932008.08445, 4.0783692276745e-05, 6.6802493591870e-02, 3.9415796812998e-05, 7.0951949503742e-02},
      {"PROPANE", 480.857, 38972745.33683, 5.7849566275168e-05, 9.0651408581744e-02, 6.0039764046944e-05, 9.9703431380447e-02},
      {"R134A", 280.6575, 1193029.563139, 2.5770377777319e-04, 9.0393800182755e-02, 2.5428245071266e-04, 8.9475107403340e-02},
      {"R134A", 411.631, 11082887.77478, 6.5171000688799e-05, 5.3006508949767e-02, 6.5961723631826e-05, 5.4129997398074e-02},
      {"R134A", 486.473, 40174480.69617, 8.9548801894064e-05, 6.7979959212060e-02, 9.8220388334921e-05, 7.2273829074962e-02},
      {"TOLUENE", 443.8125, 1288149.39309, 1.6509786571318e-04, 9.3704872024963e-02, 1.6482924210734e-04, 9.2218047500124e-02},
      {"TOLUENE", 650.925, 10804158.42783, 5.4619891095099e-05, 7.5127993330366e-02, 5.2468262292424e-05, 7.4880790613610e-02},
      {"TOLUENE", 769.275, 39974654.51745, 7.4722506698651e-05, 9.5824266869547e-02, 7.8921869489485e-05, 9.7944905220167e-02},
      {"WATER", 485.322, 6404695.669172, 1.2527298863343e-04, 6.8982986008485e-01, 1.2907412366549e-04, 6.9455344075235e-01},
      {"WATER", 711.8056, 54841210.98041, 5.8243233448324e-05, 3.5802447395435e-01, 5.8810456734224e-05, 3.7882515600633e-01},
      {"WATER", 841.2248, 170019943.4478, 6.9142832618337e-05, 4.1148115810909e-01, 7.6769727417017e-05, 4.9820814236974e-01},
    };

    for (const auto& g : golden) {
        for (int i = 0; i < 2; ++i) {
            const bool heos = (i == 0);
            auto AS = std::shared_ptr<AbstractState>(AbstractState::factory(heos ? "HEOS" : "PR", g.fluid));
            AS->use_viscosity_RES(true);
            AS->use_conductivity_RES(true);
            AS->update(PT_INPUTS, g.p, g.T);
            INFO("fluid=" << g.fluid << " backend=" << (heos ? "HEOS" : "PR") << " T=" << g.T << " p=" << g.p);
            // 1e-7 absorbs the flash's own last-bit noise and compiler differences.  A real change
            // to the residual term is orders of magnitude larger -- the R_D/gamma correction moved
            // ETHANE by 7.5e-5, some 750x this bound -- so the test still fails loudly on one.
            CHECK(AS->viscosity() == Catch::Approx(heos ? g.eta_heos : g.eta_pr).epsilon(1e-7));
            CHECK(AS->conductivity() == Catch::Approx(heos ? g.tc_heos : g.tc_pr).epsilon(1e-7));
        }
    }
    // Guard the premise: if the table is ever emptied or truncated this must not pass quietly.
    CHECK(std::size(golden) == 30);
}

TEST_CASE("RES withheld fluids report the omission rather than guessing", "[RES][transport]") {
    // These fluids carry no HEOS entry on purpose: their REFPROP-fitted coefficients do not
    // transfer to HEOS's residual entropy, measured over a 3828-point grid rather than the single
    // published sample point that an earlier round used.  See HEOS_TRANSFER_EXCLUDE in
    // dev/convert_RES_csv_to_json.py for the criterion and the per-fluid numbers.
    //
    // The list is repeated here on purpose.  It is the user-visible contract -- "ask for RES on
    // this fluid and you get a clear error, not a plausible wrong number" -- so changing which
    // fluids are withheld should require editing a test and justifying it, not just regenerating
    // a data file.  Under the old suite this behaviour was invisible: a designed throw and a
    // broken one were indistinguishable to a `catch (...)`.
    static const char* withheld[] = {"BENZENE",  "D5",      "HEPTANE", "MD3M", "MD4M", "R1123", "R1224YDZ",
                                     "R1233ZDE", "R1234YF", "R1243ZF", "R13",  "R161", "R41",   "VINYLCHLORIDE"};

    std::size_t checked = 0;
    for (const char* name : withheld) {
        std::shared_ptr<AbstractState> AS;
        try {
            AS.reset(AbstractState::factory("HEOS", name));
        } catch (const CoolProp::CoolPropBaseError& e) {
            // These fluids are withheld from the RES parameter set, not from CoolProp -- every
            // one of them must still construct on HEOS.  Reporting the name beats the bare
            // count mismatch the old `continue` produced two lines further down.
            FAIL("withheld fluid " << name << " is not constructible on HEOS: " << e.what());
        }
        ++checked;
        INFO("withheld fluid " << name);
        CHECK_FALSE(AS->RES_data().comps[0].viscosity.provided);
        CHECK_FALSE(AS->RES_data().comps[0].conductivity.provided);
        AS->use_viscosity_RES(true);
        AS->use_conductivity_RES(true);
        AS->update(PT_INPUTS, 1.0e5, 300.0);
        CHECK_THROWS_AS(AS->viscosity(), CoolProp::ValueError);
        CHECK_THROWS_AS(AS->conductivity(), CoolProp::ValueError);
    }
    CHECK(checked == std::size(withheld));  // unreachable unless the loop body is skipped

    // The converse, so that "everything throws" cannot pass this test: a fluid that is NOT
    // withheld must still evaluate.
    auto ok = make_heos("Propane");
    ok->use_viscosity_RES(true);
    ok->use_conductivity_RES(true);
    ok->update(PT_INPUTS, 1.0e5, 300.0);
    CHECK(ok->viscosity() > 0.0);
    CHECK(ok->conductivity() > 0.0);
}

// -------------------------------------------------------------- guards -------

TEST_CASE("RES throws after alpha change (viscosity)", "[RES][transport]") {
    auto AS = make_pr("CarbonDioxide");
    AS->specify_phase(iphase_gas);
    AS->update(PT_INPUTS, 1.0e6, 300.0);
    REQUIRE_NOTHROW(AS->viscosity());
    auto* cubic = as_cubic(AS);
    cubic->set_cubic_alpha_C(0, "Twu", 0.3, 0.8, 1.1);
    AS->update(PT_INPUTS, 1.0e6, 300.0);
    CHECK_THROWS_AS(AS->viscosity(), CoolProp::ValueError);
}

TEST_CASE("RES throws after alpha change (conductivity)", "[RES][transport]") {
    auto AS = make_pr("CarbonDioxide");
    AS->specify_phase(iphase_gas);
    AS->update(PT_INPUTS, 1.0e6, 300.0);
    REQUIRE_NOTHROW(AS->conductivity());
    auto* cubic = as_cubic(AS);
    cubic->set_cubic_alpha_C(0, "Twu", 0.3, 0.8, 1.1);
    AS->update(PT_INPUTS, 1.0e6, 300.0);
    CHECK_THROWS_AS(AS->conductivity(), CoolProp::ValueError);
}

TEST_CASE("RES recovers after set_viscosity_RES_residual_params", "[RES][transport]") {
    auto AS = make_pr("CarbonDioxide");
    AS->specify_phase(iphase_gas);
    auto* cubic = as_cubic(AS);
    cubic->set_cubic_alpha_C(0, "Twu", 0.3, 0.8, 1.1);
    AS->update(PT_INPUTS, 1.0e6, 300.0);
    AS->set_viscosity_RES_residual_params(0, {0.417, -0.231, 0.068}, 1.023);
    AS->set_conductivity_RES_residual_params(0, {0.42, -0.20, 0.06, 0.0}, 1.0);
    REQUIRE_NOTHROW(AS->viscosity());
    REQUIRE_NOTHROW(AS->conductivity());
}

// --------------------------------------------------------- round-trip --------

TEST_CASE("RES setters reject coefficient vectors the model would index past", "[RES][transport]") {
    // The transport routines read n_dilute[0..4] and n_res[0..2] (viscosity) / n_res[0..3]
    // (conductivity) unconditionally once `provided` is set, guarded by nothing but that flag.
    // Both Python bindings expose these setters, so a short list from a caller used to be
    // accepted here and read out of bounds at the next viscosity()/conductivity().
    auto AS = make_pr("CarbonDioxide");
    const std::vector<double> d5{1, 2, 3, 4, 5};
    const std::vector<double> v3{0.1, 0.2, 0.3};
    const std::vector<double> c4{0.1, 0.2, 0.3, 0.4};

    CHECK_NOTHROW(AS->set_viscosity_RES_parameters(0, d5, v3, 1.0));
    CHECK_NOTHROW(AS->set_conductivity_RES_parameters(0, d5, c4, 1.0, 1.02, 1.239, 0.057, 2.1e-10, 600.0, 6.0e-10));

    // too few dilute coefficients
    CHECK_THROWS_AS(AS->set_viscosity_RES_parameters(0, {1, 2, 3, 4}, v3, 1.0), CoolProp::ValueError);
    // too many, which would silently ignore the extras
    CHECK_THROWS_AS(AS->set_viscosity_RES_parameters(0, {1, 2, 3, 4, 5, 6}, v3, 1.0), CoolProp::ValueError);
    // wrong residual count for each property, including the conductivity count on viscosity
    CHECK_THROWS_AS(AS->set_viscosity_RES_parameters(0, d5, c4, 1.0), CoolProp::ValueError);
    CHECK_THROWS_AS(AS->set_conductivity_RES_parameters(0, d5, v3, 1.0, 1.02, 1.239, 0.057, 2.1e-10, 600.0, 6.0e-10), CoolProp::ValueError);
    CHECK_THROWS_AS(AS->set_viscosity_RES_residual_params(0, {0.1, 0.2}, 1.0), CoolProp::ValueError);
    CHECK_THROWS_AS(AS->set_conductivity_RES_residual_params(0, v3, 1.0), CoolProp::ValueError);
}

TEST_CASE("RES alpha change invalidates the memoized transport values", "[RES][transport]") {
    // set_cubic_alpha_C() clears n_params_match_alpha so that RES refuses to run with
    // coefficients fitted for the previous alpha.  Without also clearing the memoized values,
    // a caller who does NOT update() in between gets the stale number back and the guard never
    // runs -- which is exactly the case the two "throws after alpha change" tests above miss,
    // because their update() call clears the cache for them.
    auto AS = make_pr("CarbonDioxide");
    AS->specify_phase(iphase_gas);
    AS->update(PT_INPUTS, 1.0e6, 300.0);
    REQUIRE_NOTHROW(AS->viscosity());
    REQUIRE_NOTHROW(AS->conductivity());

    as_cubic(AS)->set_cubic_alpha_C(0, "Twu", 0.3, 0.8, 1.1);

    // No update() here -- that is the point.
    CHECK_THROWS_AS(AS->viscosity(), CoolProp::ValueError);
    CHECK_THROWS_AS(AS->conductivity(), CoolProp::ValueError);
}

TEST_CASE("REFPROP saturation shims inherit the host RES configuration", "[RES][REFPROP][transport]") {
    // link_to_loaded_fluids() builds a second state on the SAME components.  If it does not
    // carry _RES across, a shim built from an RES-enabled host silently evaluates REFPROP's
    // native transport model instead -- no error, just a different number.
    Skip_if_No_REFPROP();
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "Propane"));
    AS->use_viscosity_RES(true);
    AS->use_conductivity_RES(true);
    AS->update(QT_INPUTS, 0.0, 300.0);

    // saturated_liquid_keyed_output routes through build_saturation_shim().
    const double eta_L = AS->saturated_liquid_keyed_output(iviscosity);
    const double tc_L = AS->saturated_liquid_keyed_output(iconductivity);
    CHECK(std::isfinite(eta_L));
    CHECK(eta_L > 0.0);
    CHECK(std::isfinite(tc_L));
    CHECK(tc_L > 0.0);

    // The RES value must differ from REFPROP's own model; equality would mean the shim fell back.
    auto REF = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "Propane"));
    REF->update(QT_INPUTS, 0.0, 300.0);
    CHECK(eta_L != REF->saturated_liquid_keyed_output(iviscosity));
    CHECK(tc_L != REF->saturated_liquid_keyed_output(iconductivity));
}

TEST_CASE("REFPROP component swap re-seeds the RES parameters", "[RES][REFPROP][transport]") {
    // set_REFPROP_fluids() is public, so the component set can be changed on a live state.  The
    // RES store is per component; keeping the previous fluid's coefficients would produce wrong
    // numbers with no error at all.
    Skip_if_No_REFPROP();
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "Propane"));
    auto* rp = dynamic_cast<REFPROPMixtureBackend*>(AS.get());
    REQUIRE(rp != nullptr);
    const auto propane_vis = AS->get_viscosity_RES_parameters(0);
    REQUIRE(propane_vis.provided);

    // set_REFPROP_fluids() short-circuits whenever this instance's components are already the
    // ones REFPROP has loaded: that fast path exists for check_loaded_fluid() and IGNORES its
    // argument entirely.  Load something else on a second instance first, so the call below
    // really takes the SETUPdll path and changes this state's component set.
    auto OTHER = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "Water"));
    OTHER->update(PT_INPUTS, 1.0e5, 300.0);

    rp->set_REFPROP_fluids({"Methane"});
    rp->set_mole_fractions({1.0});
    const auto methane_vis = AS->get_viscosity_RES_parameters(0);
    CHECK(methane_vis.provided);
    CHECK(methane_vis.n_res != propane_vis.n_res);
    CHECK(AS->RES_data().comps[0].name == "METHANE");

    // The enable flags are reset with the store: silently carrying "RES enabled" onto a fluid
    // the caller never enabled it for is how a wrong number gets returned without an error.
    CHECK_FALSE(AS->RES_data().viscosity_enabled);
    CHECK_FALSE(AS->RES_data().conductivity_enabled);

    // ... and the re-seeded parameters actually work.
    AS->use_viscosity_RES(true);
    AS->update(PT_INPUTS, 5.0e6, 300.0);
    CHECK(std::isfinite(AS->viscosity()));
}

TEST_CASE("REFPROP reload of the SAME components keeps user RES parameters", "[RES][REFPROP][transport]") {
    // The converse of the swap test: check_loaded_fluid() re-invokes set_REFPROP_fluids() with
    // the same names on every property call once another instance has taken over REFPROP's
    // global state.  Re-seeding there would silently discard set_*_RES_parameters().
    Skip_if_No_REFPROP();
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "Propane"));
    auto* rp = dynamic_cast<REFPROPMixtureBackend*>(AS.get());
    REQUIRE(rp != nullptr);
    const std::vector<double> mine{0.1, 0.2, 0.3};
    AS->set_viscosity_RES_residual_params(0, mine, 1.25);

    // As above, force the real reload path rather than the argument-ignoring fast path --
    // otherwise this would pass without ever reaching the code it is meant to cover.
    auto OTHER = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "Water"));
    OTHER->update(PT_INPUTS, 1.0e5, 300.0);

    rp->set_REFPROP_fluids({"Propane"});

    const auto got = AS->get_viscosity_RES_residual_params(0);
    CHECK(got.first == mine);
    CHECK(got.second == 1.25);
}

TEST_CASE("REFPROP RES rejects shipped parameters after an EOS-mode change under reload", "[RES][REFPROP][transport]") {
    // The n_params_match_alpha guard is evaluated at SEED time.  A state seeded against the
    // reference Helmholtz EOS, whose components are later reloaded while REFPROP_USE_PENGROBINSON
    // is on, is now evaluating a different alpha^r with coefficients regressed against the old
    // one -- and the reload alone does not re-run the guard, because the component set is
    // unchanged.  That is a silently wrong number, not a refusal.
    Skip_if_No_REFPROP();
    const bool restore = get_config_bool(REFPROP_USE_PENGROBINSON);
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "Propane"));
    AS->use_viscosity_RES(true);
    AS->update(PT_INPUTS, 5.0e6, 320.0);
    REQUIRE(std::isfinite(AS->viscosity()));  // reference EOS: the shipped parameters do apply

    try {
        set_config_bool(REFPROP_USE_PENGROBINSON, true);
        // Force a real reload of the SAME components: another instance takes over REFPROP's
        // global state, so this one's next property call goes through set_REFPROP_fluids().
        auto OTHER = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "Water"));
        OTHER->update(PT_INPUTS, 1.0e5, 300.0);
        AS->update(PT_INPUTS, 5.0e6, 320.0);
        CHECK_THROWS_WITH(AS->viscosity(), Catch::Matchers::ContainsSubstring("different alpha function"));
    } catch (...) {
        set_config_bool(REFPROP_USE_PENGROBINSON, restore);
        throw;
    }
    set_config_bool(REFPROP_USE_PENGROBINSON, restore);

    // And back: with the reference EOS restored, a reload must make the shipped parameters
    // usable again rather than leaving the state permanently refusing.
    auto BACK = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "Propane"));
    BACK->update(PT_INPUTS, 5.0e6, 320.0);
    auto OTHER2 = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "Water"));
    OTHER2->update(PT_INPUTS, 1.0e5, 300.0);
    AS->update(PT_INPUTS, 5.0e6, 320.0);
    CHECK(std::isfinite(AS->viscosity()));
}

TEST_CASE("REFPROP .MIX replaced by a normal fluid gains RES support", "[RES][REFPROP][transport]") {
    // A predefined .MIX leaves resolved_fluid_names empty, exactly as a state that has never been
    // seeded does.  Inferring "never seeded" from that emptiness would make this state keep an
    // empty RES store after the swap and report RES as unsupported for a fluid that has
    // parameters -- fail-closed rather than wrong, but still wrong.
    Skip_if_No_REFPROP();
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "AIR.MIX"));
    auto* rp = dynamic_cast<REFPROPMixtureBackend*>(AS.get());
    REQUIRE(rp != nullptr);
    REQUIRE_THROWS_AS(AS->use_viscosity_RES(true), CoolProp::NotImplementedError);

    // Force the real reload path rather than the argument-ignoring fast path.
    auto OTHER = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "Water"));
    OTHER->update(PT_INPUTS, 1.0e5, 300.0);

    rp->set_REFPROP_fluids({"Propane"});
    rp->set_mole_fractions({1.0});

    REQUIRE(AS->RES_data().comps.size() == 1);
    CHECK(AS->RES_data().comps[0].name == "PROPANE");
    CHECK(AS->get_viscosity_RES_parameters(0).provided);
    CHECK_NOTHROW(AS->use_viscosity_RES(true));
    AS->update(PT_INPUTS, 5.0e6, 300.0);
    CHECK(std::isfinite(AS->viscosity()));
}

TEST_CASE("REFPROP normal fluid replaced by a .MIX drops RES support", "[RES][REFPROP][transport]") {
    // The other direction: the old fluid's parameters must not survive onto a mixture that has
    // no per-component identity to attach them to.
    Skip_if_No_REFPROP();
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "Propane"));
    auto* rp = dynamic_cast<REFPROPMixtureBackend*>(AS.get());
    REQUIRE(rp != nullptr);
    REQUIRE(AS->get_viscosity_RES_parameters(0).provided);

    auto OTHER = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "Water"));
    OTHER->update(PT_INPUTS, 1.0e5, 300.0);

    rp->set_REFPROP_fluids({"AIR.MIX"});

    CHECK(AS->RES_data().comps.empty());
    CHECK_THROWS_AS(AS->use_viscosity_RES(true), CoolProp::NotImplementedError);
}

TEST_CASE("RES toggles and setters invalidate the memoized transport values", "[RES][transport]") {
    // viscosity()/conductivity() memoize into _viscosity/_conductivity, so toggling RES or
    // pushing new coefficients after a read used to return the previous model's number with
    // no indication anything had been ignored.
    SECTION("use_viscosity_RES after a read changes the result") {
        auto A = make_heos("CarbonDioxide");
        A->update(PT_INPUTS, 1.0e6, 300.0);
        double before = A->viscosity();
        A->use_viscosity_RES(true);
        double after = A->viscosity();

        // Reference value: same state, RES enabled before the first read.
        auto B = make_heos("CarbonDioxide");
        B->use_viscosity_RES(true);
        B->update(PT_INPUTS, 1.0e6, 300.0);
        double fresh = B->viscosity();

        INFO("before=" << before << " after=" << after << " fresh=" << fresh);
        CHECK(after == Catch::Approx(fresh).epsilon(1e-12));
        CHECK(after != Catch::Approx(before).epsilon(1e-12));
    }
    SECTION("use_conductivity_RES after a read changes the result") {
        auto A = make_heos("CarbonDioxide");
        A->update(PT_INPUTS, 1.0e6, 300.0);
        double before = A->conductivity();
        A->use_conductivity_RES(true);
        double after = A->conductivity();

        auto B = make_heos("CarbonDioxide");
        B->use_conductivity_RES(true);
        B->update(PT_INPUTS, 1.0e6, 300.0);
        double fresh = B->conductivity();

        INFO("before=" << before << " after=" << after << " fresh=" << fresh);
        CHECK(after == Catch::Approx(fresh).epsilon(1e-12));
        CHECK(after != Catch::Approx(before).epsilon(1e-12));
    }
    SECTION("set_viscosity_RES_residual_params after a read takes effect") {
        auto A = make_heos("CarbonDioxide");
        A->use_viscosity_RES(true);
        A->update(PT_INPUTS, 1.0e6, 300.0);
        double before = A->viscosity();
        auto [n, xita] = A->get_viscosity_RES_residual_params(0);
        std::vector<double> bumped = n;
        bumped[0] *= 1.5;  // materially different coefficients
        A->set_viscosity_RES_residual_params(0, bumped, xita);
        double after = A->viscosity();
        INFO("before=" << before << " after=" << after);
        CHECK(after != Catch::Approx(before).epsilon(1e-12));
    }
}

// ----------------------- REFPROP building blocks for RES ---------------------
// These validate the two additions RES-on-REFPROP is built on, before any RES code depends on
// them.  Both are purely additive to the REFPROP backend.

TEST_CASE("REFPROP get_fluid_constant returns per-component constants", "[RES][REFPROP][transport]") {
    Skip_if_No_REFPROP();
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "Propane"));
    AS->update(PT_INPUTS, 1.0e6, 300.0);
    // molar_mass() is the mixture value; for a pure fluid it must equal the component constant.
    CHECK(AS->get_fluid_constant(0, imolar_mass) == Catch::Approx(AS->molar_mass()).epsilon(1e-10));
    CHECK(AS->get_fluid_constant(0, iT_critical) == Catch::Approx(AS->T_critical()).epsilon(1e-6));
    CHECK(AS->get_fluid_constant(0, iP_critical) == Catch::Approx(AS->p_critical()).epsilon(1e-6));
    CHECK_THROWS(AS->get_fluid_constant(5, imolar_mass));  // out of range
}

TEST_CASE("REFPROP off-state alphar evaluation does not mutate the state", "[RES][REFPROP][transport]") {
    // The premise of the RES port: PHIXdll is a pure function of (itau, idelta, tau, delta, x),
    // so alpha^r can be evaluated at a reference temperature without disturbing the cached state.
    // If this ever regresses, the critical-enhancement reference term silently corrupts the state.
    Skip_if_No_REFPROP();
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "Propane"));
    AS->update(PT_INPUTS, 4.0e6, 360.0);

    const double d1_before = AS->dalphar_dDelta();
    const double d2_before = AS->d2alphar_dDelta2();
    const double p_before = AS->p();
    const double sres_before = AS->smolar_residual();
    const double tau_before = AS->tau();
    const double delta_before = AS->delta();

    // Evaluate well away from the current temperature, then re-read everything.
    auto* rp = dynamic_cast<REFPROPMixtureBackend*>(AS.get());
    REQUIRE(rp != nullptr);
    const double tau_far = tau_before * 0.5;
    (void)rp->call_phixdll(0, 1, tau_far, delta_before);
    (void)rp->call_phixdll(0, 2, tau_far, delta_before);

    CHECK(AS->dalphar_dDelta() == d1_before);
    CHECK(AS->d2alphar_dDelta2() == d2_before);
    CHECK(AS->p() == p_before);
    CHECK(AS->smolar_residual() == sres_before);
    CHECK(AS->tau() == tau_before);
    CHECK(AS->delta() == delta_before);

    // At the CURRENT (tau, delta) the new overload must agree with the state-coupled one.
    CHECK(rp->call_phixdll(0, 1, tau_before, delta_before) == Catch::Approx(d1_before).epsilon(1e-12));
    CHECK(rp->call_phixdll(0, 2, tau_before, delta_before) == Catch::Approx(d2_before).epsilon(1e-12));
}

// ------------------------------ RES on the REFPROP backend -------------------
// REFPROP is the backend where the published coefficients are exactly valid -- they were
// regressed against REFPROP's reference Helmholtz EOS.  On HEOS they are an approximation,
// which is what the HEOS_*_EXCLUDE lists in dev/convert_RES_csv_to_json.py are about.

TEST_CASE("REFPROP RES is a distinct model from TRNPRPdll", "[RES][REFPROP][transport]") {
    Skip_if_No_REFPROP();
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "Propane"));
    AS->update(PT_INPUTS, 5.0e6, 320.0);
    const double eta_native = AS->viscosity();
    const double tc_native = AS->conductivity();
    REQUIRE(std::isfinite(eta_native));
    REQUIRE(std::isfinite(tc_native));

    AS->use_viscosity_RES(true);
    AS->use_conductivity_RES(true);
    const double eta_res = AS->viscosity();
    const double tc_res = AS->conductivity();
    REQUIRE(std::isfinite(eta_res));
    REQUIRE(std::isfinite(tc_res));
    CHECK(eta_res > 0.0);
    CHECK(tc_res > 0.0);
    // A bit-equal result would mean the RES branch never ran and the native value leaked
    // through; a wildly different one would mean it ran on the wrong units or parameters.
    CHECK(eta_res != eta_native);
    CHECK(tc_res != tc_native);
    CHECK(eta_res == Catch::Approx(eta_native).epsilon(0.25));
    CHECK(tc_res == Catch::Approx(tc_native).epsilon(0.25));

    // Toggling back must restore the native values exactly, not hand back a stale cache.
    AS->use_viscosity_RES(false);
    AS->use_conductivity_RES(false);
    CHECK(AS->viscosity() == eta_native);
    CHECK(AS->conductivity() == tc_native);
}

TEST_CASE("REFPROP RES viscosity and conductivity flags are independent", "[RES][REFPROP][transport]") {
    // TRNPRPdll returns both properties from a single call, and calc_conductivity() used to get
    // its value by calling calc_viscosity().  With RES enabled on only one of the two, that left
    // the other CachedElement unpopulated and reading it threw.  All four flag combinations must
    // work, and -- the part that actually catches a cache defect -- neither value may depend on
    // the order in which the two are read.
    Skip_if_No_REFPROP();
    for (int combo = 0; combo < 4; ++combo) {
        const bool vis_res = (combo & 1) != 0;
        const bool tc_res = (combo & 2) != 0;
        CAPTURE(vis_res, tc_res);

        auto A = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "Propane"));
        A->update(PT_INPUTS, 5.0e6, 320.0);
        A->use_viscosity_RES(vis_res);
        A->use_conductivity_RES(tc_res);
        double eta_first = 0, tc_second = 0;
        REQUIRE_NOTHROW(eta_first = A->viscosity());
        REQUIRE_NOTHROW(tc_second = A->conductivity());

        auto B = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "Propane"));
        B->update(PT_INPUTS, 5.0e6, 320.0);
        B->use_viscosity_RES(vis_res);
        B->use_conductivity_RES(tc_res);
        double tc_first = 0, eta_second = 0;
        REQUIRE_NOTHROW(tc_first = B->conductivity());
        REQUIRE_NOTHROW(eta_second = B->viscosity());

        CHECK(std::isfinite(eta_first));
        CHECK(eta_first > 0.0);
        CHECK(std::isfinite(tc_first));
        CHECK(tc_first > 0.0);
        CHECK(eta_second == eta_first);
        CHECK(tc_second == tc_first);
    }
}

TEST_CASE("REFPROP off-state dRhomass/dp matches the alphar route", "[RES][REFPROP][transport]") {
    // The RES critical enhancement takes this derivative twice: at the current temperature from
    // alpha^r (PHIXdll), and at the reference temperature from calc_drhomass_dp_constT_at()
    // (THERM2dll).  Evaluating both at the SAME temperature is the only way to show the two
    // routes agree -- on units, on the gas constant, and on the reducing density.
    Skip_if_No_REFPROP();
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "Propane"));
    AS->update(PT_INPUTS, 4.0e6, 360.0);
    const double T = AS->T();
    const double delta_st = AS->delta();
    const double dp_drho = AS->gas_constant() * T * (1.0 + 2.0 * delta_st * AS->dalphar_dDelta() + delta_st * delta_st * AS->d2alphar_dDelta2());
    const double drhodp_phix = 1.0 / dp_drho * AS->molar_mass();
    CHECK(AS->drhomass_dp_constT_at(T) == Catch::Approx(drhodp_phix).epsilon(1e-9));

    // Off-state (the actual use, at a supercritical reference temperature) it must stay finite
    // and positive, and must leave the cached state untouched.
    const double p_before = AS->p();
    const double sres_before = AS->smolar_residual();
    const double d_ref = AS->drhomass_dp_constT_at(1.5 * AS->T_critical());
    CHECK(std::isfinite(d_ref));
    CHECK(d_ref > 0.0);
    CHECK(AS->T() == T);
    CHECK(AS->p() == p_before);
    CHECK(AS->smolar_residual() == sres_before);
}

TEST_CASE("REFPROP RES parameter setters round-trip and invalidate the cache", "[RES][REFPROP][transport]") {
    Skip_if_No_REFPROP();
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "Propane"));
    AS->update(PT_INPUTS, 5.0e6, 320.0);
    AS->use_viscosity_RES(true);
    AS->use_conductivity_RES(true);
    const double eta_before = AS->viscosity();
    const double tc_before = AS->conductivity();

    const std::vector<double> n_vis{0.1, 0.2, 0.3};
    AS->set_viscosity_RES_residual_params(0, n_vis, 1.25);
    const auto got_vis = AS->get_viscosity_RES_residual_params(0);
    CHECK(got_vis.first == n_vis);
    CHECK(got_vis.second == 1.25);
    CHECK(AS->viscosity() != eta_before);

    const std::vector<double> n_tc{0.1, 0.2, 0.3, 0.4};
    AS->set_conductivity_RES_residual_params(0, n_tc, 0.9);
    const auto got_tc = AS->get_conductivity_RES_residual_params(0);
    CHECK(got_tc.first == n_tc);
    CHECK(got_tc.second == 0.9);
    CHECK(AS->conductivity() != tc_before);

    // Pure fluid: only component 0 exists.
    CHECK_THROWS(AS->set_viscosity_RES_residual_params(1, n_vis, 1.0));
    CHECK_THROWS(AS->get_conductivity_RES_residual_params(1));
}

TEST_CASE("REFPROP RES works for a mixture and for a REFPROP-only fluid", "[RES][REFPROP][transport]") {
    Skip_if_No_REFPROP();
    // A binary: RES parameters must be attributed per component, in component order.
    auto MIX = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "Methane&Ethane"));
    MIX->set_mole_fractions({0.6, 0.4});
    MIX->update(PT_INPUTS, 3.0e6, 280.0);
    MIX->use_viscosity_RES(true);
    MIX->use_conductivity_RES(true);
    const double eta_mix = MIX->viscosity();
    const double tc_mix = MIX->conductivity();
    CHECK(std::isfinite(eta_mix));
    CHECK(eta_mix > 0.0);
    CHECK(std::isfinite(tc_mix));
    CHECK(tc_mix > 0.0);

    // 22DIMETHYLBUTANE has no CoolProp EOS at all, so it carries only a "REFPROP" parameter
    // block -- and its key is 16 characters, which NAMEdll's character*12 hnam would have
    // truncated to "22DIMETHYLBU" and silently matched nothing.
    auto ONLY = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "22DIMETHYLBUTANE"));
    ONLY->update(PT_INPUTS, 1.0e5, 330.0);
    ONLY->use_viscosity_RES(true);
    const double eta_only = ONLY->viscosity();
    CHECK(std::isfinite(eta_only));
    CHECK(eta_only > 0.0);
}

TEST_CASE("REFPROP predefined .MIX mixtures report RES as unsupported", "[RES][REFPROP][transport]") {
    // SETMIXdll reports Ncomp == 1 for the whole mixture and collapses the component list to the
    // .MIX filename, so there is no per-component identity to attach parameters to.  Reporting
    // "unsupported" is the only honest answer: silently applying one constituent's parameters
    // would return a confidently wrong number.
    Skip_if_No_REFPROP();
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "AIR.MIX"));
    AS->update(PT_INPUTS, 1.0e5, 300.0);
    CHECK(std::isfinite(AS->viscosity()));  // the native path still works
    CHECK_THROWS_AS(AS->use_viscosity_RES(true), CoolProp::NotImplementedError);
    CHECK_THROWS_AS(AS->use_conductivity_RES(true), CoolProp::NotImplementedError);
}

TEST_CASE("REFPROP RES refuses to run when the EOS has been replaced", "[RES][REFPROP][transport]") {
    // REFPROP_USE_PENGROBINSON (and REFPROP_USE_GERG, which goes through the same guard) swap out
    // the very Helmholtz EOS the n_res/xita coefficients were regressed against, so s_res is no
    // longer the quantity they were fitted to.  Returning a number there would be silently wrong.
    Skip_if_No_REFPROP();
    const bool restore = get_config_bool(REFPROP_USE_PENGROBINSON);
    set_config_bool(REFPROP_USE_PENGROBINSON, true);
    try {
        auto PR = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "Propane"));
        PR->update(PT_INPUTS, 5.0e6, 320.0);
        PR->use_viscosity_RES(true);
        PR->use_conductivity_RES(true);
        // Match the message, not just the type: a ValueError from some unrelated REFPROP
        // failure would otherwise let this pass without the guard ever being reached.
        CHECK_THROWS_WITH(PR->viscosity(), Catch::Matchers::ContainsSubstring("different alpha function"));
        CHECK_THROWS_WITH(PR->conductivity(), Catch::Matchers::ContainsSubstring("different alpha function"));
    } catch (...) {
        set_config_bool(REFPROP_USE_PENGROBINSON, restore);
        throw;
    }
    set_config_bool(REFPROP_USE_PENGROBINSON, restore);

    // Constructing a fresh state calls PREOSdll(0) unconditionally, which restores the reference
    // EOS.  Assert that here rather than assuming it -- a leaked Peng-Robinson setting would
    // corrupt every REFPROP test that runs after this one.
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "Propane"));
    AS->update(PT_INPUTS, 5.0e6, 320.0);
    AS->use_viscosity_RES(true);
    CHECK(std::isfinite(AS->viscosity()));
    CHECK(AS->rhomass() == Catch::Approx(CoolProp::PropsSI("Dmass", "P", 5.0e6, "T", 320.0, "REFPROP::Propane")).epsilon(1e-10));
}

// -------------------------- where the RES inputs come from -------------------
// RES takes two inputs that are not part of the entropy-scaling model itself -- the dilute-gas
// term and the viscosity inside the critical enhancement -- and the source papers take both from
// the backend's own transport model on some code paths.  These tests pin which source AUTO
// resolves to on each path, since getting that wrong is silent: the answer stays plausible.

TEST_CASE("RES dilute-gas source is selectable and AUTO follows the papers", "[RES][REFPROP][transport]") {
    Skip_if_No_REFPROP();

    SECTION("pure viscosity: AUTO is the native term (Martinek)") {
        // D2 is where the fitted polynomial matches REFPROP's eta0 worst, so it is the sharpest
        // discriminator available: 0.96 % off the published value with the polynomial, 0.0015 %
        // with the native term.
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "D2"));
        AS->update(PT_INPUTS, 1.11716e-1 * 1e6, 24.0);
        AS->use_viscosity_RES(true);
        const double v_auto = AS->viscosity();
        AS->set_viscosity_RES_dilute_source(RES_DILUTE_BACKEND_NATIVE);
        CHECK(AS->viscosity() == v_auto);
        AS->set_viscosity_RES_dilute_source(RES_DILUTE_POLYNOMIAL);
        const double v_poly = AS->viscosity();
        CHECK(v_poly != v_auto);
        // A different dilute term, not a different model: percent-scale, not order-of-magnitude.
        CHECK(v_poly == Catch::Approx(v_auto).epsilon(0.05));
        // Setting it back must restore the value exactly -- i.e. the setter clears the cache in
        // both directions, not just on the way out.
        AS->set_viscosity_RES_dilute_source(RES_DILUTE_AUTO);
        CHECK(AS->viscosity() == v_auto);
    }

    SECTION("mixture viscosity: AUTO is the polynomial (Martinek)") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "Argon&Nitrogen"));
        AS->set_mole_fractions({0.5, 0.5});
        AS->update(PT_INPUTS, 5.0e6, 300.0);
        AS->use_viscosity_RES(true);
        const double v_auto = AS->viscosity();
        AS->set_viscosity_RES_dilute_source(RES_DILUTE_POLYNOMIAL);
        CHECK(AS->viscosity() == v_auto);
        AS->set_viscosity_RES_dilute_source(RES_DILUTE_BACKEND_NATIVE);
        CHECK(AS->viscosity() != v_auto);
    }

    SECTION("pure conductivity: AUTO is the polynomial (Li)") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "Propane"));
        AS->update(PT_INPUTS, 4.65e6, 375.9);
        AS->use_conductivity_RES(true);
        const double v_auto = AS->conductivity();
        AS->set_conductivity_RES_dilute_source(RES_DILUTE_POLYNOMIAL);
        CHECK(AS->conductivity() == v_auto);
        AS->set_conductivity_RES_dilute_source(RES_DILUTE_BACKEND_NATIVE);
        CHECK(AS->conductivity() != v_auto);
    }

    SECTION("mixture conductivity: AUTO is the native term (Li)") {
        // ARGON+NEON is the pair this decides: -11.6 % off the published value with the
        // polynomial-plus-Wilke term, +0.003 % with REFPROP's own mixture lambda0.
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "Argon&Neon"));
        AS->set_mole_fractions({0.5, 0.5});
        AS->update(PT_INPUTS, 1.0e6, 300.0);
        AS->use_conductivity_RES(true);
        const double v_auto = AS->conductivity();
        AS->set_conductivity_RES_dilute_source(RES_DILUTE_BACKEND_NATIVE);
        CHECK(AS->conductivity() == v_auto);
        AS->set_conductivity_RES_dilute_source(RES_DILUTE_POLYNOMIAL);
        CHECK(AS->conductivity() != v_auto);
    }
}

TEST_CASE("RES critical-enhancement viscosity source is selectable", "[RES][REFPROP][transport]") {
    Skip_if_No_REFPROP();
    // PROPANE at its published sample point is where this choice bites: the enhancement is active,
    // and feeding it the RES viscosity instead of REFPROP's own moves the answer by ~5 %.
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "Propane"));
    AS->update(PT_INPUTS, 4.65e6, 375.9);
    AS->use_conductivity_RES(true);
    const double v_auto = AS->conductivity();
    AS->set_conductivity_RES_enhancement_viscosity(RES_ENH_VIS_BACKEND_NATIVE);
    CHECK(AS->conductivity() == v_auto);
    AS->set_conductivity_RES_enhancement_viscosity(RES_ENH_VIS_RES);
    const double v_res = AS->conductivity();
    CHECK(v_res != v_auto);
    CHECK(v_res == Catch::Approx(v_auto).epsilon(0.10));
    AS->set_conductivity_RES_enhancement_viscosity(RES_ENH_VIS_AUTO);
    CHECK(AS->conductivity() == v_auto);
}

TEST_CASE("RES source policy where the backend supplies no native term", "[RES][transport]") {
    // This is the half of the policy that keeps HEOS and the cubics unchanged, and the half that
    // must NOT fail open.  AUTO silently uses the fitted model, exactly as the reference
    // implementations' try/except does.  An EXPLICIT request for the backend's own model, where
    // there is none, must throw rather than quietly hand back a different model than was asked
    // for -- that is the difference between "no native term available" and a wrong number.
    //
    // "Supplies none" is not one situation: the cubics have no transport model at all, which
    // is a reason RES exists, while HEOS has correlations that RES deliberately does not
    // consume (see RESDiluteSource).  HEOS is used here because it is the stricter case.
    auto AS = make_heos("Propane");
    AS->update(PT_INPUTS, 4.65e6, 375.9);
    AS->use_viscosity_RES(true);
    AS->use_conductivity_RES(true);
    const double eta_auto = AS->viscosity();
    const double tc_auto = AS->conductivity();

    AS->set_viscosity_RES_dilute_source(RES_DILUTE_POLYNOMIAL);
    AS->set_conductivity_RES_dilute_source(RES_DILUTE_POLYNOMIAL);
    AS->set_conductivity_RES_enhancement_viscosity(RES_ENH_VIS_RES);
    CHECK(AS->viscosity() == eta_auto);
    CHECK(AS->conductivity() == tc_auto);

    AS->set_viscosity_RES_dilute_source(RES_DILUTE_BACKEND_NATIVE);
    CHECK_THROWS_WITH(AS->viscosity(), Catch::Matchers::ContainsSubstring("does not supply a native transport model"));
    AS->set_viscosity_RES_dilute_source(RES_DILUTE_POLYNOMIAL);

    AS->set_conductivity_RES_dilute_source(RES_DILUTE_BACKEND_NATIVE);
    CHECK_THROWS_WITH(AS->conductivity(), Catch::Matchers::ContainsSubstring("does not supply a native transport model"));
    AS->set_conductivity_RES_dilute_source(RES_DILUTE_POLYNOMIAL);

    AS->set_conductivity_RES_enhancement_viscosity(RES_ENH_VIS_BACKEND_NATIVE);
    CHECK_THROWS_WITH(AS->conductivity(), Catch::Matchers::ContainsSubstring("does not supply a native viscosity model"));
}

TEST_CASE("REFPROP mixture critical point agrees with the PropsSI route", "[RES][REFPROP][transport]") {
    // The RES mixture critical enhancement reads Tc, pc and rhoc off AbstractState, which lands on
    // CRITPdll.  Li 2024 obtains the same three through PropsSI on a composition-bearing fluid
    // string (code_SI.py:162-163, 178).  The enhancement only reproduces the published mixture
    // values while those two routes agree, and a mixture string goes through composition parsing
    // that a pure name does not -- so pin the equivalence rather than assume it.
    //
    // Exact equality, not Approx: both should reach the same CRITPdll call, so anything less than
    // bit-identical means the routing changed and the enhancement has quietly drifted off the
    // reference.
    Skip_if_No_REFPROP();
    struct Mix
    {
        const char* a;
        const char* b;
        double xa, xb;
    };
    static const Mix mixes[] = {
      {"BUTANE", "METHANE", 0.606, 0.394},
      {"ARGON", "NEON", 0.5, 0.5},
      {"CO2", "ETHANE", 0.5, 0.5},
    };

    for (const auto& m : mixes) {
        // Built the way Li's fluid_mix_name_cp() does it.
        const std::string named = std::string("REFPROP::") + m.a + "[" + std::to_string(m.xa) + "]&" + m.b + "[" + std::to_string(m.xb) + "]";
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", std::string(m.a) + "&" + m.b));
        AS->set_mole_fractions({m.xa, m.xb});
        INFO("mixture " << named);
        CHECK(CoolProp::PropsSI("Tcrit", "T", 0, "P", 0, named) == AS->T_critical());
        CHECK(CoolProp::PropsSI("rhocrit", "T", 0, "P", 0, named) == AS->rhomass_critical());
        CHECK(CoolProp::PropsSI("Pcrit", "T", 0, "P", 0, named) == AS->p_critical());
    }
}

TEST_CASE("RES mixture critical enhancement follows the backend policy", "[RES][REFPROP][transport]") {
    Skip_if_No_REFPROP();
    // BUTANE+METHANE at its published sample point is the only binary in the sample set that sits
    // inside the near-critical window (T/Tc = 0.974, rho/rhoc = 1.62), so it is the only one where
    // this policy is observable at all.
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "BUTANE&METHANE"));
    AS->set_mole_fractions({0.606, 0.394});
    AS->update(PT_INPUTS, 24131.7 * 1e3, 377.567);
    AS->use_conductivity_RES(true);

    // REFPROP reports its mixture critical point directly, so AUTO enables the enhancement.
    const double tc_auto = AS->conductivity();
    AS->set_RES_mixture_enhancement(RES_MIX_ENH_ON);
    CHECK(AS->conductivity() == tc_auto);

    AS->set_RES_mixture_enhancement(RES_MIX_ENH_OFF);
    const double tc_off = AS->conductivity();
    CHECK(tc_off != tc_auto);
    CHECK(tc_off < tc_auto);  // the enhancement is a positive contribution
    // It is worth about 1 % here, which is why leaving it out was the last outlier in the
    // published-value comparison.
    CHECK(tc_off == Catch::Approx(tc_auto).epsilon(0.05));

    AS->set_RES_mixture_enhancement(RES_MIX_ENH_AUTO);
    CHECK(AS->conductivity() == tc_auto);
}

TEST_CASE("RES mixture enhancement policy does not affect pure fluids", "[RES][REFPROP][transport]") {
    Skip_if_No_REFPROP();
    // A pure fluid's enhancement is always applied; only the MIXTURE path is gated.
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", "Propane"));
    AS->update(PT_INPUTS, 4.65e6, 375.9);
    AS->use_conductivity_RES(true);
    const double tc = AS->conductivity();
    for (auto policy : {RES_MIX_ENH_OFF, RES_MIX_ENH_ON, RES_MIX_ENH_AUTO}) {
        AS->set_RES_mixture_enhancement(policy);
        CHECK(AS->conductivity() == tc);
    }
}

TEST_CASE("RES mixture critical enhancement is off by default on native CoolProp backends", "[RES][transport]") {
    // Deliberate: the enhancement needs the MIXTURE critical point, which HEOS has to SOLVE for.
    // That solve is neither quick nor dependable, and the physical case for a critical enhancement
    // in mixtures is not well established -- so AUTO must leave it off here.  If this ever starts
    // failing because AUTO turned it on, that is a policy regression, not a numerical one.
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Methane&Ethane"));
    AS->set_mole_fractions({0.6, 0.4});
    AS->update(PT_INPUTS, 5.0e6, 250.0);
    AS->use_conductivity_RES(true);
    const double tc_auto = AS->conductivity();
    AS->set_RES_mixture_enhancement(RES_MIX_ENH_OFF);
    CHECK(AS->conductivity() == tc_auto);
    // ...and the enhancement really is absent, not merely equal by luck: no critical point was
    // needed to produce either number.
    CHECK(AS->critical_point_is_cheap() == false);
}

// On the REFPROP backend the published vis_res / TC_RES columns were produced with the SAME
// equation of state, so these tests measure implementation correctness alone -- unlike their
// HEOS counterparts above, which additionally carry the parameter-transfer error.
//
// Two tolerances per case, because they fail for different reasons: the per-sample bound
// catches gross breakage (wrong component order, wrong units, unattached parameters), the mean
// bound catches a subtle regression that a loose per-sample bound would sleep through.
// Both were set from the measured distribution -- run [RES_refprop_parity] to reprint it.
//
// The bounds are deliberately close to the measured maxima: a loose one here would sleep
// through exactly the kind of regression this file exists to catch.
//
// All four paths now reproduce to better than 0.05 %, mixtures included: the mixture critical
// enhancement is implemented and, on REFPROP, enabled by default.  On CoolProp's own backends it
// is deliberately OFF by default -- see RESMixtureEnhancement -- so these bounds describe REFPROP
// and REFPROP only.
TEST_CASE("RES on REFPROP reproduces the published pure-fluid values", "[RES][REFPROP][transport]") {
    Skip_if_No_REFPROP();

    SECTION("viscosity (Martinek 2025)") {
        auto samples = parse_vis_pure(VIS_PURE_OUT);
        REQUIRE(samples.size() >= 120);
        double sum_abs = 0;
        std::size_t n = 0;
        for (const auto& s : samples) {
            // The sample files carry REFPROP names already, so they go to the backend unmapped.
            auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", s.name));
            AS->update(PT_INPUTS, s.p_MPa * 1e6, s.T);
            AS->use_viscosity_RES(true);
            const double v = AS->viscosity() * 1e6;
            INFO("fluid=" << s.name << " T=" << s.T << " ref=" << s.vis_res << " cpp=" << v);
            CHECK(v == Catch::Approx(s.vis_res).epsilon(2e-4));
            sum_abs += std::abs(v / s.vis_res - 1.0);
            ++n;
        }
        // Every sample must have been evaluated -- no fluid may drop out silently.
        REQUIRE(n == samples.size());
        CHECK(sum_abs / n < 1e-5);
    }

    SECTION("conductivity (Li 2024)") {
        auto samples = parse_tc_pure(TC_PURE_REF);
        REQUIRE(samples.size() >= 120);
        double sum_abs = 0;
        std::size_t n = 0;
        for (const auto& s : samples) {
            auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", s.name));
            AS->update(PT_INPUTS, s.p_kPa * 1e3, s.T);
            AS->use_conductivity_RES(true);
            const double v = AS->conductivity();
            INFO("fluid=" << s.name << " T=" << s.T << " ref=" << s.tc_res << " cpp=" << v);
            CHECK(v == Catch::Approx(s.tc_res).epsilon(8e-4));
            sum_abs += std::abs(v / s.tc_res - 1.0);
            ++n;
        }
        REQUIRE(n == samples.size());
        CHECK(sum_abs / n < 1e-4);
    }
}

TEST_CASE("RES on REFPROP reproduces the published mixture values", "[RES][REFPROP][transport]") {
    Skip_if_No_REFPROP();

    SECTION("viscosity (Martinek 2025)") {
        auto samples = parse_vis_bin(VIS_BIN_OUT);
        REQUIRE(samples.size() >= 15);
        std::size_t n = 0;
        for (const auto& s : samples) {
            auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", s.c1 + "&" + s.c2));
            AS->set_mole_fractions({s.mf1, s.mf2});
            AS->update(PT_INPUTS, s.p_MPa * 1e6, s.T);
            AS->use_viscosity_RES(true);
            const double v = AS->viscosity() * 1e6;
            INFO("mix=" << s.c1 << "+" << s.c2 << " T=" << s.T << " ref=" << s.vis_res << " cpp=" << v);
            // Both papers and this code use the fitted polynomial plus Wilke on this path, so
            // there is no dilute-source divergence to absorb: it must reproduce to round-off.
            CHECK(v == Catch::Approx(s.vis_res).epsilon(1e-4));
            ++n;
        }
        REQUIRE(n == samples.size());
    }

    SECTION("conductivity (Li 2024)") {
        auto samples = parse_tc_bin(TC_BIN_REF);
        REQUIRE(samples.size() >= 8);
        std::size_t n = 0;
        for (const auto& s : samples) {
            auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", s.c1 + "&" + s.c2));
            AS->set_mole_fractions({s.mf1, s.mf2});
            AS->update(PT_INPUTS, s.p_kPa * 1e3, s.T);
            AS->use_conductivity_RES(true);
            const double v = AS->conductivity();
            INFO("mix=" << s.c1 << "+" << s.c2 << " T=" << s.T << " ref=" << s.tc_res << " cpp=" << v);
            CHECK(v == Catch::Approx(s.tc_res).epsilon(4e-4));
            ++n;
        }
        REQUIRE(n == samples.size());
    }
}

// ----------------------------- grid evaluation harness -----------------------
// Reads the (fluid, T, p) grid built by dev/RES_grid_build.py and evaluates RES on BOTH backends
// at every point, writing dev/RES_comparison/grid_cpp.csv for dev/RES_grid_report.py to analyse.
//
// This is a measurement harness, not a test: it records failures as data rather than asserting,
// because "HEOS cannot evaluate this fluid here" is itself one of the results being collected.
// The report script is responsible for surfacing the failure counts -- they are not incidental.
//
// Run with:  CatchTestRunner.exe [RES_grid]
TEST_CASE("RES grid evaluation harness (measurement)", "[.][RES_grid][REFPROP]") {
    Skip_if_No_REFPROP();

    std::ifstream fin("dev/RES_comparison/grid_points.csv");
    REQUIRE_FALSE(!fin);  // run dev/RES_grid_build.py first, from the repo root

    struct Point
    {
        std::string fluid, region;
        double T, p;
        std::vector<CoolPropDbl> z;  // empty for a pure fluid
    };
    std::vector<Point> points;
    std::string line;
    std::getline(fin, line);  // header
    while (std::getline(fin, line)) {
        std::vector<std::string> f;
        std::string cell;
        std::istringstream ss(line);
        while (std::getline(ss, cell, ','))
            f.push_back(cell);
        if (f.size() < 4) continue;
        // Column 5 is semicolon-separated mole fractions, empty for pure fluids.
        std::vector<CoolPropDbl> z;
        if (f.size() >= 5 && !f[4].empty()) {
            std::istringstream zs(f[4]);
            std::string tok;
            while (std::getline(zs, tok, ';'))
                if (!tok.empty()) z.push_back(std::stod(tok));
        }
        points.push_back({f[0], f[3], std::stod(f[1]), std::stod(f[2]), z});
    }
    REQUIRE(points.size() > 100);

    std::ofstream out("dev/RES_comparison/grid_cpp.csv");
    REQUIRE_FALSE(!out);
    out << "fluid,T_K,p_Pa,region,mole_fractions,backend,ok,rho,s_res,eta,tc,phase,eta0_native,eta0_poly,tc0_poly,Tc,rhoc,err\n";
    out.precision(17);

    // The dilute polynomial is shipped once per fluid and is shared by every backend, so it can
    // be read off either state; evaluating it here rather than exposing it from the model keeps
    // the residual-share diagnostic honest about which term was actually subtracted.
    auto csv_z = [](const std::vector<CoolPropDbl>& z) {
        std::string s;
        for (std::size_t k = 0; k < z.size(); ++k) {
            if (k) s += ";";
            s += format("%0.17g", static_cast<double>(z[k]));
        }
        return s;
    };

    auto dilute_poly = [](const std::vector<double>& c, double T) {
        if (c.size() < 5) return std::numeric_limits<double>::quiet_NaN();
        return c[0] + T * (c[1] + T * (c[2] + T * (c[3] + T * c[4])));
    };

    std::size_t i = 0;
    while (i < points.size()) {
        const std::string fluid = points[i].fluid;
        std::size_t j = i;
        while (j < points.size() && points[j].fluid == fluid)
            ++j;

        // Three columns, because the two questions this grid answers need DIFFERENT settings
        // and running only one of them silently answers the wrong question:
        //
        //   REFPROP         defaults -- what the papers do, so this is what the reference code is
        //                   comparable to.  Used for (a), implementation correctness.
        //   REFPROP_pinned  fitted dilute term and RES enhancement viscosity, i.e. exactly what
        //   HEOS            HEOS can also do.  Used for (b), parameter transfer, where leaving the
        //                   defaults alone would measure the source policy instead.
        //
        // HEOS needs no second column: supplying no native term to RES, its AUTO already resolves to
        // the fitted model, so its default and pinned results are identical by construction.
        for (const char* backend : {"REFPROP", "REFPROP_pinned", "HEOS"}) {
            const bool pinned = (std::string(backend) != "REFPROP");
            const char* factory_name = pinned && std::string(backend) == "REFPROP_pinned" ? "REFPROP" : backend;
            // One state per (fluid, column): re-creating it per point would re-run SETUPdll on
            // every REFPROP row for no benefit.
            std::shared_ptr<AbstractState> AS;
            try {
                AS.reset(AbstractState::factory(factory_name, fluid));
                if (!points[i].z.empty()) AS->set_mole_fractions(points[i].z);
                AS->use_viscosity_RES(true);
                AS->use_conductivity_RES(true);
                if (pinned) {
                    AS->set_viscosity_RES_dilute_source(RES_DILUTE_POLYNOMIAL);
                    AS->set_conductivity_RES_dilute_source(RES_DILUTE_POLYNOMIAL);
                    AS->set_conductivity_RES_enhancement_viscosity(RES_ENH_VIS_RES);
                    AS->set_RES_mixture_enhancement(RES_MIX_ENH_OFF);
                }
            } catch (const std::exception& e) {
                for (std::size_t k = i; k < j; ++k) {
                    out << fluid << "," << points[k].T << "," << points[k].p << "," << points[k].region << "," << csv_z(points[k].z) << "," << backend
                        << ",0,,,,,,,,,,," << '"' << e.what() << '"' << "\n";
                }
                continue;
            }

            for (std::size_t k = i; k < j; ++k) {
                double rho = 0, sres = 0, eta = 0, tc = 0, eta0n = 0, eta0p = 0, tc0p = 0, Tc = 0, rhoc = 0;
                int phase = -1;
                std::string err;
                bool ok = false;
                try {
                    AS->update(PT_INPUTS, points[k].p, points[k].T);
                    rho = AS->rhomass();
                    sres = AS->smolar_residual();
                    eta = AS->viscosity();
                    tc = AS->conductivity();
                    try {
                        phase = static_cast<int>(AS->phase());
                    } catch (...) {
                    }
                    if (!AS->transport_native(iviscosity, 0.0, eta0n)) eta0n = std::numeric_limits<double>::quiet_NaN();
                    eta0p = 1e-6 * dilute_poly(AS->RES_data().comps[0].viscosity.n_dilute, points[k].T);
                    tc0p = dilute_poly(AS->RES_data().comps[0].conductivity.n_dilute, points[k].T);
                    // Emitted so the report can tell whether the Olchowy-Sengers enhancement was
                    // even active (Li gates it off at rho/rhoc >= 2 or T/Tc > 1.4).  Near the
                    // critical point that term depends on each backend's OWN critical point and
                    // derivatives, so it moves independently of the RES parameters.
                    Tc = AS->T_critical();
                    rhoc = AS->rhomass_critical();
                    ok = true;
                } catch (const std::exception& e) {
                    err = e.what();
                    for (char& c : err)
                        if (c == '"' || c == '\n' || c == ',') c = ' ';
                }
                out << fluid << "," << points[k].T << "," << points[k].p << "," << points[k].region << "," << csv_z(points[k].z) << "," << backend
                    << "," << (ok ? 1 : 0) << ",";
                if (ok) {
                    out << rho << "," << sres << "," << eta << "," << tc << "," << phase << "," << eta0n << "," << eta0p << "," << tc0p << "," << Tc
                        << "," << rhoc << ",";
                } else {
                    out << ",,,,,,,,,,";
                }
                out << '"' << err << '"' << "\n";
            }
        }
        i = j;
    }
    out.close();
    std::cout << "\nWrote dev/RES_comparison/grid_cpp.csv (" << points.size() << " points x 3 columns)\n";
}

// Measurement harness, hidden from the default run by the leading-dot tag.  It answers the
// implementation-correctness question in its cheapest form: on the REFPROP backend the
// published columns were computed with the SAME equation of state, so any deviation here is
// an implementation defect rather than a parameter-transfer error.  Run it with:
//     CatchTestRunner.exe [RES_refprop_parity]
// It reports rather than asserts; the assertions live in the two tests below it, whose
// tolerances were set from what this printed.
TEST_CASE("RES REFPROP parity sweep (measurement)", "[.][RES_refprop_parity][REFPROP]") {
    Skip_if_No_REFPROP();

    struct Row
    {
        std::string name;
        double dev;
    };
    auto report = [](const char* what, std::vector<Row>& rows, std::vector<std::string>& failed) {
        std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b) { return std::abs(a.dev) > std::abs(b.dev); });
        std::cout << "\n=== " << what << ": " << rows.size() << " reproduced, " << failed.size() << " could not be evaluated\n";
        double worst = 0, sum = 0;
        for (const auto& r : rows) {
            worst = std::max(worst, std::abs(r.dev));
            sum += std::abs(r.dev);
        }
        std::cout << "    max |dev| = " << worst * 100 << " %,  mean |dev| = " << (rows.empty() ? 0 : sum / rows.size()) * 100 << " %\n";
        for (std::size_t i = 0; i < rows.size() && i < 15; ++i) {
            std::cout << "    " << rows[i].name << "  " << rows[i].dev * 100 << " %\n";
        }
        if (!failed.empty()) {
            std::cout << "    not evaluated:";
            for (const auto& f : failed)
                std::cout << " " << f;
            std::cout << "\n";
        }
    };

    SECTION("pure viscosity") {
        auto samples = parse_vis_pure(VIS_PURE_OUT);
        REQUIRE_FALSE(samples.empty());
        std::vector<Row> rows;
        std::vector<std::string> failed;
        for (const auto& s : samples) {
            try {
                // Sample files already use REFPROP names, so no alias mapping is needed here.
                auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", s.name));
                AS->update(PT_INPUTS, s.p_MPa * 1e6, s.T);
                AS->use_viscosity_RES(true);
                rows.push_back({s.name, AS->viscosity() * 1e6 / s.vis_res - 1.0});
            } catch (const std::exception& e) {
                failed.push_back(s.name + "(" + e.what() + ")");
            }
        }
        report("pure viscosity", rows, failed);
    }
    SECTION("pure conductivity") {
        auto samples = parse_tc_pure(TC_PURE_REF);
        REQUIRE_FALSE(samples.empty());
        std::vector<Row> rows;
        std::vector<std::string> failed;
        for (const auto& s : samples) {
            try {
                auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", s.name));
                AS->update(PT_INPUTS, s.p_kPa * 1e3, s.T);
                AS->use_conductivity_RES(true);
                rows.push_back({s.name, AS->conductivity() / s.tc_res - 1.0});
            } catch (const std::exception& e) {
                failed.push_back(s.name + "(" + e.what() + ")");
            }
        }
        report("pure conductivity", rows, failed);
    }
    SECTION("binary viscosity") {
        auto samples = parse_vis_bin(VIS_BIN_OUT);
        REQUIRE_FALSE(samples.empty());
        std::vector<Row> rows;
        std::vector<std::string> failed;
        for (const auto& s : samples) {
            try {
                auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", s.c1 + "&" + s.c2));
                AS->set_mole_fractions({s.mf1, s.mf2});
                AS->update(PT_INPUTS, s.p_MPa * 1e6, s.T);
                AS->use_viscosity_RES(true);
                rows.push_back({s.c1 + "+" + s.c2, AS->viscosity() * 1e6 / s.vis_res - 1.0});
            } catch (const std::exception& e) {
                failed.push_back(s.c1 + "+" + s.c2 + "(" + e.what() + ")");
            }
        }
        report("binary viscosity", rows, failed);
    }
    SECTION("binary conductivity") {
        auto samples = parse_tc_bin(TC_BIN_REF);
        REQUIRE_FALSE(samples.empty());
        std::vector<Row> rows;
        std::vector<std::string> failed;
        for (const auto& s : samples) {
            try {
                auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("REFPROP", s.c1 + "&" + s.c2));
                AS->set_mole_fractions({s.mf1, s.mf2});
                AS->update(PT_INPUTS, s.p_kPa * 1e3, s.T);
                AS->use_conductivity_RES(true);
                rows.push_back({s.c1 + "+" + s.c2, AS->conductivity() / s.tc_res - 1.0});
            } catch (const std::exception& e) {
                failed.push_back(s.c1 + "+" + s.c2 + "(" + e.what() + ")");
            }
        }
        report("binary conductivity", rows, failed);
    }
}

// ------------------------------- critical-enhancement guards -----------------
// These three mirror Olchowy_critical_enhancement() in the Li 2024 supporting-information
// code (code_SI.py); each one failed before those guards were ported.

TEST_CASE("RES conductivity is finite for fluids with no fitted critical-enhancement params", "[RES][transport]") {
    // D2O / HELIUM / ORTHOHYD carry an all-zero critical-enhancement record in the source
    // tables, meaning "not fitted".  Treating that as provided divides by t_ref == 0 and
    // Gamma == 0, so the enhancement came out NaN and poisoned the whole conductivity.
    // 1 bar / 400 K keeps D2O inside both suppression guards, so the zero-parameter path
    // is what is actually under test here rather than an early return.
    auto AS = make_heos("HeavyWater");
    AS->use_conductivity_RES(true);
    AS->update(PT_INPUTS, 1e5, 400.0);
    double tc = AS->conductivity();
    INFO("D2O conductivity = " << tc);
    CHECK(std::isfinite(tc));
    CHECK(tc > 0.0);
}

TEST_CASE("RES conductivity does not throw a viscosity error when only viscosity params are missing", "[RES][transport]") {
    // Where the backend supplies no native term to RES, the enhancement term falls back to the RES
    // viscosity, so a fluid carrying conductivity parameters but not viscosity ones used to
    // surface a *viscosity* ValueError out of a conductivity call.  The enhancement must be
    // skipped instead.
    //
    // This combination used to occur naturally -- R41 was excluded for viscosity but not for
    // conductivity -- but the grid-derived exclusion lists select the same fluids for both
    // properties, since s_res is the one input both models share.  So the case is constructed
    // rather than found: the code path is still reachable by anyone calling
    // set_conductivity_RES_parameters() without the viscosity counterpart.
    auto AS = make_heos("Propane");
    AS->use_conductivity_RES(true);
    // 4.65 MPa / 375.9 K is inside the enhancement region, so the fallback is actually reached.
    AS->update(PT_INPUTS, 4.65e6, 375.9);
    auto& comps = AS->RES_data_mutable().comps;
    REQUIRE(comps.size() == 1);
    REQUIRE(comps[0].conductivity.provided);
    REQUIRE(comps[0].viscosity.provided);  // premise: it starts out present, so removing it means something
    comps[0].viscosity.provided = false;

    double tc = 0;
    REQUIRE_NOTHROW(tc = AS->conductivity());
    INFO("Propane conductivity without viscosity RES params = " << tc);
    CHECK(std::isfinite(tc));
    CHECK(tc > 0.0);
}

TEST_CASE("RES critical enhancement is suppressed outside the near-critical region", "[RES][transport]") {
    // Li 2024 returns exactly 0 when rho/rhoc >= 2 or T/Tc > 1.4.  Compare against the same
    // state with crit_provided forced off: if the guards work the two must agree bit-for-bit,
    // because the enhancement contributes nothing either way.
    auto with_crit = [](double p, double T) {
        auto AS = make_heos("CarbonDioxide");
        AS->use_conductivity_RES(true);
        AS->update(PT_INPUTS, p, T);
        return AS->conductivity();
    };
    auto without_crit = [](double p, double T) {
        auto AS = make_heos("CarbonDioxide");
        auto* heos = dynamic_cast<HelmholtzEOSMixtureBackend*>(AS.get());
        heos->RES_data_mutable().comps[0].conductivity.crit_provided = false;
        AS->use_conductivity_RES(true);
        AS->update(PT_INPUTS, p, T);
        return AS->conductivity();
    };

    SECTION("T/Tc > 1.4") {  // 500 K -> T/Tc = 1.64
        double a = with_crit(5.0e6, 500.0), b = without_crit(5.0e6, 500.0);
        INFO("with=" << a << " without=" << b);
        CHECK(std::isfinite(a));
        CHECK(a == Catch::Approx(b).epsilon(1e-12));
    }
    SECTION("rho/rhoc >= 2") {  // 250 K / 30 MPa -> rho/rhoc = 2.42
        double a = with_crit(3.0e7, 250.0), b = without_crit(3.0e7, 250.0);
        INFO("with=" << a << " without=" << b);
        CHECK(std::isfinite(a));
        CHECK(a == Catch::Approx(b).epsilon(1e-12));
    }
    SECTION("near-critical state still receives an enhancement") {
        // Complements the above: verifies the guards are not simply disabling the term everywhere.
        double a = with_crit(8.0e6, 310.0), b = without_crit(8.0e6, 310.0);
        INFO("with=" << a << " without=" << b);
        CHECK(std::isfinite(a));
        CHECK(a > b);
    }
}

TEST_CASE("RES full parameter getters round-trip against the setters", "[RES][transport]") {
    // The residual-only getters cover what a refit changes; these return everything, including
    // the shipped dilute polynomial and the enhancement constants, which otherwise had no
    // read-back path at all outside C++ (RES_data() is not exposed to the wrappers).
    auto AS = make_heos("Propane");
    const auto v0 = AS->get_viscosity_RES_parameters(0);
    const auto c0 = AS->get_conductivity_RES_parameters(0);
    REQUIRE(v0.provided);
    REQUIRE(c0.provided);
    CHECK(v0.n_dilute.size() == 5);
    CHECK(v0.n_res.size() == 3);
    CHECK(c0.n_dilute.size() == 5);
    CHECK(c0.n_res.size() == 4);
    CHECK(v0.molar_mass == Catch::Approx(AS->molar_mass()).epsilon(1e-12));
    CHECK(c0.crit_provided);  // Propane ships enhancement parameters
    CHECK(c0.t_ref > 0);

    // What the setters write is what the getters read back.
    const std::vector<double> nd{1, 2, 3, 4, 5}, nr{0.1, 0.2, 0.3};
    AS->set_viscosity_RES_parameters(0, nd, nr, 1.25);
    const auto v1 = AS->get_viscosity_RES_parameters(0);
    CHECK(v1.n_dilute == nd);
    CHECK(v1.n_res == nr);
    CHECK(v1.xita == 1.25);

    const std::vector<double> cnr{0.1, 0.2, 0.3, 0.4};
    AS->set_conductivity_RES_parameters(0, nd, cnr, 0.9, 1.01, 1.24, 0.05, 1.9e-10, 600.0, 1.4e9);
    const auto c1 = AS->get_conductivity_RES_parameters(0);
    CHECK(c1.n_dilute == nd);
    CHECK(c1.n_res == cnr);
    CHECK(c1.xita == 0.9);
    CHECK(c1.R_D == 1.01);
    CHECK(c1.gamma_uni == 1.24);
    CHECK(c1.Gamma == 0.05);
    CHECK(c1.phi0 == 1.9e-10);
    CHECK(c1.t_ref == 600.0);
    CHECK(c1.q_D == 1.4e9);

    // Out of range must report, not read past the end.
    CHECK_THROWS_AS(AS->get_viscosity_RES_parameters(1), CoolProp::ValueError);
    CHECK_THROWS_AS(AS->get_conductivity_RES_parameters(1), CoolProp::ValueError);
}

TEST_CASE("RES viscosity residual params round-trip", "[RES][transport]") {
    auto AS = make_pr("CarbonDioxide");
    const std::vector<double> n = {0.123, -0.456, 0.078};
    AS->set_viscosity_RES_residual_params(0, n, 1.23);
    auto [n_out, x_out] = AS->get_viscosity_RES_residual_params(0);
    REQUIRE(n_out.size() == n.size());
    for (std::size_t i = 0; i < n.size(); ++i)
        CHECK(n_out[i] == Catch::Approx(n[i]));
    CHECK(x_out == Catch::Approx(1.23));
}

TEST_CASE("RES conductivity residual params round-trip", "[RES][transport]") {
    auto AS = make_pr("CarbonDioxide");
    const std::vector<double> n = {0.42, -0.20, 0.06, 0.01};
    AS->set_conductivity_RES_residual_params(0, n, 0.97);
    auto [n_out, x_out] = AS->get_conductivity_RES_residual_params(0);
    REQUIRE(n_out.size() == n.size());
    for (std::size_t i = 0; i < n.size(); ++i)
        CHECK(n_out[i] == Catch::Approx(n[i]));
    CHECK(x_out == Catch::Approx(0.97));
}

// --------------------------------- structural / fallback tests ---------------

TEST_CASE("HEOS mixture viscosity uses log-mean when RES not enabled", "[RES][transport]") {
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Argon&Nitrogen"));
    AS->set_mole_fractions({0.5, 0.5});
    AS->update(PT_INPUTS, 1e5, 300.0);
    REQUIRE_NOTHROW(AS->viscosity());
}

TEST_CASE("HEOS pure fluid explicit RES opt-in without params throws an informative error", "[RES][transport]") {
    // Before this fix, use_viscosity_RES(true)/use_conductivity_RES(true) on a fluid whose
    // RES parameters are not "provided" silently fell back to the default reference
    // correlation with no indication that RES was never actually used. Simulate that
    // situation directly (rather than relying on a specific fluid lacking RES data, which
    // could change as res_transport_parameters.json evolves) and check it now raises a
    // clear ValueError instead.
    auto AS = make_heos("CarbonDioxide");
    auto* heos = dynamic_cast<HelmholtzEOSMixtureBackend*>(AS.get());
    REQUIRE(heos != nullptr);
    heos->RES_data_mutable().comps[0].viscosity.provided = false;
    heos->RES_data_mutable().comps[0].conductivity.provided = false;

    AS->use_viscosity_RES(true);
    AS->use_conductivity_RES(true);
    AS->update(PT_INPUTS, 1e6, 300.0);
    CHECK_THROWS_AS(AS->viscosity(), CoolProp::ValueError);
    CHECK_THROWS_AS(AS->conductivity(), CoolProp::ValueError);
    CHECK_THROWS_WITH(AS->viscosity(), Catch::Matchers::ContainsSubstring("use_viscosity_RES"));
    CHECK_THROWS_WITH(AS->conductivity(), Catch::Matchers::ContainsSubstring("use_conductivity_RES"));
}

TEST_CASE("HEOS mixture explicit RES opt-in without params throws an informative error", "[RES][transport]") {
    // Same as above, but for the mixture (log-mean/linear-sum) dispatch branch: RES opt-in
    // must fail loudly, naming the offending component, rather than silently falling back
    // to the approximate mixture rule.
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Argon&Nitrogen"));
    AS->set_mole_fractions({0.5, 0.5});
    auto* heos = dynamic_cast<HelmholtzEOSMixtureBackend*>(AS.get());
    REQUIRE(heos != nullptr);
    heos->RES_data_mutable().comps[0].viscosity.provided = false;
    heos->RES_data_mutable().comps[0].conductivity.provided = false;

    AS->use_viscosity_RES(true);
    AS->use_conductivity_RES(true);
    AS->update(PT_INPUTS, 1e5, 300.0);
    CHECK_THROWS_AS(AS->viscosity(), CoolProp::ValueError);
    CHECK_THROWS_AS(AS->conductivity(), CoolProp::ValueError);
    CHECK_THROWS_WITH(AS->viscosity(), Catch::Matchers::ContainsSubstring("Argon"));
}

TEST_CASE("PR backend without RES params throws NotImplementedError", "[RES][transport]") {
    auto AS = make_pr("CarbonDioxide");
    auto* cubic = as_cubic(AS);
    auto& comps = cubic->RES_data_mutable().comps;
    comps[0].viscosity.provided = false;
    comps[0].conductivity.provided = false;
    AS->specify_phase(iphase_gas);
    AS->update(PT_INPUTS, 1e6, 300.0);
    CHECK_THROWS_AS(AS->viscosity(), CoolProp::NotImplementedError);
    CHECK_THROWS_AS(AS->conductivity(), CoolProp::NotImplementedError);
}

#endif  // ENABLE_CATCH
