// Catch2 tests that number parsing outside the expression DSL ignores the C
// locale (tag [locale]).  A host program (Python, EES, Mathcad, ...) that calls
// setlocale(LC_NUMERIC, "de_DE") made strtod, atof and std::stod stop at the
// '.'; the DSL's own locale tests are in CoolProp-Tests-Expression.cpp.

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include "LocaleGuard.h"

#    include <cmath>
#    include <cstdlib>
#    include <string>
#    include <system_error>

#    include "CoolProp/CoolProp.h"
#    include "CoolProp/Configuration.h"
#    include "CoolProp/Exceptions.h"
#    include "CoolProp/detail/strings.h"

namespace {

struct Parsed
{
    std::errc ec;
    double value;
    std::size_t length;
};

Parsed parse(const std::string& s) {
    double v = -99.0;
    const char* end = nullptr;
    const std::errc ec = parse_double_C(s.data(), s.data() + s.size(), v, end);
    return {ec, v, static_cast<std::size_t>(end - s.data())};
}

// Sets an environment variable for the lifetime of the guard.
class EnvGuard
{
    std::string name_;

   public:
    EnvGuard(const char* name, const char* value) : name_(name) {
#    if defined(_WIN32)
        (void)_putenv_s(name, value);
#    else
        (void)setenv(name, value, 1);
#    endif
    }
    ~EnvGuard() {
#    if defined(_WIN32)
        (void)_putenv_s(name_.c_str(), "");
#    else
        (void)unsetenv(name_.c_str());
#    endif
    }
    EnvGuard(const EnvGuard&) = delete;
    EnvGuard& operator=(const EnvGuard&) = delete;
    EnvGuard(EnvGuard&&) = delete;
    EnvGuard& operator=(EnvGuard&&) = delete;
};

}  // namespace

TEST_CASE("parse_double_C ignores the C locale", "[locale]") {
    NumericLocaleGuard guard("de_DE.UTF-8");
    require_decimal_comma_locale(guard);

    Parsed p = parse("8.4e3");
    CHECK(p.ec == std::errc());
    CHECK(p.value == 8400.0);
    CHECK(p.length == 5);

    p = parse("+1.5%");
    CHECK(p.ec == std::errc());
    CHECK(p.value == 1.5);
    CHECK(p.length == 4);  // stops before the '%'

    p = parse("-2.5e-1");
    CHECK(p.ec == std::errc());
    CHECK(p.value == -0.25);

    p = parse("0x10");  // hex is not accepted: reads the 0 only
    CHECK(p.ec == std::errc());
    CHECK(p.value == 0.0);
    CHECK(p.length == 1);

    for (const char* bad : {"", "+-1", "+", ".", "abc", " 1"}) {
        INFO(bad);
        p = parse(bad);
        CHECK(p.ec == std::errc::invalid_argument);
        CHECK(p.length == 0);
        CHECK(p.value == -99.0);
    }

    p = parse("1e400");
    CHECK(p.ec == std::errc::result_out_of_range);
    CHECK(p.length == 5);
    CHECK(p.value == -99.0);
}

TEST_CASE("string2double ignores the C locale", "[locale]") {
    NumericLocaleGuard guard("de_DE.UTF-8");
    require_decimal_comma_locale(guard);

    CHECK(string2double("1.5") == 1.5);
    CHECK(string2double("1.5D-3") == 1.5e-3);  // FORTRAN exponent, as in HMX.BNC
    CHECK(string2double("-2.25d2") == -225.0);
    CHECK(string2double("  0.125") == 0.125);  // leading whitespace, as strtod allowed
    CHECK(string2double("+7") == 7.0);
    CHECK_THROWS_AS(string2double("1,5"), CoolProp::ValueError);
    CHECK_THROWS_AS(string2double("1.5x"), CoolProp::ValueError);
    CHECK_THROWS_AS(string2double(""), CoolProp::ValueError);
    CHECK_THROWS_AS(string2double("1e400"), CoolProp::ValueError);
}

// Under a decimal-comma strtod, "MEG-20.5%" read as 20 then ".5%", which is not
// "%", so the mass fraction silently became 20.
TEST_CASE("incompressible solution concentration ignores the C locale", "[locale]") {
    const double rho_ref = CoolProp::PropsSI("D", "T", 300.0, "P", 101325.0, "INCOMP::MEG-20.5%");
    REQUIRE(std::isfinite(rho_ref));
    CHECK(rho_ref != CoolProp::PropsSI("D", "T", 300.0, "P", 101325.0, "INCOMP::MEG-20%"));

    double rho_de = 0.0;
    {
        NumericLocaleGuard guard("de_DE.UTF-8");
        require_decimal_comma_locale(guard);
        rho_de = CoolProp::PropsSI("D", "T", 300.0, "P", 101325.0, "INCOMP::MEG-20.5%");
    }
    CHECK(rho_de == rho_ref);

    // Malformed concentrations are errors rather than a silent 0 or an unscaled percentage.
    CHECK_FALSE(std::isfinite(CoolProp::PropsSI("D", "T", 300.0, "P", 101325.0, "INCOMP::MEG-abc%")));
    CHECK_FALSE(std::isfinite(CoolProp::PropsSI("D", "T", 300.0, "P", 101325.0, "INCOMP::MEG-20 %")));
}

TEST_CASE("double-valued COOLPROP_* environment variables ignore the C locale", "[locale]") {
    NumericLocaleGuard guard("de_DE.UTF-8");
    require_decimal_comma_locale(guard);
    {
        EnvGuard env("COOLPROP_SPINODAL_MINIMUM_DELTA", "0.25");
        CoolProp::Configuration config;  // reads COOLPROP_* in its constructor
        CHECK(static_cast<double>(config.get_item(SPINODAL_MINIMUM_DELTA)) == 0.25);
    }
    {
        EnvGuard env("COOLPROP_SPINODAL_MINIMUM_DELTA", "0.25abc");
        CHECK_THROWS_AS(CoolProp::Configuration(), CoolProp::ValueError);
    }
}

#endif
