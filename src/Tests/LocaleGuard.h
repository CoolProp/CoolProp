// Shared helper for Catch2 tests that switch to a decimal-comma C locale
// (tag [locale]).  Header-only; include after catch2/catch_all.hpp.

#ifndef COOLPROP_TESTS_LOCALE_GUARD_H
#define COOLPROP_TESTS_LOCALE_GUARD_H

#include <clocale>
#include <cstdlib>
#include <string>

// Switches LC_NUMERIC for the lifetime of the guard and restores the previous
// setting on scope exit, including when a REQUIRE throws.  `active()` is false
// when the host has no such locale installed; callers SKIP in that case.
class NumericLocaleGuard
{
    static std::string current() {
        const char* cur = std::setlocale(LC_NUMERIC, nullptr);
        return (cur != nullptr) ? cur : "C";
    }
    // Declaration order matters: saved_ is captured before active_ switches.
    std::string saved_;
    bool active_;

   public:
    explicit NumericLocaleGuard(const char* name) : saved_(current()), active_(std::setlocale(LC_NUMERIC, name) != nullptr) {}
    ~NumericLocaleGuard() {
        (void)std::setlocale(LC_NUMERIC, saved_.c_str());
    }
    NumericLocaleGuard(const NumericLocaleGuard&) = delete;
    NumericLocaleGuard& operator=(const NumericLocaleGuard&) = delete;
    NumericLocaleGuard(NumericLocaleGuard&&) = delete;
    NumericLocaleGuard& operator=(NumericLocaleGuard&&) = delete;
    [[nodiscard]] bool active() const {
        return active_;
    }
};

// CI installs de_DE.UTF-8 and sets COOLPROP_REQUIRE_LOCALE_TESTS=1, so a missing
// locale there is a failure, not a silent SKIP that would let a regression through.
// Catch2's SKIP/FAIL/REQUIRE throw, so they end the calling test from here too.
inline void require_decimal_comma_locale(const NumericLocaleGuard& guard) {
    if (!guard.active()) {
        if (std::getenv("COOLPROP_REQUIRE_LOCALE_TESTS") != nullptr) {
            FAIL("de_DE.UTF-8 locale required (COOLPROP_REQUIRE_LOCALE_TESTS) but absent");
        }
        SKIP("de_DE.UTF-8 locale not installed");
    }
    REQUIRE(std::localeconv()->decimal_point[0] == ',');
}

#endif
