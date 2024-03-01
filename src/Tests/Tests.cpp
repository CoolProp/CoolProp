/**
This file includes some testing functions that will get built
into the program.  Otherwise CTest can be used by removing this file from
the build to avoid double declaration of the main function and
Catch clashing
*/
#include "Tests.h"
#include <ctime>

#if defined ENABLE_CATCH
#    include <catch2/catch_all.hpp>

static Catch::Session session;  // There must be exactly one instance

#endif  // ENABLE_CATCH

int run_fast_tests() {
#ifdef ENABLE_CATCH
    Catch::ConfigData& config = session.configData();
    config.testsOrTags.clear();
    config.testsOrTags.emplace_back("[fast]");
    session.useConfigData(config);
    return session.run();
#else
    return 0;
#endif
}

int run_not_slow_tests() {
#ifdef ENABLE_CATCH
    Catch::ConfigData& config = session.configData();
    config.testsOrTags.clear();
    config.testsOrTags.emplace_back("~[slow]");
    session.useConfigData(config);

    time_t t1 = clock();
    session.run();
    time_t t2 = clock();
    printf("Elapsed time for not slow tests: %g s", (double)(t2 - t1) / CLOCKS_PER_SEC);

    return 1;
#else
    return 0;
#endif
}

int run_user_defined_tests(const std::vector<std::string>& tests_or_tags) {
#ifdef ENABLE_CATCH
    Catch::ConfigData& config = session.configData();
    config.testsOrTags = tests_or_tags;
    session.useConfigData(config);

    time_t t1 = clock();
    session.run();
    time_t t2 = clock();
    printf("Elapsed time for user defined tests: %g s", (double)(t2 - t1) / CLOCKS_PER_SEC);

    return 1;
#else
    return 0;
#endif
}

void run_tests() {
#ifdef ENABLE_CATCH
    Catch::ConfigData& config = session.configData();
    config.testsOrTags.clear();
    //config.shouldDebugBreak = true;
    session.useConfigData(config);
    session.run();
#endif
}
