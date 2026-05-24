// Shared helpers for Catch2 test TUs.  Header-only.

#ifndef COOLPROP_TESTS_TEST_UTILS_H
#define COOLPROP_TESTS_TEST_UTILS_H

#if defined(_WIN32)
#    include <process.h>  // _getpid
#else
#    include <unistd.h>  // getpid
#endif

namespace CoolProp {
namespace tests {

// Portable process ID for tests that stage scratch files under
// fs::temp_directory_path() — suffix the path with this so two
// concurrent CatchTestRunner instances (preflight in two worktrees,
// dev + CI, etc.) don't collide on the same hardcoded path.  See
// CoolProp-8ft.
inline int test_pid() {
#if defined(_WIN32)
    return ::_getpid();
#else
    return ::getpid();
#endif
}

}  // namespace tests
}  // namespace CoolProp

#endif  // COOLPROP_TESTS_TEST_UTILS_H
