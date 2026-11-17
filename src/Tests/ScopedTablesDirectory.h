// Test helper: point the tabular backends at a fresh temporary table directory
// for the lifetime of the object.
//
// The on-disk table path is built from the backend and the composition only,
// not the grid size (TabularBackend::path_to_tables).  A test that builds a
// non-default grid in the shared directory therefore overwrites the
// default-grid tables there, and concurrent test processes can write the same
// files.  A unique directory also gives the dataset a unique in-memory cache
// key, so the tables are built cold.
//
// Catch2 re-runs a TEST_CASE body for each SECTION, so an instance at the top
// of a test case gives every section its own directory.

#ifndef COOLPROP_TESTS_SCOPED_TABLES_DIRECTORY_H
#define COOLPROP_TESTS_SCOPED_TABLES_DIRECTORY_H

#include <atomic>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <string>
#include <system_error>

#include "CoolProp/Configuration.h"

namespace CoolProp {
namespace testing {

class ScopedTablesDirectory
{
   public:
    explicit ScopedTablesDirectory(const std::string& tag) : prev_(get_config_string(ALTERNATIVE_TABLES_DIRECTORY)) {
        static std::atomic<unsigned long> counter{0};
        const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
        dir_ = (std::filesystem::temp_directory_path() / ("CoolProp-" + tag + "-" + std::to_string(stamp) + "-" + std::to_string(counter++))).string()
               + "/";
        set_config_string(ALTERNATIVE_TABLES_DIRECTORY, dir_);
    }
    ~ScopedTablesDirectory() {
        // A destructor must not throw (bad_alloc in either call would terminate)
        try {
            set_config_string(ALTERNATIVE_TABLES_DIRECTORY, prev_);
            std::error_code ec;
            std::filesystem::remove_all(dir_, ec);
        } catch (...) {
            // Cannot rethrow from a destructor; say so instead of terminating
            (void)std::fputs("ScopedTablesDirectory: could not restore ALTERNATIVE_TABLES_DIRECTORY or remove the temp directory\n", stderr);
        }
    }
    ScopedTablesDirectory(const ScopedTablesDirectory&) = delete;
    ScopedTablesDirectory& operator=(const ScopedTablesDirectory&) = delete;
    ScopedTablesDirectory(ScopedTablesDirectory&&) = delete;
    ScopedTablesDirectory& operator=(ScopedTablesDirectory&&) = delete;

   private:
    std::string prev_;
    std::string dir_;
};

}  // namespace testing
}  // namespace CoolProp

#endif
