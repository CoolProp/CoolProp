#ifndef COOLPROP_FILE_PATH_H
#define COOLPROP_FILE_PATH_H

#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

/// Get directory separator
std::string get_separator();

/// Get the user's home directory;  It is believed that is is always a place that files can be written
std::string get_home_dir();

/// Return true if path exists
bool path_exists(const std::string& path);

/// Return merged path, append separator if string two is empty
std::string join_path(const std::string& one, const std::string& two);

/// Make directory and all required intermediate directories
void make_dirs(std::string file_path);

/// Get the size of a directory in bytes.
/// Uses std::filesystem on modern platforms (Windows, Linux, macOS).
/// Legacy Android/powerpc stub returns 0.
unsigned long long CalculateDirSize(const std::string& path);

// Get all the contents of a file and dump into a STL string
// Thanks to http://stackoverflow.com/questions/2602013/read-whole-ascii-file-into-c-stdstring
std::string get_file_contents(const char* filename);

/// Get all the contents of a binary file
std::vector<char> get_binary_file_contents(const char* filename);

/// Atomically write `size` bytes from `bytes` to `target`.
///
/// Writes to a sibling `<target>.tmp.<process-salt>.<seq>` first, then
/// `std::filesystem::rename`s onto the target.  The rename is atomic on
/// POSIX same-filesystem and on Win32 (replace-existing semantics), so
/// concurrent writers in different processes / threads either see the
/// previous complete file or one complete writer's payload — never a
/// partial-write file on the visible path.
///
/// The temp-path suffix combines a process-unique 64-bit salt (drawn
/// once per process from `std::random_device`) with a process-local
/// atomic counter, so concurrent writers never collide on their
/// intermediates.
///
/// When `restrict_perms` is true, the temp file is chmod'd to
/// owner-only (0600) *before* the rename, so the post-rename file
/// inherits the tight permissions without a brief world-readable
/// window.  No-op on Windows where the POSIX permission model doesn't
/// apply.
///
/// Throws `std::runtime_error` on open/write/rename failure; the temp
/// file is unlinked before throwing.  Assumes `target.parent_path()`
/// already exists.
void write_bytes_atomic(const std::filesystem::path& target, const void* bytes, std::size_t size, bool restrict_perms = false);

#endif