#ifndef COOLPROP_DETAIL_ATOMIC_WRITE_H
#define COOLPROP_DETAIL_ATOMIC_WRITE_H

// Narrow home of write_bytes_atomic, split out of detail/filepaths.h (GH:
// de-leak <filesystem>).  filepaths.h is pulled into the public surface via
// detail/tools.h -> ... -> AbstractState.h; keeping the only std::filesystem
// signature here means including AbstractState.h no longer requires <filesystem>
// (and thus C++17 for that reason).  Include this header where you actually
// call write_bytes_atomic.

#include <cstddef>
#include <filesystem>

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

#endif  // COOLPROP_DETAIL_ATOMIC_WRITE_H
