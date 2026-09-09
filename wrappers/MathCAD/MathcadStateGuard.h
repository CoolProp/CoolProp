// MathcadStateGuard.h : a Mathcad-wrapper-specific convenience layered on top
// of CoolProp's handle-based low-level C API (include/CoolProp/CoolPropLib.h).
//
// Deliberately free of any Mathcad SDK (mcadincl.h) dependency, even though
// CoolPropMathcad.cpp is currently its only caller -- keeps this class's
// logic decoupled from Mathcad-specific types and buildable/testable on its
// own if that's ever useful, without requiring the Mathcad Prime SDK.

#ifndef MATHCAD_STATE_GUARD_H
#define MATHCAD_STATE_GUARD_H

#include <map>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

// Gives the Mathcad wrapper's AS_factory() function "get-or-create" (memoized)
// semantics on top of AbstractState_factory(): calling get_or_create() again
// with the same (backend, fluids) pair returns the SAME live handle that
// pair returned last time, instead of destroying and rebuilding the backend
// -- rebuilding on every Mathcad recalculation (every worksheet edit/
// recalculation re-executes AS_factory()'s cell) would pay the backend's
// full construction cost -- 80-140 ms for tabular backends (BICUBIC/TTSE)
// -- every single time, defeating the entire reason the Low-Level API exists
// (build once, reuse for many flashes). A *different* (backend, fluids) pair
// gets its own independent entry, so multi-fluid worksheets (e.g. two
// independent states for two sides of a heat exchanger) are unaffected, and
// nothing is ever destroyed just because a *different* key was requested --
// only a genuinely dead handle under the SAME key (see get_or_create()'s own
// comment) triggers creating a replacement.
//
// This must NOT change AbstractStateLibrary's shared semantics in
// CoolPropLib.cpp -- other callers of the low-level C API (Fortran, Julia,
// ...) legitimately want multiple concurrent handles for the same
// backend+fluids, so the memoization policy lives here, one layer up, rather
// than in the shared handle table itself.
class MathcadStateGuard
{
   public:
    // Returns the live handle already registered for this (backend, fluids)
    // key, or creates and registers a new one via AbstractState_factory() if
    // there isn't one (or the registered one turned out to be dead -- e.g.
    // freed directly via AS_free() outside this registry). On reuse, also
    // clears any phase constraint a *previous* worksheet state may have
    // imposed via AbstractState_specify_phase(), so removing/changing an
    // AS_specify_phase() call in the worksheet can't leave a stale
    // constraint silently in effect on the reused handle -- this is the one
    // piece of AbstractState configuration state cheap enough to reset
    // without a full rebuild. (Mixture fractions are NOT reset: there is no
    // cheap "clear fractions" call, and AS_set_fractions() is always
    // re-chained immediately after AS_factory() in normal use, so it
    // naturally re-applies on every recalculation regardless.)
    //
    // Same errcode/message_buffer/buffer_length contract as the
    // AbstractState_* functions in CoolPropLib.h: *errcode == 0 on success;
    // message_buffer/buffer_length follow the same "caller-owned buffer,
    // truncate if it doesn't fit" convention used throughout that API.
    long get_or_create(const std::string& backend, const std::string& fluids, long* errcode, char* message_buffer, long buffer_length);

    // Returns a point-in-time copy of every currently-registered (key, handle)
    // pair whose handle is still actually alive, in key order -- i.e. the
    // same order std::map<std::string, long> iterates in, which depends only
    // on the CURRENT set of keys, not on when the copy is taken. That
    // property is what lets two independent callers (AS_list_handles() and
    // AS_list_states() in MathcadLowLevel.h, each calling snapshot()
    // separately since a Mathcad function can only return one of
    // COMPLEX_ARRAY/MC_STRING, never both) agree on ordering as long as no
    // get_or_create()/erase happens between the two calls -- the normal case
    // for two list-cells on the same worksheet. See those functions'
    // comments for the one edge case where that assumption doesn't hold.
    //
    // A handle registered here can go dead without this class knowing --
    // AS_free() in MathcadLowLevel.h releases handles directly via
    // AbstractState_free(), deliberately bypassing this registry (see that
    // function's own comment), so it has no way to remove the now-stale
    // entry at free time. snapshot() is where that staleness actually gets
    // noticed and cleaned up: like get_or_create(), it probes each handle
    // before reporting it, and silently drops (erases) any that are dead.
    // Not const, since it can mutate `live` this way.
    std::vector<std::pair<std::string, long>> snapshot();

   private:
    std::map<std::string, long> live;  // key: backend + "|" + fluids -> the live handle for that key, if any
    mutable std::mutex mtx;            // mutable: snapshot() is logically read-only but still needs to lock
};

#endif
