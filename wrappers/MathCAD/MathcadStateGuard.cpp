#include "MathcadStateGuard.h"

#include "CoolProp/CoolPropLib.h"

#include <cstring>

long MathcadStateGuard::get_or_create(const std::string& backend, const std::string& fluids, long* errcode, char* message_buffer,
                                      long buffer_length) {
    std::scoped_lock guard(mtx);

    const std::string key = backend + "|" + fluids;
    auto it = live.find(key);
    if (it != live.end()) {
        // Probe-and-reset in one call: AbstractState_unspecify_phase() both
        // clears any stale phase constraint AND tells us, via the errcode,
        // whether the handle we're tracking is still actually alive --
        // internally it does handle_manager.get(handle) before touching the
        // object, so a dead handle (freed directly via AS_free() outside
        // this registry) reports the same "HandleError:"-prefixed message
        // AbstractState_free()/AbstractState_factory() use elsewhere in this
        // API for that case.
        //
        // A live handle whose backend simply doesn't implement phase
        // specification (NotImplementedError, not a HandleError) is NOT a
        // reason to discard it -- there's nothing to clear, not a failure.
        // Zero-initialized: HandleException() (src/CoolPropLib.cpp) only
        // writes message_buffer when the formatted error text fits in it --
        // on its "didn't fit" path (errcode==2) the buffer is left as-is, so
        // an uninitialized array here could leave strncmp() below reading
        // garbage stack memory instead of a real (or empty) message.
        long probe_errcode = 0;
        char probe_message[256] = {};
        AbstractState_unspecify_phase(it->second, &probe_errcode, probe_message, static_cast<long>(sizeof(probe_message)));
        if (probe_errcode != 0 && std::strncmp(probe_message, "HandleError:", 12) == 0) {
            live.erase(it);  // dead -- fall through to creating a fresh one below
        } else {
            *errcode = 0;
            return it->second;
        }
    }

    long handle = AbstractState_factory(backend.c_str(), fluids.c_str(), errcode, message_buffer, buffer_length);
    if (*errcode == 0) {
        live[key] = handle;
    }
    return handle;
}

std::vector<std::pair<std::string, long>> MathcadStateGuard::snapshot() {
    std::scoped_lock guard(mtx);

    std::vector<std::pair<std::string, long>> result;
    result.reserve(live.size());

    for (auto it = live.begin(); it != live.end();) {
        // Read-only aliveness probe. Unlike get_or_create()'s reuse-time
        // probe (AbstractState_unspecify_phase(), which deliberately also
        // clears a stale phase constraint as part of reusing the handle),
        // this call must not have any side effect on a state it is merely
        // listing -- AbstractState_backend_name() only reads the object
        // (handle_manager.get(handle) then AS->backend_name()), so it's safe
        // to use purely as a liveness check.
        // See the matching comment in get_or_create() above -- same
        // zero-init reasoning applies to this probe's message buffer.
        long probe_errcode = 0;
        char probe_message[256] = {};
        char backend_buf[256] = {};
        AbstractState_backend_name(it->second, backend_buf, &probe_errcode, probe_message, static_cast<long>(sizeof(probe_message)));
        if (probe_errcode != 0 && std::strncmp(probe_message, "HandleError:", 12) == 0) {
            it = live.erase(it);  // dead -- freed directly via AS_free() outside this registry
        } else {
            result.emplace_back(it->first, it->second);
            ++it;
        }
    }

    return result;
}
