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
        long probe_errcode = 0;
        char probe_message[256];
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
