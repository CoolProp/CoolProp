
#include <vector>
#include <string>
#include <chrono>
#include <memory>

#include "CoolProp.h"
#include "AbstractState.h"


int main(int argc, const char* argv[]) {

    CoolProp::set_debug_level(1000);
    std::shared_ptr<CoolProp::AbstractState> ptr;

    std::string backend;
    std::vector<std::string> fluids;
    double Q, T, p, s, h;

    // Test as described in https://github.com/CoolProp/CoolProp/issues/1611
    backend = "HEOS";
    fluids = {"R407C"};
    ptr.reset(CoolProp::AbstractState::factory(backend, fluids));
    p = 4863285.0;
    Q = 0;
    ptr->update(CoolProp::PQ_INPUTS, p, Q);

    // test as described in https://github.com/CoolProp/CoolProp/issues/1678
    backend = "HEOS";
    fluids = {"Water"};
    ptr.reset(CoolProp::AbstractState::factory(backend, fluids));
    p = ptr->p_critical();
    Q = 0;
    ptr->update(CoolProp::PQ_INPUTS, p, Q);
    s = 4000;
    ptr->update(CoolProp::PSmass_INPUTS, p, s);


    return 0;
}
