

#include "IF97Backend.h"
#include "AbstractState.h"
#include "DataStructures.h"

static class IF97BackendGenerator : public CoolProp::AbstractStateGenerator{
public:
    IF97BackendGenerator(){
        register_backend(CoolProp::IF97_BACKEND_FAMILY, shared_ptr<AbstractStateGenerator>(this));
    }
    CoolProp::AbstractState * get_AbstractState(const std::vector<std::string> &fluid_names){
        return new CoolProp::IF97Backend();
    };
} if97_gen; // This static initialization will cause the generator to register