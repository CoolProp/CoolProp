

#include "IF97Backend.h"
#include "AbstractState.h"
#include "DataStructures.h"

class IF97BackendGenerator : public CoolProp::AbstractStateGenerator{
public:
    CoolProp::AbstractState * get_AbstractState(const std::vector<std::string> &fluid_names){
        return new CoolProp::IF97Backend();
    };
} ;
// This static initialization will cause the generator to register
static CoolProp::GeneratorInitializer<CoolProp::IF97_BACKEND_FAMILY, IF97BackendGenerator> if97_gen;