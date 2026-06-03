// DEPRECATED forwarding shim (GH #1280). Remove at v9.
// The canonical location is <CoolProp/detail/rapidjson.h>.
#ifndef COOLPROP_NO_DEPRECATED_HEADER_WARNINGS
#    pragma message( \
      "rapidjson_include.h is deprecated; include \"CoolProp/detail/rapidjson.h\". Define COOLPROP_NO_DEPRECATED_HEADER_WARNINGS to silence.")
#endif
#include "CoolProp/detail/rapidjson.h"
