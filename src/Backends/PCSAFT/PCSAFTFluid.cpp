#include <string>
#include <vector>
#include <map>

#include "PCSAFTFluid.h"
#include "rapidjson_include.h"

namespace CoolProp {

PCSAFTFluid::PCSAFTFluid(rapidjson::Value::ValueIterator itr) {
    name = cpjson::get_string(*itr, "name");
    CAS = cpjson::get_string(*itr, "CAS");
    params.m = cpjson::get_double(*itr, "m");
    params.sigma = cpjson::get_double(*itr, "sigma");
    params.u = cpjson::get_double(*itr, "u");
    if (itr->HasMember("uAB") && (*itr)["uAB"].IsNumber()){
        params.uAB = cpjson::get_double(*itr, "uAB");
    }
    if (itr->HasMember("volA") && (*itr)["volA"].IsNumber()){
        params.volA = cpjson::get_double(*itr, "volA");
    }
    if (itr->HasMember("dipm") && (*itr)["dipm"].IsNumber()){
        params.dipm = cpjson::get_double(*itr, "dipm");
    }
    if (itr->HasMember("dipnum") && (*itr)["dipnum"].IsNumber()){
        params.dipnum = cpjson::get_double(*itr, "dipnum");
    }
    molemass = cpjson::get_double(*itr, "molemass");
    aliases = cpjson::get_string_array(*itr, "aliases");
}

} /* namespace CoolProp */
