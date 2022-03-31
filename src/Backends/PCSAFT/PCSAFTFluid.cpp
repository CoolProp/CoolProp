#include <string>
#include <vector>
#include <map>
#include <math.h>

#include "PCSAFTFluid.h"
#include "rapidjson_include.h"

namespace CoolProp {

PCSAFTFluid::PCSAFTFluid(rapidjson::Value::ValueIterator itr) {
    name = cpjson::get_string(*itr, "name");
    CAS = cpjson::get_string(*itr, "CAS");
    params.m = cpjson::get_double(*itr, "m");
    params.sigma = cpjson::get_double(*itr, "sigma");
    params.u = cpjson::get_double(*itr, "u");

    if (itr->HasMember("uAB") && (*itr)["uAB"].IsNumber()) {
        params.uAB = cpjson::get_double(*itr, "uAB");
    } else {
        params.uAB = 0.;
    }

    if (itr->HasMember("volA") && (*itr)["volA"].IsNumber()) {
        params.volA = cpjson::get_double(*itr, "volA");
    } else {
        params.volA = 0.;
    }

    if (itr->HasMember("dipm") && (*itr)["dipm"].IsNumber()) {
        params.dipm = cpjson::get_double(*itr, "dipm");
    } else {
        params.dipm = 0.;
    }

    if (itr->HasMember("dipnum") && (*itr)["dipnum"].IsNumber()) {
        params.dipnum = cpjson::get_double(*itr, "dipnum");
    } else {
        params.dipnum = 0.;
    }

    if (itr->HasMember("charge") && (*itr)["charge"].IsNumber()) {
        params.z = cpjson::get_double(*itr, "charge");
    } else {
        params.z = 0.;
    }

    molemass = cpjson::get_double(*itr, "molemass");
    aliases = cpjson::get_string_array(*itr, "aliases");
}

void PCSAFTFluid::calc_water_sigma(double t) {
    if (t > 473.16) {
        throw ValueError("The current function for sigma for water is only valid for temperatures below 473.15 K.");
    } else if (t < 273) {
        throw ValueError("The current function for sigma for water is only valid for temperatures above 273.15 K.");
    }

    params.sigma = 3.8395 + 1.2828 * exp(-0.0074944 * t) - 1.3939 * exp(-0.00056029 * t);
}

} /* namespace CoolProp */
