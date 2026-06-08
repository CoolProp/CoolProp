// Byte-equivalence gate: the embedded CBOR blob
// must decode to exactly the source dev/all_fluids.json (value-equality). Catches
// encoder/decoder drift (cbor2/nlohmann version bumps) and staleness.
#if defined(ENABLE_CATCH)
#    include <catch2/catch_all.hpp>
#    include <fstream>
#    include <sstream>
#    include <string>
#    include "CoolProp/detail/json.h"

extern "C"
{
    extern unsigned char gall_fluids_CBORData[];
    extern unsigned int gall_fluids_CBORSize;
}

TEST_CASE("embedded CBOR decodes to the source all_fluids.json", "[cbor]") {
    nlohmann::json from_blob = cpjson::from_cbor(gall_fluids_CBORData, gall_fluids_CBORSize);
    std::ifstream f(COOLPROP_ALL_FLUIDS_JSON_PATH, std::ios::binary);
    REQUIRE(f.good());
    std::ostringstream ss;
    ss << f.rdbuf();
    nlohmann::json from_src = cpjson::parse(ss.str());
    REQUIRE(from_blob == from_src);
}
#endif
