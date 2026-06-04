// Interop proof: Python cbor2 -> nlohmann from_cbor must reproduce the source
// JSON exactly, including full-precision doubles. This is the logic a CI
// byte-equivalence test would assert.
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>

using nlohmann::json;

static void must_open(const std::ifstream& f, const char* p) {
    if (!f) {
        std::fprintf(stderr, "cannot open %s\n", p);
        std::exit(1);
    }
}
static std::string read_text(const char* p) {
    std::ifstream f(p, std::ios::binary);
    must_open(f, p);
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}
static std::vector<std::uint8_t> read_bytes(const char* p) {
    std::ifstream f(p, std::ios::binary);
    must_open(f, p);
    return std::vector<std::uint8_t>((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}

int main(int argc, char** argv) {
    const char* json_path = (argc > 1) ? argv[1] : "dev/all_fluids.json";
    const char* cbor_path = (argc > 2) ? argv[2] : "/tmp/all_fluids.cbor";

    json ref = json::parse(read_text(json_path));        // source of truth
    json from = json::from_cbor(read_bytes(cbor_path));  // python cbor2 -> nlohmann

    // 1. Deep value equality — fails if ANY number/string/structure differs,
    //    including a single-ULP double mismatch.
    bool values_equal = (ref == from);

    // 2. Byte-equality of the canonical JSON serialization (the form a CI test
    //    would compare). dump() uses shortest round-trippable float repr, so
    //    equal doubles produce identical bytes.
    std::string ref_dump = ref.dump();
    std::string from_dump = from.dump();
    bool bytes_equal = (ref_dump == from_dump);

    std::printf("fluids in array      : %zu\n", ref.is_array() ? ref.size() : 0);
    std::printf("deep value-equal     : %s\n", values_equal ? "YES" : "NO");
    std::printf("canonical-dump equal : %s  (%zu vs %zu bytes)\n", bytes_equal ? "YES" : "NO", ref_dump.size(), from_dump.size());

    // 3. Full-precision double spot-check on a known deep path.
    try {
        double a = ref.at(0).at("ANCILLARIES").at("hL").at("A").at(0).get<double>();
        double b = from.at(0).at("ANCILLARIES").at("hL").at("A").at(0).get<double>();
        std::printf("sample double  json  : %.17g\n", a);
        std::printf("sample double  cbor  : %.17g\n", b);
        std::printf("bit-identical        : %s\n", (a == b) ? "YES" : "NO");
    } catch (const std::exception& e) {
        std::printf("(sample-path probe skipped: %s)\n", e.what());
    }

    bool ok = values_equal && bytes_equal;
    std::printf("\nROUND-TRIP %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
