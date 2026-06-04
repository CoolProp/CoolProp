// Standalone parse benchmark: RapidJSON vs nlohmann/json on the real
// CoolProp all_fluids.json. De-risking gate for the migration — measures
// full-DOM parse time of the embedded fluid blob with each library.
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

#include <rapidjson/document.h>
#include <nlohmann/json.hpp>

static std::string read_file(const char* path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) {
        std::fprintf(stderr, "cannot open %s\n", path);
        std::exit(1);
    }
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

template <typename F>
static void run(const char* name, int iters, F&& fn) {
    std::vector<double> ms;
    ms.reserve(iters);
    for (int i = 0; i < iters; ++i) {
        auto t0 = std::chrono::steady_clock::now();
        fn();
        auto t1 = std::chrono::steady_clock::now();
        ms.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());
    }
    std::sort(ms.begin(), ms.end());
    double mn = ms.front();
    double md = ms[ms.size() / 2];
    double mx = ms.back();
    std::printf("%-14s  min %7.1f ms   median %7.1f ms   max %7.1f ms   (%d runs)\n", name, mn, md, mx, iters);
}

int main(int argc, char** argv) {
    const char* path = (argc > 1) ? argv[1] : "dev/all_fluids.json";
    std::string s = read_file(path);
    std::printf("file: %s   size: %.2f MB\n\n", path, s.size() / (1024.0 * 1024.0));

    const int iters = 15;

    // RapidJSON full DOM parse (matches the loaders' d.Parse<0>(...) usage).
    run("rapidjson", iters, [&] {
        rapidjson::Document d;
        d.Parse<0>(s.c_str());
        if (d.HasParseError()) {
            std::fprintf(stderr, "rapidjson parse error\n");
            std::abort();
        }
    });

    // nlohmann full DOM parse.
    run("nlohmann-json", iters, [&] {
        auto j = nlohmann::json::parse(s);
        if (!j.is_array() && !j.is_object()) {
            std::fprintf(stderr, "nlohmann parse error\n");
            std::abort();
        }
    });

    // Mitigation: pre-parsed binary blobs (nlohmann reads these natively and
    // much faster than text JSON). Generate once from the parsed doc, then
    // time deserialization — this is what a CBOR/MessagePack-embedded build
    // would pay at first load.
    nlohmann::json doc = nlohmann::json::parse(s);
    std::vector<std::uint8_t> mp = nlohmann::json::to_msgpack(doc);
    std::vector<std::uint8_t> cb = nlohmann::json::to_cbor(doc);
    std::printf("\nbinary blob sizes:  msgpack %.2f MB   cbor %.2f MB   (json text %.2f MB)\n\n", mp.size() / (1024.0 * 1024.0),
                cb.size() / (1024.0 * 1024.0), s.size() / (1024.0 * 1024.0));

    run("nlohmann-mpack", iters, [&] {
        auto j = nlohmann::json::from_msgpack(mp);
        if (!j.is_array() && !j.is_object()) {
            std::fprintf(stderr, "msgpack error\n");
            std::abort();
        }
    });
    run("nlohmann-cbor", iters, [&] {
        auto j = nlohmann::json::from_cbor(cb);
        if (!j.is_array() && !j.is_object()) {
            std::fprintf(stderr, "cbor error\n");
            std::abort();
        }
    });

    return 0;
}
