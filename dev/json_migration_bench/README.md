# JSON migration benchmarks & interop proof

Throwaway, reproducible programs used as an **up-front de-risking gate** for the
RapidJSON → nlohmann/json migration (see
`docs/superpowers/specs/2026-06-04-rapidjson-to-nlohmann-migration-design.md`
§6). They are standalone — not part of the CMake build — and read the real
`dev/all_fluids.json`. Run from the repo root after a `build_catch` configure
(so the CPM-fetched headers exist under `build_catch/_deps/`).

## Programs

| File | Measures |
|------|----------|
| `parse_bench.cpp` | Full-DOM parse time of `all_fluids.json`: RapidJSON vs nlohmann JSON vs nlohmann MessagePack/CBOR. |
| `zlib_uncompress_bench.c` | miniz `uncompress()` time for the embedded zlib blob (the step at `FluidLibrary.cpp:49`). |
| `cbor_roundtrip.cpp` | Interop proof: Python `cbor2` → nlohmann `from_cbor` reproduces the source JSON exactly (deep value-equality + full-precision doubles). The logic the Phase-1 CI `[cbor]` byte-equivalence test will assert. |

## Build & run

```bash
RJ=build_catch/_deps/rapidjson-src/include
NL=build_catch/_deps/nlohmann_json-src/include

# Parse benchmark (rapidjson vs nlohmann json/msgpack/cbor)
c++ -std=c++17 -O2 -DNDEBUG -I "$RJ" -I "$NL" dev/json_migration_bench/parse_bench.cpp -o /tmp/parse_bench
/tmp/parse_bench dev/all_fluids.json

# zlib uncompress benchmark (miniz)
cc -O2 -DNDEBUG -I externals/miniz-3.1.1 \
   dev/json_migration_bench/zlib_uncompress_bench.c externals/miniz-3.1.1/miniz.c -o /tmp/zlib_bench
/tmp/zlib_bench dev/all_fluids.json

# CBOR interop proof: encode with Python cbor2, decode + compare with nlohmann
uvx --from cbor2 python -c "import json,cbor2; open('/tmp/all_fluids.cbor','wb').write(cbor2.dumps(json.load(open('dev/all_fluids.json'))))"
c++ -std=c++17 -O2 -DNDEBUG -I "$NL" dev/json_migration_bench/cbor_roundtrip.cpp -o /tmp/cbor_roundtrip
/tmp/cbor_roundtrip dev/all_fluids.json /tmp/all_fluids.cbor
```

## Recorded results (9.71 MB all_fluids.json, -O2, single dev machine, 2026-06-04)

```
rapidjson  JSON      median  12.8 ms     <- current parse
nlohmann   JSON      median  82.4 ms     <- ~6.4x slower (naive swap, +70 ms)
nlohmann   msgpack   median  36.6 ms
nlohmann   cbor      median  35.1 ms
miniz uncompress     mean    24.6 ms     <- dominant cost in today's load path

binary blob: ~4.3 MB (cbor/msgpack) vs 9.71 MB json text

full first-load (decompress + parse):
  current  (zlib + rapidjson json)      ~37 ms
  naive    (zlib + nlohmann json)       ~107 ms
  cbor     (no zlib + nlohmann cbor)    ~35 ms   <- parity with today, no zlib

cbor interop (python cbor2 -> nlohmann from_cbor): ROUND-TRIP PASS
  136 fluids, deep value-equal YES, canonical dumps byte-identical,
  sample double -246119.46072729214 bit-identical YES
```

**Conclusion:** embed the fluid data as **uncompressed CBOR** in Phase 1 — it
matches today's first-load (~35 ms), removes zlib/miniz from the fluid-load
path, and round-trips losslessly.
