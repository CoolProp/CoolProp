#ifndef COOLPROP_HASH_H
#define COOLPROP_HASH_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>

namespace CoolProp {

// FNV-1a 64-bit hash over a raw byte buffer.
//
// Same hash family that CoolProp already uses for stamping
// `source_eos_hash` on superancillaries (see the `Superancillary
// source_eos_hash matches current EOS at bit level` test in
// src/Tests/CoolProp-Tests.cpp).  Deterministic across compilers
// and platforms.  Used here for short cache-key prefixes computed
// over canonical-JSON option blobs.
//
// References:
//   - Fowler / Noll / Vo 1991 — http://www.isthe.com/chongo/tech/comp/fnv/
//   - constants per the canonical specification (offset basis +
//     prime are 64-bit)
inline std::uint64_t fnv1a_64(const void* data, std::size_t n) noexcept {
    constexpr std::uint64_t kFnv64OffsetBasis = 0xCBF29CE484222325ULL;
    constexpr std::uint64_t kFnv64Prime = 0x100000001B3ULL;
    std::uint64_t h = kFnv64OffsetBasis;
    const auto* p = static_cast<const std::uint8_t*>(data);
    for (std::size_t i = 0; i < n; ++i) {
        h ^= static_cast<std::uint64_t>(p[i]);
        h *= kFnv64Prime;
    }
    return h;
}

inline std::uint64_t fnv1a_64(std::string_view s) noexcept {
    return fnv1a_64(s.data(), s.size());
}

// Lowercase 16-char hex representation of a 64-bit hash.  No "0x"
// prefix.  Suitable as a filename-safe cache-key fragment.
inline std::string to_hex16(std::uint64_t h) {
    static constexpr char kHex[] = "0123456789abcdef";
    std::string out(16, '0');
    for (int i = 15; i >= 0; --i) {
        out[i] = kHex[h & 0xF];
        h >>= 4;
    }
    return out;
}

}  // namespace CoolProp

#endif  // COOLPROP_HASH_H
