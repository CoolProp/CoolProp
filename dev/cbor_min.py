"""Minimal CBOR (RFC 8949) codec for the JSON data model — standard library only.

Vendored to keep CoolProp's build path free of a third-party Python dependency:
`dev/generate_headers.py` runs during every C++ build, in ~12 heterogeneous CI
environments (manylinux, emscripten, MSVC, conda, …), and getting a compiled
package (cbor2) installed into the exact interpreter CMake invokes proved
fragile. CBOR for the JSON value model is small and fully specified, so we
encode it directly here.

Supports exactly the JSON types: None, bool, int, float, str, list, dict (str
keys). Floats are always written as 64-bit doubles to preserve full precision.
JSON has no NaN/Infinity, so non-finite floats are rejected loudly rather than
encoded (cbor2 would shorten them to half-floats, which would also break the
byte-identity gate). The authoritative correctness gate is the C++ ``[cbor]``
byte-equivalence test (nlohmann ``from_cbor`` of the embedded blob == the source
``all_fluids.json``); ``loads`` here exists only for the generation-time
round-trip self-check.
"""
import math
import struct


def _encode_head(major, n, out):
    """Write a CBOR head: 3-bit major type + definite length/value n."""
    mt = major << 5
    if n < 24:
        out.append(mt | n)
    elif n < 0x100:
        out.append(mt | 24)
        out.append(n)
    elif n < 0x10000:
        out.append(mt | 25)
        out += struct.pack(">H", n)
    elif n < 0x100000000:
        out.append(mt | 26)
        out += struct.pack(">I", n)
    elif n < 0x10000000000000000:
        out.append(mt | 27)
        out += struct.pack(">Q", n)
    else:
        raise ValueError("integer out of range for CBOR encoding: %d" % n)


def _encode(obj, out):
    # Order matters: bool is a subclass of int, so test it first.
    if obj is None:
        out.append(0xF6)
    elif obj is True:
        out.append(0xF5)
    elif obj is False:
        out.append(0xF4)
    elif isinstance(obj, int):
        if obj >= 0:
            _encode_head(0, obj, out)
        else:
            _encode_head(1, -1 - obj, out)
    elif isinstance(obj, float):
        if not math.isfinite(obj):
            raise ValueError("non-finite float is not representable in JSON/CBOR: %r" % obj)
        out.append(0xFB)
        out += struct.pack(">d", obj)
    elif isinstance(obj, str):
        b = obj.encode("utf-8")
        _encode_head(3, len(b), out)
        out += b
    elif isinstance(obj, (bytes, bytearray)):
        _encode_head(2, len(obj), out)
        out += bytes(obj)
    elif isinstance(obj, (list, tuple)):
        _encode_head(4, len(obj), out)
        for item in obj:
            _encode(item, out)
    elif isinstance(obj, dict):
        _encode_head(5, len(obj), out)
        for k, v in obj.items():
            if not isinstance(k, str):
                raise ValueError("only string keys are supported, got %r" % type(k))
            kb = k.encode("utf-8")
            _encode_head(3, len(kb), out)
            out += kb
            _encode(v, out)
    else:
        raise ValueError("unsupported type for CBOR encoding: %r" % type(obj))


def dumps(obj):
    """Encode a JSON-model object to CBOR bytes."""
    out = bytearray()
    _encode(obj, out)
    return bytes(out)


def _decode(data, i):
    b = data[i]
    i += 1
    major = b >> 5
    ai = b & 0x1F
    if major == 7:  # simple values / floats — ai is the value, not a length
        if ai == 20:
            return False, i
        if ai == 21:
            return True, i
        if ai == 22:
            return None, i
        if ai == 26:
            return struct.unpack_from(">f", data, i)[0], i + 4
        if ai == 27:
            return struct.unpack_from(">d", data, i)[0], i + 8
        raise ValueError("unsupported simple/float value %d" % ai)
    if ai < 24:
        n = ai
    elif ai == 24:
        n = data[i]
        i += 1
    elif ai == 25:
        n = struct.unpack_from(">H", data, i)[0]
        i += 2
    elif ai == 26:
        n = struct.unpack_from(">I", data, i)[0]
        i += 4
    elif ai == 27:
        n = struct.unpack_from(">Q", data, i)[0]
        i += 8
    else:
        raise ValueError("unsupported additional info %d" % ai)
    if major == 0:
        return n, i
    if major == 1:
        return -1 - n, i
    if major == 2:
        return bytes(data[i:i + n]), i + n
    if major == 3:
        return data[i:i + n].decode("utf-8"), i + n
    if major == 4:
        arr = []
        for _ in range(n):
            v, i = _decode(data, i)
            arr.append(v)
        return arr, i
    if major == 5:
        d = {}
        for _ in range(n):
            k, i = _decode(data, i)
            v, i = _decode(data, i)
            d[k] = v
        return d, i
    raise ValueError("unsupported major type %d" % major)


def loads(data):
    """Decode CBOR bytes back to a JSON-model object (used for the self-check)."""
    obj, _ = _decode(bytes(data), 0)
    return obj
