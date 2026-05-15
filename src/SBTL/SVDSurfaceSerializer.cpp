#include "CoolProp/sbtl/SVDSurfaceSerializer.h"

#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include <msgpack.hpp>

#include "CPfilepaths.h"
#include "CoolProp/region/ConstantCurve.h"
#include "CoolProp/region/CubicSplineCurve.h"
#include "CoolProp/region/PiecewiseChebyshevCurve.h"
#include "CoolProp/region/Region.h"
#include "CoolProp/svd/SVDDecomposition.h"
#include "miniz.h"

namespace CoolProp {
namespace sbtl {

namespace {

// Discriminator for the BoundaryCurve subclass stored in a packed
// region blob.  Adding a new subclass: bump revision, give it a new
// kind id, write/read its State.
enum class CurveKind : std::uint8_t
{
    CONSTANT = 0,
    CUBIC_SPLINE = 1,
    PIECEWISE_CHEBYSHEV = 2,
};

// ---------- packing helpers -------------------------------------------

template <typename Packer>
void pack_curve(Packer& pk, const region::BoundaryCurve& curve) {
    // dynamic_cast probes for each known concrete type.  Order is
    // chosen so the cheapest cast (ConstantCurve) is tried first.
    if (const auto* c = dynamic_cast<const region::ConstantCurve*>(&curve)) {
        const auto s = c->state();
        pk.pack_array(4);
        pk.pack(static_cast<std::uint8_t>(CurveKind::CONSTANT));
        pk.pack(s.a_lo);
        pk.pack(s.a_hi);
        pk.pack(s.b);
        return;
    }
    if (const auto* c = dynamic_cast<const region::CubicSplineCurve*>(&curve)) {
        const auto s = c->state();
        pk.pack_array(6);
        pk.pack(static_cast<std::uint8_t>(CurveKind::CUBIC_SPLINE));
        pk.pack(s.a);
        pk.pack(s.b);
        pk.pack(s.M);
        pk.pack(s.b_min);
        pk.pack(s.b_max);
        return;
    }
    if (const auto* c = dynamic_cast<const region::PiecewiseChebyshevCurve*>(&curve)) {
        const auto s = c->state();
        pk.pack_array(7);
        pk.pack(static_cast<std::uint8_t>(CurveKind::PIECEWISE_CHEBYSHEV));
        pk.pack(s.a_lo);
        pk.pack(s.a_hi);
        pk.pack(static_cast<std::uint8_t>(s.scale));
        pk.pack_array(s.pieces.size());
        for (const auto& p : s.pieces) {
            pk.pack_array(6);
            pk.pack(p.t_lo);
            pk.pack(p.t_hi);
            pk.pack(p.inv_half_span);
            pk.pack(p.t_mid);
            pk.pack(p.coeffs);
            pk.pack(p.deriv_coeffs);
        }
        pk.pack(s.b_min);
        pk.pack(s.b_max);
        return;
    }
    throw std::runtime_error("SVDSurfaceSerializer: unknown BoundaryCurve subclass");
}

template <typename Packer>
void pack_decomp(Packer& pk, std::size_t region_idx, ::CoolProp::parameters prop_key, const svd::SVDDecomposition& d) {
    pk.pack_array(14);
    pk.pack(static_cast<std::uint32_t>(region_idx));
    pk.pack(static_cast<std::int32_t>(prop_key));
    pk.pack(d.NX);
    pk.pack(d.NY);
    pk.pack(d.rank);
    pk.pack(static_cast<std::uint8_t>(d.out_transform));
    pk.pack(static_cast<std::uint8_t>(d.slope_source));
    pk.pack(d.x_grid);
    pk.pack(d.y_grid);
    pk.pack(d.U);
    pk.pack(d.dU_dx);
    pk.pack(d.V_S);
    pk.pack(d.dV_S_dy);
    pk.pack(d.S);
}

template <typename Packer>
void pack_region(Packer& pk, const region::Region& r) {
    const auto& axis = r.primary();
    pk.pack_array(8);
    pk.pack(static_cast<std::uint8_t>(axis.scale));
    pk.pack(axis.a_lo);
    pk.pack(axis.a_hi);
    pk.pack(axis.a_lo_t);
    pk.pack(axis.a_hi_t);
    pk.pack(axis.inv_span_t);
    pack_curve(pk, r.b_lo());
    pack_curve(pk, r.b_hi());
}

template <typename Packer>
void pack_surface(Packer& pk, const SVDSurface& s) {
    pk.pack_array(4);
    pk.pack(static_cast<std::int32_t>(s.input_pair()));

    const auto& props = s.properties();
    pk.pack_array(props.size());
    for (const auto& p : props) {
        // Property entry is [key, transform].  Transform is currently
        // baked into the SVDDecomposition itself, but the surface-level
        // record is the spec — kept here so a deserializer can verify
        // the transform per-property matches what the decomp expects.
        pk.pack_array(1);
        pk.pack(static_cast<std::int32_t>(p));
    }

    const std::size_t n_regions = s.region_count();
    pk.pack_array(n_regions);
    for (std::size_t r = 0; r < n_regions; ++r) {
        pack_region(pk, s.atlas().region(r));
    }

    // Decomps are flattened: region_idx ranges first, then properties.
    pk.pack_array(n_regions * props.size());
    for (std::size_t r = 0; r < n_regions; ++r) {
        for (const auto& p : props) {
            pack_decomp(pk, r, p, s.decomposition(r, p));
        }
    }
}

// ---------- unpacking helpers ----------------------------------------

// Convenience wrappers around msgpack::object → typed conversion.
template <typename T>
T as(const msgpack::object& o) {
    return o.as<T>();
}

void check_array(const msgpack::object& o, std::size_t expected_min, const char* what) {
    if (o.type != msgpack::type::ARRAY) {
        throw std::runtime_error(std::string("SVDSurfaceSerializer: expected array for ") + what);
    }
    if (o.via.array.size < expected_min) {
        throw std::runtime_error(std::string("SVDSurfaceSerializer: array for ") + what + " too short");
    }
}

std::unique_ptr<region::BoundaryCurve> unpack_curve(const msgpack::object& o) {
    check_array(o, 1, "curve");
    const auto kind = static_cast<CurveKind>(as<std::uint8_t>(o.via.array.ptr[0]));
    switch (kind) {
        case CurveKind::CONSTANT: {
            check_array(o, 4, "constant curve");
            const auto a_lo = as<double>(o.via.array.ptr[1]);
            const auto a_hi = as<double>(o.via.array.ptr[2]);
            const auto b = as<double>(o.via.array.ptr[3]);
            return std::make_unique<region::ConstantCurve>(a_lo, a_hi, b);
        }
        case CurveKind::CUBIC_SPLINE: {
            check_array(o, 6, "cubic spline curve");
            region::CubicSplineCurve::State s;
            s.a = as<std::vector<double>>(o.via.array.ptr[1]);
            s.b = as<std::vector<double>>(o.via.array.ptr[2]);
            s.M = as<std::vector<double>>(o.via.array.ptr[3]);
            s.b_min = as<double>(o.via.array.ptr[4]);
            s.b_max = as<double>(o.via.array.ptr[5]);
            return region::CubicSplineCurve::from_state(std::move(s));
        }
        case CurveKind::PIECEWISE_CHEBYSHEV: {
            check_array(o, 7, "piecewise chebyshev curve");
            region::PiecewiseChebyshevCurve::State s;
            s.a_lo = as<double>(o.via.array.ptr[1]);
            s.a_hi = as<double>(o.via.array.ptr[2]);
            s.scale = static_cast<region::PiecewiseChebyshevCurve::ParamScale>(as<std::uint8_t>(o.via.array.ptr[3]));
            const auto& pieces_obj = o.via.array.ptr[4];
            check_array(pieces_obj, 1, "chebyshev pieces");
            s.pieces.reserve(pieces_obj.via.array.size);
            for (std::uint32_t i = 0; i < pieces_obj.via.array.size; ++i) {
                const auto& po = pieces_obj.via.array.ptr[i];
                check_array(po, 6, "chebyshev piece");
                region::PiecewiseChebyshevCurve::PieceState ps;
                ps.t_lo = as<double>(po.via.array.ptr[0]);
                ps.t_hi = as<double>(po.via.array.ptr[1]);
                ps.inv_half_span = as<double>(po.via.array.ptr[2]);
                ps.t_mid = as<double>(po.via.array.ptr[3]);
                ps.coeffs = as<std::vector<double>>(po.via.array.ptr[4]);
                ps.deriv_coeffs = as<std::vector<double>>(po.via.array.ptr[5]);
                s.pieces.push_back(std::move(ps));
            }
            s.b_min = as<double>(o.via.array.ptr[5]);
            s.b_max = as<double>(o.via.array.ptr[6]);
            return region::PiecewiseChebyshevCurve::from_state(std::move(s));
        }
    }
    throw std::runtime_error("SVDSurfaceSerializer: unknown curve kind");
}

region::Region unpack_region(const msgpack::object& o) {
    check_array(o, 8, "region");
    const auto scale = static_cast<region::AxisScale>(as<std::uint8_t>(o.via.array.ptr[0]));
    const auto a_lo = as<double>(o.via.array.ptr[1]);
    const auto a_hi = as<double>(o.via.array.ptr[2]);
    auto axis = region::AxisTransform::make(scale, a_lo, a_hi);
    // a_lo_t/a_hi_t/inv_span_t are derived; we re-derive via make()
    // and don't read the stream values — they're stored for round-
    // trip diagnostics only.
    auto b_lo = unpack_curve(o.via.array.ptr[6]);
    auto b_hi = unpack_curve(o.via.array.ptr[7]);
    return region::Region(axis, std::move(b_lo), std::move(b_hi));
}

svd::SVDDecomposition unpack_decomp(const msgpack::object& o, std::size_t* out_region_idx, ::CoolProp::parameters* out_prop) {
    check_array(o, 14, "decomp");
    *out_region_idx = static_cast<std::size_t>(as<std::uint32_t>(o.via.array.ptr[0]));
    *out_prop = static_cast<::CoolProp::parameters>(as<std::int32_t>(o.via.array.ptr[1]));
    svd::SVDDecomposition d;
    d.NX = as<std::int32_t>(o.via.array.ptr[2]);
    d.NY = as<std::int32_t>(o.via.array.ptr[3]);
    d.rank = as<std::int32_t>(o.via.array.ptr[4]);
    d.out_transform = static_cast<svd::OutputTransform>(as<std::uint8_t>(o.via.array.ptr[5]));
    d.slope_source = static_cast<svd::SlopeSource>(as<std::uint8_t>(o.via.array.ptr[6]));
    d.x_grid = as<std::vector<double>>(o.via.array.ptr[7]);
    d.y_grid = as<std::vector<double>>(o.via.array.ptr[8]);
    d.U = as<std::vector<double>>(o.via.array.ptr[9]);
    d.dU_dx = as<std::vector<double>>(o.via.array.ptr[10]);
    d.V_S = as<std::vector<double>>(o.via.array.ptr[11]);
    d.dV_S_dy = as<std::vector<double>>(o.via.array.ptr[12]);
    d.S = as<std::vector<double>>(o.via.array.ptr[13]);
    return d;
}

SVDSurface unpack_surface(const std::string& fluid_name, const msgpack::object& o) {
    check_array(o, 4, "surface");
    const auto input_pair = static_cast<::CoolProp::input_pairs>(as<std::int32_t>(o.via.array.ptr[0]));

    const auto& props_obj = o.via.array.ptr[1];
    check_array(props_obj, 0, "properties");
    std::vector<::CoolProp::parameters> properties;
    properties.reserve(props_obj.via.array.size);
    for (std::uint32_t i = 0; i < props_obj.via.array.size; ++i) {
        const auto& pe = props_obj.via.array.ptr[i];
        check_array(pe, 1, "property entry");
        properties.push_back(static_cast<::CoolProp::parameters>(as<std::int32_t>(pe.via.array.ptr[0])));
    }

    SVDSurface surface(fluid_name, input_pair, properties);

    const auto& regions_obj = o.via.array.ptr[2];
    check_array(regions_obj, 0, "regions");
    for (std::uint32_t i = 0; i < regions_obj.via.array.size; ++i) {
        surface.add_region(unpack_region(regions_obj.via.array.ptr[i]));
    }

    const auto& decomps_obj = o.via.array.ptr[3];
    check_array(decomps_obj, 0, "decomps");
    for (std::uint32_t i = 0; i < decomps_obj.via.array.size; ++i) {
        std::size_t r_idx = 0;
        ::CoolProp::parameters prop{};
        auto decomp = unpack_decomp(decomps_obj.via.array.ptr[i], &r_idx, &prop);
        surface.add_region_property_svd(r_idx, prop, std::move(decomp));
    }

    surface.seal();
    return surface;
}

// ---------- compression / decompression -------------------------------

std::vector<char> zlib_compress(const msgpack::sbuffer& sbuf) {
    // Match TabularBackends.cpp's allocation strategy: start with the
    // raw size and let zlib write into that buffer.  miniz returns
    // Z_OK when it fits.
    std::vector<char> out(sbuf.size() + (sbuf.size() / 1000) + 128);
    mz_ulong out_size = static_cast<mz_ulong>(out.size());
    const int code = compress((unsigned char*)out.data(), &out_size, (const unsigned char*)sbuf.data(), static_cast<mz_ulong>(sbuf.size()));
    if (code != Z_OK) {
        throw std::runtime_error(std::string("SVDSurfaceSerializer: zlib compress failed (code ") + std::to_string(code) + ")");
    }
    out.resize(out_size);
    return out;
}

std::vector<char> zlib_uncompress(const std::vector<char>& compressed) {
    // Mirrors TabularBackends.cpp:52-69: start at 5x the compressed
    // size, double on Z_BUF_ERROR until it fits or we exhaust patience.
    std::vector<char> out(compressed.size() * 5);
    mz_ulong out_size = static_cast<mz_ulong>(out.size());
    mz_ulong in_size = static_cast<mz_ulong>(compressed.size());
    int code = 0;
    int retries = 0;
    do {
        code = uncompress((unsigned char*)out.data(), &out_size, (const unsigned char*)compressed.data(), in_size);
        if (code == Z_BUF_ERROR) {
            if (++retries > 8) {
                throw std::runtime_error("SVDSurfaceSerializer: zlib uncompress would not fit after 8 retries");
            }
            out.resize(out.size() * 2);
            out_size = static_cast<mz_ulong>(out.size());
        } else if (code != Z_OK) {
            throw std::runtime_error(std::string("SVDSurfaceSerializer: zlib uncompress failed (code ") + std::to_string(code) + ")");
        }
    } while (code != Z_OK);
    out.resize(out_size);
    return out;
}

}  // namespace

// ---------- public API ------------------------------------------------

std::vector<char> SVDSurfaceSerializer::save(const SVDSurface& surface) {
    if (!surface.sealed()) {
        throw std::logic_error("SVDSurfaceSerializer::save: surface must be sealed");
    }
    msgpack::sbuffer sbuf;
    msgpack::packer<msgpack::sbuffer> pk(sbuf);
    pk.pack_array(4);
    pk.pack(std::string("SVDS"));
    pk.pack(static_cast<std::int32_t>(kRevision));
    pk.pack(surface.fluid_name());
    // Single-surface file for now — but the format puts surfaces in
    // an array so a single fluid's PH + PT surfaces could share one
    // file in the future without a revision bump.
    pk.pack_array(1);
    pack_surface(pk, surface);
    return zlib_compress(sbuf);
}

SVDSurface SVDSurfaceSerializer::load(const std::vector<char>& compressed) {
    const auto plain = zlib_uncompress(compressed);

    msgpack::object_handle oh;
    msgpack::unpack(oh, plain.data(), plain.size());
    const msgpack::object& root = oh.get();

    check_array(root, 4, "root");
    const auto magic = as<std::string>(root.via.array.ptr[0]);
    if (magic != "SVDS") {
        throw std::runtime_error("SVDSurfaceSerializer::load: bad magic (expected 'SVDS', got '" + magic + "')");
    }
    const auto revision = as<std::int32_t>(root.via.array.ptr[1]);
    if (revision != kRevision) {
        throw std::runtime_error("SVDSurfaceSerializer::load: revision mismatch (expected " + std::to_string(kRevision) + ", got "
                                 + std::to_string(revision) + ")");
    }
    const auto fluid_name = as<std::string>(root.via.array.ptr[2]);

    const auto& surfaces_obj = root.via.array.ptr[3];
    check_array(surfaces_obj, 1, "surfaces");
    if (surfaces_obj.via.array.size != 1) {
        throw std::runtime_error("SVDSurfaceSerializer::load: expected exactly 1 surface (got " + std::to_string(surfaces_obj.via.array.size) + ")");
    }
    return unpack_surface(fluid_name, surfaces_obj.via.array.ptr[0]);
}

void SVDSurfaceSerializer::save_to_file(const SVDSurface& surface, const std::string& path) {
    const auto compressed = save(surface);
    std::ofstream out(path, std::ios::binary);
    if (!out) {
        throw std::runtime_error("SVDSurfaceSerializer::save_to_file: cannot open " + path);
    }
    out.write(compressed.data(), static_cast<std::streamsize>(compressed.size()));
    if (!out) {
        throw std::runtime_error("SVDSurfaceSerializer::save_to_file: write failed for " + path);
    }
}

SVDSurface SVDSurfaceSerializer::load_from_file(const std::string& path) {
    std::ifstream in(path, std::ios::binary | std::ios::ate);
    if (!in) {
        throw std::runtime_error("SVDSurfaceSerializer::load_from_file: cannot open " + path);
    }
    const auto size = in.tellg();
    in.seekg(0, std::ios::beg);
    std::vector<char> buf(static_cast<std::size_t>(size));
    if (!in.read(buf.data(), size)) {
        throw std::runtime_error("SVDSurfaceSerializer::load_from_file: read failed for " + path);
    }
    return load(buf);
}

std::string SVDSurfaceSerializer::default_cache_dir() {
    const std::string dir = ::get_home_dir() + "/.CoolProp/SVDTables";
    std::error_code ec;
    std::filesystem::create_directories(dir, ec);
    return dir + "/";
}

std::string SVDSurfaceSerializer::default_cache_path(const std::string& fluid_name, ::CoolProp::input_pairs input_pair) {
    return default_cache_dir() + fluid_name + "." + std::to_string(static_cast<int>(input_pair)) + ".svd.bin.z";
}

}  // namespace sbtl
}  // namespace CoolProp
