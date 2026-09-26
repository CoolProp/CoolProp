"""Extract the Bollengier et al. (2019) water Gibbs-energy spline coefficients.

ONE-TIME DEVELOPER STEP.  The generated header is committed, so a normal
build never runs this and CoolProp gains no dependency from it.  Re-run it
only to regenerate or re-verify the coefficients.

Source
------
Bollengier, Brown & Shaw, "Thermodynamics of pure liquid water: Sound speed
measurements to 700 MPa down to the freezing point, and an equation of state
to 2300 MPa from 240 to 500 K", J. Chem. Phys. 151, 054501 (2019).
doi:10.1063/1.5097179

The LBF parameters are in the published Supplementary Material, file
  SM_D_data_eos_scripts.zip -> WaterEOS.mat
in the HDF5 group `G_H2O_2GPa_500K`.  (The paper's body refers to
"Supplementary Material C" for these; that does not match the published
layout -- SM_C is a prose discussion and SM_D holds the data.)

The supplementary archive is paywalled, so it is NOT vendored here.  Only
the coefficients -- which are numerical facts, not a creative work -- are
committed, together with the hashes below so the derivation is auditable by
anyone who has the supplement.

  supplement zip  sha256 955dfd39800c28209f00b8e8d6337de02d46b7a5e6145bf97b409cb3f59dcdff
  WaterEOS.mat    sha256 7f630a868db76a8e29f3ed911be93e69f8ee0500c75a04ab5a95dad479a52b27

None of the authors' MATLAB code (fnGval.m, IAPWS95.m, therm_surf.m,
DemoScript) is used or reproduced; the evaluator is CoolProp's own de Boor
implementation in include/CoolProp/spline/TensorBSpline2D.h.

Usage
-----
    python3 dev/scripts/extract_bollengier_coefficients.py /path/to/WaterEOS.mat

WaterEOS.mat is a MATLAB v7.3 file, i.e. HDF5.  scipy.io.loadmat cannot read
that format.  This script uses h5py when available and otherwise shells out to
the `h5dump` CLI, so it needs neither as a hard dependency of CoolProp itself.
"""
import hashlib
import re
import subprocess
import sys
from pathlib import Path

GROUP = "G_H2O_2GPa_500K"
EXPECTED_MAT_SHA256 = "7f630a868db76a8e29f3ed911be93e69f8ee0500c75a04ab5a95dad479a52b27"

# Expected SUPPORT of the fitted surface -- the knot span, used to verify
# the extraction rather than to describe the model's validity.  If a future
# supplement revision moves these, the mismatch should be loud.
#
# NOTE these are deliberately NOT the paper's stated validity range.  The
# paper's title and abstract give 240-500 K; the fitted knots run slightly
# wider.  The backend advertises the PAPER's range (see
# BollengierBackend::kPaperTminK); these values describe the data file.
SUPPORT_P_MIN_MPA, SUPPORT_P_MAX_MPA = 0.0, 2300.6
SUPPORT_T_MIN_K, SUPPORT_T_MAX_K = 239.0, 501.0
PAPER_ORDER = (6, 6)
PAPER_SHAPE = (80, 40)


def _h5dump_numbers(mat, path):
    """Read a numeric dataset via the h5dump CLI (no h5py needed).

    `-m %.17g` is NOT optional.  h5dump's default float format is %g, i.e.
    SIX significant figures.  Without this flag every coefficient and knot
    comes back truncated -- and truncated values still round-trip through
    `%.17g` on output, so the generated header looks like full-precision
    doubles while carrying ~1e-6 relative error.  That is worth ~0.2% in
    sound speed and ~1e-4 in density, far outside the paper's uncertainty
    and outside this project's own test tolerances.
    """
    out = subprocess.run(["h5dump", "-m", "%.17g", "-d", path, str(mat)], capture_output=True, text=True, check=True).stdout
    idx = re.compile(r"^\((\d+)(?:,(\d+))?\):\s*(.*)$")
    flt = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")
    vals = []
    for line in out.splitlines():
        m = idx.match(line.strip())
        if not m:
            continue
        body = m.group(3)
        if '"' in body:  # H5PATH / MATLAB_class attribute strings
            continue
        vals += [float(v) for v in flt.findall(body)]
    return vals


def _h5dump_ref_targets(mat, path):
    out = subprocess.run(["h5dump", "-d", path, str(mat)], capture_output=True, text=True, check=True).stdout
    return re.findall(r'DATASET "[^"]*#refs#/([A-Za-z])"', out)


def read_spline(mat):
    """Return (order_x, order_y, knots_x, knots_y, coefs_row_major)."""
    try:
        import h5py  # noqa: F401  (optional fast path)
    except ImportError:
        order = [int(v) for v in _h5dump_numbers(mat, f"/{GROUP}/order")]
        number = [int(v) for v in _h5dump_numbers(mat, f"/{GROUP}/number")]
        ka, kb = _h5dump_ref_targets(mat, f"/{GROUP}/knots")
        kx = _h5dump_numbers(mat, f"/#refs#/{ka}")
        ky = _h5dump_numbers(mat, f"/#refs#/{kb}")
        flat = _h5dump_numbers(mat, f"/{GROUP}/coefs")
    else:
        import h5py
        import numpy as np

        with h5py.File(mat, "r") as f:
            g = f[GROUP]
            order = [int(v) for v in np.array(g["order"]).ravel()]
            number = [int(v) for v in np.array(g["number"]).ravel()]
            kx = [float(v) for v in np.array(f[g["knots"][0, 0]]).ravel()]
            ky = [float(v) for v in np.array(f[g["knots"][1, 0]]).ravel()]
            flat = [float(v) for v in np.array(g["coefs"]).ravel()]

    nx, ny = number
    # Check the RAW read before reshaping.  The comprehension below always
    # produces exactly nx*ny elements regardless of how many were read, so a
    # short or long `flat` would otherwise be silently masked (or raise an
    # opaque IndexError) instead of reported.
    if len(flat) != nx * ny:
        raise ValueError(f"read {len(flat)} coefficients from {GROUP}/coefs, expected {nx} x {ny} = {nx * ny}")
    # HDF5 stores the (nx, ny) MATLAB array transposed as (ny, nx).
    coefs = [flat[b * nx + a] for a in range(nx) for b in range(ny)]
    return order[0], order[1], kx, ky, coefs


def main(argv):
    if len(argv) != 2:
        print(__doc__)
        return 2
    mat = Path(argv[1])
    digest = hashlib.sha256(mat.read_bytes()).hexdigest()
    if digest != EXPECTED_MAT_SHA256:
        print(f"ERROR: {mat} sha256 {digest}\n       expected {EXPECTED_MAT_SHA256}", file=sys.stderr)
        print("Refusing to emit coefficients from an unrecognised file.", file=sys.stderr)
        return 1

    ox, oy, kx, ky, coefs = read_spline(mat)

    # Cross-check the extraction against the paper's own stated values.  These
    # are literals above, not values read back from the file, so a mis-parse
    # or a silently different revision fails here instead of downstream.
    problems = []
    if (ox, oy) != PAPER_ORDER:
        problems.append(f"order {(ox, oy)} != paper {PAPER_ORDER}")
    nx, ny = len(kx) - ox, len(ky) - oy
    if (nx, ny) != PAPER_SHAPE:
        problems.append(f"shape {(nx, ny)} != paper {PAPER_SHAPE}")
    if len(coefs) != nx * ny:
        problems.append(f"{len(coefs)} coefficients for a {nx}x{ny} grid")
    # Compared with a tolerance, NOT exact equality.  The fitted knots carry
    # ULP-class noise -- the true upper P knot is 2300.5999999999995, not
    # 2300.6 -- so an exact comparison against the paper's rounded literals
    # rejects the correct data and accepts a truncated copy.  (These are the
    # fit's own rounded bounds, not a validity claim.)  An earlier
    # version of this script did exactly that: it passed only because the
    # values had been degraded to six significant figures, and refused to
    # run at all on full precision.  The sha256 gate above already pins the
    # file's content; these checks exist to catch a mis-parse or a different
    # revision, which a 1e-9 relative tolerance detects perfectly well.
    def _near(a, b):
        return abs(a - b) <= 1e-9 * max(1.0, abs(b))

    if not (_near(kx[ox - 1], SUPPORT_P_MIN_MPA) and _near(kx[nx], SUPPORT_P_MAX_MPA)):
        problems.append(f"P support {(kx[ox-1], kx[nx])} != expected {(SUPPORT_P_MIN_MPA, SUPPORT_P_MAX_MPA)}")
    if not (_near(ky[oy - 1], SUPPORT_T_MIN_K) and _near(ky[ny], SUPPORT_T_MAX_K)):
        problems.append(f"T support {(ky[oy-1], ky[ny])} != expected {(SUPPORT_T_MIN_K, SUPPORT_T_MAX_K)}")
    if problems:
        for p in problems:
            print("ERROR: " + p, file=sys.stderr)
        return 1

    def arr(name, vals, per_line=6):
        lines = [f"inline constexpr double {name}[] = {{"]
        for i in range(0, len(vals), per_line):
            lines.append("    " + ", ".join(f"{v:.17g}" for v in vals[i : i + per_line]) + ",")
        lines.append("};")
        return "\n".join(lines)

    print(f"""// GENERATED FILE -- do not edit by hand.
//
// Regenerate with:
//   python3 dev/scripts/extract_bollengier_coefficients.py <path>/WaterEOS.mat
//
// Gibbs-energy spline coefficients for liquid water from
//   Bollengier, Brown & Shaw, J. Chem. Phys. 151, 054501 (2019),
//   doi:10.1063/1.5097179
// published Supplementary Material, SM_D_data_eos_scripts.zip -> WaterEOS.mat,
// HDF5 group `{GROUP}`.
//
//   WaterEOS.mat sha256 {digest}
//
// The supplementary archive is paywalled and is deliberately NOT vendored.
// These coefficients are numerical facts extracted from it; none of the
// authors' MATLAB code is used or reproduced.  See the extraction script for
// the full provenance note.
//
// G is in J/kg, P in MPa, T in K.  Order is the B-spline ORDER (degree + 1),
// matching CoolProp::spline::TensorBSpline2D.
//
// Knot support: P in [{kx[ox-1]!r}, {kx[nx]!r}] MPa,
// T in [{ky[oy-1]!r}, {ky[ny]!r}] K -- the ACTUAL fitted bounds, not the
// paper's rounded statement of them.  This is the span of the FIT, not a
// validity claim: the paper states 240-500 K, which is what the backend
// advertises.  Liquid only: no vapour branch, no saturation curve.

#ifndef COOLPROP_BOLLENGIER_WATER_COEFFICIENTS_H
#define COOLPROP_BOLLENGIER_WATER_COEFFICIENTS_H

#include <cstddef>

// clang-format off
//
// Formatting is disabled for this file on purpose.  It is generated and
// committed (the source supplement is paywalled and cannot be vendored), so
// the generator's output has to stay byte-stable under the formatter --
// otherwise every regeneration produces a spurious diff, or worse, the
// committed file silently stops matching what the script emits.

namespace CoolProp {{
namespace Bollengier {{

inline constexpr std::size_t kOrderP = {ox};
inline constexpr std::size_t kOrderT = {oy};
inline constexpr std::size_t kNP = {nx};
inline constexpr std::size_t kNT = {ny};

{arr("kKnotsP", kx)}

{arr("kKnotsT", ky)}

// Row-major, shape (kNP, kNT).
{arr("kCoefs", coefs)}

}}  // namespace Bollengier
}}  // namespace CoolProp

// clang-format on

#endif  // COOLPROP_BOLLENGIER_WATER_COEFFICIENTS_H""")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
