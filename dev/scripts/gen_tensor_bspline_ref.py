"""Generate scipy reference values for the TensorBSpline2D Catch2 tests.

Emits C++ literals for src/Tests/CoolProp-Tests-TensorBSpline.cpp.  The
fixture uses a DIFFERENT order per axis on purpose so that a bug which
swaps or conflates the two axes cannot pass, and non-uniform interior
knots so the interval search is exercised.

Order convention here and in TensorBSpline2D matches MATLAB B-form and
the Bollengier et al. (2019) coefficient file: order = degree + 1.
scipy.interpolate.BSpline takes the DEGREE, hence the `- 1` below.

Usage:
    python3 dev/scripts/gen_tensor_bspline_ref.py
"""
import numpy as np
from scipy.interpolate import BSpline

# Deliberately DIFFERENT order per axis, so an axis-mixing bug cannot hide,
# and non-uniform interior knots so the interval search is exercised.
kx = np.array([0.0,0.0,0.0,0.0, 0.3, 0.7, 1.0,1.0,1.0,1.0])   # order 4 (cubic)
ky = np.array([0.0,0.0,0.0, 1.5, 2.2, 4.0,4.0,4.0])           # order 3 (quadratic)
ox, oy = 4, 3
nx, ny = len(kx) - ox, len(ky) - oy
assert (nx, ny) == (6, 5), (nx, ny)

# Deterministic, non-trivial, mixed-sign coefficients.
rng = np.random.default_rng(20260926)
C = np.round(rng.uniform(-3.0, 3.0, size=(nx, ny)), 6)

def ev(x, y, dx=0, dy=0):
    Bx = BSpline(kx, np.eye(nx), ox - 1, extrapolate=False)
    By = BSpline(ky, np.eye(ny), oy - 1, extrapolate=False)
    vx = (Bx.derivative(dx) if dx else Bx)(x)
    vy = (By.derivative(dy) if dy else By)(y)
    return float(vx @ C @ vy)

print("// knots_x")
print("    const std::vector<double> kx{%s};" % ", ".join(f"{v:.1f}" for v in kx))
print("// knots_y")
print("    const std::vector<double> ky{%s};" % ", ".join(f"{v:.1f}" for v in ky))
print("// coefs row-major (nx=%d, ny=%d)" % (nx, ny))
rows = [", ".join(f"{v:.6f}" for v in C[i]) for i in range(nx)]
print("    const std::vector<double> coefs{\n      " + ",\n      ".join(rows) + "};")

# Probe points: interior, plus x exactly ON an interior knot (0.3) and
# y exactly on an interior knot (2.2), plus both domain corners.
pts = [(0.15, 0.9), (0.3, 2.2), (0.55, 3.1), (0.92, 1.2), (0.0, 0.0), (1.0, 4.0)]
derivs = [(0,0),(1,0),(0,1),(2,0),(0,2),(1,1)]
print()
print("    // {x, y, f, f_x, f_y, f_xx, f_yy, f_xy}  (scipy.interpolate.BSpline)")
print("    const double ref[][8] = {")
for (x,y) in pts:
    vals = [ev(x,y,dx,dy) for (dx,dy) in derivs]
    print("      {%s, %s}," % (f"{x:.2f}, {y:.2f}", ", ".join(f"{v:.17g}" for v in vals)))
print("    };")

# ---------------------------------------------------------------------------
# Fixture 2: interior knot of MULTIPLICITY 2.  Exercises a zero-width span.
# (Previously generated ad hoc; emitted here so the test file's "re-run this
# script to regenerate" claim is actually true.)
# ---------------------------------------------------------------------------
kx2 = np.array([0.0,0.0,0.0,0.0, 0.5,0.5, 1.0,1.0,1.0,1.0]); ox2 = 4
ky2 = np.array([0.0,0.0,0.0, 2.0, 4.0,4.0,4.0]);             oy2 = 3
nx2, ny2 = len(kx2)-ox2, len(ky2)-oy2
C2 = np.round(np.random.default_rng(7).uniform(-2, 2, size=(nx2, ny2)), 6)

# ---------------------------------------------------------------------------
# Fixture 3: UNCLAMPED uniform knots.  The support is [knots[order-1],
# knots[n]], strictly inside the knot vector -- which is what check_in_domain
# must use.  Using knots[0]/knots[-1] instead would be wrong here but
# indistinguishable on a clamped vector.
# ---------------------------------------------------------------------------
kx3 = np.arange(0.0, 10.0); ox3 = 4          # n = 6, support [3, 6]
ky3 = np.arange(0.0, 8.0);  oy3 = 3          # n = 5, support [2, 5]
nx3, ny3 = len(kx3)-ox3, len(ky3)-oy3
C3 = np.round(np.random.default_rng(11).uniform(-2, 2, size=(nx3, ny3)), 6)

# ---------------------------------------------------------------------------
# Fixture 4: HIGH ORDER.  The falling-factorial scaling in basis_ders is
# p*(p-1)*...; at order 15 it exceeds 2^31, so computing it in `int` silently
# returns wrong derivatives.  kMaxOrder advertises 16, so this range must work.
# ---------------------------------------------------------------------------
ox4 = 15
kx4 = np.concatenate([np.zeros(ox4), np.ones(ox4)])   # clamped, n = 15
ky4 = np.array([0.0,0.0, 1.0,1.0]); oy4 = 2           # n = 2
nx4, ny4 = len(kx4)-ox4, len(ky4)-oy4
C4 = np.round(np.random.default_rng(13).uniform(-1, 1, size=(nx4, ny4)), 6)

def emit(name, kxx, kyy, oxx, oyy, CC, pts, derivs):
    nxx, nyy = len(kxx)-oxx, len(kyy)-oyy
    print(f"\n// ---- {name}: order {oxx}/{oyy}, n = {nxx}x{nyy} ----")
    print("    const std::vector<double> %s_kx{%s};" % (name, ", ".join(f"{v:g}" for v in kxx)))
    print("    const std::vector<double> %s_ky{%s};" % (name, ", ".join(f"{v:g}" for v in kyy)))
    print("    const std::vector<double> %s_c{%s};" % (name, ", ".join(f"{v:.6f}" for v in CC.ravel())))
    Bx = BSpline(kxx, np.eye(nxx), oxx-1, extrapolate=False)
    By = BSpline(kyy, np.eye(nyy), oyy-1, extrapolate=False)
    print("    // {x, y, " + ", ".join(f"d{dx}{dy}" for dx,dy in derivs) + "}")
    for (x, y) in pts:
        vals = []
        for (dx, dy) in derivs:
            vx = (Bx.derivative(dx) if dx else Bx)(x)
            vy = (By.derivative(dy) if dy else By)(y)
            vals.append(float(vx @ CC @ vy))
        print("      {%g, %g, %s}," % (x, y, ", ".join(f"{v:.17g}" for v in vals)))

emit("mult", kx2, ky2, ox2, oy2, C2, [(0.25,1.0),(0.5,2.0),(0.75,3.0),(0.5,4.0)], [(0,0),(1,0),(0,1)])
emit("unclamped", kx3, ky3, ox3, oy3, C3, [(3.0,2.0),(4.25,3.5),(5.5,4.0),(6.0,5.0)], [(0,0),(1,0),(0,1)])
emit("highorder", kx4, ky4, ox4, oy4, C4, [(0.3,0.5),(0.7,0.25)], [(0,0),(10,0),(12,0),(14,0)])

