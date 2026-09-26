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
