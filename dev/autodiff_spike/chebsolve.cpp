// Experiment 5: speed of a Chebyshev all-roots density solver for PC-SAFT (hard chain +
// dispersion) built from *universal* coefficient tables.
//
// Rootfinding function on a fixed partition of eta in [0, 0.74] (see doc/):
//     G(eta) = eta * Q(eta) - t * P(eta)^2,   Q = P^2 Z,   t = p q(T,x) / (R T)
// Q is linear in ~30 (T,x)-dependent scalar weights times fluid-independent eta-functions,
// plus, per component, a one-parameter (eta, a) family for the chain term.  Three tiers:
//   offline      : Chebyshev coefficients of every universal function on every piece
//   per (T, x)   : scalar weights, then G's coefficients as a linear combination
//   per p        : subtract t * P^2, exclude root-free pieces, find roots on the rest
//
// Build: c++ -std=c++20 -O3 -march=native -I<eigen> chebsolve.cpp -o chebsolve
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <functional>
#include <vector>

namespace {
constexpr double PI = 3.14159265358979323846, NA = 6.02214076e23, R = 8.31446261815324;
constexpr int NQ = 12;      // degree of Q per piece
constexpr int NG = NQ + 1;  // degree of G (eta * Q)
constexpr int NAD = 10;     // degree in a for the chain table
constexpr double A_LO = 0.75, A_HI = 2.5;
constexpr std::array<double, 9> EDGES = {0.0, 0.140867, 0.281733, 0.382913, 0.484093, 0.556767, 0.629441, 0.681641, 0.74};
constexpr int NP = EDGES.size() - 1;
using Vec = std::array<double, NG + 1>;

constexpr double GA[3][7] = {{0.9105631445, 0.6361281449, 2.6861347891, -26.547362491, 97.759208784, -159.59154087, 91.297774084},
                             {-0.3084016918, 0.1860531159, -2.5030047259, 21.419793629, -65.255885330, 83.318680481, -33.746922930},
                             {-0.0906148351, 0.4527842806, 0.5962700728, -1.7241829131, -4.1302112531, 13.776631870, -8.6728470368}};
constexpr double GB[3][7] = {{0.7240946941, 2.2382791861, -4.0025849485, -21.003576815, 26.855641363, 206.55133841, -355.60235612},
                             {-0.5755498075, 0.6995095521, 3.8925673390, -17.215471648, 192.67226447, -161.82646165, -165.20769346},
                             {0.0976883116, -0.2557574982, -9.1558561530, 20.642075974, -38.804430052, 93.626774077, -29.666905585}};

// ------------------------------------------------------------------ universal eta-functions
// value and derivative of a few polynomials via a tiny dual type (offline only)
struct D
{
    double v, d;
};
D operator+(D a, D b) {
    return {a.v + b.v, a.d + b.d};
}
D operator-(D a, D b) {
    return {a.v - b.v, a.d - b.d};
}
D operator*(D a, D b) {
    return {a.v * b.v, a.d * b.v + a.v * b.d};
}
D operator*(double s, D a) {
    return {s * a.v, s * a.d};
}
D cst(double c) {
    return {c, 0};
}

D Kf(D e) {
    const D o = cst(1) - e, t = cst(2) - e, o2 = o * o;
    return o2 * o2 * t * t;
}
D Mf(D e) {
    const D t = cst(2) - e;
    return (8.0 * e - 2.0 * e * e) * t * t;
}
D Lf(D e) {
    const D o = cst(1) - e, e2 = e * e;
    return (20.0 * e - 27.0 * e2 + 12.0 * e2 * e - 2.0 * e2 * e2) * o * o;
}
D P0f(D e) {
    return Kf(e) + Lf(e);
}
D P1f(D e) {
    return Mf(e) - Lf(e);
}
D Isf(const double (&c)[3][7], int s, D e) {
    D r = cst(0), ek = cst(1);
    for (int k = 0; k < 7; ++k) {
        r = r + c[s][k] * ek;
        ek = ek * e;
    }
    return r;
}
double Pr(int r, double e) {
    return (r == 0 ? P0f(cst(e)) : P1f(cst(e))).v;
}
double Pi_r(int r, double e) {  // P^2 = Pi0 + mbar Pi1 + mbar^2 Pi2
    const double p0 = Pr(0, e), p1 = Pr(1, e);
    return r == 0 ? p0 * p0 : r == 1 ? 2 * p0 * p1 : p1 * p1;
}

// The 27 universal Q-basis functions, in a fixed order that the weights below mirror.
std::vector<std::function<double(double)>> make_Qbasis() {
    std::vector<std::function<double(double)>> B;
    for (int r = 0; r < 3; ++r)
        B.push_back([r](double e) { return Pi_r(r, e); });  // 1
    for (int r = 0; r < 3; ++r)
        B.push_back([r](double e) { return Pi_r(r, e) * e / ((1 - e) * (1 - e)); });  // h1
    for (int r = 0; r < 3; ++r)
        B.push_back([r](double e) { return Pi_r(r, e) * e * (1 + e) / ((1 - e) * (1 - e) * (1 - e)); });  // h2
    for (int r = 0; r < 3; ++r)
        B.push_back([r](double e) { return Pi_r(r, e) * e / (1 - e); });  // h3
    for (int s = 0; s < 3; ++s)
        for (int r = 0; r < 3; ++r)
            B.push_back([r, s](double e) {  // eta d(eta I1_s)/d eta
                const D E{e, 1.0}, f = E * Isf(GA, s, E);
                return Pi_r(r, e) * e * f.d;
            });
    for (int s = 0; s < 3; ++s)
        for (int r = 0; r < 2; ++r)
            B.push_back([r, s](double e) {  // eta [P_r (eta K I2_s)' - eta K I2_s P_r']
                const D E{e, 1.0}, f = E * Kf(E) * Isf(GB, s, E), P = r == 0 ? P0f(E) : P1f(E);
                return e * (P.v * f.d - f.v * P.d);
            });
    return B;
}
// chain family: Pi_r * eta * chi(eta; a),  chi = -alpha/(1-alpha eta) + beta/(1+beta eta)
double chain_fn(int r, double e, double a) {
    const double al = (3 - a) / 3, be = (2 * a - 3) / 3;
    return Pi_r(r, e) * e * (-al / (1 - al * e) + be / (1 + be * e));
}

// ------------------------------------------------------------------ Chebyshev utilities
template <int N, class F>
std::array<double, N + 1> chebfit(F f, double lo, double hi) {  // interpolation at first-kind nodes
    std::array<double, N + 1> fv{}, c{};
    for (int j = 0; j <= N; ++j)
        fv[j] = f(lo + (hi - lo) * (std::cos(PI * (j + 0.5) / (N + 1)) + 1) / 2);
    for (int k = 0; k <= N; ++k) {
        double s = 0;
        for (int j = 0; j <= N; ++j)
            s += fv[j] * std::cos(PI * k * (j + 0.5) / (N + 1));
        c[k] = s * (k == 0 ? 1.0 : 2.0) / (N + 1);
    }
    return c;
}
// eta * Q exactly: eta = m + h u, u T_k = (T_{k+1} + T_{|k-1|}) / 2
Vec times_eta(const std::array<double, NQ + 1>& c, double lo, double hi) {
    const double m = 0.5 * (lo + hi), h = 0.5 * (hi - lo);
    Vec r{};
    for (int k = 0; k <= NQ; ++k) {
        r[k] += m * c[k];
        if (k == 0)
            r[1] += h * c[0];
        else {
            r[k + 1] += 0.5 * h * c[k];
            r[k - 1] += 0.5 * h * c[k];
        }
    }
    return r;
}
double clenshaw(const Vec& c, double u, int n = NG) {
    double b1 = 0, b2 = 0;
    for (int k = n; k >= 1; --k) {
        const double b0 = c[k] + 2 * u * b1 - b2;
        b2 = b1;
        b1 = b0;
    }
    return c[0] + u * b1 - b2;
}

// ------------------------------------------------------------------ offline tables
struct Tables
{
    std::vector<std::array<Vec, NP>> U;                          // eta * Q-basis, [j][piece]
    std::array<std::array<Vec, NP>, 3> P2;                       // Pi_r (no eta factor), [r][piece]
    std::array<std::array<std::array<Vec, NAD + 1>, NP>, 3> Ch;  // eta * chain table, [r][piece][ka]
};
Tables build_tables() {
    Tables T;
    const auto B = make_Qbasis();
    T.U.resize(B.size());
    for (std::size_t j = 0; j < B.size(); ++j)
        for (int p = 0; p < NP; ++p)
            T.U[j][p] = times_eta(chebfit<NQ>(B[j], EDGES[p], EDGES[p + 1]), EDGES[p], EDGES[p + 1]);
    for (int r = 0; r < 3; ++r)
        for (int p = 0; p < NP; ++p) {
            const auto c = chebfit<NQ>([r](double e) { return Pi_r(r, e); }, EDGES[p], EDGES[p + 1]);  // degree 12: exact
            T.P2[r][p] = {};
            for (int k = 0; k <= NQ; ++k)
                T.P2[r][p][k] = c[k];
        }
    // 2-D table: Chebyshev in eta (per piece) x Chebyshev in a on [A_LO, A_HI]
    for (int r = 0; r < 3; ++r)
        for (int p = 0; p < NP; ++p) {
            std::array<std::array<double, NQ + 1>, NAD + 1> ce;  // eta-coeffs at each a-node
            for (int j = 0; j <= NAD; ++j) {
                const double a = A_LO + (A_HI - A_LO) * (std::cos(PI * (j + 0.5) / (NAD + 1)) + 1) / 2;
                ce[j] = chebfit<NQ>([r, a](double e) { return chain_fn(r, e, a); }, EDGES[p], EDGES[p + 1]);
            }
            for (int ka = 0; ka <= NAD; ++ka) {
                std::array<double, NQ + 1> c{};
                for (int j = 0; j <= NAD; ++j) {
                    const double w = std::cos(PI * ka * (j + 0.5) / (NAD + 1)) * (ka == 0 ? 1.0 : 2.0) / (NAD + 1);
                    for (int k = 0; k <= NQ; ++k)
                        c[k] += w * ce[j][k];
                }
                T.Ch[r][p][ka] = times_eta(c, EDGES[p], EDGES[p + 1]);
            }
        }
    return T;
}

// ------------------------------------------------------------------ per (T, x)
constexpr int NCMAX = 3;
struct Fluid
{
    int N;
    double m[NCMAX], sig[NCMAX], eps[NCMAX];
};
struct Scalars
{
    double mbar, q, A1, A2, A3, B1, B2, c1, c2, a[NCMAX], w[NCMAX];
};
Scalars scalars(const Fluid& F, double T, const double* x) {
    Scalars S{};
    double d[NCMAX], s[4] = {0, 0, 0, 0}, E1 = 0, E2 = 0;
    for (int i = 0; i < F.N; ++i) {
        d[i] = F.sig[i] * (1 - 0.12 * std::exp(-3 * F.eps[i] / T));
        const double xm = x[i] * F.m[i];
        s[0] += xm;
        s[1] += xm * d[i];
        s[2] += xm * d[i] * d[i];
        s[3] += xm * d[i] * d[i] * d[i];
        S.mbar += xm;
        for (int j = 0; j < F.N; ++j) {
            const double sij = 0.5 * (F.sig[i] + F.sig[j]), eij = std::sqrt(F.eps[i] * F.eps[j]);
            const double pre = x[i] * x[j] * F.m[i] * F.m[j] * sij * sij * sij;
            E1 += pre * eij;
            E2 += pre * eij * eij;
        }
    }
    const double r0 = s[0] / s[3], r1 = s[1] / s[3], r2 = s[2] / s[3];
    S.q = PI / 6 * NA * 1e-30 * s[3];
    S.A1 = 3 * r1 * r2 / r0;
    S.A2 = r2 * r2 * r2 / r0;
    S.A3 = S.A2 - 1;
    S.B1 = -12 * E1 / (T * s[3]);
    S.B2 = -6 * S.mbar * E2 / (T * T * s[3]);
    S.c1 = (S.mbar - 1) / S.mbar;
    S.c2 = S.c1 * (S.mbar - 2) / S.mbar;
    for (int i = 0; i < F.N; ++i) {
        S.a[i] = 1.5 * d[i] * r2;
        S.w[i] = x[i] * (F.m[i] - 1);
    }
    return S;
}

struct Assembled
{
    std::array<Vec, NP> G;   // eta Q
    std::array<Vec, NP> P2;  // P^2
    double q;
    bool ok;
};
// Tier 2: weights + linear combination of the universal tables
Assembled assemble(const Tables& Tb, const Fluid& F, double T, const double* x) {
    const Scalars S = scalars(F, T, x);
    Assembled A{};
    A.q = S.q;
    A.ok = true;
    const double mr[3] = {1, S.mbar, S.mbar * S.mbar};
    double sw = 0;
    for (int i = 0; i < F.N; ++i)
        sw += S.w[i];
    double W[27];
    int j = 0;
    for (int r = 0; r < 3; ++r)
        W[j++] = mr[r];
    for (int r = 0; r < 3; ++r)
        W[j++] = mr[r] * S.mbar * S.A1;
    for (int r = 0; r < 3; ++r)
        W[j++] = mr[r] * S.mbar * S.A2;
    for (int r = 0; r < 3; ++r)
        W[j++] = mr[r] * (-S.mbar * S.A3 - 3 * sw);
    const double cs[3] = {1, S.c1, S.c2};
    for (int s = 0; s < 3; ++s)
        for (int r = 0; r < 3; ++r)
            W[j++] = mr[r] * S.B1 * cs[s];
    for (int s = 0; s < 3; ++s)
        for (int r = 0; r < 2; ++r)
            W[j++] = mr[r] * S.B2 * cs[s];
    for (int p = 0; p < NP; ++p) {
        Vec g{}, pp{};
        for (int jj = 0; jj < 27; ++jj)
            for (int k = 0; k <= NG; ++k)
                g[k] += W[jj] * Tb.U[jj][p][k];
        for (int r = 0; r < 3; ++r)
            for (int k = 0; k <= NG; ++k)
                pp[k] += mr[r] * Tb.P2[r][p][k];
        A.G[p] = g;
        A.P2[p] = pp;
    }
    // chain term, per component: contract the a-direction, weight -w_i mbar^r
    for (int i = 0; i < F.N; ++i) {
        if (F.m[i] == 1.0) continue;  // w_i = 0
        if (S.a[i] < A_LO || S.a[i] > A_HI) {
            A.ok = false;  // outside the table rectangle: reject, caller falls back
            return A;
        }
        const double ua = 2 * (S.a[i] - A_LO) / (A_HI - A_LO) - 1;
        double Ta[NAD + 1];
        Ta[0] = 1;
        Ta[1] = ua;
        for (int k = 2; k <= NAD; ++k)
            Ta[k] = 2 * ua * Ta[k - 1] - Ta[k - 2];
        for (int r = 0; r < 3; ++r) {
            const double wr = -S.w[i] * mr[r];
            for (int p = 0; p < NP; ++p)
                for (int ka = 0; ka <= NAD; ++ka) {
                    const double f = wr * Ta[ka];
                    for (int k = 0; k <= NG; ++k)
                        A.G[p][k] += f * Tb.Ch[r][p][ka][k];
                }
        }
    }
    return A;
}

// ------------------------------------------------------------------ per p: roots
// colleague-matrix eigenvalues (numpy chebcompanion form), real roots in [-1, 1]
int colleague_roots(const Vec& c, int n, double* out) {
    while (n > 0 && std::abs(c[n]) < 1e-14 * std::abs(c[0]))
        --n;
    if (n == 0) return 0;
    if (n == 1) {
        const double u = -c[0] / c[1];
        if (std::abs(u) <= 1 + 1e-10) {
            out[0] = u;
            return 1;
        }
        return 0;
    }
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, NG, NG> M;
    M.setZero(n, n);
    const double s = std::sqrt(0.5);
    M(0, 1) = M(1, 0) = s;
    for (int k = 1; k < n - 1; ++k)
        M(k, k + 1) = M(k + 1, k) = 0.5;
    for (int k = 0; k < n; ++k) {
        const double scl = (k == 0 ? 1.0 : s) / s;
        M(k, n - 1) -= c[k] / c[n] * scl * 0.5;
    }
    Eigen::EigenSolver<decltype(M)> es(M, false);
    int nr = 0;
    for (int k = 0; k < n; ++k) {
        const auto z = es.eigenvalues()[k];
        if (std::abs(z.imag()) < 1e-8 && std::abs(z.real()) <= 1 + 1e-10) out[nr++] = z.real();
    }
    std::sort(out, out + nr);
    return nr;
}
// cheap alternative: sign changes on a Chebyshev-Lobatto grid + Illinois (not a guarantee)
int grid_roots(const Vec& c, double* out, int M = 24) {
    int nr = 0;
    double u0 = -1, g0 = clenshaw(c, -1);
    for (int j = 1; j <= M; ++j) {
        const double u1 = -std::cos(PI * j / M), g1 = clenshaw(c, u1);
        if ((g0 < 0) != (g1 < 0)) {
            double a = u0, b = u1, fa = g0, fb = g1;
            int side = 0;
            for (int it = 0; it < 60 && b - a > 1e-15; ++it) {
                const double m = (a * fb - b * fa) / (fb - fa), fm = clenshaw(c, m);
                if ((fm < 0) == (fb < 0)) {
                    b = m;
                    fb = fm;
                    if (side == -1) fa *= 0.5;
                    side = -1;
                } else {
                    a = m;
                    fa = fm;
                    if (side == 1) fb *= 0.5;
                    side = 1;
                }
                if (fm == 0) break;
            }
            out[nr++] = 0.5 * (a + b);
        }
        u0 = u1;
        g0 = g1;
    }
    return nr;
}
// Certified subdivision (no eigensolve): on a subinterval [ua, ub] of the piece variable,
//   |c0| > sum|ck|        -> provably no root
//   same test on c'       -> monotone: at most one root, bracketed by endpoint signs
//   otherwise             -> re-expand on both halves and recurse
// A root is only missed if it is a tangency (double root) finer than the depth cap.
bool excluded(const Vec& c, int n) {
    double t = 0;
    for (int k = 1; k <= n; ++k)
        t += std::abs(c[k]);
    return std::abs(c[0]) > t;
}
Vec chebder(const Vec& c, int n) {  // derivative coefficients, degree n-1
    Vec d{};
    if (n >= 1) {
        d[n - 1] = 2 * n * c[n];
        if (n >= 2) d[n - 2] = 2 * (n - 1) * c[n - 1];
        for (int k = n - 3; k >= 0; --k)
            d[k] = d[k + 2] + 2 * (k + 1) * c[k + 1];
        d[0] *= 0.5;
    }
    return d;
}
// Re-expansion onto the left/right half is a fixed linear map: precompute both 14x14 matrices.
struct HalfMaps
{
    double L[NG + 1][NG + 1], Rm[NG + 1][NG + 1];
    HalfMaps() {
        for (int side = 0; side < 2; ++side)
            for (int col = 0; col <= NG; ++col) {  // image of T_col
                double f[NG + 1];
                for (int j = 0; j <= NG; ++j) {
                    const double x = std::cos(PI * (j + 0.5) / (NG + 1));
                    const double u = side == 0 ? (x - 1) / 2 : (x + 1) / 2;  // node in the parent variable
                    f[j] = std::cos(col * std::acos(u));
                }
                for (int k = 0; k <= NG; ++k) {
                    double s = 0;
                    for (int j = 0; j <= NG; ++j)
                        s += f[j] * std::cos(PI * k * (j + 0.5) / (NG + 1));
                    (side == 0 ? L : Rm)[k][col] = s * (k == 0 ? 1.0 : 2.0) / (NG + 1);
                }
            }
    }
};
const HalfMaps HM;
Vec reexpand_half(const Vec& c, bool right) {
    const auto& A = right ? HM.Rm : HM.L;
    Vec r{};
    for (int k = 0; k <= NG; ++k) {
        double s = 0;
        for (int j = 0; j <= NG; ++j)
            s += A[k][j] * c[j];
        r[k] = s;
    }
    return r;
}
double illinois(const Vec& c, double a, double b, double fa, double fb) {
    int side = 0;
    for (int it = 0; it < 80 && b - a > 1e-15; ++it) {
        const double m = (a * fb - b * fa) / (fb - fa), fm = clenshaw(c, m);
        if (fm == 0) return m;
        if ((fm < 0) == (fb < 0)) {
            b = m;
            fb = fm;
            if (side == -1) fa *= 0.5;
            side = -1;
        } else {
            a = m;
            fa = fm;
            if (side == 1) fb *= 0.5;
            side = 1;
        }
    }
    return 0.5 * (a + b);
}
// local: coefficients on the current subinterval; [ua, ub] its position in the piece variable
void cert_rec(const Vec& local, double ua, double ub, int depth, double* out, int& nr) {
    if (excluded(local, NG)) return;
    const Vec d = chebder(local, NG);
    if (excluded(d, NG - 1) || depth >= 12) {
        const double fa = clenshaw(local, -1), fb = clenshaw(local, 1);
        if ((fa < 0) != (fb < 0) || fa == 0) {
            const double v = fa == 0 ? -1.0 : illinois(local, -1, 1, fa, fb);
            out[nr++] = ua + (ub - ua) * (v + 1) / 2;
        }
        return;
    }
    const double um = 0.5 * (ua + ub);
    cert_rec(reexpand_half(local, false), ua, um, depth + 1, out, nr);
    cert_rec(reexpand_half(local, true), um, ub, depth + 1, out, nr);
}
int cert_roots(const Vec& c, double* out) {
    int nr = 0;
    cert_rec(c, -1, 1, 0, out, nr);
    return nr;
}
template <int METHOD>  // 0 eigen, 1 grid, 2 certified
int solve(const Assembled& A, double T, double p, double* rho) {
    const double t = p * A.q / (R * T);
    int n = 0;
    for (int pc = 0; pc < NP; ++pc) {
        Vec g;
        for (int k = 0; k <= NG; ++k)
            g[k] = A.G[pc][k] - t * A.P2[pc][k];
        double tail = 0;
        for (int k = 1; k <= NG; ++k)
            tail += std::abs(g[k]);
        if (std::abs(g[0]) > tail) continue;  // provably root-free piece
        double u[64];
        const int nr = METHOD == 0 ? colleague_roots(g, NG, u) : METHOD == 1 ? grid_roots(g, u) : cert_roots(g, u);
        for (int k = 0; k < nr; ++k) {
            const double e = EDGES[pc] + (EDGES[pc + 1] - EDGES[pc]) * (u[k] + 1) / 2;
            if (n > 0 && std::abs(e - rho[n - 1] * A.q) < 1e-10 * e) continue;  // shared edge
            rho[n++] = e / A.q;
        }
    }
    return n;
}

// ------------------------------------------------------------------ reference (direct Z)
double Z_direct(const Fluid& F, double T, double eta, const double* x) {
    const Scalars S = scalars(F, T, x);
    const double e = eta, om = 1 - e;
    double Z = 1 + S.mbar * (S.A1 * e / (om * om) + S.A2 * e * (1 + e) / (om * om * om) - S.A3 * e / om);
    for (int i = 0; i < F.N; ++i) {
        const double al = (3 - S.a[i]) / 3, be = (2 * S.a[i] - 3) / 3;
        Z -= S.w[i] * e * (3 / om - al / (1 - al * e) + be / (1 + be * e));
    }
    const D E{e, 1.0};
    const double cs[3] = {1, S.c1, S.c2};
    D I1 = cst(0), I2 = cst(0);
    for (int s = 0; s < 3; ++s) {
        I1 = I1 + cs[s] * Isf(GA, s, E);
        I2 = I2 + cs[s] * Isf(GB, s, E);
    }
    const D P = P0f(E) + S.mbar * P1f(E), K = Kf(E);
    const D f1 = E * I1;
    // eta C1 I2 = eta K I2 / P
    const D num = E * K * I2;
    const double dC = (num.d * P.v - num.v * P.d) / (P.v * P.v);
    Z += S.B1 * e * f1.d + S.B2 * e * dC;
    return Z;
}
int ref_roots(const Fluid& F, double T, const double* x, double p, double* rho) {
    const double q = scalars(F, T, x).q, t = p * q / (R * T);
    auto g = [&](double e) { return e * Z_direct(F, T, e, x) - t; };
    const int M = 200000;
    int n = 0;
    double e0 = 1e-12, g0 = g(e0);
    for (int j = 1; j <= M; ++j) {
        const double e1 = 0.74 * j / M, g1 = g(e1);
        if ((g0 < 0) != (g1 < 0)) {
            double a = e0, b = e1, fa = g0;
            for (int it = 0; it < 80; ++it) {
                const double m = 0.5 * (a + b), fm = g(m);
                if ((fm < 0) == (fa < 0)) {
                    a = m;
                    fa = fm;
                } else
                    b = m;
            }
            rho[n++] = 0.5 * (a + b) / q;
        }
        e0 = e1;
        g0 = g1;
    }
    return n;
}

// direct-fit comparison: evaluate Q = P^2 Z at the nodes of every piece, per (T, x)
Assembled assemble_direct(const Fluid& F, double T, const double* x) {
    const Scalars S = scalars(F, T, x);
    Assembled A{};
    A.q = S.q;
    A.ok = true;
    for (int p = 0; p < NP; ++p) {
        const auto c = chebfit<NQ>(
          [&](double e) {
              const double P = Pr(0, e) + S.mbar * Pr(1, e);
              return P * P * Z_direct(F, T, e, x);
          },
          EDGES[p], EDGES[p + 1]);
        A.G[p] = times_eta(c, EDGES[p], EDGES[p + 1]);
        const auto c2 = chebfit<NQ>(
          [&](double e) {
              const double P = Pr(0, e) + S.mbar * Pr(1, e);
              return P * P;
          },
          EDGES[p], EDGES[p + 1]);
        A.P2[p] = {};
        for (int k = 0; k <= NQ; ++k)
            A.P2[p][k] = c2[k];
    }
    return A;
}
}  // namespace

int main() {
    auto t0 = std::chrono::steady_clock::now();
    const Tables Tb = build_tables();
    const double t_build = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    std::printf("offline tables: %.1f ms, %zu doubles (%.0f kB)\n", t_build, (Tb.U.size() * NP + 3 * NP + 3 * NP * (NAD + 1)) * (NG + 1),
                (Tb.U.size() * NP + 3 * NP + 3 * NP * (NAD + 1)) * (NG + 1) * 8 / 1024.0);

    struct Case
    {
        const char* name;
        Fluid F;
        double x[NCMAX];
    };
    auto check_assembly = [&](const char* name, const Fluid& F, const double* x) {
        double worst = 0, worst_T = 0, worst_e = 0;
        for (double T : {150.0, 250.0, 400.0}) {
            const Assembled A = assemble(Tb, F, T, x);
            const Scalars S = scalars(F, T, x);
            for (int p = 0; p < NP; ++p)
                for (int j = 0; j <= 200; ++j) {
                    const double u = -1 + 2.0 * j / 200, e = EDGES[p] + (EDGES[p + 1] - EDGES[p]) * (u + 1) / 2;
                    if (e <= 0) continue;
                    const double P = Pr(0, e) + S.mbar * Pr(1, e), ref = e * P * P * Z_direct(F, T, e, x);
                    const double scale = e * P * P * (1 + std::abs(Z_direct(F, T, e, x)));
                    const double err = std::abs(clenshaw(A.G[p], u) - ref) / scale;
                    if (err > worst) {
                        worst = err;
                        worst_T = T;
                        worst_e = e;
                    }
                }
        }
        std::printf("assembly error %-16s max |G_asm - eta P^2 Z| / (eta P^2 (1+|Z|)) = %.1e  (T=%g, eta=%.4f)\n", name, worst, worst_T, worst_e);
    };
    const Case cases[] = {
      {"propane", {1, {2.0020}, {3.6184}, {208.11}}, {1.0}},
      {"C1-nC10 70/30", {2, {1.0, 4.6627}, {3.7039, 3.8384}, {150.03, 243.87}}, {0.7, 0.3}},
      {"C1C2C3 50/30/20", {3, {1.0, 1.6069, 2.0020}, {3.7039, 3.5206, 3.6184}, {150.03, 191.42, 208.11}}, {0.5, 0.3, 0.2}},
    };
    for (const auto& cs : cases)
        check_assembly(cs.name, cs.F, cs.x);
    const double Ts[] = {150, 200, 250, 300, 350, 400, 500};
    const double ps[] = {1e3, 1e4, 1e5, 1e6, 3e6, 1e7, 1e8};

    std::printf("\n%-16s %8s %9s %9s %9s %9s %9s %10s %9s %9s %9s\n", "case", "roots", "err_eig", "err_grid", "err_cert", "asm_us", "eig_us",
                "grid_us", "cert_us", "direct_us", "count_ok");
    for (const auto& cs : cases) {
        // accuracy + root-count check against a 200k-point scan of the direct model
        double worst = 0, worstg = 0, worstc = 0, wT = 0, wp = 0, wr = 0;
        int wn = 0;
        int nroots = 0, bad = 0;
        for (double T : Ts) {
            const Assembled A = assemble(Tb, cs.F, T, cs.x);
            if (!A.ok) {
                ++bad;
                continue;
            }
            for (double p : ps) {
                double r1[32], r2[32], r3[32], r4[32];
                const int n1 = solve<0>(A, T, p, r1), n2 = ref_roots(cs.F, T, cs.x, p, r2), n3 = solve<1>(A, T, p, r3), n4 = solve<2>(A, T, p, r4);
                if (n4 != n2) {
                    ++bad;
                    std::printf("   certified count mismatch %s T=%g p=%g: %d vs %d\n", cs.name, T, p, n4, n2);
                } else
                    for (int k = 0; k < n4; ++k) {
                        const double e = std::abs(r4[k] - r2[k]) / r2[k];
                        if (e > worstc) {
                            worstc = e;
                            wT = T;
                            wp = p;
                            wr = r2[k] * A.q;
                            wn = n2;
                        }
                    }
                nroots += n2;
                if (n1 != n2) {
                    ++bad;
                    std::printf("   count mismatch %s T=%g p=%g: eig %d ref %d\n", cs.name, T, p, n1, n2);
                    continue;
                }
                for (int k = 0; k < n1; ++k)
                    worst = std::max(worst, std::abs(r1[k] - r2[k]) / r2[k]);
                if (n3 == n2)
                    for (int k = 0; k < n3; ++k)
                        worstg = std::max(worstg, std::abs(r3[k] - r2[k]) / r2[k]);
                else
                    worstg = INFINITY;
            }
        }
        // timing
        volatile double sink = 0;
        const int REP = 20000;
        auto time = [&](auto&& fn) {
            double best = 1e300;
            for (int rep = 0; rep < 5; ++rep) {
                const auto a = std::chrono::steady_clock::now();
                for (int it = 0; it < REP; ++it)
                    fn(it);
                best = std::min(best, std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - a).count() / REP);
            }
            return best;
        };
        const double t_asm = time([&](int it) {
            const Assembled A = assemble(Tb, cs.F, Ts[it % 7] + 1e-9 * it, cs.x);
            sink = sink + A.G[0][0];
        });
        std::vector<Assembled> pre;
        for (double T : Ts)
            pre.push_back(assemble(Tb, cs.F, T, cs.x));
        const double t_eig = time([&](int it) {
            double r[32];
            sink = sink + solve<0>(pre[it % 7], Ts[it % 7], ps[(it / 7) % 7], r);
        });
        const double t_grid = time([&](int it) {
            double r[32];
            sink = sink + solve<1>(pre[it % 7], Ts[it % 7], ps[(it / 7) % 7], r);
        });
        const double t_cert = time([&](int it) {
            double r[32];
            sink = sink + solve<2>(pre[it % 7], Ts[it % 7], ps[(it / 7) % 7], r);
        });
        const double t_dir = time([&](int it) {
            const Assembled A = assemble_direct(cs.F, Ts[it % 7] + 1e-9 * it, cs.x);
            sink = sink + A.G[0][0];
        });
        std::printf("%-16s %8d %9.1e %9.1e %9.1e %9.2f %9.2f %10.2f %9.2f %9.2f %9s\n", cs.name, nroots, worst, worstg, worstc, t_asm, t_eig, t_grid,
                    t_cert, t_dir, bad ? "NO" : "yes");
        std::printf("   worst certified root: T=%g K p=%g Pa eta=%.6f (%d roots at this state)\n", wT, wp, wr, wn);
    }
}
