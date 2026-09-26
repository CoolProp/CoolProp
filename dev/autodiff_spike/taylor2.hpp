// Bivariate truncated Taylor polynomial of total degree N: carries every normalized
// coefficient c_ij = f_ij / (i! j!) with i + j <= N, i.e. the whole "half matrix" of
// derivatives, in one forward pass.  (N+1)(N+2)/2 coefficients; a product costs
// C(N+4, 4) multiply-adds (70 for N = 4).
//
// Eager, no expression templates, same design philosophy as numdual.hpp.
#pragma once
#include <array>
#include <cmath>

namespace nd {

template <int N>
struct Taylor2
{
    static constexpr int M = (N + 1) * (N + 2) / 2;
    // packed by total degree k = i + j, then by j: (0,0) (1,0) (0,1) (2,0) (1,1) (0,2) ...
    static constexpr int idx(int i, int j) {
        return (i + j) * (i + j + 1) / 2 + j;
    }
    std::array<double, M> c{};
    Taylor2() = default;
    Taylor2(double v) {
        c[0] = v;
    }
    double operator()(int i, int j) const {
        return c[idx(i, j)];
    }
};

// Flattened product table: for output m, the (a, b) index pairs with a + b = m in
// multi-index sense.  Built at compile time so the product is a straight run of FMAs.
template <int N>
struct ProdTable
{
    static constexpr int M = Taylor2<N>::M;
    static constexpr int P = (N + 1) * (N + 2) * (N + 3) * (N + 4) / 24;  // C(N+4, 4)
    std::array<int, M + 1> start{};
    std::array<int, P> a{}, b{};
    constexpr ProdTable() {
        int n = 0;
        for (int k = 0; k <= N; ++k)
            for (int j = 0; j <= k; ++j) {
                const int i = k - j;
                start[Taylor2<N>::idx(i, j)] = n;
                // (0,0) first so division can skip it by starting at start[m] + 1
                for (int ia = 0; ia <= i; ++ia)
                    for (int jb = 0; jb <= j; ++jb) {
                        a[n] = Taylor2<N>::idx(ia, jb);
                        b[n] = Taylor2<N>::idx(i - ia, j - jb);
                        ++n;
                    }
            }
        start[M] = n;
    }
};
template <int N>
inline constexpr ProdTable<N> prod_table{};

template <int N>
Taylor2<N> operator*(const Taylor2<N>& x, const Taylor2<N>& y) {
    constexpr auto& t = prod_table<N>;
    Taylor2<N> r;
    // Full unrolling lets the table indices constant-fold; without it clang keeps an
    // indexed, latency-bound loop from N = 3 upward (~10x slower).
#pragma clang loop unroll(full)
    for (int m = 0; m < Taylor2<N>::M; ++m) {
        double s = 0;
#pragma clang loop unroll(full)
        for (int n = t.start[m]; n < t.start[m + 1]; ++n)
            s += x.c[t.a[n]] * y.c[t.b[n]];
        r.c[m] = s;
    }
    return r;
}
// q = x / y from y*q = x, solved coefficient by coefficient in increasing degree
template <int N>
Taylor2<N> operator/(const Taylor2<N>& x, const Taylor2<N>& y) {
    constexpr auto& t = prod_table<N>;
    Taylor2<N> q;
    const double inv0 = 1.0 / y.c[0];
#pragma clang loop unroll(full)
    for (int m = 0; m < Taylor2<N>::M; ++m) {
        double s = x.c[m];
#pragma clang loop unroll(full)
        for (int n = t.start[m] + 1; n < t.start[m + 1]; ++n)
            s -= y.c[t.a[n]] * q.c[t.b[n]];
        q.c[m] = s * inv0;
    }
    return q;
}
// f(g) = sum_k fk[k] h^k with h = g - g0 and fk[k] = f^(k)(g0)/k!, by Horner in h
template <int N>
Taylor2<N> compose(const Taylor2<N>& g, const std::array<double, N + 1>& fk) {
    Taylor2<N> h = g;
    h.c[0] = 0.0;
    Taylor2<N> r(fk[N]);
    for (int k = N - 1; k >= 0; --k) {
        r = r * h;
        r.c[0] += fk[k];
    }
    return r;
}
template <int N>
Taylor2<N> exp(const Taylor2<N>& g) {
    std::array<double, N + 1> fk;
    fk[0] = std::exp(g.c[0]);
    for (int k = 1; k <= N; ++k)
        fk[k] = fk[k - 1] / k;
    return compose(g, fk);
}
template <int N>
Taylor2<N> log(const Taylor2<N>& g) {
    std::array<double, N + 1> fk;
    const double r = 1.0 / g.c[0];
    fk[0] = std::log(g.c[0]);
    double rk = 1.0;
    for (int k = 1; k <= N; ++k) {
        rk *= r;
        fk[k] = ((k % 2) ? 1.0 : -1.0) * rk / k;
    }
    return compose(g, fk);
}

template <int N>
Taylor2<N> operator+(Taylor2<N> a, const Taylor2<N>& b) {
    for (int m = 0; m < Taylor2<N>::M; ++m)
        a.c[m] += b.c[m];
    return a;
}
template <int N>
Taylor2<N> operator-(Taylor2<N> a, const Taylor2<N>& b) {
    for (int m = 0; m < Taylor2<N>::M; ++m)
        a.c[m] -= b.c[m];
    return a;
}
template <int N>
Taylor2<N> operator-(Taylor2<N> a) {
    for (auto& v : a.c)
        v = -v;
    return a;
}
template <int N>
Taylor2<N> operator*(Taylor2<N> a, double s) {
    for (auto& v : a.c)
        v *= s;
    return a;
}
template <int N>
Taylor2<N> operator*(double s, Taylor2<N> a) {
    return a * s;
}
template <int N>
Taylor2<N> operator/(Taylor2<N> a, double s) {
    return a * (1.0 / s);
}
template <int N>
Taylor2<N> operator/(double s, const Taylor2<N>& a) {
    return Taylor2<N>(s) / a;
}
template <int N>
Taylor2<N> operator+(Taylor2<N> a, double s) {
    a.c[0] += s;
    return a;
}
template <int N>
Taylor2<N> operator+(double s, Taylor2<N> a) {
    return a + s;
}
template <int N>
Taylor2<N> operator-(Taylor2<N> a, double s) {
    a.c[0] -= s;
    return a;
}
template <int N>
Taylor2<N> operator-(double s, const Taylor2<N>& a) {
    return (-a) + s;
}
template <int N>
Taylor2<N>& operator+=(Taylor2<N>& a, const Taylor2<N>& b) {
    return a = a + b;
}
template <int N>
Taylor2<N>& operator-=(Taylor2<N>& a, const Taylor2<N>& b) {
    return a = a - b;
}
template <int N>
Taylor2<N>& operator*=(Taylor2<N>& a, const Taylor2<N>& b) {
    return a = a * b;
}

}  // namespace nd
