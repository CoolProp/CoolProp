// Univariate truncated Taylor polynomial, normalized coefficients t[k] = f^(k)/k!.
// Product: (N+1)(N+2)/2 FMAs (15 at N = 4); exp/log/division by the standard O(N^2)
// recurrences.  Plus the glue for building a bivariate Taylor2 from univariate pieces:
// embed along u or v, and compose a univariate expansion with a bivariate argument.
#pragma once
#include <array>
#include <cmath>
#include "taylor2.hpp"

namespace nd {

template <int N>
struct Taylor1
{
    std::array<double, N + 1> t{};
    Taylor1() = default;
    Taylor1(double v) {
        t[0] = v;
    }
    static Taylor1 variable(double x0, double slope = 1.0) {
        Taylor1 r(x0);
        if constexpr (N >= 1) r.t[1] = slope;
        return r;
    }
};

template <int N>
Taylor1<N> operator*(const Taylor1<N>& a, const Taylor1<N>& b) {
    Taylor1<N> r;
    for (int k = 0; k <= N; ++k) {
        double s = 0;
        for (int j = 0; j <= k; ++j)
            s += a.t[j] * b.t[k - j];
        r.t[k] = s;
    }
    return r;
}
template <int N>
Taylor1<N> operator/(const Taylor1<N>& a, const Taylor1<N>& b) {
    Taylor1<N> q;
    const double inv0 = 1.0 / b.t[0];
    for (int k = 0; k <= N; ++k) {
        double s = a.t[k];
        for (int j = 1; j <= k; ++j)
            s -= b.t[j] * q.t[k - j];
        q.t[k] = s * inv0;
    }
    return q;
}
// e = exp(g): e_k = (1/k) sum_{j=1..k} j g_j e_{k-j}
template <int N>
Taylor1<N> exp(const Taylor1<N>& g) {
    Taylor1<N> e;
    e.t[0] = std::exp(g.t[0]);
    for (int k = 1; k <= N; ++k) {
        double s = 0;
        for (int j = 1; j <= k; ++j)
            s += j * g.t[j] * e.t[k - j];
        e.t[k] = s / k;
    }
    return e;
}
// l = log(g): l_k = (g_k - (1/k) sum_{j=1..k-1} j l_j g_{k-j}) / g_0
template <int N>
Taylor1<N> log(const Taylor1<N>& g) {
    Taylor1<N> l;
    l.t[0] = std::log(g.t[0]);
    const double inv0 = 1.0 / g.t[0];
    for (int k = 1; k <= N; ++k) {
        double s = 0;
        for (int j = 1; j < k; ++j)
            s += j * l.t[j] * g.t[k - j];
        l.t[k] = (g.t[k] - s / k) * inv0;
    }
    return l;
}
template <int N>
Taylor1<N> operator+(Taylor1<N> a, const Taylor1<N>& b) {
    for (int k = 0; k <= N; ++k)
        a.t[k] += b.t[k];
    return a;
}
template <int N>
Taylor1<N> operator-(Taylor1<N> a, const Taylor1<N>& b) {
    for (int k = 0; k <= N; ++k)
        a.t[k] -= b.t[k];
    return a;
}
template <int N>
Taylor1<N> operator*(Taylor1<N> a, double s) {
    for (auto& v : a.t)
        v *= s;
    return a;
}
template <int N>
Taylor1<N> operator*(double s, Taylor1<N> a) {
    return a * s;
}
template <int N>
Taylor1<N> operator+(Taylor1<N> a, double s) {
    a.t[0] += s;
    return a;
}
template <int N>
Taylor1<N> operator+(double s, Taylor1<N> a) {
    return a + s;
}
template <int N>
Taylor1<N> operator-(double s, Taylor1<N> a) {
    for (auto& v : a.t)
        v = -v;
    a.t[0] += s;
    return a;
}
template <int N>
Taylor1<N> operator/(double s, const Taylor1<N>& a) {
    return Taylor1<N>(s) / a;
}
template <int N>
Taylor1<N>& operator+=(Taylor1<N>& a, const Taylor1<N>& b) {
    return a = a + b;
}
template <int N>
Taylor1<N>& operator*=(Taylor1<N>& a, const Taylor1<N>& b) {
    return a = a * b;
}

// ---- bivariate glue
template <int N>
Taylor2<N> embed_u(const Taylor1<N>& a) {  // a(u) as a function of (u, v)
    Taylor2<N> r;
    for (int i = 0; i <= N; ++i)
        r.c[Taylor2<N>::idx(i, 0)] = a.t[i];
    return r;
}
// (c(u) * B(u,v)): univariate-in-u times bivariate, sum over a <= i only
template <int N>
Taylor2<N> mul_u(const Taylor1<N>& c, const Taylor2<N>& B) {
    Taylor2<N> r;
    for (int k = 0; k <= N; ++k)
        for (int j = 0; j <= k; ++j) {
            const int i = k - j;
            double s = 0;
            for (int a = 0; a <= i; ++a)
                s += c.t[a] * B.c[Taylor2<N>::idx(i - a, j)];
            r.c[Taylor2<N>::idx(i, j)] = s;
        }
    return r;
}
// Powers H^1..H^N of a bivariate increment H (H(0,0) = 0), computed once and shared by
// every composition phi(eta(u,v)) = sum_k phi_k H^k.
template <int N>
struct PowerBasis
{
    std::array<Taylor2<N>, N + 1> H;
    explicit PowerBasis(const Taylor2<N>& h) {
        H[0] = Taylor2<N>(1.0);
        if constexpr (N >= 1) H[1] = h;
        for (int k = 2; k <= N; ++k)
            H[k] = H[k - 1] * h;
    }
    Taylor2<N> compose(const Taylor1<N>& phi) const {
        Taylor2<N> r(phi.t[0]);
        for (int k = 1; k <= N; ++k)
            for (int m = 0; m < Taylor2<N>::M; ++m)
                r.c[m] += phi.t[k] * H[k].c[m];
        return r;
    }
};

}  // namespace nd
