// A deliberately small C++ analogue of the Rust num-dual crate used by FeOs
// (Rehner & Bauer, doi:10.3389/fceng.2021.758090): plain structs, eager evaluation,
// no expression templates.  Each type stores exactly the derivative coefficients it
// needs and every elementary function goes through one chain-rule helper, so adding a
// function means supplying f, f', f'', f''' at the real part.
//
//   Dual<T>         f, f'                       (nests: Dual<Dual<T>> etc.)
//   Dual2<T>        f, f', f''                  univariate, 2nd order in one pass
//   Dual3<T>        f, f', f'', f'''            univariate, 3rd order in one pass
//   HyperDual<T>    f, f_x, f_y, f_xy           mixed 2nd derivative in one pass
//   DualVec<T, N>   f, grad f (N directions)    whole gradient in one pass
//
// Spike simplification: the scalar is templated for shape but the elementary-function
// kernels assume T = double (num-dual proper is generic, which is what allows nesting).
#pragma once
#include <array>
#include <cmath>
#include <type_traits>

namespace nd {

// ---------------------------------------------------------------- Dual
template <class T>
struct Dual
{
    T re{}, eps{};
    Dual() = default;
    Dual(double r) : re(r) {}
    Dual(T r, T e) : re(r), eps(e) {}
};
template <class T>
Dual<T> operator*(const Dual<T>& a, const Dual<T>& b) {
    return {a.re * b.re, a.eps * b.re + a.re * b.eps};
}
template <class T>
Dual<T> chain(const Dual<T>& g, T f0, T f1, T, T) {
    return {f0, f1 * g.eps};
}
template <class T>
Dual<T> scale(const Dual<T>& a, double s) {
    return {a.re * s, a.eps * s};
}
template <class T>
Dual<T> add(const Dual<T>& a, const Dual<T>& b) {
    return {a.re + b.re, a.eps + b.eps};
}
template <class T>
Dual<T> addre(const Dual<T>& a, double s) {
    return {a.re + s, a.eps};
}

// ---------------------------------------------------------------- Dual2
template <class T>
struct Dual2
{
    T re{}, v1{}, v2{};
    Dual2() = default;
    Dual2(double r) : re(r) {}
    Dual2(T r, T a, T b) : re(r), v1(a), v2(b) {}
};
template <class T>
Dual2<T> operator*(const Dual2<T>& a, const Dual2<T>& b) {
    return {a.re * b.re, a.v1 * b.re + a.re * b.v1, a.v2 * b.re + 2.0 * a.v1 * b.v1 + a.re * b.v2};
}
template <class T>
Dual2<T> chain(const Dual2<T>& g, T f0, T f1, T f2, T) {
    return {f0, f1 * g.v1, f2 * g.v1 * g.v1 + f1 * g.v2};
}
template <class T>
Dual2<T> scale(const Dual2<T>& a, double s) {
    return {a.re * s, a.v1 * s, a.v2 * s};
}
template <class T>
Dual2<T> add(const Dual2<T>& a, const Dual2<T>& b) {
    return {a.re + b.re, a.v1 + b.v1, a.v2 + b.v2};
}
template <class T>
Dual2<T> addre(const Dual2<T>& a, double s) {
    return {a.re + s, a.v1, a.v2};
}

// ---------------------------------------------------------------- Dual3
template <class T>
struct Dual3
{
    T re{}, v1{}, v2{}, v3{};
    Dual3() = default;
    Dual3(double r) : re(r) {}
    Dual3(T r, T a, T b, T c) : re(r), v1(a), v2(b), v3(c) {}
};
template <class T>
Dual3<T> operator*(const Dual3<T>& a, const Dual3<T>& b) {
    return {a.re * b.re, a.v1 * b.re + a.re * b.v1, a.v2 * b.re + 2.0 * a.v1 * b.v1 + a.re * b.v2,
            a.v3 * b.re + 3.0 * (a.v2 * b.v1 + a.v1 * b.v2) + a.re * b.v3};
}
template <class T>
Dual3<T> chain(const Dual3<T>& g, T f0, T f1, T f2, T f3) {
    return {f0, f1 * g.v1, f2 * g.v1 * g.v1 + f1 * g.v2, f3 * g.v1 * g.v1 * g.v1 + 3.0 * f2 * g.v1 * g.v2 + f1 * g.v3};
}
template <class T>
Dual3<T> scale(const Dual3<T>& a, double s) {
    return {a.re * s, a.v1 * s, a.v2 * s, a.v3 * s};
}
template <class T>
Dual3<T> add(const Dual3<T>& a, const Dual3<T>& b) {
    return {a.re + b.re, a.v1 + b.v1, a.v2 + b.v2, a.v3 + b.v3};
}
template <class T>
Dual3<T> addre(const Dual3<T>& a, double s) {
    return {a.re + s, a.v1, a.v2, a.v3};
}

// ---------------------------------------------------------------- HyperDual
template <class T>
struct HyperDual
{
    T re{}, e1{}, e2{}, e12{};
    HyperDual() = default;
    HyperDual(double r) : re(r) {}
    HyperDual(T r, T a, T b, T c) : re(r), e1(a), e2(b), e12(c) {}
};
template <class T>
HyperDual<T> operator*(const HyperDual<T>& a, const HyperDual<T>& b) {
    return {a.re * b.re, a.e1 * b.re + a.re * b.e1, a.e2 * b.re + a.re * b.e2, a.e12 * b.re + a.e1 * b.e2 + a.e2 * b.e1 + a.re * b.e12};
}
template <class T>
HyperDual<T> chain(const HyperDual<T>& g, T f0, T f1, T f2, T) {
    return {f0, f1 * g.e1, f1 * g.e2, f1 * g.e12 + f2 * g.e1 * g.e2};
}
template <class T>
HyperDual<T> scale(const HyperDual<T>& a, double s) {
    return {a.re * s, a.e1 * s, a.e2 * s, a.e12 * s};
}
template <class T>
HyperDual<T> add(const HyperDual<T>& a, const HyperDual<T>& b) {
    return {a.re + b.re, a.e1 + b.e1, a.e2 + b.e2, a.e12 + b.e12};
}
template <class T>
HyperDual<T> addre(const HyperDual<T>& a, double s) {
    return {a.re + s, a.e1, a.e2, a.e12};
}

// ---------------------------------------------------------------- DualVec
template <class T, int N>
struct DualVec
{
    T re{};
    std::array<T, N> eps{};
    DualVec() = default;
    DualVec(double r) : re(r) {}
    DualVec(T r, const std::array<T, N>& e) : re(r), eps(e) {}
};
template <class T, int N>
DualVec<T, N> operator*(const DualVec<T, N>& a, const DualVec<T, N>& b) {
    DualVec<T, N> r;
    r.re = a.re * b.re;
    for (int i = 0; i < N; ++i)
        r.eps[i] = a.eps[i] * b.re + a.re * b.eps[i];
    return r;
}
template <class T, int N>
DualVec<T, N> chain(const DualVec<T, N>& g, T f0, T f1, T, T) {
    DualVec<T, N> r;
    r.re = f0;
    for (int i = 0; i < N; ++i)
        r.eps[i] = f1 * g.eps[i];
    return r;
}
template <class T, int N>
DualVec<T, N> scale(const DualVec<T, N>& a, double s) {
    DualVec<T, N> r;
    r.re = a.re * s;
    for (int i = 0; i < N; ++i)
        r.eps[i] = a.eps[i] * s;
    return r;
}
template <class T, int N>
DualVec<T, N> add(const DualVec<T, N>& a, const DualVec<T, N>& b) {
    DualVec<T, N> r;
    r.re = a.re + b.re;
    for (int i = 0; i < N; ++i)
        r.eps[i] = a.eps[i] + b.eps[i];
    return r;
}
template <class T, int N>
DualVec<T, N> addre(const DualVec<T, N>& a, double s) {
    DualVec<T, N> r = a;
    r.re += s;
    return r;
}

// ---------------------------------------------------------------- generic operators
// Everything below is written once against the (chain, scale, add, addre, *) kernel.
template <class D>
struct is_nd : std::false_type
{
};
template <class T>
struct is_nd<Dual<T>> : std::true_type
{
};
template <class T>
struct is_nd<Dual2<T>> : std::true_type
{
};
template <class T>
struct is_nd<Dual3<T>> : std::true_type
{
};
template <class T>
struct is_nd<HyperDual<T>> : std::true_type
{
};
template <class T, int N>
struct is_nd<DualVec<T, N>> : std::true_type
{
};
template <class D>
using if_nd = std::enable_if_t<is_nd<D>::value, int>;

template <class D, if_nd<D> = 0>
D operator+(const D& a, const D& b) {
    return add(a, b);
}
template <class D, if_nd<D> = 0>
D operator+(const D& a, double s) {
    return addre(a, s);
}
template <class D, if_nd<D> = 0>
D operator+(double s, const D& a) {
    return addre(a, s);
}
template <class D, if_nd<D> = 0>
D operator-(const D& a) {
    return scale(a, -1.0);
}
template <class D, if_nd<D> = 0>
D operator-(const D& a, const D& b) {
    return add(a, scale(b, -1.0));
}
template <class D, if_nd<D> = 0>
D operator-(const D& a, double s) {
    return addre(a, -s);
}
template <class D, if_nd<D> = 0>
D operator-(double s, const D& a) {
    return addre(scale(a, -1.0), s);
}
template <class D, if_nd<D> = 0>
D operator*(const D& a, double s) {
    return scale(a, s);
}
template <class D, if_nd<D> = 0>
D operator*(double s, const D& a) {
    return scale(a, s);
}
template <class D, if_nd<D> = 0>
D recip(const D& a) {
    const double x = a.re, r = 1.0 / x;
    return chain(a, r, -r * r, 2.0 * r * r * r, -6.0 * r * r * r * r);
}
template <class D, if_nd<D> = 0>
D operator/(const D& a, const D& b) {
    return a * recip(b);
}
template <class D, if_nd<D> = 0>
D operator/(const D& a, double s) {
    return scale(a, 1.0 / s);
}
template <class D, if_nd<D> = 0>
D operator/(double s, const D& a) {
    return scale(recip(a), s);
}
template <class D, if_nd<D> = 0>
D& operator+=(D& a, const D& b) {
    return a = a + b;
}
template <class D, if_nd<D> = 0>
D& operator-=(D& a, const D& b) {
    return a = a - b;
}
template <class D, if_nd<D> = 0>
D& operator*=(D& a, const D& b) {
    return a = a * b;
}
template <class D, if_nd<D> = 0>
D exp(const D& a) {
    const double e = std::exp(a.re);
    return chain(a, e, e, e, e);
}
template <class D, if_nd<D> = 0>
D log(const D& a) {
    const double r = 1.0 / a.re;
    return chain(a, std::log(a.re), r, -r * r, 2.0 * r * r * r);
}

}  // namespace nd
