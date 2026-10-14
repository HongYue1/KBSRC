// Double-double arithmetic (~106-bit significand) built on error-free transforms
#pragma once
#include <cmath>
#include "isa.hpp"

namespace kbsrc::KBSRC_ISA::dsp {

// Exact error a*b - p of the rounded product p = fl(a*b). Hardware FMA when available, otherwise Dekker's
// split plus a scaled path for tiny products, which returns the same bits as the FMA in round-to-nearest (the
// mode Exact always runs in, see dispatch.cpp) whenever the split cannot overflow (|a|, |b| up to about 2^995,
// |p| up to about 2^1000; larger values behave as before). This is the only arithmetic that differs between
// the engines, so it is what keeps the generic and AVX2 engines bit-identical.
//   prod_err_fast  the same without the tiny-product path, for vectorised loops: equal to prod_err unless
//                  prod_err_tiny(a, b, p), which is always false with FMA.
#if defined(__FMA__) || defined(__ARM_FEATURE_FMA) || defined(__aarch64__)
// __builtin_fma, not std::fma: unoptimised builds compile std::fma to a call into the C runtime, and MinGW-w64's
// fma is not a correctly rounded fused multiply-add (it differs from the instruction on about one random product
// in ten), which made the AVX2 engine's Exact output wrong in Debug builds. The builtin is always the
// instruction; optimised code is unchanged.
static inline double prod_err(double a, double b, double p) {
    return __builtin_fma(a, b, -p);
}
static inline double prod_err_fast(double a, double b, double p) { return prod_err(a, b, p); }
static inline bool prod_err_tiny(double, double, double) { return false; }
#else
// Dekker's split. Exact, hence equal to the FMA result, while no partial product underflows: guaranteed for
// |p| >= 2^-968, since then e_a + e_b >= -970 and every partial product is a multiple of
// ulp(a) * ulp(b) >= 2^-1074.
static inline double prod_err_fast(double a, double b, double p) {
    const double S = 134217729.0;   // 2^27 + 1
    double t = S * a, ah = t - (t - a), al = a - ah;
    t = S * b; double bh = t - (t - b), bl = b - bh;
    return ((ah * bh - p) + ah * bl + al * bh) + al * bl;
}
// The only inputs where the split can be inexact: a tiny product of nonzero factors (subnormal operands or
// results, e.g. denormal samples). Zero factors give +0 on both paths; a NaN or infinite p fails the compare.
// Bitwise & instead of && keeps loops that call this branch-free (vectorisable).
static inline bool prod_err_tiny(double a, double b, double p) {
    return (std::fabs(p) < 0x1p-968) & (a != 0) & (b != 0);
}
// Cold path for prod_err_tiny. Here |a|, |b| < 2^106 (the other factor is at least 2^-1074), so scaling both by
// 2^589 is exact and cannot overflow, and it lifts e_a + e_b back to >= -970: the split of the scaled factors
// is exact. P - ps is exact (Sterbenz); where p was rounded more coarsely than P (|a*b - p| <= 2^-1075) every
// path gives a signed zero, and otherwise P == ps and the error is scaled back with a single final rounding.
// Either way the result is fma(a, b, -p) bit for bit.
[[gnu::noinline, gnu::cold]] static inline double prod_err_scaled(double a, double b, double p) {
    const double up = 0x1p589, down = 0x1p-589;
    double as = a * up, bs = b * up, ps = p * up * up, P = as * bs;
    return ((P - ps) + prod_err_fast(as, bs, P)) * down * down;
}
static inline double prod_err(double a, double b, double p) {
    if (std::fabs(p) < 0x1p-968) [[unlikely]] {
        if (a != 0 && b != 0) return prod_err_scaled(a, b, p);
    }
    return prod_err_fast(a, b, p);
}
#endif

struct dd {
    double hi, lo;
    dd() : hi(0), lo(0) {}
    dd(double h) : hi(h), lo(0) {}
    dd(double h, double l) : hi(h), lo(l) {}
};
static inline dd qts(double a, double b) { double s = a + b; return dd(s, b - (s - a)); }
static inline dd ts(double a, double b) { double s = a + b, bb = s - a; return dd(s, (a - (s - bb)) + (b - bb)); }
static inline dd operator+(const dd& a, const dd& b) {
    dd s = ts(a.hi, b.hi), t = ts(a.lo, b.lo);
    s.lo += t.hi; s = qts(s.hi, s.lo); s.lo += t.lo; return qts(s.hi, s.lo);
}
static inline dd operator-(const dd& a) { return dd(-a.hi, -a.lo); }
static inline dd operator-(const dd& a, const dd& b) { return a + (-b); }
static inline dd operator*(const dd& a, const dd& b) {
    double p = a.hi * b.hi, e = prod_err(a.hi, b.hi, p);
    e += a.hi * b.lo + a.lo * b.hi; return qts(p, e);
}
static inline dd operator*(const dd& a, double b) {
    double p = a.hi * b, e = prod_err(a.hi, b, p); e += a.lo * b; return qts(p, e);
}
static inline dd operator/(const dd& a, double b) {
    double q1 = a.hi / b; dd r = a - dd(b) * q1;
    double q2 = r.hi / b; r = r - dd(b) * q2;
    double q3 = r.hi / b; return qts(q1, q2) + dd(q3);
}
static inline dd& operator+=(dd& a, const dd& b) { a = a + b; return a; }
static inline dd& operator-=(dd& a, const dd& b) { a = a - b; return a; }
static inline dd& operator*=(dd& a, const dd& b) { a = a * b; return a; }
static inline dd operator/(const dd& a, const dd& b) {  // long division, 3 steps
    double q1 = a.hi / b.hi; dd r = a - b * q1;
    double q2 = r.hi / b.hi; r = r - b * q2;
    double q3 = r.hi / b.hi; return qts(q1, q2) + dd(q3);
}
// sqrt: one Newton step on the double estimate (exact residual via prod_err)
static inline dd dd_sqrt(const dd& a) {
    if (a.hi <= 0) return dd(0.0);
    double y = std::sqrt(a.hi); dd r = a - dd(y) * y;
    return qts(y, r.hi / (2 * y));
}
// modified Bessel function I0, power series in double-double
static inline dd dd_i0(const dd& x) {
    dd q = x * x / 4.0, t(1.0), s(1.0);
    for (int k = 1; k < 2000; k++) {
        t = t * q / (double(k) * double(k)); s += t;
        if (t.hi < s.hi * 1e-34) break;
    }
    return s;
}
static const dd DD_PI(3.141592653589793116e+00, 1.224646799147353207e-16);

// cos/sin of 2*pi*k/N, exact integer argument reduction, Taylor in dd on [0, pi/4]
static inline void dd_sincos_small(const dd& x, dd& c, dd& s) {
    dd x2 = x * x, t = x, ss = x, cc(1.0), tc(1.0);
    for (int k = 1; k < 40; k++) {
        tc = -(tc * x2) / double((2 * k - 1) * (2 * k)); cc += tc;
        t = -(t * x2) / double((2 * k) * (2 * k + 1)); ss += t;
        if (std::fabs(t.hi) < 1e-40 && std::fabs(tc.hi) < 1e-40) break;
    }
    c = cc; s = ss;
}
static inline void dd_cis(long long k, long long N, dd& c, dd& s) {  // e^{+i 2 pi k/N} = c + i s
    k %= N; if (k < 0) k += N;
    long long m = 8 * k, oct = m / N, rem = m - oct * N;
    bool odd = oct & 1; long long j = oct / 2;
    dd phi = (dd(DD_PI.hi / 4, DD_PI.lo / 4) * double(odd ? (N - rem) : rem)) / double(N);
    dd cc, sn; dd_sincos_small(phi, cc, sn);
    if (odd) { sn = -sn; j += 1; }
    switch (j & 3) {
        case 0: c = cc; s = sn; break;
        case 1: c = -sn; s = cc; break;
        case 2: c = -cc; s = -sn; break;
        default: c = sn; s = -cc; break;
    }
}

}  // namespace kbsrc::KBSRC_ISA::dsp
