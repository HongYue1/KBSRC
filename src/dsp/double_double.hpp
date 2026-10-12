// Double-double arithmetic (~106-bit significand) built on error-free transforms
#pragma once
#include <cmath>
#include "isa.hpp"

namespace kbsrc::KBSRC_ISA::dsp {

// Exact error a*b - p of the rounded product p = fl(a*b). Hardware FMA when available, otherwise
// Dekker's split (also exact while no partial product underflows, which audio-range values never do).
static inline double prod_err(double a, double b, double p) {
#if defined(__FMA__) || defined(__ARM_FEATURE_FMA) || defined(__aarch64__)
    return std::fma(a, b, -p);
#else
    const double S = 134217729.0;   // 2^27 + 1
    double t = S * a, ah = t - (t - a), al = a - ah;
    t = S * b; double bh = t - (t - b), bl = b - bh;
    return ((ah * bh - p) + ah * bl + al * bh) + al * bl;
#endif
}

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
