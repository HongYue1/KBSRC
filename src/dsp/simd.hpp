// 4-lane types: V4 (4 x double) and DD4 (4 x double-double).
// AVX2+FMA intrinsics when the target has them, otherwise plain arrays (auto-vectorised: SSE2, NEON, ...).
// Both backends perform the same IEEE operations in the same order, and prod_err returns the FMA bits on both
// (round-to-nearest, see double_double.hpp), so results are bit-identical.
#pragma once
#include "dsp/double_double.hpp"
#include "isa.hpp"

#if defined(__AVX2__) && defined(__FMA__)
#define KBSRC_SIMD_AVX2 1
#include <immintrin.h>
#else
#define KBSRC_SIMD_AVX2 0
#endif

namespace kbsrc::KBSRC_ISA::dsp {

// scratch storage unit (32-byte aligned) for FFT work buffers
struct alignas(32) Scratch32 { double d[4]; };

#if KBSRC_SIMD_AVX2
inline constexpr const char* simd_name = "avx2";
struct V4 {
    __m256d v;
    V4() : v(_mm256_setzero_pd()) {}
    V4(__m256d x) : v(x) {}
    V4(double x) : v(_mm256_set1_pd(x)) {}
};
static inline V4 operator+(V4 a, V4 b) { return _mm256_add_pd(a.v, b.v); }
static inline V4 operator-(V4 a, V4 b) { return _mm256_sub_pd(a.v, b.v); }
static inline V4 operator*(V4 a, V4 b) { return _mm256_mul_pd(a.v, b.v); }
static inline V4 operator-(V4 a) { return _mm256_xor_pd(a.v, _mm256_set1_pd(-0.0)); }

struct DD4 {
    __m256d hi, lo;
    DD4() : hi(_mm256_setzero_pd()), lo(_mm256_setzero_pd()) {}
    DD4(__m256d h, __m256d l) : hi(h), lo(l) {}
    DD4(double x) : hi(_mm256_set1_pd(x)), lo(_mm256_setzero_pd()) {}
    DD4(const dd& x) : hi(_mm256_set1_pd(x.hi)), lo(_mm256_set1_pd(x.lo)) {}
};
static inline DD4 qts4(__m256d a, __m256d b) {
    __m256d s = _mm256_add_pd(a, b);
    return DD4(s, _mm256_sub_pd(b, _mm256_sub_pd(s, a)));
}
static inline DD4 operator+(const DD4& a, const DD4& b) {
    __m256d s = _mm256_add_pd(a.hi, b.hi), bb = _mm256_sub_pd(s, a.hi);
    __m256d e = _mm256_add_pd(_mm256_sub_pd(a.hi, _mm256_sub_pd(s, bb)), _mm256_sub_pd(b.hi, bb));
    e = _mm256_add_pd(e, _mm256_add_pd(a.lo, b.lo));
    return qts4(s, e);
}
static inline DD4 operator-(const DD4& a) {
    __m256d m = _mm256_set1_pd(-0.0);
    return DD4(_mm256_xor_pd(a.hi, m), _mm256_xor_pd(a.lo, m));
}
static inline DD4 operator-(const DD4& a, const DD4& b) { return a + (-b); }
static inline DD4 operator*(const DD4& a, const DD4& b) {
    __m256d p = _mm256_mul_pd(a.hi, b.hi), e = _mm256_fmsub_pd(a.hi, b.hi, p);
    e = _mm256_add_pd(e, _mm256_add_pd(_mm256_mul_pd(a.hi, b.lo), _mm256_mul_pd(a.lo, b.hi)));
    return qts4(p, e);
}
// lane I/O: d[0..3] <-> lanes 0..3
static inline V4 v4_load(const double* d) { return _mm256_loadu_pd(d); }
static inline void v4_store(const V4& x, double* d) { _mm256_storeu_pd(d, x.v); }
static inline DD4 dd4_load(const double* h, const double* l) { return DD4(_mm256_loadu_pd(h), _mm256_loadu_pd(l)); }
static inline void dd4_store(const DD4& x, double* h, double* l) {
    _mm256_storeu_pd(h, x.hi); _mm256_storeu_pd(l, x.lo);
}

#else
inline constexpr const char* simd_name = "generic";
struct V4 {
    double v[4];
    V4() : v{0, 0, 0, 0} {}
    V4(double x) : v{x, x, x, x} {}
};
#define KBSRC_L4(expr) for (int i = 0; i < 4; i++) { expr; }
static inline V4 operator+(V4 a, V4 b) { V4 r; KBSRC_L4(r.v[i] = a.v[i] + b.v[i]) return r; }
static inline V4 operator-(V4 a, V4 b) { V4 r; KBSRC_L4(r.v[i] = a.v[i] - b.v[i]) return r; }
static inline V4 operator*(V4 a, V4 b) { V4 r; KBSRC_L4(r.v[i] = a.v[i] * b.v[i]) return r; }
static inline V4 operator-(V4 a) { V4 r; KBSRC_L4(r.v[i] = -a.v[i]) return r; }

struct DD4 {
    double hi[4], lo[4];
    DD4() : hi{0, 0, 0, 0}, lo{0, 0, 0, 0} {}
    DD4(double x) : hi{x, x, x, x}, lo{0, 0, 0, 0} {}
    DD4(const dd& x) : hi{x.hi, x.hi, x.hi, x.hi}, lo{x.lo, x.lo, x.lo, x.lo} {}
};
static inline DD4 operator+(const DD4& a, const DD4& b) {
    DD4 r;
    KBSRC_L4(double s = a.hi[i] + b.hi[i]; double bb = s - a.hi[i];
             double e = (a.hi[i] - (s - bb)) + (b.hi[i] - bb); e = e + (a.lo[i] + b.lo[i]);
             double t = s + e; r.hi[i] = t; r.lo[i] = e - (t - s))
    return r;
}
static inline DD4 operator-(const DD4& a) { DD4 r; KBSRC_L4(r.hi[i] = -a.hi[i]; r.lo[i] = -a.lo[i]) return r; }
static inline DD4 operator-(const DD4& a, const DD4& b) { return a + (-b); }
// Product of all four lanes with prod_err (tiny-product path included). Out of line: never needed with FMA,
// and only for subnormal-range products without it.
[[gnu::noinline, gnu::cold]] static inline DD4 dd4_mul_exact(const DD4& a, const DD4& b) {
    DD4 r;
    KBSRC_L4(double p = a.hi[i] * b.hi[i]; double e = prod_err(a.hi[i], b.hi[i], p);
             e = e + (a.hi[i] * b.lo[i] + a.lo[i] * b.hi[i]);
             double t = p + e; r.hi[i] = t; r.lo[i] = e - (t - p))
    return r;
}
// Cheap screen: true when some lane's |hi * hi| is below the prod_err_tiny threshold or NaN. !(m >= T) is
// also true for NaN, so a NaN lane cannot hide a tiny one in the pairwise minimum.
static inline bool dd4_small(const DD4& a, const DD4& b) {
    double m0 = std::fabs(a.hi[0] * b.hi[0]), m2 = std::fabs(a.hi[2] * b.hi[2]);
    double m1 = std::fabs(a.hi[1] * b.hi[1]), m3 = std::fabs(a.hi[3] * b.hi[3]);
    m0 = m2 < m0 ? m2 : m0; m1 = m3 < m1 ? m3 : m1;
    return !(m0 >= 0x1p-968) | !(m1 >= 0x1p-968);
}
// Exact test behind the screen: does some lane need prod_err's tiny-product path?
static inline bool dd4_tiny(const DD4& a, const DD4& b) {
    int tiny = 0;
    KBSRC_L4(tiny |= prod_err_tiny(a.hi[i], b.hi[i], a.hi[i] * b.hi[i]))
    return tiny;
}
// The hot loop uses the branch-free prod_err_fast (vectorised); a product with a tiny lane is redone exactly.
static inline DD4 operator*(const DD4& a, const DD4& b) {
    DD4 r;
    KBSRC_L4(double p = a.hi[i] * b.hi[i]; double e = prod_err_fast(a.hi[i], b.hi[i], p);
             e = e + (a.hi[i] * b.lo[i] + a.lo[i] * b.hi[i]);
             double t = p + e; r.hi[i] = t; r.lo[i] = e - (t - p))
    if (dd4_small(a, b) && dd4_tiny(a, b)) [[unlikely]] r = dd4_mul_exact(a, b);
    return r;
}
#undef KBSRC_L4
static inline V4 v4_load(const double* d) { V4 r; for (int i = 0; i < 4; i++) r.v[i] = d[i]; return r; }
static inline void v4_store(const V4& x, double* d) { for (int i = 0; i < 4; i++) d[i] = x.v[i]; }
static inline DD4 dd4_load(const double* h, const double* l) {
    DD4 r; for (int i = 0; i < 4; i++) { r.hi[i] = h[i]; r.lo[i] = l[i]; }
    return r;
}
static inline void dd4_store(const DD4& x, double* h, double* l) {
    for (int i = 0; i < 4; i++) { h[i] = x.hi[i]; l[i] = x.lo[i]; }
}
#endif

}  // namespace kbsrc::KBSRC_ISA::dsp
