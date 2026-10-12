// 4-lane AVX2 types: V4 (4 x double) and DD4 (4 x double-double)
#pragma once
#include <immintrin.h>
#include "dsp/double_double.hpp"

namespace kbsrc::dsp {
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
static inline DD4 qts4(__m256d a, __m256d b) { __m256d s = _mm256_add_pd(a, b); return DD4(s, _mm256_sub_pd(b, _mm256_sub_pd(s, a))); }
static inline DD4 operator+(const DD4& a, const DD4& b) {
    __m256d s = _mm256_add_pd(a.hi, b.hi), bb = _mm256_sub_pd(s, a.hi);
    __m256d e = _mm256_add_pd(_mm256_sub_pd(a.hi, _mm256_sub_pd(s, bb)), _mm256_sub_pd(b.hi, bb));
    e = _mm256_add_pd(e, _mm256_add_pd(a.lo, b.lo));
    return qts4(s, e);
}
static inline DD4 operator-(const DD4& a) { __m256d m = _mm256_set1_pd(-0.0); return DD4(_mm256_xor_pd(a.hi, m), _mm256_xor_pd(a.lo, m)); }
static inline DD4 operator-(const DD4& a, const DD4& b) { return a + (-b); }
static inline DD4 operator*(const DD4& a, const DD4& b) {
    __m256d p = _mm256_mul_pd(a.hi, b.hi), e = _mm256_fmsub_pd(a.hi, b.hi, p);
    e = _mm256_add_pd(e, _mm256_add_pd(_mm256_mul_pd(a.hi, b.lo), _mm256_mul_pd(a.lo, b.hi)));
    return qts4(p, e);
}
// lane access
static inline double lane(const V4& x, int l) { alignas(32) double t[4]; _mm256_store_pd(t, x.v); return t[l]; }
static inline double lane(const DD4& x, int l) { alignas(32) double h[4], o[4]; _mm256_store_pd(h, x.hi); _mm256_store_pd(o, x.lo); return h[l] + o[l]; }
static inline void lanes(const V4& x, double* d) { _mm256_storeu_pd(d, x.v); }
static inline void lanes(const DD4& x, double* d) { _mm256_storeu_pd(d, _mm256_add_pd(x.hi, x.lo)); }
static inline V4 pack4(double a, double b, double c, double d, V4*) { return _mm256_set_pd(d, c, b, a); }
static inline DD4 pack4(double a, double b, double c, double d, DD4*) { return DD4(_mm256_set_pd(d, c, b, a), _mm256_setzero_pd()); }

}  // namespace kbsrc::dsp
