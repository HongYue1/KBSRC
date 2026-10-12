// Mixed-radix (2,3,4,5,7) Stockham FFT, four-step decomposition with SIMD lanes, Bluestein for other lengths
#pragma once
#include <vector>
#include <stdexcept>
#include <type_traits>
#include <algorithm>
#include <memory>
#include <cstdint>
#include "dsp/double_double.hpp"
#include "dsp/simd.hpp"

namespace kbsrc::dsp {

template <class T> static inline T from_dd(const dd& v);
template <> inline double from_dd<double>(const dd& v) { return v.hi; }
template <> inline long double from_dd<long double>(const dd& v) { return (long double)v.hi + (long double)v.lo; }
template <> inline dd from_dd<dd>(const dd& v) { return v; }
static inline double to_double(double v) { return v; }
static inline double to_double(long double v) { return (double)v; }
static inline double to_double(const dd& v) { return v.hi + v.lo; }
static inline dd to_dd(double v) { return dd(v); }
static inline dd to_dd(long double v) { double h = (double)v; return dd(h, (double)(v - (long double)h)); }
static inline dd to_dd(const dd& v) { return v; }

template <class T> struct cx { T re, im; };

// fast exact-ish twiddles: cis(2 pi k/N) = cis(2 pi jB/N) * cis(2 pi l/N), tables of ~2 sqrt(N) Taylor values
struct CisTable {
    long long N = 0, B = 1; std::vector<dd> ac, as, bc, bs;
    explicit CisTable(long long n) : N(n) {
        B = (long long)std::ceil(std::sqrt((double)n)); long long J = n / B + 1;
        ac.resize(J); as.resize(J); bc.resize(B); bs.resize(B);
        for (long long j = 0; j < J; j++) dd_cis(j * B, n, ac[j], as[j]);
        for (long long l = 0; l < B; l++) dd_cis(l, n, bc[l], bs[l]);
    }
    void get(long long k, dd& c, dd& s) const {
        k %= N; if (k < 0) k += N; long long j = k / B, l = k % B;
        if (l == 0) { c = ac[j]; s = as[j]; return; }
        if (j == 0) { c = bc[l]; s = bs[l]; return; }
        c = ac[j] * bc[l] - as[j] * bs[l]; s = ac[j] * bs[l] + as[j] * bc[l];
    }
};
template <class T> static inline cx<T> operator+(const cx<T>& a, const cx<T>& b) { return {a.re + b.re, a.im + b.im}; }
template <class T> static inline cx<T> operator-(const cx<T>& a, const cx<T>& b) { return {a.re - b.re, a.im - b.im}; }
template <class T> static inline cx<T> operator*(const cx<T>& a, const cx<T>& b) {
    return {a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re};
}
template <class T, class W> static inline cx<T> mulw(const cx<T>& a, const cx<W>& b) {
    T br(b.re), bi(b.im); return {a.re * br - a.im * bi, a.re * bi + a.im * br};
}
template <class T> static inline cx<T> mul_mi(const cx<T>& a) { return {a.im, -a.re}; }  // a * (-i)
template <class T> static inline cx<T> conj(const cx<T>& a) { return {a.re, -a.im}; }

static inline bool smooth7(long long n) {
    if (n < 1) return false;
    for (int p : {2, 3, 5, 7}) while (n % p == 0) n /= p;
    return n == 1;
}

template <class T, class W, int R> static inline void odd_bfly(const cx<T>* src, size_t sm, cx<T>* dst, int s, int q, const W* cr, const W* sr, const cx<W>* tw, bool p0) {
    constexpr int H = (R - 1) / 2;
    cx<T> a0 = src[q], S[H], D[H];
    for (int t = 1; t <= H; t++) { cx<T> x = src[q + t * sm], y = src[q + (R - t) * sm]; S[t - 1] = x + y; D[t - 1] = x - y; }
    cx<T> sum = a0; for (int t = 0; t < H; t++) sum = sum + S[t];
    dst[q] = sum;
    for (int u = 1; u <= H; u++) {
        cx<T> A = a0, B = {T(0), T(0)};
        for (int t = 1; t <= H; t++) {
            int m = (t * u) % R; T c(cr[m]), sn(sr[m]);
            A = {A.re + S[t - 1].re * c, A.im + S[t - 1].im * c};
            B = {B.re + D[t - 1].re * sn, B.im + D[t - 1].im * sn};
        }
        cx<T> xu = {A.re + B.im, A.im - B.re}, xr = {A.re - B.im, A.im + B.re};   // A -/+ iB
        if (p0) { dst[q + u * s] = xu; dst[q + (R - u) * s] = xr; }
        else { dst[q + u * s] = mulw(xu, tw[u - 1]); dst[q + (R - u) * s] = mulw(xr, tw[R - u - 1]); }
    }
}

template <class T, class W = T> struct SFFT {  // Stockham forward complex DFT, exponent sign -1; data T, twiddles W
    int n = 0;
    struct Stage { int r, m, s; std::vector<cx<W>> tw; };
    std::vector<Stage> st;
    W cr[3][7], sr[3][7];
    void init(int n_) {
        n = n_; st.clear();
        std::vector<int> rs; int t = n;
        while (t % 4 == 0) { rs.push_back(4); t /= 4; }
        while (t % 2 == 0) { rs.push_back(2); t /= 2; }
        for (int p : {3, 5, 7}) while (t % p == 0) { rs.push_back(p); t /= p; }
        if (t != 1) throw std::runtime_error("FFT length not 7-smooth");
        int ns = n, s = 1; CisTable ct(n);
        for (int r : rs) {
            Stage g; g.r = r; g.m = ns / r; g.s = s; g.tw.resize((size_t)g.m * (r - 1));
            for (int p = 0; p < g.m; p++)
                for (int u = 1; u < r; u++) {
                    dd c, sn; ct.get(-(long long)s * p * u, c, sn);
                    g.tw[(size_t)p * (r - 1) + u - 1] = {from_dd<W>(c), from_dd<W>(sn)};
                }
            st.push_back(std::move(g)); ns /= r; s *= r;
        }
        int rr[3] = {3, 5, 7};
        for (int i = 0; i < 3; i++) for (int k = 0; k < rr[i]; k++) { dd c, sn; dd_cis(k, rr[i], c, sn); cr[i][k] = from_dd<W>(c); sr[i][k] = from_dd<W>(sn); }
    }
    // in: x (n), work: y (n). returns pointer to the result (x or y).
    cx<T>* run(cx<T>* x, cx<T>* y) const {
        cx<T>* a = x; cx<T>* b = y;
        for (const Stage& g : st) {
            const int r = g.r, m = g.m, s = g.s;
            for (int p = 0; p < m; p++) {
                const cx<W>* tw = &g.tw[(size_t)p * (r - 1)];
                const cx<T>* src = a + (size_t)s * p; cx<T>* dst = b + (size_t)s * r * p;
                const size_t sm = (size_t)s * m;
                if (r == 4) {
                    for (int q = 0; q < s; q++) {
                        cx<T> a0 = src[q], a1 = src[q + sm], a2 = src[q + 2 * sm], a3 = src[q + 3 * sm];
                        cx<T> t0 = a0 + a2, t1 = a0 - a2, t2 = a1 + a3, t3 = mul_mi(a1 - a3);
                        dst[q] = t0 + t2;
                        if (p == 0) { dst[q + s] = t1 + t3; dst[q + 2 * s] = t0 - t2; dst[q + 3 * s] = t1 - t3; }
                        else { dst[q + s] = mulw(t1 + t3, tw[0]); dst[q + 2 * s] = mulw(t0 - t2, tw[1]); dst[q + 3 * s] = mulw(t1 - t3, tw[2]); }
                    }
                } else if (r == 2) {
                    for (int q = 0; q < s; q++) {
                        cx<T> a0 = src[q], a1 = src[q + sm];
                        dst[q] = a0 + a1; dst[q + s] = p ? mulw(a0 - a1, tw[0]) : a0 - a1;
                    }
                } else if (r == 3) { for (int q = 0; q < s; q++) odd_bfly<T, W, 3>(src, sm, dst, s, q, cr[0], sr[0], tw, p == 0); }
                else if (r == 5) { for (int q = 0; q < s; q++) odd_bfly<T, W, 5>(src, sm, dst, s, q, cr[1], sr[1], tw, p == 0); }
                else { for (int q = 0; q < s; q++) odd_bfly<T, W, 7>(src, sm, dst, s, q, cr[2], sr[2], tw, p == 0); }
            }
            std::swap(a, b);
        }
        return a;
    }
};

// lane packing between scalar storage S and SIMD compute type T (4 independent sub-FFTs per vector)
template <class S, class T> struct Lane;
template <class S> struct Lane<S, S> {
    static constexpr int N = 1;
    static cx<S> pack(const cx<S>* const* p, size_t) { return *p[0]; }
    static void unpack(const cx<S>& v, cx<S>* const* o, int) { *o[0] = v; }
    static cx<S> packw(const cx<S>* const* p) { return *p[0]; }
};
template <> struct Lane<double, V4> {
    static constexpr int N = 4;
    static inline cx<V4> pack4(const cx<double>& a, const cx<double>& b, const cx<double>& c, const cx<double>& d) {
        __m256d l0 = _mm256_loadu2_m128d((const double*)&c, (const double*)&a);   // a.re a.im c.re c.im
        __m256d l1 = _mm256_loadu2_m128d((const double*)&d, (const double*)&b);   // b.re b.im d.re d.im
        return {V4(_mm256_unpacklo_pd(l0, l1)), V4(_mm256_unpackhi_pd(l0, l1))};   // a b c d order
    }
    static cx<V4> packw(const cx<double>* const* p) { return pack4(*p[0], *p[1], *p[2], *p[3]); }
    static void unpack(const cx<V4>& v, cx<double>* const* o, int n) {
        __m256d lo = _mm256_unpacklo_pd(v.re.v, v.im.v), hi = _mm256_unpackhi_pd(v.re.v, v.im.v);  // (a,c) (b,d)
        alignas(32) double t0[4], t1[4]; _mm256_store_pd(t0, lo); _mm256_store_pd(t1, hi);
        o[0]->re = t0[0]; o[0]->im = t0[1]; if (n > 1) { o[1]->re = t1[0]; o[1]->im = t1[1]; }
        if (n > 2) { o[2]->re = t0[2]; o[2]->im = t0[3]; } if (n > 3) { o[3]->re = t1[2]; o[3]->im = t1[3]; }
    }
};
template <> struct Lane<dd, DD4> {
    static constexpr int N = 4;
    static cx<DD4> packw(const cx<dd>* const* p) {
        return {DD4(_mm256_set_pd(p[3]->re.hi, p[2]->re.hi, p[1]->re.hi, p[0]->re.hi), _mm256_set_pd(p[3]->re.lo, p[2]->re.lo, p[1]->re.lo, p[0]->re.lo)),
                DD4(_mm256_set_pd(p[3]->im.hi, p[2]->im.hi, p[1]->im.hi, p[0]->im.hi), _mm256_set_pd(p[3]->im.lo, p[2]->im.lo, p[1]->im.lo, p[0]->im.lo))};
    }
    static void unpack(const cx<DD4>& v, cx<dd>* const* o, int n) {
        alignas(32) double a[4], b[4], c[4], d[4];
        _mm256_store_pd(a, v.re.hi); _mm256_store_pd(b, v.re.lo); _mm256_store_pd(c, v.im.hi); _mm256_store_pd(d, v.im.lo);
        for (int l = 0; l < n; l++) { o[l]->re = qts(a[l], b[l]); o[l]->im = qts(c[l], d[l]); }
    }
};
template <class S, class T> static inline cx<T> lane_pack(const cx<S>* const* p) {
    if constexpr (std::is_same_v<S, T>) return *p[0];
    else if constexpr (std::is_same_v<T, V4>) return Lane<double, V4>::pack4(*p[0], *p[1], *p[2], *p[3]);
    else return Lane<S, T>::packw(p);
}

// Four-step (Bailey) FFT: n = N1*N2. Column FFTs (length N1) then twiddle, row FFTs (length N2) written
// transposed into natural order. Sub-FFTs run in SIMD type T on 4 columns / rows at once (L1/L2 resident),
// data is stored in scalar type S -> 2 streaming passes over memory, 4x less traffic than block lanes.
static inline long long next_smooth7(long long n) { while (!smooth7(n)) n++; return n; }
template <class S, class W = S, class T = S> struct CFFT {
    using LN = Lane<S, T>; static constexpr int NL = LN::N;
    int n = 0, N1 = 0, N2 = 0; bool four = false;
    SFFT<T, W> s1, s2; SFFT<S, W> sd; std::vector<cx<W>> tw;   // tw[k1*N2+n2] = w_n^{k1 n2}
    // Bluestein (chirp-z) for lengths that are not 7-smooth: b_t = e^{+i pi t^2/n} (exact, t^2 mod 2n)
    int m = 0; std::shared_ptr<CFFT> inner; std::vector<cx<W>> bc, Vf;
    size_t scratch_bytes() const {
        if (m) return (size_t)2 * m * sizeof(cx<S>) + inner->scratch_bytes() + 128;
        return four ? (size_t)2 * std::max(N1, N2) * sizeof(cx<T>) + 64 : 0;
    }
    void init(int n_) {
        n = n_; m = 0; inner.reset();
        if (!smooth7(n)) {
            m = (int)next_smooth7(2LL * n - 1); inner = std::make_shared<CFFT>(); inner->init(m);
            bc.resize(n); std::vector<cx<dd>> v(m, cx<dd>{dd(0.0), dd(0.0)}), w(m);
            for (long long t = 0; t < n; t++) { long long e = (long long)((__int128)t * t % (2LL * n)); dd c, sn; dd_cis(e, 2LL * n, c, sn); bc[t] = {from_dd<W>(c), from_dd<W>(sn)};
                v[t] = {c, sn}; if (t) v[m - t] = {c, sn}; }
            CFFT<dd, dd, dd> fd; fd.init(m); std::vector<__m256d> sc(fd.scratch_bytes() / 32 + 2);
            cx<dd>* V = fd.run(v.data(), w.data(), sc.data()); Vf.resize(m);
            for (int k = 0; k < m; k++) Vf[k] = {from_dd<W>(V[k].re / (double)m), from_dd<W>(V[k].im / (double)m)};
            return;
        }
        four = n > 2048;
        if (!four) { sd.init(n); return; }
        int best = 1; for (int d = 1; (long long)d * d <= n; d++) if (n % d == 0) best = d;
        N1 = best; N2 = n / best; if (N1 < 8) { four = false; sd.init(n); return; }
        s1.init(N1); s2.init(N2); tw.resize((size_t)n);
        CisTable ct(n);
        for (int k1 = 0; k1 < N1; k1++) for (int n2 = 0; n2 < N2; n2++) { dd c, sn; ct.get(-(long long)k1 * n2, c, sn); tw[(size_t)k1 * N2 + n2] = {from_dd<W>(c), from_dd<W>(sn)}; }
    }
    // x: input (n, destroyed), y: n, scratch: scratch_bytes(). returns result pointer (x or y).
    cx<S>* run(cx<S>* x, cx<S>* y, void* scratch) const {
        if (m) {   // X_k = conj(b_k) * sum_j (x_j conj(b_j)) b_{k-j}
            cx<S>* u = (cx<S>*)scratch; cx<S>* u2 = u + m; void* isc = (void*)(((uintptr_t)(u2 + m) + 63) & ~(uintptr_t)63);
            for (int j = 0; j < n; j++) u[j] = mulw(x[j], conj(bc[j]));
            for (int j = n; j < m; j++) u[j] = {S(0.0), S(0.0)};
            cx<S>* U = inner->run(u, u2, isc);
            for (int k = 0; k < m; k++) U[k] = conj(mulw(U[k], Vf[k]));          // conj -> inverse via forward
            cx<S>* other = U == u ? u2 : u; cx<S>* C = inner->run(U, other, isc);
            for (int k = 0; k < n; k++) y[k] = mulw(conj(C[k]), conj(bc[k]));
            return y;
        }
        if (!four) return sd.run(x, y);
        const int M = std::max(N1, N2);
        cx<T>* tile = (cx<T>*)scratch; cx<T>* wk = tile + M;
        const cx<S>* ip[4]; cx<S>* op[4]; const cx<W>* wp[4];
        for (int c0 = 0; c0 < N2; c0 += NL) {
            int g = std::min(NL, N2 - c0);
            for (int r = 0; r < N1; r++) { for (int j = 0; j < NL; j++) ip[j] = x + (size_t)r * N2 + c0 + std::min(j, g - 1); tile[r] = lane_pack<S, T>(ip); }
            cx<T>* res = s1.run(tile, wk);
            for (int k1 = 0; k1 < N1; k1++) {
                cx<T> v = res[k1];
                if (k1) {
                    if constexpr (NL == 1) { if (c0) v = mulw(v, tw[(size_t)k1 * N2 + c0]); }
                    else { for (int j = 0; j < NL; j++) wp[j] = &tw[(size_t)k1 * N2 + c0 + std::min(j, g - 1)]; v = v * lane_pack<W, T>(wp); }
                }
                for (int j = 0; j < g; j++) op[j] = x + (size_t)k1 * N2 + c0 + j;
                LN::unpack(v, op, g);
            }
        }
        for (int r0 = 0; r0 < N1; r0 += NL) {
            int g = std::min(NL, N1 - r0);
            for (int i = 0; i < N2; i++) { for (int j = 0; j < NL; j++) ip[j] = x + (size_t)(r0 + std::min(j, g - 1)) * N2 + i; tile[i] = lane_pack<S, T>(ip); }
            cx<T>* res = s2.run(tile, wk);
            for (int k2 = 0; k2 < N2; k2++) { for (int j = 0; j < g; j++) op[j] = y + (size_t)N1 * k2 + r0 + j; LN::unpack(res[k2], op, g); }
        }
        return y;
    }
};

template <class T, class W = T, class C = T> struct RFFT {  // real transform of even length N via complex N/2 (C: SIMD compute type)
    int N = 0, h = 0; CFFT<T, W, C> c; std::vector<cx<W>> w;  // w[k] = e^{-2 pi i k/N}
    void init(int N_) {
        N = N_; h = N / 2; if (N % 2) throw std::runtime_error("RFFT length must be even");
        c.init(h); w.resize(h + 1); CisTable ct(N);
        for (int k = 0; k <= h; k++) { dd cc, s; ct.get(-k, cc, s); w[k] = {from_dd<W>(cc), from_dd<W>(s)}; }
    }
    // r: N reals -> X: h+1 bins (unnormalized). z, y: h scratch each.
    size_t scratch_bytes() const { return c.scratch_bytes(); }
    void forward(const T* r, cx<T>* X, cx<T>* z, cx<T>* y, void* sc) const {
        for (int i = 0; i < h; i++) z[i] = {r[2 * i], r[2 * i + 1]};
        z = c.run(z, y, sc);
        X[0] = {z[0].re + z[0].im, T(0.0)}; X[h] = {z[0].re - z[0].im, T(0.0)};
        for (int k = 1; k < h; k++) {
            cx<T> A = z[k], B = conj(z[h - k]);
            cx<T> E = A + B, O = mul_mi(A - B);   // 2E, 2O
            cx<T> v = E + mulw(O, w[k]);
            X[k] = {v.re * T(0.5), v.im * T(0.5)};
        }
    }
    // X: h+1 hermitian bins -> r: N reals, unnormalized (= N * true inverse)
    void inverse(const cx<T>* X, T* r, cx<T>* z, cx<T>* y, void* sc) const {
        for (int k = 0; k < h; k++) {
            cx<T> A = X[k], B = conj(X[h - k]);
            cx<T> E = A + B, D = mulw(A - B, conj(w[k]));
            cx<T> v = {E.re - D.im, E.im + D.re};   // E + i D
            z[k] = conj(v);
        }
        z = c.run(z, y, sc);
        for (int i = 0; i < h; i++) { r[2 * i] = z[i].re; r[2 * i + 1] = -z[i].im; }
    }
};

}  // namespace kbsrc::dsp
