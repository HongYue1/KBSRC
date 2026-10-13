#include "io/wav.hpp"

#include <algorithm>
#include <cfenv>
#include <cmath>
#include <cstdlib>
#if defined(__x86_64__) || defined(__i386__)
#include <smmintrin.h>
#endif
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <random>
#include <stdexcept>

#ifdef _WIN32
#include <windows.h>
#endif

namespace kbsrc::io {
namespace {

[[noreturn]] void fail(const char* m) { throw std::runtime_error(m); }

#ifdef _WIN32
FILE* uopen(const char* p, const char* m) {   // UTF-8 path -> wide
    wchar_t wp[4096], wm[8]; MultiByteToWideChar(CP_UTF8, 0, p, -1, wp, 4096); MultiByteToWideChar(CP_UTF8, 0, m, -1, wm, 8);
    return _wfopen(wp, wm);
}
#define FSEEK _fseeki64
#define FTELL _ftelli64
#else
FILE* uopen(const char* p, const char* m) { return fopen(p, m); }
#define FSEEK fseeko
#define FTELL ftello
#endif

struct File {  // closes on scope exit so errors thrown mid-way do not leak the handle
    FILE* f;
    ~File() { if (f) fclose(f); }
};


using PackFn = void (*)(uint8_t*, const std::vector<std::vector<double>>&, int, long long, long long,
                       std::mt19937_64&, std::uniform_real_distribution<double>&, long long&);
template<int B, bool D>
void pack_scalar(uint8_t* p, const std::vector<std::vector<double>>& y, int ch, long long first, long long last,
                 std::mt19937_64& rng, std::uniform_real_distribution<double>& U, long long& clipped) {
    constexpr double sc = B == 2 ? 32768.0 : B == 3 ? 8388608.0 : 2147483648.0;
    for (long long i = first; i < last; ++i) for (int c = 0; c < ch; ++c, p += B) {
        double v = y[c][i] * sc;
        if constexpr (D) v += U(rng) + U(rng);
        v = std::nearbyint(v);
        if (v > sc-1) { v = sc-1; ++clipped; } else if (v < -sc) { v = -sc; ++clipped; }
        // The release serializes NaN as zero; avoid its undefined cast.
        int64_t iv = std::isnan(v) ? 0 : (int64_t)v;
        for (int t = 0; t < B; ++t) p[t] = (iv >> (8*t)) & 255;
    }
}

#if defined(__x86_64__) || defined(__i386__)
template<int B, bool D, int R>
__attribute__((target("sse4.1")))
void pack_sse(uint8_t* p, const std::vector<std::vector<double>>& y, int ch, long long first, long long last,
              std::mt19937_64& rng, std::uniform_real_distribution<double>& U, long long& clipped) {
    constexpr double sc = B == 2 ? 32768.0 : B == 3 ? 8388608.0 : 2147483648.0;
    for (long long i = first; i < last; ++i) for (int c = 0; c < ch; ++c, p += B) {
        double v = y[c][i] * sc;
        if constexpr (D) v += U(rng) + U(rng);
        v = _mm_cvtsd_f64(_mm_round_sd(_mm_setzero_pd(), _mm_set_sd(v), R | _MM_FROUND_NO_EXC));
        if (v > sc-1) { v = sc-1; ++clipped; } else if (v < -sc) { v = -sc; ++clipped; }
        // The release serializes NaN as zero; avoid its undefined cast.
        int64_t iv = std::isnan(v) ? 0 : (int64_t)v;
        for (int t = 0; t < B; ++t) p[t] = (iv >> (8*t)) & 255;
    }
}
#endif

template<int B, bool D>
PackFn select_pack(bool use_sse, int rounding) {
#if defined(__x86_64__) || defined(__i386__)
    if (use_sse) {
        switch (rounding) {
        case FE_TONEAREST: return &pack_sse<B,D,_MM_FROUND_TO_NEAREST_INT>;
        case FE_DOWNWARD: return &pack_sse<B,D,_MM_FROUND_TO_NEG_INF>;
        case FE_UPWARD: return &pack_sse<B,D,_MM_FROUND_TO_POS_INF>;
        case FE_TOWARDZERO: return &pack_sse<B,D,_MM_FROUND_TO_ZERO>;
        }
    }
#else
    (void)use_sse; (void)rounding;
#endif
    return &pack_scalar<B,D>;
}

}  // namespace

Audio read_wav(const char* path) {
    File fh{uopen(path, "rb")}; FILE* f = fh.f; if (!f) fail("cannot open input");
    FSEEK(f, 0, SEEK_END); long long sz = FTELL(f); FSEEK(f, 0, SEEK_SET);
    uint8_t hd[12]; if (sz < 12 || fread(hd, 1, 12, f) != 12 || memcmp(hd, "RIFF", 4) || memcmp(hd + 8, "WAVE", 4)) fail("not a RIFF/WAVE file");
    auto u16 = [](const uint8_t* p) { return (unsigned)p[0] | (unsigned)p[1] << 8; };
    auto u32 = [&](const uint8_t* p) { return (uint32_t)u16(p) | (uint32_t)u16(p + 2) << 16; };
    long long pos = 12, dpos = -1, dlen = 0; int tag = 0, ch = 0, rate = 0, bits = 0, align = 0;
    while (pos + 8 <= sz) {
        uint8_t c[8]; FSEEK(f, pos, SEEK_SET); if (fread(c, 1, 8, f) != 8) break;
        uint32_t cl = u32(c + 4); long long body = pos + 8;
        if (!memcmp(c, "fmt ", 4)) {
            uint8_t fm[40] = {}; size_t n = fread(fm, 1, std::min<uint32_t>(cl, 40), f); if (n < 16) fail("bad fmt chunk");
            tag = u16(fm); ch = u16(fm + 2); rate = u32(fm + 4); align = u16(fm + 12); bits = u16(fm + 14);
            if (tag == 0xFFFE && n >= 26) tag = u16(fm + 24);
        } else if (!memcmp(c, "data", 4)) {
            dpos = body; dlen = std::min<long long>(cl, sz - body);
            if (cl == 0xFFFFFFFFu || cl == 0) dlen = sz - body;
            break;
        }
        pos = body + cl + (cl & 1);
    }
    if (dpos < 0 || !ch || !rate) fail("missing fmt/data chunk");
    if (align <= 0 || align % ch) fail("invalid WAV block alignment");
    if (!(tag == 1 || (tag == 3 && (bits == 32 || bits == 64)))) fail("unsupported WAV format (need PCM or IEEE float)");
    int bps = align / ch; if (tag == 1 && (bps < 1 || bps > 4)) fail("unsupported PCM width");
    Audio a; a.rate = rate; a.ch = ch; a.fmt_bits = bits; a.fmt_float = tag == 3;
    a.frames = dlen / align; size_t n = (size_t)a.frames * ch; a.d.resize(n);
    FSEEK(f, dpos, SEEK_SET);
    const size_t CS = (size_t)1 << 20; std::vector<uint8_t> buf(CS * bps);
    for (size_t k0 = 0; k0 < n; k0 += CS) {
        size_t m = std::min(CS, n - k0); if (fread(buf.data(), (size_t)bps, m, f) != m) fail("read error");
        const uint8_t* data = buf.data(); double* o = a.d.data() + k0;
        if (tag == 3 && bits == 64) memcpy(o, data, 8 * m);
        else if (tag == 3) for (size_t k = 0; k < m; k++) { float v; memcpy(&v, data + 4 * k, 4); o[k] = v; }
        else if (bps == 1) for (size_t k = 0; k < m; k++) o[k] = ((int)data[k] - 128) / 128.0;
        else if (bps == 2) for (size_t k = 0; k < m; k++) o[k] = (int16_t)(data[2 * k] | data[2 * k + 1] << 8) / 32768.0;
        else if (bps == 3) for (size_t k = 0; k < m; k++) { const uint8_t* p = data + 3 * k; int32_t v = p[0] | p[1] << 8 | p[2] << 16; if (v & 0x800000) v -= 1 << 24; o[k] = v / 8388608.0; }
        else for (size_t k = 0; k < m; k++) { const uint8_t* p = data + 4 * k; o[k] = (int32_t)(p[0] | p[1] << 8 | p[2] << 16 | (uint32_t)p[3] << 24) / 2147483648.0; }
    }
    return a;
}

void write_wav(const char* path, const std::vector<std::vector<double>>& y, int ch, long long frames, int rate, const std::string& fmt, bool dither, long long& clipped) {
    bool fl = fmt == "f64" || fmt == "f32"; int bits = fmt == "f64" ? 64 : fmt == "f32" ? 32 : fmt == "s32" ? 32 : fmt == "s24" ? 24 : 16;
    int kind = fmt == "f64" ? 0 : fmt == "f32" ? 1 : 2;
    int bps = bits / 8; size_t n = (size_t)frames * ch, dl = n * bps;
    if (dl + 80 > 0xFFFFFFFFull) fail("output too large for WAV (4 GB limit)");
    std::vector<uint8_t> out; out.reserve(96);
    auto p16 = [&](unsigned v) { out.push_back(v & 255); out.push_back(v >> 8 & 255); };
    auto p32 = [&](uint32_t v) { p16(v & 0xFFFF); p16(v >> 16); };
    auto tag4 = [&](const char* s) { out.insert(out.end(), s, s + 4); };
    bool ext = ch > 2 || (!fl && bits > 16);
    unsigned fmtlen = ext ? 40 : (fl ? 18 : 16);
    uint32_t riff = 4 + 8 + fmtlen + (fl ? 12 : 0) + 8 + dl + (dl & 1);
    tag4("RIFF"); p32(riff); tag4("WAVE"); tag4("fmt "); p32(fmtlen);
    p16(ext ? 0xFFFE : fl ? 3 : 1); p16(ch); p32(rate); p32(rate * ch * bps); p16(ch * bps); p16(bits);
    if (fmtlen >= 18) p16(ext ? 22 : 0);
    if (ext) {
        p16(bits); p32(ch == 1 ? 4 : ch == 2 ? 3 : 0);
        static const uint8_t guid[14] = {0x00,0x00,0x00,0x00,0x10,0x00,0x80,0x00,0x00,0xAA,0x00,0x38,0x9B,0x71};
        p16(fl ? 3 : 1); out.insert(out.end(), guid, guid + 14);
    }
    if (fl) { tag4("fact"); p32(4); p32((uint32_t)frames); }
    tag4("data"); p32((uint32_t)dl);
    File fh{uopen(path, "wb")}; FILE* f = fh.f; if (!f) fail("cannot open output");
    if (fwrite(out.data(), 1, out.size(), f) != out.size()) fail("write error");
    std::mt19937_64 rng(12345); std::uniform_real_distribution<double> U(-0.5, 0.5);
    bool use_sse = false;
    int rounding = std::fegetround();
#if defined(_WIN32) && (defined(__x86_64__) || defined(__i386__))
    // MinGW nearbyint uses x87, whereas fegetround may report MXCSR.
    unsigned short cw;
    __asm__ volatile("fnstcw %0" : "=m"(cw));
    const int x87_modes[] = {FE_TONEAREST, FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO};
    rounding = x87_modes[(cw >> 10) & 3];
#endif
#if defined(__x86_64__) || defined(__i386__)
    const char* cpu = std::getenv("KBSRC_CPU");
    use_sse = __builtin_cpu_supports("sse4.1") && !(cpu && std::strcmp(cpu,"generic") == 0);
#if !defined(_WIN32)
    // If the FP units disagree, retain the platform library behavior.
    const int mxcsr_modes[] = {FE_TONEAREST, FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO};
    if (mxcsr_modes[(_mm_getcsr() >> 13) & 3] != rounding) use_sse = false;
#endif
#endif
    PackFn pack;
    if (bps == 2) pack = dither ? select_pack<2,true>(use_sse,rounding) : select_pack<2,false>(use_sse,rounding);
    else if (bps == 3) pack = dither ? select_pack<3,true>(use_sse,rounding) : select_pack<3,false>(use_sse,rounding);
    else pack = dither ? select_pack<4,true>(use_sse,rounding) : select_pack<4,false>(use_sse,rounding);
    clipped = 0;
    const long long CF = 1 << 16; std::vector<uint8_t> buf((size_t)CF * ch * bps);
    for (long long f0 = 0; f0 < frames; f0 += CF) {
        long long nf = std::min(CF, frames - f0); uint8_t* p = buf.data();
        if (kind == 2) {
            pack(buf.data(), y, ch, f0, f0+nf, rng, U, clipped);
            p += (size_t)nf * ch * bps;
        } else {
            for (long long i = f0; i < f0+nf; ++i) for (int c = 0; c < ch; ++c, p += bps) {
                double x = y[c][i];
                if (kind == 0) memcpy(p, &x, 8);
                else { float v = (float)x; memcpy(p, &v, 4); }
            }
        }
        size_t nb = (size_t)(p - buf.data()); if (fwrite(buf.data(), 1, nb, f) != nb) fail("write error");
    }
    if (dl & 1) fputc(0, f);
    fh.f = nullptr; if (fclose(f)) fail("write error");
}

}  // namespace kbsrc::io
