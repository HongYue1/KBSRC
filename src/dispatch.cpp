// Public entry points: picks the engine compiled for the best instruction set this CPU supports.
#include <cstdlib>
#include <cstring>
#if defined(__SSE2_MATH__)
#include <xmmintrin.h>
#else
#include <cfenv>
#endif

#include "kbsrc/kbsrc.hpp"

#ifndef KBSRC_VERSION
#define KBSRC_VERSION "dev"
#endif

namespace kbsrc {

#define KBSRC_DECLARE_ENGINE(ns)                                                                              \
    namespace ns {                                                                                            \
    std::vector<std::vector<double>> resample_impl(const double*, long long, int, int, int, const Options&, \
                                                   Stats*);                                                   \
    }
KBSRC_DECLARE_ENGINE(base)
#ifdef KBSRC_ENGINE_AVX2
KBSRC_DECLARE_ENGINE(avx2)
#endif

const char* version() { return KBSRC_VERSION; }

namespace {
// Exact is computed in round-to-nearest whatever the caller's rounding mode. Its error-free transforms (TwoSum,
// the product error) are only exact there; under a directed mode the FMA and non-FMA engines would round
// differently, and so would worker threads, which may start in the default mode (Windows). The caller's mode
// is restored on return or throw. Only the rounding-control bits are touched: exception flags raised meanwhile
// stay raised, flush-to-zero/denormals-are-zero are left as they are (they must be off, see kbsrc.hpp), and so
// is the x87 control word (LongDouble, not used by Exact).
class NearestRounding {
public:
    explicit NearestRounding(bool active) {
        if (!active) return;
#if defined(__SSE2_MATH__)
        saved_ = _mm_getcsr() & kRoundingBits;
        if (saved_) _mm_setcsr(_mm_getcsr() & ~kRoundingBits);
#elif defined(FE_TONEAREST)
        int mode = std::fegetround();
        if (mode >= 0 && mode != FE_TONEAREST && std::fesetround(FE_TONEAREST) == 0) saved_ = mode;
#endif
    }
    ~NearestRounding() {
#if defined(__SSE2_MATH__)
        if (saved_) _mm_setcsr((_mm_getcsr() & ~kRoundingBits) | saved_);
#elif defined(FE_TONEAREST)
        if (saved_ != kNone) std::fesetround(saved_);
#endif
    }
    NearestRounding(const NearestRounding&) = delete;
    NearestRounding& operator=(const NearestRounding&) = delete;

private:
#if defined(__SSE2_MATH__)
    static constexpr unsigned kRoundingBits = 0x6000;  // MXCSR RC field; 0 = round-to-nearest
    unsigned saved_ = 0;                               // caller's RC bits when they had to be cleared
#else
    static constexpr int kNone = -1;
    int saved_ = kNone;                                // caller's mode when it had to be changed
#endif
};
}  // namespace

std::vector<std::vector<double>> resample(const double* in, long long frames, int channels, int rate_in, int rate_out,
                                          const Options& opt, Stats* stats) {
    // resample_impl runs the double-double (Exact) engine for every precision except Fast and LongDouble
    NearestRounding nearest(opt.precision != Precision::Fast && opt.precision != Precision::LongDouble);
#ifdef KBSRC_ENGINE_AVX2
    const char* force = std::getenv("KBSRC_CPU");
    bool generic = force && std::strcmp(force, "generic") == 0;
    if (!generic && __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma"))
        return avx2::resample_impl(in, frames, channels, rate_in, rate_out, opt, stats);
#endif
    return base::resample_impl(in, frames, channels, rate_in, rate_out, opt, stats);
}

}  // namespace kbsrc
