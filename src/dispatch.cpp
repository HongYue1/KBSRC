// Public entry points: picks the engine compiled for the best instruction set this CPU supports.
#include <cstdlib>
#include <cstring>

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

std::vector<std::vector<double>> resample(const double* in, long long frames, int channels, int rate_in, int rate_out,
                                          const Options& opt, Stats* stats) {
#ifdef KBSRC_ENGINE_AVX2
    const char* force = std::getenv("KBSRC_CPU");
    bool generic = force && std::strcmp(force, "generic") == 0;
    if (!generic && __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma"))
        return avx2::resample_impl(in, frames, channels, rate_in, rate_out, opt, stats);
#endif
    return base::resample_impl(in, frames, channels, rate_in, rate_out, opt, stats);
}

}  // namespace kbsrc
