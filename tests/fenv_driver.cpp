// Test driver for tests/test_engines.py (not installed): runs kbsrc::resample with the caller in a given
// floating-point rounding mode and reports whether that mode survived the call.
//
// usage: kbsrc_fenv_driver <nearest|down|up|zero> <exact|fast|ld> <rate_in> <rate_out> <channels> <length>
//                          <beta> <threads> <in.f64> <out.f64>
//   in.f64   interleaved little-endian doubles (frames = size / 8 / channels)
//   out.f64  planar output, channel after channel
// Prints "isa=<engine> mode_kept=<0|1> threw=<0|1> frames=<n>". Exit 0; 3 when resample threw
// std::invalid_argument (the mode must survive that too); 2 on usage or I/O errors.
#include <cfenv>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>
#if defined(__x86_64__) || defined(__i386__)
#include <xmmintrin.h>
#endif

#include "kbsrc/kbsrc.hpp"

namespace {
// The rounding state as the hardware holds it: fegetround() plus, on x86, the MXCSR (SSE) and x87 rounding
// controls separately (MinGW's fegetround reports only one of them).
struct FpMode {
    int fe = -1;
    unsigned sse = 0, x87 = 0;
    static FpMode now() {
        FpMode m;
        m.fe = std::fegetround();
#if defined(__x86_64__) || defined(__i386__)
        unsigned short cw = 0;
        __asm__ volatile("fnstcw %0" : "=m"(cw));
        m.sse = _mm_getcsr() & 0x6000u;
        m.x87 = cw & 0x0C00u;
#endif
        return m;
    }
    bool operator==(const FpMode& o) const { return fe == o.fe && sse == o.sse && x87 == o.x87; }
};

int usage() {
    std::fprintf(stderr, "usage: kbsrc_fenv_driver <nearest|down|up|zero> <exact|fast|ld> <rate_in> <rate_out> "
                         "<channels> <length> <beta> <threads> <in.f64> <out.f64>\n");
    return 2;
}
}  // namespace

int main(int argc, char** argv) {
    if (argc != 11) return usage();
    // Parse everything before changing the rounding mode: strtod honours it.
    const std::string mode = argv[1], prec = argv[2];
    const int fe = mode == "nearest" ? FE_TONEAREST : mode == "down" ? FE_DOWNWARD
                 : mode == "up" ? FE_UPWARD : mode == "zero" ? FE_TOWARDZERO : -1;
    kbsrc::Options opt;
    if (prec == "exact") opt.precision = kbsrc::Precision::Exact;
    else if (prec == "fast") opt.precision = kbsrc::Precision::Fast;
    else if (prec == "ld") opt.precision = kbsrc::Precision::LongDouble;
    else return usage();
    const int rate_in = std::atoi(argv[3]), rate_out = std::atoi(argv[4]), channels = std::atoi(argv[5]);
    opt.length = std::strtod(argv[6], nullptr);
    opt.beta = std::strtod(argv[7], nullptr);
    opt.threads = std::atoi(argv[8]);
    if (fe < 0 || channels < 1) return usage();

    std::vector<double> in;
    if (FILE* f = std::fopen(argv[9], "rb")) {
        double buf[512];
        for (size_t n; (n = std::fread(buf, sizeof(double), 512, f)) > 0;) in.insert(in.end(), buf, buf + n);
        std::fclose(f);
    } else {
        std::perror(argv[9]);
        return 2;
    }
    const long long frames = (long long)(in.size() / channels);

    if (std::fesetround(fe) != 0) {
        std::fprintf(stderr, "cannot set rounding mode %s\n", argv[1]);
        return 2;
    }
    const FpMode before = FpMode::now();
    std::vector<std::vector<double>> out;
    kbsrc::Stats stats;
    bool threw = false;
    try {
        out = kbsrc::resample(in.data(), frames, channels, rate_in, rate_out, opt, &stats);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    const FpMode after = FpMode::now();
    std::fesetround(FE_TONEAREST);

    const bool kept = before == after && before.fe == fe;
    std::printf("isa=%s mode_kept=%d threw=%d frames=%zu\n", threw ? "none" : stats.isa, kept ? 1 : 0,
                threw ? 1 : 0, out.empty() ? (size_t)0 : out[0].size());
    FILE* f = std::fopen(argv[10], "wb");
    if (!f) {
        std::perror(argv[10]);
        return 2;
    }
    bool written = true;
    for (const auto& ch : out) written &= std::fwrite(ch.data(), sizeof(double), ch.size(), f) == ch.size();
    if (std::fclose(f) != 0 || !written) {
        std::fprintf(stderr, "%s: write error\n", argv[10]);
        return 2;
    }
    return threw ? 3 : 0;
}
