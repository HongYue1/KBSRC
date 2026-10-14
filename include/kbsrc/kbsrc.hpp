// KBSRC - Kaiser-Bessel fold-back sample-rate converter
#pragma once

#include <vector>

namespace kbsrc {

enum class Precision {
    Fast,        // float64 pipeline
    Exact,       // double-double (~106-bit) pipeline
    LongDouble,  // x87 80-bit reference
};

struct Options {
    Precision precision = Precision::Exact;
    double length = 0.1;    // one-sided kernel length, seconds
    double beta = 40.0;     // Kaiser window beta
    int threads = 0;        // 0 = automatic
    int block = 0;          // block size factor, 0 = automatic
    bool lpc_edges = true;  // Burg LPC extrapolation at the signal edges
};

struct Stats {
    long long p = 0, q = 0;     // rate ratio P/Q
    long long taps = 0;         // kernel taps at the kernel rate
    int block = 0, threads = 0;
    bool general_fold = false;  // near-unity ratio path
    double design_seconds = 0, process_seconds = 0;
    const char* isa = "";       // engine used: "avx2" or "generic" (identical output)
};

// Resamples interleaved input (frames x channels) from rate_in to rate_out.
// Returns planar output: one vector per channel, ceil(frames * rate_out / rate_in) samples each.
// Throws std::invalid_argument on bad arguments. Portable builds pick the AVX2 engine when the CPU
// supports it; setting the environment variable KBSRC_CPU=generic forces the generic one.
// Floating-point environment: Exact always computes in round-to-nearest, so its output is bit-identical on
// both engines and does not depend on the caller's rounding mode; the caller's mode is restored on return or
// throw. Flush-to-zero and denormals-are-zero must be off. Fast and LongDouble run in the caller's mode and
// are specified (and engine-identical) only for round-to-nearest.
std::vector<std::vector<double>> resample(const double* interleaved, long long frames, int channels,
                                          int rate_in, int rate_out, const Options& options = {},
                                          Stats* stats = nullptr);

const char* version();

}  // namespace kbsrc
