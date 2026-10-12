// KBSRC - Kaiser-Bessel fold-back sample-rate converter
#pragma once

#include <vector>

namespace kbsrc {

enum class Precision {
    Fast,        // float64 pipeline, AVX2
    Exact,       // double-double (~106-bit) pipeline, AVX2
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
};

// Resamples interleaved input (frames x channels) from rate_in to rate_out.
// Returns planar output: one vector per channel, ceil(frames * rate_out / rate_in) samples each.
// Throws std::invalid_argument on bad arguments.
std::vector<std::vector<double>> resample(const double* interleaved, long long frames, int channels,
                                          int rate_in, int rate_out, const Options& options = {},
                                          Stats* stats = nullptr);

const char* version();

}  // namespace kbsrc
