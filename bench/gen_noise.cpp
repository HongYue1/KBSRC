// Writes the timing input: 180 s of stereo pink-ish noise at 96 kHz, -20 dBFS RMS, 64-bit float.
// Fixed seed. build_tools.sh compiles it with -ffp-contract=off, so the -march flags cannot change the
// samples through FMA contraction (std::normal_distribution and std::log still depend on the C++ and
// C libraries). A busy full-band signal avoids the denormal slowdowns that long digital silence causes
// in some time-domain resamplers.
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "io/wav.hpp"

int main(int argc, char** argv) {
    if (argc != 2) { std::fprintf(stderr, "usage: gen_noise out.wav\n"); return 2; }
    const int rate = 96000, ch = 2; const long long n = 180LL * rate;
    std::vector<std::vector<double>> y(ch, std::vector<double>((size_t)n));
    std::mt19937_64 rng(20261012); std::normal_distribution<double> g(0.0, 1.0);
    for (int c = 0; c < ch; c++) {
        double b0 = 0, b1 = 0, b2 = 0, ss = 0;   // Paul Kellet's economy pink filter
        for (long long i = 0; i < n; i++) {
            double w = g(rng);
            b0 = 0.99765 * b0 + w * 0.0990460; b1 = 0.96300 * b1 + w * 0.2965164; b2 = 0.57000 * b2 + w * 1.0526913;
            double v = b0 + b1 + b2 + w * 0.1848; y[c][i] = v; ss += v * v;
        }
        double k = 0.1 / std::sqrt(ss / n);
        for (double& v : y[c]) v *= k;
    }
    long long clipped = 0;
    kbsrc::io::write_wav(argv[1], y, ch, n, rate, "f64", false, clipped);
    return 0;
}
