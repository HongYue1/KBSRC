#include "design.hpp"

#include <algorithm>
#include <cmath>

#include "dsp/fft.hpp"

namespace kbsrc::KBSRC_ISA::detail {

using namespace dsp;

namespace {

using LD = long double;
LD i0l(LD x) { LD q = x * x / 4, s = 1, t = 1; for (int k = 1; k < 1000; k++) { t *= q / ((LD)k * k); s += t; if (t < s * 1e-22L) break; } return s; }
const LD PIL = 3.141592653589793238462643383279502884L;

}  // namespace

// one-sided taps of the Kaiser-windowed sinc sampled at fk, cutoff fc = (a/b)*fk/2 (a/b reduced)
std::vector<dd> kernel_taps(int L, long long a, long long b, double beta) {
    std::vector<dd> h(L + 1); LD ib = i0l(beta);
    for (int n = 0; n <= L; n++) {
        LD s;
        if (n == 0) s = 1;
        else { long long m = (long long)(((__int128)n * a) % (2 * (__int128)b)); s = sinl(PIL * (LD)m / (LD)b) / (PIL * (LD)n * (LD)a / (LD)b); }
        LD u = (LD)n / L, w = i0l((LD)beta * sqrtl(std::max<LD>(0, 1 - u * u))) / ib;
        LD v = (LD)a / (LD)b * s * w; double hi = (double)v; h[n] = dd(hi, (double)(v - (LD)hi));
    }
    return h;
}
std::vector<dd> kernel(long long fk, long long a, long long b, double T, double beta) { return kernel_taps((int)std::floor(T * (double)fk), a, b, beta); }

// spectrum of the circularly placed symmetric kernel, length Nlen, times 1/scale_den (dd)
std::vector<dd> kernel_spectrum(const std::vector<dd>& h, int Nlen, double scale_den) {
    RFFT<dd, dd, dd> R; R.init(Nlen); int Lk = (int)h.size() - 1;
    std::vector<dd> r(Nlen); r[0] = h[0];
    for (int n = 1; n <= Lk; n++) { r[n] = h[n]; r[Nlen - n] = h[n]; }
    std::vector<cx<dd>> X(Nlen / 2 + 1), z(Nlen / 2), y(Nlen / 2); std::vector<Scratch32> sc(R.scratch_bytes() / 32 + 2);
    R.forward(r.data(), X.data(), z.data(), y.data(), sc.data());
    std::vector<dd> H(Nlen / 2 + 1); for (int k = 0; k <= Nlen / 2; k++) H[k] = X[k].re / scale_den;
    return H;
}

}  // namespace kbsrc::KBSRC_ISA::detail
