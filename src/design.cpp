#include "design.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <thread>

#include "dsp/fft.hpp"
#include "parallel.hpp"

namespace kbsrc::KBSRC_ISA::detail {

using namespace dsp;

// one-sided taps of the Kaiser-windowed sinc sampled at fk, cutoff fc = (a/b)*fk/2 (a/b reduced).
// Everything is double-double built from IEEE basic operations, so the taps are identical on every platform.
std::vector<dd> kernel_taps(int L, long long a, long long b, double beta) {
    std::vector<dd> h(L + 1); const dd ib = dd_i0(dd(beta)), ab = dd((double)a) / (double)b;
    auto tap = [&](int n) {
        dd s(1.0);
        if (n) {   // sin(pi n a/b) / (pi n a/b), argument reduced exactly: n a mod 2b
            long long m = (long long)(((__int128)n * a) % (2 * (__int128)b));
            dd c, sn; dd_cis(m, 2 * b, c, sn);
            s = sn / (DD_PI * (dd((double)n) * (double)a) / (double)b);
        }
        dd u = dd((double)n) / (double)L, w = dd_i0(dd(beta) * dd_sqrt(dd(1.0) - u * u)) / ib;
        h[n] = ab * s * w;
    };
    // taps are independent: split them over threads (the result does not depend on the split)
    int nt = (int)std::min<long long>(std::max(1u, std::thread::hardware_concurrency()), std::max(1, (L + 1) / 512));
    std::atomic<int> next(0);
    auto work = [&](const std::atomic<bool>& stop) {
        for (int c; !stop && (c = next.fetch_add(256)) <= L;) for (int n = c; n <= std::min(L, c + 255); n++) tap(n);
    };
    parallel::run(nt, work);
    return h;
}

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
