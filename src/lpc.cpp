#include "lpc.hpp"

#include <algorithm>
#include <cmath>

#include "dsp/double_double.hpp"

namespace kbsrc::KBSRC_ISA::detail {

using dsp::dd;

// Burg recursion in double-double: float64 coefficients break gapless output
std::vector<double> burg_extrap(const double* x, long long len, long long n, int order, int train) {
    std::vector<double> out(n, 0.0);
    long long m = std::min<long long>(len, train); const double* seg = x + len - m;
    if (m < 2 * order + 8) return out;
    bool nz = false; double mx = 0; for (long long i = 0; i < m; i++) { if (seg[i] != 0) nz = true; mx = std::max(mx, std::fabs(seg[i])); }
    if (!nz) return out;
    std::vector<dd> f(seg + 1, seg + m), bk(seg, seg + m - 1), a{dd(1.0)};
    for (int it = 0; it < order; it++) {
        dd num(0.0), den(0.0); for (size_t i = 0; i < f.size(); i++) { num += f[i] * bk[i]; den += f[i] * f[i] + bk[i] * bk[i]; }
        if (den.hi <= 0) break;
        dd k = -(num * 2.0) / den; a.push_back(dd(0.0)); std::vector<dd> ar(a.rbegin(), a.rend());
        for (size_t i = 0; i < a.size(); i++) a[i] += k * ar[i];
        std::vector<dd> f2(f.size() - 1), b2(f.size() - 1);
        for (size_t i = 0; i + 1 < f.size(); i++) { f2[i] = f[i + 1] + k * bk[i + 1]; b2[i] = bk[i] + k * f[i]; }
        f.swap(f2); bk.swap(b2);
    }
    int p = (int)a.size() - 1; std::vector<double> c(p);
    for (int i = 0; i < p; i++) c[i] = -(a[p - i].hi + a[p - i].lo);      // c = -a[1:][::-1]
    std::vector<double> buf(p + n); for (int i = 0; i < p; i++) buf[i] = seg[m - p + i];
    for (long long i = 0; i < n; i++) { double s = 0; for (int j = 0; j < p; j++) s += c[j] * buf[i + j]; buf[p + i] = s; }
    double lim = 4 * mx;
    for (long long i = 0; i < n; i++) if (!std::isfinite(buf[p + i]) || std::fabs(buf[p + i]) > lim) return out;
    for (long long i = 0; i < n; i++) out[i] = buf[p + i];
    return out;
}

}  // namespace kbsrc::KBSRC_ISA::detail
