#include "kbsrc/kbsrc.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <climits>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <thread>
#include <vector>

#include "design.hpp"
#include "dsp/fft.hpp"
#include "lpc.hpp"

// One copy of this file is compiled per instruction set (see isa.hpp); src/dispatch.cpp picks one at run time.
namespace kbsrc::KBSRC_ISA {

using namespace dsp;
using detail::burg_extrap;
using detail::kernel;
using detail::kernel_spectrum;
using detail::kernel_taps;

namespace {

using LD = long double;

double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

// Overlap-save block engine. T: arithmetic type, W: coefficient type, C: SIMD type used inside the FFT.
template <class T, class W, class C>
void run_engine(const std::vector<std::vector<double>>& xx, std::vector<std::vector<double>>& out, long long P,
                long long Q, bool up, int nb, int marg, int hop, int nblk, long long o, long long nout,
                const std::vector<dd>& Hdd, int threads, int Jg) {
    const int Nb = nb * Q, Mb = nb * P, K = Mb / 2 + 1;
    RFFT<T, W, C> Fw, Iv; Fw.init(Nb); Iv.init(Mb);
    std::vector<W> H(Hdd.size()); for (size_t k = 0; k < Hdd.size(); k++) H[k] = from_dd<W>(Hdd[k]);
    const int ch = (int)xx.size(); const long long ntask = (long long)ch * nblk;
    std::atomic<long long> next(0);
    auto work = [&]() {
        const int mx = std::max(Nb, Mb);
        std::vector<T> r(Mb); std::vector<cx<T>> X(Nb / 2 + 1), Y(K), z(mx / 2), y(mx / 2);
        std::vector<Scratch32> sc(std::max(Fw.scratch_bytes(), Iv.scratch_bytes()) / 32 + 2);
        for (;;) {
            long long t = next++; if (t >= ntask) break;
            int c = (int)(t / nblk); long long b = t % nblk;
            const double* s = xx[c].data() + b * hop * Q;
            // Pack f64 input directly into the complex FFT buffer; no expanded staging pass.
            Fw.forward(s, X.data(), z.data(), y.data(), sc.data());
            if (Jg) {   // general path: true (unaliased) kernel response G[j] on the common grid, fold every bin mod Mb
                for (int k = 0; k < K; k++) Y[k] = {T(0.0), T(0.0)};
                for (int j = 0; j <= Jg; j++) {
                    int jj = j % Nb; cx<T> x = jj <= Nb / 2 ? X[jj] : conj(X[Nb - jj]); T g(H[j]);
                    cx<T> v = {x.re * g, x.im * g}; int k = j % Mb;
                    if (j == 0) Y[0] = Y[0] + v;
                    else if (k == 0 || 2 * k == Mb) Y[k].re = Y[k].re + T(2.0) * v.re;
                    else if (2 * k < Mb) Y[k] = Y[k] + v;
                    else Y[Mb - k] = Y[Mb - k] + conj(v);
                }
                Y[0].im = T(0.0); Y[Mb / 2].im = T(0.0);
            } else if (!up) {
                for (int k = 0; k < K; k++) { T hk(H[k]); Y[k] = {X[k].re * hk, X[k].im * hk}; }
                int jmax = std::min(Nb / 2, Mb - 1);
                for (int j = K; j <= jmax; j++) { T hj(H[j]); cx<T> v = {X[j].re * hj, X[j].im * hj}; Y[Mb - j] = Y[Mb - j] + conj(v); }
                Y[Mb / 2] = {T(2.0) * X[Mb / 2].re * T(H[Mb / 2]), T(0.0)};
            } else {
                for (int k = 0; k < K; k++) {
                    int kk = k % Nb; cx<T> v = kk <= Nb / 2 ? X[kk] : conj(X[Nb - kk]); T hk(H[k]);
                    Y[k] = {v.re * hk, v.im * hk};
                }
                Y[Mb / 2].im = T(0.0);
            }
            Iv.inverse(Y.data(), r.data(), z.data(), y.data(), sc.data());
            double* dst = out[c].data();
            for (long long i = 0; i < (long long)hop * P; i++) {
                long long d = b * hop * P + i - o;
                if (d >= 0 && d < nout) dst[d] = to_double(r[marg * P + i]);
            }
        }
    };
    std::vector<std::thread> th; for (int i = 1; i < threads; i++) th.emplace_back(work);
    work(); for (auto& t : th) t.join();
}

}  // namespace

std::vector<std::vector<double>> resample_impl(const double* in, long long frames, int channels, int rate_in, int rate_out,
                                          const Options& opt, Stats* stats) {
    if (channels < 1 || frames < 0 || (frames > 0 && !in)) throw std::invalid_argument("bad input buffer");
    if (rate_in < 1 || rate_out < 1) throw std::invalid_argument("sample rates must be positive");
    if (!(opt.length > 0 && opt.length <= 10) || !(opt.beta > 0 && opt.beta <= 200) || opt.threads < 0 || opt.block < 0)
        throw std::invalid_argument("option out of range");
    const int fo = rate_out;
    long long fi = rate_in, g = std::gcd(fi, (long long)fo), P = fo / g, Q = fi / g; bool up = fo > fi;
    long long nin = frames, nout = (nin * P + Q - 1) / Q; const int ch = channels;
    std::vector<std::vector<double>> out(ch, std::vector<double>(nout));
    if (fo == fi) { for (int c = 0; c < ch; c++) for (long long i = 0; i < nin; i++) out[c][i] = in[i * ch + c]; if (stats) { stats->p = P; stats->q = Q; } return out; }
    const double t0 = now();
    long long fk = up ? fo : fi;                 // kernel sampling rate; cutoff = min(fi,fo)/2
    long long mn = std::min(fi, (long long)fo), gg = std::gcd(mn, fk);
    // Near-unity ratios need a finer-grid kernel below. Compute the original reach
    // without designing taps that would never be used in that path.
    // Lk = one-sided kernel length in kernel samples. A length shorter than one kernel sample would collapse the
    // kernel (L = 0) and return nonfinite output as success; it is rejected, as is one too long for an int.
    const double lk = std::floor(opt.length * (double)fk);
    if (!(lk >= 1 && lk <= INT_MAX)) throw std::invalid_argument("option out of range");
    long long Lk = (long long)lk;
    long long L = (Lk * fi + fk - 1) / fk + 1;   // reach in input samples
    long long margl = std::max((L + Q + Q - 1) / Q, up ? (Lk + P + P - 1) / P : 0LL);
    // Transform sizes block * Q and block * P are ints, and the finer kernel grid needs 5/4 of the larger one:
    // refuse what would overflow them (huge explicit blocks, extreme rates or lengths) instead of misbehaving.
    if ((opt.block ? opt.block : 10 * margl) > INT_MAX / 2 / std::max(P, Q))
        throw std::invalid_argument("block, rates or length too large");
    // default thread cap: the pipeline is memory-bandwidth bound beyond a few cores
    int threads = opt.threads ? opt.threads : (int)std::min(std::max(1u, std::thread::hardware_concurrency()), opt.precision == Precision::Fast ? 4u : 8u);
    long long E = L + Q, pre = (E + 64 + Q - 1) / Q * Q, span = pre + nin + E + 64;
    auto fix = [&](int v) { while (v % 2 || !smooth7(v / 2)) v++; return v; };   // non-7-smooth P/Q -> Bluestein
    int nb = fix(opt.block ? opt.block : (int)(10 * margl));
    // Each block must leave a positive hop after its two margins. Automatic blocks always do (10 * marg, never
    // shrunk below 4 * marg); an explicit block at or below 2 * marg used to crash the process.
    if (nb <= 2 * margl) throw std::invalid_argument("option out of range");
    const int marg = (int)margl;
    // Short inputs: shrink blocks so up to 8 threads get work. Block size affects rounding, so this
    // heuristic is part of the output definition: it must not depend on the machine or on --threads.
    const long long shrink_threads = 8;
    while (!opt.block && nb > 4 * marg) {
        long long nbk = (span + (long long)(nb - 2 * marg) * Q - 1) / ((long long)(nb - 2 * marg) * Q);
        if ((long long)ch * ((nbk + 3) / 4) >= 2 * shrink_threads) break;
        int nn = fix(std::max(4 * marg, nb * 3 / 4)); if (nn >= nb) break; nb = nn;
    }
    int Nb = nb * Q, Mb = nb * P, hop = nb - 2 * marg;
    long long nblk = (span + (long long)hop * Q - 1) / ((long long)hop * Q);
    long long base = (long long)marg * Q, tot = base + nblk * hop * Q + Nb;
    std::vector<std::vector<double>> xx(ch);
    const int train = 4096;
    for (int c = 0; c < ch; c++) {
        xx[c].resize(tot);                       // value-initialised; pads must be zero anyway
        double* X = xx[c].data() + base + pre;
        for (long long i = 0; i < nin; i++) X[i] = in[i * ch + c];
        if (opt.lpc_edges) {
            long long ext = E + 32; double fd[32]; for (int i = 0; i < 32; i++) { dd c, sn; dd_cis(i, 64, c, sn); fd[i] = ((dd(1.0) + c) * 0.5).hi; }   // raised cosine
            std::vector<double> tl = burg_extrap(X, nin, ext); for (int i = 0; i < 32; i++) tl[ext - 32 + i] *= fd[i];
            long long m = std::min<long long>(nin, train); std::vector<double> xr(m); for (long long i = 0; i < m; i++) xr[i] = X[m - 1 - i];
            std::vector<double> hd = burg_extrap(xr.data(), m, ext);
            for (int i = 0; i < 32; i++) hd[ext - 32 + i] *= fd[i];   // = head[:32] *= fade[::-1] after reversal
            for (long long i = 0; i < ext; i++) { X[nin + i] = tl[i]; X[-1 - i] = hd[i]; }
        }
    }
    // If the transition (fc +- main-lobe half-width) comes close to the kernel's own Nyquist (near-unity ratios),
    // the fi/fo-sampled kernel spectrum is aliased onto itself. Then sample the kernel on a finer grid Nk (rate Nk*df)
    // so G is the true continuous response, and fold every bin modulo Mb in the engine (exact for any ratio).
    double fml = std::sqrt(opt.beta * opt.beta + 9.8696044010893586) / (2 * 3.14159265358979323846 * opt.length);
    bool gen = (double)(std::max(fi, (long long)fo) - mn) / 2 < 3 * fml;
    int Jg = 0; std::vector<dd> H;
    if (gen) {
        long long Nk = std::max(Nb, Mb) * 5LL / 4; Nk += Nk & 1; while (!smooth7(Nk / 2)) Nk += 2;
        // kernel rate fk2 = Nk*df, df = fi/Nb; relative cutoff (mn/2)/(fk2/2) = mn*Nb/(Nk*fi)
        __int128 A = (__int128)mn * Nb, B = (__int128)Nk * fi;
        long long g2 = (long long)std::gcd((unsigned long long)mn * (unsigned long long)Nb, (unsigned long long)Nk * (unsigned long long)fi);
        double fk2 = (double)Nk * fi / Nb;
        std::vector<dd> h2 = kernel_taps((int)std::floor(opt.length * fk2), (long long)(A / g2), (long long)(B / g2), opt.beta);
        H = kernel_spectrum(h2, (int)Nk, (double)Nb);
        Jg = (int)std::min<long long>(Nk / 2, (long long)std::ceil((mn / 2.0 + 4 * fml) * Nb / fi));
    } else {
        std::vector<dd> h = kernel(fk, mn / gg, fk / gg, opt.length, opt.beta);
        H = kernel_spectrum(h, up ? Mb : Nb, (double)Nb);
    }
    const double t1 = now();
    threads = (int)std::max<long long>(1, std::min<long long>(threads, (long long)ch * ((nblk + 3) / 4)));
    long long o = pre / Q * P;
    switch (opt.precision) {
        case Precision::Fast: run_engine<double, double, V4>(xx, out, P, Q, up, nb, marg, hop, (int)nblk, o, nout, H, threads, Jg); break;
        case Precision::LongDouble: run_engine<LD, LD, LD>(xx, out, P, Q, up, nb, marg, hop, (int)nblk, o, nout, H, threads, Jg); break;
        default: run_engine<dd, dd, DD4>(xx, out, P, Q, up, nb, marg, hop, (int)nblk, o, nout, H, threads, Jg); break;
    }
    if (stats) {
        stats->isa = simd_name;
        stats->p = P; stats->q = Q; stats->taps = 2 * Lk + 1; stats->block = nb; stats->threads = threads;
        stats->general_fold = Jg != 0; stats->design_seconds = t1 - t0; stats->process_seconds = now() - t1;
    }
    return out;
}

}  // namespace kbsrc::KBSRC_ISA
