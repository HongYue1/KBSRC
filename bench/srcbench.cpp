// Benchmark harness for library resamplers (soxr, libsamplerate, speexdsp).
// Uses the same WAV reader/writer as kbsrc so file I/O cost is identical.
// usage: srcbench soxr|soxr-vhq|lsr|speex in.wav out.wav rate
#include <cmath>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <exception>
#include <stdexcept>
#include <string>
#include <vector>

#include <samplerate.h>
#include <soxr.h>
#include <speex_resampler.h>

#include "io/wav.hpp"

namespace {

using Interleaved = std::vector<double>;

// soxr 0.1.3: given preset, linear phase, 99.5% bandwidth, float64 I/O, 1 thread
Interleaved run_soxr(const Interleaved& x, long long n, int ch, int fi, int fo, long long nout, unsigned long recipe) {
    Interleaved y((size_t)nout * ch);
    soxr_io_spec_t io = soxr_io_spec(SOXR_FLOAT64_I, SOXR_FLOAT64_I);
    soxr_quality_spec_t q = soxr_quality_spec(recipe, SOXR_LINEAR_PHASE);
    q.passband_end = 0.995;
    soxr_runtime_spec_t rt = soxr_runtime_spec(1);
    soxr_error_t e = nullptr;
    soxr_t r = soxr_create(fi, fo, ch, &e, &io, &q, &rt);
    if (!r) throw std::runtime_error(std::string("soxr: ") + e);
    // stream in 64k-frame blocks (one-shot on the whole file is much slower)
    const size_t B = 65536;
    size_t ipos = 0, opos = 0;
    while (opos < (size_t)nout) {
        size_t ilen = ipos < (size_t)n ? std::min(B, (size_t)n - ipos) : 0;
        size_t idone = 0, odone = 0;
        e = soxr_process(r, ilen ? x.data() + ipos * ch : nullptr, ilen, &idone,
                         y.data() + opos * ch, (size_t)nout - opos, &odone);
        if (e) { soxr_delete(r); throw std::runtime_error(std::string("soxr: ") + e); }
        ipos += idone; opos += odone;
        if (!ilen && !odone) break;
    }
    soxr_delete(r);
    return y;
}

// libsamplerate 0.2.2: SRC_SINC_BEST_QUALITY (float32 API)
Interleaved run_lsr(const Interleaved& x, long long n, int ch, int fi, int fo, long long nout) {
    std::vector<float> in(x.begin(), x.end()), out(((size_t)nout + 64) * ch);
    SRC_DATA d{};
    d.data_in = in.data(); d.input_frames = n;
    d.data_out = out.data(); d.output_frames = nout + 64;
    d.src_ratio = (double)fo / fi; d.end_of_input = 1;
    if (int e = src_simple(&d, SRC_SINC_BEST_QUALITY, ch)) throw std::runtime_error(std::string("libsamplerate: ") + src_strerror(e));
    Interleaved y((size_t)nout * ch, 0.0);
    for (long long i = 0; i < std::min<long long>(nout, d.output_frames_gen) * ch; i++) y[i] = out[i];
    return y;
}

// speexdsp 1.2.1: quality 10 (float32 API), leading zeros skipped and tail flushed so output is aligned
Interleaved run_speex(const Interleaved& x, long long n, int ch, int fi, int fo, long long nout) {
    int err = 0;
    SpeexResamplerState* st = speex_resampler_init(ch, fi, fo, 10, &err);
    if (!st || err) throw std::runtime_error("speex: init failed");
    speex_resampler_skip_zeros(st);
    const long long lat = speex_resampler_get_input_latency(st);
    std::vector<float> in((size_t)(n + lat + 1) * ch, 0.0f), out((size_t)nout * ch);
    for (size_t i = 0; i < x.size(); i++) in[i] = (float)x[i];
    long long ip = 0, op = 0, total = n + lat + 1;
    while (op < nout && ip < total) {
        spx_uint32_t il = (spx_uint32_t)std::min<long long>(total - ip, 1 << 20);
        spx_uint32_t ol = (spx_uint32_t)std::min<long long>(nout - op, 1 << 20);
        speex_resampler_process_interleaved_float(st, in.data() + ip * ch, &il, out.data() + op * ch, &ol);
        ip += il; op += ol;
        if (il == 0 && ol == 0) break;
    }
    speex_resampler_destroy(st);
    return Interleaved(out.begin(), out.end());
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 5) { std::fprintf(stderr, "usage: srcbench soxr|soxr-vhq|lsr|speex in.wav out.wav rate\n"); return 2; }
        const std::string lib = argv[1];
        const int fo = std::stoi(argv[4]);
        kbsrc::io::Audio a = kbsrc::io::read_wav(argv[2]);
        const long long nout = (a.frames * fo + a.rate - 1) / a.rate;
        Interleaved y;
        if (lib == "soxr") y = run_soxr(a.d, a.frames, a.ch, a.rate, fo, nout, SOXR_32_BITQ);
        else if (lib == "soxr-vhq") y = run_soxr(a.d, a.frames, a.ch, a.rate, fo, nout, SOXR_VHQ);
        else if (lib == "lsr") y = run_lsr(a.d, a.frames, a.ch, a.rate, fo, nout);
        else if (lib == "speex") y = run_speex(a.d, a.frames, a.ch, a.rate, fo, nout);
        else { std::fprintf(stderr, "unknown library %s\n", lib.c_str()); return 2; }
        std::vector<std::vector<double>> planar(a.ch, std::vector<double>((size_t)nout));
        for (long long i = 0; i < nout; i++)
            for (int c = 0; c < a.ch; c++) planar[c][i] = y[(size_t)i * a.ch + c];
        long long clipped = 0;
        kbsrc::io::write_wav(argv[3], planar, a.ch, nout, fo, "f64", false, clipped);
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "srcbench: %s\n", e.what());
        return 1;
    }
}
