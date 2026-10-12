// kbsrc command-line tool
#include <chrono>
#include <cmath>
#include <cstdio>
#include <exception>
#include <string>
#include <vector>

#include "io/wav.hpp"
#include "kbsrc/kbsrc.hpp"

#ifdef _WIN32
#include <windows.h>
#include <shellapi.h>
#endif

namespace {

struct UsageError : std::exception {
    std::string msg;
    explicit UsageError(std::string m) : msg(std::move(m)) {}
    const char* what() const noexcept override { return msg.c_str(); }
};

double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

int to_int(const std::string& v, const char* name) {
    try { size_t p; int r = std::stoi(v, &p); if (p == v.size()) return r; } catch (...) {}
    throw UsageError(std::string("invalid value for ") + name + ": " + v);
}

double to_double(const std::string& v, const char* name) {
    try { size_t p; double r = std::stod(v, &p); if (p == v.size() && std::isfinite(r)) return r; } catch (...) {}
    throw UsageError(std::string("invalid value for ") + name + ": " + v);
}

void print_help() {
    std::printf(
        "kbsrc %s - Kaiser-Bessel fold-back sample-rate converter\n"
        "\n"
        "usage: kbsrc [options] input.wav output.wav\n"
        "\n"
        "  -r, --rate N        output sample rate in Hz (default 44100)\n"
        "  -m, --mode M        exact (default), fast or ld\n"
        "  -f, --format F      f64 (default), f32, s32, s24 or s16\n"
        "      --dither        TPDF dither for integer formats\n"
        "  -t, --threads N     worker threads (default: 8 for exact, 4 for fast)\n"
        "      --len S         one-sided kernel length in seconds (default 0.1)\n"
        "      --beta B        Kaiser beta (default 40)\n"
        "      --block N       block size factor (default: automatic)\n"
        "      --no-lpc        zero-pad edges instead of LPC extrapolation\n"
        "  -v, --verbose       print plan and timing\n"
        "  -h, --help          show this help\n"
        "      --version       show version\n",
        kbsrc::version());
}

// Windows: re-read the command line as UTF-16 and convert to UTF-8 so non-ASCII paths work.
std::vector<std::string> utf8_args(int argc, char** argv) {
    std::vector<std::string> args;
#ifdef _WIN32
    int wc = 0;
    if (LPWSTR* wv = CommandLineToArgvW(GetCommandLineW(), &wc)) {
        for (int i = 0; i < wc; i++) {
            int n = WideCharToMultiByte(CP_UTF8, 0, wv[i], -1, nullptr, 0, nullptr, nullptr);
            std::string u(n, '\0');
            WideCharToMultiByte(CP_UTF8, 0, wv[i], -1, u.data(), n, nullptr, nullptr);
            u.resize(n - 1);
            args.push_back(std::move(u));
        }
        LocalFree(wv);
        return args;
    }
#endif
    for (int i = 0; i < argc; i++) args.emplace_back(argv[i]);
    return args;
}

int run(const std::vector<std::string>& args) {
    kbsrc::Options opt;
    int rate = 44100;
    std::string format = "f64", in_path, out_path;
    bool dither = false, verbose = false;

    for (size_t i = 1; i < args.size(); i++) {
        const std::string& a = args[i];
        auto value = [&]() -> const std::string& {
            if (i + 1 >= args.size()) throw UsageError("missing value for " + a);
            return args[++i];
        };
        if (a == "-r" || a == "--rate") rate = to_int(value(), "--rate");
        else if (a == "-m" || a == "--mode") {
            const std::string& m = value();
            if (m == "exact") opt.precision = kbsrc::Precision::Exact;
            else if (m == "fast") opt.precision = kbsrc::Precision::Fast;
            else if (m == "ld") opt.precision = kbsrc::Precision::LongDouble;
            else throw UsageError("mode must be exact, fast or ld");
        }
        else if (a == "-f" || a == "--format") format = value();
        else if (a == "-t" || a == "--threads") opt.threads = to_int(value(), "--threads");
        else if (a == "--block") opt.block = to_int(value(), "--block");
        else if (a == "--len") opt.length = to_double(value(), "--len");
        else if (a == "--beta") opt.beta = to_double(value(), "--beta");
        else if (a == "--dither") dither = true;
        else if (a == "--no-lpc") opt.lpc_edges = false;
        else if (a == "-v" || a == "--verbose") verbose = true;
        else if (a == "-h" || a == "--help") { print_help(); return 0; }
        else if (a == "--version") { std::printf("kbsrc %s\n", kbsrc::version()); return 0; }
        else if (a.size() > 1 && a[0] == '-') throw UsageError("unknown option " + a);
        else if (in_path.empty()) in_path = a;
        else if (out_path.empty()) out_path = a;
        else throw UsageError("too many arguments");
    }
    if (in_path.empty() || out_path.empty()) throw UsageError("expected input and output paths");
    if (rate < 1 || rate > 10000000) throw UsageError("--rate must be between 1 and 10000000");
    if (format != "f64" && format != "f32" && format != "s32" && format != "s24" && format != "s16")
        throw UsageError("format must be f64, f32, s32, s24 or s16");

    const double t0 = now();
    kbsrc::io::Audio a = kbsrc::io::read_wav(in_path.c_str());
    const double t1 = now();
    kbsrc::Stats st;
    auto y = kbsrc::resample(a.d.data(), a.frames, a.ch, a.rate, rate, opt, &st);
    const double t2 = now();
    long long nout = y.empty() ? 0 : (long long)y[0].size(), clipped = 0;
    kbsrc::io::write_wav(out_path.c_str(), y, a.ch, nout, rate, format, dither, clipped);
    const double t3 = now();

    if (clipped) std::fprintf(stderr, "kbsrc: warning: %lld samples clipped\n", clipped);
    if (verbose) {
        std::fprintf(stderr, "%d Hz -> %d Hz, %d ch, %lld -> %lld frames\n", a.rate, rate, a.ch, a.frames, nout);
        std::fprintf(stderr, "P/Q=%lld/%lld taps=%lld block=%d threads=%d engine=%s%s\n", st.p, st.q, st.taps, st.block,
                     st.threads, st.isa, st.general_fold ? " (general fold)" : "");
        std::fprintf(stderr, "read %.3fs  design %.3fs  resample %.3fs  write %.3fs  total %.3fs\n", t1 - t0,
                     st.design_seconds, t2 - t1 - st.design_seconds, t3 - t2, t3 - t0);
    }
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return run(utf8_args(argc, argv));
    } catch (const UsageError& e) {
        std::fprintf(stderr, "kbsrc: %s (see --help)\n", e.what());
        return 2;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "kbsrc: %s\n", e.what());
        return 1;
    }
}
