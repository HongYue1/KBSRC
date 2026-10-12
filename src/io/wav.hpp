// WAV file reading and writing (PCM 8/16/24/32-bit, IEEE float 32/64, WAVE_FORMAT_EXTENSIBLE)
#pragma once

#include <string>
#include <vector>

namespace kbsrc::io {

struct Audio {
    int rate = 0, ch = 0;
    long long frames = 0;
    std::vector<double> d;  // interleaved
    int fmt_bits = 0;
    bool fmt_float = false;
};

// format: "f64", "f32", "s32", "s24" or "s16". Integer formats round (optionally with TPDF dither) and clip.
// Paths are UTF-8 on every platform. Both functions throw std::runtime_error on failure.
Audio read_wav(const char* path);
void write_wav(const char* path, const std::vector<std::vector<double>>& y, int ch, long long frames, int rate,
               const std::string& format, bool dither, long long& clipped);

}  // namespace kbsrc::io
