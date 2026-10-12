# KBSRC

An FFT sample-rate converter for WAV files. It filters with a Kaiser-windowed sinc whose transition band is centred on the output Nyquist frequency, and it runs in double-double (~106-bit) arithmetic by default.

The approach follows [ARDFTSRC](https://github.com/mrspoonsi/ARDFTSRC): DFT-domain rational resampling with a Kaiser taper that straddles the output Nyquist. KBSRC is an independent C++ implementation focused on speed: double-double arithmetic instead of 128-bit floats gives the same top score in a fraction of the time, and LPC edge extension avoids pre-echo at the start and end of a file.

## Building

You need CMake 3.20+, Ninja, and GCC 12+ or Clang 15+. MSVC isn't supported. It builds on x86-64 and ARM64.

```
cmake --preset release
cmake --build --preset release
```

The binary is written to `build/release/kbsrc`. The `release` preset tunes for the build machine. Use `portable` for a binary that runs on any CPU of the same architecture: on x86-64 it contains an AVX2 engine and a generic one and picks at run time. The generic engine is about 1.4x slower in exact mode and 1.15x slower in fast mode (1.85 s vs 1.34 s and 0.75 s vs 0.65 s on the benchmark below). Set `KBSRC_CPU=generic` to force it. Both engines give bit-identical output. Other compilers or platforms (such as ARM64) can differ by an ULP or two, below -300 dBFS, because `long double` and the maths library differ. On Windows, [w64devkit](https://github.com/skeeto/w64devkit) works out of the box and produces a static `kbsrc.exe`.

## Usage

```
kbsrc input.wav output.wav --rate 48000
```

| Option | Default | |
|---|---|---|
| `-r, --rate N` | 44100 | output sample rate |
| `-m, --mode M` | `exact` | `exact` (double-double), `fast` (float64) or `ld` (x87 80-bit) |
| `-f, --format F` | `f64` | `f64`, `f32`, `s32`, `s24` or `s16` |
| `--dither` | off | TPDF dither for integer output |
| `-t, --threads N` | 8 (exact), 4 (fast) | worker threads |
| `--len S`, `--beta B` | 0.1, 40 | kernel half-length in seconds, Kaiser beta |
| `--no-lpc` | | zero-pad the edges instead of extrapolating them |
| `-v` | | print the plan and timings |

Input can be 8/16/24/32-bit PCM or 32/64-bit float, with any number of channels. Any pair of integer sample rates works.

## Comparison

All runs were done on one machine (Ryzen 7 5800H, Windows 11). Quality is the Hydrogenaudio balanced score, computed locally with the official test scripts in GNU Octave 11.3. Time is the median of 5 runs converting a 180-second stereo 96 kHz noise file to 44.1 kHz, 64-bit float in and out, wall-clock including file I/O.

| Resampler | Version | Settings | Threads | Score | Time |
|---|---|---|---|---|---|
| kbsrc | 1.2.0 | `--mode exact` | 8 | 99.93 | 1.40 s |
| kbsrc | 1.2.0 | `--mode exact -t 1` | 1 | 99.93 | 4.21 s |
| kbsrc | 1.2.0 | `--mode fast` | 4 | 99.87 | 0.67 s |
| kbsrc | 1.2.0 | `--mode fast -t 1` | 1 | 99.87 | 0.83 s |
| ardftsrc-rs | 0.3.1 | Insane, f64 | default | 99.74 | 2.46 s |
| ardftsrc-rs | 0.3.1 | Insane, double-double | default | 99.92 | 31.1 s |
| soxr | 0.1.3 | `SOXR_32_BITQ`, linear phase, 99.5% | 1 | 82.56 | 0.78 s |
| soxr | 0.1.3 | `SOXR_VHQ`, linear phase, 99.5% | 1 | 81.88 | 0.74 s |
| SoX | 14.4.2 | `rate -v -L -b 99.5` | 1 | 72.71 | 1.30 s |
| libsamplerate | 0.2.2 | `SRC_SINC_BEST_QUALITY` | 1 | 68.23 | 14.95 s |
| speexdsp | 1.2.1 | quality 10, float | 1 | 62.11 | 9.46 s |

"Insane" is `--quality 262144 --bandwidth 0.995 --taper-type kaiser --beta 32 --alias-floor-db 0`, the settings behind dBpoweramp's ARDFTSRC Insane entry. ardftsrc-rs is the official Windows release. Everything else was built from source with GCC 16.1 at `-O3 -march=x86-64-v3`: soxr from tag `0.1.3` (945b592), SoX from tag `sox-14.4.2` (45b161d), and libsamplerate and speexdsp from their release tarballs. The libraries are called through `bench/srcbench.cpp`. Most of the gap in the lower rows comes from the aliasing, impulse, intermodulation and gapless tests.

KBSRC achieves 99.93% on the Hydrogenaudio SRC comparison ([ranking](https://src.hydrogenaudio.org/), [result](https://src.hydrogenaudio.org/compareresults?id1=b7d5b47e-dd06-41cd-892a-b0b73d6e433d&id2=0)).

To reproduce, put the Hydrogenaudio test pack in `tests/data/`, then run `bench/build_tools.sh`, `bench/run.sh` and `bench/score.sh`. Results go to `bench/out/`.

## How it works

- Each block is processed with overlap-save in the rational domain. Nb = nb·Q input bins map to Mb = nb·P output bins, where P/Q is the reduced rate ratio.
- When downsampling, spectrum bins above the output Nyquist are folded back onto it, so the transition band is power-complementary. When upsampling, the input spectrum is extended periodically.
- Rates that are very close together, like 44100 and 44101, put the transition band on the kernel's own Nyquist. In that case the kernel is sampled on a finer grid and every bin is folded modulo Mb.
- The signal edges are extended with Burg linear prediction (order 32), which avoids the pre-echo that zero padding causes.
- The FFT is custom. It uses mixed-radix Stockham passes (2, 3, 4, 5, 7) with a four-step split that runs four columns per vector (AVX2, or plain arrays the compiler vectorises). Other lengths go through a Bluestein transform. Twiddle factors are exact in double-double.
- When the input and output rates are equal, the samples are copied unchanged.

## Tests

```
ctest --preset release
```

The tests need Python 3 with NumPy and SciPy. `test_reference` compares the output against the Python implementation in `reference/python`. It needs the 64-bit float WAVs from the Hydrogenaudio test pack in `tests/data/`, and it is skipped if they're missing.

## Layout

```
include/kbsrc/   public API (kbsrc::resample)
src/             library: design, LPC, engine
src/dsp/         double-double, SIMD and FFT
src/io/          WAV reader and writer
app/             command-line tool
tests/           CTest scripts
reference/       Python reference implementation
```

## License

MIT, see [LICENSE](LICENSE).
