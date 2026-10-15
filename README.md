# KBSRC

An FFT sample-rate converter for WAV files. It filters with a Kaiser-windowed sinc whose transition band ends at the lower of the two Nyquist frequencies, so nothing above it aliases or images, and it runs in double-double (~106-bit) arithmetic by default.

The approach follows [ARDFTSRC](https://github.com/mrspoonsi/ARDFTSRC): rational resampling in the DFT domain. KBSRC is an independent C++ implementation focused on speed: double-double arithmetic instead of 128-bit floats costs a fraction of the time, and LPC edge extension avoids pre-echo at the start and end of a file. Its filter does not straddle the Nyquist frequency: the whole transition band lies below it.

## Building

You need CMake 3.20+, Ninja, and GCC 12+ or Clang 15+. MSVC isn't supported. It builds on x86-64 and ARM64.

```
cmake --preset release
cmake --build --preset release
```

The binary is written to `build/release/kbsrc`. The `release` preset tunes for the build machine. Use `portable` for a binary that runs on any CPU of the same architecture: on x86-64 it contains an AVX2 engine and a generic one and picks at run time. Set `KBSRC_CPU=generic` to force the generic engine. Output is bit-identical across engines, thread counts, compilers and platforms (checked with MinGW GCC on Windows, GCC and Clang on Linux, and ARM64), because every step uses only IEEE basic operations. The exception is `ld`, which uses the x87 80-bit format and so matches only on x86-64. On Windows, [w64devkit](https://github.com/skeeto/w64devkit) works out of the box and produces a static `kbsrc.exe`.

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
| `-t, --threads N` | up to 8 (4 for `fast`) | worker threads |
| `--len S`, `--beta B` | 0.2, 40 | kernel half-length in seconds, Kaiser beta |
| `--no-lpc` | | zero-pad the edges instead of extrapolating them |
| `--block N` | automatic | block size nb (see [How it works](#how-it-works)); it changes the rounding of the output |
| `-v` | | print the plan and timings |

Input can be 8/16/24/32-bit PCM or 32/64-bit float, with any number of channels. Any pair of integer sample rates from 1 Hz to 10 MHz works if the transition band fits below the lower Nyquist frequency with room for a passband: the band is 2·√(β² + π²) / (2π·S) Hz wide, 64 Hz at the defaults, so the lower rate must be at least 128 Hz; a longer `--len` allows lower rates. Converting 96 kHz to 44.1 kHz with the defaults, the response is flat to 21986 Hz and in the stopband from 22050 Hz.

## Comparison

All runs were done on one machine (Ryzen 7 5800H, Windows 11). Time is the median of 5 runs after warm-up, converting a 180-second stereo 96 kHz noise file to 44.1 kHz, 64-bit float in and out, wall-clock including file I/O. The KBSRC rows were measured on 2026-10-15 at a fixed 3600 MHz using `bench/run.sh` and the portable GCC 16.1 build of 1.3.0. Competitor timings are from a batch on 2026-10-13 on the same machine and clock.

Quality is the Hydrogenaudio balanced score, computed locally with the official scripts in GNU Octave 11.3. KBSRC keeps its transition band below the output Nyquist frequency (see [How it works](#how-it-works)), so nothing above that frequency folds back into the output. Letting the band straddle Nyquist scores higher on bandwidth but folds content back: release 1.2.4 did that and scored 99.93 (exact) and 99.87 (fast).

| Resampler | Version | Settings | Threads | Score | Time |
|---|---|---|---|---|---|
| kbsrc | 1.3.0 | `--mode exact` | 8 | 99.85 | 1.331 s |
| kbsrc | 1.3.0 | `--mode exact -t 1` | 1 | 99.85 | 3.033 s |
| kbsrc | 1.3.0 | `--mode fast` | 4 | 99.79 | 0.702 s |
| kbsrc | 1.3.0 | `--mode fast -t 1` | 1 | 99.79 | 0.965 s |
| ardftsrc-rs | 0.3.1 | Insane, f64 | default | 99.74 | 2.250 s |
| ardftsrc-rs | 0.3.1 | Insane, double-double | default | 99.92 | 30.985 s |
| soxr | 0.1.3 | `SOXR_32_BITQ`, linear phase, 99.5% | 1 | 82.56 | 0.697 s |
| soxr | 0.1.3 | `SOXR_VHQ`, linear phase, 99.5% | 1 | 81.88 | 0.670 s |
| SoX | 14.4.2 | `rate -v -L -b 99.5` | 1 | 72.71 | 1.289 s |
| libsamplerate | 0.2.2 | `SRC_SINC_BEST_QUALITY` | 1 | 68.23 | 14.778 s |
| speexdsp | 1.2.1 | quality 10, float | 1 | 62.11 | 9.330 s |

"Insane" is `--quality 262144 --bandwidth 0.995 --taper-type kaiser --beta 32 --alias-floor-db 0`, the settings behind dBpoweramp's ARDFTSRC Insane entry. ardftsrc-rs is the official Windows release. Everything else was built from source with GCC 16.1 at `-O3 -march=x86-64-v3`: soxr from tag `0.1.3` (945b592), SoX from tag `sox-14.4.2` (45b161d), and libsamplerate and speexdsp from their release tarballs. The libraries are called through `bench/srcbench.cpp`. Most of the gap in the lower rows comes from the aliasing, impulse, intermodulation and gapless tests.

The [Hydrogenaudio SRC comparison](https://src.hydrogenaudio.org/) ranks resamplers with the same test scripts.

To reproduce, put the Hydrogenaudio test pack in `tests/data/`, then run `bench/build_tools.sh`, `bench/run.sh` and `bench/score.sh`. Results go to `bench/out/`.

### AVX2 vs generic

KBSRC 1.3.0 on the same Ryzen 7 5800H at 3600 MHz, using the portable build and the same 180-second f64 noise conversion. Each engine was warmed up, then timed in five alternating trials with Python subprocess wall-clock timing, including setup and file I/O.

| Mode | Threads | Generic | AVX2 | AVX2 speedup |
|---|---:|---:|---:|---:|
| exact | 1 | 9.975 s | 2.957 s | 3.37× |
| exact | 8 | 2.192 s | 1.298 s | 1.69× |
| fast | 1 | 1.724 s | 0.909 s | 1.90× |
| fast | 4 | 0.870 s | 0.602 s | 1.44× |

AVX2 was faster in every paired trial. Outputs matched byte-for-byte in all four comparisons. Set `KBSRC_CPU=generic` to force the portable array backend; otherwise the binary selects AVX2 when supported. Raw trials: `bench/out/_speed/engines-1.3.0-3600mhz-20261015.json`.

## How it works

- Each block is processed with overlap-save in the rational domain. Nb = nb·Q input bins map to Mb = nb·P output bins, where P/Q is the reduced rate ratio. By default a block spans ten kernel margins in `exact` and `ld`. `fast` is limited by memory traffic, so it aims for transforms of about 51200 points, with four to ten margins per block. nb is rounded up to an even number whose half is 7-smooth (5-smooth for automatic `fast` blocks).
- The cutoff sits half a transition band below the lower Nyquist frequency, so the response is in the stopband (below -300 dB with the default beta) from that frequency up. When downsampling, the spectrum bins above the output Nyquist are folded onto the output bins as decimation requires, but they hold only stopband residue. When upsampling, the input spectrum is extended periodically and the images lie in the stopband.
- Rates that are very close together, like 44100 and 44101, put the transition band on the kernel's own Nyquist. In that case the kernel is sampled on a finer grid and every bin is folded modulo Mb.
- The kernel taps (sinc and Kaiser window) are computed in double-double, with exact integer reduction of the sinc argument.
- The signal edges are extended with Burg linear prediction (order 32), which avoids the pre-echo that zero padding causes. The Burg recursion also runs in double-double.
- The FFT is custom. It uses mixed-radix Stockham passes (2, 3, 4, 5, 7) with a four-step split that runs four columns per vector (AVX2, or plain arrays the compiler vectorises). Other lengths go through a Bluestein transform. Twiddle factors are exact in double-double.
- When the input and output rates are equal, the samples are copied unchanged.

## Tests

```
ctest --preset release
```

The tests need Python 3 with NumPy and SciPy. `test_reference` compares the output against the Python implementation in `reference/python`. It needs the 64-bit float WAVs from the Hydrogenaudio test pack in `tests/data/`, and CTest reports it as skipped if any is missing.

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
