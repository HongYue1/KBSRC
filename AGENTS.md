# AGENTS.md

Rules for anyone changing KBSRC, human or automated. README.md describes what the converter does and how it
is used; this file covers how to build, test and change it without breaking its guarantees.

## Layout

- `include/kbsrc/kbsrc.hpp`: public API (`kbsrc::resample`, `Options`, `Stats`) and its floating-point contract.
- `src/dispatch.cpp`: picks the engine at run time and runs Exact in round-to-nearest.
- `src/resampler.cpp`, `src/design.cpp`, `src/lpc.cpp`, `src/dsp/`: the engine (block planning, kernel design,
  Burg edge extension, double-double arithmetic, SIMD types, FFT). It is compiled once per instruction set.
- `src/io/`: WAV reader and writer, used by the command-line tool and the bench tools, not by the library.
- `app/main.cpp`: command-line tool. `tests/`: CTest scripts and `fenv_driver.cpp`, a test-only executable.
- `reference/python/`: NumPy reference implementation used by `tests/test_reference.py`.
- `bench/`: quality and speed comparison against other resamplers.
- `build/`, `out/`, `tests/data/`, `docs/`, `bench/tools/` and `bench/out/` are git-ignored.

## Build

CMake 3.20+, Ninja, and GCC 12+ or Clang 15+. MSVC is rejected (the code needs `__int128` and x87
`long double`). On Windows use w64devkit (MinGW-w64 GCC); its executables are linked statically.

```
cmake --preset portable && cmake --build --preset portable   # build/portable: generic + AVX2 engines
cmake --preset release && cmake --build --preset release     # build/release: tuned with -march=native
cmake --preset debug && cmake --build --preset debug         # build/debug: Debug, KBSRC_ARCH=native
```

- `portable` is the configuration that ships and the one to test and time. On x86-64 it compiles the generic
  engine for baseline x86-64 and a second copy with `-march=x86-64-v3`, chosen at run time;
  `KBSRC_CPU=generic` forces the generic engine. `native` builds have a single engine.
- For a Debug build with both engines, configure a separate directory:
  `cmake -S . -B build/debug-portable -G Ninja -DCMAKE_BUILD_TYPE=Debug -DKBSRC_ARCH=portable`.
- Every change must build warning-free (`-Wall -Wextra`) in Release and Debug.

## Test

```
ctest --test-dir build/portable --output-on-failure      # or: ctest --preset release
```

| Test | Needs | Checks |
|---|---|---|
| `robust` | Python 3 standard library | rejected options and WAV headers fail cleanly; valid neighbours convert |
| `engines` | Python 3 standard library, `kbsrc_fenv_driver` | pinned outputs; engine, rounding-mode, thread identity |
| `formats` | NumPy | WAV formats, edge cases, command-line errors |
| `rates` | NumPy | sine accuracy over common, odd and near-unity ratios |
| `alias` | NumPy | nothing at or above the lower Nyquist frequency passes the filter (no aliasing or imaging) |
| `reference` | NumPy, SciPy, `tests/data/` | agreement with `reference/python` on the Hydrogenaudio test WAVs |
| `isa` | Python 3 standard library, objdump | portable x86-64 only: no AVX instruction outside the AVX2 engine |

- `test_robust.py`, `test_engines.py` and `test_isa.py` must stay standard-library only, so they run on any
  machine. They do not import `tests/common.py`, which needs NumPy.
- When the Python that CMake found has no NumPy/SciPy, run those suites directly, for example
  `uv run --no-project --with numpy,scipy python tests/test_rates.py --bin build/portable/kbsrc.exe`.
- `reference` needs the 64-bit float WAVs of the Hydrogenaudio test pack in `tests/data/` (`KBSRC_TEST_DATA`
  overrides the path) and reports a skip without them; read its output, a passing CTest row may be a skip.
  With GCC 16.1 on Windows four checks currently exceed their limits: gaplesstest_s exact (5.4e-15),
  sweep_1_to_44KHz fast and exact (5.7e-6) and aliasing150db exact (8.4e-15). Treat any other failure as a
  regression, and a change in these values as a change of output that needs an explanation.
- Run every suite on Release and on a Debug build with both engines before committing a code change.

## Filter

- The kernel is a Kaiser-windowed sinc whose transition band is 2 * fml wide, fml = sqrt(beta^2 + pi^2) /
  (2 pi len). It must end at or below the lower Nyquist frequency, so nothing above that frequency is folded
  back into the output (downsampling) or passed as an image (upsampling) above the stopband level. Do not move
  the cutoff up to gain bandwidth; `alias` measures the response the way the Hydrogenaudio impulse test does.
- The defaults (`--len 0.2`, `--beta 40`) give the best Hydrogenaudio balanced score under that rule: shorter
  kernels lose bandwidth, longer ones ring longer after a sweep crosses the cutoff (spectrogram test), and
  beyond about 0.29 s a kernel overlaps the neighbouring impulses of the impulse test, which are 0.6 s apart.
  A lower beta narrows the transition band but raises the stopband near the Nyquist frequency. Change the
  defaults only with scores measured by `bench/score.sh`.

## Numerical rules

The output is defined to the bit. These rules protect that:

- Every target compiles with `-ffp-contract=off` through the `kbsrc_options` interface library; link new
  targets against it. Never add `-ffast-math`, `-Ofast` or any flag that reassociates, contracts, assumes
  finite values or flushes denormals.
- The double-double code depends on IEEE rounding of each basic operation in round-to-nearest. Do not
  reorder or fuse its operations, and do not replace error-free transforms with "equivalent" algebra.
- Hardware FMA is used only in engine copies compiled for FMA hardware, through `__builtin_fma`. Do not call
  `std::fma`/`fma` from libm in hot code: unoptimised builds call the C runtime, and MinGW-w64's `fma` is not
  correctly rounded.
- The generic engine must stay baseline x86-64 in portable builds: no `-march`, target attributes or
  intrinsics beyond SSE2 in code compiled into it. Instruction-set-specific code belongs in the per-ISA
  engine copies (`src/isa.hpp`, `kbsrc_engine()` in CMakeLists.txt) or behind a run-time CPU check.
- A portable build must start on any x86-64 CPU. Objects at namespace scope in engine code must be
  constant-initialised (`constexpr`): a dynamic initialiser in the AVX2 copy is AVX code that runs before
  `main` on every CPU. The `isa` test disassembles the executable and fails for any AVX instruction outside the
  AVX2 engine, including start-up code and target-attribute functions such as the SSE4.1 WAV writer.
- Output must be bit-identical across engines, thread counts, compilers and platforms, and for Exact also
  across the caller's rounding modes (`resample()` switches to round-to-nearest and restores the caller's
  mode). Block size affects rounding, so block planning may depend only on the input and the options, never
  on the machine or the thread count. Flush-to-zero and denormals-are-zero must be off.
- `tests/test_engines.py` pins SHA-256 hashes and samples of known outputs. A change that is not meant to
  alter the output must leave every pin unchanged; one that is meant to must update the pins in the same
  commit and explain why.
- Speed claims need measurements: portable build, both engines, median of at least 5 runs alternating old and
  new binaries on a fixed CPU clock. Report slowdowns as well as speed-ups, and leave the README tables alone
  unless they are re-measured the way they describe.

## Dependencies and pins

- The library and the command-line tool use only the C++20 standard library and threads (plus `shell32` on
  Windows). Do not add libraries, submodules, `FetchContent` or new Python packages; the test suites use only
  NumPy and SciPy beyond the standard library.
- `bench/build_tools.sh` fetches the comparison resamplers at pinned versions (git tags, release tarballs, and
  the ardftsrc-rs release zip checked by SHA-256) into `bench/tools/`. Change a pin only deliberately, with the
  new tag or checksum in the script and re-measured results.

## Code and commits

- C++20. 4-space indentation (2 in CMake, JSON and YAML), UTF-8, LF line endings, final newline, no trailing
  whitespace, lines up to 120 columns (`.editorconfig`, `.gitattributes`).
- Bad arguments to the library throw `std::invalid_argument`. The command-line tool prints `kbsrc: <reason>`,
  exits with status 1 (2 for usage errors) and writes no output file for a rejected conversion.
- Keep comments and README statements true after a change; update them in the same commit.
- One logical change per commit. Subject in the imperative, at most 72 characters, no trailing period; then a
  blank line and a body wrapped at 72 columns that explains the problem, the change and how it was verified.
- Describe behaviour. Commit messages, comments, test names and test labels must not refer to review
  findings, tickets or task lists by code or number; say what was wrong and what changed instead.
- Never commit build outputs, archives, WAV files, test data, local notes or editor settings. `docs/` holds the
  maintainer's local notes and journal of experiments that were tried and rejected; it is not part of the
  repository. Add rejected experiments there instead of retrying them.
- Do not push, tag or publish releases unless the maintainer asks. The version is `project(VERSION)` in
  CMakeLists.txt and changes only with a release.
