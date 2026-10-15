#!/usr/bin/env bash
# Fetches and builds the resamplers used by bench/run.sh into bench/tools (git-ignored).
# Needs: git, curl, cmake, ninja, gcc/g++ (w64devkit on Windows).
set -euo pipefail
cd "$(dirname "$0")"
T=$PWD/tools S=$PWD/tools/src P=$PWD/tools/prefix
mkdir -p "$S"
FLAGS="-O3 -march=x86-64-v3"   # same ISA as kbsrc's portable preset
CM=(-G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -Wno-dev -DCMAKE_INSTALL_PREFIX="$P")

# fetch_git dir url commit; fetch_tar dir url sha256 (the tarball unpacks to dir)
fetch_git() { [ -d "$S/$1" ] || git clone -q "$2" "$S/$1"; git -C "$S/$1" checkout -q "$3"; }
fetch_tar() {
  [ -d "$S/$1" ] && return
  curl -sSLf -o "$S/$1.tar" "$2"
  echo "$3  $S/$1.tar" | sha256sum -c -
  tar xf "$S/$1.tar" -C "$S"; rm "$S/$1.tar"
}

# ardftsrc-rs v0.3.1, official Windows release binary
if [ ! -x "$T/ardftsrc-rs/ardftsrc-rs.exe" ] && [ "${OS:-}" = Windows_NT ]; then
  mkdir -p "$T/ardftsrc-rs"
  curl -sSLf -o "$T/ardftsrc-rs/a.zip" \
    https://github.com/phayes/ardftsrc-rs/releases/download/v0.3.1/ardftsrc-rs-x86_64-pc-windows-msvc.zip
  echo "f2c51a801225879ed2a0f5881e119182368da83ab70d94d5edded6ebd434c6a5  $T/ardftsrc-rs/a.zip" | sha256sum -c -
  (cd "$T/ardftsrc-rs" && unzip -o -q a.zip && rm a.zip)
fi

# soxr 0.1.3
fetch_git soxr https://git.code.sf.net/p/soxr/code 945b592b70470e29f917f4de89b4281fbbd540c0
cmake -S "$S/soxr" -B "$S/soxr/build" "${CM[@]}" -DCMAKE_C_FLAGS="$FLAGS" -DBUILD_SHARED_LIBS=OFF \
  -DWITH_OPENMP=OFF -DBUILD_TESTS=OFF -DBUILD_EXAMPLES=OFF -DWITH_LSR_BINDINGS=OFF >/dev/null
cmake --build "$S/soxr/build" >/dev/null && cmake --install "$S/soxr/build" >/dev/null

# libsamplerate 0.2.2
fetch_tar libsamplerate-0.2.2 \
  https://github.com/libsndfile/libsamplerate/releases/download/0.2.2/libsamplerate-0.2.2.tar.xz \
  3258da280511d24b49d6b08615bbe824d0cacc9842b0e4caf11c52cf2b043893
cmake -S "$S/libsamplerate-0.2.2" -B "$S/libsamplerate-0.2.2/build" "${CM[@]}" -DCMAKE_C_FLAGS="$FLAGS" \
  -DBUILD_SHARED_LIBS=OFF -DBUILD_TESTING=OFF -DLIBSAMPLERATE_EXAMPLES=OFF -DLIBSAMPLERATE_INSTALL=ON >/dev/null
cmake --build "$S/libsamplerate-0.2.2/build" >/dev/null && cmake --install "$S/libsamplerate-0.2.2/build" >/dev/null

# speexdsp 1.2.1 (resampler only, built into srcbench)
fetch_tar speexdsp-1.2.1 https://downloads.xiph.org/releases/speex/speexdsp-1.2.1.tar.gz \
  8c777343e4a6399569c72abc38a95b24db56882c83dbdb6c6424a5f4aeb54d3d
mkdir -p "$T/gen"
printf '%s\n' '#pragma once' '#include <stdint.h>' 'typedef int16_t spx_int16_t;' 'typedef uint16_t spx_uint16_t;' \
  'typedef int32_t spx_int32_t;' 'typedef uint32_t spx_uint32_t;' > "$T/gen/speexdsp_config_types.h"

# SoX 14.4.2 (last release). GCC 14+ rejects some of its old C, so those errors are downgraded to warnings.
fetch_git sox https://git.code.sf.net/p/sox/code 45b161d73ec087a8e003747b1aed07cd33589bca
SOX_WNO="-Wno-error=implicit-function-declaration -Wno-error=incompatible-pointer-types -Wno-error=int-conversion"
cmake -S "$S/sox" -B "$S/sox/build" "${CM[@]}" -DCMAKE_C_FLAGS="$FLAGS -std=gnu17 $SOX_WNO" >/dev/null
cmake --build "$S/sox/build" --target sox >/dev/null

# srcbench harness
SPX=$S/speexdsp-1.2.1
gcc $FLAGS -c -DFLOATING_POINT -DUSE_SSE -DUSE_SSE2 -DEXPORT= -I"$SPX/include" -I"$SPX/include/speex" -I"$T/gen" \
  "$SPX/libspeexdsp/resample.c" -o "$T/speex_resample.o"
g++ $FLAGS -std=c++20 -ffp-contract=off -I"$P/include" -I"$SPX/include/speex" -I"$T/gen" -I../src \
  srcbench.cpp ../src/io/wav.cpp "$T/speex_resample.o" -L"$P/lib" -lsoxr -lsamplerate -o "$T/srcbench" -static
g++ $FLAGS -std=c++20 -ffp-contract=off -I../src gen_noise.cpp ../src/io/wav.cpp -o "$T/gen_noise" -static
echo "tools ready in $T"
