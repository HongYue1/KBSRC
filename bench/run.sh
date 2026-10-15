#!/usr/bin/env bash
# Converts the Hydrogenaudio 96 kHz test signals to 44.1 kHz with every resampler (for scoring),
# then times a 180 s stereo 96 kHz noise file (bench/gen_noise.cpp).
# Output: bench/out/<id>/*.wav and bench/out/timing.tsv
# usage: bench/run.sh [runs=5] [ids...]   (SKIP_CONVERT=1 to only re-time)
set -euo pipefail
cd "$(dirname "$0")/.."
RUNS=${1:-5}; shift || true
DATA=tests/data T=bench/tools O=bench/out
KB=build/portable/kbsrc.exe; [ -x "$KB" ] || KB=build/portable/kbsrc
ARD=$T/ardftsrc-rs/ardftsrc-rs.exe
SOX=$T/src/sox/build/src/sox.exe
SB=$T/srcbench.exe
ARD_INSANE="--quality 262144 --bandwidth 0.995 --taper-type kaiser --beta 32 --alias-floor-db 0 --out-format f64"

# id | how to run it ($1 = input, $2 = output)
convert() {
  local id=$1 in=$2 out=$3
  case $id in
    kbsrc-exact)    $KB "$in" "$out" --rate 44100 ;;
    kbsrc-exact-1t) $KB "$in" "$out" --rate 44100 -t 1 ;;
    kbsrc-fast)     $KB "$in" "$out" --rate 44100 --mode fast ;;
    kbsrc-fast-1t)  $KB "$in" "$out" --rate 44100 --mode fast -t 1 ;;
    ardftsrc-f64)   $ARD --input "$in" --output "$out" --output-rate 44100 $ARD_INSANE ;;
    ardftsrc-dd)    $ARD --input "$in" --output "$out" --output-rate 44100 $ARD_INSANE --high-precision double-double ;;
    sox)            $SOX -D "$in" -e floating-point -b 64 "$out" rate -v -L -b 99.5 44100 ;;
    soxr)           $SB soxr "$in" "$out" 44100 ;;
    soxr-vhq)       $SB soxr-vhq "$in" "$out" 44100 ;;
    libsamplerate)  $SB lsr "$in" "$out" 44100 ;;
    speexdsp)       $SB speex "$in" "$out" 44100 ;;
    *) echo "unknown id $id" >&2; return 2 ;;
  esac
}

IDS=("$@")
[ ${#IDS[@]} -gt 0 ] || IDS=(kbsrc-exact kbsrc-exact-1t kbsrc-fast kbsrc-fast-1t ardftsrc-f64 ardftsrc-dd
                             sox soxr-vhq soxr libsamplerate speexdsp)
FILES=(impulse sweep_1_to_44KHz aliasing150db intermodulation_sine gaplesstest_s gaplesstest_m)
mkdir -p "$O"

for id in "${IDS[@]}"; do
  [ "${SKIP_CONVERT:-0}" = 1 ] && continue
  mkdir -p "$O/$id"
  for f in "${FILES[@]}"; do
    convert "$id" "$DATA/$f-64bitfloat.wav" "$O/$id/$f-64bitfloat.wav" 2>>"$O/$id/stderr.log"
  done
  echo "converted $id"
done

ms() { echo $(( $(date +%s%N) / 1000000 )); }
TSV=$O/timing.tsv
[ -f "$TSV" ] || printf "id\truns_ms\tmedian_ms\n" > "$TSV"
IN=$O/_input/noise-96k-stereo.wav TMP=$O/_timing.wav
[ -f "$IN" ] || { mkdir -p "$O/_input"; $T/gen_noise "$IN"; }
for id in "${IDS[@]}"; do
  convert "$id" "$IN" "$TMP" 2>/dev/null   # warm-up (file cache)
  t=()
  for ((r = 0; r < RUNS; r++)); do s=$(ms); convert "$id" "$IN" "$TMP" 2>/dev/null; t+=($(( $(ms) - s ))); done
  med=$(printf "%s\n" "${t[@]}" | sort -n | sed -n "$(( (RUNS + 1) / 2 ))p")
  { grep -v "^$id"$'\t' "$TSV" || true; } > "$TSV.tmp"; mv "$TSV.tmp" "$TSV"
  printf "%s\t%s\t%s\n" "$id" "${t[*]}" "$med" | tee -a "$TSV"
done
rm -f "$TMP"
