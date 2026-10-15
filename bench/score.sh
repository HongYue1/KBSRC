#!/usr/bin/env bash
# Scores bench/out/<id> with the official Hydrogenaudio SRC Octave scripts.
# Needs GNU Octave with the signal and image packages, and tests/data/TestScripts.zip.
# usage: bench/score.sh [ids...]   (default: every converted id)
set -euo pipefail
cd "$(dirname "$0")/.."
ROOT=$PWD O=$ROOT/bench/out D=$ROOT/tests/data
TS=$D/TestScripts/TestScripts
[ -d "$TS" ] || { mkdir -p "$D/TestScripts"; unzip -q -o "$D/TestScripts.zip" -d "$D/TestScripts"; }
if [ -z "${OCTAVE:-}" ]; then
  OCTAVE=$(command -v octave || true)
  [ -n "$OCTAVE" ] ||
    OCTAVE=$(ls -d /c/Program\ Files/GNU\ Octave/Octave-*/mingw64/bin/octave.exe 2>/dev/null | tail -1)
fi
[ -x "$OCTAVE" ] || { echo "GNU Octave not found; set OCTAVE=/path/to/octave" >&2; exit 1; }
oct() {
  "$OCTAVE" --no-gui -q --eval \
    "pkg load signal; pkg load image; graphics_toolkit('qt'); set(0,'defaultfigurevisible','off'); $1"
}

REF=$O/_refs
if [ ! -f "$REF/intermodulation_sine-InternalUseAt44kHz.wav" ]; then
  mkdir -p "$REF"; cp bench/make_refs.m "$REF/"; (cd "$REF" && oct make_refs)
fi

IDS=("$@")
[ ${#IDS[@]} -gt 0 ] || IDS=($(cd "$O" && ls -d */ | tr -d / | grep -v '^_'))
S=$O/scores.tsv   # one row per id; rows for ids not scored in this run are kept
[ -f "$S" ] ||
  printf "id\tspectrogram\tbandwidth\timpulse_freq\taliasing\tintermod\tgapless\tdelay_samples\tpreringing\n" > "$S"
for id in "${IDS[@]}"; do
  W=$O/$id/_eval; rm -rf "$W"; mkdir -p "$W"
  cp "$TS"/TestScripts/*.m "$TS"/TestScripts/*REF.png "$W/"
  cp "$TS/TestSignals96KHzto44KHz/gaplesstest_c-64bitfloat-InternalUseAt44kHz.wav" \
    "$REF/intermodulation_sine-InternalUseAt44kHz.wav" "$W/"
  cp "$D/impulse-64bitfloat.wav" "$W/impulse-64bitfloat-InternalUse.wav"
  for f in aliasing150db impulse intermodulation_sine gaplesstest_s gaplesstest_m sweep_1_to_44KHz; do
    cp "$O/$id/$f-64bitfloat.wav" "$W/$f.wav"
  done
  for s in impulse aliasing150db intermodulation_harmonic_distortion gaplesstest \
           spectogram_sweep_1_to_44kHz_96kHz_srcto_44kHz; do
    (cd "$W" && oct "$s" > "log_$s.txt" 2>&1) || echo "$id: $s failed (see $W/log_$s.txt)" >&2
  done
  q() { [ -f "$W/$1" ] && tr -d ' \r' < "$W/$1" | grep . | tail -1 || echo NA; }  # last value
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$id" "$(q quality-spectrogram.txt)" "$(q quality-bandwidth.txt)" \
    "$(q quality-impulse_freq.txt)" "$(q quality-alias.txt)" "$(q quality-intermoddiff.txt)" \
    "$(q quality-gapless.txt)" "$(q calculateddelay.txt)" "$(q quality-preringing.txt)" > "$W/row.tsv"
  cat "$W/row.tsv"; { grep -v "^$id	" "$S" || true; cat "$W/row.tsv"; } > "$S.tmp" && mv "$S.tmp" "$S"
  rm -f "$W"/*.wav
done
