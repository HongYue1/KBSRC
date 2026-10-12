"""Sine-wave accuracy across common, odd and near-unity rate ratios."""
import numpy as np

from common import Checker, kbsrc, parse_args, tmp, wavio

PAIRS = [(96000, 44100), (44100, 48000), (48000, 44100), (96000, 48000), (44100, 96000), (192000, 44100),
         (44100, 88200), (8000, 44100), (44100, 8000), (48000, 192000),
         (44100, 44101), (44101, 44100), (48000, 44101), (96000, 32001), (22050, 48001), (48000, 48001),
         (44100, 44200)]
LIMIT_DB = -200.0  # the float64 sine reference itself limits this to about -215 dB

args = parse_args()
c = Checker()
rng = np.random.default_rng(1)
for fi, fo in PAIRS:
    n = fi * 2
    freqs, phases = [100, 1000, 0.4 * min(fi, fo)], rng.uniform(0, 6, 3)

    def sig(t):
        return sum(0.3 * np.sin(2 * np.pi * f * t + p) for f, p in zip(freqs, phases))

    x = sig(np.arange(n) / fi)
    wavio.write(tmp("in.wav"), np.stack([x, -x], 1), fi, 64)
    for mode in ["fast", "exact"]:
        kbsrc(args.bin, tmp("in.wav"), tmp("out.wav"), "--rate", fo, "--mode", mode)
        y, sr, _ = wavio.read(tmp("out.wav"))
        want = sig(np.arange(len(y)) / fo)
        mid = slice(len(y) // 4, 3 * len(y) // 4)
        err = max(np.abs(y[mid, 0] - want[mid]).max(), np.abs(y[mid, 1] + want[mid]).max())
        db = 20 * np.log10(err + 1e-300)
        ok = sr == fo and len(y) == -(-n * fo // fi) and db < LIMIT_DB
        c.check(ok, f"{fi} -> {fo} {mode}: {db:.1f} dB")

x = rng.standard_normal((4410, 2)) * 0.1
wavio.write(tmp("in.wav"), x, 44100, 64)
kbsrc(args.bin, tmp("in.wav"), tmp("out.wav"), "--rate", 44100)
c.check(np.array_equal(wavio.read(tmp("out.wav"))[0], x), "equal rates: bit-exact copy")
c.done()
