"""WAV input/output formats, edge cases and command-line error handling."""
import struct
import wave

import numpy as np

from common import Checker, kbsrc, parse_args, tmp, wavio

args = parse_args()
c = Checker()
t = np.arange(48000) / 48000
s = 0.5 * np.sin(2 * np.pi * 440 * t)

# PCM input decoding (identity rate, so output equals the decoded input)
for width in [1, 2, 3, 4]:
    with wave.open(tmp("pcm.wav"), "wb") as w:
        w.setnchannels(2), w.setsampwidth(width), w.setframerate(48000)
        if width == 1:
            d = np.round(s * 127 + 128).astype(np.uint8).reshape(-1, 1)
        else:
            v = np.round(s * (2 ** (8 * width - 1) - 1)).astype(np.int64)
            d = np.stack([(v >> (8 * i)) & 255 for i in range(width)], 1).astype(np.uint8)
        w.writeframes(np.stack([d, d], 1).tobytes())
    kbsrc(args.bin, tmp("pcm.wav"), tmp("out.wav"), "--rate", 48000)
    y = wavio.read(tmp("out.wav"))[0]
    lsb = 2.0 ** -(8 * width - 1)
    c.check(np.abs(y[:, 0] - s).max() <= 1.5 * lsb and np.array_equal(y[:, 0], y[:, 1]), f"{8 * width}-bit PCM input")

# Output formats
wavio.write(tmp("in.wav"), s[:, None], 48000, 64)
kbsrc(args.bin, tmp("in.wav"), tmp("ref.wav"), "--rate", 44100)
ref = wavio.read(tmp("ref.wav"))[0]
for fmt, bits, tol in [("f32", 32, 2 ** -24), ("s32", 32, 2 ** -31), ("s24", 24, 2 ** -23), ("s16", 16, 2 ** -15)]:
    kbsrc(args.bin, tmp("in.wav"), tmp("out.wav"), "--rate", 44100, "--format", fmt)
    y, sr, (tag, b) = wavio.read(tmp("out.wav"))
    c.check(b == bits and sr == 44100 and np.abs(y - ref).max() <= tol, f"{fmt} output")

# Multichannel and tiny inputs
wavio.write(tmp("in.wav"), np.stack([s * (i + 1) / 7 for i in range(6)], 1), 48000, 64)
kbsrc(args.bin, tmp("in.wav"), tmp("out.wav"), "--rate", 96000, "--mode", "fast")
c.check(wavio.read(tmp("out.wav"))[0].shape == (96000, 6), "6 channels")
for n in [1, 10]:
    wavio.write(tmp("in.wav"), s[:n, None], 48000, 64)
    kbsrc(args.bin, tmp("in.wav"), tmp("out.wav"), "--rate", 44100)
    c.check(wavio.read(tmp("out.wav"))[0].shape == (-(-n * 44100 // 48000), 1), f"{n}-sample input")
hdr = b"WAVEfmt " + struct.pack("<IHHIIHH", 16, 1, 1, 48000, 96000, 2, 16) + b"data" + struct.pack("<I", 0)
with open(tmp("empty.wav"), "wb") as f:
    f.write(b"RIFF" + struct.pack("<I", len(hdr)) + hdr)
kbsrc(args.bin, tmp("empty.wav"), tmp("out.wav"), "--rate", 44100)
c.check(wavio.read(tmp("out.wav"))[0].shape == (0, 1), "empty input")

# Errors: bad usage exits 2, bad files exit 1, never a crash
with open(tmp("ref.wav"), "rb") as f, open(tmp("trunc.wav"), "wb") as g:
    g.write(f.read(30))
for argv, code in [([tmp("ref.wav"), tmp("o.wav"), "--rate", "0"], 2),
                   ([tmp("ref.wav"), tmp("o.wav"), "--rate", "abc"], 2),
                   ([tmp("ref.wav"), tmp("o.wav"), "--mode", "nope"], 2),
                   ([tmp("ref.wav"), tmp("o.wav"), "--len", "-1"], 1),
                   ([tmp("ref.wav")], 2),
                   ([tmp("missing.wav"), tmp("o.wav")], 1),
                   ([tmp("trunc.wav"), tmp("o.wav")], 1)]:
    r = kbsrc(args.bin, *argv, check=False)
    c.check(r.returncode == code and r.stderr.startswith("kbsrc: "), f"error exit {code}: {' '.join(argv[2:]) or argv[0][-12:]}")
c.done()
