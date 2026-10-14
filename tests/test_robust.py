"""Rejected options and malformed WAV headers must fail cleanly; valid neighbours must still convert.

Covers block sizes (an explicit --block at or below twice the kernel margin crashed the process, and huge blocks
overflowed the int transform sizes), kernel lengths (a --len shorter than one kernel sample was accepted and gave
nonfinite output) and the WAV header checks (a PCM header with 0 bits per sample was accepted, and so was an IEEE
float header whose sample size disagrees with the block alignment). A rejection must exit with status 1, print
"kbsrc: <reason>" and write no output file. Stdlib only: it does not import common.py, which needs numpy.

Usage: test_robust.py --bin <kbsrc>
"""
import argparse
import math
import os
import random
import re
import struct
import subprocess
import tempfile

ENGINES = {"default": None, "generic": "generic"}  # engine label -> KBSRC_CPU
MODES = ["exact", "fast", "ld"]
OUT_OF_RANGE = "option out of range"
TOO_LARGE = "block, rates or length too large"
PCM_WIDTH = "unsupported PCM width"
FLOAT_SIZE = "IEEE float sample size does not match the block alignment"


class Checker:
    def __init__(self):
        self.failed = 0

    def check(self, ok, label):
        print(("ok    " if ok else "FAIL  ") + label, flush=True)
        self.failed += not ok

    def done(self):
        if self.failed:
            raise SystemExit(f"{self.failed} check(s) failed")
        print("all checks passed")


def wav(tag, ch, rate, bits, align, data):
    """RIFF/WAVE bytes with a 16-byte fmt chunk; the header fields are written as given, valid or not."""
    fmt = struct.pack("<HHIIHH", tag, ch, rate, rate * align & 0xFFFFFFFF, align, bits)
    body = b"WAVE" + b"fmt " + struct.pack("<I", len(fmt)) + fmt + b"data" + struct.pack("<I", len(data)) + data
    return b"RIFF" + struct.pack("<I", len(body)) + body


def pcm(width, values):
    """Little-endian PCM samples in containers of width bytes (8-bit PCM is unsigned)."""
    if width == 1:
        return bytes(v + 128 for v in values)
    return b"".join((v & (1 << 8 * width) - 1).to_bytes(width, "little") for v in values)


def wav_samples(path):
    """The samples of a 64-bit float WAV written by kbsrc, or None."""
    with open(path, "rb") as f:
        blob = f.read()
    pos = 12
    while pos + 8 <= len(blob):
        tag, size = blob[pos:pos + 4], struct.unpack("<I", blob[pos + 4:pos + 8])[0]
        if tag == b"data":
            data = blob[pos + 8:pos + 8 + size]
            return struct.unpack("<%dd" % (len(data) // 8), data)
        pos += 8 + size + (size & 1)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True, help="path to the kbsrc executable")
    args = ap.parse_args()
    c = Checker()
    rng = random.Random(0x524F42)

    with tempfile.TemporaryDirectory(prefix="kbsrc-robust-") as tmp:
        out = os.path.join(tmp, "out.wav")

        def write(name, blob):
            path = os.path.join(tmp, name + ".wav")
            with open(path, "wb") as f:
                f.write(blob)
            return path

        def run(src, *opts, engine="default", timeout=120):
            """Convert src to out; returns (exit status or None on timeout, stderr)."""
            if os.path.exists(out):
                os.remove(out)
            env = dict(os.environ)
            env.pop("KBSRC_CPU", None)
            if ENGINES[engine]:
                env["KBSRC_CPU"] = ENGINES[engine]
            try:
                r = subprocess.run([args.bin, *map(str, opts), src, out], capture_output=True, text=True, env=env,
                                   timeout=timeout)
            except subprocess.TimeoutExpired:
                return None, f"timed out after {timeout} s"
            return r.returncode, r.stderr

        def rejected(label, reason, src, *opts, **kw):
            rc, err = run(src, *opts, **kw)
            ok = rc == 1 and err.startswith("kbsrc: ") and reason in err and not os.path.exists(out)
            c.check(ok, f"{label}: rejected, {reason}"
                    + ("" if ok else f" [rc={rc}, output {'written' if os.path.exists(out) else 'absent'}, "
                                     f"{err.strip()[:200]!r}]"))

        def converts(label, src, *opts, plan=None, **kw):
            """Accepted with finite output; plan = fragments the -v plan line must contain. Returns the samples."""
            rc, err = run(src, "-v", *opts, **kw)
            y = wav_samples(out) if rc == 0 and os.path.exists(out) else None
            ok = y is not None and len(y) > 0 and all(map(math.isfinite, y))
            missing = [p for p in plan or [] if not re.search(r"(^|\s)" + re.escape(p) + r"(\s|$)", err)]
            c.check(ok and not missing, f"{label}: converts, finite output" + (f", {' '.join(plan)}" if plan else "")
                    + ("" if ok and not missing else f" [rc={rc}, missing {missing}, {err.strip()[-200:]!r}]"))
            return y

        x = [rng.uniform(-0.5, 0.5) for _ in range(2 * 64)]
        src = write("f64_48k", wav(3, 2, 48000, 64, 16, struct.pack("<%dd" % len(x), *x)))

        # --block: at 48 kHz -> 96 kHz with the default --len 0.1 the kernel margin is 4802 frames. Blocks up to
        # 2 * 4802 = 9604 leave no hop; 9605 is the smallest accepted value and is rounded up to 9720 (7-smooth).
        for engine in ENGINES:
            for mode in MODES:
                for block in (1, 1000, 9604):
                    rejected(f"--block {block} {mode} [{engine}]", OUT_OF_RANGE, src, "--rate", 96000,
                             "--mode", mode, "--block", block, engine=engine)
                converts(f"--block 9605 {mode} [{engine}]", src, "--rate", 96000, "--mode", mode,
                         "--block", 9605, plan=["block=9720"], engine=engine)
        converts("automatic block", src, "--rate", 96000)
        converts("--block 20000", src, "--rate", 96000, "--block", 20000, plan=["block=20000"])
        # Transform sizes are block * P and block * Q (here P/Q = 2/1) and must stay well inside int: such blocks
        # used to run for minutes, fail to allocate or overflow.
        for block in (536870912, 2147483647):
            rejected(f"--block {block}", TOO_LARGE, src, "--rate", 96000, "--block", block, timeout=60)

        # --len is the one-sided kernel length in seconds, sampled at the higher rate (here 96 kHz):
        # 1e-5 s is 0.96 samples and leaves no kernel, 1.1e-5 s is one sample (3 taps).
        for engine in ENGINES:
            for mode in MODES:
                rejected(f"--len 1e-5 {mode} [{engine}]", OUT_OF_RANGE, src, "--rate", 96000, "--mode", mode,
                         "--len", "1e-5", engine=engine)
                converts(f"--len 1.1e-5 {mode} [{engine}]", src, "--rate", 96000, "--mode", mode,
                         "--len", "1.1e-5", plan=["taps=3"], engine=engine)
        for length in ("1e-8", "0", "-0.1", "10.5"):
            rejected(f"--len {length}", OUT_OF_RANGE, src, "--rate", 96000, "--len", length)
        # downsampling samples the kernel at the input rate (48 kHz): 2e-5 s is 0.96 samples
        rejected("--len 2e-5 48k -> 44.1k", OUT_OF_RANGE, src, "--rate", 44100, "--len", "2e-5")
        converts("--len 2.1e-5 48k -> 44.1k", src, "--rate", 44100, "--len", "2.1e-5", plan=["taps=3"])
        converts("--len 10 48k -> 16k (longest kernel)", src, "--rate", 16000, "--len", 10, plan=["taps=960001"])
        # A WAV header can carry any rate: at 2 GHz a 10 s kernel has too many samples for an int, and at the
        # prime 2147483647 Hz the transforms of P/Q = 48000/2147483647 cannot fit.
        hi = write("f64_2ghz", wav(3, 2, 2000000000, 64, 16, struct.pack("<%dd" % len(x), *x)))
        rejected("2 GHz input, --len 10", OUT_OF_RANGE, hi, "--rate", 48000, "--len", 10, timeout=60)
        prime = write("f64_prime", wav(3, 2, 2147483647, 64, 16, struct.pack("<%dd" % len(x), *x)))
        rejected("2147483647 Hz input", TOO_LARGE, prime, "--rate", 48000, timeout=60)

        # WAV headers. PCM is decoded by container size (block alignment / channels), so the declared bits per
        # sample must need exactly that many bytes; float samples are 4 or 8 bytes.
        for name, tag, bits, align, reason in (
                ("PCM 0 bits in 2 bytes", 1, 0, 2, PCM_WIDTH),
                ("PCM 8 bits in 2 bytes", 1, 8, 2, PCM_WIDTH),
                ("PCM 17 bits in 2 bytes", 1, 17, 2, PCM_WIDTH),
                ("PCM 33 bits in 4 bytes", 1, 33, 4, PCM_WIDTH),
                ("PCM 40 bits in 5 bytes", 1, 40, 5, PCM_WIDTH),
                ("float 64 bits in 4 bytes", 3, 64, 4, FLOAT_SIZE),
                ("float 32 bits in 8 bytes", 3, 32, 8, FLOAT_SIZE),
                ("float 32 bits in 2 bytes", 3, 32, 2, FLOAT_SIZE),
                ("float 16 bits in 2 bytes", 3, 16, 2, "unsupported WAV format")):
            bad = write("bad", wav(tag, 1, 48000, bits, align, bytes(align * 4096)))
            rejected(f"WAV {name}", reason, bad, "--rate", 96000)

        ints = [rng.randint(-(1 << 31), (1 << 31) - 1) for _ in range(2 * 64)]
        for name, tag, bits, width in (("u8", 1, 8, 1), ("s16", 1, 16, 2), ("s24", 1, 24, 3), ("s32", 1, 32, 4),
                                       ("f32", 3, 32, 4), ("f64", 3, 64, 8)):
            if tag == 3:
                data = struct.pack("<%d%s" % (len(x), "f" if width == 4 else "d"), *x)
            else:
                data = pcm(width, [v >> 32 - bits for v in ints])
            converts(f"WAV {name} stereo (control)", write(name, wav(tag, 2, 48000, bits, 2 * width, data)),
                     "--rate", 96000)
        # Fewer significant bits than the container are valid (left-justified, low bits zero) and decode exactly
        # like the full container width.
        for bits, width in ((12, 2), (20, 3)):
            data = pcm(width, [(v >> 32 - bits) << 8 * width - bits for v in ints])
            full = write("full", wav(1, 2, 48000, 8 * width, 2 * width, data))
            part = write("part", wav(1, 2, 48000, bits, 2 * width, data))
            y_full = converts(f"WAV PCM {8 * width} bits in {width} bytes", full, "--rate", 96000)
            y_part = converts(f"WAV PCM {bits} bits in {width} bytes", part, "--rate", 96000)
            c.check(y_full is not None and y_part == y_full,
                    f"WAV PCM {bits} bits in {width} bytes decodes like {8 * width} bits")
    c.done()


if __name__ == "__main__":
    main()
