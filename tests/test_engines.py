"""Engine identity and rounding-mode independence of Exact.

The generic engine (Dekker's split) and the AVX2 engine (hardware FMA) must produce the same bits, also on
subnormal input, where the split's partial products can underflow, and Exact must not depend on the caller's
floating-point rounding mode (TwoSum and the product error are exact only in round-to-nearest).
KBSRC_CPU=generic selects the generic engine; without an AVX2 engine (other CPUs or builds) the cross-engine
comparisons are reported as skipped and the pinned values are still checked. Stdlib only: it does not import
common.py, which needs numpy.

Usage: test_engines.py --bin <kbsrc> [--driver <kbsrc_fenv_driver>]
"""
import argparse
import hashlib
import math
import os
import platform
import random
import re
import struct
import subprocess
import sys
import tempfile

# Pinned results, measured on the AVX2 engine: sha256 of the generated input WAVs, of the CLI output WAVs and of
# the driver's planar output, and the exact output samples of the one-frame case.
PRIMARY_WAV = "59c2fbaef9de2e52113f766058c3cc8d9cdac07fab7839c3666cf7600a34d5fe"
SUBNORMAL_WAV = "83d27bd8fb3438a12503cd1c0a9e78a1c980917ab4c808d226f14aca1b5a8f12"
CLI_PINS = [  # input, rate, --mode, sha256 of the output WAV, pinned only on x86-64 (x87 long double)
    ("primary", 96000, "exact", "355d557004d979d59358252264c63e689a136b51103c4390b2985a1c85f8398b", False),
    ("primary", 96000, "ld", "c55bae01612386ce25ad884afcc9f1223b19237d13debdeea508a0bc078edfdd", True),
    ("primary", 96000, "fast", "5c71eba51cfc4958af6736b2c4657dd3cb6bf291f643b23e22f3abe799d5cb07", False),
    ("primary", 44100, "exact", "752c3fe855e795ecbba79e83aaf00cfe47fa09c9d4905d2e3e61235197dc1ccc", False),
    ("subnormal", 96000, "exact", "90ae0178d5ae1fea6505f0c03c130a6ddb3a880d2b74ad21dd25cf53dc5e7fea", False),
    ("subnormal", 44100, "exact", "0a03683c6809085e3d1d2e926d21aa3faec9aecc3db66590f501c68440b7bef2", False),
]
SUBNORMAL_96K_PLANAR = "b1a721cde705af8f2a055dc45dc4575a46b29999f48da18ae4da590a04d8c7c3"
FRAME1 = (-0.5955104827880859, 0.11259841918945312)
FRAME1_OUT = ["bfe003180f1a9fbe", "bfd77c39344e9ed3", "3fb8386a3126e979", "3fb1c322eee114b9"]
MODES = ["nearest", "down", "up", "zero"]
X86_64 = platform.machine().lower() in ("x86_64", "amd64", "x64")


def primary97():
    """97 stereo frames of dyadic values in [-1, 1]."""
    rng = random.Random(0x4B42535243)
    return [rng.randint(-1048576, 1048576) / 1048576 for _ in range(194)]


def subnormal97():
    """97 stereo frames of subnormal doubles, |x| <= 2^-1024."""
    rng = random.Random(12345)
    return [math.ldexp(rng.randint(-1024, 1024) / 1024, -1024) for _ in range(194)]


def f64(samples):
    return struct.pack("<%dd" % len(samples), *samples)


def wav_f64(samples, ch=2, rate=48000):
    """64-bit float WAV with an 18-byte fmt chunk (cbSize = 0) and a fact chunk."""
    data = f64(samples)
    fmt = struct.pack("<HHIIHHH", 3, ch, rate, rate * 8 * ch, 8 * ch, 64, 0)
    body = (b"WAVE" + b"fmt " + struct.pack("<I", len(fmt)) + fmt + b"fact" + struct.pack("<II", 4, len(samples) // ch)
            + b"data" + struct.pack("<I", len(data)) + data)
    return b"RIFF" + struct.pack("<I", len(body)) + body


def sha(b):
    return hashlib.sha256(b).hexdigest()


def hex64(b):
    return [b[i:i + 8][::-1].hex() for i in range(0, len(b), 8)]


def read(path):
    """The bytes of path, or b"" if it does not exist."""
    if not os.path.exists(path):
        return b""
    with open(path, "rb") as f:
        return f.read()


class Checker:
    def __init__(self):
        self.failed = 0

    def check(self, ok, label):
        print(("ok    " if ok else "FAIL  ") + label, flush=True)
        self.failed += not ok

    def skip(self, label):
        print("skip  " + label, flush=True)

    def done(self):
        if self.failed:
            raise SystemExit(f"{self.failed} check(s) failed")
        print("all checks passed")


ENGINES = {"avx2": None, "generic": "generic"}  # requested engine -> KBSRC_CPU


def env_for(engine):
    env = dict(os.environ)
    env.pop("KBSRC_CPU", None)
    if ENGINES[engine]:
        env["KBSRC_CPU"] = ENGINES[engine]
    return env


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True, help="path to the kbsrc executable")
    ap.add_argument("--driver", help="path to kbsrc_fenv_driver (default: next to --bin)")
    args = ap.parse_args()
    driver = args.driver
    if driver is None:
        base = os.path.join(os.path.dirname(os.path.abspath(args.bin)), "kbsrc_fenv_driver")
        driver = next((p for p in (base, base + ".exe") if os.path.exists(p)), None)
    c = Checker()

    with tempfile.TemporaryDirectory(prefix="kbsrc-engines-") as tmp:
        inputs = {"primary": primary97(), "subnormal": subnormal97()}
        wavs = {}
        for name, want in (("primary", PRIMARY_WAV), ("subnormal", SUBNORMAL_WAV)):
            blob = wav_f64(inputs[name])
            c.check(sha(blob) == want, f"input {name}97 == {want[:16]}")
            wavs[name] = os.path.join(tmp, name + ".wav")
            with open(wavs[name], "wb") as f:
                f.write(blob)

        # CLI (production path, round-to-nearest): both engines give the pinned bytes
        seen = set()
        out = os.path.join(tmp, "out.wav")
        for name, rate, mode, want, x86_only in CLI_PINS:
            if x86_only and not X86_64:
                c.skip(f"cli {name} {rate} {mode}: pinned for x86-64 only")
                continue
            for engine in ENGINES:
                if os.path.exists(out):
                    os.remove(out)
                r = subprocess.run([args.bin, "--rate", str(rate), "--mode", mode, "-v", wavs[name], out],
                                   capture_output=True, text=True, env=env_for(engine))
                m = re.search(r"engine=(\w+)", r.stderr)
                isa = m.group(1) if m else "?"
                seen.add(isa)
                got = sha(read(out)) if r.returncode == 0 and os.path.exists(out) else f"rc={r.returncode}"
                c.check(got == want, f"cli {name} {rate} {mode} [{isa}] == {want[:16]}")
        cross = len(seen) > 1
        if not cross:
            c.skip(f"cross-engine identity: only the {', '.join(sorted(seen))} engine ran (no AVX2 engine here)")

        if not driver or not os.path.exists(driver):
            if args.driver:
                c.check(False, f"driver {args.driver} exists")
            else:
                c.skip("rounding-mode checks: kbsrc_fenv_driver not found next to --bin (pass --driver)")
            c.done()
            return

        raw = {}
        for name, samples in list(inputs.items()) + [("frame1", list(FRAME1))]:
            raw[name] = os.path.join(tmp, name + ".f64")
            with open(raw[name], "wb") as f:
                f.write(f64(samples))
        pout = os.path.join(tmp, "out.f64")

        def drive(engine, mode, name, rate_in, rate_out, prec="exact", length="0.1", threads=1):
            if os.path.exists(pout):
                os.remove(pout)
            r = subprocess.run([driver, mode, prec, str(rate_in), str(rate_out), "2", length, "40", str(threads),
                                raw[name], pout], capture_output=True, text=True, env=env_for(engine))
            info = dict(kv.split("=", 1) for kv in r.stdout.split() if "=" in kv)
            data = read(pout)
            return r.returncode, info, data

        # One stereo frame 8 kHz -> 16 kHz, sensitive to the rounding mode: the same samples in every mode
        for engine in ENGINES:
            for mode in MODES:
                rc, info, data = drive(engine, mode, "frame1", 8000, 16000, length="0.01")
                c.check(rc == 0 and hex64(data) == FRAME1_OUT and info.get("mode_kept") == "1",
                        f"exact frame1 8k->16k {mode} [{info.get('isa', '?')}] == pinned, caller's mode kept")

        # Subnormal input, 48k -> 96k: every mode and engine == the pinned round-to-nearest output
        for engine in ENGINES:
            for mode in MODES:
                rc, info, data = drive(engine, mode, "subnormal", 48000, 96000)
                c.check(rc == 0 and sha(data) == SUBNORMAL_96K_PLANAR and info.get("mode_kept") == "1",
                        f"exact subnormal97 48k->96k {mode} [{info.get('isa', '?')}] == pinned, caller's mode kept")

        # Normal input, up- and downsampling: every mode and engine == the round-to-nearest output
        for rate_out in (96000, 44100):
            rc0, _, ref = drive("avx2", "nearest", "primary", 48000, rate_out)
            for engine in ENGINES:
                for mode in MODES:
                    rc, info, data = drive(engine, mode, "primary", 48000, rate_out)
                    c.check(rc0 == 0 and rc == 0 and data == ref and info.get("mode_kept") == "1",
                            f"exact primary97 48k->{rate_out // 1000}k {mode} [{info.get('isa', '?')}] == nearest, "
                            "caller's mode kept")

        # Worker threads may start in round-to-nearest while the caller is not: 1 and 8 threads == nearest
        for engine in ENGINES:
            rc0, info0, ref = drive(engine, "nearest", "primary", 48000, 96000)
            for threads in (1, 8):
                rc, info, data = drive(engine, "down", "primary", 48000, 96000, threads=threads)
                c.check(rc0 == 0 and rc == 0 and data == ref and info.get("mode_kept") == "1",
                        f"exact primary97 48k->96k down, {threads} thread(s) [{info.get('isa', '?')}] == nearest")

        # the caller's mode also survives a rejected call, and Fast/LongDouble leave it alone
        rc, info, _ = drive("avx2", "down", "frame1", 0, 16000)
        c.check(rc == 3 and info.get("threw") == "1" and info.get("mode_kept") == "1",
                "exact rejected argument under down: invalid_argument, caller's mode kept")
        for prec in ("fast", "ld"):
            rc, info, _ = drive("avx2", "up", "frame1", 8000, 16000, prec=prec, length="0.01")
            c.check(rc == 0 and info.get("mode_kept") == "1", f"{prec} under up: caller's mode kept")
    c.done()


if __name__ == "__main__":
    main()
