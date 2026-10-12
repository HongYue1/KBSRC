"""Compares kbsrc against the Python reference on the Hydrogenaudio test files."""
import os
import sys

import numpy as np

from common import Checker, kbsrc, parse_args, tmp, wavio

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "reference", "python"))
import kbsrc as ref  # noqa: E402

FILES = ["gaplesstest_m", "gaplesstest_s", "sweep_1_to_44KHz", "aliasing150db"]
TOL = {"fast": 1e-14, "exact": 1e-15}

args = parse_args(need_data=True)
if not os.path.isdir(args.data):
    print(f"skipped: test data not found in {args.data}")
    sys.exit(0)
c = Checker()
for name in FILES:
    path = os.path.join(args.data, f"{name}-64bitfloat.wav")
    x, sr, _ = wavio.read(path)
    for mode, tol in TOL.items():
        expected = ref.resample(x, sr, 44100, mode)
        kbsrc(args.bin, path, tmp("out.wav"), "--mode", mode)
        y = wavio.read(tmp("out.wav"))[0]
        err = np.abs(y - expected).max() if y.shape == expected.shape else np.inf
        c.check(err < tol, f"{name} {mode}: max diff {err:.2e} (limit {tol:.0e})")
c.done()
