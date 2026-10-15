"""Compares kbsrc against the Python reference on the Hydrogenaudio test files.

Exits with 77 (skipped) if any of the files is missing from --data; the ones present are still compared.
Prints a note where NumPy's long double is only float64 (as on Windows): the reference is then less precise.
"""
import os
import sys

import numpy as np

from common import Checker, kbsrc, parse_args, tmp, wavio
import kbsrc as ref  # reference/python, which common puts on sys.path

FILES = ["gaplesstest_m", "gaplesstest_s", "sweep_1_to_44KHz", "aliasing150db"]
TOL = {"fast": 1e-14, "exact": 1e-15}
SKIP = 77

args = parse_args(need_data=True)
c = Checker()
if np.finfo(np.longdouble).eps >= np.finfo(np.float64).eps:
    print("note  NumPy's long double is float64 here, so the reference runs its exact FFT and its Burg edge"
          " prediction in float64; AGENTS.md lists the checks that then exceed their limits", flush=True)
missing = []
for name in FILES:
    path = os.path.join(args.data, f"{name}-64bitfloat.wav")
    if not os.path.isfile(path):
        missing.append(name)
        print(f"skip  {name}: {path} not found", flush=True)
        continue
    x, sr, _ = wavio.read(path)
    for mode, tol in TOL.items():
        expected = ref.resample(x, sr, 44100, mode)
        kbsrc(args.bin, path, tmp("out.wav"), "--mode", mode)
        y = wavio.read(tmp("out.wav"))[0]
        err = np.abs(y - expected).max() if y.shape == expected.shape else np.inf
        c.check(err < tol, f"{name} {mode}: max diff {err:.2e} (limit {tol:.0e})")
if missing and not c.failed:
    print(f"skipped: {len(missing)} of {len(FILES)} test files not found in {args.data}")
    sys.exit(SKIP)
c.done()
