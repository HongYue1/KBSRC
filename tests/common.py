"""Shared helpers for the kbsrc test scripts."""
import argparse
import atexit
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "reference", "python"))
import wavio  # noqa: E402


def parse_args(need_data=False):
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True, help="path to the kbsrc executable")
    if need_data:
        ap.add_argument("--data", default=os.path.join(ROOT, "tests", "data"),
                        help="directory with the 64-bit float test WAVs")
    return ap.parse_args()


TMP = tempfile.mkdtemp(prefix="kbsrc-test-")
atexit.register(shutil.rmtree, TMP, ignore_errors=True)


def tmp(name):
    return os.path.join(TMP, name)


def kbsrc(binary, *args, check=True):
    r = subprocess.run([binary, *map(str, args)], capture_output=True, text=True)
    if check and r.returncode:
        raise SystemExit(f"kbsrc failed ({r.returncode}): {' '.join(map(str, args))}\n{r.stderr}")
    return r


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
