"""No aliasing or imaging: nothing at or above the lower Nyquist frequency passes the filter.

Measures the converter's continuous impulse response the way the Hydrogenaudio impulse test does: Q unit impulses
(P/Q the reduced rate ratio), one at each input position modulo Q, cover every phase between the input and the
output grid, so their responses interleave into the response sampled at lcm(rate_in, rate_out). Its spectrum must
be flat in the passband and below LIMIT_DB from the lower Nyquist frequency up: when downsampling such content
would fold back below the output Nyquist frequency, when upsampling it would be an image of the input.
"""
import math
import re

import numpy as np

from common import Checker, kbsrc, parse_args, tmp, wavio

# fi, fo: down- and upsampling by the Hydrogenaudio ratio, and near-unity ratios (general fold path)
PAIRS = [(96000, 44100), (44100, 96000), (1001, 1000), (1000, 1001)]
LENGTH = 0.2          # kernel half-length, seconds (the default)
LIMIT_DB = -280.0     # the stopband of the default filter is near -345 dB, float64 output noise near -320 dB
FLAT = 1e-12          # passband deviation limit, relative

args = parse_args()
c = Checker()
for fi, fo in PAIRS:
    g = math.gcd(fi, fo)
    P, Q, R = fo // g, fi // g, fi * (fo // g)       # R = lcm(fi, fo): rate of the interleaved response
    reach = math.ceil(LENGTH * fi) + 1                 # kernel reach in input samples
    D = Q * math.ceil((2 * reach + 64) / Q)            # impulse spacing, a multiple of Q
    m = D * np.arange(1, Q + 1) + np.arange(Q)         # positions: m mod Q = 0, 1, ..., Q - 1
    x = np.zeros(D * (Q + 2))
    x[m] = 1.0
    wavio.write(tmp("in.wav"), x, fi, 64)
    W = reach * P                                      # half-width of the response on the R grid
    for mode in ["fast", "exact"]:
        r = kbsrc(args.bin, tmp("in.wav"), tmp("out.wav"), "--rate", fo, "--mode", mode, "--len", LENGTH, "--no-lpc",
                  "-v")
        y = wavio.read(tmp("out.wav"))[0][:, 0]
        res = np.zeros(2 * W + 1)
        for mi in m:
            k = np.arange(-(-(mi * P - W) // Q), (mi * P + W) // Q + 1)
            res[k * Q - mi * P + W] += y[k]
        H = np.abs(np.fft.rfft(res))
        H /= H[0]
        f = np.arange(len(H)) * R / len(res)
        edge = re.search(r"passband 0-([0-9.]+) Hz, stopband from ([0-9.]+) Hz", r.stderr)
        passband, stopband = (float(edge.group(1)), float(edge.group(2))) if edge else (0.0, math.inf)
        stop_db = 20 * np.log10(H[f >= min(fi, fo) / 2].max() + 1e-300)
        flat = np.abs(H[f <= passband] - 1).max()
        label = f"{fi} -> {fo} {mode}"
        c.check(edge is not None and stopband <= min(fi, fo) / 2 and passband > 0,
                f"{label}: passband 0-{passband} Hz, stopband from {stopband} Hz <= {min(fi, fo) / 2} Hz")
        c.check(stop_db < LIMIT_DB, f"{label}: from {min(fi, fo) / 2} Hz up {stop_db:.1f} dB (limit {LIMIT_DB:.0f})")
        c.check(flat < FLAT, f"{label}: passband flat to {flat:.1e} (limit {FLAT:.0e})")
c.done()
