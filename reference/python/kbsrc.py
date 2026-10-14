#!/usr/bin/env python3
"""KBSRC - Kaiser-Bessel sample-rate converter, NumPy reference implementation.

Design (see README):
  * Kaiser-windowed sinc kernel, sampled at the higher of the two rates. The window's main lobe makes the
    transition band 2*w Hz wide, w = sqrt(beta^2 + pi^2) / (2 pi T). The cutoff sits ceil(2*w)/2 Hz below
    the lower Nyquist frequency, so the transition band ends at or just below it.
  * Block overlap-save DFT resampling with block lengths that are integer multiples of the reduced ratio
    (exact phase, zero delay). Downsampling folds the bins above the output Nyquist back, exactly as sampling
    at the lower rate does; they lie in the stopband. Upsampling extends the input spectrum periodically.
  * Precision modes: 'fast' (float64 FFT) or 'exact' (80-bit extended FFT; output limited only by float64
    output rounding).
  * Burg-LPC edge extrapolation at file ends (no edge transients, gapless-friendly).
"""
import numpy as np, scipy.fft as F, math, sys, os, argparse
from fractions import Fraction
LD = np.longdouble
PI = LD('3.14159265358979323846264338327950288')
_KC, _HC = {}, {}

def _i0(x):
    x = np.asarray(x, LD); s = np.ones_like(x); t = np.ones_like(x); q = (x/2)**2
    for k in range(1, 500):
        t = t*q/(k*k); s += t
        if np.all(t < s*LD(1e-22)): break
    return s

def kernel_kw(fs, fc, T, beta):
    """classic Kaiser-windowed sinc, long double: one-sided taps"""
    key = ('kw', fs, fc, T, beta)
    if key in _KC: return _KC[key]
    L = int(math.floor(T*fs)); n = np.arange(0, L+1).astype(LD); x = 2*LD(fc)*n/LD(fs)
    with np.errstate(all='ignore'): s = np.where(n == 0, LD(1), np.sin(PI*x)/(PI*x))
    u = n/LD(L); w = _i0(LD(beta)*np.sqrt(np.clip(1-u*u, 0, None)))/_i0(LD(beta))
    h = 2*LD(fc)*s*w/LD(fs); _KC[key] = h; return h

def _spec(h, Nb, scale, T):
    key = (id(h), Nb, scale, T)
    if key not in _HC:
        L = len(h)-1; hk = np.zeros(Nb, LD); hk[:L+1] = h; hk[Nb-L:] = h[:0:-1]
        _HC[key] = (F.rfft(hk).real*LD(scale)).astype(T)
    return _HC[key]

def _burg(x, p):
    x = np.asarray(x, LD); f = x[1:].copy(); b = x[:-1].copy(); a = np.array([LD(1)])
    for _ in range(p):
        den = np.dot(f, f)+np.dot(b, b)
        if den <= 0: break
        k = -2*np.dot(f, b)/den; a = np.concatenate([a, [LD(0)]]); a = a+k*a[::-1]
        f, b = f[1:]+k*b[1:], b[:-1]+k*f[:-1]
    return a

def _extrap(x, n, order=32, train=4096):
    """Burg-LPC forward prediction of n samples after x"""
    seg = np.asarray(x[-train:], np.float64)
    if len(seg) < 2*order+8 or not np.any(seg): return np.zeros(n)
    c = -np.asarray(_burg(seg, order), np.float64)[1:][::-1]; p = len(c)
    buf = np.concatenate([seg[-p:], np.zeros(n)])
    for i in range(n): buf[p+i] = c @ buf[i:i+p]
    out = buf[p:]
    if not np.all(np.isfinite(out)) or np.max(np.abs(out)) > 4*np.max(np.abs(seg)): return np.zeros(n)
    return out

def resample(x, fi, fo, mode='exact', T=0.2, tb=40.0, nb=None, lpc=True, batch=6, workers=-1):
    """x: (n,) or (n, ch) float array -> (ceil(n*fo/fi), ch) float64"""
    x = np.asarray(x, np.float64); one = x.ndim == 1
    if one: x = x[:, None]
    r = Fraction(fo, fi); P, Q = r.numerator, r.denominator
    DT = LD if mode == 'exact' else np.float64; CT = np.clongdouble if mode == 'exact' else np.complex128
    up = fo > fi
    fk = fo if up else fi                                  # rate at which the kernel is sampled
    w = math.sqrt(tb*tb + 9.8696044010893586)/(2*3.14159265358979323846*T)   # main-lobe half-width, Hz
    drop = math.ceil(2*w)                                  # cutoff = (min(fi, fo) - drop) / 2, as in the C++ code
    if not min(fi, fo) - drop > 2*w: raise ValueError('length too short for these sample rates')
    h = kernel_kw(fk, (min(fi, fo) - drop)/2, T, tb); Lk = len(h)-1
    L = -(-Lk*fi//fk)+1                                    # kernel reach in input samples
    marg = max(-(-(L+Q)//Q), -(-(Lk+P)//P) if up else 0)   # reach in units of Q input / P output samples
    if nb is None: nb = (10 if mode == 'fast' else 20)*marg
    while F.next_fast_len(nb*Q) != nb*Q or F.next_fast_len(nb*P) != nb*P: nb += 1
    Nb, Mb = nb*Q, nb*P; hop = nb-2*marg
    nin, ch = x.shape; nout = -(-nin*P//Q)
    E = L+Q; pre = -(-(E+64)//Q)*Q
    nblk = -(-(pre+nin+E+64)//(hop*Q))
    xx = np.zeros((marg*Q+nblk*hop*Q+Nb, ch)); base = marg*Q
    xx[base+pre:base+pre+nin] = x
    if lpc:
        ext = E+32; fd = 0.5*(1+np.cos(np.pi*np.arange(32)/32))
        for c in range(ch):
            t = _extrap(x[:, c], ext); t[-32:] *= fd; hd = _extrap(x[::-1, c], ext)[::-1]; hd[:32] *= fd[::-1]
            xx[base+pre+nin:base+pre+nin+ext, c] = t; xx[base+pre-ext:base+pre, c] = hd
    K = Mb//2+1
    if up:   # kernel spectrum on the output grid; input spectrum extended periodically (exact images)
        H = _spec(h, Mb, Fraction(Mb, Nb), DT); kk = np.arange(K) % Nb
        src = np.where(kk <= Nb//2, kk, Nb-kk); cj = kk > Nb//2
    else:
        H = _spec(h, Nb, Fraction(Mb, Nb), DT)
        j = np.arange(K, min(len(H)-1, Mb-1)+1); j = j[H[j] != 0]
    out = np.empty((nblk*hop*P, ch))
    for c in range(ch):
        xc = np.ascontiguousarray(xx[:, c])
        for b0 in range(0, nblk, batch):
            bs = np.arange(b0, min(nblk, b0+batch))
            seg = np.lib.stride_tricks.sliding_window_view(xc, Nb)[b0*hop*Q:(b0+len(bs)-1)*hop*Q+1:hop*Q].astype(DT)
            Z = F.rfft(seg, axis=1, workers=workers)
            if up:
                Y = Z[:, src]; Y[:, cj] = np.conj(Y[:, cj]); Y *= H
                if Mb % 2 == 0: Y[:, -1] = Y[:, -1].real
            else:
                Z *= H
                Y = Z[:, :K].copy(); Y[:, Mb-j] += np.conj(Z[:, j])
                if Mb % 2 == 0: Y[:, -1] = 2*Z[:, Mb//2].real
            y = F.irfft(Y, Mb, axis=1, workers=workers)
            out[b0*hop*P:(b0+len(bs))*hop*P, c] = y[:, marg*P:marg*P+hop*P].reshape(-1)
    o = pre//Q*P; y = out[o:o+nout]
    return y[:, 0] if one else y

if __name__ == '__main__':
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import wavio
    ap = argparse.ArgumentParser(description='KBSRC resampler')
    ap.add_argument('inp'); ap.add_argument('out'); ap.add_argument('--rate', type=int, default=44100)
    ap.add_argument('--mode', choices=['exact', 'fast'], default='exact')
    ap.add_argument('--bits', type=int, default=64, choices=[32, 64])
    a = ap.parse_args()
    x, sr, _ = wavio.read(a.inp)
    wavio.write(a.out, resample(x, sr, a.rate, a.mode), a.rate, a.bits)
