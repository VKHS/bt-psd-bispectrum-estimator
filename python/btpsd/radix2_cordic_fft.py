"""Radix-2 DIF FFT whose twiddle multiplications are integer CORDIC rotations (educational baseline).

This is the "single CORDIC processor" style of FFT that the paper's CRU-based
radix-8 design is compared against.  Every butterfly output is scaled by 1/2
(so the transform is scaled by 1/N and cannot overflow); the twiddle
``W_N^e`` is applied with :func:`btpsd.cordic.cordic_rotate_int` which reduces
the angle to the quadrant first (plain CORDIC diverges beyond ~99.9 degrees,
which was a defect of the original script).
"""

from __future__ import annotations

import numpy as np

from .cordic import cordic_rotate_int, twiddle_q
from .fixedpoint import rshift, saturate
from .r8fft import to_fixed
from .reference import frame_signal, get_window

__all__ = ["fft_cordic_fixed", "bt_psd_cordic"]


def _bitrev(n: int) -> np.ndarray:
    bits = int(np.log2(n))
    idx = np.arange(n)
    rev = np.zeros(n, dtype=np.int64)
    for b in range(bits):
        rev |= ((idx >> b) & 1) << (bits - 1 - b)
    return rev


def fft_cordic_fixed(re, im, wl: int = 16, n_iter: int = 14, frac: int = 14, rounding: str = "round"):
    """``FFT_N / N`` of integer words (last axis, ``N`` a power of two <= 1024); natural-order output."""
    xr = np.array(re, dtype=np.int64, copy=True)
    xi = np.array(im, dtype=np.int64, copy=True)
    n = xr.shape[-1]
    if n & (n - 1) or n > 1024:
        raise ValueError("N must be a power of two <= 1024")
    half = n // 2
    stage = 0
    while half >= 1:
        for blk in range(0, n, 2 * half):
            i = np.arange(half)
            a, b = blk + i, blk + i + half
            ur, ui, vr, vi = xr[..., a], xi[..., a], xr[..., b], xi[..., b]
            sr, si = rshift(ur + vr, 1, rounding), rshift(ui + vi, 1, rounding)
            dr, di = rshift(ur - vr, 1, rounding), rshift(ui - vi, 1, rounding)
            q = twiddle_q(i * (1 << stage), n)
            dr, di = cordic_rotate_int(dr, di, q, wl, n_iter, frac, rounding)
            xr[..., a], xi[..., a] = saturate(sr, wl), saturate(si, wl)
            xr[..., b], xi[..., b] = dr, di
        half //= 2
        stage += 1
    rev = _bitrev(n)
    return xr[..., rev], xi[..., rev]


def bt_psd_cordic(x, n: int = 256, wl: int = 16, window: str = "hann", overlap: float = 0.5):
    """Welch-style PSD with the CORDIC FFT, in the scale of ``reference.welch_twosided(scaling='paper')``."""
    x = np.asarray(x, dtype=float)
    hop = int(round(n * (1 - overlap)))
    frames = frame_signal(x, n, hop) * get_window(window, n)
    re, _ = to_fixed(frames, wl)
    Fr, Fi = fft_cordic_fixed(re, np.zeros_like(re), wl)
    P = np.mean(Fr.astype(float) ** 2 + Fi.astype(float) ** 2, axis=0)
    return P * n * n / float(1 << (2 * wl - 2)) / n
