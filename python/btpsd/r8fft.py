"""Radix-8 512-point FFT models: floating point and bit-accurate fixed point.

Both follow :mod:`btpsd.schedule` (in place, 3 passes of 64 radix-8
butterflies).  A radix-8 butterfly is the 8-point DFT built from three
radix-2 DIF stages (outputs are produced in bit-reversed order and re-ordered
by fixed wiring), followed by one constant-rotator twiddle multiplication per
output ``k = 1..7``.

Fixed-point contract (shared with the RTL)
------------------------------------------
* words are ``wl``-bit two's complement, real and imaginary parts separately;
* every radix-2 add/subtract is followed by ``>> 1`` (rounding mode of the
  configuration), so a radix-8 butterfly scales by 1/8 and the whole FFT by
  ``1/512``; no overflow can occur in the add/subtract network;
* multiplications by ``W8^1 = exp(-j pi/4)`` and ``W8^3`` go through the same
  constant rotator (CRU ROM entry ``r = 128``); ``-j`` is a free swap/negate;
* the twiddle multiplication after each butterfly is
  :func:`btpsd.cordic.rotate_fixed`.
"""

from __future__ import annotations

import numpy as np

from .cordic import RotRom, rotate_fixed, rotate_float
from .fixedpoint import rshift, saturate
from .schedule import N, NPASS, output_index, pass_addresses, pass_twiddle_q

__all__ = ["dft8_float", "r8fft_float", "dft8_fixed", "r8fft_fixed", "to_fixed", "fft512_reference_scaled"]


# ---------------------------------------------------------------------------
# 8-point DFT (three radix-2 DIF stages, natural-order output)
# ---------------------------------------------------------------------------
_BITREV3 = np.array([0, 4, 2, 6, 1, 5, 3, 7])  # output slot j of the network holds X[_BITREV3[j]]


def dft8_float(x):
    """8-point DFT over the last axis (exact twiddles), natural-order output."""
    x = np.asarray(x, dtype=complex)
    s = x[..., :4] + x[..., 4:]
    d = (x[..., :4] - x[..., 4:]) * np.exp(-1j * np.pi / 4 * np.arange(4))
    out = np.empty_like(x)
    for grp, base in ((s, 0), (d, 1)):  # X[2m + base]
        c0, c1 = grp[..., 0] + grp[..., 2], grp[..., 1] + grp[..., 3]
        c2, c3 = grp[..., 0] - grp[..., 2], (grp[..., 1] - grp[..., 3]) * (-1j)
        out[..., 0 + base] = c0 + c1
        out[..., 4 + base] = c0 - c1
        out[..., 2 + base] = c2 + c3
        out[..., 6 + base] = c2 - c3
    return out


def r8fft_float(x):
    """512-point FFT via the in-place radix-8 schedule (float, exact twiddles); last axis."""
    x = np.asarray(x, dtype=complex)
    if x.shape[-1] != N:
        raise ValueError("last axis must have 512 samples")
    m = x.copy()
    for p in range(NPASS):
        adr = pass_addresses(p)  # (64, 8)
        tw = pass_twiddle_q(p)
        y = dft8_float(m[..., adr])
        m[..., adr] = rotate_float(y, tw)
    out = np.empty_like(m)
    out[..., output_index(np.arange(N))] = m
    return out


# ---------------------------------------------------------------------------
# bit-accurate fixed point
# ---------------------------------------------------------------------------
def to_fixed(x, wl: int, rounding: str = "round"):
    """Complex floats in [-1, 1) -> ``(re, im)`` int64 codes (``wl`` bits, ``wl-1`` fraction bits)."""
    x = np.asarray(x)
    sc = float(1 << (wl - 1))
    f = np.floor if rounding == "round" else (lambda v: np.floor(v))
    re = f(np.real(x) * sc + (0.5 if rounding == "round" else 0.0)).astype(np.int64)
    im = f(np.imag(x) * sc + (0.5 if rounding == "round" else 0.0)).astype(np.int64)
    return saturate(re, wl), saturate(im, wl)


def _addsub_half(ar, ai, br, bi, rounding, wl=None):
    """((a + b) >> 1, (a - b) >> 1) for complex integer pairs, saturated to ``wl`` bits when given.

    With round-half-up the difference of the two extreme codes, ``(2**(wl-1) - 1) - (-2**(wl-1))``, halves to
    ``+2**(wl-1)``, one LSB above the largest code; the hardware clamps that single case to ``2**(wl-1) - 1``.
    """
    sr = rshift(ar + br, 1, rounding)
    si = rshift(ai + bi, 1, rounding)
    dr = rshift(ar - br, 1, rounding)
    di = rshift(ai - bi, 1, rounding)
    if wl is not None:
        sr, si, dr, di = (saturate(v, wl) for v in (sr, si, dr, di))
    return sr, si, dr, di


def dft8_fixed(re, im, rom: RotRom, wl: int, rounding: str = "round"):
    """Scaled (1/8) 8-point DFT over the last axis; natural-order output ``(re, im)``."""
    re = np.asarray(re, dtype=np.int64)
    im = np.asarray(im, dtype=np.int64)
    # stage 1: distance 4
    sr, si, dr, di = _addsub_half(re[..., :4], im[..., :4], re[..., 4:], im[..., 4:], rounding, wl)
    # W8^i on the difference branch: i=0 identity, 1 -> q=128, 2 -> -j, 3 -> q=384
    q = np.array([0, 128, 256, 384])
    dr, di = rotate_fixed(dr, di, q, rom, wl, rounding)
    outr = np.empty_like(re)
    outi = np.empty_like(im)
    for (gr, gi), base in (((sr, si), 0), ((dr, di), 1)):
        # stage 2: distance 2 (c0,c1 sums; c2,c3 differences, c3 *= -j)
        c0r, c0i, c2r, c2i = _addsub_half(gr[..., 0], gi[..., 0], gr[..., 2], gi[..., 2], rounding, wl)
        c1r, c1i, e3r, e3i = _addsub_half(gr[..., 1], gi[..., 1], gr[..., 3], gi[..., 3], rounding, wl)
        c3r, c3i = saturate(e3i, wl), saturate(-e3r, wl)  # (-j)(a + jb) = b - ja
        # stage 3: distance 1
        for slot, (ur, ui, vr, vi) in ((0, (c0r, c0i, c1r, c1i)), (2, (c2r, c2i, c3r, c3i))):
            pr, pi_, mr, mi = _addsub_half(ur, ui, vr, vi, rounding, wl)
            outr[..., slot + base] = pr
            outi[..., slot + base] = pi_
            outr[..., slot + 4 + base] = mr
            outi[..., slot + 4 + base] = mi
    return outr, outi


def r8fft_fixed(re, im, rom: RotRom, wl: int, rounding: str = "round"):
    """Bit-accurate 512-point FFT scaled by 1/512 (natural-order output ``(re, im)`` int64).

    ``re, im``: integer arrays with 512 samples on the last axis (leading
    axes are independent chunks).
    """
    m_re = np.array(re, dtype=np.int64, copy=True)
    m_im = np.array(im, dtype=np.int64, copy=True)
    if m_re.shape[-1] != N:
        raise ValueError("last axis must have 512 samples")
    for p in range(NPASS):
        adr = pass_addresses(p)
        tw = pass_twiddle_q(p)
        yr, yi = dft8_fixed(m_re[..., adr], m_im[..., adr], rom, wl, rounding)
        yr, yi = rotate_fixed(yr, yi, tw, rom, wl, rounding)
        m_re[..., adr] = yr
        m_im[..., adr] = yi
    out_re = np.empty_like(m_re)
    out_im = np.empty_like(m_im)
    idx = output_index(np.arange(N))
    out_re[..., idx] = m_re
    out_im[..., idx] = m_im
    return out_re, out_im


def fft512_reference_scaled(re, im):
    """Exact ``FFT512(re + j im) / 512`` (float) - the target of :func:`r8fft_fixed`."""
    z = np.asarray(re, dtype=float) + 1j * np.asarray(im, dtype=float)
    return np.fft.fft(z, axis=-1) / 512.0
