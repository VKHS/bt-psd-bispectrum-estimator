"""Bispectrum estimators: direct (FFT triple product), indirect (third-order cumulant), bicoherence,
and a bit-accurate fixed-point direct estimator that reuses the radix-8 memory FFT.

Definitions (two-sided frequency indices ``f1, f2`` in DFT bins)::

    direct      B(f1, f2) = 1/S  sum_s  F_s(f1) F_s(f2) F_s*(f1 + f2) / L          (Eq. 28)
    indirect    B(f1, f2) = sum_{t1,t2} c3(t1, t2) w(t1) w(t2) w(t1 - t2) exp(-j 2 pi (f1 t1 + f2 t2)/nfft)
    bicoherence b^2(f1, f2) = |sum_s F1 F2 F3*|^2 / ( sum_s |F1 F2|^2 * sum_s |F3|^2 )   in [0, 1]

``F_s`` is the DFT of the windowed segment ``s``; ``c3(t1, t2) = E[x(n) x(n+t1) x(n+t2)]``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .cordic import PAPER_TOL, default_rom
from .fixedpoint import rshift, saturate
from .r8fft import r8fft_fixed, to_fixed
from .reference import frame_signal, get_window
from .windows import spectral_kernel

__all__ = [
    "bispectrum_direct",
    "third_order_cumulant",
    "bispectrum_indirect",
    "bicoherence",
    "BispecConfig",
    "bispectrum_fixed",
    "bispectrum_hw_to_float",
]


def _segment_spectra(x, L, hop, window, nfft):
    x = np.asarray(x, dtype=float)
    frames = frame_signal(x, L, hop)
    w = get_window(window, L)
    return np.fft.fft(frames * w, n=nfft, axis=-1)


def _triple(F, max_bin, nfft):
    f = np.arange(max_bin)
    F1 = F[:, f][:, :, None]
    F2 = F[:, f][:, None, :]
    F3 = F[:, (f[:, None] + f[None, :]) % nfft]
    return F1 * F2 * np.conj(F3)  # (S, max_bin, max_bin)


def bispectrum_direct(x, L: int = 256, hop: int | None = None, window: str = "hann",
                      nfft: int | None = None, max_bin: int | None = None) -> np.ndarray:
    """Direct bispectrum ``B[f1, f2]``, ``f1, f2 = 0..max_bin-1`` (``max_bin`` defaults to ``nfft//4``)."""
    hop = L // 2 if hop is None else hop
    nfft = nfft or L
    max_bin = nfft // 4 if max_bin is None else max_bin
    F = _segment_spectra(x, L, hop, window, nfft)
    return np.mean(_triple(F, max_bin, nfft), axis=0) / L


def bicoherence(x, L: int = 256, hop: int | None = None, window: str = "hann",
                nfft: int | None = None, max_bin: int | None = None) -> np.ndarray:
    """Squared bicoherence (Kim & Powers normalisation) on ``f1, f2 < max_bin``."""
    hop = L // 2 if hop is None else hop
    nfft = nfft or L
    max_bin = nfft // 4 if max_bin is None else max_bin
    F = _segment_spectra(x, L, hop, window, nfft)
    f = np.arange(max_bin)
    T = _triple(F, max_bin, nfft)
    num = np.abs(np.sum(T, axis=0)) ** 2
    d12 = np.sum(np.abs(F[:, f][:, :, None] * F[:, f][:, None, :]) ** 2, axis=0)
    d3 = np.sum(np.abs(F[:, (f[:, None] + f[None, :]) % nfft]) ** 2, axis=0)
    return num / (d12 * d3)


def third_order_cumulant(x, maxlag: int = 24, L: int = 256, hop: int | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Segment-averaged third-order moment ``c3[t1, t2] = E[x(n) x(n+t1) x(n+t2)]`` of the mean-removed segments.

    Returns ``(lags, c3)`` with ``lags = -maxlag .. maxlag`` and ``c3`` of shape ``(2 maxlag + 1, 2 maxlag + 1)``.
    """
    hop = L // 2 if hop is None else hop
    frames = frame_signal(np.asarray(x, dtype=float), L, hop)
    frames = frames - frames.mean(axis=1, keepdims=True)
    M = maxlag
    lags = np.arange(-M, M + 1)
    c3 = np.zeros((2 * M + 1, 2 * M + 1))
    for i1, t1 in enumerate(lags):
        for i2, t2 in enumerate(lags):
            lo = max(0, -t1, -t2)
            hi = min(L, L - t1, L - t2)
            if hi <= lo:
                continue
            seg = frames[:, lo:hi] * frames[:, lo + t1:hi + t1] * frames[:, lo + t2:hi + t2]
            c3[i1, i2] = seg.mean()
    return lags, c3


def bispectrum_indirect(x, maxlag: int = 24, L: int = 256, hop: int | None = None, nfft: int = 256,
                        max_bin: int | None = None) -> np.ndarray:
    """Indirect bispectrum: third-order cumulant -> lag window -> 2-D DFT (Parzen lag window)."""
    from scipy import signal as sg

    max_bin = nfft // 4 if max_bin is None else max_bin
    M = maxlag
    lags, c3 = third_order_cumulant(x, maxlag, L, hop)
    pw = sg.get_window("parzen", 2 * M + 1)  # symmetric lag window, w(0) = 1
    t1g, t2g = np.meshgrid(np.arange(2 * M + 1), np.arange(2 * M + 1), indexing="ij")
    idx_d = np.clip(t1g - t2g + M, 0, 2 * M)
    W = pw[t1g] * pw[t2g] * np.where(np.abs(t1g - t2g) <= M, pw[idx_d], 0.0)
    c = c3 * W
    # place lag (t1, t2) at array index (t1 mod nfft, t2 mod nfft) and take the 2-D DFT
    big = np.zeros((nfft, nfft))
    big[np.ix_(lags % nfft, lags % nfft)] = c
    B = np.fft.fft2(big)
    return B[:max_bin, :max_bin]


# ---------------------------------------------------------------------------
# bit-accurate direct estimator
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class BispecConfig:
    wl: int = 16
    rounding: str = "round"
    window: str = "hann"
    hop: int = 256  # 512-sample segments, 50 % overlap
    win_frac: int = 10
    max_bin: int = 128
    rom_tol: float = PAPER_TOL
    avg_shift: int | None = None


def bispectrum_fixed(x, cfg: BispecConfig | None = None):
    """Bit-accurate ``B[f1, f2]`` on 512-sample segments.

    Per segment: R8 memory FFT (scale 1/512), frequency-domain window (integer
    taps), ``P12 = F1 F2`` rounded back to ``wl`` bits (fractional multiply),
    ``F1 F2 F3*`` exact, accumulated over segments in a wide register and
    shifted by ``avg_shift``.  Returns ``(B_re, B_im, S, avg_shift)``.
    """
    cfg = cfg or BispecConfig()
    wl = cfg.wl
    re, _ = to_fixed(np.asarray(x, dtype=float), wl, cfg.rounding)
    frames = frame_signal(re, 512, cfg.hop)
    S = frames.shape[0]
    rom = default_rom(cfg.rom_tol)
    Fr, Fi = r8fft_fixed(frames, np.zeros_like(frames), rom, wl, cfg.rounding)
    if cfg.window.lower() not in ("rect", "rectangular", "boxcar"):
        accr = np.zeros_like(Fr)
        acci = np.zeros_like(Fi)
        for d, c in spectral_kernel(cfg.window).items():
            ci = int(np.floor(c * (1 << cfg.win_frac) + 0.5))
            accr = accr + ci * np.roll(Fr, d, axis=-1)
            acci = acci + ci * np.roll(Fi, d, axis=-1)
        Fr = saturate(rshift(accr, cfg.win_frac, cfg.rounding), wl)
        Fi = saturate(rshift(acci, cfg.win_frac, cfg.rounding), wl)
    f = np.arange(cfg.max_bin)
    s = (f[:, None] + f[None, :]) % 512
    ar, ai = Fr[:, f][:, :, None], Fi[:, f][:, :, None]
    br, bi = Fr[:, f][:, None, :], Fi[:, f][:, None, :]
    pr = rshift(ar * br - ai * bi, wl - 1, cfg.rounding)
    pi_ = rshift(ar * bi + ai * br, wl - 1, cfg.rounding)
    cr, ci_ = Fr[:, s], Fi[:, s]
    tr = pr * cr + pi_ * ci_  # P12 * conj(F3)
    ti = pi_ * cr - pr * ci_
    accr, acci = tr.sum(axis=0), ti.sum(axis=0)
    shift = int(np.ceil(np.log2(S))) if cfg.avg_shift is None else cfg.avg_shift
    return rshift(accr, shift, cfg.rounding), rshift(acci, shift, cfg.rounding), S, shift


def bispectrum_hw_to_float(B_re, B_im, wl: int, S: int, avg_shift: int) -> np.ndarray:
    """Analytic scale of :func:`bispectrum_fixed` to ``bispectrum_direct(L=512, nfft=512)``.

    ``B_ref = 512**2 * mean_s(F1u F2u F3u*)`` with ``F_u = F_int / 2**(wl-1)`` and
    ``mean = B_int * 2**avg_shift / (S * 2**(2 wl - 2))``.
    """
    return (B_re + 1j * B_im) * (512.0**2) * (2.0**avg_shift) / (S * 2.0 ** (2 * wl - 2))
