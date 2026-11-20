"""Chunk-merge PSD estimator of the paper: float and bit-accurate models.

Datapath (paper Figs. 3-5)
--------------------------
::

    x --> [512-pt chunks c_0, c_1, ...] --> R8 memory FFT --> A_s (scale 1/512)
    segment s = (c_s, c_{s+1}):   even bins  X(2m)   = (A_s(m) + A_{s+1}(m)) / 2     exact
                                  odd  bins  X(2m+1) = FD-filter{ A_s - A_{s+1} } / 2  (half-bin interpolation)
    --> frequency-domain window  Xw(k) = sum_d W(d) X(k-d)          (1/3/5 taps, Eq. 13)
    --> |Xw|^2 (sum-of-odds squarer) --> Save/Add accumulation over segments --> >> avg_shift

Consecutive segments share a chunk (50 % overlap of the 1024-sample
segments); every chunk is transformed once and its spectrum is reused by the
two segments that contain it (the ninth memory M_9 of the paper).

The quantity computed is the Welch average of windowed periodograms
``P(k) = (1/S) sum_s |DFT_1024{w x_s}(k)|^2 / 1024``.  With an ideal FD filter
(exact half-bin interpolation) and exact arithmetic, :func:`psd_float_fd`
reproduces :func:`btpsd.reference.welch_twosided` to rounding error; the
differences measured here are therefore *only* due to the 6-tap
interpolation (``mode='fd'``) and to fixed-point arithmetic (``psd_fixed``).
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np

from .cordic import PAPER_TOL, default_rom
from .fdfilter import PAPER_FD6, PAPER_FD6_START, merge_chunk_spectra
from .fixedpoint import best_csd_int, rshift, saturate
from .r8fft import r8fft_fixed, to_fixed
from .reference import welch_twosided
from .squarer import mag_sq
from .windows import apply_window_freq, spectral_kernel, approx_spectral_kernel, RAISED_COSINE

__all__ = [
    "HwConfig",
    "PsdResult",
    "psd_exact_welch",
    "psd_float_fd",
    "psd_fixed",
    "hw_power_scale",
    "chunk_view",
]

CH = 512  # chunk length
SEG = 1024  # segment length


@dataclass(frozen=True)
class HwConfig:
    """Parameters of the bit-accurate datapath (all constants are integers in the RTL)."""

    wl: int = 16  # data word length (memory words, merge and window outputs)
    rounding: str = "round"  # 'round' (half-up) or 'trunc'
    window: str = "hann"  # rect | hann | hamming | blackman
    fd_h: tuple = tuple(float(v) for v in PAPER_FD6)
    fd_start: int = PAPER_FD6_START
    fd_frac: int = 12  # fraction bits of the FD tap constants
    fd_max_terms: int | None = None  # CSD digits per tap (None: plain rounding)
    win_frac: int = 10  # fraction bits of the window-kernel constants
    win_max_terms: int | None = None
    rom_tol: float = PAPER_TOL  # CRU angle tolerance (see btpsd.cordic.generate_rom)
    avg_shift: int | None = None  # None: ceil(log2(S))
    squarer_units: int = 2

    # ---- integer constants ------------------------------------------------
    def fd_int(self) -> np.ndarray:
        h = np.asarray(self.fd_h, dtype=float) * (1 << self.fd_frac)
        if self.fd_max_terms is None:
            return np.floor(h + 0.5).astype(np.int64)
        return np.array([best_csd_int(v, self.fd_max_terms) for v in h], dtype=np.int64)

    def win_kernel(self) -> dict[int, float]:
        return spectral_kernel(self.window)

    def win_int(self) -> dict[int, int]:
        out = {}
        for d, c in self.win_kernel().items():
            v = c * (1 << self.win_frac)
            out[d] = int(np.floor(v + 0.5)) if self.win_max_terms is None else best_csd_int(v, self.win_max_terms)
        return out


@dataclass
class PsdResult:
    P_int: np.ndarray  # (1024,) accumulated, shifted power (integers)
    acc: np.ndarray  # (1024,) wide accumulator before the final shift
    S: int  # number of segments
    avg_shift: int
    cfg: HwConfig
    X_re: np.ndarray  # (S, 1024) windowed segment spectra (ints, scale 1/1024)
    X_im: np.ndarray
    squarer_cycles: int  # total sum-of-odds cycles over all bins and segments

    def to_float(self) -> np.ndarray:
        """Analytic conversion to the scale of ``welch_twosided(scaling='paper')`` (no fitting)."""
        return self.P_int.astype(float) * hw_power_scale(self.cfg.wl, self.S, self.avg_shift)


def hw_power_scale(wl: int, S: int, avg_shift: int) -> float:
    """Factor ``P_ref ~= P_int * factor`` between the integer datapath and the Welch float reference.

    ``|Xw_int| = |DFT(w x)| * 2**(wl-1) / 1024`` and ``P_ref = mean|DFT(w x)|^2 / 1024``
    give ``P_ref = 1024 * mean|Xw_int|^2 / 2**(2 wl - 2)`` with
    ``mean|Xw_int|^2 = P_int * 2**avg_shift / S``.
    """
    return 1024.0 * (2.0**avg_shift) / (S * 2.0 ** (2 * wl - 2))


def chunk_view(x, n_chunks: int | None = None) -> np.ndarray:
    """``(C, 512)`` view of the first ``C`` whole chunks of ``x``."""
    x = np.asarray(x)
    c = len(x) // CH if n_chunks is None else n_chunks
    return x[: c * CH].reshape(c, CH)


# ---------------------------------------------------------------------------
# float references
# ---------------------------------------------------------------------------
def psd_exact_welch(x, window: str = "hann") -> np.ndarray:
    """Reference: Welch average of windowed 1024-pt periodograms, hop 512, ``/L`` scaling."""
    x = np.asarray(x, dtype=float)
    x = x[: (len(x) // CH) * CH]
    return welch_twosided(x, L=SEG, hop=CH, window=window, periodic=True, scaling="paper")[0]


def psd_float_fd(
    x,
    window: str = "hann",
    h=None,
    start: int | None = None,
    mode: str = "fd",
    return_segments: bool = False,
):
    """Chunk-merge PSD in floating point.

    ``mode='exact'`` uses the exact half-bin interpolation (must equal
    :func:`psd_exact_welch`), ``mode='fd'`` the filter ``(h, start)`` (default:
    paper Eq. 6).  Returns ``P`` or ``(P, X_segments)``.
    """
    x = np.asarray(x, dtype=float)
    A = np.fft.fft(chunk_view(x), axis=-1)
    X = merge_chunk_spectra(A[:-1], A[1:], h=h, start=start, mode=mode)  # (S, 1024)
    if window.lower() not in ("rect", "rectangular", "boxcar"):
        Xw = apply_window_freq(X, spectral_kernel(window), axis=-1)
    else:
        Xw = X
    P = np.mean(np.abs(Xw) ** 2, axis=0) / SEG
    return (P, Xw) if return_segments else P


# ---------------------------------------------------------------------------
# bit-accurate model
# ---------------------------------------------------------------------------
def _merge_fixed(Ar, Ai, Br, Bi, cfg: HwConfig):
    """Segment spectrum (scale 1/1024) from two chunk spectra (scale 1/512); arrays (S, 512)."""
    rnd = cfg.rounding
    wl = cfg.wl
    # even bins
    Er = saturate(rshift(Ar + Br, 1, rnd), wl)
    Ei = saturate(rshift(Ai + Bi, 1, rnd), wl)
    # odd bins: FD filter on the sign-alternated difference (wide accumulator)
    Dr, Di = Ar - Br, Ai - Bi
    m = np.arange(CH)
    s = np.where(m % 2 == 0, 1, -1)
    Dpr, Dpi = Dr * s, Di * s
    hq = cfg.fd_int()
    accr = np.zeros_like(Dpr)
    acci = np.zeros_like(Dpi)
    for i, hi in enumerate(hq):
        sh = -(cfg.fd_start + i)
        accr = accr + int(hi) * np.roll(Dpr, sh, axis=-1)
        acci = acci + int(hi) * np.roll(Dpi, sh, axis=-1)
    # O = -j (-1)^m acc = (s*acci, -s*accr)
    Or_w, Oi_w = s * acci, -s * accr
    Or = saturate(rshift(Or_w, cfg.fd_frac + 1, rnd), wl)
    Oi = saturate(rshift(Oi_w, cfg.fd_frac + 1, rnd), wl)
    S = Ar.shape[0]
    Xr = np.empty((S, SEG), dtype=np.int64)
    Xi = np.empty((S, SEG), dtype=np.int64)
    Xr[:, 0::2], Xi[:, 0::2] = Er, Ei
    Xr[:, 1::2], Xi[:, 1::2] = Or, Oi
    return Xr, Xi


def _window_fixed(Xr, Xi, cfg: HwConfig):
    """Frequency-domain window ``Xw[k] = sum_d W[d] X[k-d]`` with integer taps, scale unchanged."""
    ker = cfg.win_int()
    if cfg.window.lower() in ("rect", "rectangular", "boxcar"):
        return Xr, Xi
    accr = np.zeros_like(Xr)
    acci = np.zeros_like(Xi)
    for d, c in ker.items():
        accr = accr + c * np.roll(Xr, d, axis=-1)
        acci = acci + c * np.roll(Xi, d, axis=-1)
    return (saturate(rshift(accr, cfg.win_frac, cfg.rounding), cfg.wl),
            saturate(rshift(acci, cfg.win_frac, cfg.rounding), cfg.wl))


def psd_fixed(x, cfg: HwConfig | None = None, quantize: bool = True) -> PsdResult:
    """Bit-accurate PSD of a real signal.

    ``x``: floats in ``[-1, 1)`` (quantised to ``wl`` bits, round-half-up) or,
    with ``quantize=False``, integer codes already in ``wl``-bit range.
    """
    cfg = cfg or HwConfig()
    x = np.asarray(x)
    if quantize:
        re, _ = to_fixed(x.astype(float), cfg.wl, cfg.rounding)
    else:
        re = x.astype(np.int64)
    ch = chunk_view(re)
    C = ch.shape[0]
    if C < 2:
        raise ValueError("need at least two 512-sample chunks")
    rom = default_rom(cfg.rom_tol)
    Ar, Ai = r8fft_fixed(ch, np.zeros_like(ch), rom, cfg.wl, cfg.rounding)  # (C, 512), scale 1/512
    Xr, Xi = _merge_fixed(Ar[:-1], Ai[:-1], Ar[1:], Ai[1:], cfg)
    Wr, Wi = _window_fixed(Xr, Xi, cfg)
    P, cyc = mag_sq(Wr, Wi, cfg.squarer_units)
    S = P.shape[0]
    shift = int(np.ceil(np.log2(S))) if cfg.avg_shift is None else cfg.avg_shift
    acc = P.sum(axis=0)
    return PsdResult(
        P_int=rshift(acc, shift, cfg.rounding),
        acc=acc,
        S=S,
        avg_shift=shift,
        cfg=cfg,
        X_re=Wr,
        X_im=Wi,
        squarer_cycles=int(cyc.sum()),
    )
