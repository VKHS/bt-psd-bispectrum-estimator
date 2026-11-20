"""Floating-point reference estimators.

Two spectral estimators appear in the paper and both are provided here:

``welch_psd``  (modified-periodogram average, Eqs. 7-11)
    ``P = 1/K * sum_k (1/L) |FFT{w[n] x_k[n]}|^2`` over 50 %-overlapped
    segments.  *This is what the hardware datapath of Figs. 3-5 computes* when
    all arithmetic is exact, so it is the natural reference for SQNR studies.

``bt_correlogram_psd``  (Blackman-Tukey, Eq. 1)
    ``P_BT(w) = sum_{|k|<=M} r(k) w(k) exp(-jwk) = 1/(2 pi) P_per(w) * W(w)``,
    the Fourier transform of the lag-windowed autocorrelation, i.e. the
    periodogram smoothed in frequency by the lag-window transform.
    :func:`bt_via_periodogram_smoothing` evaluates the right-hand side of
    Eq. (1) and the tests verify the identity numerically.

The two coincide only in the loose sense that both smooth the periodogram; the
repository keeps them apart so that a result is never compared against the
wrong "ground truth".
"""

from __future__ import annotations

import numpy as np

from .windows import window as _window

__all__ = [
    "get_window",
    "frame_signal",
    "periodogram",
    "welch_twosided",
    "welch_psd",
    "autocorrelation_biased",
    "bt_correlogram_psd",
    "bt_via_periodogram_smoothing",
    "bt_psd_reference",
    "bispectrum_reference",
]


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def get_window(window_name: str, L: int, periodic: bool = True) -> np.ndarray:
    """Window of length ``L`` (periodic/DFT-even by default, see :mod:`btpsd.windows`)."""
    return _window(window_name, L, periodic=periodic)


def frame_signal(x: np.ndarray, frame_length: int, hop_length: int) -> np.ndarray:
    """Slice ``x`` into overlapping frames, shape ``(n_frames, frame_length)``.

    A signal shorter than one frame is zero-padded to one frame; trailing
    samples that do not fill a whole frame are dropped.
    """
    x = np.asarray(x)
    if x.ndim != 1:
        raise ValueError("x must be a 1-D array")
    if frame_length <= 0 or hop_length <= 0:
        raise ValueError("frame_length and hop_length must be > 0")
    if len(x) < frame_length:
        x = np.concatenate([x, np.zeros(frame_length - len(x), dtype=x.dtype)])
    n_frames = 1 + (len(x) - frame_length) // hop_length
    idx = np.arange(frame_length)[None, :] + hop_length * np.arange(n_frames)[:, None]
    return x[idx]


# ---------------------------------------------------------------------------
# periodogram family
# ---------------------------------------------------------------------------
def periodogram(x: np.ndarray, nfft: int | None = None, window: str = "rect") -> np.ndarray:
    """Two-sided modified periodogram ``|FFT{w x}|^2 / L`` of a single block."""
    x = np.asarray(x, dtype=float)
    L = len(x)
    w = get_window(window, L)
    return np.abs(np.fft.fft(x * w, n=nfft or L)) ** 2 / L


def welch_twosided(
    x: np.ndarray,
    L: int = 1024,
    hop: int | None = None,
    window: str = "hann",
    periodic: bool = True,
    scaling: str = "paper",
    fs: float = 1.0,
) -> tuple[np.ndarray, int]:
    """Averaged modified periodogram, full ``L`` bins; returns ``(P, K)``.

    scaling
        ``'paper'``    ``(1/K) sum |X|^2 / L``            (Eq. 8, window included)
        ``'density'``  ``(1/K) sum |X|^2 / (fs * sum w^2)``  (two-sided PSD)
    """
    x = np.asarray(x, dtype=float)
    hop = L // 2 if hop is None else hop
    frames = frame_signal(x, L, hop)
    w = get_window(window, L, periodic=periodic)
    X = np.fft.fft(frames * w, axis=-1)
    P = np.mean(np.abs(X) ** 2, axis=0)
    if scaling == "paper":
        P = P / L
    elif scaling == "density":
        P = P / (fs * np.sum(w**2))
    else:
        raise ValueError("scaling must be 'paper' or 'density'")
    return P, frames.shape[0]


def welch_psd(
    x: np.ndarray,
    fs: float = 1.0,
    N_fft: int = 1024,
    segment_length: int | None = None,
    overlap: float = 0.5,
    window: str = "hann",
    periodic: bool = True,
    scaling: str = "paper",
) -> tuple[np.ndarray, np.ndarray]:
    """One-sided ``(freqs, psd)`` of the averaged modified periodogram.

    ``segment_length`` defaults to ``N_fft``; shorter segments are zero-padded
    to ``N_fft`` before the FFT.
    """
    x = np.asarray(x, dtype=float)
    seg = N_fft if segment_length is None else segment_length
    hop = int(round(seg * (1.0 - overlap)))
    if hop <= 0:
        raise ValueError("overlap too large: hop length <= 0")
    frames = frame_signal(x, seg, hop)
    w = get_window(window, seg, periodic=periodic)
    X = np.fft.fft(frames * w, n=N_fft, axis=-1)
    P = np.mean(np.abs(X) ** 2, axis=0)
    if scaling == "paper":
        P = P / seg
    elif scaling == "density":
        P = P / (fs * np.sum(w**2))
    else:
        raise ValueError("scaling must be 'paper' or 'density'")
    half = N_fft // 2
    return np.arange(half + 1) * fs / N_fft, P[: half + 1]


# ---------------------------------------------------------------------------
# Blackman-Tukey (correlogram) estimator, Eq. (1)
# ---------------------------------------------------------------------------
def autocorrelation_biased(x: np.ndarray, max_lag: int) -> np.ndarray:
    """Biased autocorrelation ``r(k) = 1/N sum_n x[n] x[n+k]`` for ``k = 0..max_lag``."""
    x = np.asarray(x, dtype=float)
    N = len(x)
    nfft = 1 << int(np.ceil(np.log2(2 * N - 1)))
    X = np.fft.rfft(x, nfft)
    r = np.fft.irfft(np.abs(X) ** 2, nfft)[: max_lag + 1] / N
    return r


def _lag_window(name: str, M: int) -> np.ndarray:
    """Lag window ``w(k)``, ``k = -M..M`` (length ``2M+1``, w(0) = 1)."""
    n = name.lower()
    k = np.arange(-M, M + 1)
    if n in ("bartlett", "triangular"):
        return 1.0 - np.abs(k) / (M + 1.0)
    if n in ("hann", "hanning"):
        return 0.5 + 0.5 * np.cos(np.pi * k / (M + 1.0))
    if n == "hamming":
        return 0.54 + 0.46 * np.cos(np.pi * k / M)
    if n in ("rect", "rectangular", "boxcar"):
        return np.ones_like(k, dtype=float)
    raise ValueError(f"unsupported lag window {name!r}")


def bt_correlogram_psd(
    x: np.ndarray, M: int, lag_window: str = "bartlett", nfft: int = 1024
) -> np.ndarray:
    """Blackman-Tukey PSD ``sum_{k=-M}^{M} r(k) w(k) exp(-jwk)`` on ``nfft`` bins (two-sided)."""
    if nfft < 2 * M + 1:
        raise ValueError("nfft must be >= 2*M+1")
    r = autocorrelation_biased(x, M)
    w = _lag_window(lag_window, M)
    full = np.concatenate([r[:0:-1], r])  # r(-M..M), r(-k) = r(k) for real x
    seq = np.zeros(nfft)
    seq[: M + 1] = (full * w)[M:]  # k = 0..M
    seq[nfft - M :] = (full * w)[:M]  # k = -M..-1
    return np.real(np.fft.fft(seq))


def bt_via_periodogram_smoothing(
    x: np.ndarray, M: int, lag_window: str = "bartlett", nfft: int = 1024
) -> np.ndarray:
    """Right-hand side of Eq. (1): periodogram circularly convolved with ``W``.

    ``nfft`` must be >= ``2N-1`` for the identity with the correlogram to hold
    exactly (otherwise the periodogram itself is time-aliased).
    """
    x = np.asarray(x, dtype=float)
    N = len(x)
    if nfft < 2 * N - 1:
        raise ValueError("need nfft >= 2N-1 for the exact identity")
    Pper = np.abs(np.fft.fft(x, nfft)) ** 2 / N
    w = _lag_window(lag_window, M)
    seq = np.zeros(nfft)
    seq[: M + 1] = w[M:]
    seq[nfft - M :] = w[:M]
    W = np.real(np.fft.fft(seq))  # transform of the lag window
    # (1/2pi) P_per * W  ==  (1/nfft) circular convolution of the sampled spectra
    return np.real(np.fft.ifft(np.fft.fft(Pper) * np.fft.fft(W))) / nfft


# ---------------------------------------------------------------------------
# names used by the original scripts (kept for compatibility)
# ---------------------------------------------------------------------------
def bt_psd_reference(
    x: np.ndarray,
    fs: float = 1.0,
    N_fft: int = 1024,
    segment_length: int | None = None,
    overlap: float = 0.5,
    window: str = "hann",
) -> tuple[np.ndarray, np.ndarray]:
    """Alias of :func:`welch_psd` (the form computed by the paper's datapath)."""
    return welch_psd(
        x, fs=fs, N_fft=N_fft, segment_length=segment_length, overlap=overlap, window=window
    )


def bispectrum_reference(
    x: np.ndarray,
    fs: float = 1.0,
    N_fft: int = 256,
    segment_length: int | None = None,
    overlap: float = 0.5,
    window: str = "hann",
    max_freq_bins: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Direct triple-product bispectrum (Eq. 28); see :mod:`btpsd.bispectrum`."""
    from .bispectrum import bispectrum_direct

    seg = N_fft if segment_length is None else segment_length
    hop = int(round(seg * (1.0 - overlap)))
    max_freq_bins = N_fft // 4 if max_freq_bins is None else max_freq_bins
    B = bispectrum_direct(
        x, L=seg, hop=hop, window=window, nfft=N_fft, max_bin=max_freq_bins
    )
    freqs = np.fft.fftfreq(N_fft, d=1.0 / fs)[:max_freq_bins]
    return freqs, B
