"""Synthetic test signals (all *synthetic*; none is a recording and none is the paper's data).

Every generator returns a float array scaled to a peak of ``peak`` (default
0.9 of full scale) so that it can be quantised directly to a Q1.(wl-1) word.
"""

from __future__ import annotations

import numpy as np

__all__ = ["white", "pink", "ar2_lowpass", "sines_in_noise", "qpc", "eeg_like", "scale_peak", "SIGNALS"]


def scale_peak(x: np.ndarray, peak: float = 0.9) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    return x * (peak / np.max(np.abs(x)))


def white(n: int, seed: int = 0, peak: float = 0.9) -> np.ndarray:
    return scale_peak(np.random.default_rng(seed).standard_normal(n), peak)


def pink(n: int, seed: int = 0, peak: float = 0.9) -> np.ndarray:
    """1/f noise by spectral shaping of white noise."""
    w = np.random.default_rng(seed).standard_normal(n)
    W = np.fft.rfft(w)
    f = np.arange(len(W), dtype=float)
    f[0] = 1.0
    return scale_peak(np.fft.irfft(W / np.sqrt(f), n), peak)


def ar2_lowpass(n: int, seed: int = 0, peak: float = 0.9, pole_r: float = 0.97, pole_w: float = 0.35) -> np.ndarray:
    """Narrow-band AR(2) process (steep spectrum: strong dynamic range)."""
    from scipy import signal

    w = np.random.default_rng(seed).standard_normal(n + 512)
    a = [1.0, -2 * pole_r * np.cos(pole_w), pole_r**2]
    return scale_peak(signal.lfilter([1.0], a, w)[512:], peak)


def sines_in_noise(n: int, freqs=(0.05, 0.12, 0.31), amps=(1.0, 0.5, 0.25), snr_db: float = 20.0,
                   seed: int = 0, peak: float = 0.9) -> np.ndarray:
    """Sum of sinusoids (normalised frequencies, cycles/sample) plus white noise."""
    t = np.arange(n)
    rng = np.random.default_rng(seed)
    s = sum(a * np.cos(2 * np.pi * f * t + rng.uniform(0, 2 * np.pi)) for f, a in zip(freqs, amps))
    p = np.mean(s**2)
    noise = rng.standard_normal(n) * np.sqrt(p / 10 ** (snr_db / 10))
    return scale_peak(s + noise, peak)


def qpc(n: int, f1: float = 0.0625, f2: float = 0.203125, coupled: bool = True, snr_db: float = 25.0,
        seed: int = 0, peak: float = 0.9, block: int = 512):
    """Quadratic-phase-coupling test signal (synthetic).

    ``x = cos(2 pi f1 t + p1) + cos(2 pi f2 t + p2) + cos(2 pi (f1+f2) t + p3)`` plus white noise.
    The phases are redrawn independently every ``block`` samples (piecewise constant), so that
    non-overlapping segments of length ``block`` are independent realisations: for ``coupled=True``
    ``p3 = p1 + p2`` (phase-locked), otherwise ``p3`` is independent.  A stationary signal with fixed
    phases would give a bicoherence of 1 in both cases.  Defaults are bin-centred for
    FFT sizes that divide 512.  Returns ``(x, (f1, f2, f1 + f2))``.
    """
    rng = np.random.default_rng(seed)
    t = np.arange(n)
    nb = -(-n // block)
    ph = rng.uniform(0, 2 * np.pi, (nb, 3))
    if coupled:
        ph[:, 2] = ph[:, 0] + ph[:, 1]
    P = np.repeat(ph, block, axis=0)[:n]
    s = (np.cos(2 * np.pi * f1 * t + P[:, 0]) + np.cos(2 * np.pi * f2 * t + P[:, 1])
         + np.cos(2 * np.pi * (f1 + f2) * t + P[:, 2]))
    noise = rng.standard_normal(n) * np.sqrt(np.mean(s**2) / 10 ** (snr_db / 10))
    return scale_peak(s + noise, peak), (f1, f2, f1 + f2)


def eeg_like(n: int, fs: float = 256.0, seed: int = 0, peak: float = 0.9) -> np.ndarray:
    """Synthetic EEG-flavoured signal: 1/f background + 10 Hz alpha + 20 Hz beta, sampled at ``fs``."""
    t = np.arange(n) / fs
    bg = pink(n, seed, peak=1.0)
    rng = np.random.default_rng(seed + 1)
    x = bg + 0.6 * np.sin(2 * np.pi * 10 * t + rng.uniform(0, 6.28)) + 0.25 * np.sin(2 * np.pi * 20 * t + rng.uniform(0, 6.28))
    return scale_peak(x, peak)


def eeg_coupled(n: int, fs: float = 256.0, seed: int = 0, peak: float = 0.9, block: int = 512, f_alpha: float = 10.0) -> np.ndarray:
    """Synthetic EEG-flavoured signal with harmonic phase coupling: 1/f background + alpha tone + its phase-locked
    second harmonic (``p2 = 2 p1``) + a weaker independent beta tone.

    The alpha phase is redrawn every ``block`` samples (so non-overlapping 512-sample segments are independent
    realisations); within a block the harmonic stays locked, which produces a bispectrum peak at
    ``(f_alpha, f_alpha)`` in bins.  Synthetic, not a recording.
    """
    rng = np.random.default_rng(seed + 7)
    t = np.arange(n) / fs
    nb = -(-n // block)
    p1 = np.repeat(rng.uniform(0, 2 * np.pi, nb), block)[:n]
    p3 = np.repeat(rng.uniform(0, 2 * np.pi, nb), block)[:n]
    x = (pink(n, seed, peak=1.0) + 0.7 * np.cos(2 * np.pi * f_alpha * t + p1) + 0.5 * np.cos(2 * np.pi * 2 * f_alpha * t + 2 * p1)
         + 0.25 * np.cos(2 * np.pi * 1.5 * f_alpha * t + p3))
    return scale_peak(x, peak)


SIGNALS = {
    "white": white,
    "pink": pink,
    "ar2": ar2_lowpass,
    "sines": sines_in_noise,
    "eeg_like": eeg_like,
    "eeg_coupled": eeg_coupled,
}
