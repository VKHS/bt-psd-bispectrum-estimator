# bt_reference.py
#
# Floating‑point reference implementations of:
#  - Blackman–Tukey style PSD (1024‑point)
#  - Bispectrum via triple‑product definition
#
# These follow the equations and concepts from:
#   "High‑performance power spectral/bispectral estimator ..." (Integration, 2024)
#   especially Eqs. (1), (8–13), (28–32). :contentReference[oaicite:2]{index=2}

import numpy as np


# ---------------------------
# Windowing & framing helpers
# ---------------------------

def get_window(window_name: str, L: int) -> np.ndarray:
    """
    Return a 1‑D window of length L.

    Supported names: 'rect', 'rectangular', 'hann', 'hamming', 'blackman'.
    """
    window_name = window_name.lower()
    if window_name in ("rect", "rectangular", "boxcar"):
        return np.ones(L, dtype=float)
    if window_name in ("hann", "hanning"):
        return np.hanning(L)
    if window_name in ("hamming",):
        return np.hamming(L)
    if window_name in ("blackman",):
        return np.blackman(L)
    raise ValueError(f"Unsupported window: {window_name}")


def frame_signal(x: np.ndarray, frame_length: int, hop_length: int) -> np.ndarray:
    """
    Slice a 1‑D signal into overlapping frames using striding.
    Returns an array of shape (n_frames, frame_length).

    frame_length: number of samples per frame
    hop_length:   step between consecutive frames
    """
    x = np.asarray(x)
    if x.ndim != 1:
        raise ValueError("x must be a 1‑D array")
    if frame_length <= 0 or hop_length <= 0:
        raise ValueError("frame_length and hop_length must be > 0")

    n = len(x)
    if n < frame_length:
        # Zero‑pad to at least one frame
        pad = frame_length - n
        x = np.concatenate([x, np.zeros(pad, dtype=x.dtype)])
        n = len(x)

    n_frames = 1 + (n - frame_length) // hop_length
    shape = (n_frames, frame_length)
    strides = (hop_length * x.strides[0], x.strides[0])
    frames = np.lib.stride_tricks.as_strided(
        x, shape=shape, strides=strides
    ).copy()
    return frames


# ---------------------------
# Blackman–Tukey PSD reference
# ---------------------------

def bt_psd_reference(
    x: np.ndarray,
    fs: float = 1.0,
    N_fft: int = 1024,
    segment_length: int | None = None,
    overlap: float = 0.5,
    window: str = "hann",
) -> tuple[np.ndarray, np.ndarray]:
    """
    Floating‑point Blackman–Tukey‑style PSD estimator.

    Conceptually implements:
        P_x(ω) = (1/L) * average_k |FFT{ w[n] x_k[n] }|^2
    using overlapping segments, Eq. (8)–(11) in the paper. :contentReference[oaicite:3]{index=3}

    Parameters
    ----------
    x : array_like
        Input 1‑D signal.
    fs : float
        Sampling frequency.
    N_fft : int
        FFT size for the PSD.
    segment_length : int or None
        Length L of each segment. If None, uses N_fft.
    overlap : float in [0,1)
        Fractional overlap between segments (0.5 → 50% overlap).
    window : str
        Window type ('rect', 'hann', 'hamming', 'blackman').

    Returns
    -------
    freqs : ndarray
        Frequency vector (one‑sided).
    psd : ndarray
        One‑sided PSD estimate.
    """
    x = np.asarray(x, dtype=float)
    if segment_length is None:
        segment_length = N_fft

    hop = int(round(segment_length * (1.0 - overlap)))
    if hop <= 0:
        raise ValueError("overlap too large: hop length <= 0")

    frames = frame_signal(x, segment_length, hop)
    w = get_window(window, segment_length)

    psd_accum = np.zeros(N_fft, dtype=float)
    for frame in frames:
        xf = frame * w
        X = np.fft.fft(xf, n=N_fft)
        # Eq. (8): (1/L) |X_i(e^{jω})|^2
        psd_seg = np.abs(X) ** 2 / segment_length
        psd_accum += psd_seg

    psd = psd_accum / frames.shape[0]

    freqs = np.fft.fftfreq(N_fft, d=1.0 / fs)
    # Return one‑sided spectrum
    half = N_fft // 2
    if N_fft % 2 == 0:
        idx = slice(0, half + 1)
    else:
        idx = slice(0, half + 1)
    return freqs[idx], psd[idx]


# ---------------------------
# Bispectrum reference
# ---------------------------

def bispectrum_reference(
    x: np.ndarray,
    fs: float = 1.0,
    N_fft: int = 256,
    segment_length: int | None = None,
    overlap: float = 0.5,
    window: str = "hann",
    max_freq_bins: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Simple third‑order bispectrum estimator using the triple‑product definition
    (Eq. (28) in the paper): B(f1,f2) = F(f1) F(f2) F*(f1+f2). :contentReference[oaicite:4]{index=4}

    This is O(max_freq_bins^2 * n_frames), used as a reference
    (not a high‑performance implementation).

    Parameters
    ----------
    x : array_like
        Input 1‑D signal.
    fs : float
        Sampling frequency.
    N_fft : int
        FFT size used per segment.
    segment_length : int or None
        Segment length; if None, uses N_fft.
    overlap : float
        Segment overlap (0.5 → 50%).
    window : str
        Window applied in time domain before FFT.
    max_freq_bins : int or None
        Number of low‑frequency bins (0..max_freq_bins‑1)
        for which to compute B(f1,f2). If None, defaults to N_fft//4.

    Returns
    -------
    freqs : ndarray
        Frequency vector for the f1,f2 axes (0..max_freq_bins‑1).
    B : ndarray, shape (max_freq_bins, max_freq_bins)
        Complex bispectrum estimate.
    """
    x = np.asarray(x, dtype=float)
    if segment_length is None:
        segment_length = N_fft

    hop = int(round(segment_length * (1.0 - overlap)))
    if hop <= 0:
        raise ValueError("overlap too large: hop length <= 0")

    frames = frame_signal(x, segment_length, hop)
    w = get_window(window, segment_length)

    if max_freq_bins is None:
        max_freq_bins = N_fft // 4
    max_freq_bins = min(max_freq_bins, N_fft)

    B = np.zeros((max_freq_bins, max_freq_bins), dtype=complex)

    for frame in frames:
        xf = frame * w
        X = np.fft.fft(xf, n=N_fft)

        for f1 in range(max_freq_bins):
            for f2 in range(max_freq_bins):
                k3 = f1 + f2
                if k3 >= N_fft:
                    continue
                B[f1, f2] += X[f1] * X[f2] * np.conj(X[k3])

    B /= frames.shape[0]
    freqs = np.fft.fftfreq(N_fft, d=1.0 / fs)[:max_freq_bins]
    return freqs, B


# Example quick check (remove or guard with if __name__ == "__main__": in real code)
if __name__ == "__main__":
    fs = 1000.0
    t = np.arange(0, 1.0, 1 / fs)
    x = np.sin(2 * np.pi * 50 * t) + 0.1 * np.random.randn(len(t))

    f_psd, P = bt_psd_reference(x, fs=fs, N_fft=1024, window="hann")
    f_bisp, B = bispectrum_reference(x, fs=fs, N_fft=128, max_freq_bins=16)
    print("PSD shape:", P.shape)
    print("Bispectrum shape:", B.shape)
