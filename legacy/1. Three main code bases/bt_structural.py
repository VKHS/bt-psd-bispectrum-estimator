# bt_structural.py
#
# Structural Blackman–Tukey PSD model inspired by the hardware architecture:
#  - 512‑sample chunks with 50% overlap
#  - Bidirectional 6‑tap fractional delay filter (Eq. (6))
#  - Even / odd combination of two chunks (Eqs. (3–4))
#  - 1024‑point FFT for PSD estimation
#
# This is a high‑level *software* model; actual hardware does the
# combination after a memory‑based 512‑point FFT, but the structure
# (chunks, FD, even/odd interleave, windowing) matches the paper. :contentReference[oaicite:6]{index=6}

import numpy as np
from bt_reference import get_window, frame_signal


# 6‑tap fractional delay coefficients from Eq. (6) and Fig. 9. :contentReference[oaicite:7]{index=7}
FD_COEFFS_6TAP = np.array(
    [0.125, -0.2122, 0.6366, 0.6366, -0.2122, 0.125],
    dtype=float,
)


def bidirectional_fd_filter(x: np.ndarray, h: np.ndarray = FD_COEFFS_6TAP) -> np.ndarray:
    """
    Apply a symmetric bidirectional FIR fractional‑delay filter.

    This models the bidirectional FD filter chosen in the paper (Table 3,
    Fig. 9), which gives nearly constant group delay and better SNR than
    unidirectional variants. :contentReference[oaicite:8]{index=8}
    """
    x = np.asarray(x, dtype=float)

    # Forward pass
    y_fwd = np.convolve(x, h, mode="same")
    # Backward pass for bidirectional (zero‑phase) response
    y_bwd = np.convolve(x[::-1], h, mode="same")[::-1]

    return 0.5 * (y_fwd + y_bwd)


def bt_psd_structural_2x512(
    x: np.ndarray,
    fs: float = 1.0,
    window: str = "hann",
) -> tuple[np.ndarray, np.ndarray]:
    """
    Structural 1024‑point BT‑style PSD estimator using:
      - 512‑sample chunks with 50% overlap,
      - 6‑tap bidirectional FD filter (Eq. (6)), and
      - even/odd interleaving of two chunks, following Eqs. (3–4):

          X(2n)   = x1(n) + x2(n)
          X(2n+1) = x1(n+1/2) − x2(n+1/2)     (informally)

    We approximate this by:
      even[n] = x1[n] + x2[n]
      odd[n]  = FD{x1}[n] − FD{x2}[n]
      seg[2n]   = even[n]
      seg[2n+1] = odd[n]

    Then apply a 1024‑point window and FFT to compute PSD.

    Parameters
    ----------
    x : array_like
        Input 1‑D signal.
    fs : float
        Sampling frequency.
    window : str
        Window type ('rect', 'hann', 'hamming', 'blackman').

    Returns
    -------
    freqs : ndarray
        One‑sided frequency vector.
    psd : ndarray
        One‑sided PSD estimate (structural BT).
    """
    x = np.asarray(x, dtype=float)

    chunk_len = 512
    hop = chunk_len // 2  # 50% overlap as in Fig. 4a. :contentReference[oaicite:9]{index=9}

    chunks = frame_signal(x, chunk_len, hop)  # shape: (n_chunks, 512)
    n_chunks = chunks.shape[0]
    if n_chunks < 2:
        raise ValueError("Signal too short: need at least 2 chunks for 2×512 structural BT.")

    N_fft = 1024
    w1024 = get_window(window, N_fft)

    psd_accum = np.zeros(N_fft, dtype=float)
    n_segments = 0

    for i in range(n_chunks - 1):
        x1 = chunks[i]
        x2 = chunks[i + 1]

        # Even path: direct sum of current and next chunk
        even = x1 + x2

        # Odd path: half‑sample delayed difference using bidirectional FD filter
        x1_fd = bidirectional_fd_filter(x1)
        x2_fd = bidirectional_fd_filter(x2)
        odd = x1_fd - x2_fd

        # Interleave even and odd to form a 1024‑sample segment:
        # seg[0] = even[0], seg[1] = odd[0], seg[2] = even[1], seg[3] = odd[1], ...
        seg = np.empty(2 * len(even), dtype=float)
        seg[0::2] = even
        seg[1::2] = odd

        # Apply 1024‑point window in time domain
        seg_win = seg * w1024

        # FFT + periodogram (scaled by 1/N_fft)
        X = np.fft.fft(seg_win, n=N_fft)
        psd_seg = np.abs(X) ** 2 / N_fft

        psd_accum += psd_seg
        n_segments += 1

    psd = psd_accum / n_segments

    freqs = np.fft.fftfreq(N_fft, d=1.0 / fs)
    half = N_fft // 2
    idx = slice(0, half + 1)
    return freqs[idx], psd[idx]


# Quick check (optional)
if __name__ == "__main__":
    fs = 1000.0
    t = np.arange(0, 1.5, 1 / fs)
    x = np.sin(2 * np.pi * 50 * t) + 0.1 * np.random.randn(len(t))

    f, P = bt_psd_structural_2x512(x, fs=fs, window="hann")
    print("Structural PSD shape:", P.shape)
