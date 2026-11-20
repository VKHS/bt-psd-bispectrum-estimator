"""
fft_cordic_fixed.py

Small CORDIC-based fixed-point FFT + BT-style PSD.

This is the "next step" on top of:
  - a floating-point BT PSD reference
  - fixed-point quantization tools

It gives you a 1D radix-2 FFT that:
  * uses CORDIC-generated twiddles (dyadic-angle rotations as in the CRU) 
  * uses quantized complex add / sub / mul after every butterfly
  * can be dropped into a BT-style PSD estimator

NOTE:
  - This is an *algorithmic* model of the CRU-based FFT, not a cycle-accurate R8MDC.
  - Word-length / frac-bits are configurable via FixedPointConfig.
"""

from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from typing import Tuple


# ---------------------------------------------------------------------
# Fixed-point configuration and helpers
# ---------------------------------------------------------------------

@dataclass
class FixedPointConfig:
    wl: int        # total word length (including sign)
    frac: int      # number of fractional bits
    signed: bool = True

    @property
    def scale(self) -> int:
        return 1 << self.frac

    @property
    def min_int(self) -> int:
        if self.signed:
            return -(1 << (self.wl - 1))
        return 0

    @property
    def max_int(self) -> int:
        if self.signed:
            return (1 << (self.wl - 1)) - 1
        return (1 << self.wl) - 1


def _saturate_int(x: np.ndarray | int, cfg: FixedPointConfig) -> np.ndarray:
    x = np.asarray(x, dtype=np.int64)
    return np.clip(x, cfg.min_int, cfg.max_int).astype(np.int64)


def quantize_real(x: np.ndarray | float, cfg: FixedPointConfig) -> np.ndarray:
    """
    Quantize real float(s) to fixed-point integer representation using
    round-to-nearest and saturation.
    """
    x = np.asarray(x, dtype=float)
    scaled = np.round(x * cfg.scale).astype(np.int64)
    return _saturate_int(scaled, cfg)


def dequantize_real(q: np.ndarray | int, cfg: FixedPointConfig) -> np.ndarray:
    """
    Convert fixed-point integer(s) back to float.
    """
    q = np.asarray(q, dtype=np.int64)
    return q.astype(float) / cfg.scale


def quantize_complex(z: np.ndarray | complex, cfg: FixedPointConfig) -> np.ndarray:
    """
    Quantize complex float(s) to complex fixed-point (two ints).
    """
    z = np.asarray(z, dtype=complex)
    rq = quantize_real(np.real(z), cfg)
    iq = quantize_real(np.imag(z), cfg)
    return rq + 1j * iq


def dequantize_complex(zq: np.ndarray | complex, cfg: FixedPointConfig) -> np.ndarray:
    """
    Dequantize complex fixed-point back to complex float.
    """
    zq = np.asarray(zq)
    r = dequantize_real(np.real(zq), cfg)
    i = dequantize_real(np.imag(zq), cfg)
    return r + 1j * i


# ---------------------------------------------------------------------
# CORDIC rotation (for twiddle generation)
# ---------------------------------------------------------------------

def _cordic_gain(n_iter: int) -> float:
    """Compute CORDIC gain K = Π_i 1/sqrt(1+2^-2i)."""
    return np.prod([1.0 / np.sqrt(1.0 + 2.0 ** (-2 * i)) for i in range(n_iter)])


def cordic_rotate_unit_vector(
    angle_rad: float,
    cfg: FixedPointConfig,
    n_iter: int = 14,
) -> Tuple[float, float]:
    """
    CORDIC rotation of (1,0) by 'angle_rad' using dyadic angles atan(2^-i).

    This is a *fixed-point* algorithm at the integer level, dequantized
    back to float at the end. It is a software stand-in for the CRU
    "twiddle generation" using dyadic tangents as in eqs. (20)–(23). 
    """
    # angle accumulator in float (we don't quantize the angle here)
    z = angle_rad
    atan_lut = [np.arctan(2.0 ** (-i)) for i in range(n_iter)]
    K = _cordic_gain(n_iter)

    # Start from fixed-point representation of (1,0)
    xi = quantize_real(1.0, cfg)[()]
    yi = quantize_real(0.0, cfg)[()]

    for i in range(n_iter):
        if z >= 0:
            # rotate negative (clockwise) or positive depending on sign convention
            x_new = xi - (yi >> i)
            y_new = yi + (xi >> i)
            z -= atan_lut[i]
        else:
            x_new = xi + (yi >> i)
            y_new = yi - (xi >> i)
            z += atan_lut[i]

        xi = _saturate_int(x_new, cfg)
        yi = _saturate_int(y_new, cfg)

    # Dequantize and apply gain compensation
    xr = (xi / cfg.scale) * K
    yr = (yi / cfg.scale) * K
    return float(xr), float(yr)


def cordic_twiddle(angle_rad: float, cfg: FixedPointConfig, n_iter: int = 14) -> complex:
    """
    Approximate e^{-j * angle_rad} using fixed-point CORDIC.

    Note: FFT uses e^{-j 2π k n / N} twiddles; we pass -angle to get
    the correct sign. This returns a complex float that includes the
    finite-precision CORDIC error (like the approximated angles in
    Table 6). 
    """
    xr, yr = cordic_rotate_unit_vector(-angle_rad, cfg, n_iter=n_iter)
    return xr + 1j * yr


# ---------------------------------------------------------------------
# Quantized complex arithmetic
# ---------------------------------------------------------------------

def q_add(a: complex, b: complex, cfg: FixedPointConfig) -> complex:
    """
    Quantized complex addition: (a + b) rounded and saturated.
    """
    a = complex(a)
    b = complex(b)
    r = quantize_real(a.real + b.real, cfg)
    i = quantize_real(a.imag + b.imag, cfg)
    return dequantize_real(r, cfg) + 1j * dequantize_real(i, cfg)


def q_sub(a: complex, b: complex, cfg: FixedPointConfig) -> complex:
    """
    Quantized complex subtraction: (a - b) rounded and saturated.
    """
    a = complex(a)
    b = complex(b)
    r = quantize_real(a.real - b.real, cfg)
    i = quantize_real(a.imag - b.imag, cfg)
    return dequantize_real(r, cfg) + 1j * dequantize_real(i, cfg)


def q_mul(a: complex, b: complex, cfg: FixedPointConfig) -> complex:
    """
    Quantized complex multiplication using integer arithmetic:

      (a_r + j a_i)*(b_r + j b_i)
      = (a_r b_r - a_i b_i) + j (a_r b_i + a_i b_r)

    All intermediate products are scaled and saturated back to cfg.frac.
    """
    a = complex(a)
    b = complex(b)

    ar_q = quantize_real(a.real, cfg)
    ai_q = quantize_real(a.imag, cfg)
    br_q = quantize_real(b.real, cfg)
    bi_q = quantize_real(b.imag, cfg)

    ar = ar_q.astype(np.int64)
    ai = ai_q.astype(np.int64)
    br = br_q.astype(np.int64)
    bi = bi_q.astype(np.int64)

    # Compute real and imag parts in integer domain
    real_int = ar * br - ai * bi
    imag_int = ar * bi + ai * br

    # Rescale back to Q(frac)
    real_q = _saturate_int(real_int >> cfg.frac, cfg)
    imag_q = _saturate_int(imag_int >> cfg.frac, cfg)

    real_f = dequantize_real(real_q, cfg)
    imag_f = dequantize_real(imag_q, cfg)
    return complex(real_f, imag_f)


# ---------------------------------------------------------------------
# Radix-2 FFT using CORDIC twiddles + quantized complex ops
# ---------------------------------------------------------------------

def _bit_reverse_indices(N: int) -> np.ndarray:
    """
    Compute bit-reversal permutation indices for length-N (power of 2).
    """
    m = int(np.log2(N))
    idx = np.arange(N, dtype=int)
    rev = np.zeros(N, dtype=int)
    for i in range(N):
        b = format(i, f"0{m}b")[::-1]
        rev[i] = int(b, 2)
    return rev


def fft_fixed_radix2_cordic(
    x: np.ndarray,
    cfg: FixedPointConfig,
) -> np.ndarray:
    """
    In-place radix-2 Cooley-Tukey FFT using:
      - CORDIC-based twiddles (cordic_twiddle)
      - quantized complex add/sub/mul after each butterfly

    This mimics the core numerical behavior of the CRU-based FFT, but
    in a simple algorithmic form (not the memory-based R8MDC topology). 
    """
    x = np.asarray(x, dtype=complex)
    N = x.shape[0]
    if N & (N - 1) != 0:
        raise ValueError("N must be a power of 2 for radix-2 FFT.")

    # Bit-reversal reordering
    br = _bit_reverse_indices(N)
    X = x[br].copy()

    stages = int(np.log2(N))

    for s in range(1, stages + 1):
        m = 1 << s
        half_m = m >> 1
        angle_step = 2.0 * np.pi / m

        for k in range(0, N, m):
            for j in range(half_m):
                angle = angle_step * j
                tw = cordic_twiddle(angle, cfg)  # approximate e^{-j angle}
                u = X[k + j]
                v = X[k + j + half_m]

                t = q_mul(v, tw, cfg)
                X[k + j] = q_add(u, t, cfg)
                X[k + j + half_m] = q_sub(u, t, cfg)

        # Optional: emulate "safe scaling" (divide by 2 each stage).
        # Uncomment if you want that behavior:
        # X = dequantize_complex(quantize_complex(X / 2.0, cfg), cfg)

    return X


# ---------------------------------------------------------------------
# BT-style PSD using CORDIC-fixed FFT
# ---------------------------------------------------------------------

def _window(window_type: str, L: int) -> np.ndarray:
    wt = window_type.lower()
    if wt in ("rect", "rectangular", "boxcar"):
        return np.ones(L, dtype=float)
    elif wt in ("hann", "hanning"):
        return np.hanning(L).astype(float)
    elif wt == "hamming":
        return np.hamming(L).astype(float)
    elif wt == "blackman":
        return np.blackman(L).astype(float)
    else:
        raise ValueError(f"Unsupported window type: {window_type}")


def _frame_signal(x: np.ndarray, L: int, hop: int) -> np.ndarray:
    """
    Simple overlapping framing (no fancy stride tricks, to keep it clear).
    """
    x = np.asarray(x, dtype=float)
    N = len(x)
    if N < L:
        raise ValueError("Signal shorter than frame length.")

    frames = []
    start = 0
    while start + L <= N:
        frames.append(x[start:start + L])
        start += hop
    return np.stack(frames, axis=0)


def bt_psd_cordic_fixed(
    x: np.ndarray,
    fs: float,
    cfg: FixedPointConfig,
    L: int = 64,
    overlap: float = 0.5,
    window_type: str = "hann",
) -> Tuple[np.ndarray, np.ndarray]:
    """
    BT-style PSD using CORDIC-based fixed-point FFT (length L):

      - segment x into overlapping frames (L, hop=L*(1-overlap))
      - window in time domain
      - run fft_fixed_radix2_cordic on each frame
      - embed safe scaling by dividing FFT output by L (eq. (12)) 
      - form |X|^2, accumulate, then average

    Parameters
    ----------
    x : np.ndarray
        Real-valued input signal.
    fs : float
        Sampling frequency [Hz].
    cfg : FixedPointConfig
        Fixed-point configuration.
    L : int
        FFT / frame length. Must be a power of 2.
    overlap : float
        Fractional overlap (0.5 -> 50%).
    window_type : str
        Time-domain window to apply.

    Returns
    -------
    f : np.ndarray
        One-sided frequency axis.
    Pxx : np.ndarray
        One-sided PSD estimate.
    """
    x = np.asarray(x, dtype=float).flatten()
    if L & (L - 1) != 0:
        raise ValueError("L must be a power of 2 for radix-2 FFT.")
    if not (0 <= overlap < 1):
        raise ValueError("overlap must be in [0, 1).")

    hop = int(L * (1.0 - overlap))
    if hop <= 0:
        raise ValueError("Overlap too large: hop length <= 0.")

    frames = _frame_signal(x, L, hop)
    n_frames = frames.shape[0]

    w = _window(window_type, L)
    U = np.sum(w ** 2)

    P_accum = np.zeros(L, dtype=float)

    for n in range(n_frames):
        seg = frames[n] * w
        # CORDIC-based fixed-point FFT
        X = fft_fixed_radix2_cordic(seg, cfg)
        # Safe scaling: embed 1/L factor into FFT result (as in eq. (12)) 
        X_scaled = X / float(L)
        X_scaled_q = quantize_complex(X_scaled, cfg)
        X_scaled_f = dequantize_complex(X_scaled_q, cfg)

        P_seg = (1.0 / (fs * U)) * np.abs(X_scaled_f) ** 2
        P_accum += P_seg

    P_avg = P_accum / n_frames

    freqs = np.fft.fftfreq(L, d=1.0 / fs)
    half = L // 2
    return freqs[:half], P_avg[:half]


# ---------------------------------------------------------------------
# Simple demo / sanity check
# ---------------------------------------------------------------------

if __name__ == "__main__":
    # Example: compare CORDIC-fixed PSD vs floating BT PSD
    from bt_reference_float import blackman_tukey_psd  # adjust path if needed

    fs = 1000.0
    t = np.arange(0, 1.0, 1.0 / fs)
    # 50 Hz tone + small noise
    x = np.sin(2 * np.pi * 50 * t) + 0.05 * np.random.randn(len(t))

    # Floating reference PSD (BT-style)
    f_ref, P_ref = blackman_tukey_psd(x, fs=fs, L=64, overlap=0.5, window_type="hann", nfft=64)

    # CORDIC-based fixed-point PSD
    cfg = FixedPointConfig(wl=12, frac=10, signed=True)
    f_fx, P_fx = bt_psd_cordic_fixed(x, fs=fs, cfg=cfg, L=64, overlap=0.5, window_type="hann")

    # Interpolate reference onto same freq grid if necessary (here L==nfft so equal)
    # Compute SQNR (dB) between |P_fx| and |P_ref|
    eps = 1e-20
    noise_power = np.mean((P_ref - P_fx) ** 2)
    signal_power = np.mean(P_ref ** 2) + eps
    sqnr_db = 10 * np.log10(signal_power / (noise_power + eps))

    print("CORDIC-fixed BT PSD vs float reference:")
    print("SQNR ≈ %.2f dB" % sqnr_db)
