# fixed_point_model.py
#
# Fixed‑point utilities and CORDIC‑style rotation tools that reflect
# the dyadic‑angle, shift‑and‑add philosophy used in the paper’s CRU
# (Table 5 & Table 6, Fig. 17). :contentReference[oaicite:14]{index=14}
#
# This code is intended for algorithm‑level and SQNR studies, not as
# an exact micro‑architecture replica of the FPGA design.

from __future__ import annotations
import numpy as np
from dataclasses import dataclass


# ---------------------------
# Fixed‑point configuration & quantization
# ---------------------------

@dataclass
class FixedPointConfig:
    wl: int         # total word length in bits (including sign)
    frac: int       # fractional bits
    signed: bool = True

    @property
    def scale(self) -> int:
        return 1 << self.frac

    @property
    def min_int(self) -> int:
        if self.signed:
            return -(1 << (self.wl - 1))
        else:
            return 0

    @property
    def max_int(self) -> int:
        if self.signed:
            return (1 << (self.wl - 1)) - 1
        else:
            return (1 << self.wl) - 1


def quantize_real(x: np.ndarray | float, cfg: FixedPointConfig) -> np.ndarray:
    """
    Quantize real values to fixed‑point integers using round‑to‑nearest
    and saturation.
    """
    x = np.asarray(x, dtype=float)
    val = np.round(x * cfg.scale).astype(np.int64)
    val = np.clip(val, cfg.min_int, cfg.max_int)
    return val


def dequantize_real(q: np.ndarray | int, cfg: FixedPointConfig) -> np.ndarray:
    """
    Convert fixed‑point integers back to float using cfg.frac.
    """
    q = np.asarray(q, dtype=np.int64)
    return q.astype(float) / cfg.scale


# ---------------------------
# CORDIC‑style rotation (vector rotation mode)
# ---------------------------

def _cordic_k_gain(n_iter: int) -> float:
    """Compute the CORDIC gain K for n_iter steps."""
    return np.prod([1.0 / np.sqrt(1.0 + 2.0 ** (-2 * i)) for i in range(n_iter)])


def cordic_rotate_float(
    x: float,
    y: float,
    angle_rad: float,
    n_iter: int = 20,
) -> tuple[float, float]:
    """
    Floating‑point CORDIC rotation (rotation mode) for reference.

    Starting from (x,y), apply a rotation of 'angle_rad' using
    dyadic steps atan(2^-i). Returns (x', y') in float.
    """
    z = angle_rad
    atan_lut = [np.arctan(2.0 ** (-i)) for i in range(n_iter)]
    K = _cordic_k_gain(n_iter)

    for i in range(n_iter):
        if z >= 0:
            x_new = x - y * (2.0 ** (-i))
            y_new = y + x * (2.0 ** (-i))
            z -= atan_lut[i]
        else:
            x_new = x + y * (2.0 ** (-i))
            y_new = y - x * (2.0 ** (-i))
            z += atan_lut[i]
        x, y = x_new, y_new

    # Compensate for CORDIC gain
    return x * K, y * K


def cordic_rotate_fixed_int(
    angle_rad: float,
    n_iter: int,
    cfg: FixedPointConfig,
) -> tuple[float, float]:
    """
    Fixed‑point CORDIC rotation of the unit vector (1,0) by angle_rad.

    - Uses integer shifts and adds on a fixed‑point representation.
    - Keeps the angle accumulator z in float for simplicity.
    - Returns (cos(angle), sin(angle)) *approximately* as floats.

    This models the idea of Table 5 & 6 (dyadic tangents, small shift
    depths) without encoding every micro‑rotation case separately. :contentReference[oaicite:15]{index=15}
    """
    z = angle_rad
    atan_lut = [np.arctan(2.0 ** (-i)) for i in range(n_iter)]
    K = _cordic_k_gain(n_iter)

    # Start with (1,0) in fixed‑point
    xi = int(round(1.0 * cfg.scale))
    yi = 0

    for i in range(n_iter):
        if z >= 0:
            xi_new = xi - (yi >> i)
            yi_new = yi + (xi >> i)
            z -= atan_lut[i]
        else:
            xi_new = xi + (yi >> i)
            yi_new = yi - (xi >> i)
            z += atan_lut[i]

        # Saturation
        xi = int(np.clip(xi_new, cfg.min_int, cfg.max_int))
        yi = int(np.clip(yi_new, cfg.min_int, cfg.max_int))

    # Dequantize and apply gain compensation
    x_float = (xi / cfg.scale) * K
    y_float = (yi / cfg.scale) * K
    return x_float, y_float


def cordic_twiddle(
    angle_rad: float,
    cfg: FixedPointConfig,
    n_iter: int = 14,
) -> complex:
    """
    Use integer CORDIC to approximate e^{j * angle_rad} in fixed‑point.

    This plays the role of a "CORDIC‑generated twiddle factor" as in the
    central rotational unit (CRU), though here it is just a software model. :contentReference[oaicite:16]{index=16}
    """
    xr, yr = cordic_rotate_fixed_int(angle_rad, n_iter=n_iter, cfg=cfg)
    return xr + 1j * yr


# ---------------------------
# Quantized complex multiplication
# ---------------------------

def quantized_complex_mult(
    a: complex,
    b: complex,
    cfg: FixedPointConfig,
) -> complex:
    """
    Multiply two complex numbers in fixed‑point:

      (a_r + j a_i) * (b_r + j b_i)

    using integer arithmetic with saturation and scaling by cfg.frac.
    This lets you study the effect of finite precision on complex
    multiplications (like FFT twiddle rotations, FD/window taps, etc.).
    """
    a = complex(a)
    b = complex(b)

    # Quantize real and imaginary parts separately
    ar_q = quantize_real(a.real, cfg)
    ai_q = quantize_real(a.imag, cfg)
    br_q = quantize_real(b.real, cfg)
    bi_q = quantize_real(b.imag, cfg)

    ar = ar_q.astype(np.int64)
    ai = ai_q.astype(np.int64)
    br = br_q.astype(np.int64)
    bi = bi_q.astype(np.int64)

    # Integer complex multiply:
    # (ar + j ai)*(br + j bi) = (ar*br - ai*bi) + j (ar*bi + ai*br)
    real_int = ar * br - ai * bi
    imag_int = ar * bi + ai * br

    # Scale back from scale^2 to scale^1
    real_q = np.clip(real_int // cfg.scale, cfg.min_int, cfg.max_int)
    imag_q = np.clip(imag_int // cfg.scale, cfg.min_int, cfg.max_int)

    real = dequantize_real(real_q, cfg)
    imag = dequantize_real(imag_q, cfg)
    return real + 1j * imag


# Example quick checks
if __name__ == "__main__":
    cfg = FixedPointConfig(wl=16, frac=14, signed=True)

    # Test CORDIC rotation vs float
    for ang in [0.0, np.pi / 4, np.pi / 2]:
        xf, yf = cordic_rotate_float(1.0, 0.0, ang, n_iter=20)
        xi, yi = cordic_rotate_fixed_int(ang, n_iter=14, cfg=cfg)
        print(f"angle={ang:.3f} float=({xf:.6f},{yf:.6f})  fixed=({xi:.6f},{yi:.6f})")

    # Test twiddle + quantized complex multiplication
    angle = 2 * np.pi * 5 / 64
    tw = cordic_twiddle(angle, cfg)
    z = 1.0 + 0.5j
    zq = quantized_complex_mult(z, tw, cfg)
    print("Twiddle approx:", tw)
    print("Quantized product:", zq)
