"""Rotators for the central rotational unit (CRU): constant-coefficient rotators and classical CORDIC.

Constant rotators (the CRU of the paper, Table 5/6, Fig. 17)
-------------------------------------------------------------
A twiddle multiplication ``y = x * exp(-j q pi/512)`` is split into

1. a *free* quadrant operation ``(-j)^quad`` (swap / negate), and
2. a rotation by ``r * pi/512`` with ``0 <= r <= 128`` (at most 45 degrees)
   in either direction (:func:`octant_reduce`),

and the remaining rotation is a multiplication by the integer pair ``(a, b)``
followed by a gain correction::

    y = (x * (a -/+ j b)) * inv / 2**inv_frac,         inv/2**inv_frac ~ 1/|a + j b|

``a``, ``b`` and ``inv`` are constants with few canonical-signed-digit (CSD)
terms, so the datapath is shift-and-add only.  :func:`generate_rom` searches,
for every ``r``, the cheapest ``(a, b, inv)`` whose angle error stays within a
tolerance ``tol`` (the paper uses ``0.01 rad``) and whose gain error stays
below ``gain_tol``.  The printed Table 6 of the paper is kept untouched as
reference data (``data/table6_printed.csv``); this ROM is *generated and
verified* rather than transcribed (see ``docs/PAPER_NOTES.md``).

The classical iterative CORDIC (micro-rotations by ``atan 2**-i``) is provided
as well, with the octant pre-rotation that the original repository script was
missing (plain CORDIC only converges for ``|angle| <= ~99.9 degrees``).
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from .fixedpoint import naf_weight, best_csd_int, rshift, saturate

__all__ = [
    "STEP",
    "PAPER_TOL",
    "RotRom",
    "generate_rom",
    "default_rom",
    "octant_reduce",
    "rotate_fixed",
    "rotate_float",
    "twiddle_q",
    "cordic_gain",
    "cordic_rotate_float",
    "cordic_rotate_int",
]

STEP = np.pi / 512.0  # angle unit of the twiddle index q
PAPER_TOL = 0.01  # rad; angle tolerance quoted for Table 6


# ---------------------------------------------------------------------------
# ROM
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class RotRom:
    """Coefficients for ``r = 0..128`` (units of pi/512)."""

    a: np.ndarray  # int64[129]
    b: np.ndarray  # int64[129]
    inv: np.ndarray  # int64[129]
    inv_frac: int
    tol: float
    gain_tol: float

    @property
    def n(self) -> int:
        return len(self.a)

    def angle_error(self) -> np.ndarray:
        """Implemented minus ideal angle in rad, per ``r``."""
        r = np.arange(self.n)
        return np.arctan2(self.b, self.a) - r * STEP

    def gain_error(self) -> np.ndarray:
        """Relative gain error (gain - 1), per ``r``."""
        rad = np.hypot(self.a, self.b)
        return self.inv * rad / float(1 << self.inv_frac) - 1.0

    def weights(self) -> np.ndarray:
        """Per-``r`` CSD digit counts ``(w(a), w(b), w(inv))``."""
        return np.array([[naf_weight(int(x)) for x in (a, b, i)]
                         for a, b, i in zip(self.a, self.b, self.inv)])

    def adders(self) -> np.ndarray:
        """Adder count of one complex rotation (model: constant multipliers + combine + scale)."""
        w = self.weights()
        wa, wb, wi = w[:, 0], w[:, 1], w[:, 2]
        out = np.where(self.b == 0, 0, 2 * ((wa - 1) + (wb - 1) + 1 + (wi - 1)))
        return out.astype(int)


def generate_rom(
    tol: float = PAPER_TOL,
    gain_tol: float = 2.0**-10,
    inv_frac: int = 14,
    a_max: int = 255,
    k_max: int = 128,
) -> RotRom:
    """Cheapest ``(a, b, inv)`` per angle subject to ``tol`` / ``gain_tol``.

    Cost is the total CSD digit count ``w(a) + w(b) + w(inv)``; ties are broken
    by the smaller angle error, then the smaller ``a``.  For ``r >= 1`` the
    pair must be a genuine rotation (``b >= 1``).  ``r = 0`` is the identity.
    """
    A = [1]
    B = [0]
    INV = [1 << inv_frac]
    for k in range(1, k_max + 1):
        theta = k * STEP
        tn = np.tan(theta)
        best_key = None
        best_val = None
        for a in range(1, a_max + 1):
            wa = naf_weight(a)
            b0 = a * tn
            for b in {int(np.floor(b0)), int(np.ceil(b0))}:
                if b < 1 or b > a:
                    continue
                err = np.arctan2(b, a) - theta
                if abs(err) > tol:
                    continue
                cab = wa + naf_weight(b)
                if best_key is not None and cab + 1 > best_key[0]:
                    continue
                rad = np.hypot(a, b)
                inv = None
                for t in range(1, 7):
                    m = best_csd_int((1 << inv_frac) / rad, t)
                    if m > 0 and abs(m * rad / (1 << inv_frac) - 1.0) <= gain_tol:
                        inv = m
                        break
                if inv is None:
                    continue
                key = (cab + naf_weight(inv), abs(err), a)
                if best_key is None or key < best_key:
                    best_key, best_val = key, (a, b, inv)
        if best_val is None:
            raise ValueError(f"no coefficient pair for k={k}; relax tol or raise a_max")
        a, b, inv = best_val
        A.append(a)
        B.append(b)
        INV.append(inv)
    return RotRom(
        np.array(A, dtype=np.int64),
        np.array(B, dtype=np.int64),
        np.array(INV, dtype=np.int64),
        inv_frac,
        tol,
        gain_tol,
    )


@lru_cache(maxsize=8)
def _cached_rom(tol: float, gain_tol: float, inv_frac: int, a_max: int) -> RotRom:
    return generate_rom(tol, gain_tol, inv_frac, a_max)


def default_rom(tol: float = PAPER_TOL, gain_tol: float = 2.0**-10, inv_frac: int = 14,
                a_max: int = 255) -> RotRom:
    """Cached :func:`generate_rom` (the ROM used by the bit-accurate models)."""
    return _cached_rom(float(tol), float(gain_tol), int(inv_frac), int(a_max))


# ---------------------------------------------------------------------------
# Octant reduction and rotation
# ---------------------------------------------------------------------------
def twiddle_q(index, n: int = 512) -> np.ndarray:
    """Exponent ``e`` of ``W_n^e = exp(-j 2 pi e / n)`` -> ``q`` in units of pi/512 (mod 1024).

    ``2 pi e / n = e * (1024 / n) * pi/512``; ``n`` must divide ``1024 * e``.
    """
    return (np.asarray(index, dtype=np.int64) * 1024 // n) % 1024


def octant_reduce(q):
    """Split ``q`` (units of pi/512, any integer) into ``(quad, r, ccw)``.

    ``exp(-j q pi/512) = (-j)^quad * exp(-/+ j r pi/512)`` with ``0 <= r <= 128``;
    ``ccw`` is True when the residual rotation is counter-clockwise (``+j``).
    """
    q = np.asarray(q, dtype=np.int64) % 1024
    quad = q // 256
    r = q % 256
    flip = r > 128
    quad = (quad + flip.astype(np.int64)) % 4
    r = np.where(flip, 256 - r, r)
    return quad, r, flip


def _quadrant(zr, zi, quad):
    """Multiply by ``(-j)^quad`` (swap / negate)."""
    outr = np.select([quad == 0, quad == 1, quad == 2, quad == 3], [zr, zi, -zr, -zi])
    outi = np.select([quad == 0, quad == 1, quad == 2, quad == 3], [zi, -zr, -zi, zr])
    return outr, outi


def rotate_fixed(xr, xi, q, rom: RotRom, wl: int, rounding: str = "round"):
    """Bit-accurate ``y = x * exp(-j q pi/512)`` through the CRU datapath.

    ``xr, xi`` are integer arrays (``wl``-bit words); ``q`` broadcasts against
    them.  Products are exact, one rounding after the gain correction, then
    the free quadrant operation and saturation to ``wl`` bits.
    """
    if wl > 28:
        raise ValueError("int64 model supports wl <= 28")
    xr, xi, q = np.broadcast_arrays(np.asarray(xr, dtype=np.int64), np.asarray(xi, dtype=np.int64),
                                    np.asarray(q, dtype=np.int64))
    quad, r, ccw = octant_reduce(q)
    a = rom.a[r]
    b = rom.b[r]
    inv = rom.inv[r]
    s = np.where(ccw, -1, 1)
    yr = xr * a + s * xi * b
    yi = xi * a - s * xr * b
    zr = rshift(yr * inv, rom.inv_frac, rounding)
    zi = rshift(yi * inv, rom.inv_frac, rounding)
    zr, zi = _quadrant(zr, zi, quad)
    return saturate(zr, wl), saturate(zi, wl)


def rotate_float(x, q, rom: RotRom | None = None):
    """Float ``x * exp(-j q pi/512)``; with ``rom`` uses the ROM's (a, b, inv) exactly
    (isolating the coefficient error from rounding noise)."""
    x, q = np.broadcast_arrays(np.asarray(x, dtype=complex), np.asarray(q, dtype=np.int64))
    if rom is None:
        return x * np.exp(-1j * q * STEP)
    quad, r, ccw = octant_reduce(q)
    a = rom.a[r].astype(float)
    b = rom.b[r].astype(float)
    gain = rom.inv[r] / float(1 << rom.inv_frac)
    rot = np.where(ccw, a + 1j * b, a - 1j * b) * gain
    return x * rot * (-1j) ** quad


# ---------------------------------------------------------------------------
# Classical CORDIC (rotation mode)
# ---------------------------------------------------------------------------
def cordic_gain(n_iter: int) -> float:
    """Aggregate gain ``prod sqrt(1 + 2**-2i)`` of ``n_iter`` micro-rotations."""
    i = np.arange(n_iter)
    return float(np.prod(np.sqrt(1.0 + 4.0 ** (-i))))


def cordic_rotate_float(x: complex, angle: float, n_iter: int = 16) -> complex:
    """Rotate ``x`` by ``angle`` rad (any value) with ``n_iter`` micro-rotations (gain removed).

    The angle is first reduced to ``[-pi/4, pi/4]`` by an exact quadrant
    rotation; plain CORDIC alone only converges for ``|angle| <~ 99.9 deg``.
    """
    quad = int(np.round(angle / (np.pi / 2)))
    ang = angle - quad * np.pi / 2
    z = complex(x) * (1j) ** quad
    xr, xi = z.real, z.imag
    for i in range(n_iter):
        d = 1.0 if ang >= 0 else -1.0
        xr, xi = xr - d * xi * 2.0**-i, xi + d * xr * 2.0**-i
        ang -= d * np.arctan(2.0**-i)
    g = cordic_gain(n_iter)
    return complex(xr / g, xi / g)


def cordic_rotate_int(xr, xi, q, wl: int, n_iter: int = 14, frac: int = 14, rounding: str = "round"):
    """Integer CORDIC twiddle ``x * exp(-j q pi/512)`` (``q`` in units of pi/512).

    Quadrant handled exactly, residual angle (|.| <= pi/4) by ``n_iter``
    micro-rotations with ``frac``-bit angle constants; the aggregate gain is
    removed by a final multiplication with a rounded constant.
    """
    xr, xi, q = np.broadcast_arrays(np.asarray(xr, dtype=np.int64), np.asarray(xi, dtype=np.int64),
                                    np.asarray(q, dtype=np.int64) % 1024)
    xr, xi = xr.copy(), xi.copy()
    quad = (q + 128) // 256  # nearest quadrant, residual in [-128, 127] steps
    res = q - quad * 256
    # exp(-j q pi/512) = (-j)^quad * exp(-j res pi/512): rotate by -res steps
    zr, zi = _quadrant(xr, xi, quad % 4)
    ang = -res * STEP
    atan_tab = np.round(np.arctan(2.0 ** -np.arange(n_iter)) * (1 << frac)).astype(np.int64)
    zang = np.round(ang * (1 << frac)).astype(np.int64)
    for i in range(n_iter):
        d = np.where(zang >= 0, 1, -1)
        zr, zi = zr - d * (zi >> i), zi + d * (zr >> i)
        zang = zang - d * atan_tab[i]
    k = int(round((1 << frac) / cordic_gain(n_iter)))
    zr = rshift(zr * k, frac, rounding)
    zi = rshift(zi * k, frac, rounding)
    return saturate(zr, wl), saturate(zi, wl)
