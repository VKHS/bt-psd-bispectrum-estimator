"""Fixed-point arithmetic primitives shared by every bit-accurate model.

Everything here is *integer* arithmetic on ``numpy.int64`` arrays (or Python
ints).  The same conventions are used by the RTL in ``rtl/`` so that the
Python models can serve as golden references for the testbenches:

* two's-complement words of ``wl`` bits, ``frac`` fractional bits
  (``Q1.(wl-1)`` for normalised data in ``[-1, 1)``),
* right shifts are *arithmetic*; ``rounding='round'`` adds half an LSB before
  the shift (round-half-up), ``rounding='trunc'`` simply floors,
* overflow either saturates (default, what the RTL does) or wraps.

The module also provides canonical-signed-digit (CSD / non-adjacent form)
helpers.  Multiplication by a constant realised as a handful of shifted
additions is the arithmetic primitive of the paper's multiplierless central
rotational unit (CRU), of the folded fractional-delay filter and of the folded
window filter.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "QFormat",
    "FixedPointConfig",
    "saturate",
    "wrap",
    "rshift",
    "quantize_real",
    "dequantize_real",
    "naf",
    "naf_weight",
    "naf_weight_array",
    "csd_value",
    "best_csd_int",
    "csd_from_real",
]

ROUNDING_MODES = ("trunc", "round")


# ---------------------------------------------------------------------------
# Word formats
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class QFormat:
    """Signed (or unsigned) fixed-point word: ``wl`` bits, ``frac`` fraction bits."""

    wl: int
    frac: int
    signed: bool = True

    def __post_init__(self) -> None:
        if self.wl < 2:
            raise ValueError("wl must be >= 2")
        if not 0 <= self.frac < self.wl + 32:
            raise ValueError("frac out of range")

    @property
    def scale(self) -> int:
        return 1 << self.frac

    @property
    def min_int(self) -> int:
        return -(1 << (self.wl - 1)) if self.signed else 0

    @property
    def max_int(self) -> int:
        return (1 << (self.wl - 1)) - 1 if self.signed else (1 << self.wl) - 1

    def quantize(self, x, rounding: str = "round", sat: bool = True) -> np.ndarray:
        """Real values -> integer codes."""
        return quantize_real(x, self, rounding=rounding, sat=sat)

    def to_float(self, q) -> np.ndarray:
        """Integer codes -> real values."""
        return dequantize_real(q, self)


# Name used by the original repository scripts.
FixedPointConfig = QFormat


# ---------------------------------------------------------------------------
# Elementary integer operations
# ---------------------------------------------------------------------------
def saturate(x, wl: int, signed: bool = True) -> np.ndarray:
    """Clip integers to the representable range of a ``wl``-bit word."""
    x = np.asarray(x, dtype=np.int64)
    if signed:
        lo, hi = -(1 << (wl - 1)), (1 << (wl - 1)) - 1
    else:
        lo, hi = 0, (1 << wl) - 1
    return np.clip(x, lo, hi)


def wrap(x, wl: int) -> np.ndarray:
    """Two's-complement wrap-around to ``wl`` bits."""
    x = np.asarray(x, dtype=np.int64)
    m = (1 << wl) - 1
    x = x & m
    return np.where(x >= (1 << (wl - 1)), x - (1 << wl), x)


def rshift(x, s: int, rounding: str = "round") -> np.ndarray:
    """Arithmetic right shift by ``s >= 0`` with the chosen rounding.

    ``'trunc'`` is ``x >> s`` (floor).  ``'round'`` is round-half-up,
    ``(x + 2**(s-1)) >> s`` - one adder in hardware.
    """
    if rounding not in ROUNDING_MODES:
        raise ValueError(f"rounding must be one of {ROUNDING_MODES}")
    x = np.asarray(x, dtype=np.int64)
    if s < 0:
        raise ValueError("shift must be >= 0")
    if s == 0:
        return x
    if rounding == "trunc":
        return x >> s
    return (x + (1 << (s - 1))) >> s


def quantize_real(x, cfg: QFormat, rounding: str = "round", sat: bool = True) -> np.ndarray:
    """Quantise real numbers to integer codes of ``cfg``."""
    x = np.asarray(x, dtype=float)
    scaled = x * cfg.scale
    q = np.floor(scaled + 0.5) if rounding == "round" else np.floor(scaled)
    q = q.astype(np.int64)
    return saturate(q, cfg.wl, cfg.signed) if sat else wrap(q, cfg.wl)


def dequantize_real(q, cfg: QFormat) -> np.ndarray:
    """Integer codes -> real numbers."""
    return np.asarray(q, dtype=np.int64).astype(float) / cfg.scale


# ---------------------------------------------------------------------------
# Canonical signed digit (non-adjacent form)
# ---------------------------------------------------------------------------
def naf(n: int) -> list[tuple[int, int]]:
    """Non-adjacent form of the integer ``n``.

    Returns ``[(sign, shift), ...]`` (ascending shifts) with
    ``n == sum(sign * 2**shift)`` and the minimum possible number of non-zero
    digits.  ``naf(0) == []``.
    """
    n = int(n)
    neg = n < 0
    n = -n if neg else n
    digits: list[tuple[int, int]] = []
    shift = 0
    while n:
        if n & 1:
            d = 2 - (n & 3)  # +1 if n % 4 == 1, -1 if n % 4 == 3
            digits.append((-d if neg else d, shift))
            n -= d
        n >>= 1
        shift += 1
    return digits


def _popcount_u64(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x).astype(np.uint64)
    with np.errstate(over="ignore"):  # SWAR popcount relies on modular uint64 arithmetic
        x = x - ((x >> np.uint64(1)) & np.uint64(0x5555555555555555))
        x = (x & np.uint64(0x3333333333333333)) + ((x >> np.uint64(2)) & np.uint64(0x3333333333333333))
        x = (x + (x >> np.uint64(4))) & np.uint64(0x0F0F0F0F0F0F0F0F)
        return ((x * np.uint64(0x0101010101010101)) >> np.uint64(56)).astype(np.int64)


def naf_weight_array(n) -> np.ndarray:
    """Vectorised NAF weight: ``popcount((|n| ^ 3|n|) >> 1)`` (Reitwiesner)."""
    a = np.abs(np.asarray(n, dtype=np.int64)).astype(np.uint64)
    with np.errstate(over="ignore"):
        return _popcount_u64((a ^ (a * np.uint64(3))) >> np.uint64(1))


def naf_weight(n: int) -> int:
    """Number of non-zero CSD digits of ``n`` (= adders needed + 1)."""
    return int(naf_weight_array(int(n)))


def csd_value(digits: list[tuple[int, int]]) -> int:
    """Inverse of :func:`naf`."""
    return sum(s * (1 << e) for s, e in digits)


def best_csd_int(target: float, max_terms: int) -> int:
    """Integer closest to ``target`` whose CSD weight is ``<= max_terms``.

    Exact (the search window doubles until no farther integer can win).
    Ties go to the smaller magnitude.  ``max_terms >= 1``; a target that
    rounds to 0 returns 0.
    """
    if max_terms < 1:
        raise ValueError("max_terms must be >= 1")
    sign = -1 if target < 0 else 1
    t = abs(float(target))
    base = int(np.floor(t))
    w = 16
    while True:
        lo, hi = max(0, base - w), base + 1 + w
        cand = np.arange(lo, hi + 1, dtype=np.int64)
        ok = naf_weight_array(cand) <= max_terms
        if ok.any():
            c = cand[ok]
            dist = np.abs(c - t)
            j = int(np.argmin(dist))  # ascending candidates -> ties keep the smaller one
            if dist[j] <= w:  # every integer outside the window is farther than w
                return sign * int(c[j])
        w *= 2


def csd_from_real(value: float, frac: int, max_terms: int) -> tuple[int, list[tuple[int, int]]]:
    """Approximate ``value`` by ``m / 2**frac`` with at most ``max_terms`` CSD digits.

    Returns ``(m, digits)`` where ``m`` is the integer numerator and ``digits``
    its non-adjacent form.  The implied constant is ``m / 2**frac``.
    """
    m = best_csd_int(value * (1 << frac), max_terms)
    return m, naf(m)
