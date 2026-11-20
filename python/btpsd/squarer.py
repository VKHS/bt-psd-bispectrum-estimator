"""Sum-of-odd-numbers squarer and magnitude-squared unit (paper Sec. on |X|^2, Fig. 22-23 style).

``n**2 = 1 + 3 + 5 + ... + (2|n| - 1)``: a squarer built from one adder and a
counter needs ``|n|`` cycles and no multiplier.  The value is exact; the
*cycle count* is data dependent, which matters for throughput (see
:func:`spectrum_power_cycles`).
"""

from __future__ import annotations

import numpy as np

__all__ = ["square_odds", "square_serial", "mag_sq", "spectrum_power_cycles"]


def square_odds(n: int) -> int:
    """Literal serial algorithm (for tests): accumulate the first ``|n|`` odd numbers."""
    acc, odd = 0, 1
    for _ in range(abs(int(n))):
        acc += odd
        odd += 2
    return acc


def square_serial(n):
    """``(n**2, cycles)`` of the sum-of-odds squarer (vectorised; cycles = ``|n|``)."""
    n = np.asarray(n, dtype=np.int64)
    return n * n, np.abs(n)


def mag_sq(re, im, units: int = 2):
    """``(re**2 + im**2, cycles)``.

    ``units = 2`` squarers in parallel -> ``max(|re|, |im|)`` cycles,
    ``units = 1`` shared squarer -> ``|re| + |im|`` cycles.
    """
    sr, cr = square_serial(re)
    si, ci = square_serial(im)
    cyc = np.maximum(cr, ci) if units == 2 else cr + ci
    return sr + si, cyc


def spectrum_power_cycles(re, im, units: int = 2) -> dict:
    """Cycles needed to square all bins of one segment spectrum (sequential bins)."""
    _, cyc = mag_sq(re, im, units)
    return {"total": int(cyc.sum()), "max_per_bin": int(cyc.max()), "mean_per_bin": float(cyc.mean())}
