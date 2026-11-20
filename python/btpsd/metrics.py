"""Accuracy metrics.

``sqnr_db`` is defined in the *linear power domain*::

    SQNR = 10 log10( sum P_ref^2 / sum (P_ref - P_hw)^2 )

The paper does not state its SQNR definition; this one is used throughout the
repository and every number produced here says so.  No fitting or rescaling
of the hardware output is done: the scale between the integer datapath and
the float reference is analytic (see :mod:`btpsd.estimator`).
"""

from __future__ import annotations

import numpy as np

__all__ = ["sqnr_db", "snr_complex_db", "log_spectral_error_db", "rel_rms"]


def sqnr_db(ref, test) -> float:
    ref = np.asarray(ref, dtype=float)
    test = np.asarray(test, dtype=float)
    num = float(np.sum(ref**2))
    den = float(np.sum((ref - test) ** 2))
    if den == 0.0:
        return float("inf")
    return 10.0 * np.log10(num / den)


def snr_complex_db(ref, test) -> float:
    """``10 log10( sum |ref|^2 / sum |ref - test|^2 )`` for complex spectra."""
    ref = np.asarray(ref, dtype=complex)
    test = np.asarray(test, dtype=complex)
    den = float(np.sum(np.abs(ref - test) ** 2))
    if den == 0.0:
        return float("inf")
    return 10.0 * np.log10(float(np.sum(np.abs(ref) ** 2)) / den)


def log_spectral_error_db(ref, test, floor: float = 1e-30) -> np.ndarray:
    """Per-bin ``10 log10(test/ref)`` (dB)."""
    ref = np.maximum(np.asarray(ref, dtype=float), floor)
    test = np.maximum(np.asarray(test, dtype=float), floor)
    return 10.0 * np.log10(test / ref)


def rel_rms(ref, test) -> float:
    ref = np.asarray(ref, dtype=float)
    test = np.asarray(test, dtype=float)
    return float(np.sqrt(np.sum((ref - test) ** 2) / np.sum(ref**2)))
