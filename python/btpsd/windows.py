"""Data windows and their *frequency-domain* kernels (paper Table 4, Eq. (13)).

The estimator of the paper applies the window **after** the FFT, as a short
convolution on the spectrum instead of a multiplication in time::

    F[x(n) w(n)] = 1/(2*pi) * X(k) * W(k)                    (Eq. 13)

For the raised-cosine family

    w[n] = a0 - a1*cos(2*pi*n/L) + a2*cos(4*pi*n/L)     (periodic form)

the DFT of the window has only 1, 3 or 5 non-zero bins, so the convolution is
an exact 1/3/5-tap filter - the numbers quoted as "significant taps" in
Table 4 (Rectangular 1, Hann 3, Hamming 3, Blackman 5).  This only holds for
the *periodic* (DFT-even) window, which is therefore the default here.  The
symmetric windows returned by ``numpy.hanning`` & co. are not exactly sparse
in frequency.

For windows that are not raised cosines (Triangular, Parzen, Kaiser, ...)
:func:`approx_spectral_kernel` truncates the DFT to its central taps and
reports how much window energy the truncation keeps.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "RAISED_COSINE",
    "window",
    "spectral_kernel",
    "approx_spectral_kernel",
    "apply_window_freq",
    "TABLE4_SIGNIFICANT_TAPS",
]

# name -> (a0, a1, a2) of  w[n] = a0 - a1 cos(2 pi n / L) + a2 cos(4 pi n / L)
RAISED_COSINE: dict[str, tuple[float, float, float]] = {
    "rect": (1.0, 0.0, 0.0),
    "hann": (0.5, 0.5, 0.0),
    "hamming": (0.54, 0.46, 0.0),
    "blackman": (0.42, 0.5, 0.08),
}

_ALIASES = {
    "rectangular": "rect",
    "boxcar": "rect",
    "hanning": "hann",
}

# "F. response significant taps" as printed in Table 4 of the paper.
TABLE4_SIGNIFICANT_TAPS = {
    "rect": 1,
    "triangular": 5,
    "parzen": 3,
    "hann": 3,
    "hamming": 3,
    "blackman": 5,
    "kaiser": 5,
}


def _canon(name: str) -> str:
    n = name.lower()
    return _ALIASES.get(n, n)


def window(name: str, L: int, periodic: bool = True) -> np.ndarray:
    """Window of length ``L``.

    ``periodic=True`` (default) gives the DFT-even form whose spectrum is
    exactly 1/3/5 bins wide for the raised-cosine windows.
    ``periodic=False`` gives the symmetric form (``numpy.hanning`` etc.).
    Also supports ``'triangular'``/``'bartlett'``, ``'parzen'`` and
    ``'kaiser'`` (beta=8.6) via :mod:`scipy.signal`.
    """
    n = _canon(name)
    if n in RAISED_COSINE:
        a0, a1, a2 = RAISED_COSINE[n]
        N = L if periodic else L - 1
        k = np.arange(L)
        return a0 - a1 * np.cos(2 * np.pi * k / N) + a2 * np.cos(4 * np.pi * k / N)
    from scipy import signal as sg

    if n in ("triangular", "bartlett"):
        return sg.get_window("bartlett", L, fftbins=periodic)
    if n == "parzen":
        return sg.get_window("parzen", L, fftbins=periodic)
    if n == "kaiser":
        return sg.get_window(("kaiser", 8.6), L, fftbins=periodic)
    raise ValueError(f"unsupported window: {name!r}")


def spectral_kernel(name: str) -> dict[int, float]:
    """Exact DFT of a periodic raised-cosine window as ``{offset: coefficient}``.

    With ``w[n] = sum_d W[d] exp(j 2 pi d n / L)`` the windowed spectrum is
    ``Xw[k] = sum_d W[d] X[k - d]``.
    """
    n = _canon(name)
    if n not in RAISED_COSINE:
        raise ValueError(f"{name!r} is not a raised-cosine window; use approx_spectral_kernel")
    a0, a1, a2 = RAISED_COSINE[n]
    ker = {0: a0}
    if a1:
        ker[-1] = ker[1] = -a1 / 2.0
    if a2:
        ker[-2] = ker[2] = a2 / 2.0
    return dict(sorted(ker.items()))


def approx_spectral_kernel(name: str, ntaps: int, L: int = 1024) -> tuple[dict[int, float], float]:
    """Central ``ntaps`` (odd) DFT coefficients of an arbitrary window.

    Returns ``(kernel, kept_energy_fraction)`` where the second item is the
    fraction of ``sum |W[d]|^2`` retained by the truncation.
    """
    if ntaps % 2 != 1:
        raise ValueError("ntaps must be odd")
    w = window(name, L, periodic=True)
    W = np.fft.fft(w) / L  # w[n] = sum_d W[d] e^{+j2 pi d n/L}
    half = ntaps // 2
    ker = {d: float(np.real(W[d % L])) for d in range(-half, half + 1)}
    # window is real & even about n=0 (periodic form) -> W is real
    kept = sum(v * v for v in ker.values()) / float(np.sum(np.abs(W) ** 2))
    return ker, float(kept)


def apply_window_freq(X, kernel: dict[int, float], axis: int = -1) -> np.ndarray:
    """Circular convolution ``Xw[k] = sum_d W[d] X[k-d]`` along ``axis``."""
    X = np.asarray(X)
    out = np.zeros_like(X, dtype=np.result_type(X.dtype, np.float64))
    for d, c in kernel.items():
        out = out + c * np.roll(X, d, axis=axis)
    return out
