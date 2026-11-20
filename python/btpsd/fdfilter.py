"""Fractional-delay (half-bin) merge of two 512-point chunk spectra (paper Eqs. 3-6, Figs. 4-5).

Setting
-------
Two consecutive 512-sample chunks ``a`` and ``b`` form one 1024-sample segment
``[a b]``.  With ``A = FFT512(a)`` and ``B = FFT512(b)`` the 1024-point DFT of
the segment is, *exactly*::

    X(2m)   = A(m) + B(m)                                       (Eq. 3)
    X(2m+1) = A(m + 1/2) - B(m + 1/2)                           (Eq. 4)

where ``A(m + 1/2) = sum_n a[n] exp(-j 2 pi (m + 1/2) n / 512)`` is the DTFT of
the chunk sampled half-way between two FFT bins.  Eq. (4) therefore needs a
*half-bin shift of the spectrum*, i.e. a half-sample fractional delay along the
frequency axis.  The paper replaces the 512 radix-2 butterflies of a classic
last FFT stage by a short fractional-delay (FD) FIR filter applied to the
difference ``D = A - B`` (Fig. 5: the half-sample delay ``Z^{1/2}`` sits behind
the subtractor).

Derivation of the interpolation kernel (verified numerically in the tests)
-------------------------------------------------------------------------
The DTFT of a length-``N`` sequence is interpolated exactly by the Dirichlet
kernel::

    A(m + 1/2) = sum_d K(d) A(m + d),
    K(d) = j exp(-j pi (d - 1/2)/N) / (N sin(pi (d - 1/2)/N))
         ~ j / (pi (d - 1/2))                      for |d| << N.

Inserting the sign alternation ``A'(k) = (-1)^k A(k)`` turns the (anti-
symmetric) kernel into a real, *symmetric* half-sample interpolator::

    A(m + 1/2) = -j (-1)^m  sum_d h[d] A'(m + d),   h[d] = sinc(d - 1/2)

For taps ``d = -2..3`` this is ``{0.1273, -0.2122, 0.6366, 0.6366, -0.2122,
0.1273}`` - exactly the 6-tap filter of Eq. (6) (the paper rounds the outer
taps to ``1/8`` so that they are a single shift).  The factors ``(-1)^k`` and
``-j`` are sign changes / a re-im swap, i.e. free in hardware.  The filter
needs the three bins ``m+1 .. m+3`` ahead of the output bin, which is the
"D = 3" delay of Fig. 5.

Odd-length and one-sided ("unidirectional") kernels are available for the
group-delay / SNR study of Table 3 and Figs. 6-11.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "PAPER_FD6",
    "PAPER_FD6_START",
    "fd_kernel",
    "exact_half_bin",
    "interp_half_bin",
    "merge_chunk_spectra",
    "fd_frequency_response",
    "fd_group_delay",
]

#: Eq. (6) of the paper, taps at offsets d = -2 .. +3 around the output bin.
PAPER_FD6 = np.array([0.125, -0.2122, 0.6366, 0.6366, -0.2122, 0.125])
PAPER_FD6_START = -2


def fd_kernel(L: int, kind: str = "sinc") -> tuple[np.ndarray, int]:
    """Half-sample interpolation kernel ``(h, start)``.

    ``h[i]`` multiplies ``A'(m + start + i)``.

    kind
        ``'paper6'``  Eq. (6) (``L`` must be 6)
        ``'sinc'``    truncated ideal interpolator, *symmetric* for even ``L``
                      (constant group delay); for odd ``L`` it is truncated
                      asymmetrically (one extra tap on the left)
        ``'unidirectional'``  all taps on the already-available side (causal
                      only, ``d = -(L-1)..0``)
        ``'lagrange'``  maximally-flat (Lagrange) half-sample interpolator
    """
    if kind == "paper6":
        if L != 6:
            raise ValueError("paper6 kernel has exactly 6 taps")
        return PAPER_FD6.copy(), PAPER_FD6_START
    if L < 1:
        raise ValueError("L must be >= 1")
    if kind in ("sinc", "lagrange"):
        start = -(L // 2 - 1) if L % 2 == 0 else -((L - 1) // 2)
    elif kind == "unidirectional":
        start = -(L - 1)
    else:
        raise ValueError(f"unknown kernel kind {kind!r}")
    d = start + np.arange(L)
    if kind in ("sinc", "unidirectional"):
        h = np.sinc(d - 0.5)
    else:  # Lagrange basis polynomials, interpolation point at d = 1/2
        tau = 0.5 - start
        p = np.arange(L, dtype=float)
        h = np.empty(L)
        for i in range(L):
            num = np.prod([tau - p[j] for j in range(L) if j != i])
            den = np.prod([p[i] - p[j] for j in range(L) if j != i])
            h[i] = num / den
    return h, start


def exact_half_bin(A: np.ndarray) -> np.ndarray:
    """Exact ``A(m + 1/2)`` from the 512-point DFT ``A`` (last axis).

    Uses the inverse DFT to recover the chunk and re-evaluates it with a
    half-bin phase ramp, i.e. no truncation of the Dirichlet kernel.
    """
    A = np.asarray(A, dtype=complex)
    N = A.shape[-1]
    a = np.fft.ifft(A, axis=-1)
    ramp = np.exp(-1j * np.pi * np.arange(N) / N)
    return np.fft.fft(a * ramp, axis=-1)


def interp_half_bin(D: np.ndarray, h: np.ndarray, start: int) -> np.ndarray:
    """FD-filter estimate of ``D(m + 1/2)`` from the DFT bins ``D(m)`` (last axis).

    Implements ``-j (-1)^m  sum_i h[i] (-1)^(m+d_i) D(m + d_i)`` with
    ``d_i = start + i`` and circular bin indexing.
    """
    D = np.asarray(D, dtype=complex)
    N = D.shape[-1]
    s = (-1.0) ** np.arange(N)
    Dp = D * s
    acc = np.zeros_like(Dp)
    for i, hi in enumerate(np.asarray(h, dtype=float)):
        acc = acc + hi * np.roll(Dp, -(start + i), axis=-1)
    return -1j * s * acc


def merge_chunk_spectra(
    A: np.ndarray,
    B: np.ndarray,
    h: np.ndarray | None = None,
    start: int | None = None,
    mode: str = "fd",
) -> np.ndarray:
    """1024-bin segment spectrum from two 512-point chunk spectra.

    ``mode='exact'`` evaluates Eq. (4) exactly (reference; infeasible as a
    cheap datapath), ``mode='fd'`` uses the FD filter ``(h, start)``
    (default: the paper's 6-tap filter).  Output bins are interleaved
    ``X[0::2] = even``, ``X[1::2] = odd``; arrays may carry leading batch axes.
    """
    A = np.asarray(A, dtype=complex)
    B = np.asarray(B, dtype=complex)
    if A.shape != B.shape:
        raise ValueError("A and B must have the same shape")
    N = A.shape[-1]
    X = np.empty(A.shape[:-1] + (2 * N,), dtype=complex)
    X[..., 0::2] = A + B
    if mode == "exact":
        X[..., 1::2] = exact_half_bin(A - B)
    elif mode == "fd":
        if h is None:
            h, start = PAPER_FD6, PAPER_FD6_START
        if start is None:
            raise ValueError("start must be given together with h")
        X[..., 1::2] = interp_half_bin(A - B, h, start)
    else:
        raise ValueError("mode must be 'fd' or 'exact'")
    return X


def fd_frequency_response(h: np.ndarray, n: int = 1024) -> tuple[np.ndarray, np.ndarray]:
    """``(omega, H)`` of the FIR ``h`` over ``[0, pi)`` (rad/sample)."""
    w = np.linspace(0.0, np.pi, n, endpoint=False)
    H = np.exp(-1j * np.outer(w, np.arange(len(h)))) @ np.asarray(h, dtype=float)
    return w, H


def fd_group_delay(h: np.ndarray, n: int = 1024) -> tuple[np.ndarray, np.ndarray]:
    """``(omega, group delay in samples)`` of the FIR ``h`` over ``[0, pi)``.

    Computed analytically as ``Re{ sum(i h_i e^{-jwi}) / sum(h_i e^{-jwi}) }``.
    Even-length symmetric kernels give the constant ``(L-1)/2``; the other
    kinds give a frequency-dependent delay (Figs. 6-11).
    """
    h = np.asarray(h, dtype=float)
    w = np.linspace(0.0, np.pi, n, endpoint=False)
    E = np.exp(-1j * np.outer(w, np.arange(len(h))))
    num = E @ (np.arange(len(h)) * h)
    den = E @ h
    with np.errstate(divide="ignore", invalid="ignore"):
        gd = np.real(num / den)
    return w, gd
