"""Memory map and butterfly schedule of the radix-8 memory-based 512-point FFT.

This is the repository's own in-place scheme; the paper's organisation (register sequences, Fig. 3 and Tables 1-2) is
different, see docs/ARCHITECTURE.md section 6.

Data layout
-----------
A 512-word chunk lives in 8 banks of 64 words (one read and one write port each).  With the address
split ``a = 64*a2 + 8*a1 + a0`` (digits ``a2, a1, a0`` in 0..7)::

    bank(a) = (a2 + a1 + a0) mod 8        row(a) = 8*a2 + a1

A radix-8 butterfly of pass ``p`` combines the 8 words that differ in exactly
one digit (pass 0: ``a2``, pass 1: ``a1``, pass 2: ``a0``).  Because the bank
index is the digit sum modulo 8, those 8 words always sit in 8 *different*
banks: every butterfly reads and writes all banks in parallel, one word each,
with no conflict (checked in the tests).  The butterfly is *in place*.

Algorithm (in-place radix-8 DIF, natural-order input)
-----------------------------------------------------
With ``n = 64 n2 + 8 n1 + n0`` and ``k = k0 + 8 k1 + 64 k2``::

    pass 0   8-pt DFT over n2 -> digit k0,  twiddle W512^{(8 n1 + n0) k0}
    pass 1   8-pt DFT over n1 -> digit k1,  twiddle W512^{8 n0 k1}
    pass 2   8-pt DFT over n0 -> digit k2,  no twiddle

After pass 2 the word at address ``(a2, a1, a0)`` is ``X[a2 + 8 a1 + 64 a0]``
(digit reversal on the read-out address).

Twiddle exponents are returned in units of pi/512 (``W512^e = exp(-j 2e * pi/512)``),
i.e. the index ``q`` of :mod:`btpsd.cordic`.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "N",
    "NBANK",
    "ROWS",
    "NPASS",
    "split",
    "bank_of",
    "row_of",
    "addr_of",
    "pass_addresses",
    "pass_twiddle_q",
    "output_index",
    "check_conflict_free",
    "MemoryFFT",
    "schedule_cycles",
]

N = 512
NBANK = 8
ROWS = 64
NPASS = 3


def split(a):
    """Digits ``(a2, a1, a0)`` of an address."""
    a = np.asarray(a, dtype=np.int64)
    return (a >> 6) & 7, (a >> 3) & 7, a & 7


def bank_of(a):
    a2, a1, a0 = split(a)
    return (a2 + a1 + a0) % 8


def row_of(a):
    a2, a1, _ = split(a)
    return 8 * a2 + a1


def addr_of(bank, row):
    """Inverse of ``(bank_of, row_of)``."""
    bank = np.asarray(bank, dtype=np.int64)
    row = np.asarray(row, dtype=np.int64)
    a2, a1 = row >> 3, row & 7
    a0 = (bank - a2 - a1) % 8
    return 64 * a2 + 8 * a1 + a0


def pass_addresses(p: int) -> np.ndarray:
    """``(64, 8)`` addresses: butterfly ``b`` (natural order) of pass ``p``, input slot ``m``."""
    b = np.arange(64)[:, None]
    m = np.arange(8)[None, :]
    hi, lo = b >> 3, b & 7
    if p == 0:  # fixed (a1, a0) = (hi, lo), vary a2
        return 64 * m + 8 * hi + lo
    if p == 1:  # fixed (a2, a0) = (hi, lo), vary a1
        return 64 * hi + 8 * m + lo
    if p == 2:  # fixed (a2, a1) = (hi, lo), vary a0
        return 64 * hi + 8 * lo + m
    raise ValueError("pass must be 0, 1 or 2")


def pass_twiddle_q(p: int) -> np.ndarray:
    """``(64, 8)`` twiddle index ``q`` (units of pi/512) applied to butterfly output ``k``."""
    b = np.arange(64)[:, None]
    k = np.arange(8)[None, :]
    if p == 0:  # W512^{b k}  with b = 8 n1 + n0
        return (2 * b * k) % 1024
    if p == 1:  # W512^{8 n0 k}, n0 = b mod 8
        return (16 * (b & 7) * k) % 1024
    if p == 2:
        return np.zeros((64, 8), dtype=np.int64)
    raise ValueError("pass must be 0, 1 or 2")


def output_index(a):
    """Frequency index held at address ``a`` after pass 2: ``a2 + 8 a1 + 64 a0``."""
    a2, a1, a0 = split(a)
    return a2 + 8 * a1 + 64 * a0


def check_conflict_free() -> bool:
    """True if every butterfly of every pass touches 8 distinct banks and 8 distinct rows/banks pairs."""
    for p in range(NPASS):
        adr = pass_addresses(p)
        if sorted(adr.ravel().tolist()) != list(range(N)):
            return False
        bk = bank_of(adr)
        if any(len(set(r.tolist())) != 8 for r in bk):
            return False
    return True


# ---------------------------------------------------------------------------
# Explicit-bank execution (golden model of the memory organisation)
# ---------------------------------------------------------------------------
class MemoryFFT:
    """Executes the in-place radix-8 FFT on explicit 8x64 banks (float, exact twiddles).

    Used to verify that the bank/row mapping realises the algorithm and to
    produce the address traces that the RTL address generators must match.
    """

    def __init__(self) -> None:
        self.banks = np.zeros((NBANK, ROWS), dtype=complex)
        self.trace: list[tuple[int, int, list[int], list[int]]] = []

    def load(self, x) -> None:
        x = np.asarray(x, dtype=complex)
        a = np.arange(N)
        self.banks[bank_of(a), row_of(a)] = x

    def read(self, a):
        return self.banks[bank_of(a), row_of(a)]

    def write(self, a, v) -> None:
        self.banks[bank_of(a), row_of(a)] = v

    def unload(self) -> np.ndarray:
        """Natural-order spectrum (read-out address bit-reversed in digits)."""
        out = np.empty(N, dtype=complex)
        a = np.arange(N)
        out[output_index(a)] = self.read(a)
        return out

    def run(self) -> np.ndarray:
        from .r8fft import dft8_float
        from .cordic import rotate_float

        self.trace.clear()
        for p in range(NPASS):
            adr = pass_addresses(p)
            tw = pass_twiddle_q(p)
            for b in range(64):
                a = adr[b]
                banks = bank_of(a)
                rows = row_of(a)
                self.trace.append((p, b, banks.tolist(), rows.tolist()))
                y = dft8_float(self.read(a))
                self.write(a, rotate_float(y, tw[b]))
        return self.unload()


# ---------------------------------------------------------------------------
# Throughput model (RAW hazards between passes)
# ---------------------------------------------------------------------------
def schedule_cycles(latency: int) -> dict:
    """Issue schedule of the 3*64 butterflies in counter (natural) order.

    One butterfly issues per cycle (all 8 banks read in parallel); its 8
    results are written ``latency`` cycles after issue.  A butterfly of
    pass p+1 may only issue once the producers of its 8 words have written.
    Returns the total cycle count (first issue to last write) and the stall
    cycles inserted per pass; the stalls come out as ``max(0, latency - 8)``
    at the pass 0 -> 1 boundary and none elsewhere.
    """
    if latency < 1:
        raise ValueError("latency >= 1")
    write_done: dict[int, int] = {}
    t = 0
    stalls = []
    for p in range(NPASS):
        adr = pass_addresses(p)
        pass_stall = 0
        for b in range(64):
            ready = max((write_done.get(int(a), 0) for a in adr[b]), default=0)
            issue = max(t, ready)
            pass_stall += issue - t
            for a in adr[b]:
                write_done[int(a)] = issue + latency
            t = issue + 1
        stalls.append(pass_stall)
    return {"latency": latency, "total_cycles": int(max(write_done.values())),
            "issue_cycles": int(t), "stalls_per_pass": stalls}
