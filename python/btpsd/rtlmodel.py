"""Python transliteration of the RTL control logic.

``rtl/fft512_mem.v`` and ``rtl/psd_engine.v`` are written from the formulas implemented here, line by line:
the issue schedule with its stall, the bank rotation, the row functions, the twiddle exponents, the write-back
timing and the unload address map (:class:`Fft512Rtl`), and the sequential merge / window / accumulate loops of the
segment engine (:func:`psd_segment_loops`).  The tests compare them with the bit-accurate models, so a mistake in the
*design* of the RTL (an address, a hazard, a tap order, a sign) shows up without an HDL simulator.  A typing mistake
in the Verilog itself is only caught by the testbenches in ``tb/``.

The RTL pipeline is modelled at the level that matters for correctness: a butterfly reads the banks in its issue cycle
(read-first RAM: words written at the end of cycle ``t`` are visible to reads from cycle ``t + 1``), and its results
are written ``LAT_BFLY + 1`` cycles later.
"""

from __future__ import annotations

import numpy as np

from .cordic import default_rom, rotate_fixed
from .estimator import HwConfig
from .fixedpoint import rshift, saturate
from .r8fft import dft8_fixed

__all__ = ["LAT_CRU", "LAT_BFLY", "LAT", "stall_cycles", "row_fn", "twq", "Fft512Rtl", "psd_segment_loops",
           "bispectrum_frame_loops", "PSD_STATES", "BSP_STATES", "psd_engine_state_cycles", "bispectrum_engine_state_cycles"]

LAT_CRU = 4  # register stages of rtl/cru_rotate.v
LAT_BFLY = 1 + LAT_CRU + 1 + LAT_CRU  # rtl/r8_butterfly.v
LAT = LAT_BFLY + 2  # issue -> cycle in which the written words can be read again


def stall_cycles(lat: int = LAT) -> int:
    """Idle issue cycles after pass 0 (``STALL`` of rtl/fft512_mem.v)."""
    return lat - 8 if lat > 8 else 0


def row_fn(p: int, m: int, hi: int, lo: int) -> int:
    """Row of slot ``m`` of butterfly ``(hi, lo)`` in pass ``p`` (``row_fn`` of rtl/fft512_mem.v)."""
    if p == 0:
        return (m << 3) | hi
    if p == 1:
        return (hi << 3) | m
    return (hi << 3) | lo


def twq(p: int, b: int, k: int) -> int:
    """Twiddle exponent (units of pi/512) of butterfly output ``k`` (``twq`` of rtl/fft512_mem.v)."""
    if p == 0:
        return ((b * k) << 1) & 1023
    if p == 1:
        return (((b & 7) * k) << 4) & 1023
    return 0


class Fft512Rtl:
    """Cycle-level model of ``fft512_mem`` (8 banks of 64 words, one butterfly issue per cycle)."""

    def __init__(self, wl: int = 16, rounding: str = "round", tol: float = 0.01) -> None:
        self.wl, self.rounding = wl, rounding
        self.rom = default_rom(tol)
        self.mem = np.zeros((8, 64, 2), dtype=np.int64)
        self.done_cycle = 0

    # -- load / unload ---------------------------------------------------------
    def load(self, re, im) -> None:
        for n in range(512):
            a2, a1, a0 = n >> 6, (n >> 3) & 7, n & 7
            self.mem[(a2 + a1 + a0) & 7, (a2 << 3) | a1] = (re[n], im[n])

    def read_bin(self, k: int):
        a2, a1, a0 = k & 7, (k >> 3) & 7, k >> 6
        w = self.mem[(a2 + a1 + a0) & 7, (a2 << 3) | a1]
        return int(w[0]), int(w[1])

    # -- run ---------------------------------------------------------------------
    def issue_schedule(self, stall: int):
        t, sched = 0, []
        for p in range(3):
            for b in range(64):
                sched.append((t, p, b))
                t += 1
            if p == 0:
                t += stall
        return sched

    def run(self, stall: int | None = None, lat: int = LAT) -> int:
        """Execute the three passes; returns the cycle in which ``done`` is visible (relative to the first issue)."""
        stall = stall_cycles(lat) if stall is None else stall
        sched = self.issue_schedule(stall)
        writes: dict[int, list] = {}
        last = sched[-1][0]
        for cyc in range(last + lat + 1):
            # reads of the issue in this cycle see the memory before this cycle's writes
            issue = [s for s in sched if s[0] == cyc]
            for _, p, b in issue:
                hi, lo = b >> 3, b & 7
                s = (hi + lo) & 7
                words = np.zeros((8, 2), dtype=np.int64)  # words[j]: what bank j delivers
                for j in range(8):
                    m_j = (j - s) & 7
                    words[j] = self.mem[j, row_fn(p, m_j, hi, lo)]
                x = np.array([words[(m + s) & 7] for m in range(8)])  # slot m <- bank (m + s) mod 8
                yr, yi = dft8_fixed(x[:, 0], x[:, 1], self.rom, self.wl, self.rounding)
                tq = np.array([0] + [twq(p, b, k) for k in range(1, 8)])
                yr, yi = rotate_fixed(yr, yi, tq, self.rom, self.wl, self.rounding)
                writes.setdefault(cyc + lat - 1, []).append((p, b, yr.copy(), yi.copy()))
            for p, b, yr, yi in writes.pop(cyc, []):
                hi, lo = b >> 3, b & 7
                s = (hi + lo) & 7
                for j in range(8):
                    m_w = (j - s) & 7
                    self.mem[j, row_fn(p, m_w, hi, lo)] = (yr[m_w], yi[m_w])
        self.done_cycle = last + lat
        return self.done_cycle

    def spectrum(self):
        out = [self.read_bin(k) for k in range(512)]
        return np.array([o[0] for o in out], dtype=np.int64), np.array([o[1] for o in out], dtype=np.int64)


def psd_segment_loops(Ar, Ai, Br, Bi, cfg: HwConfig):
    """The loops of ``psd_engine.v`` for one segment: returns the windowed spectrum ``(Wr, Wi)`` (1024 ints each)."""
    wl, rnd = cfg.wl, cfg.rounding
    h = [int(v) for v in cfg.fd_int()]
    win = cfg.win_int()
    wcoef = [int(win.get(d, 0)) for d in (-2, -1, 0, 1, 2)]
    fd_sh = cfg.fd_frac + 1
    xr = np.zeros(1024, dtype=np.int64)
    xi = np.zeros(1024, dtype=np.int64)
    # merge: six taps per m, tap t reads index (m + t - 2) mod 512
    for m in range(512):
        facc_r = facc_i = 0
        for t in range(6):
            idx = (m + t - 2) & 511
            dre, dim = int(Ar[idx]) - int(Br[idx]), int(Ai[idx]) - int(Bi[idx])
            if idx & 1:
                dre, dim = -dre, -dim
            facc_r += dre * h[t]
            facc_i += dim * h[t]
            if t == 2:  # d = 0
                xr[2 * m] = int(rshift(int(Ar[m]) + int(Br[m]), 1, rnd))
                xi[2 * m] = int(rshift(int(Ai[m]) + int(Bi[m]), 1, rnd))
        ow_r, ow_i = (-facc_i, facc_r) if (m & 1) else (facc_i, -facc_r)
        xr[2 * m + 1] = int(saturate(rshift(ow_r, fd_sh, rnd), wl))
        xi[2 * m + 1] = int(saturate(rshift(ow_i, fd_sh, rnd), wl))
    # window: five taps per k, tap t reads index (k + 2 - t) mod 1024 with coefficient W(t - 2)
    wr = np.zeros(1024, dtype=np.int64)
    wi = np.zeros(1024, dtype=np.int64)
    for k in range(1024):
        ar = ai = 0
        for t in range(5):
            idx = (k + 2 - t) & 1023
            ar += int(xr[idx]) * wcoef[t]
            ai += int(xi[idx]) * wcoef[t]
        wr[k] = int(saturate(rshift(ar, cfg.win_frac, rnd), wl))
        wi[k] = int(saturate(rshift(ai, cfg.win_frac, rnd), wl))
    return wr, wi


def bispectrum_frame_loops(Fr, Fi, cfg, logb: int):
    """The loops of ``bispectrum_engine.v`` for one frame spectrum ``(Fr, Fi)`` (512 ints each, scale 1/512).

    Returns the contribution ``(T_re, T_im)`` (``(2**logb, 2**logb)`` int64) that is added to the accumulators.
    """
    maxb, nb, wl, rnd = 1 << logb, 2 << logb, cfg.wl, cfg.rounding
    win = HwConfig(window=cfg.window, win_frac=cfg.win_frac).win_int()
    wcoef = [int(win.get(d, 0)) for d in (-2, -1, 0, 1, 2)]
    fw_r = np.zeros(nb, dtype=np.int64)
    fw_i = np.zeros(nb, dtype=np.int64)
    for k in range(nb):
        ar = ai = 0
        for t in range(5):
            idx = (k + 2 - t) & 511
            ar += int(Fr[idx]) * wcoef[t]
            ai += int(Fi[idx]) * wcoef[t]
        fw_r[k] = int(saturate(rshift(ar, cfg.win_frac, rnd), wl))
        fw_i[k] = int(saturate(rshift(ai, cfg.win_frac, rnd), wl))
    tre = np.zeros((maxb, maxb), dtype=np.int64)
    tim = np.zeros((maxb, maxb), dtype=np.int64)
    for f1 in range(maxb):
        for f2 in range(maxb):
            ar, ai, br, bi = int(fw_r[f1]), int(fw_i[f1]), int(fw_r[f2]), int(fw_i[f2])
            cr, ci = int(fw_r[f1 + f2]), int(fw_i[f1 + f2])
            pr = int(rshift(ar * br - ai * bi, wl - 1, rnd))
            pi_ = int(rshift(ar * bi + ai * br, wl - 1, rnd))
            tre[f1, f2] = pr * cr + pi_ * ci
            tim[f1, f2] = pi_ * cr - pr * ci
    return tre, tim


# ---- cycles per FSM state (formulas read off the case arms of rtl/psd_engine.v and rtl/bispectrum_engine.v) -------------
PSD_STATES = ("S_IDLE", "S_CLR", "M_A", "M_B", "M_C", "M_F", "W_A", "W_B", "W_C", "W_F", "W_G", "W_S", "W_ACC", "S_END")
BSP_STATES = ("S_IDLE", "S_CLR", "W_A", "W_B", "W_C", "W_F", "R_A", "R_B", "R_C", "T_A", "T_B", "T_C", "T_D", "S_END")


def psd_engine_state_cycles(segments: int, squarer_cycles: int, clears: int = 1, serial: bool = True) -> dict[str, int]:
    """Cycles spent in each busy state of ``psd_engine`` for ``clears`` clear operations and ``segments`` segments.

    ``squarer_cycles`` is the sum over all bins of all segments of ``max(|re|, |im|)`` (:func:`btpsd.squarer.mag_sq`).
    Per segment: merge ``512 x (6 taps x 3 cycles + M_F)``, window ``1024 x (5 taps x 3 cycles + W_F + W_G + W_S + W_ACC)``.
    ``W_S`` is entered with the cycle in which ``sq_start`` is high, so it lasts ``n + 3`` cycles for a serial squarer
    (``done`` is visible ``n + 2`` cycles after the start cycle) and 2 cycles for the multiplier variant.  ``S_IDLE`` is
    not modelled (it is the time spent waiting).
    """
    wsum = (3 * 1024 * segments + squarer_cycles) if serial else 2 * 1024 * segments
    return {
        "S_CLR": 1024 * clears,
        "M_A": 3072 * segments, "M_B": 3072 * segments, "M_C": 3072 * segments, "M_F": 512 * segments,
        "W_A": 5120 * segments, "W_B": 5120 * segments, "W_C": 5120 * segments, "W_F": 1024 * segments,
        "W_G": 1024 * segments, "W_S": wsum, "W_ACC": 1024 * segments,
        "S_END": clears + segments,
    }


def bispectrum_engine_state_cycles(frames: int, logb: int, clears: int = 1) -> dict[str, int]:
    """Cycles spent in each busy state of ``bispectrum_engine`` for ``clears`` clears and ``frames`` frames.

    ``MAXB = 2**logb``.  Per frame: window ``2 MAXB x (5 taps x 3 cycles + W_F)``, ``MAXB`` visits of ``R_A .. R_C`` and
    ``MAXB**2`` visits of ``T_A .. T_D`` (``5216`` cycles for ``MAXB = 32``).
    """
    maxb = 1 << logb
    nb = 2 * maxb
    return {
        "S_CLR": maxb * maxb * clears,
        "W_A": 5 * nb * frames, "W_B": 5 * nb * frames, "W_C": 5 * nb * frames, "W_F": nb * frames,
        "R_A": maxb * frames, "R_B": maxb * frames, "R_C": maxb * frames,
        "T_A": maxb * maxb * frames, "T_B": maxb * maxb * frames, "T_C": maxb * maxb * frames, "T_D": maxb * maxb * frames,
        "S_END": clears + frames,
    }
