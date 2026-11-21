"""Numbers behind the FFT figure sheets (fig24 index cube, fig25 pass snapshots, fig27 twiddle atlas, fig28 accuracy).

Everything is computed from the models in this package on synthetic signals (``btpsd.r8fft``, ``btpsd.schedule``,
``btpsd.cordic``); nothing is measured hardware.  The functions are separate from the drawing code so that
``python/tests/test_fft_study.py`` can pin the numbers down: the snapshots end in the model's FFT, the operation counts
equal an explicit enumeration, and the analytic noise prediction is compared with the simulated error.

Fixed-point conventions are those of :mod:`btpsd.r8fft`: ``wl``-bit words, every radix-2 add/subtract followed by a
rounded ``>> 1`` (so each radix-8 butterfly scales by 1/8 and the whole transform by 1/512).
"""

from __future__ import annotations

import numpy as np

from .cordic import RotRom, default_rom, octant_reduce, rotate_float
from .r8fft import dft8_fixed, dft8_float, to_fixed
from .cordic import rotate_fixed
from .rtlmodel import LAT, LAT_BFLY, LAT_CRU, row_fn, stall_cycles, twq
from .schedule import N, NPASS, bank_of, output_index, pass_addresses, pass_twiddle_q, row_of

__all__ = [
    "pass_snapshots_fixed", "pass_snapshots_float", "fft_float_rom", "natural_order", "digit_table",
    "rotation_counts", "operation_counts", "fft_inputs", "fft_sqnr", "predicted_error_variance", "predicted_bias_power",
    "predicted_error_power", "predicted_sqnr", "rom_only_sqnr", "combined_sqnr", "error_floor_per_bin", "pass_error_power",
    "rom_error_vectors", "twiddle_exponents_used", "tone_addresses", "ACCURACY_HEADER", "ACCURACY_INPUTS", "ACCURACY_WLS",
    "accuracy_rows", "TWO_TONE_BINS", "two_tone_record", "mean_error_per_bin", "sqnr_vs_level", "LEVELS_DB",
    "butterfly_signals", "fft_cycle_budget",
]

DIGIT_LABELS = [("n2", "n1", "n0"), ("k0", "n1", "n0"), ("k0", "k1", "n0"), ("k0", "k1", "k2")]


# ---------------------------------------------------------------------------------------------------------------
# memory snapshots after every pass
# ---------------------------------------------------------------------------------------------------------------
def pass_snapshots_fixed(re, im, rom: RotRom, wl: int, rounding: str = "round") -> list[np.ndarray]:
    """Memory content (complex integer codes, indexed by address ``a``) after load and after each of the three passes.

    The last entry holds the digit-reversed spectrum of :func:`btpsd.r8fft.r8fft_fixed`; :func:`natural_order` maps it
    to frequency order.
    """
    m_re = np.array(re, dtype=np.int64, copy=True)
    m_im = np.array(im, dtype=np.int64, copy=True)
    snaps = [m_re + 1j * m_im]
    for p in range(NPASS):
        adr = pass_addresses(p)
        yr, yi = dft8_fixed(m_re[adr], m_im[adr], rom, wl, rounding)
        yr, yi = rotate_fixed(yr, yi, pass_twiddle_q(p), rom, wl, rounding)
        m_re[adr] = yr
        m_im[adr] = yi
        snaps.append(m_re + 1j * m_im)
    return snaps


def pass_snapshots_float(x) -> list[np.ndarray]:
    """The same algorithm in floating point with exact twiddles and the 1/8 scaling of every pass."""
    m = np.array(x, dtype=complex, copy=True)
    snaps = [m.copy()]
    for p in range(NPASS):
        adr = pass_addresses(p)
        y = dft8_float(m[adr]) / 8.0
        m[adr] = rotate_float(y, pass_twiddle_q(p))
        snaps.append(m.copy())
    return snaps


def _dft8_float_rom(x, rom: RotRom):
    """Scaled (1/8) 8-point DFT in floating point whose three ``W8^i`` rotations use the ROM coefficients."""
    x = np.asarray(x, dtype=complex)
    s = (x[..., :4] + x[..., 4:]) / 2.0
    d = rotate_float((x[..., :4] - x[..., 4:]) / 2.0, np.array([0, 128, 256, 384]), rom)
    out = np.empty_like(x)
    for grp, base in ((s, 0), (d, 1)):
        c0, c1 = (grp[..., 0] + grp[..., 2]) / 2, (grp[..., 1] + grp[..., 3]) / 2
        c2, c3 = (grp[..., 0] - grp[..., 2]) / 2, (grp[..., 1] - grp[..., 3]) / 2 * (-1j)
        out[..., 0 + base] = (c0 + c1) / 2
        out[..., 4 + base] = (c0 - c1) / 2
        out[..., 2 + base] = (c2 + c3) / 2
        out[..., 6 + base] = (c2 - c3) / 2
    return out


def fft_float_rom(x, rom: RotRom | None = None) -> list[np.ndarray]:
    """Floating-point snapshots (as :func:`pass_snapshots_float`) with the ROM coefficients in every rotation.

    No rounding anywhere: the difference to the exact-twiddle snapshots is the error of the rotator coefficients alone.
    """
    rom = rom or default_rom()
    m = np.array(x, dtype=complex, copy=True)
    snaps = [m.copy()]
    for p in range(NPASS):
        adr = pass_addresses(p)
        m[adr] = rotate_float(_dft8_float_rom(m[adr], rom), pass_twiddle_q(p), rom)
        snaps.append(m.copy())
    return snaps


def natural_order(mem) -> np.ndarray:
    """Spectrum in frequency order from the memory content after pass 2 (address ``a`` holds ``X[output_index(a)]``)."""
    mem = np.asarray(mem)
    out = np.empty_like(mem)
    out[output_index(np.arange(N))] = mem
    return out


def digit_table() -> list[tuple[str, tuple[str, str, str]]]:
    """Meaning of the address digits ``(a2, a1, a0)`` after load and after each pass."""
    names = ["after load", "after pass 0", "after pass 1", "after pass 2"]
    return list(zip(names, DIGIT_LABELS))


def tone_addresses(bins) -> list[dict]:
    """Where the words of the given frequency bins sit after pass 2: digits, address, bank and row."""
    rows = []
    for k in bins:
        k0, k1, k2 = k & 7, (k >> 3) & 7, (k >> 6) & 7
        a = 64 * k0 + 8 * k1 + k2  # (a2, a1, a0) = (k0, k1, k2)
        assert int(output_index(a)) == k
        rows.append({"k": int(k), "k0": k0, "k1": k1, "k2": k2, "addr": int(a), "bank": int(bank_of(a)), "row": int(row_of(a))})
    return rows


# ---------------------------------------------------------------------------------------------------------------
# operation counts
# ---------------------------------------------------------------------------------------------------------------
def rotation_counts() -> dict:
    """Rotations that need more than a swap/negate (``q`` not a multiple of 256 in units of pi/512).

    ``direct``: the N^2 products of the DFT; ``radix2``: the twiddles of an in-place radix-2 DIF FFT;
    ``radix8``: the ``W8^i`` and output twiddles of this implementation (pass 2 has no twiddle).
    """
    n = np.arange(N)
    e = (n[:, None] * n[None, :]) % N
    direct = int(np.count_nonzero(e % (N // 4)))
    radix2 = 0
    for s in range(1, 10):
        j = np.arange(N // 2)
        exp = (j % (N >> s)) << (s - 1)  # twiddle exponent of butterfly j in stage s (W_N^exp)
        radix2 += int(np.count_nonzero(exp % (N // 4)))
    per_pass = [int(np.count_nonzero(pass_twiddle_q(p) % 256)) for p in range(NPASS)]
    w8 = 2 * 64 * NPASS  # W8^1 and W8^3 of every butterfly; W8^2 = -j is free
    return {"direct": direct, "radix2": radix2, "twiddle_per_pass": per_pass, "w8": w8, "radix8": sum(per_pass) + w8}


def operation_counts() -> dict:
    """Arithmetic and memory traffic of one 512-point transform: this implementation against an in-place radix-2 FFT.

    ``complex_addsub``: complex additions/subtractions (24 per radix-8 butterfly: 3 stages of 4 radix-2 butterflies, each
    one addition and one subtraction); ``words_*``: words moved to and from the memory; ``passes``: sweeps over the whole memory.
    """
    rc = rotation_counts()
    bfly = 64 * NPASS
    r2 = 256 * 9
    return {**rc, "butterflies": bfly, "complex_addsub": 24 * bfly, "words_read": 8 * bfly,
            "words_written": 8 * bfly, "passes": NPASS, "radix2_butterflies": r2, "radix2_addsub": 2 * r2,
            "radix2_words_read": 2 * r2, "radix2_words_written": 2 * r2, "radix2_passes": 9}


# ---------------------------------------------------------------------------------------------------------------
# accuracy
# ---------------------------------------------------------------------------------------------------------------
def fft_inputs(kind: str, wl: int, seed: int = 0):
    """Real test records quantised to ``wl``-bit codes (peak 0.9 of full scale), returned as ``(re, im)`` integer arrays."""
    from . import signals as sg

    n = np.arange(N)
    if kind == "white noise":
        x = sg.white(N, seed=seed)
    elif kind == "full-scale tone":
        x = 0.9 * np.cos(2 * np.pi * 37 * n / N + 0.3)
    elif kind == "six tones":
        rng = np.random.default_rng(seed)
        bins = rng.choice(np.arange(3, 250), 6, replace=False)
        x = sum(np.cos(2 * np.pi * b * n / N + rng.uniform(0, 6.28)) for b in bins)
        x = sg.scale_peak(x, 0.9)
    elif kind == "EEG-like":
        x = sg.eeg_like(N, seed=seed)
    else:
        raise ValueError(kind)
    re, im = to_fixed(x + 0j, wl)
    return re, im


def fft_sqnr(re, im, wl: int, rom: RotRom | None = None, rounding: str = "round") -> float:
    """SQNR (dB) of the fixed-point FFT against the exact transform of the same integer input, scaled by 1/512."""
    rom = rom or default_rom()
    snap = pass_snapshots_fixed(re, im, rom, wl, rounding)[-1]
    got = natural_order(snap)
    ref = np.fft.fft(np.asarray(re, float) + 1j * np.asarray(im, float)) / N
    return 10 * np.log10(np.sum(np.abs(ref) ** 2) / np.sum(np.abs(got - ref) ** 2))


def _nontrivial_fraction(p: int) -> float:
    return float(np.count_nonzero(pass_twiddle_q(p) % 256)) / pass_twiddle_q(p).size


def predicted_error_variance() -> float:
    """Predicted error variance per real component of every output word, in LSB^2 (random part, rounding to nearest).

    * the nine rounded halvings add a half-LSB error with probability 1/2 (variance 1/16) after stage ``s``; the
      ``9 - s`` radix-2 stages that follow average it down by 1/2 each;
    * a rotator rounds once (uniform, variance 1/12) on the elements that pass through a non-trivial rotation:
      a quarter of the elements after stages 1, 4 and 7 (the two 45-degree ``W8^i``) and the fractions of
      ``pass_twiddle_q`` after stages 3 and 6.
    """
    halving = sum((1.0 / 16.0) * 2.0 ** -(9 - s) for s in range(1, 10))
    w8 = sum(0.25 * (1.0 / 12.0) * 2.0 ** -m for m in (8, 5, 2))
    tw = (_nontrivial_fraction(0) * (1.0 / 12.0) * 2.0**-6) + (_nontrivial_fraction(1) * (1.0 / 12.0) * 2.0**-3)
    return halving + w8 + tw


def predicted_bias_power() -> float:
    """Predicted power of the systematic error per bin (complex, LSB^2) for rounding to nearest.

    Round-half-up of ``(a +- b) / 2`` has a mean error of +1/4 LSB per real component.  The offset stays on the sum
    branches of the following stages and cancels on the difference branches, so after the last stage a bin collects
    ``n`` offsets, ``n = 1`` for every bin, ``n >= 2`` for half of them, and so on (``P(n >= j) = 2^-(j-1)``).  Adding the
    offsets in phase, the power is ``2 (1/4)^2 E[n^2]`` per complex bin, ``E[n^2] = sum (2j - 1) 2^-(j-1)``, j = 1..9.
    The twiddle rotations turn the offsets of different stages against each other, so this is an upper estimate.
    """
    en2 = sum((2 * j - 1) * 2.0 ** -(j - 1) for j in range(1, 10))
    return 2.0 * (1.0 / 16.0) * en2


def predicted_error_power() -> float:
    """Predicted mean-square error per complex bin (LSB^2), random part plus systematic part, rounding to nearest."""
    return 2.0 * predicted_error_variance() + predicted_bias_power()


def predicted_sqnr(re, im) -> float:
    """SQNR (dB) predicted from :func:`predicted_error_power` for the integer input ``re + j im`` (quantisation only).

    The estimate assumes that the roundings behave like random errors, which holds for noise-like records; for a
    record with few non-zero spectral lines the measured SQNR is higher (see ``tests/test_fft_study.py``).
    """
    ref = np.fft.fft(np.asarray(re, float) + 1j * np.asarray(im, float)) / N
    return 10 * np.log10(np.sum(np.abs(ref) ** 2) / (N * predicted_error_power()))


def rom_only_sqnr(re, im, rom: RotRom | None = None) -> float:
    """SQNR (dB) of the floating-point FFT with the ROM coefficients (no rounding anywhere) against the exact transform."""
    rom = rom or default_rom()
    x = np.asarray(re, float) + 1j * np.asarray(im, float)
    ref = np.fft.fft(x) / N
    got = natural_order(fft_float_rom(x, rom)[-1])
    return 10 * np.log10(np.sum(np.abs(ref) ** 2) / np.sum(np.abs(got - ref) ** 2))


def combined_sqnr(quant_db: float, rom_db: float) -> float:
    """Two independent error powers add: ``1 / (10^(-q/10) + 10^(-r/10))`` in dB."""
    return -10 * np.log10(10 ** (-quant_db / 10) + 10 ** (-rom_db / 10))


def error_floor_per_bin(re, im, wl: int, rom: RotRom | None = None, rounding: str = "round") -> np.ndarray:
    """``|X_fixed - X_exact|^2`` per bin (LSB^2), natural order."""
    rom = rom or default_rom()
    got = natural_order(pass_snapshots_fixed(re, im, rom, wl, rounding)[-1])
    ref = np.fft.fft(np.asarray(re, float) + 1j * np.asarray(im, float)) / N
    return np.abs(got - ref) ** 2


def pass_error_power(re, im, wl: int, rom: RotRom | None = None, rounding: str = "round") -> dict:
    """Mean-square error per word (LSB^2) after load, pass 0, 1, 2 of the fixed-point FFT.

    ``total``: fixed point against the exact floating-point snapshot of the same pass; ``rom``: the floating-point
    snapshot with the ROM coefficients against the exact one (coefficient error alone); ``arith``: fixed point against the
    ROM snapshot (rounding of the datapath alone).
    """
    rom = rom or default_rom()
    x = np.asarray(re, float) + 1j * np.asarray(im, float)
    fx = pass_snapshots_fixed(re, im, rom, wl, rounding)
    ex = pass_snapshots_float(x)
    fr = fft_float_rom(x, rom)
    ms = lambda a, b: float(np.mean(np.abs(a - b) ** 2))  # noqa: E731
    return {"total": [ms(f, e) for f, e in zip(fx, ex)], "rom": [ms(r, e) for r, e in zip(fr, ex)],
            "arith": [ms(f, r) for f, r in zip(fx, fr)], "signal": [float(np.mean(np.abs(e) ** 2)) for e in ex]}


# ---------------------------------------------------------------------------------------------------------------
# twiddle atlas
# ---------------------------------------------------------------------------------------------------------------
def rom_error_vectors(rom: RotRom | None = None) -> dict:
    """Error of the ROM rotation against the exact ``exp(-j q pi/512)`` for every ``q = 0..1023``."""
    rom = rom or default_rom()
    q = np.arange(1024)
    got = rotate_float(np.ones(1024, dtype=complex), q, rom)
    exact = np.exp(-1j * q * np.pi / 512)
    err = got - exact
    _, r, _ = octant_reduce(q)
    return {"q": q, "r": r, "err": err, "abs": np.abs(err), "angle": np.angle(got / exact), "gain": np.abs(got) - 1.0}


def twiddle_exponents_used() -> dict:
    """Exponents ``e`` of ``W512^e`` (``q = 2e`` in units of pi/512) that the output twiddles of each pass use, and those of ``W8^i``."""
    out = {}
    for p in range(2):
        q = pass_twiddle_q(p)
        out[f"pass {p}"] = sorted({int(v) // 2 for v in q.ravel()})
    out["W8"] = [64 * i for i in range(8)]
    return out


# ---------------------------------------------------------------------------------------------------------------
# table behind the accuracy sheet (also stored in results/model_fft_accuracy.csv)
# ---------------------------------------------------------------------------------------------------------------
ACCURACY_HEADER = ["input", "wl", "rounding", "sqnr_fixed_db", "sqnr_predicted_db", "sqnr_quantisation_only_predicted_db",
                   "sqnr_rom_only_db"]
ACCURACY_INPUTS = ("white noise", "full-scale tone", "six tones", "EEG-like")
ACCURACY_WLS = (8, 10, 12, 14, 16, 18, 20, 24)


def accuracy_rows(seed: int = 1) -> list[list]:
    """Rows ``[input, wl, rounding, fixed SQNR, predicted SQNR, quantisation-only prediction, ROM-only SQNR]`` (dB, 0.01).

    The predictions are for rounding to nearest and are left empty for truncation.
    """
    rom = default_rom()
    rows = []
    for kind in ACCURACY_INPUTS:
        for wl in ACCURACY_WLS:
            re, im = fft_inputs(kind, wl, seed)
            r_db = rom_only_sqnr(re, im, rom)
            q_db = predicted_sqnr(re, im)
            for rounding in ("round", "trunc"):
                pred = [round(float(combined_sqnr(q_db, r_db)), 2), round(float(q_db), 2)] if rounding == "round" else ["", ""]
                rows.append([kind, wl, rounding, round(float(fft_sqnr(re, im, wl, rom, rounding)), 2), *pred, round(float(r_db), 2)])
    return rows


# ---------------------------------------------------------------------------------------------------------------
# records used by the sheets
# ---------------------------------------------------------------------------------------------------------------
TWO_TONE_BINS = (37, 150)


def two_tone_record(wl: int = 12, seed: int = 1, snr_db: float = 40.0):
    """Two tones (bins 37 and 150, amplitudes 1 and 1/2) in white noise, quantised to ``wl`` bits; returns ``(re, im)``.

    The SNR is the power of the two tones over the power of the noise before scaling to a peak of 0.9 of full scale.
    """
    from . import signals as sg

    rng = np.random.default_rng(seed)
    n = np.arange(N)
    x = sum(a * np.cos(2 * np.pi * b * n / N + rng.uniform(0, 2 * np.pi)) for b, a in zip(TWO_TONE_BINS, (1.0, 0.5)))
    x = x + rng.standard_normal(N) * np.sqrt(np.mean(x**2) / 10 ** (snr_db / 10))
    return to_fixed(sg.scale_peak(x, 0.9) + 0j, wl)


def mean_error_per_bin(wl: int = 12, rounding: str = "round", records: int = 60, kind: str = "white noise", rom: RotRom | None = None):
    """Error of the fixed-point FFT per bin over ``records`` records: ``(mean error, mean-square error)`` as complex / real arrays (LSB)."""
    rom = rom or default_rom()
    errs = []
    for seed in range(records):
        re, im = fft_inputs(kind, wl, seed)
        got = natural_order(pass_snapshots_fixed(re, im, rom, wl, rounding)[-1])
        errs.append(got - np.fft.fft(np.asarray(re, float) + 1j * np.asarray(im, float)) / N)
    e = np.array(errs)
    return e.mean(axis=0), np.mean(np.abs(e) ** 2, axis=0)


LEVELS_DB = tuple(range(-66, -11, 6))


def sqnr_vs_level(wl: int, levels_db=LEVELS_DB, rounding: str = "round", seed: int = 1, rom: RotRom | None = None) -> list[float]:
    """SQNR (dB) of the fixed-point FFT for white noise whose rms is ``level`` dB below full scale (``2**(wl-1)`` codes).

    The error of the arithmetic does not depend on the level, so the SQNR falls by one decibel per decibel of level until
    the coefficient error of the rotator ROM takes over at the top.
    """
    rom = rom or default_rom()
    rng = np.random.default_rng(seed)
    g = rng.standard_normal(N)
    out = []
    for lev in levels_db:
        re, im = to_fixed(g * 10 ** (lev / 20) + 0j, wl)
        out.append(float(fft_sqnr(re, im, wl, rom, rounding)))
    return out


# ---------------------------------------------------------------------------------------------------------------
# the signals of one butterfly (rtl/fft512_mem.v)
# ---------------------------------------------------------------------------------------------------------------
def butterfly_signals(p: int, b: int) -> dict:
    """What ``fft512_mem`` derives for butterfly ``b = 8 hi + lo`` of pass ``p``.

    ``s`` is the rotation ``(hi + lo) mod 8`` that routes bank ``(m + s) mod 8`` to butterfly slot ``m``.  ``slots`` lists,
    for every slot ``m``, the bank, the row (``row_fn`` of the RTL), the address digits ``(a2, a1, a0)`` recovered from
    ``bank = (a2 + a1 + a0) mod 8`` and ``row = 8 a2 + a1``, the address ``a``, and the twiddle exponent ``tq`` (units of
    pi/512, ``twq`` of the RTL) with its angle in degrees.  Slot 0 never gets a twiddle.
    """
    if not (0 <= p < NPASS and 0 <= b < 64):
        raise ValueError("pass 0..2, butterfly 0..63")
    hi, lo = b >> 3, b & 7
    s = (hi + lo) & 7
    slots = []
    for m in range(8):
        bank, row = (m + s) & 7, row_fn(p, m, hi, lo)
        a2, a1 = row >> 3, row & 7
        a0 = (bank - a2 - a1) & 7
        t = twq(p, b, m)
        slots.append(dict(m=m, bank=bank, row=row, a2=a2, a1=a1, a0=a0, addr=64 * a2 + 8 * a1 + a0, tq=t,
                          angle_deg=t * 180.0 / 512.0))
    return dict(p=p, b=b, hi=hi, lo=lo, s=s, slots=slots)


def fft_cycle_budget() -> dict:
    """Latency constants of the core and the transform time.

    ``total_cycles`` counts from the cycle of the first issue to the cycle in which ``done`` becomes visible (the
    ``cycles`` of ``tb/tb_fft512_mem.v``, ``schedule_cycles`` of the model, ``Fft512Rtl.run``): the 196 issue slots
    (3 x 64 butterflies plus the stall after pass 0) and the ``LAT - 1`` cycles until the last write is readable.
    In the CI trace ``results/rtl_trace_fft512.csv`` (cycle 0 = ``start`` high, signals sampled at the rising edge) the
    first issue is sampled in cycle 1, so ``done`` is sampled in cycle ``1 + total_cycles``.
    """
    stall = stall_cycles()
    issue_slots = NPASS * 64 + stall
    total = issue_slots + LAT - 1
    return dict(lat_cru=LAT_CRU, lat_bfly=LAT_BFLY, lat=LAT, stall=stall, issue_slots=issue_slots, total_cycles=total,
                trace_first_issue=1, trace_done_cycle=1 + total)
