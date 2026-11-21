"""Spectral behaviour of the whole estimator chain on synthetic signals: the data behind fig29.

The chain under test is the one of the RTL: two 512-point FFTs per segment, the six-tap fractional-delay merge that
produces the odd bins of the 1 024-point segment spectrum, the frequency-domain window, the squarer and the averaging
shift.  Every study runs the bit-accurate model (``psd_fixed``, word length 16 unless stated) next to the exact Welch
average of the same record (``psd_exact_welch``) and, to separate causes, the floating-point six-tap merge
(``psd_float_fd``).  The odd-bin merge error itself is quantified in ``docs/PAPER_NOTES.md`` section 2.2; the studies here
show what it does to line spectra, to the choice of window and to the averaging.

Nothing here is measured hardware; every record is synthetic and seeded.  All spectra are on the Welch scale of
``psd_exact_welch`` and are compared as dB relative to the peak of the *exact* spectrum of the same record.
"""
from __future__ import annotations

import numpy as np

from .estimator import CH, SEG, HwConfig, psd_exact_welch, psd_fixed, psd_float_fd
from .windows import window as make_window

__all__ = [
    "WINDOWS", "WL", "NCHUNKS", "TONE_AMP", "TONE_PHASE", "EXCLUDE", "BAND", "FLOOR_DB", "TONE_BIN", "OFFSETS",
    "SEPARATIONS", "VALLEY_F1", "VALLEY_RATIOS_DB", "AVG_SEGMENTS", "AVG_NREAL", "AVG_SEED0", "AVG_SIGMA", "AVG_BAND",
    "CSV_HEADER", "MODEL_LABELS",
    "tone", "two_tones", "psds", "db_re", "leakage_floor_db", "valley_depth_db", "floor_vs_offset", "valley_vs_separation",
    "welch_overlap_rho", "welch_overlap_variance", "noise_record", "averaging_study", "tone_panel", "run_study", "study_rows",
]

WINDOWS = ("rect", "hann", "hamming", "blackman")
WL = 16  # word length of the hardware model: long enough that quantisation does not show
NCHUNKS = 9  # nine 512-sample chunks = eight overlapping 1 024-sample segments
TONE_AMP = 0.9  # tone amplitude, fraction of full scale
TONE_PHASE = 0.3  # rad
EXCLUDE = 12  # bins around a tone that are left out of a leakage floor
BAND = (4, 508)  # one-sided analysis band (bins of the 1 024-point grid)
FLOOR_DB = -100.0  # clamp for bins whose power is zero or below the clamp
TONE_BIN = 300.5  # tone of panel (a): half way between two bins of the 1 024-point grid
OFFSETS = tuple(np.arange(0, 17) / 8.0)  # tone at 300 + offset, offset 0..2 bins in steps of 1/8
SEPARATIONS = tuple(range(2, 11))  # bins between two tones
VALLEY_F1 = 200.3  # first tone of the two-tone study (the second is at VALLEY_F1 + separation)
VALLEY_RATIOS_DB = (0.0, -20.0)  # amplitude of the second tone relative to the first
AVG_SEGMENTS = (1, 2, 4, 8, 16, 32, 64)
AVG_NREAL = 48  # records per averaging point
AVG_SEED0 = 1000
AVG_SIGMA = 0.12  # Gaussian record, standard deviation as a fraction of full scale (no clipping)
AVG_BAND = (8, 504)
MODEL_LABELS = {"ref": "exact_welch", "fd": "float_six_tap", "hw": f"hw_wl{WL}"}
CSV_HEADER = ["study", "window", "model", "x", "y"]


# ---------------------------------------------------------------------------------------------------------------
# records
# ---------------------------------------------------------------------------------------------------------------
def tone(f_bins: float, amp: float = TONE_AMP, phase: float = TONE_PHASE, nchunks: int = NCHUNKS) -> np.ndarray:
    """Cosine of ``f_bins`` bins of the 1 024-point grid (cycles per 1 024 samples)."""
    n = np.arange(nchunks * CH)
    return amp * np.cos(2 * np.pi * f_bins / SEG * n + phase)


def two_tones(f1: float, f2: float, ratio_db: float = 0.0, amp: float = 0.45, nchunks: int = NCHUNKS) -> np.ndarray:
    """Two cosines; the second is ``ratio_db`` dB relative to the first (amplitude ``amp`` for the first)."""
    n = np.arange(nchunks * CH)
    a2 = 10.0 ** (ratio_db / 20.0)
    return amp * (np.cos(2 * np.pi * f1 / SEG * n + 0.4) + a2 * np.cos(2 * np.pi * f2 / SEG * n + 1.1))


def noise_record(seed: int, nchunks: int, sigma: float = AVG_SIGMA) -> np.ndarray:
    """Gaussian white record (seeded), ``sigma`` as a fraction of full scale."""
    return sigma * np.random.default_rng(seed).standard_normal(nchunks * CH)


# ---------------------------------------------------------------------------------------------------------------
# spectra and figures of merit
# ---------------------------------------------------------------------------------------------------------------
def psds(x, window: str = "hann", wl: int = WL) -> dict[str, np.ndarray]:
    """Three estimates of the same record, ``(1024,)`` each on the Welch scale.

    ``ref``: exact Welch average; ``fd``: floating-point six-tap merge (the error of the merge, no quantisation);
    ``hw``: bit-accurate hardware model at word length ``wl``.
    """
    return {
        "ref": psd_exact_welch(x, window),
        "fd": psd_float_fd(x, window),
        "hw": psd_fixed(x, HwConfig(wl=wl, window=window)).to_float(),
    }


def db_re(p, peak: float, floor_db: float = FLOOR_DB):
    """``10 log10(p / peak)``, clamped at ``floor_db`` (zero power becomes the clamp)."""
    with np.errstate(divide="ignore"):
        v = 10.0 * np.log10(np.asarray(p, dtype=float) / peak)
    return np.maximum(v, floor_db)


def leakage_floor_db(p, f_bins: float, peak: float, exclude: int = EXCLUDE, band=BAND) -> float:
    """Highest bin more than ``exclude`` bins from the tone inside the one-sided ``band``, in dB relative to ``peak``."""
    p = np.asarray(p, dtype=float)
    k = np.arange(len(p))
    sel = (np.abs(k - f_bins) > exclude) & (k >= band[0]) & (k <= band[1])
    return float(db_re(p[sel].max(), peak))


def valley_depth_db(p, f1: float, f2: float) -> float:
    """Depth of the valley between two peaks in dB: the lower of the two peak values (maximum over the bin of the tone and its
    two neighbours) over the lowest bin strictly between the two peak bins.  NaN if no bin lies between them."""
    p = np.asarray(p, dtype=float)
    b1, b2 = int(round(f1)), int(round(f2))
    if b2 - b1 < 2:
        return float("nan")
    top = min(p[b1 - 1:b1 + 2].max(), p[b2 - 1:b2 + 2].max())
    mid = p[b1 + 1:b2].min()
    return float(min(10.0 * np.log10(top / max(mid, 1e-300)), -FLOOR_DB))


def floor_vs_offset(window: str, offsets=OFFSETS, base: float = 300.0, wl: int = WL) -> dict[str, np.ndarray]:
    """Leakage floor (dB re the exact peak) of a tone at ``base + offset`` for the three estimates (keys ``ref``, ``fd``, ``hw``).

    ``base`` is an even bin: offset 0 is a tone on an even bin of the 1 024-point grid, offset 1 on an odd bin, offset 1/2 half way
    between two bins.
    """
    out: dict[str, list[float]] = {"ref": [], "fd": [], "hw": []}
    for off in offsets:
        f = base + float(off)
        P = psds(tone(f), window, wl)
        peak = P["ref"][: SEG // 2].max()
        for key in out:
            out[key].append(leakage_floor_db(P[key], f, peak))
    return {k: np.array(v) for k, v in out.items()}


def valley_vs_separation(window: str, separations=SEPARATIONS, f1: float = VALLEY_F1, ratio_db: float = 0.0,
                         wl: int = WL) -> dict[str, np.ndarray]:
    """Valley depth (dB) between two tones ``separation`` bins apart for the three estimates (keys ``ref``, ``fd``, ``hw``)."""
    out: dict[str, list[float]] = {"ref": [], "fd": [], "hw": []}
    for d in separations:
        f2 = f1 + float(d)
        P = psds(two_tones(f1, f2, ratio_db), window, wl)
        for key in out:
            out[key].append(valley_depth_db(P[key][: SEG // 2], f1, f2))
    return {k: np.array(v) for k, v in out.items()}


# ---------------------------------------------------------------------------------------------------------------
# averaging
# ---------------------------------------------------------------------------------------------------------------
def welch_overlap_rho(window: str, hop: int = CH) -> float:
    """Correlation coefficient of the windows of two neighbouring segments: ``|sum w[n] w[n + hop]| / sum w[n]^2``."""
    w = make_window(window, SEG)
    return float(abs(np.sum(w[: SEG - hop] * w[hop:])) / np.sum(w * w))


def welch_overlap_variance(window: str, segments: int) -> float:
    """Relative variance (variance over squared mean) of a Welch average of ``segments`` segments of Gaussian white noise
    at 50 % overlap: ``(1 + 2 (1 - 1/S) rho^2) / S`` with ``rho`` of :func:`welch_overlap_rho`.  Only neighbouring segments
    overlap at hop 512, so a single correlation term remains; valid for bins away from DC and Nyquist."""
    rho = welch_overlap_rho(window)
    return (1.0 + 2.0 * (1.0 - 1.0 / segments) * rho * rho) / segments


def averaging_study(window: str = "hann", segments=AVG_SEGMENTS, nreal: int = AVG_NREAL, wl: int = WL,
                    seed0: int = AVG_SEED0, band=AVG_BAND) -> dict[str, np.ndarray]:
    """Relative standard deviation of the PSD estimate against the number of segments ``S``.

    ``nreal`` seeded Gaussian records of ``max(segments) + 1`` chunks; the first ``S + 1`` chunks of each give an ``S``-segment
    estimate.  For every bin in ``band`` the standard deviation over the records is divided by the mean over the records, and
    the root mean square over the bins is reported.  Keys: ``S``, ``hw``, ``ref`` (measured, hardware model and exact average),
    ``analytic`` (square root of :func:`welch_overlap_variance`), ``inv_sqrt`` (``1 / sqrt(S)``), ``bias_db`` (mean over records
    and bins of the hardware estimate against the exact average, in dB).
    """
    segs = tuple(int(s) for s in segments)
    k0, k1 = band
    hw = {s: [] for s in segs}
    ref = {s: [] for s in segs}
    top = max(segs) + 1
    for r in range(nreal):
        x = noise_record(seed0 + r, top)
        for s in segs:
            xs = x[: (s + 1) * CH]
            hw[s].append(psd_fixed(xs, HwConfig(wl=wl, window=window)).to_float()[k0:k1 + 1])
            ref[s].append(psd_exact_welch(xs, window)[k0:k1 + 1])

    def rel_std(a):
        a = np.asarray(a)
        return float(np.sqrt(np.mean((a.std(axis=0, ddof=1) / a.mean(axis=0)) ** 2)))

    return {
        "S": np.array(segs),
        "hw": np.array([rel_std(hw[s]) for s in segs]),
        "ref": np.array([rel_std(ref[s]) for s in segs]),
        "analytic": np.array([np.sqrt(welch_overlap_variance(window, s)) for s in segs]),
        "inv_sqrt": 1.0 / np.sqrt(np.array(segs, dtype=float)),
        "bias_db": np.array([10.0 * np.log10(np.mean(hw[s]) / np.mean(ref[s])) for s in segs]),
    }


# ---------------------------------------------------------------------------------------------------------------
# all studies and the table behind the sheet
# ---------------------------------------------------------------------------------------------------------------
def tone_panel(f_bins: float = TONE_BIN, window: str = "hann", wl: int = WL) -> dict:
    """Panel (a): the three estimates of one tone, in dB relative to the exact peak, for bins ``0..512``."""
    P = psds(tone(f_bins), window, wl)
    peak = float(P["ref"][: SEG // 2].max())
    return {
        "f": f_bins, "window": window, "peak": peak,
        "bins": np.arange(SEG // 2 + 1),
        **{k: db_re(P[k][: SEG // 2 + 1], peak) for k in P},
        "floor": {k: leakage_floor_db(P[k], f_bins, peak) for k in P},
    }


def run_study(nreal: int = AVG_NREAL) -> dict:
    """Run the four studies of fig29 with the module's parameters."""
    return {
        "tone": tone_panel(),
        "floor": {w: floor_vs_offset(w) for w in WINDOWS},
        "valley": {r: {w: valley_vs_separation(w, ratio_db=r) for w in WINDOWS} for r in VALLEY_RATIOS_DB},
        "averaging": averaging_study("hann", nreal=nreal),
    }


def study_rows(res: dict) -> list[list]:
    """Long-format rows ``[study, window, model, x, y]`` for ``results/model_spectral_study.csv`` from :func:`run_study`."""
    rows: list[list] = []
    t = res["tone"]
    for key in ("ref", "fd", "hw"):
        for b, v in zip(t["bins"], t[key]):
            rows.append(["tone_psd_db_re_exact_peak", t["window"], MODEL_LABELS[key], int(b), round(float(v), 4)])
    for w in WINDOWS:
        for key in ("ref", "fd", "hw"):
            for off, v in zip(OFFSETS, res["floor"][w][key]):
                rows.append(["floor_db_re_exact_peak_vs_offset", w, MODEL_LABELS[key], float(off), round(float(v), 4)])
    for r in VALLEY_RATIOS_DB:
        for w in WINDOWS:
            for key in ("ref", "fd", "hw"):
                for d, v in zip(SEPARATIONS, res["valley"][r][w][key]):
                    rows.append([f"valley_db_vs_separation_ratio{int(r)}db", w, MODEL_LABELS[key], int(d), round(float(v), 4)])
    a = res["averaging"]
    for name, key in (("hw", "hw"), ("exact_welch", "ref"), ("analytic_overlap", "analytic"), ("inv_sqrt_S", "inv_sqrt"),
                      ("hw_bias_db", "bias_db")):
        label = MODEL_LABELS["hw"] if name == "hw" else name
        for s, v in zip(a["S"], a[key]):
            rows.append(["averaging_rel_std" if name != "hw_bias_db" else "averaging_bias_db", "hann", label, int(s), round(float(v), 5)])
    return rows
