"""Data helpers shared by the figure scripts (model outputs on synthetic signals; nothing here is measured hardware)."""
import numpy as np

from _common import RESULTS
from btpsd.estimator import HwConfig, hw_power_scale, psd_exact_welch, psd_fixed
from btpsd.fdfilter import PAPER_FD6, PAPER_FD6_START, fd_kernel, merge_chunk_spectra
from btpsd.metrics import sqnr_db
from btpsd.reference import periodogram
from btpsd.signals import ar2_lowpass, eeg_like, qpc, white


def read_csv(name):
    import csv
    with open(RESULTS / name) as f:
        return list(csv.DictReader(f))


def fd_gain(h, start, n=512):
    """Complex gain G(n0) of the odd-bin interpolator for a unit impulse at sample n0 of the chunk difference.

    The merge computes O(m) = -j (-1)^m sum_i h_i (-1)^(m+d_i) D(m+d_i); for d[n] = delta[n - n0] this is G(n0) D(m).
    The ideal half-bin value is D(m + 1/2) = exp(-j pi n0 / n) D(m).
    """
    n0 = np.arange(n)
    d = start + np.arange(len(h))
    ph = np.exp(-2j * np.pi * np.outer(n0, d) / n)
    G = -1j * (ph * (np.asarray(h)[None, :] * ((-1.0) ** d)[None, :])).sum(axis=1)
    return n0, G, np.exp(-1j * np.pi * n0 / n)


def check_fd_gain():
    """The closed form must equal the merge code applied to an impulse (float, mode='fd')."""
    n0 = 100
    a = np.zeros(512)
    a[n0] = 1.0
    A = np.fft.fft(a)[None, :]
    X = merge_chunk_spectra(A, np.zeros_like(A), h=PAPER_FD6, start=PAPER_FD6_START, mode="fd")[0]
    r = X[1::2] / A[0]
    _, G, _ = fd_gain(PAPER_FD6, PAPER_FD6_START)
    return float(np.max(np.abs(r - G[n0])))


def exact_kernel(dmin=-9, dmax=10):
    """Ideal (infinite) half-bin kernel after sign alternation: h_d = -(-1)^d / (pi (d - 1/2))."""
    d = np.arange(dmin, dmax + 1)
    return d, -((-1.0) ** d) / (np.pi * (d - 0.5))


def psd_example(wl=12, nchunks=9, seed=3, window="hann"):
    x = eeg_like(512 * nchunks, fs=256.0, seed=seed)
    ref = psd_exact_welch(x, window)
    res = psd_fixed(x, HwConfig(wl=wl, window=window))
    hw = res.to_float()
    return x, ref, hw, sqnr_db(ref, hw)


def bispec_pair(signal="qpc_coupled", wl=14, window="hann", nchunks=60, seed=0):
    """``(reference, hardware model)`` bispectra on f1, f2 < 128 for a synthetic signal.

    The reference is the floating-point direct estimator (512-point segments, hop 256); the model is
    ``bispectrum_fixed`` scaled back by ``bispectrum_hw_to_float``.
    """
    from btpsd import bispectrum as B
    from btpsd import signals as sg

    n = 512 * nchunks
    gens = {"qpc_coupled": lambda: sg.qpc(n, coupled=True, seed=seed)[0],
            "qpc_uncoupled": lambda: sg.qpc(n, coupled=False, seed=seed)[0],
            "pink": lambda: sg.pink(n, seed=seed), "eeg_like": lambda: sg.eeg_like(n, seed=seed),
            "eeg_coupled": lambda: sg.eeg_coupled(n, seed=seed)}
    x = gens[signal]()
    ref = B.bispectrum_direct(x, L=512, hop=256, nfft=512, max_bin=128, window=window)
    re, im, S, sh = B.bispectrum_fixed(x, B.BispecConfig(wl=wl, window=window))
    return ref, B.bispectrum_hw_to_float(re, im, wl, S, sh)


def squarer_cycles(wl, nchunks=9, seed=3, window="hann"):
    """Cycles of the sum-of-odds squarer pair (re, im) for every bin 0..512 of every segment (EEG-like signal)."""
    from btpsd.squarer import mag_sq

    x = eeg_like(512 * nchunks, fs=256.0, seed=seed)
    res = psd_fixed(x, HwConfig(wl=wl, window=window))
    _, cyc = mag_sq(res.X_re[:, :513], res.X_im[:, :513], units=2)
    return np.asarray(cyc).ravel()


def psd_examples(wl: int = 10, window_name: str = "hann", nchunks: int = 10, seed: int = 3, fs: float = 256.0):
    """Synthetic records with their float and hardware-model PSDs (``wl``-bit datapath); sheet ``fig22_psd_examples``."""
    n = 512 * nchunks
    sigs = {"EEG-like": eeg_like(n, fs=fs, seed=seed), "AR(2) narrow-band": ar2_lowpass(n, seed=seed), "white noise": white(n, seed=seed)}
    out = {}
    for name, x in sigs.items():
        ref = psd_exact_welch(x, window_name)
        res = psd_fixed(x, HwConfig(wl=wl, window=window_name))
        hw = res.to_float()
        out[name] = {"x": x, "ref": ref, "hw": hw, "P_int": res.P_int, "lsb": hw_power_scale(wl, res.S, res.avg_shift),
                     "blk": periodogram(x[:1024], window=window_name), "sqnr": sqnr_db(ref, hw)}
    return out


PSD_EXAMPLE_HEADER = ["signal", "window", "wl", "record_samples", "seed", "sqnr_vs_exact_welch_db", "psd_lsb_db", "zero_bins_of_513"]


def psd_example_rows(wl: int = 10, window_name: str = "hann", nchunks: int = 10, seed: int = 3):
    """Numbers printed on ``fig22_psd_examples`` (rows of ``results/model_psd_examples_wl10.csv``)."""
    rows = []
    for name, e in psd_examples(wl, window_name, nchunks, seed).items():
        rows.append([name, window_name, wl, 512 * nchunks, seed, round(float(e["sqnr"]), 3), round(float(10 * np.log10(e["lsb"])), 3),
                     int((e["P_int"][:513] == 0).sum())])
    return rows
