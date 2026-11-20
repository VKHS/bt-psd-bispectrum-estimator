"""Model experiments -> ``results/*.csv`` (+ ``results/provenance.json``).

Everything here is produced by *this repository's models on synthetic signals*.
Nothing is copied from the paper; paper-reported values live in ``data/paper_*.csv``.
SQNR = 10 log10(sum P_ref^2 / sum (P_ref - P_hw)^2), linear power domain, analytic scale (no fitting).
"""
import csv
import json
import platform
import time

import numpy as np
import scipy

from _common import RESULTS
from btpsd import bispectrum as B
from btpsd import estimator as E
from btpsd import signals as sg
from btpsd.cordic import default_rom
from btpsd.fdfilter import fd_kernel, merge_chunk_spectra
from btpsd.metrics import snr_complex_db, sqnr_db
from btpsd.r8fft import fft512_reference_scaled, r8fft_fixed, to_fixed
from btpsd.schedule import schedule_cycles

NCH = 40  # 512-sample chunks per record
N = 512 * NCH
SEEDS = (0, 1, 2, 3, 4)


def brown(n, seed=0, peak=0.9):
    w = np.random.default_rng(seed).standard_normal(n + 2000)
    x = np.cumsum(w)[2000:]
    x = x - np.linspace(x[0], x[-1], n)  # remove the linear trend
    return sg.scale_peak(x, peak)


SIG = {
    "white": sg.white,
    "pink": sg.pink,
    "brown": brown,
    "eeg_like": sg.eeg_like,
    "ar2_narrowband": sg.ar2_lowpass,
}


def write(name, header, rows):
    with open(RESULTS / name, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
    print(f"  wrote results/{name} ({len(rows)} rows)")


def merge_snr(x, h, st):
    A = np.fft.fft(x[: (len(x) // 512) * 512].reshape(-1, 512), axis=-1)
    ref = merge_chunk_spectra(A[:-1], A[1:], mode="exact")[:, 1::2]
    out = merge_chunk_spectra(A[:-1], A[1:], h=h, start=st, mode="fd")[:, 1::2]
    return snr_complex_db(ref, out)


def exp_table3_analog():
    rows = []
    for name, gen in SIG.items():
        for kind in ("sinc", "lagrange"):
            for L in range(3, 11):
                h, st = fd_kernel(L, kind)
                s_c, s_p = [], []
                for sd in SEEDS:
                    x = gen(N, seed=sd)
                    s_c.append(merge_snr(x, h, st))
                    s_p.append(sqnr_db(E.psd_exact_welch(x, "hann"), E.psd_float_fd(x, "hann", h=h, start=st)))
                rows.append([name, kind, L, round(float(np.median(s_c)), 3), round(float(np.median(s_p)), 3)])
        h, st = fd_kernel(6, "paper6")
        s_c, s_p = [], []
        for sd in SEEDS:
            x = gen(N, seed=sd)
            s_c.append(merge_snr(x, h, st))
            s_p.append(sqnr_db(E.psd_exact_welch(x, "hann"), E.psd_float_fd(x, "hann", h=h, start=st)))
        rows.append([name, "paper6", 6, round(float(np.median(s_c)), 3), round(float(np.median(s_p)), 3)])
    write("fd_merge_snr_vs_taps.csv",
          ["signal", "kernel", "taps", "snr_complex_spectrum_db", "sqnr_psd_hann_db"], rows)


def exp_sqnr_vs_wl():
    rows = []
    for tol in (0.01, 0.0025):
        for name in ("white", "pink", "eeg_like", "ar2_narrowband"):
            for wl in (10, 12, 14, 16, 18, 20, 24):
                tot, arith, merge_lim = [], [], []
                for sd in SEEDS:
                    x = SIG[name](N, seed=sd)
                    ref = E.psd_exact_welch(x, "hann")
                    fd = E.psd_float_fd(x, "hann")
                    hw = E.psd_fixed(x, E.HwConfig(wl=wl, rom_tol=tol)).to_float()
                    tot.append(sqnr_db(ref, hw))
                    arith.append(sqnr_db(fd, hw))
                    merge_lim.append(sqnr_db(ref, fd))
                rows.append([name, tol, wl, round(float(np.median(tot)), 3), round(float(np.median(arith)), 3),
                             round(float(np.median(merge_lim)), 3)])
    write("sqnr_vs_wordlength.csv",
          ["signal", "rom_angle_tol_rad", "wl", "sqnr_vs_exact_welch_db", "sqnr_vs_float_fd6_db",
           "sqnr_float_fd6_vs_exact_db"], rows)


def exp_rom_tradeoff():
    rng = np.random.default_rng(0)
    rows = []
    for tol in (0.02, 0.01, 0.005, 0.0025, 0.001):
        rom = default_rom(tol=tol)
        ae = np.abs(rom.angle_error()[1:])
        snr = {}
        for wl in (16, 24):
            x = 0.5 * (rng.uniform(-1, 1, (8, 512)) + 1j * rng.uniform(-1, 1, (8, 512)))
            re, im = to_fixed(x, wl)
            fr, fi = r8fft_fixed(re, im, rom, wl)
            snr[wl] = snr_complex_db(fft512_reference_scaled(re, im), fr + 1j * fi)
        rows.append([tol, round(float(rom.adders()[1:].mean()), 3), round(float(ae.max()), 6),
                     round(float(np.sqrt((ae**2).mean())), 6), round(float(snr[16]), 2), round(float(snr[24]), 2)])
    write("cru_rom_tradeoff.csv",
          ["angle_tol_rad", "mean_adders_per_rotation", "max_angle_err_rad", "rms_angle_err_rad",
           "fft512_snr_wl16_db", "fft512_snr_wl24_db"], rows)


def exp_bispectrum():
    rows = []
    nb = 512 * 60
    sigs = {"qpc_coupled": sg.qpc(nb, coupled=True, seed=0)[0], "qpc_uncoupled": sg.qpc(nb, coupled=False, seed=0)[0],
            "pink": sg.pink(nb, seed=0)}
    for name, x in sigs.items():
        ref = B.bispectrum_direct(x, L=512, hop=256, nfft=512, max_bin=128)
        for wl in (10, 12, 14, 16, 18):
            re, im, S, sh = B.bispectrum_fixed(x, B.BispecConfig(wl=wl))
            hw = B.bispectrum_hw_to_float(re, im, wl, S, sh)
            rows.append([name, wl, round(float(snr_complex_db(ref, hw)), 3),
                         round(float(sqnr_db(np.abs(ref), np.abs(hw))), 3)])
    write("bispectrum_sqnr_vs_wordlength.csv", ["signal", "wl", "snr_complex_db", "sqnr_magnitude_db"], rows)


def exp_cycle_model():
    rows = [[lat, *[schedule_cycles(lat)[k] for k in ("issue_cycles", "total_cycles")],
             sum(schedule_cycles(lat)["stalls_per_pass"])] for lat in range(4, 17)]
    write("fft512_cycles_vs_latency.csv", ["pipeline_latency_cycles", "issue_cycles", "total_cycles", "stall_cycles"], rows)
    rows = []
    cyc_budget = 512 / 0.88  # cycles per chunk at the paper's Table 9 'proposed' ratio R = 0.88 samples/cycle
    for name in ("white", "pink", "eeg_like", "ar2_narrowband"):
        for wl in (10, 16):
            res = E.psd_fixed(SIG[name](N, seed=0), E.HwConfig(wl=wl, window="hann"))
            from btpsd.squarer import mag_sq
            _, cyc = mag_sq(res.X_re[:, :513], res.X_im[:, :513], units=2)  # real input: bins 0..512 suffice
            per_seg = cyc.sum(axis=1)
            rows.append([name, wl, round(float(per_seg.mean()), 1), round(float(per_seg.max()), 1),
                         round(float(cyc.max()), 1), round(float(per_seg.mean() / cyc_budget), 2)])
    write("squarer_cycle_model.csv",
          ["signal", "wl", "mean_cycles_per_segment_513_bins", "max_cycles_per_segment", "max_cycles_single_bin",
           "serial_squarer_units_needed_for_R0.88"], rows)


def exp_wl10_windows():
    rows = []
    for name in ("white", "pink", "eeg_like", "ar2_narrowband"):
        for win in ("rect", "hamming", "hann", "blackman"):
            vals = []
            for sd in SEEDS:
                x = SIG[name](N, seed=sd)
                vals.append(sqnr_db(E.psd_exact_welch(x, win), E.psd_fixed(x, E.HwConfig(wl=10, window=win)).to_float()))
            rows.append([name, win, round(float(np.median(vals)), 3)])
    write("model_sqnr_wl10_synthetic.csv", ["signal", "window", "sqnr_vs_exact_welch_db"], rows)


if __name__ == "__main__":
    t0 = time.time()
    RESULTS.mkdir(exist_ok=True)
    for f in (exp_table3_analog, exp_sqnr_vs_wl, exp_rom_tradeoff, exp_bispectrum, exp_cycle_model, exp_wl10_windows):
        print(f.__name__)
        f()
    with open(RESULTS / "provenance.json", "w") as f:
        json.dump({
            "generated_by": "python/scripts/run_experiments.py",
            "numpy": np.__version__, "scipy": scipy.__version__, "python": platform.python_version(),
            "record_length_samples": N, "seeds": list(SEEDS), "signals": "synthetic (btpsd.signals); not the paper's data",
            "sqnr_definition": "10*log10(sum(P_ref^2)/sum((P_ref-P_hw)^2)), linear power domain, analytic scale",
            "note": "These numbers are model results on synthetic signals. Paper-reported values are in data/paper_*.csv.",
            "runtime_s": round(time.time() - t0, 1),
        }, f, indent=2)
    print(f"done in {time.time() - t0:.1f}s")
