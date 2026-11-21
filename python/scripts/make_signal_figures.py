"""Figure sheets 21-22: window responses and signal-plus-PSD examples (counterparts of paper Figs. 12-15 and 18-20).

Everything is computed live from the model on synthetic signals (``btpsd.windows``, ``btpsd.estimator``); nothing is
measured hardware and nothing is taken from the paper's recordings.  The paper's own panels are kept unchanged in
``docs/paper_figures``.
"""
from __future__ import annotations

import csv

import numpy as np

from _common import FIGS, RESULTS
from _figdata import PSD_EXAMPLE_HEADER, psd_example_rows, psd_examples
from _sheet import PAPER, card, labels, legend
from btpsd.estimator import HwConfig
from btpsd.viz import GRID, INK, INK2, MUTED, PAL, Diagram, tint
from btpsd.windows import TABLE4_SIGNIFICANT_TAPS, kernel_exactness, window, window_metrics

BLUE, ORANGE, AQUA, VIOLET = PAL["blue"], PAL["orange"], PAL["aqua"], PAL["violet"]
WINS = ("rect", "hann", "hamming", "blackman")
WINCOL = {"rect": MUTED, "hann": BLUE, "hamming": ORANGE, "blackman": AQUA}
PAPER_MARK = dict(PAPER, ms=3.2, mew=0.9)
WINNAME = {"rect": "Rectangular", "hann": "Hann", "hamming": "Hamming", "blackman": "Blackman"}
COEF = {"rect": "a0 = 1", "hann": "a0 = a1 = 0.5", "hamming": "a0 = 0.54, a1 = 0.46", "blackman": "a0 = 0.42, a1 = 0.5, a2 = 0.08"}


def _fmt(v: float) -> str:
    return f"{v:.1f}".rstrip("0").rstrip(".")


def _title(ax, s):
    ax.set_title(s, fontsize=5.0, color=INK2, loc="left", pad=2.2)


def _sci(v: float) -> str:
    m, e = f"{v:.0e}".split("e")
    return f"{m}e{int(e)}"


# ---------------------------------------------------------------------------------------------------------------
# fig21: window responses
# ---------------------------------------------------------------------------------------------------------------
def fig21_window_responses():
    W, H = 150, 136
    D = Diagram(W, H, "Data windows of the estimator: impulse, spectrum and kernel (model)",
                subtitle="w[n] = a0 - a1 cos(2 pi n / L) + a2 cos(4 pi n / L), periodic form. Counterpart of the paper's Figs. 12-15, which show "
                         "128-point windows; the windows below are computed here.")
    L = 128
    n = np.arange(L)
    for i, name in enumerate(WINS):
        col = WINCOL[name]
        x0 = 2 + i * 37
        card(D, x0, 38, 35, 88, "(" + "abcd"[i] + ")", WINNAME[name], COEF[name], plot=False, fs=7.4)
        # ---- impulse response ---------------------------------------------------------------------------------------
        w = window(name, L)
        ax = D.inset(x0 + 7.0, 96.6, 26.0, 16.0, fs=4.4)
        ax.fill_between(n, w, color=tint(col, 0.82), lw=0)
        ax.plot(n, w, color=col, lw=1.0)
        ax.set_xlim(0, L - 1)
        ax.set_ylim(0, 1.1)
        ax.set_xticks([0, 64, 127])
        ax.set_yticks([0, 0.5, 1])
        _title(ax, "impulse response w[n], L = 128")
        labels(ax, "n", None, fs=4.8)
        # ---- spectrum: the few non-zero DFT bins -----------------------------------------------------------------------
        F = np.abs(np.fft.fft(w))
        d = np.arange(-3, 4)
        v = F[d % L]
        ax = D.inset(x0 + 7.0, 71.1, 26.0, 16.0, fs=4.4)
        ax.vlines(d, 0, v, color=col, lw=1.3)
        nz = v > 1e-9
        ax.plot(d[nz], v[nz], "s", mfc="#ffffff", mec=col, mew=1.0, ms=3.4)
        ax.plot(d[~nz], v[~nz], ".", color=MUTED, ms=2.4)
        for dd, vv in zip(d[nz], v[nz]):
            ax.text(dd, vv + 7.0, _fmt(vv), ha="center", va="bottom", fontsize=4.2, color=INK)
        ax.set_xlim(-3.6, 3.6)
        ax.set_ylim(0, 160)
        ax.set_xticks(d)
        ax.set_yticks([0, 64, 128])
        _title(ax, "|DFT of w[n]|, bins around DC")
        labels(ax, "bin offset d", None, fs=4.8)
        # ---- frequency response in dB ----------------------------------------------------------------------------------
        Lw = 1024
        pad = 64
        wl_ = window(name, Lw)
        R = np.abs(np.fft.fft(wl_, Lw * pad))
        R = R / R[0]
        f = np.fft.fftshift(np.fft.fftfreq(Lw * pad, d=1.0 / Lw))
        R = np.fft.fftshift(R)
        sel = np.abs(f) <= 8
        db = np.maximum(20 * np.log10(R[sel] + 1e-12), -90)
        ax = D.inset(x0 + 7.0, 45.6, 26.0, 16.0, fs=4.4)
        ax.plot(f[sel], db, color=col, lw=0.9)
        met = window_metrics(name)
        ax.axhline(met["psl_db"], color=INK2, lw=0.6, ls=(0, (2.5, 1.5)))
        ax.plot([met["psl_at"], -met["psl_at"]], [met["psl_db"]] * 2, "o", mfc="#ffffff", mec=INK, mew=0.8, ms=2.8, zorder=5)
        ax.text(7.6, met["psl_db"] + 3.0, f"{met['psl_db']:.1f} dB", ha="right", va="bottom", fontsize=4.4, color=INK)
        ax.set_xlim(-8, 8)
        ax.set_ylim(-90, 0)
        ax.set_xticks([-8, -4, 0, 4, 8])
        ax.set_yticks([-80, -60, -40, -20, 0])
        _title(ax, "response of the L = 1024 window (dB)")
        labels(ax, "frequency offset (bins)", None, fs=4.8)

    # ---- table ----------------------------------------------------------------------------------------------------------
    D.card(2, 2, 146, 34)
    D.panel(3.8, 32.7, "(e)", "Figures of merit and kernel constants of the four windows", fs=7.4)
    cols = [("window", 5.0, "left"), ("taps", 27.0, "center"), ("taps, paper (Table 4)", 45.0, "center"), ("kernel constants / 2^10", 66.0, "center"),
            ("coherent gain", 85.0, "center"), ("ENBW (bins)", 98.0, "center"), ("-3 dB width", 110.0, "center"),
            ("main lobe (+-bins)", 123.0, "center"), ("peak sidelobe", 137.0, "center"), ("peak sidelobe", 137.0, "center")][:-1]
    yh = 28.0
    for lab, xx, ha in cols:
        D.text(xx, yh, lab, fs=4.9, ha=ha, color=INK2, bold=True)
    D.ax.plot([4.0, 146.0], [yh - 1.4, yh - 1.4], color=GRID, lw=0.7, zorder=2)
    for r, name in enumerate(WINS):
        y = yh - 4.6 - 3.9 * r
        m = window_metrics(name)
        ints = HwConfig(window=name).win_int()
        kc = "  ".join(str(ints[k]) for k in sorted(ints))
        D.text(5.0, y, WINNAME[name], fs=5.4, ha="left", color=INK, bold=True)
        D.ax.plot([3.3], [y + 0.25], "s", color=WINCOL[name], ms=3.0, zorder=3)
        D.text(27.0, y, str(len(ints)), fs=5.4, ha="center", color=INK)
        D.text(45.0, y, str(TABLE4_SIGNIFICANT_TAPS[name]), fs=5.4, ha="center", color=INK)
        D.ax.plot([41.6], [y + 0.25], **PAPER_MARK)
        D.text(66.0, y, kc, fs=5.4, ha="center", color=INK)
        D.text(85.0, y, f"{m['cg']:.2f}", fs=5.4, ha="center", color=INK)
        D.text(98.0, y, f"{m['enbw']:.3f}", fs=5.4, ha="center", color=INK)
        D.text(110.0, y, f"{m['bw3']:.3f}", fs=5.4, ha="center", color=INK)
        D.text(123.0, y, f"{m['null']:.0f}", fs=5.4, ha="center", color=INK)
        D.text(137.0, y, f"{m['psl_db']:.1f} dB", fs=5.4, ha="center", color=INK)
    D.text(5.0, 6.4, "Computed here from the 1024-point periodic windows (64 x zero-padded DFT). The DFT of a periodic raised-cosine window has exactly "
           "1, 3, 3 or 5 non-zero bins, so the window is a short filter on the spectrum;", fs=4.6, ha="left", color=INK2)
    errs = {nm: kernel_exactness(nm, periodic=False) for nm in WINS[1:]}
    D.text(5.0, 3.8, "applying that kernel to the symmetric window (numpy.hanning form) leaves a relative error of "
           + ", ".join(f"{_sci(errs[nm])} ({nm})" for nm in WINS[1:]) + "; for the periodic window the difference is round-off only (below 1e-12).",
           fs=4.6, ha="left", color=INK2)
    D.save(FIGS / "fig21_window_responses")


# ---------------------------------------------------------------------------------------------------------------
# fig22: signal and PSD examples
# ---------------------------------------------------------------------------------------------------------------
def fig22_psd_examples():
    wl, win_name, fs = 10, "hann", 256.0
    ex = psd_examples(wl, win_name, fs=fs)
    W, H = 150, 128
    D = Diagram(W, H, "Signal and estimated PSD, word length 10 (model, synthetic signals)",
                subtitle=f"Counterpart of the paper's Figs. 18-20 (recorded EEG, EOG and EMG; its third curve, an ARMA(18,20) estimator, is not part of this "
                         f"repository). {win_name.capitalize()} window, fs = {fs:.0f} Hz nominal, 9 segments of 1 024 samples.")
    k = np.arange(513)
    fr = k * fs / 1024
    for r, (name, e) in enumerate(ex.items()):
        y0 = 2 + (2 - r) * 39
        D.card(2, y0, 146, 36)
        D.panel(3.8, y0 + 32.7, "(" + "abc"[r] + ")", name, fs=7.6)
        # ---- record --------------------------------------------------------------------------------------------------------
        ax = D.inset(10.0, y0 + 5.0, 46.0, 22.0, fs=4.6)
        ax.plot(np.arange(len(e["x"])), e["x"], color=BLUE, lw=0.35)
        ax.set_xlim(0, len(e["x"]))
        ax.set_ylim(-1, 1)
        ax.set_yticks([-1, 0, 1])
        _title(ax, "record (scaled to 0.9 of full scale)")
        labels(ax, "sample", None, fs=4.8)
        # ---- PSD -------------------------------------------------------------------------------------------------------------
        ax = D.inset(70.0, y0 + 5.0, 74.0, 22.0, fs=4.6)
        lsb_db = 10 * np.log10(e["lsb"])
        ref_db = 10 * np.log10(np.maximum(e["ref"][:513], 1e-300))
        top = np.ceil(max(ref_db.max(), lsb_db) / 10.0) * 10.0
        bot = np.floor((min(ref_db.min(), lsb_db) - 6.0) / 10.0) * 10.0
        db = lambda p: 10 * np.log10(np.maximum(p, 1e-300))  # noqa: E731
        nz = e["P_int"][:513] > 0
        ax.plot(fr, db(e["blk"][:513]), color="#b8b8b2", lw=0.6, label="float periodogram, one 1 024-sample block")
        ax.plot(fr, ref_db, color=INK, lw=0.8, label="float Welch average, 9 segments")
        ax.plot(fr, np.where(nz, db(e["hw"][:513]), np.nan), color=BLUE, lw=0.8, label=f"hardware model, WL = {wl}")
        ax.plot(fr[~nz], np.full(int((~nz).sum()), bot + 0.012 * (top - bot)), "|", color=BLUE, ms=2.6, mew=0.5,
                label="bins whose integer result is 0")
        ax.axhline(lsb_db, color=BLUE, lw=0.6, ls=(0, (1.2, 1.8)))
        ax.text(fs / 2 - 0.5, lsb_db + 0.012 * (top - bot), f"1 LSB of the PSD word ({lsb_db:.1f} dB)", ha="right", va="bottom", fontsize=4.3, color=BLUE,
                bbox=dict(fc="#ffffff", ec="none", pad=0.8, alpha=0.85), zorder=6)
        ax.set_xlim(0, fs / 2)
        ax.set_ylim(bot, top)
        ax.set_yticks(np.arange(bot, top + 1, 10.0 if top - bot <= 50 else 20.0))
        _title(ax, f"power (dB);  SQNR against the float Welch average: {e['sqnr']:.1f} dB")
        labels(ax, "frequency (Hz)", None, fs=4.8)
        if r == 0:
            legend(ax, "upper right", fs=4.4)
    D.save(FIGS / "fig22_psd_examples")
    with open(RESULTS / "model_psd_examples_wl10.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(PSD_EXAMPLE_HEADER)
        w.writerows(psd_example_rows(wl, win_name))


FIGURES = [fig21_window_responses, fig22_psd_examples]


def main():
    for fn in FIGURES:
        fn()
        print("wrote", fn.__name__)


if __name__ == "__main__":
    main()
