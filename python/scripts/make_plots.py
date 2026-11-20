"""Figure sheets 12-14: word-length / error budget, 2-D bispectrum atlas, bispectrum views.

All data come from the Python model on synthetic signals (``results/*.csv`` and live model calls); nothing here is
measured hardware and nothing is taken from the paper's recordings.  The paper's own bispectrum figures
(Figs. 24-32) are kept unchanged in docs/paper_figures."""
import matplotlib.pyplot as plt
import numpy as np

from _common import FIGS
from _figdata import bispec_pair, read_csv
from _sheet import card, cmap_blue, cmap_div, cmap_orange, labels, legend
from btpsd import bispectrum as B
from btpsd import signals as sg
from btpsd.metrics import sqnr_db
from btpsd.viz import INK, INK2, MUTED, PAL, Diagram

BLUE, ORANGE, AQUA, VIOLET, YELLOW = PAL["blue"], PAL["orange"], PAL["aqua"], PAL["violet"], PAL["yellow"]
DASH = (0, (3.0, 1.6))
SIGCOL = {"white": BLUE, "pink": ORANGE, "eeg_like": AQUA, "ar2_narrowband": VIOLET}
SIGNAME = {"white": "white", "pink": "pink", "eeg_like": "EEG-like", "ar2_narrowband": "AR(2)"}
WINCOL = {"rect": MUTED, "hann": BLUE, "hamming": ORANGE, "blackman": AQUA}


def fig12_wordlength():
    rows = read_csv("sqnr_vs_wordlength.csv")

    def get(sig, tol):
        r = sorted((int(q["wl"]), float(q["sqnr_vs_exact_welch_db"]), float(q["sqnr_vs_float_fd6_db"]), float(q["sqnr_float_fd6_vs_exact_db"]))
                   for q in rows if q["signal"] == sig and abs(float(q["rom_angle_tol_rad"]) - tol) < 1e-12)
        return np.array(r)

    W, H = 130, 97
    D = Diagram(W, H, "Word-length study of the PSD datapath (model, synthetic signals)",
                subtitle="SQNR = 10 log10(sum P_ref^2 / sum (P_ref - P_hw)^2), no fitting.  Counterpart of the paper's SQNR-versus-word-length "
                         "plots (Figs. 21-23, which use recorded EEG, EOG and EMG).")
    # ---- (a) total SQNR and the 6-tap merge limit ---------------------------------------------------------------------
    ax = card(D, 2, 43.5, 62, 44, "(a)", "SQNR against the exact Welch PSD", "ROM tolerance 0.01 rad; dashed: limit of the 6-tap merge in floating point")
    for sig in SIGCOL:
        r = get(sig, 0.01)
        ax.axhline(r[0, 3], color=SIGCOL[sig], lw=0.7, ls=DASH)
        ax.plot(r[:, 0], r[:, 1], "-o", color=SIGCOL[sig], lw=1.0, ms=2.8, label=SIGNAME[sig])
    ax.set_xlim(9.5, 24.5)
    ax.set_ylim(10, 50)
    ax.set_xticks([10, 12, 14, 16, 18, 20, 22, 24])
    labels(ax, "word length WL (bits)", "SQNR (dB)")
    legend(ax, "upper left", ncol=2, fs=5.2, columnspacing=1.0)
    # ---- (b) arithmetic error alone ---------------------------------------------------------------------------------------
    ax = card(D, 66, 43.5, 62, 44, "(b)", "Arithmetic error alone", "SQNR against the floating-point 6-tap merge; dotted: ROM tolerance 0.0025 rad")
    for sig in SIGCOL:
        r = get(sig, 0.01)
        r2 = get(sig, 0.0025)
        ax.plot(r2[:, 0], r2[:, 2], color=SIGCOL[sig], lw=0.8, ls=":")
        ax.plot(r[:, 0], r[:, 2], "-o", color=SIGCOL[sig], lw=1.0, ms=2.8)
    r = get("white", 0.01)
    x0 = np.array([10.0, 16.0])
    ax.plot(x0, r[0, 2] + 6.02 * (x0 - 10), color=INK2, lw=0.7, ls=DASH)
    ax.text(15.3, 40.0, "6.02 dB / bit", fontsize=5.2, color=INK2, rotation=0, ha="left")
    ax.set_xlim(9.5, 24.5)
    ax.set_ylim(10, 90)
    ax.set_xticks([10, 12, 14, 16, 18, 20, 22, 24])
    labels(ax, "word length WL (bits)", "SQNR (dB)")
    # ---- (c) error budget ----------------------------------------------------------------------------------------------------
    ax = card(D, 2, 2, 62, 40, "(c)", "Error budget at WL = 12 and 16", "share of the error power, assuming uncorrelated errors")
    ys, lab = [], []
    k = 0
    for sig in SIGCOL:
        for wl in (12, 16):
            r = get(sig, 0.01)
            row = r[r[:, 0] == wl][0]
            p_fd, p_ar = 10 ** (-row[3] / 10), 10 ** (-row[2] / 10)
            share = p_fd / (p_fd + p_ar)
            y = -k
            ax.barh(y, share * 100, height=0.72, color=ORANGE, ec="#ffffff", lw=0.8, label="6-tap truncation" if k == 0 else None)
            ax.barh(y, (1 - share) * 100, left=share * 100, height=0.72, color=BLUE, ec="#ffffff", lw=0.8,
                    label="arithmetic (rounding, ROM)" if k == 0 else None)
            ax.text(101.5, y, f"{row[1]:.1f} dB", va="center", fontsize=5.2, color=INK)
            ys.append(y)
            lab.append(f"{SIGNAME[sig]}, {wl}")
            k += 1
        k += 0.0
    ax.set_yticks(ys)
    ax.set_yticklabels(lab, fontsize=5.0)
    ax.set_xlim(0, 100)
    ax.set_ylim(-7.6, 0.6)
    ax.set_xticks([0, 25, 50, 75, 100])
    labels(ax, "share of error power (%);  total SQNR at the right")
    ax.grid(axis="y", visible=False)
    ax.legend(fontsize=5.0, frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, handlelength=1.2, borderpad=0.1,
              columnspacing=1.0)
    # ---- (d) window choice at WL = 10 -------------------------------------------------------------------------------------
    rows = read_csv("model_sqnr_wl10_synthetic.csv")
    sigs = ["white", "pink", "eeg_like", "ar2_narrowband"]
    wins = ["rect", "hamming", "hann", "blackman"]
    val = {(q["signal"], q["window"]): float(q["sqnr_vs_exact_welch_db"]) for q in rows}
    ax = card(D, 66, 2, 62, 40, "(d)", "Window choice at WL = 10", "SQNR against the exact Welch PSD with the same window")
    xs = np.arange(len(sigs))
    for i, w in enumerate(wins):
        ax.bar(xs + (i - 1.5) * 0.2, [val[(s, w)] for s in sigs], 0.18, color=WINCOL[w], label={"rect": "rectangular"}.get(w, w))
    ax.set_xticks(xs)
    ax.set_xticklabels([SIGNAME[s] for s in sigs], fontsize=5.4)
    ax.set_ylim(0, 40)
    ax.set_xlim(-0.55, len(sigs) - 0.45)
    labels(ax, None, "SQNR (dB)")
    legend(ax, "upper left", ncol=4, fs=5.0, columnspacing=0.9, handlelength=1.0)
    D.save(FIGS / "fig12_wordlength_study")


def _map_axes(D, x, y, size, data, cmap, lo, hi, ticks_left=True, ticks_bottom=True, extent=(0, 128, 0, 128)):
    ax = D.inset(x, y, size, size, grid=False)
    im = ax.imshow(np.clip(data, lo, hi).T, origin="lower", extent=extent, cmap=cmap, vmin=lo, vmax=hi, interpolation="nearest", aspect="equal")
    ax.set_xticks([0, 64, 128] if extent[1] == 128 else [extent[0], 0, extent[1]])
    ax.set_yticks([0, 64, 128] if extent[3] == 128 else [extent[2], 0, extent[3]])
    if not ticks_left:
        ax.set_yticklabels([])
    if not ticks_bottom:
        ax.set_xticklabels([])
    return ax, im


def _cbar(D, im, x, y, h, ticks=None, labs=None, w=1.5):
    cax = D.inset(x, y, w, h, grid=False)
    cb = plt.colorbar(im, cax=cax)
    cb.outline.set_linewidth(0.4)
    cb.ax.tick_params(labelsize=5.0, length=1.5, pad=1.0, width=0.5)
    if ticks is not None:
        cb.set_ticks(ticks)
        if labs:
            cb.set_ticklabels(labs)
    return cb


def fig13_bispectrum_atlas():
    sigs = [("eeg_coupled", "EEG-like with a phase-locked harmonic", "alpha 10 Hz and its second harmonic on a 1/f background, fs = 256 Hz",
             [(20, 20)], (4, 36, 4, 36)),
            ("qpc_coupled", "Quadratic phase coupling", "tones at f1 = 32, f2 = 104 and f1 + f2 (bins), phase-locked",
             [(32, 104), (104, 32)], (16, 48, 88, 120)),
            ("pink", "Pink noise, no coupling", "Gaussian, so the true bispectrum is zero; the corner shows the 1/f energy",
             [], (0, 32, 0, 32))]
    wins = [("rect", "rectangular"), ("hann", "Hann"), ("hamming", "Hamming")]
    W, H = 130, 124
    D = Diagram(W, H, "2-D bispectrum atlas of the fixed-point model: signal class x window",
                subtitle="|B(f1, f2)| in dB relative to the peak of the floating-point direct estimator with the same window; 512-point segments, hop 256, WL = 14.  "
                         "Synthetic signals.")
    s = 24.0
    xs = [9.0 + k * 29.4 for k in range(4)]
    heads = ["overview, Hann window", "zoom, rectangular", "zoom, Hann", "zoom, Hamming"]
    for k in range(4):
        D.text(xs[k] + s / 2, H - 12.4, heads[k], fs=7.2, color=INK, bold=True)
    last_im = None
    rowh = 34.0
    for i, (sname, title, sub, marks, (x0, x1, y0, y1)) in enumerate(sigs):
        yb = H - 15.0 - (i + 1) * rowh + 4.0
        D.text(2.0, yb + s + 5.2, title, fs=7.4, ha="left", color=INK, bold=True)
        D.text(2.0, yb + s + 2.6, sub, fs=5.6, ha="left", color=INK2)
        data = {}
        for wkey, _ in wins:
            ref, hw = bispec_pair(sname, wl=14, window=wkey, nchunks=40)
            pk = np.abs(ref).max()
            data[wkey] = (20 * np.log10(np.abs(hw) / pk + 1e-12), sqnr_db(np.abs(ref), np.abs(hw)))
        ax = D.inset(xs[0], yb, s, s, grid=False)
        last_im = ax.imshow(np.clip(data["hann"][0], -60, 0).T, origin="lower", extent=(0, 128, 0, 128), cmap=cmap_blue(), vmin=-60, vmax=0,
                            interpolation="nearest", aspect="equal")
        ax.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0, fc="none", ec=INK, lw=0.7))
        ax.set_xticks([0, 64, 128])
        ax.set_yticks([0, 64, 128])
        ax.set_ylabel("f2 (bin)", fontsize=5.6, labelpad=1.0)
        if i == 2:
            ax.set_xlabel("f1 (bin)", fontsize=5.6, labelpad=1.0)
        for j, (wkey, _) in enumerate(wins):
            mag, sq = data[wkey]
            ax = D.inset(xs[j + 1], yb, s, s, grid=False)
            sub_img = np.clip(mag[x0:x1, y0:y1], -60, 0)
            ax.imshow(sub_img.T, origin="lower", extent=(x0 - 0.5, x1 - 0.5, y0 - 0.5, y1 - 0.5), cmap=cmap_blue(), vmin=-60, vmax=0,
                      interpolation="nearest", aspect="equal")
            ticks = lambda lo, hi: [lo, (lo + hi) // 2, hi - 1]  # noqa: E731
            ax.set_xticks(ticks(x0, x1))
            ax.set_yticks(ticks(y0, y1))
            for (a, b) in marks:
                if x0 <= a < x1 and y0 <= b < y1:
                    ax.plot([a], [b], "o", ms=7.0, mfc="none", mec=INK, mew=0.7)
            ax.text(0.97, 0.04, f"SQNR {sq:.0f} dB", transform=ax.transAxes, ha="right", fontsize=5.0, color=INK,
                    bbox=dict(fc="#ffffff", ec="none", pad=0.8, alpha=0.85))
            if i == 2:
                ax.set_xlabel("f1 (bin)", fontsize=5.6, labelpad=1.0)
    _cbar(D, last_im, xs[3] + s + 2.4, H - 15.0 - 3 * rowh + 4.0, 3 * rowh - 10.0, ticks=[-60, -50, -40, -30, -20, -10, 0])
    D.text(W - 1.8, H - 15.0 - 1.5 * rowh + 4.0, "dB re peak", fs=5.6, rot=90, color=INK2)
    D.text(2.0, 2.6, "Rings: coupling bins contained in the synthetic signals by construction.  SQNR: magnitude of the model against the floating-point "
           "direct estimator (same window).  The paper's maps (Figs. 24-32) are in docs/paper_figures.", fs=5.5, ha="left", color=INK2)
    D.save(FIGS / "fig13_bispectrum_atlas")


def fig14_bispectrum_views():
    sname = "eeg_coupled"
    nch = 40
    x = sg.eeg_coupled(512 * nch)
    maxlag = 32
    lags, c3 = B.third_order_cumulant(x, maxlag=maxlag, L=512, hop=256)
    Bi = B.bispectrum_indirect(x, maxlag=maxlag, L=512, hop=256, nfft=512, max_bin=128)
    ref, hw = bispec_pair(sname, wl=14, window="hann", nchunks=nch)
    W, H = 130, 100
    D = Diagram(W, H, "Cumulant, magnitude and phase views of the bispectrum (EEG-like signal, synthetic)",
                subtitle="Top: indirect estimate (third-order cumulant, Parzen lag window, 2-D DFT).  Bottom: direct FFT estimate, floating point and the "
                         "fixed-point model (Hann, WL = 14).")
    cw, cx = 41.0, (2.0, 44.5, 87.0)
    ch = 42.0
    ytop = H - 10.0
    # (a) cumulant
    y0 = ytop - ch
    card(D, cx[0], y0, cw, ch, "(a)", "Third-order cumulant", "c3(t1, t2), lags -32 .. 32 samples", plot=False)
    cmax = np.abs(c3).max()
    ax = D.inset(cx[0] + 8.6, y0 + 7.0, 25.5, 25.5, grid=False)
    im = ax.imshow(c3.T, origin="lower", extent=(-maxlag, maxlag, -maxlag, maxlag), cmap=cmap_div(), vmin=-cmax, vmax=cmax, aspect="equal", interpolation="nearest")
    labels(ax, "t1 (samples)", "t2 (samples)")
    ax.set_xticks([-32, 0, 32])
    ax.set_yticks([-32, 0, 32])
    _cbar(D, im, cx[0] + 35.4, y0 + 7.0, 25.5)
    # (b) magnitude indirect
    pk = np.abs(Bi).max()
    magi = 20 * np.log10(np.abs(Bi) / pk + 1e-12)
    card(D, cx[1], y0, cw, ch, "(b)", "Magnitude |B|, indirect", "dB re peak", plot=False)
    ax, im = _map_axes(D, cx[1] + 8.6, y0 + 7.0, 25.5, magi, cmap_blue(), -40.0, 0.0)
    labels(ax, "f1 (bin)", "f2 (bin)")
    ax.plot([20], [20], "o", ms=6.5, mfc="none", mec=INK, mew=0.7)
    _cbar(D, im, cx[1] + 35.4, y0 + 7.0, 25.5, ticks=[-40, -30, -20, -10, 0])
    # (c) phase indirect
    ph = np.ma.masked_where(magi < -25.0, np.angle(Bi))
    card(D, cx[2], y0, cw, ch, "(c)", "Phase of B, indirect", "shown where |B| is within 25 dB of the peak", plot=False)
    ax = D.inset(cx[2] + 8.6, y0 + 7.0, 25.5, 25.5, grid=False)
    ax.set_facecolor("#ececea")
    im = ax.imshow(ph.T, origin="lower", extent=(0, 128, 0, 128), cmap="twilight", vmin=-np.pi, vmax=np.pi, aspect="equal", interpolation="nearest")
    ax.set_xticks([0, 64, 128])
    ax.set_yticks([0, 64, 128])
    labels(ax, "f1 (bin)", "f2 (bin)")
    ax.text(0.97, 0.96, "phase = 0 near the peak:\nthe harmonic is locked, p2 = 2 p1", transform=ax.transAxes, ha="right", va="top", fontsize=4.8, color=INK,
            bbox=dict(fc="#ffffff", ec="none", pad=1.0, alpha=0.9))
    _cbar(D, im, cx[2] + 35.4, y0 + 7.0, 25.5, ticks=[-np.pi, 0, np.pi], labs=["-pi", "0", "pi"])
    # bottom row
    y1 = y0 - 1.6 - ch
    pkd = np.abs(ref).max()
    magr = 20 * np.log10(np.abs(ref) / pkd + 1e-12)
    magh = 20 * np.log10(np.abs(hw) / pkd + 1e-12)
    errh = 20 * np.log10(np.abs(hw - ref) / pkd + 1e-12)
    card(D, cx[0], y1, cw, ch, "(d)", "Magnitude |B|, direct, float", "dB re peak", plot=False)
    ax, im = _map_axes(D, cx[0] + 8.6, y1 + 7.0, 25.5, magr, cmap_blue(), -70.0, 0.0)
    labels(ax, "f1 (bin)", "f2 (bin)")
    ax.plot([20], [20], "o", ms=6.5, mfc="none", mec=INK, mew=0.7)
    _cbar(D, im, cx[0] + 35.4, y1 + 7.0, 25.5, ticks=[-70, -50, -30, -10])
    card(D, cx[1], y1, cw, ch, "(e)", "Magnitude |B|, direct, model", "WL = 14, dB re peak", plot=False)
    ax, im = _map_axes(D, cx[1] + 8.6, y1 + 7.0, 25.5, magh, cmap_blue(), -70.0, 0.0)
    labels(ax, "f1 (bin)", "f2 (bin)")
    ax.plot([20], [20], "o", ms=6.5, mfc="none", mec=INK, mew=0.7)
    _cbar(D, im, cx[1] + 35.4, y1 + 7.0, 25.5, ticks=[-70, -50, -30, -10])
    card(D, cx[2], y1, cw, ch, "(f)", "Model minus float", "|B_hw - B_ref|, dB re peak", plot=False)
    ax, im = _map_axes(D, cx[2] + 8.6, y1 + 7.0, 25.5, errh, cmap_orange(), -90.0, -30.0)
    labels(ax, "f1 (bin)", "f2 (bin)")
    _cbar(D, im, cx[2] + 35.4, y1 + 7.0, 25.5, ticks=[-90, -75, -60, -45, -30])
    D.text(2.0, 1.8, "Ring: the coupling bin (20, 20) contained in the synthetic signal by construction (alpha 10 Hz, 0.5 Hz per bin).  The indirect "
           "estimate is smoother because of the lag window.", fs=5.6, ha="left", color=INK2)
    D.save(FIGS / "fig14_bispectrum_views")


if __name__ == "__main__":
    fig12_wordlength()
    fig13_bispectrum_atlas()
    fig14_bispectrum_views()
    print("figures 12-14 written")
