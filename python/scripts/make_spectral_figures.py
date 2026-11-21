"""Figure sheet 29: what the whole estimator chain does to line spectra, windows and averaging (model, synthetic signals).

Every curve and number on the sheet is computed by ``btpsd.spectralstudy`` from the bit-accurate model (word length 16)
and from the exact Welch average of the same records; nothing is measured hardware and nothing is taken from the paper.
The long-format table behind the sheet is written to ``results/model_spectral_study.csv``; ``python/tests/test_spectral_study.py``
ties both to the model.
"""
from __future__ import annotations

import csv

import numpy as np

from _common import FIGS, RESULTS
from _sheet import card, labels, legend
from btpsd import spectralstudy as ss
from btpsd.viz import INK, INK2, MUTED, PAL, Diagram, tint

BLUE, ORANGE, AQUA = PAL["blue"], PAL["orange"], PAL["aqua"]
WINCOL = {"rect": MUTED, "hann": BLUE, "hamming": ORANGE, "blackman": AQUA}
WINNAME = {"rect": "rect", "hann": "Hann", "hamming": "Hamming", "blackman": "Blackman"}
DASH = (0, (2.5, 1.5))


def _title(ax, s):
    ax.set_title(s, fontsize=5.0, color=INK2, loc="left", pad=2.2)


def _lines(D, x, y, items, fs=4.6, pitch=2.25):
    for j, (t, bold) in enumerate(items):
        D.text(x, y - pitch * j, t, fs=fs, ha="left", color=INK if bold else INK2, bold=bold)


def _signed(v: float) -> str:
    return f"{v:.1f}".replace("-", "−")


def _signed2(v: float) -> str:
    return f"{v:.2f}".replace("-", "\u2212")


def _tone_panel(D, res):
    t = res["tone"]
    f, win = t["f"], t["window"]
    card(D, 2, 66, 146, 46, "(a)", "A tone between two bins",
         f"tone at bin {f:g} of 1 024, amplitude {ss.TONE_AMP:g} of full scale, {win.capitalize()} window, eight segments; "
         f"dB relative to the peak of the exact average", plot=False, fs=7.4)
    b = t["bins"]
    # ---- zoom -------------------------------------------------------------------------------------------------------
    ax = D.inset(10.0, 72.0, 56.0, 29.5, fs=4.6)
    lo, hi = 276, 326
    sel = (b >= lo) & (b <= hi)
    ax.axvspan(f - ss.EXCLUDE, f + ss.EXCLUDE, color="#efeee9", lw=0, zorder=0)
    ax.text(f, -97, f"within {ss.EXCLUDE} bins: not part of the floor", ha="center", va="bottom", fontsize=4.2, color=INK2)
    ax.plot(b[sel], t["ref"][sel], color=INK, lw=0.9, label="exact Welch average", zorder=3)
    ax.plot(b[sel], t["fd"][sel], color=tint(ORANGE, 0.35), lw=2.2, label="float six-tap merge", zorder=2)
    ev = sel & (b % 2 == 0)
    od = sel & (b % 2 == 1)
    ax.plot(b[sel], t["hw"][sel], color=BLUE, lw=0.5, zorder=4)
    ax.plot(b[ev], t["hw"][ev], "o", color=BLUE, ms=2.2, zorder=5, label=f"hardware model, WL = {ss.WL}: even bins")
    ax.plot(b[od], t["hw"][od], "o", mfc="#ffffff", mec=BLUE, mew=0.8, ms=2.2, zorder=5, label="odd bins")
    ax.set_xlim(lo, hi)
    ax.set_ylim(-100, 3)
    ax.set_yticks(range(-100, 1, 20))
    ax.set_xticks(range(280, 326, 10))
    _title(ax, "bins 276 to 326")
    labels(ax, "bin k of the 1 024-point segment", "dB re exact peak", fs=4.8)
    legend(ax, "upper right", fs=4.0)
    # ---- full band ---------------------------------------------------------------------------------------------------
    ax = D.inset(82.0, 72.0, 62.0, 29.5, fs=4.6)
    ax.plot(b, np.where(t["ref"] <= ss.FLOOR_DB + 1e-9, np.nan, t["ref"]), color=INK, lw=0.8, label="exact Welch average")
    ax.plot(b, t["hw"], color=BLUE, lw=0.8, label=f"hardware model, WL = {ss.WL}")
    for key, col in (("hw", BLUE), ("ref", INK)):
        ax.axhline(t["floor"][key], color=col, lw=0.6, ls=DASH)
    ax.text(508, t["floor"]["hw"] + 2.0, f"highest bin more than {ss.EXCLUDE} bins from the tone: hardware {_signed(t['floor']['hw'])} dB",
            ha="right", va="bottom", fontsize=4.2, color=BLUE)
    ax.text(508, t["floor"]["ref"] + 2.0, f"exact {_signed(t['floor']['ref'])} dB", ha="right", va="bottom", fontsize=4.2, color=INK)
    ax.set_xlim(0, 512)
    ax.set_ylim(-100, 3)
    ax.set_yticks(range(-100, 1, 20))
    ax.set_xticks([0, 128, 256, 384, 512])
    _title(ax, "whole one-sided band")
    labels(ax, "bin k", None, fs=4.8)
    ax.text(14, -96, f"float six-tap merge: {_signed(t['floor']['fd'])} dB (hidden under the hardware curve)", fontsize=4.2, color=ORANGE, va="bottom")
    ax.text(338, -96, f"exact curve: below {ss.FLOOR_DB:.0f} dB away from the tone", fontsize=4.2, color=INK2, va="bottom")


def _floor_panel(D, res):
    card(D, 2, 2, 47, 62, "(b)", "Leakage floor against tone position", "highest bin more than 12 bins away; solid: hardware, dashed: exact", plot=False, fs=7.0)
    ax = D.inset(10.8, 28.0, 35.2, 26.0, fs=4.4)
    off = np.array(ss.OFFSETS)
    for w in ss.WINDOWS:
        fl = res["floor"][w]
        ref = np.where(fl["ref"] <= ss.FLOOR_DB + 1e-9, np.nan, fl["ref"])
        ax.plot(off, fl["hw"], color=WINCOL[w], lw=1.0, marker="o", ms=1.8)
        ax.plot(off, ref, color=WINCOL[w], lw=0.8, ls=DASH)
    ax.set_xlim(0, 2)
    ax.set_ylim(-100, 0)
    ax.set_xticks([0, 0.5, 1, 1.5, 2])
    ax.set_yticks(range(-100, 1, 20))
    labels(ax, "tone at 300 + offset (bins)", "dB re exact peak", fs=4.6)
    handles = [ax.plot([], [], color=WINCOL[w], lw=1.2, label=WINNAME[w])[0] for w in ss.WINDOWS]
    legend(ax, "upper center", fs=3.8, ncol=4, columnspacing=0.8, handles=handles)
    fl = {w: res["floor"][w] for w in ss.WINDOWS}
    i_half, i_odd = ss.OFFSETS.index(0.5), ss.OFFSETS.index(1.0)
    fd_dev = max(float(np.max(np.abs(fl[w]["fd"] - fl[w]["hw"])[fl[w]["fd"] > ss.FLOOR_DB])) for w in ss.WINDOWS)
    on_even = [fl[w]["hw"][0] for w in ss.WINDOWS]
    on_odd = [fl[w]["hw"][i_odd] for w in ss.WINDOWS]
    half_hw = [fl[w]["hw"][i_half] for w in ss.WINDOWS]
    half_ref = [fl[w]["ref"][i_half] for w in ss.WINDOWS]
    _lines(D, 4.0, 22.0, [
        ("Tone between two bins (offset 0.5), hardware floor:", True),
        (" / ".join(_signed(v) for v in half_hw) + " dB (rect / Hann / Hamming / Blackman)", False),
        ("exact average: " + " / ".join(_signed(v) for v in half_ref) + " dB", False),
        (f"On an even bin: {_signed(max(on_even))} to {_signed(min(on_even))} dB, set by the rotator ROM;", False),
        (f"on an odd bin: {_signed(max(on_odd))} to {_signed(min(on_odd))} dB, set by the merge.", False),
        (f"Float six-tap merge: within {fd_dev:.1f} dB of the hardware floor", False),
        ("except on even bins, where it has no floor.", False),
    ], fs=5.1, pitch=2.5)


def _valley_panel(D, res):
    card(D, 51, 2, 47, 62, "(c)", "Resolving two tones", "valley depth between the peaks", plot=False, fs=7.0)
    seps = np.array(ss.SEPARATIONS)
    for j, r in enumerate(ss.VALLEY_RATIOS_DB):
        ax = D.inset(59.0 + 20.5 * j, 28.0, 17.0, 26.0, fs=4.4)
        for w in ss.WINDOWS:
            v = res["valley"][r][w]
            ax.plot(seps, v["hw"], color=WINCOL[w], lw=1.0, marker="o", ms=1.6)
            ax.plot(seps, v["ref"], color=WINCOL[w], lw=0.8, ls=DASH)
        ax.axhline(3, color=INK2, lw=0.5, ls=(0, (1, 1.5)))
        ax.set_xlim(2, 10)
        ax.set_ylim(0, 70)
        ax.set_xticks([2, 4, 6, 8, 10])
        _title(ax, "equal amplitudes" if r == 0 else f"second tone {abs(int(r))} dB lower")
        labels(ax, "separation (bins)", "valley depth (dB)" if j == 0 else None, fs=4.6)
        if j == 1:
            ax.set_yticklabels([])
        if j == 0:
            ax.text(10, 4.0, "3 dB", fontsize=3.8, color=INK2, ha="right", va="bottom")
    i8 = ss.SEPARATIONS.index(8)
    eq, lo20 = res["valley"][0.0], res["valley"][-20.0]
    hw8 = [eq[w]["hw"][i8] for w in ss.WINDOWS]
    ref8 = [eq[w]["ref"][i8] for w in ss.WINDOWS]
    best_hw0 = max(float(np.nanmax(eq[w]["hw"])) for w in ss.WINDOWS)
    best_ref0 = max(float(np.nanmax(eq[w]["ref"])) for w in ss.WINDOWS)
    best_hw20 = max(float(np.nanmax(lo20[w]["hw"])) for w in ss.WINDOWS)
    best_ref20 = max(float(np.nanmax(lo20[w]["ref"])) for w in ss.WINDOWS)
    _lines(D, 53.0, 22.0, [
        ("Equal tones 8 bins apart, valley depth:", True),
        ("hardware " + " / ".join(f"{v:.0f}" for v in hw8) + " dB (rect / Hann / Hamming / Blackman)", False),
        ("exact " + " / ".join(f"{v:.0f}" for v in ref8) + " dB", False),
        ("Deepest valley over 2 to 10 bins:", False),
        (f"equal tones: hardware {best_hw0:.0f} dB, exact {best_ref0:.0f} dB;", False),
        (f"20 dB apart: hardware {best_hw20:.0f} dB, exact {best_ref20:.0f} dB.", False),
        ("Tones at 200.3 + separation, not on bin centres.", False),
    ], fs=5.1, pitch=2.5)


def _averaging_panel(D, res):
    a = res["averaging"]
    card(D, 100, 2, 47, 62, "(d)", "Averaging", "white noise, Hann: spread of the PSD estimate", plot=False, fs=7.0)
    ax = D.inset(108.8, 28.0, 35.2, 26.0, fs=4.4)
    S = a["S"]
    ax.plot(S, a["inv_sqrt"], color=MUTED, lw=0.9, ls=DASH, label="1 / sqrt(S)")
    ax.plot(S, a["analytic"], color=INK, lw=0.8, label="analytic, 50 % overlap")
    ax.plot(S, a["ref"], "o", mfc="#ffffff", mec=INK, mew=0.8, ms=3.2, label="exact average")
    ax.plot(S, a["hw"], "o", color=BLUE, ms=1.8, label=f"hardware model, WL = {ss.WL}")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xticks(list(S))
    ax.set_xticklabels([str(int(s)) for s in S])
    ax.set_xlim(0.8, 80)
    ax.set_ylim(0.1, 1.25)
    ax.set_yticks([0.1, 0.2, 0.5, 1.0])
    ax.set_yticklabels(["0.1", "0.2", "0.5", "1"])
    ax.minorticks_off()
    labels(ax, "segments S", "std / mean of the estimate", fs=4.6)
    legend(ax, "lower left", fs=3.8)
    rho2 = {w: ss.welch_overlap_rho(w) ** 2 for w in ss.WINDOWS}
    i64 = list(S).index(64)
    _lines(D, 102.0, 22.0, [
        (f"S = 64: hardware {a['hw'][i64]:.3f}, exact {a['ref'][i64]:.3f}, analytic {a['analytic'][i64]:.3f}", True),
        (f"{ss.AVG_NREAL} seeded records per point; rms over bins {ss.AVG_BAND[0]} to {ss.AVG_BAND[1]}.", False),
        ("Neighbouring segments overlap by half, so the", False),
        ("variance is (1 + 2 (1 - 1/S) rho^2) / S with rho^2 =", False),
        (" / ".join(f"{rho2[w]:.3f}" for w in ss.WINDOWS) + " (rect / Hann / Hamming / Blackman).", False),
        ("Mean level of the hardware estimate against the", False),
        (f"exact average: {_signed2(a['bias_db'].min())} to {_signed2(a['bias_db'].max())} dB for S = 1 to 64.", False),
    ], fs=5.1, pitch=2.5)


def fig29_spectral_study(res=None):
    res = res or ss.run_study()
    D = Diagram(150, 124, "What the whole chain does to line spectra, windows and averaging (model, synthetic signals)",
                subtitle=f"Bit-accurate model at WL = {ss.WL} against the exact Welch average of the same records. The odd bins of every segment come from the "
                         f"six-tap merge (docs/PAPER_NOTES.md section 2.2); this sheet shows its effect on spectra.")
    _tone_panel(D, res)
    _floor_panel(D, res)
    _valley_panel(D, res)
    _averaging_panel(D, res)
    D.save(FIGS / "fig29_spectral_study")
    with open(RESULTS / "model_spectral_study.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(ss.CSV_HEADER)
        w.writerows(ss.study_rows(res))
    return res


FIGURES = [fig29_spectral_study]


def main():
    for fn in FIGURES:
        fn()
        print("wrote", fn.__name__)


if __name__ == "__main__":
    main()
