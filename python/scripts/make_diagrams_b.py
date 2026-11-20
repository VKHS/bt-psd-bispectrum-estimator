"""Figure sheets 4-8: segment merge, folded FD6 filter, CRU rotator, window/squarer/accumulator, bispectrum.

Each sheet combines a block diagram (panel a) with data panels computed by the Python model on synthetic signals.
Values reported in the paper are drawn as hollow black diamonds; nothing here is measured hardware."""
import matplotlib.pyplot as plt
import numpy as np

from _common import FIGS
from _figdata import bispec_pair, exact_kernel, fd_gain, read_csv, squarer_cycles
from _sheet import PAPER, card, cmap_blue, cmap_orange, labels, legend
from btpsd.cordic import default_rom
from btpsd.estimator import HwConfig
from btpsd.fdfilter import PAPER_FD6, PAPER_FD6_START, fd_frequency_response, fd_group_delay, fd_kernel
from btpsd.fixedpoint import naf
from btpsd.viz import INK, INK2, MUTED, PAL, Diagram
from btpsd.windows import spectral_kernel

BLUE, ORANGE, AQUA, VIOLET, YELLOW = PAL["blue"], PAL["orange"], PAL["aqua"], PAL["violet"], PAL["yellow"]
DASH = (0, (3.0, 1.6))


def csd_str(n):
    t = naf(int(n))
    return " ".join(f"{'+' if s > 0 else '-'}2^{e}" for s, e in reversed(t)) or "0"


def fig04_segment_merge():
    D = Diagram(132, 96, "Segment merge: exact even bins, interpolated odd bins")
    d = D.view(0, 50)
    D.panel(1.5, 89.6, "(a)", "Merge datapath")
    d.text(2, 37.6, "Chunk spectra A (previous) and B (new), 512 bins each, give the 1024-bin spectrum of the segment [chunk A, chunk B]:",
           fs=7.2, ha="left", color=INK2)
    d.text(2, 34.0, "X(2m) = [A(m) + B(m)] / 2  (exact)            X(2m+1) = D(m + 1/2) / 2  with D = A - B:  a half-bin value, interpolated from D(m)",
           fs=7.2, ha="left", color=INK)
    yt, yb = 24.0, 10.0
    d.text(5.0, yt, "A(m)", fs=7.8, color=INK, bold=True)
    d.text(5.0, yb, "B(m)", fs=7.8, color=INK, bold=True)
    d.box(12, 4, 17, 25, "Radix-2 combine\nS = A + B\nD = A - B\n(WL+1 bits)", "compute")
    d.arrow([(8.4, yt), (12, yt)], None)
    d.arrow([(8.4, yb), (12, yb)], None)
    d.box(34, yt - 4.5, 18, 9, "Even path\n(A + B) >> 1\nround, saturate", "compute")
    d.arrow([(29, yt), (34, yt)], None)
    d.text(31.5, yt + 1.9, "S", fs=6.8, color=INK)
    d.arrow([(52, yt), (114, yt)], None)
    d.text(83, yt + 2.0, "X(2m)", fs=7.2, color=INK, bold=True)
    d.box(34, yb - 4.5, 12, 9, "(-1)^m\nsign\nalternate", "free")
    d.arrow([(29, yb), (34, yb)], None)
    d.text(31.5, yb + 1.9, "D", fs=6.8, color=INK)
    d.box(50, yb - 6, 24, 12, "6-tap symmetric FIR\nh = {h0, h1, h2, h2, h1, h0}\ntaps m-2 .. m+3\nlook-ahead of 3 bins", "compute")
    d.arrow([(46, yb), (50, yb)], None)
    d.box(77, yb - 4.5, 13, 9, "x (-j)(-1)^m\nre/im swap\nand negate", "free")
    d.arrow([(74, yb), (77, yb)], None)
    d.box(93, yb - 4.5, 12, 9, "Scale\n>> (F + 1)\nround, sat", "compute")
    d.arrow([(90, yb), (93, yb)], None)
    d.box(114, 4, 16, 25, "Interleave\nX(2m), X(2m+1)\n1024 bins\nscale 1/1024", "control")
    d.arrow([(105, yb), (114, yb)], None)
    d.text(109.5, yb + 2.0, "X(2m+1)", fs=7.2, color=INK, bold=True)
    d.text(66, 1.8, "Magenta blocks cost no arithmetic (sign flips, re/im swap).  The 6-tap truncation of the half-bin kernel, not the word length, "
           "bounds odd-bin accuracy (model results: docs/PAPER_NOTES.md).", fs=6.8, color=INK2)
    # ---- (b) gain error across the chunk --------------------------------------------------------------------
    n0, G6, Gi = fd_gain(PAPER_FD6, PAPER_FD6_START)
    hs, ss = fd_kernel(6, "sinc")
    hl, sl = fd_kernel(6, "lagrange")
    _, Gs, _ = fd_gain(hs, ss)
    _, Gl, _ = fd_gain(hl, sl)
    err = lambda G: 20 * np.log10(np.abs(G - Gi) + 1e-9)  # noqa: E731
    rms = lambda G: 10 * np.log10(np.mean(np.abs(G - Gi) ** 2))  # noqa: E731
    ax = card(D, 2, 2, 41, 44, "(b)", "Odd-bin gain error", "interpolator vs ideal half-bin value, impulse at sample n0")
    ax.plot(n0, err(Gl), color=AQUA, lw=0.8, ls=DASH, label=f"Lagrange 6 tap, rms {rms(Gl):.1f} dB")
    ax.plot(n0, err(Gs), color=ORANGE, lw=0.8, ls=DASH, label=f"sinc 6 tap, rms {rms(Gs):.1f} dB")
    ax.plot(n0, err(G6), color=BLUE, lw=1.0, label=f"Eq. (6) taps, rms {rms(G6):.1f} dB")
    ax.set_xlim(0, 511)
    ax.set_ylim(-52, 13)
    ax.set_yticks([-50, -40, -30, -20, -10, 0])
    ax.set_xticks([0, 128, 256, 384, 512])
    labels(ax, "sample n0 within the chunk", "gain error (dB)")
    legend(ax, "upper center")
    # ---- (c) exact kernel vs the six taps ---------------------------------------------------------------------
    dd, hx = exact_kernel(-9, 10)
    kept = float((PAPER_FD6 ** 2).sum())
    ax = card(D, 45, 2, 41, 44, "(c)", "Kernel truncation", "ideal kernel h(d) = -(-1)^d / (pi (d - 1/2)) and the 6 taps kept")
    ax.axhline(0, color=MUTED, lw=0.6)
    ax.vlines(dd, 0, hx, color=MUTED, lw=0.9)
    ax.plot(dd, hx, "o", ms=2.2, color=MUTED, label="ideal kernel")
    ax.plot(np.arange(-2, 4), PAPER_FD6, "o", ms=4.4, color=BLUE, mec="#ffffff", mew=0.6, label="six taps (Eq. 6)")
    ax.set_ylim(-0.46, 0.8)
    ax.set_xlim(-9.8, 10.8)
    ax.set_xticks([-8, -4, 0, 4, 8])
    labels(ax, "tap offset d (bins from the output bin)", "h(d)")
    legend(ax, "upper right")
    ax.text(0.5, 0.075, f"the six taps keep {100 * kept:.1f} % of the kernel energy;\n"
            f"the {100 * (1 - kept):.1f} % tail ({10 * np.log10(1 - kept):.1f} dB) bounds the odd-bin SNR", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=5.3, color=INK)
    # ---- (d) SNR versus tap count: model, analytic bound, paper-reported -----------------------------------------
    rows = read_csv("fd_merge_snr_vs_taps.csv")
    sel = lambda kern, col: [(int(r["taps"]), float(r[col])) for r in rows if r["signal"] == "eeg_like" and r["kernel"] == kern]  # noqa: E731
    t_s, snr_s = zip(*sel("sinc", "snr_complex_spectrum_db"))
    _, psd_s = zip(*sel("sinc", "sqnr_psd_hann_db"))
    bound = [-10 * np.log10(1 - float((fd_kernel(L, "sinc")[0] ** 2).sum())) for L in t_s]
    import csv
    from _common import DATA
    with open(DATA / "paper_table3_fd_taps_snr.csv") as f:
        pr = [(int(r["fd_taps"]), float(r["snr_db"])) for r in csv.DictReader(f)]
    ax = card(D, 88, 2, 42, 44, "(d)", "SNR versus tap count", "odd-bin model results on a synthetic EEG-like signal; Table 3 hollow")
    ax.plot(t_s, bound, color=MUTED, lw=0.9, ls=DASH, label="limit, flat spectrum")
    ax.plot(t_s, snr_s, "-o", color=BLUE, lw=0.9, ms=3.0, label="model, odd bins")
    ax.plot(t_s, psd_s, "-s", color=AQUA, lw=0.9, ms=2.8, label="model, PSD (Hann)")
    ax.plot(*zip(*pr), label="paper, Table 3", **PAPER)
    ax.set_ylim(-18, 58)
    ax.set_yticks([0, 10, 20, 30, 40, 50])
    ax.set_xticks(range(3, 11))
    labels(ax, "taps in the interpolator", "SNR (dB)")
    legend(ax, "lower center", ncol=2, fs=5.0, columnspacing=0.9)
    D.save(FIGS / "fig04_segment_merge")


def fig05_fd6_filter():
    cfg = HwConfig(fd_frac=12)
    hq = [int(v) for v in cfg.fd_int()]
    assert hq[0] == hq[5] and hq[1] == hq[4] and hq[2] == hq[3], hq
    n_const = sum(len(naf(h)) - 1 for h in hq[:3])
    n_add = 3 + n_const + 2
    F = 1 << cfg.fd_frac
    D = Diagram(121, 106, "Folded symmetric 6-tap FD filter (direct form) with shift-add constants")
    d = D.view(0, 44)
    D.panel(1.5, 100.0, "(a)", "Filter structure")
    X = [14 + 16 * i for i in range(6)]
    ys = 49.0
    names = ["m+3", "m+2", "m+1", "m", "m-1", "m-2"]
    d.arrow([(2, ys), (X[0], ys)], None, head=False)
    d.text(8, ys + 2.1, "D'(k)", fs=7.2, color=INK, bold=True)
    for i, xx in enumerate(X):
        d.junction(xx, ys, r=0.55)
        d.text(xx, ys + 3.0, f"D'({names[i]})", fs=6.9, color=INK)
        if i < 5:
            d.box(xx + 4.5, ys - 2.2, 7.0, 4.4, "z^-1", "free", fs=6.8, round_=0.6)
            d.arrow([(xx, ys), (xx + 4.5, ys)], None, head=False)
            d.arrow([(xx + 11.5, ys), (X[i + 1], ys)], None, head=False)
    rows = {2: (41.0, 50.0), 1: (35.0, 58.0), 0: (29.0, 66.0)}      # pair -> (row y, pre-adder x)
    mult_y = 21.0
    for p, (yy, ax) in rows.items():
        xl, xr = X[p], X[5 - p]
        d.arrow([(xl, ys), (xl, yy), (ax - 1.15, yy)], None, head=False)
        d.arrow([(xr, ys), (xr, yy), (ax + 1.15, yy)], None, head=False)
        d.adder(ax, yy, "+")
        d.arrow([(ax, yy - 1.15), (ax, mult_y + 1.3)], None)
        d.adder(ax, mult_y, "x", r=1.3)
        d.text(ax + 1.7, mult_y + 3.4, f"h{p}", fs=7.0, color=INK, bold=True, ha="left")
        d.arrow([(ax, mult_y - 1.3), (ax, 15.7)], None)
    d.box(44, 8, 28, 7.6, "Adder tree\n2 adders", "compute")
    d.box(77, 8, 15, 7.6, "x (-j)(-1)^m\nre/im swap", "free")
    d.box(96, 8, 12, 7.6, "Scale\n>> (F + 1)", "compute")
    d.arrow([(72, 11.8), (77, 11.8)], None)
    d.arrow([(92, 11.8), (96, 11.8)], None)
    d.arrow([(108, 11.8), (113, 11.8)], None)
    d.text(114.6, 11.8, "O(m)", fs=7.2, color=INK, bold=True, ha="left")
    lines = ["Tap constants (fd_frac = 12), CSD digits",
             f"h2 = {hq[2]:+d}/{F} = {hq[2] / F:+.4f} = {csd_str(hq[2])}",
             f"h1 = {hq[1]:+d}/{F} = {hq[1] / F:+.4f} = {csd_str(hq[1])}",
             f"h0 = {hq[0]:+d}/{F} = {hq[0] / F:+.4f} = {csd_str(hq[0])}",
             f"constant multipliers: {len(naf(hq[2])) - 1} + {len(naf(hq[1])) - 1} + {len(naf(hq[0])) - 1} = {n_const} adders"]
    d.box(2, 6, 40, 15, "\n".join(lines), "rom", fs=6.6)
    d.text(60, 1.8, f"{n_add} adders in total with these constants (3 pre-adders, {n_const} in the constant multipliers, 2 in the tree).  "
           "O(m) is the odd-bin output; the sign alternation D' = (-1)^k D is the block shown in fig04.", fs=6.6, color=INK2)
    # ---- (b) magnitude response ------------------------------------------------------------------------------
    ax = card(D, 2, 2, 38, 40, "(b)", "Magnitude response", "frequency of the FIR = position within the chunk")
    ax.axhline(0, color=MUTED, lw=0.6, ls=DASH)
    for L, col, ls in ((4, ORANGE, DASH), (8, AQUA, DASH), (6, BLUE, "-")):
        h = PAPER_FD6 if L == 6 else fd_kernel(L, "sinc")[0]
        w, Hr = fd_frequency_response(h)
        ax.plot(w / np.pi, 20 * np.log10(np.abs(Hr) + 1e-9), color=col, ls=ls, lw=1.0 if L == 6 else 0.8,
                label=f"{L} taps" + (" (Eq. 6)" if L == 6 else ""))
    ax.set_xlim(0, 1)
    ax.set_ylim(-30, 5)
    labels(ax, "frequency (x pi rad / sample)", "|H| (dB)")
    legend(ax, "lower left")
    # ---- (c) group delay ----------------------------------------------------------------------------------------
    ax = card(D, 42, 2, 38, 40, "(c)", "Group delay", "even-length symmetric kernels are linear phase")
    for L, col, ls in ((5, ORANGE, DASH), (7, AQUA, DASH), (6, BLUE, "-")):
        h = PAPER_FD6 if L == 6 else fd_kernel(L, "sinc")[0]
        w, gd = fd_group_delay(h)
        ax.plot(w / np.pi, gd, color=col, ls=ls, lw=1.0 if L == 6 else 0.8, label=f"{L} taps" + (" (Eq. 6)" if L == 6 else ""))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 6)
    labels(ax, "frequency (x pi rad / sample)", "delay (samples)")
    legend(ax, "upper right")
    # ---- (d) adder count of the folded form ---------------------------------------------------------------------
    Ls = [4, 6, 8, 10]
    pre, cst, tree = [], [], []
    for L in Ls:
        h = PAPER_FD6 if L == 6 else fd_kernel(L, "sinc")[0]
        hq = [int(round(v * (1 << cfg.fd_frac))) for v in h[: L // 2]]
        pre.append(L // 2)
        cst.append(sum(len(naf(v)) - 1 for v in hq))
        tree.append(L // 2 - 1)
    assert pre[1] + cst[1] + tree[1] == n_add, (pre, cst, tree, n_add)
    ax = card(D, 82, 2, 37, 40, "(d)", "Adder count", "folded form, taps as CSD with fd_frac = 12")
    xs = np.arange(len(Ls))
    ax.bar(xs, pre, 0.58, color=AQUA, ec="#ffffff", lw=0.8, label="pre-adders")
    ax.bar(xs, cst, 0.58, bottom=pre, color=ORANGE, ec="#ffffff", lw=0.8, label="constant multipliers")
    ax.bar(xs, tree, 0.58, bottom=np.add(pre, cst), color=BLUE, ec="#ffffff", lw=0.8, label="adder tree")
    for i in range(len(Ls)):
        ax.text(i, pre[i] + cst[i] + tree[i] + 0.5, str(pre[i] + cst[i] + tree[i]), ha="center", fontsize=6.0, color=INK)
    ax.set_xticks(xs)
    ax.set_xticklabels([str(L) for L in Ls])
    ax.set_ylim(0, 30)
    ax.set_xlim(-0.6, len(Ls) - 0.4)
    labels(ax, "taps (6 = Eq. 6 constants)", "adders")
    legend(ax, "upper left")
    D.save(FIGS / "fig05_fd6_filter")


def fig06_cru():
    rom = default_rom()
    ae = np.abs(rom.angle_error()[1:])
    D = Diagram(126, 106, "Central rotational unit (CRU): constant rotator, shift-and-add only")
    d = D.view(0, 42)
    D.panel(1.5, 100.4, "(a)", "Rotator datapath")
    d.text(6.5, 40.0, "q (10 bit)", fs=7.0, color=INK, bold=True)
    d.arrow([(11.4, 40), (16, 40)], None)
    d.box(16, 33, 40, 14, "Octant reduction\nquad = q >> 8,  r = q & 255\nif r > 128: quad += 1, r = 256 - r\ndirection = ccw if folded, else cw", "control")
    d.box(62, 33, 36, 14, "Coefficient ROM, r = 0 .. 128\n(a, b, inv) per entry, a >= b\ninv ~ 2^F / |a + j b|\nconstants of few CSD digits", "rom")
    d.arrow([(56, 40), (62, 40)], None)
    d.text(59, 41.8, "r", fs=7.0, color=INK, bold=True)
    d.text(6.0, 14.0, "x = xr + j xi", fs=6.8, color=INK, bold=True)
    d.arrow([(10.6, 14), (16, 14)], None)
    d.box(16, 8, 32, 12, "Complex rotation\nyr = xr a +/- xi b\nyi = xi a -/+ xr b\nshift-add only, exact", "compute")
    d.box(52, 9, 22, 10, "Gain correction\nz = y * inv / 2^F\nshift-add", "compute")
    d.box(78, 9, 16, 10, "Round, saturate\nhalf up\nWL bits", "compute")
    d.box(98, 9, 18, 10, "Quadrant (-j)^quad\nswap / negate\nfree", "free")
    for xa, xb in ((48, 52), (74, 78), (94, 98)):
        d.arrow([(xa, 14), (xb, 14)], None)
    d.arrow([(116, 14), (121.5, 14)], None)
    d.text(123.0, 14.0, "y", fs=7.6, color=INK, bold=True)
    d.arrow([(22, 33), (22, 20)], None, ls="--")
    d.text(23.3, 26.5, "direction", fs=6.6, ha="left", color=INK)
    d.arrow([(66, 33), (66, 27), (42, 27), (42, 20)], None, ls="--")
    d.text(54, 28.5, "a, b", fs=6.6, color=INK)
    d.arrow([(70, 33), (70, 19)], None, ls="--")
    d.text(71.3, 24.0, "inv", fs=6.6, ha="left", color=INK)
    d.arrow([(36, 47), (36, 55), (107, 55), (107, 19)], None, ls="--")
    d.text(72, 56.6, "quad", fs=6.6, color=INK)
    d.text(63, 4.8, "y = (-j)^quad * x * (a -/+ j b) * inv / 2^F  ~  x * exp(-j q pi / 512)", fs=7.8, color=INK, bold=True)
    d.text(63, 1.6, f"Generated ROM for an angle tolerance of 0.01 rad: {len(ae) + 1} entries, max angle error {ae.max():.4f} rad, "
           f"rms {np.sqrt((ae ** 2).mean()):.4f} rad, gain error <= 2^-10, {rom.adders()[1:].mean():.1f} adders per rotation on average (model).",
           fs=6.6, color=INK2)
    r = np.arange(129)
    tol = 0.01
    w4, x0 = 29.2, 2.0
    step = 30.8
    # ---- (b) angle error per ROM entry ---------------------------------------------------------------------------
    ax = card(D, x0, 2, w4, 38, "(b)", "Angle error per entry", "generated ROM, r = 0 .. 128", fs=7.8, pad=(9.2, 7.4, 2.4, 9.2))
    ax.axhspan(-tol, tol, color=BLUE, alpha=0.08, lw=0)
    ax.axhline(tol, color=MUTED, lw=0.6, ls=DASH)
    ax.axhline(-tol, color=MUTED, lw=0.6, ls=DASH)
    ax.plot(r, rom.angle_error(), color=BLUE, lw=0.8)
    ax.set_xlim(0, 128)
    ax.set_ylim(-0.0125, 0.0125)
    ax.set_xticks([0, 32, 64, 96, 128])
    ax.set_yticks([-0.01, 0, 0.01])
    ax.text(127, tol + 0.0004, "tolerance 0.01 rad", ha="right", va="bottom", fontsize=5.0, color=INK2)
    labels(ax, "ROM index r", "angle error (rad)")
    # ---- (c) adders per entry ---------------------------------------------------------------------------------------
    ad = rom.adders()
    ax = card(D, x0 + step, 2, w4, 38, "(c)", "Adders per rotation", "complex rotation + gain correction", fs=7.8, pad=(9.2, 7.4, 2.4, 9.2))
    ax.step(r, ad, where="mid", color=VIOLET, lw=0.8)
    ax.axhline(ad[1:].mean(), color=MUTED, lw=0.6, ls=DASH)
    ax.text(3, 11.7, f"dashed: mean {ad[1:].mean():.1f}", ha="left", va="center", fontsize=5.2, color=INK2)
    ax.set_xlim(0, 128)
    ax.set_ylim(0, 12.5)
    ax.set_xticks([0, 32, 64, 96, 128])
    labels(ax, "ROM index r", "adders")
    # ---- (d), (e) tolerance trade-off ----------------------------------------------------------------------------------
    rows = read_csv("cru_rom_tradeoff.csv")
    tols = [float(q["angle_tol_rad"]) for q in rows]
    mean_ad = [float(q["mean_adders_per_rotation"]) for q in rows]
    s16 = [float(q["fft512_snr_wl16_db"]) for q in rows]
    s24 = [float(q["fft512_snr_wl24_db"]) for q in rows]
    k = tols.index(0.01)
    ax = card(D, x0 + 2 * step, 2, w4, 38, "(d)", "Cost of the tolerance", "mean adders per rotation", fs=7.8, pad=(9.2, 7.4, 2.4, 9.2))
    ax.plot(tols, mean_ad, "-o", color=VIOLET, lw=0.9, ms=2.8)
    ax.plot([tols[k]], [mean_ad[k]], "o", ms=7.5, mfc="none", mec=INK, mew=0.9)
    ax.set_xscale("log")
    ticks = sorted(tols)
    ax.set_xticks(ticks)
    ax.set_xticklabels([{0.001: "0.001", 0.005: "0.005", 0.01: "0.01", 0.02: "0.02"}.get(q, "") for q in ticks], fontsize=5.0)
    ax.set_xlim(0.0008, 0.025)
    ax.minorticks_off()
    ax.set_ylim(6, 11)
    ax.text(tols[k] * 0.8, mean_ad[k] - 0.7, "paper: 0.01 rad", ha="right", va="top", fontsize=5.2, color=INK)
    labels(ax, "angle tolerance (rad)", "adders")
    ax = card(D, x0 + 3 * step, 2, w4, 38, "(e)", "512-point FFT SNR", "model, random complex input", fs=7.8, pad=(9.2, 7.4, 2.4, 9.2))
    ax.plot(tols, s24, "-s", color=AQUA, lw=0.9, ms=2.6, label="WL = 24")
    ax.plot(tols, s16, "-o", color=BLUE, lw=0.9, ms=2.8, label="WL = 16")
    ax.plot([tols[k]], [s16[k]], "o", ms=7.5, mfc="none", mec=INK, mew=0.9)
    ax.set_xscale("log")
    ticks = sorted(tols)
    ax.set_xticks(ticks)
    ax.set_xticklabels([{0.001: "0.001", 0.005: "0.005", 0.01: "0.01", 0.02: "0.02"}.get(q, "") for q in ticks], fontsize=5.0)
    ax.set_xlim(0.0008, 0.025)
    ax.minorticks_off()
    ax.set_ylim(30, 66)
    legend(ax, "lower left")
    labels(ax, "angle tolerance (rad)", "SNR (dB)")
    D.save(FIGS / "fig06_cru_rotator")


def fig07_window_power():
    D = Diagram(118, 104, "Frequency-domain window filter, squarer and segment accumulator")
    d = D.view(0, 44)
    D.panel(1.5, 97.9, "(a)", "Power-estimation datapath")
    d.text(2, 51.4, "Window as a short convolution on the spectrum, Xw(k) = sum_d W(d) X(k-d); periodic raised-cosine windows are exactly sparse (1 / 3 / 3 / 5 taps):",
           fs=7.0, ha="left", color=INK2)
    kers = {n: spectral_kernel(n) for n in ("hann", "hamming", "blackman")}
    txt = "     ".join(f"{n}: [" + ", ".join(f"{v:+.2f}" for v in k.values()) + "]" for n, k in kers.items())
    d.text(2, 48.2, txt + "     rectangular: [1]", fs=7.0, ha="left", color=INK)
    ya = 38.0
    d.text(5.2, ya, "X(k)", fs=7.4, color=INK, bold=True)
    d.arrow([(8.4, ya), (11, ya)], None)
    d.box(11, ya - 5, 21, 10, "5-word shift register\nz^-1 chain, k-2 .. k+2\ncyclic over 1024 bins", "memory")
    d.box(37, ya - 5, 22, 10, "Coefficient bank\nW(0), W(+-1), W(+-2)\nCSD shift-add", "rom")
    d.box(64, ya - 5, 14, 10, "Adder tree\nexact, wide", "compute")
    d.box(83, ya - 5, 20, 10, ">> win_frac\nround, saturate\nto WL bits", "compute")
    for xa, xb in ((32, 37), (59, 64), (78, 83)):
        d.arrow([(xa, ya), (xb, ya)], None)
    yb = 24.0
    d.arrow([(103, ya), (108, ya), (108, yb), (104, yb)], None)
    d.text(109.4, 31.0, "Xw(k)", fs=7.2, color=INK, bold=True, ha="left")
    d.box(70, yb - 5, 34, 10, "Squarer, 2 units (re and im)\nn^2 = 1 + 3 + 5 + ... + (2|n| - 1)\n|n| cycles per operand, no multiplier", "compute")
    d.box(37, yb - 5, 25, 10, "P = re^2 + im^2\n2 WL bits, unsigned\none value per bin", "compute")
    d.arrow([(70, yb), (62, yb)], None)
    yc = 8.0
    d.arrow([(49.5, yb - 5), (49.5, yc + 5)], None)
    d.box(11, yc - 5, 50, 10, "Accumulator memory, 1024 words\nsegment 1: Save (P)      segments 2 .. S: Add (acc + P)\nbins 513..1023 mirror 511..1 for a real input", "accum")
    d.box(70, yc - 5, 34, 10, "P(k) = acc >> avg_shift\navg_shift = ceil(log2 S)\n1024 bins", "io")
    d.arrow([(61, yc), (70, yc)], None)
    d.text(59, 0.9, "The sum-of-odds squarer needs |n| cycles per operand, so the power stage is data- and WL-dependent (results/squarer_cycle_model.csv).",
           fs=6.8, color=INK2)
    w3 = 36.0
    # ---- (b) window kernels ----------------------------------------------------------------------------------------
    ax = card(D, 2, 2, w3, 40, "(b)", "Window kernels", "Xw(k) = sum_d W(d) X(k - d)")
    for off, (name, col) in zip((-0.27, -0.09, 0.09, 0.27), (("rect", MUTED), ("hann", BLUE), ("hamming", ORANGE), ("blackman", AQUA))):
        ker = spectral_kernel(name)
        dd = np.array(sorted(ker))
        vv = np.array([ker[i] for i in dd])
        ax.vlines(dd + off, 0, vv, color=col, lw=1.1)
        ax.plot(dd + off, vv, "o", ms=2.6, color=col, label=name)
    ax.axhline(0, color=MUTED, lw=0.6)
    ax.set_xticks([-2, -1, 0, 1, 2])
    ax.set_xlim(-2.7, 2.7)
    ax.set_ylim(-0.4, 1.12)
    labels(ax, "offset d (bins)", "W(d)")
    legend(ax, "upper right", ncol=2, fs=5.0, columnspacing=0.8)
    # ---- (c) distribution of squarer cycles ---------------------------------------------------------------------------
    ax = card(D, 41, 2, w3, 40, "(c)", "Squarer cycles per bin", "sum of odds, EEG-like input, model")
    for wl, col in ((16, VIOLET), (10, BLUE)):
        c = np.sort(squarer_cycles(wl, nchunks=40, seed=0))
        ax.plot(np.maximum(c, 1), np.arange(1, c.size + 1) / c.size, color=col, lw=1.0,
                label=f"WL = {wl}, mean {c.mean():.1f}" if c.mean() < 10 else f"WL = {wl}, mean {c.mean():.0f}")
    ax.set_xscale("log")
    ax.set_xlim(1, 2e4)
    ax.set_ylim(0, 1.04)
    labels(ax, "cycles for one bin (max of |re|, |im|)", "fraction of bins")
    legend(ax, "lower right")
    # ---- (d) squarer pairs needed ---------------------------------------------------------------------------------------
    rows = read_csv("squarer_cycle_model.csv")
    sigs = ["white", "pink", "eeg_like", "ar2_narrowband"]
    val = {(q["signal"], int(q["wl"])): float(q["serial_squarer_units_needed_for_R0.88"]) for q in rows}
    ax = card(D, 80, 2, w3, 40, "(d)", "Squarer pairs needed", "to sustain 0.88 samples / cycle (model)")
    xs = np.arange(len(sigs))
    ax.bar(xs - 0.2, [val[(q, 10)] for q in sigs], 0.36, color=BLUE, label="WL = 10")
    ax.bar(xs + 0.2, [val[(q, 16)] for q in sigs], 0.36, color=VIOLET, label="WL = 16")
    ax.axhline(1, color=INK2, lw=0.7, ls=DASH, label="one pair")
    ax.set_yscale("log")
    ax.set_ylim(0.6, 400)
    ax.set_yticks([1, 10, 100])
    ax.set_yticklabels(["1", "10", "100"])
    ax.minorticks_off()
    ax.set_xticks(xs)
    ax.set_xticklabels(["white", "pink", "EEG", "AR(2)"], fontsize=5.4)
    ax.set_xlim(-0.6, len(sigs) - 0.4)
    labels(ax, None, "pairs (x)")
    legend(ax, "upper left")
    D.save(FIGS / "fig07_window_squarer_accumulator")


def fig08_bispectrum():
    D = Diagram(124, 100, "Bispectrum datapath: direct FFT triple product")
    d = D.view(0, 44)
    D.panel(1.5, 93.6, "(a)", "Triple-product datapath")
    d.text(2, 45.8, "B(f1, f2) = (1/S) sum_s F_s(f1) F_s(f2) F_s*(f1 + f2),   f1, f2 = 0 .. N/4 - 1  (so f1 + f2 < N/2),   N = 512, hop 256",
           fs=7.2, ha="left", color=INK2)
    yc = 31.0
    d.box(2, yc - 6, 20, 12, "Index generator\nf1, f2 loops\nf3 = f1 + f2", "control")
    d.box(28, yc - 11, 26, 22, "Spectrum memory\nF_s(f), windowed\n512 words x 2 WL\n3 read ports", "memory")
    d.arrow([(22, yc), (28, yc)], None)
    d.box(62, yc - 6, 24, 12, "Complex mult 1\nP12 = F(f1) F(f2)\nround >> (WL - 1)", "compute")
    d.box(92, yc - 6, 24, 12, "Complex mult 2\nP12 x conj(F(f3))\nexact, no rounding", "compute")
    d.arrow([(54, yc + 3.5), (62, yc + 3.5)], None)
    d.text(58, yc + 5.2, "F(f1)", fs=6.8, color=INK)
    d.arrow([(54, yc - 3.5), (62, yc - 3.5)], None)
    d.text(58, yc - 5.2, "F(f2)", fs=6.8, color=INK)
    d.arrow([(54, 22.5), (98, 22.5), (98, yc - 6)], None)
    d.text(76, 20.8, "F(f3)", fs=6.8, color=INK)
    d.arrow([(86, yc), (92, yc)], None)
    d.text(89, yc + 1.7, "P12", fs=6.6, color=INK)
    d.box(92, 4, 24, 12, "Accumulator\nover S segments\nwide register, exact", "accum")
    d.arrow([(110, yc - 6), (110, 16)], None)
    d.text(111.4, 20.8, "F1 F2 F3*", fs=6.6, color=INK, ha="left")
    d.box(56, 4, 28, 12, "B(f1, f2) memory\n>> avg_shift\n128 x 128 complex words", "io")
    d.arrow([(92, 10), (84, 10)], None)
    d.box(20, 4, 30, 12, "Bicoherence (optional)\n|B|^2 / (P(f1) P(f2) P(f3))\nreuses the PSD path", "compute", dashed=True)
    d.arrow([(56, 10), (50, 10)], None, ls="--")
    d.text(62, 1.2, "Two complex multipliers per (f1, f2) pair: P12 is rounded back to WL bits, the second product is exact and accumulated in a wide register.  "
           "Output scale: btpsd.bispectrum.bispectrum_hw_to_float.", fs=6.6, color=INK2)
    ref, hw = bispec_pair("qpc_coupled", wl=14)
    peak = np.abs(ref).max()
    # ---- (b) SQNR versus word length -------------------------------------------------------------------------------------
    rows = read_csv("bispectrum_sqnr_vs_wordlength.csv")
    ax = card(D, 2, 2, 38, 40, "(b)", "Accuracy versus word length", "magnitude SQNR against the float direct estimator")
    for name, col, mk, lab in (("pink", AQUA, "s", "pink noise"), ("qpc_coupled", BLUE, "o", "QPC, coupled"),
                               ("qpc_uncoupled", ORANGE, "^", "QPC, uncoupled")):
        pts = [(int(q["wl"]), float(q["sqnr_magnitude_db"])) for q in rows if q["signal"] == name]
        ax.plot(*zip(*pts), "-" + mk, color=col, lw=0.9, ms=3.0, label=lab)
    ax.plot([10, 18], [22.4, 22.4 + 6.02 * 8], color=MUTED, lw=0.7, ls=DASH, label="6.02 dB / bit")
    ax.set_xticks([10, 12, 14, 16, 18])
    ax.set_xlim(9.5, 18.5)
    ax.set_ylim(8, 82)
    labels(ax, "word length WL (bits)", "SQNR (dB)")
    legend(ax, "upper left")
    # ---- (c), (d) magnitude map and error map ---------------------------------------------------------------------------------
    mag_db = 20 * np.log10(np.abs(hw) / peak + 1e-12)
    err_db = 20 * np.log10(np.abs(hw - ref) / peak + 1e-12)
    for x0, tag, title, sub, data, cm, lo, hi in (
            (43, "(c)", "|B(f1, f2)|, WL = 14", "coupled QPC input, dB re peak; rings: (32, 104), (104, 32)", mag_db, cmap_blue(), -50.0, 0.0),
            (84, "(d)", "Model minus float", "|B_hw - B_ref|, dB re peak", err_db, cmap_orange(), -90.0, -30.0)):
        card(D, x0, 2, 38, 40, tag, title, sub, plot=False)
        ax = D.inset(x0 + 9.0, 8.8, 20.0, 20.0, grid=False)
        im = ax.imshow(np.clip(data, lo, hi).T, origin="lower", extent=[0, 128, 0, 128], cmap=cm, vmin=lo, vmax=hi,
                       interpolation="nearest", aspect="equal")
        ax.set_xticks([0, 32, 64, 96, 128])
        ax.set_yticks([0, 32, 64, 96, 128])
        labels(ax, "f1 (bin)", "f2 (bin)")
        cax = D.inset(x0 + 31.5, 8.8, 1.6, 20.0, grid=False)
        cb = plt.colorbar(im, cax=cax)
        cb.outline.set_linewidth(0.4)
        cb.ax.tick_params(labelsize=5.0, length=1.5, pad=1.0, width=0.5)
        if x0 == 43:
            ax.plot([32, 104], [104, 32], "o", ms=5.5, mfc="none", mec=INK, mew=0.7)
    D.save(FIGS / "fig08_bispectrum_datapath")


if __name__ == "__main__":
    fig04_segment_merge()
    fig05_fd6_filter()
    fig06_cru()
    fig07_window_power()
    fig08_bispectrum()
    print("diagrams 4-8 written")
