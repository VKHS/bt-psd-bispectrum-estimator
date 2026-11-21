"""Figure sheets 24-28: the 512-point memory-based FFT in detail (index space, memory contents, twiddles, accuracy).

Everything is computed from the models in this repository (``btpsd.fftstudy``, ``btpsd.schedule``, ``btpsd.r8fft``,
``btpsd.cordic``) on synthetic signals; nothing is measured hardware and nothing is taken from the paper.  The numbers
printed on the sheets are tied to the model by ``python/tests/test_fft_study.py`` and ``python/tests/test_fft_sheets.py``.
The RTL-level block diagram of the FFT core (fig26) is drawn from the Verilog sources via ``btpsd.rtlparse``.
"""
from __future__ import annotations

import csv
import re

import numpy as np
from matplotlib.patches import Circle, FancyBboxPatch, Polygon, Rectangle

from _common import FIGS, RESULTS
from _sheet import card, cmap_blue, cmap_orange, labels, legend
from btpsd import fftstudy as fs
from btpsd.cordic import PAPER_TOL, default_rom
from btpsd.rtlparse import RTL_DIR, define_values, parse_all
from btpsd.schedule import N, bank_of, output_index, pass_addresses, pass_twiddle_q, split
from btpsd.viz import GRID, INK, INK2, MUTED, PAL, SEQ_BLUE, SERIES, Diagram, tint

BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED = (PAL[k] for k in ("blue", "orange", "aqua", "yellow", "magenta", "green", "violet", "red"))
WL_SHEET = 12  # word length of the data sheets


def _title(ax, s):
    ax.set_title(s, fontsize=5.0, color=INK2, loc="left", pad=2.2)


def _db(v, floor=0.1):
    return 20 * np.log10(np.maximum(np.abs(v), floor))


def _gradient_bar(D, x, y, w, h, cmap, vmin, vmax, ticks, label, fs=4.4):
    """Small horizontal colour bar placed at diagram coordinates."""
    ax = D.inset(x, y, w, h, fs=fs, grid=False)
    ax.imshow(np.linspace(0, 1, 256)[None, :], cmap=cmap, aspect="auto", extent=(vmin, vmax, 0, 1))
    ax.set_yticks([])
    ax.set_xticks(ticks)
    ax.tick_params(axis="x", labelsize=fs, length=1.5, pad=1.0)
    ax.set_xlabel(label, fontsize=fs, labelpad=1.0)
    return ax


# ---------------------------------------------------------------------------------------------------------------
# fig25: memory contents after load and after every pass
# ---------------------------------------------------------------------------------------------------------------
LINE_BINS = ((37, ORANGE, "-"), (512 - 37, ORANGE, (0, (2.2, 1.6))), (150, VIOLET, "-"), (512 - 150, VIOLET, (0, (2.2, 1.6))))
PASS_TITLES = ["after load", "after pass 0", "after pass 1", "after pass 2"]
PASS_DIGITS = ["(a2, a1, a0) = (n2, n1, n0)", "(a2, a1, a0) = (k0, n1, n0)", "(a2, a1, a0) = (k0, k1, n0)", "(a2, a1, a0) = (k0, k1, k2)"]
PASS_WORDS = [None, 64, 8, 1]  # words of one spectral line after each pass


def line_cells(k: int, stage: int):
    """Cells ``(row a2, column 8 a1 + a0, width)`` that hold the spectral line ``k`` after ``stage`` passes (1..3), in the 8 x 64 layout."""
    k0, k1, k2 = k & 7, (k >> 3) & 7, (k >> 6) & 7
    if stage == 1:
        return k0, 0, 64
    if stage == 2:
        return k0, 8 * k1, 8
    return k0, 8 * k1 + k2, 1


def fig25_pass_snapshots():
    wl = WL_SHEET
    rom = default_rom()
    re, im = fs.two_tone_record(wl, seed=1)
    x = re + 1j * im
    fx = fs.pass_snapshots_fixed(re, im, rom, wl)
    ex = fs.pass_snapshots_float(x)
    sq = fs.fft_sqnr(re, im, wl, rom)
    W, H = 150, 124
    D = Diagram(W, H, f"Memory contents after load and after each pass (model, word length {wl})",
                subtitle=f"Two tones at bins 37 and 150 (amplitudes 1 and 1/2) in white noise, SNR 40 dB, peak 0.9 of full scale; the 512 words are laid out as\n"
                         f"a2 (rows) by 8 a1 + a0 (columns). Every pass resolves one digit of the frequency index, so a spectral line shrinks from 512 words to 64, 8 and finally 1.")
    cm_b, cm_o = cmap_blue(), cmap_orange()
    for r, (what, cmap, lo, hi) in enumerate((("magnitude", cm_b, -20, 70), ("error", cm_o, -20, 30))):
        y0 = 80 - 33 * r
        for i in range(4):
            if r == 1 and i == 0:
                continue
            x0 = 2 + 37 * i
            tag = "(" + "abcdefgh"[4 * r + i] + ")"
            title = PASS_TITLES[i] if r == 0 else f"error {PASS_TITLES[i]}"
            card(D, x0, y0, 35, 31, tag, title, PASS_DIGITS[i] if r == 0 else "fixed point minus exact, same pass", plot=False, fs=7.0)
            ax = D.inset(x0 + 7.0, y0 + 4.6, 26.0, 16.6, fs=4.4, grid=False)
            M = _db(fx[i]).reshape(8, 64) if r == 0 else _db(fx[i] - ex[i]).reshape(8, 64)
            ax.imshow(M, cmap=cmap, vmin=lo, vmax=hi, aspect="auto", extent=(0, 64, 8, 0), interpolation="nearest")
            for g in range(8, 64, 8):
                ax.axvline(g, color="#ffffff", lw=0.35, alpha=0.7)
            ax.set_xticks([0, 16, 32, 48, 64])
            ax.set_yticks(np.arange(8) + 0.5)
            ax.set_yticklabels([str(v) for v in range(8)])
            ax.tick_params(length=1.5, pad=1.0)
            labels(ax, "8 a1 + a0" if r == 1 else None, "a2", fs=4.6)
            if i:
                for k, col, ls in LINE_BINS:
                    row, c0, w = line_cells(k, i)
                    ax.add_patch(Rectangle((c0, row), w, 1, fill=False, ec=col, lw=0.9, ls=ls))
            if r == 0:
                D.text(x0 + 20.0, y0 + 1.8, "every line fills all 512 words" if i == 0 else f"one line: {PASS_WORDS[i]} word{'s' if PASS_WORDS[i] > 1 else ''}",
                       fs=4.8, color=INK2)
            else:
                e = np.abs(fx[i] - ex[i])
                D.text(x0 + 20.0, y0 + 1.8, f"rms {np.sqrt(np.mean(e**2)):.2f} LSB, max {e.max():.1f} LSB", fs=4.8, color=INK2)
    # ---- (e) key ------------------------------------------------------------------------------------------------------
    card(D, 2, 47, 35, 31, "(e)", "Key", plot=False, fs=7.0)
    _gradient_bar(D, 8.0, 67.0, 26.0, 1.3, cm_b, -20, 70, [-20, 0, 20, 40, 60], "|word|, dB re 1 LSB (panels a-d)")
    _gradient_bar(D, 8.0, 59.2, 26.0, 1.3, cm_o, -20, 30, [-20, -10, 0, 10, 20, 30], "|error|, dB re 1 LSB (panels f-h)")
    for j, (k, col, ls) in enumerate(LINE_BINS):
        xx = 5.0 + 8.0 * j
        D.ax.add_patch(Rectangle((xx, 52.6), 5.0, 2.0, fill=False, ec=col, lw=0.9, ls=ls))
        D.text(xx + 2.5, 50.6, f"k = {k}", fs=4.4, color=INK)
    D.text(4.0, 48.4, "outlines: words of each line (index map); dashed: mirror 512 - k", fs=3.9, ha="left", color=INK2)
    # ---- (i) spectrum after unload ------------------------------------------------------------------------------
    got = fs.natural_order(fx[-1])
    ref = np.fft.fft(x) / N
    k = np.arange(N // 2 + 1)
    card(D, 2, 2, 94, 41, "(i)", "Spectrum after the read-out, natural order", "exact transform / 512 against the fixed-point FFT and its error", plot=False, fs=7.4)
    ax = D.inset(10.0, 8.0, 83.0, 24.0, fs=4.6)
    ax.plot(k, _db(ref[k], 1e-3), color=INK, lw=0.7, label="exact FFT / 512")
    ax.plot(k, _db(got[k], 1e-3), color=BLUE, lw=0.6, label=f"fixed point, WL = {wl}")
    ax.plot(k, 20 * np.log10(np.maximum(np.abs(got - ref)[k], 1e-3)), color=ORANGE, lw=0.5, label="|error|")
    ax.axhline(10 * np.log10(fs.predicted_error_power()), color=ORANGE, lw=0.7, ls=(0, (2.5, 1.5)), label="predicted rms error per bin (noise model)")
    for kk, col in zip(fs.TWO_TONE_BINS, (ORANGE, VIOLET)):
        ax.annotate(f"k = {kk}", (kk, _db(ref[kk], 1e-3)), (kk + 12, _db(ref[kk], 1e-3) - 2), fontsize=4.6, color=INK, arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.5))
    ax.set_xlim(0, 256)
    ax.set_ylim(-20, 62)
    ax.set_xticks([0, 64, 128, 192, 256])
    labels(ax, "frequency index k", "dB re 1 LSB", fs=4.8)
    legend(ax, "upper right", fs=4.4)
    _title(ax, f"SQNR of this record: {sq:.1f} dB")
    # ---- (j) address of the line at the read-out -----------------------------------------------------------------------
    card(D, 98, 2, 50, 41, "(j)", "Read-out address, digit reversed", "frequency index k held at address a after pass 2", plot=False, fs=7.4)
    ax = D.inset(106.0, 8.0, 38.0, 24.0, fs=4.6)
    a = np.arange(N)
    ax.plot(a, output_index(a), ".", color="#b8b8b2", ms=0.9, mew=0, rasterized=False)
    for kk, col, _ in LINE_BINS:
        k0, k1, k2 = kk & 7, (kk >> 3) & 7, (kk >> 6) & 7
        ax.plot([64 * k0 + 8 * k1 + k2], [kk], "o", mfc="#ffffff", mec=col, mew=1.0, ms=3.0, zorder=5)
        ax.text(64 * k0 + 8 * k1 + k2 + 14, kk, f"{kk}", fontsize=4.4, color=INK, va="center")
    ax.set_xlim(0, 511)
    ax.set_ylim(0, 511)
    ax.set_xticks([0, 128, 256, 384, 511])
    ax.set_yticks([0, 128, 256, 384, 511])
    labels(ax, "address a = 64 a2 + 8 a1 + a0", "k = a2 + 8 a1 + 64 a0", fs=4.6)
    D.save(FIGS / "fig25_fft_pass_snapshots")



# ---------------------------------------------------------------------------------------------------------------
# fig24: index space
# ---------------------------------------------------------------------------------------------------------------
BANK_COL = SERIES  # bank j -> categorical colour j


def _lum(hex_color: str) -> float:
    h = hex_color.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def iso(x, y, z, s=1.0, cx=0.0, cy=0.0):
    """Isometric projection of the cube point ``(x, y, z) = (a0, a1, a2)`` (cell edge = ``s``): +x runs to the lower right, +y to the lower left, +z up."""
    c, h = 0.8660254037844386, 0.5
    return cx + (x - y) * c * s, cy + (-(x + y) * h + z) * s


def cube_faces(s=1.0, cx=0.0, cy=0.0):
    """The 3 x 64 visible cells of the 8 x 8 x 8 index cube: ``(polygon, (a2, a1, a0))`` for the top (a2 = 7), right (a0 = 7) and left (a1 = 7) faces."""
    out = []
    for u in range(8):
        for v in range(8):
            out.append(([iso(u, v, 8, s, cx, cy), iso(u + 1, v, 8, s, cx, cy), iso(u + 1, v + 1, 8, s, cx, cy), iso(u, v + 1, 8, s, cx, cy)], (7, v, u)))
            out.append(([iso(8, u, v, s, cx, cy), iso(8, u + 1, v, s, cx, cy), iso(8, u + 1, v + 1, s, cx, cy), iso(8, u, v + 1, s, cx, cy)], (v, u, 7)))
            out.append(([iso(u, 8, v, s, cx, cy), iso(u + 1, 8, v, s, cx, cy), iso(u + 1, 8, v + 1, s, cx, cy), iso(u, 8, v + 1, s, cx, cy)], (v, 7, u)))
    return out


def fibre(p: int, b: int):
    """The 8 words of butterfly ``b`` of pass ``p`` as ``[(a2, a1, a0), bank]`` in slot order (``pass_addresses`` re-expressed in digits)."""
    out = []
    for a in pass_addresses(p)[b]:
        a2, a1, a0 = (int(v) for v in split(a))
        out.append(((a2, a1, a0), int(bank_of(a))))
    return out


def rtl_expr(pattern: str) -> str:
    """First line of ``rtl/fft512_mem.v`` matching ``pattern``, with the declaration and the trailing comment removed."""
    import re

    from btpsd.rtlparse import RTL_DIR

    for line in (RTL_DIR / "fft512_mem.v").read_text().splitlines():
        if re.search(pattern, line):
            line = line.split("//")[0].strip().rstrip(";")
            return re.sub(r"^(wire|reg)\s*(\[[^\]]*\])?\s*", "", line)
    raise KeyError(pattern)


# lines of rtl/fft512_mem.v that fig24 (g) and fig26 (d) copy onto the sheets (regular expressions, first match per pattern)
INDEX_EXPRESSIONS = [r"wire \[2:0\] ld_bank", r"wire \[5:0\] ld_row", r"wire \[2:0\] ul_bank", r"wire \[5:0\] ul_row",
                     r"2'd0:\s+row_fn", r"2'd1:\s+row_fn", r"default: row_fn", r"2'd0:\s+t = ", r"2'd1:\s+t = "]
CORE_EXPRESSIONS = [r"wire \[2:0\] s_rd", r"wire \[2:0\]\s+m_r\b", r"wire \[5:0\]\s+raddr", r"wire \[2:0\] s1\b", r"wire \[2:0\]\s+bk\b",
                    r"assign bf_tq", r"wire \[2:0\] s_wr", r"wire \[2:0\]\s+m_w\b", r"wire \[5:0\]\s+wrow", r"wire\s+we\s+=",
                    r"wire \[DW-1:0\] rd_w"]


def _cube(D, cx, cy, s, p, b=19):
    for poly, (a2, a1, a0) in cube_faces(s, cx, cy):
        bank = (a2 + a1 + a0) % 8
        D.ax.add_patch(Polygon(poly, closed=True, fc=tint(BANK_COL[bank], 0.62), ec="#ffffff", lw=0.25, zorder=3))
    # outline of the cube
    for pts in ([iso(0, 0, 8, s, cx, cy), iso(8, 0, 8, s, cx, cy), iso(8, 8, 8, s, cx, cy), iso(0, 8, 8, s, cx, cy), iso(0, 0, 8, s, cx, cy)],
                [iso(8, 0, 8, s, cx, cy), iso(8, 0, 0, s, cx, cy), iso(8, 8, 0, s, cx, cy), iso(0, 8, 0, s, cx, cy), iso(0, 8, 8, s, cx, cy)],
                [iso(8, 8, 8, s, cx, cy), iso(8, 8, 0, s, cx, cy)]):
        xs, ys = zip(*pts)
        D.ax.plot(xs, ys, color=INK2, lw=0.5, zorder=3.5)
    fib = fibre(p, b)
    pts = [iso(a0 + 0.5, a1 + 0.5, a2 + 0.5, s, cx, cy) for (a2, a1, a0), _ in fib]
    xs, ys = zip(*pts)
    D.ax.plot(xs, ys, color=INK, lw=1.0, zorder=4)
    for (xy, (_, bank)) in zip(pts, fib):
        col = BANK_COL[bank]
        D.ax.add_patch(Circle(xy, 0.95 * s * 0.62, fc=col, ec="#ffffff", lw=0.7, zorder=5))
        D.ax.text(xy[0], xy[1] - 0.03, str(bank), fontsize=3.9, fontweight="bold", color="#ffffff" if _lum(col) < 0.55 else INK, ha="center", va="center", zorder=6)
    # axis labels along the three edges
    lx, ly = iso(4, 8.0, 0, s, cx, cy)
    D.text(lx - 3.2, ly - 1.8, "a1", fs=5.0, color=INK, bold=True)
    lx, ly = iso(8.0, 4, 0, s, cx, cy)
    D.text(lx + 3.2, ly - 1.8, "a0", fs=5.0, color=INK, bold=True)
    lx, ly = iso(0, 8, 4, s, cx, cy)
    D.text(lx - 2.6, ly, "a2", fs=5.0, color=INK, bold=True)


def fig24_index_space():
    W, H = 150, 156
    D = Diagram(W, H, "Index space of the 512-point radix-8 FFT: address digits, bank map and the three passes",
                subtitle="n = 64 n2 + 8 n1 + n0 and k = k0 + 8 k1 + 64 k2; the address a = 64 a2 + 8 a1 + a0 holds (n2, n1, n0) after the load and (k0, k1, k2) after pass 2.")
    # bank key in the header
    D.text(101.0, H - 3.0, "bank", fs=5.0, color=INK2, ha="right")
    for j in range(8):
        D.ax.add_patch(Rectangle((103.0 + 5.6 * j, H - 4.6), 4.8, 3.2, fc=tint(BANK_COL[j], 0.62), ec=BANK_COL[j], lw=0.6, zorder=3))
        D.text(105.4 + 5.6 * j, H - 3.0, str(j), fs=4.8, color=INK, bold=True)
    R1, R2, R3 = 88, 44, 2
    hi, lo = 2, 3
    b = 8 * hi + lo
    for i, (p_, dig, name) in enumerate(((0, "a2", "pass 0"), (1, "a1", "pass 1"), (2, "a0", "pass 2"))):
        x0 = 2 + 49 * i
        fib = fibre(p_, b)
        banks = " ".join(str(bk) for _, bk in fib)
        card(D, x0, R1, 47, 56, "(" + "abc"[i] + ")", f"{name}: the 8 words differ in {dig}", f"butterfly b = {b} (hi = {hi}, lo = {lo}); slot m in bank (m + {(hi + lo) % 8}) mod 8", plot=False, fs=7.0)
        _cube(D, x0 + 23.5, R1 + 27.0, 2.3, p_, b)
        D.text(x0 + 23.5, R1 + 3.4, f"banks in slot order m = 0 .. 7: {banks}", fs=5.0, color=INK)
    # ---- (d) digit flow ----------------------------------------------------------------------------------------------------
    card(D, 2, R2, 96, 42, "(d)", "What the three address digits hold", "each pass replaces one digit n by the digit k of the same weight", plot=False, fs=7.4)
    stages = [("after load", ["n2", "n1", "n0"], None), ("after pass 0", ["k0", "n1", "n0"], 0), ("after pass 1", ["k0", "k1", "n0"], 1), ("after pass 2", ["k0", "k1", "k2"], 2)]
    for j, (name, digs, changed) in enumerate(stages):
        x0 = 6.0 + 23.0 * j
        D.text(x0 + 8.6, R2 + 30.4, name, fs=5.4, color=INK, bold=True)
        for d_, dg in enumerate(digs):
            hot = changed is not None and d_ == changed
            D.box(x0 + 5.6 * d_, R2 + 22.0, 5.0, 5.2, dg, kind="compute" if hot else "memory", fs=5.4, bold=True, header=False, round_=0.7)
        D.text(x0 + 8.6, R2 + 20.3, "(a2, a1, a0)", fs=4.0, color=MUTED)
        if j < 3:
            D.arrow([(x0 + 17.4, R2 + 24.6), (x0 + 22.0, R2 + 24.6)], head=True, lw=0.9)
    notes = [("pass 0", ["8-point DFT over n2", "twiddle W512^((8 n1 + n0) k0)", "q = 2 b k, b = 8 n1 + n0"]),
             ("pass 1", ["8-point DFT over n1", "twiddle W512^(8 n0 k1)", "q = 16 lo k, lo = n0"]),
             ("pass 2", ["8-point DFT over n0", "no twiddle", ""])]
    for j, (name, rows_) in enumerate(notes):
        x0 = 6.0 + 23.0 * (j + 1)
        D.text(x0 - 5.0, R2 + 15.0, name, fs=4.8, color=INK, bold=True, ha="left")
        for r_, t in enumerate(rows_):
            D.text(x0 - 5.0, R2 + 12.4 - 2.4 * r_, t, fs=4.2, color=INK2, ha="left")
    D.text(6.0, R2 + 3.2, "read-out: the word at address a holds X[a2 + 8 a1 + 64 a0]; the frequency k is at a = 64 (k mod 8) + 8 ((k >> 3) mod 8) + (k >> 6)",
           fs=4.3, color=INK, ha="left")
    # ---- (e) worked example --------------------------------------------------------------------------------------------------
    card(D, 100, R2, 47, 42, "(e)", "Where the lines of fig25 sit", "from the index map, bins 37, 475, 150, 362", plot=False, fs=7.4)
    cols = [("k", 105.0), ("k2 k1 k0", 115.5), ("address a", 128.5), ("bank", 137.5), ("row", 143.5)]
    for lab, xx in cols:
        D.text(xx, R2 + 28.4, lab, fs=4.6, color=INK2, bold=True)
    D.arrow([(103.0, R2 + 26.8), (145.5, R2 + 26.8)], head=False, color=GRID, lw=0.7)
    for r_, row in enumerate(fs.tone_addresses([37, 475, 150, 362])):
        yy = R2 + 23.6 - 3.6 * r_
        for (lab, xx), val in zip(cols, (row["k"], f"{row['k2']} {row['k1']} {row['k0']}", row["addr"], row["bank"], row["row"])):
            D.text(xx, yy, str(val), fs=5.2, color=INK)
    D.text(104.0, R2 + 7.0, "a = 64 k0 + 8 k1 + k2,  bank = (k0 + k1 + k2) mod 8,", fs=4.2, color=INK2, ha="left")
    D.text(104.0, R2 + 4.6, "row = 8 k0 + k1  (k0 = k mod 8, k1 = (k >> 3) mod 8, k2 = k >> 6)", fs=4.2, color=INK2, ha="left")
    # ---- (f) derivation ---------------------------------------------------------------------------------------------------------
    card(D, 2, R3, 96, 40, "(f)", "Why three passes of eight-point DFTs are enough", "mixed-radix factorisation of the DFT, W = exp(-j 2 pi / 512), W8 = W^64", plot=False, fs=7.4)
    eq = [r"$W^{nk}=W_8^{\,n_2k_0}\cdot W^{(8n_1+n_0)k_0}\cdot W_8^{\,n_1k_1}\cdot W^{8n_0k_1}\cdot W_8^{\,n_0k_2}$   (the other terms contain $W^{512}=1$)",
          r"pass 0:  $A_0[k_0,n_1,n_0]=W^{(8n_1+n_0)k_0}\,\sum_{n_2}x[64n_2+8n_1+n_0]\,W_8^{\,n_2k_0}$",
          r"pass 1:  $A_1[k_0,k_1,n_0]=W^{8n_0k_1}\,\sum_{n_1}A_0[k_0,n_1,n_0]\,W_8^{\,n_1k_1}$",
          r"pass 2:  $X[k_0+8k_1+64k_2]=\sum_{n_0}A_1[k_0,k_1,n_0]\,W_8^{\,n_0k_2}$"]
    for j, t in enumerate(eq):
        D.text(5.0, R3 + 29.8 - 4.9 * j, t, fs=5.6, color=INK, ha="left")
    D.text(5.0, R3 + 10.2, "Each sum is an 8-point DFT: three radix-2 DIF stages, a rounded halving after every add/subtract, W8 and -j by the constant rotator / swap-negate.", fs=4.3, color=INK2, ha="left")
    D.text(5.0, R3 + 7.6, "A pass scales by 1/8 and the whole FFT by 1/512; the 64 butterflies of a pass read 8 words each from 8 different banks and write the same 8 addresses back.", fs=4.3, color=INK2, ha="left")
    D.text(5.0, R3 + 5.0, "In place, no second buffer. The digit reversal of the result is paid on the read-out address (unload path), not in the datapath.", fs=4.3, color=INK2, ha="left")
    # ---- (g) the same expressions in the Verilog ---------------------------------------------------------------------------------
    card(D, 100, R3, 47, 40, "(g)", "The same maps in rtl/fft512_mem.v", "copied from the source by the figure script", plot=False, fs=7.4)
    srcs = [rtl_expr(pt) for pt in INDEX_EXPRESSIONS]
    srcs = [s_.replace("2'd0:", "pass 0:").replace("2'd1:", "pass 1:").replace("default:", "pass 2:") for s_ in srcs]
    for j, t in enumerate(srcs):
        D.text(102.5, R3 + 30.4 - 2.75 * j, t, fs=3.9, color=INK, ha="left", mono=True)
    D.text(102.5, R3 + 4.8, "ld_* / ul_*: load and read-out address decode; row_fn: row of slot m;", fs=3.8, color=INK2, ha="left")
    D.text(102.5, R3 + 3.0, "twq: twiddle exponent of output m of butterfly b = {hi, lo}.", fs=3.8, color=INK2, ha="left")
    D.save(FIGS / "fig24_fft_index_space")


# ---------------------------------------------------------------------------------------------------------------
# fig27: twiddle atlas
# ---------------------------------------------------------------------------------------------------------------
def _trivial(q):
    return (np.asarray(q) % 256) == 0


def used_rom_entries():
    """ROM entries ``r`` (reduced angle, units of pi/512) that the twiddles of pass 0 and 1 and the ``W8^i`` rotations ask for."""
    from btpsd.cordic import octant_reduce

    out = {}
    for name, q in (("pass 0", pass_twiddle_q(0)), ("pass 1", pass_twiddle_q(1)), ("W8^i", np.array([128, 384]))):
        r = octant_reduce(np.unique(np.asarray(q).ravel()))[1]
        out[name] = np.unique(r[r > 0])  # r = 0 is the identity (no ROM access)
    return out


def fig27_twiddle_atlas():
    rom = default_rom()
    used = fs.twiddle_exponents_used()
    rc = fs.rotation_counts()
    oc = fs.operation_counts()
    W, H = 150, 108
    D = Diagram(W, H, "Twiddle factors of the radix-8 FFT: roots used, rotation counts and coefficient error (model)",
                subtitle="W512^e = exp(-j 2 pi e / 512); the rotator index is q = 2 e in units of pi/512. A rotation by a multiple of 90 degrees is a swap/negate; every other one "
                         "goes through the constant rotator (ROM entry, one rounding).")
    # ---- (a) the roots on the unit circle ------------------------------------------------------------------------------
    card(D, 2, 58, 52, 42, "(a)", "512th roots of unity", "which W512^e the passes use", plot=False, fs=7.4)
    ax = D.inset(10.0, 63.0, 24.0, 24.0, fs=4.4, grid=False)
    ax.set_aspect("equal")
    ang = -2 * np.pi * np.arange(512) / 512
    ax.plot(np.cos(ang), np.sin(ang), ".", color="#d6d5d0", ms=1.3, mew=0)
    p0, p1, w8 = np.array(used["pass 0"]), np.array(used["pass 1"]), np.array(used["W8"]) % 512
    quad = np.array([0, 128, 256, 384])
    sets = [(p0, dict(marker=".", ms=2.2, mfc=BLUE, mec="none", ls="none"), f"pass 0: {len(p0)} roots"),
            (p1, dict(marker="o", ms=3.6, mfc="none", mec=ORANGE, mew=0.8, ls="none"), f"pass 1: {len(p1)} roots"),
            (w8, dict(marker="D", ms=4.4, mfc="none", mec=AQUA, mew=0.9, ls="none"), "W8^i: 8 roots"),
            (quad, dict(marker="s", ms=3.0, mfc=INK, mec="none", ls="none"), "+-1, +-j: free")]
    for e, st, lab in sets:
        ax.plot(np.cos(-2 * np.pi * e / 512), np.sin(-2 * np.pi * e / 512), **st)
    ax.set_xlim(-1.15, 1.15)
    ax.set_ylim(-1.15, 1.15)
    ax.set_xticks([-1, 0, 1])
    ax.set_yticks([-1, 0, 1])
    ax.tick_params(length=1.5, pad=1.0)
    for j, (e, st, lab) in enumerate(sets):
        yy = 80.0 - 4.6 * j
        st2 = dict(st)
        st2["ms"] = min(st2["ms"], 3.4)
        D.ax.plot([37.5], [yy], **st2, zorder=6)
        D.text(40.0, yy, lab, fs=4.6, ha="left", color=INK)
    D.text(37.0, 62.5, "pass 0: e = b k <= 441, so the last\nsector of the circle is never visited", fs=4.0, ha="left", color=INK2)
    # ---- (b) exponents of pass 0 ---------------------------------------------------------------------------------------
    q0 = pass_twiddle_q(0)
    card(D, 56, 58, 46, 42, "(b)", "Pass 0: exponent e", "e = b k mod 512, b = 8 n1 + n0", plot=False, fs=7.4)
    ax = D.inset(63.0, 74.5, 36.0, 16.0, fs=4.4, grid=False)
    ax.imshow((q0.T // 2), cmap=cmap_blue(), vmin=0, vmax=511, aspect="auto", extent=(0, 64, 8, 0), interpolation="nearest")
    ty, tx = np.nonzero(_trivial(q0.T))
    ax.plot(tx + 0.5, ty + 0.5, ".", color=INK, ms=1.0, mew=0)
    ax.set_xticks([0, 16, 32, 48, 64])
    ax.set_yticks(np.arange(8) + 0.5)
    ax.set_yticklabels([str(i) for i in range(8)])
    labels(ax, "butterfly b", "output k", fs=4.6)
    ax.tick_params(length=1.5, pad=1.0)
    D.text(79.0, 66.2, f"{rc['twiddle_per_pass'][0]} of 512 need the rotator; the {512 - rc['twiddle_per_pass'][0]} dots are 1 or -j", fs=4.3, color=INK2)
    _gradient_bar(D, 67.0, 61.0, 28.0, 1.0, cmap_blue(), 0, 511, [0, 128, 256, 384, 511], "e", fs=4.0)
    # ---- (c) exponents of pass 1 ---------------------------------------------------------------------------------------
    q1 = pass_twiddle_q(1)
    card(D, 104, 58, 44, 42, "(c)", "Pass 1: exponent e", "e = 8 lo k mod 512, lo = b mod 8", plot=False, fs=7.4)
    ax = D.inset(111.0, 74.5, 34.0, 16.0, fs=4.4, grid=False)
    e1 = np.array([[(q1[lo, k] // 2) for lo in range(8)] for k in range(8)])  # rows k, columns lo (the table does not depend on hi)
    ax.imshow(e1, cmap=cmap_blue(), vmin=0, vmax=511, aspect="auto", extent=(0, 8, 8, 0), interpolation="nearest")
    for k in range(8):
        for lo in range(8):
            triv = (2 * e1[k, lo]) % 256 == 0
            ax.text(lo + 0.5, k + 0.5, str(e1[k, lo]), fontsize=3.6, ha="center", va="center", color=MUTED if triv else INK, fontweight="normal" if triv else "bold")
    ax.set_xticks(np.arange(8) + 0.5)
    ax.set_xticklabels([str(i) for i in range(8)])
    ax.set_yticks(np.arange(8) + 0.5)
    ax.set_yticklabels([str(i) for i in range(8)])
    labels(ax, "lo", "output k", fs=4.6)
    ax.tick_params(length=1.5, pad=1.0)
    D.text(126.0, 66.2, f"{rc['twiddle_per_pass'][1]} of 512 need the rotator (grey: 1 or -j)", fs=4.3, color=INK2)
    D.text(126.0, 62.2, "pass 2 has no twiddle", fs=4.3, color=INK2)
    # ---- (d) coefficient error of the ROM rotator -----------------------------------------------------------------------
    card(D, 2, 2, 78, 54, "(d)", "Error of the ROM entries", "angle and gain error of the coefficient pair (a, b, inv) of entry r (rotation by r pi/512)", plot=False, fs=7.4)
    ax = D.inset(10.0, 17.0, 66.0, 26.0, fs=4.6)
    r = np.arange(rom.n)
    ae, ge = rom.angle_error() * 1e3, rom.gain_error() * 1e3
    ax.axhspan(-PAPER_TOL * 1e3, PAPER_TOL * 1e3, color="#eef3fb", lw=0)
    ax.vlines(r, 0, ae, color=BLUE, lw=0.6)
    ax.plot(r, ae, "o", ms=1.5, mfc=BLUE, mec="none", label="angle error (mrad)")
    ax.plot(r, ge, "-", color=ORANGE, lw=0.7, label="gain error (1e-3)")
    ax.axhline(PAPER_TOL * 1e3, color=INK2, lw=0.5, ls=(0, (2.5, 1.5)))
    ax.axhline(-PAPER_TOL * 1e3, color=INK2, lw=0.5, ls=(0, (2.5, 1.5)))
    ax.text(127, PAPER_TOL * 1e3 + 0.5, f"tolerance {PAPER_TOL} rad", ha="right", va="bottom", fontsize=4.4, color=INK2)
    ax.set_xlim(-1, 129)
    ax.set_ylim(-13, 15)
    ax.set_xticks([0, 16, 32, 48, 64, 80, 96, 112, 128])
    labels(ax, None, "milli-units", fs=4.8)
    legend(ax, "lower left", fs=4.4, ncol=2)
    ang_err = rom.angle_error()
    _title(ax, f"{rom.n} entries; max angle error {np.max(np.abs(ang_err)):.4f} rad, rms {np.sqrt(np.mean(ang_err ** 2)):.4f} rad, max gain error {np.max(np.abs(rom.gain_error())):.1e}")
    rug = D.inset(10.0, 8.6, 66.0, 6.0, fs=4.2, grid=False)
    ue = used_rom_entries()
    rows_ = [("pass 0", ue["pass 0"], BLUE), ("pass 1", ue["pass 1"], ORANGE), ("W8^i", ue["W8^i"], AQUA)]
    for j, (lab, rr, col) in enumerate(rows_):
        rug.vlines(rr, j + 0.15, j + 0.85, color=col, lw=0.6)
    rug.set_xlim(-1, 129)
    rug.set_ylim(0, 3)
    rug.set_yticks([0.5, 1.5, 2.5])
    rug.set_yticklabels([row[0] for row in rows_])
    rug.set_xticks([0, 16, 32, 48, 64, 80, 96, 112, 128])
    labels(rug, "ROM entry r: which entries the three stages read (r = 0 is the identity)", None, fs=4.6)
    # ---- (e) rotation counts -----------------------------------------------------------------------------------------------
    card(D, 82, 2, 66, 54, "(e)", "Rotations and memory traffic of one FFT", "exact counts of this design against an in-place radix-2 FFT", plot=False, fs=7.4)
    ax = D.inset(98.0, 33.0, 46.0, 12.0, fs=4.6)
    ax.barh([2], [rc["direct"]], color=MUTED, height=0.6)
    ax.barh([1], [rc["radix2"]], color=BLUE, height=0.6)
    left = 0
    for val, col, lab in ((rc["twiddle_per_pass"][0], ORANGE, "pass 0 twiddles"), (rc["twiddle_per_pass"][1], AQUA, "pass 1 twiddles"), (rc["w8"], VIOLET, "W8^i")):
        ax.barh([0], [val], left=left, color=col, height=0.6, label=lab)
        left += val
    for yy, val in ((2, rc["direct"]), (1, rc["radix2"]), (0, rc["radix8"])):
        ax.text(val * 1.25, yy, f"{val:,}".replace(",", " "), va="center", fontsize=4.8, color=INK)
    ax.set_xscale("log")
    ax.set_xlim(100, 2e6)
    ax.set_yticks([0, 1, 2])
    ax.set_yticklabels(["this design", "radix-2, in place", "direct DFT"])
    ax.set_ylim(-0.5, 2.5)
    labels(ax, "rotations that need the rotator (not +-1, +-j)", None, fs=4.6)
    legend(ax, "lower right", fs=4.0, bbox_to_anchor=(1.0, 0.0))
    rows = [("", "radix-2", "this design"), ("passes over the memory", oc["radix2_passes"], oc["passes"]), ("butterflies", oc["radix2_butterflies"], oc["butterflies"]),
            ("words read", oc["radix2_words_read"], oc["words_read"]), ("words written", oc["radix2_words_written"], oc["words_written"]),
            ("complex add / sub", oc["radix2_addsub"], oc["complex_addsub"])]
    for j, (lab, a_, b_) in enumerate(rows):
        yy = 25.0 - 3.0 * j
        D.text(86.0, yy, lab, fs=5.0, ha="left", color=INK, bold=(j == 0))
        D.text(124.0, yy, f"{a_:,}".replace(",", " ") if not isinstance(a_, str) else a_, fs=5.0, ha="right", color=INK, bold=(j == 0))
        D.text(144.0, yy, f"{b_:,}".replace(",", " ") if not isinstance(b_, str) else b_, fs=5.0, ha="right", color=INK, bold=(j == 0))
    D.ax.plot([85.0, 145.0], [26.6, 26.6], color=GRID, lw=0.7, zorder=2)
    D.text(86.0, 6.4, "Same number of adders, a third of the memory accesses. The rotations drop by a fifth (output twiddles of the", fs=4.2, ha="left", color=INK2)
    D.text(86.0, 4.0, "radix-8 passes plus two 45-degree W8 rotations per butterfly); the gain is 3 sweeps instead of 9.", fs=4.2, ha="left", color=INK2)
    D.save(FIGS / "fig27_twiddle_atlas")


# ---------------------------------------------------------------------------------------------------------------
# fig28: accuracy
# ---------------------------------------------------------------------------------------------------------------
INPUT_COL = {"white noise": BLUE, "full-scale tone": ORANGE, "six tones": AQUA, "EEG-like": VIOLET}


def _acc(rows):
    out = {}
    for r in rows:
        out.setdefault((r[0], r[2]), []).append(r)
    return out


def worst_overshoot(rows):
    """Largest ``measured - predicted`` SQNR (dB) over the rounding rows of the accuracy table, with its input and word length."""
    best = max((r for r in rows if r[2] == "round"), key=lambda r: r[3] - r[4])
    return best[3] - best[4], best[0], best[1]


def fig28_accuracy():
    rom = default_rom()
    rows = fs.accuracy_rows()
    acc = _acc(rows)
    wls = list(fs.ACCURACY_WLS)
    W, H = 150, 112
    D = Diagram(W, H, "Accuracy of the fixed-point FFT: noise model, coefficient ceiling and where the error goes (model)",
                subtitle="SQNR against the exact transform of the same integer input, scaled by 1/512; synthetic records of 512 samples at peak 0.9 of full scale. "
                         "Dashed curves are the analytic prediction (noise model plus ROM floor), markers are bit-accurate simulations.")
    ch = 48  # card height
    # ---- (a) SQNR vs WL ---------------------------------------------------------------------------------------------------
    card(D, 2, 52, 48, ch, "(a)", "SQNR against word length", "rounding to nearest; dashed: prediction", plot=False, fs=7.4)
    ax = D.inset(10.0, 58.5, 36.5, 28.5, fs=4.6)
    for kind, col in INPUT_COL.items():
        rr = acc[(kind, "round")]
        ax.plot(wls, [r[3] for r in rr], "o-", color=col, lw=0.8, ms=2.4, label=kind)
        ax.plot(wls, [r[4] for r in rr], color=col, lw=0.8, ls=(0, (2.5, 1.5)))
    rom_db = [acc[(k, "round")][0][6] for k in fs.ACCURACY_INPUTS]
    ax.axhspan(min(rom_db), max(rom_db), color="#e9e8e3", lw=0, zorder=0)
    ax.text(24.2, min(rom_db) - 1.0, "ROM ceiling", fontsize=4.4, color=INK2, va="top", ha="right")
    ax.set_xlim(7.5, 24.5)
    ax.set_ylim(0, 60)
    ax.set_xticks(wls)
    labels(ax, "word length WL (bits)", "SQNR (dB)", fs=4.8)
    legend(ax, "upper left", fs=4.2)
    # ---- (b) what rounding buys -----------------------------------------------------------------------------------------------
    card(D, 52, 52, 48, ch, "(b)", "What rounding buys over truncation", "SQNR(rounding) - SQNR(truncation), measured", plot=False, fs=7.4)
    ax = D.inset(60.0, 58.5, 36.5, 28.5, fs=4.6)
    for kind, col in INPUT_COL.items():
        d_ = [a_[3] - b_[3] for a_, b_ in zip(acc[(kind, "round")], acc[(kind, "trunc")])]
        ax.plot(wls, d_, "o-", color=col, lw=0.8, ms=2.4, label=kind)
    ax.axhline(0, color=INK2, lw=0.6)
    ax.set_xlim(7.5, 24.5)
    ax.set_xticks(wls)
    labels(ax, "word length WL (bits)", "dB", fs=4.8)
    legend(ax, "upper right", fs=4.2)
    # ---- (c) error per bin --------------------------------------------------------------------------------------------------
    card(D, 102, 52, 46, ch, "(c)", "Error per bin, WL = 12", "mean square over 60 white-noise records", plot=False, fs=7.4)
    ax = D.inset(110.0, 58.5, 34.5, 28.5, fs=4.6)
    k = np.arange(N)
    for rnd, col, lab in (("round", BLUE, "rounding"), ("trunc", ORANGE, "truncation")):
        mean_e, mse = fs.mean_error_per_bin(12, rnd, 60, "white noise", rom)
        ax.plot(k, 10 * np.log10(mse), color=col, lw=0.45, label=lab)
    ax.axhline(10 * np.log10(fs.predicted_error_power()), color=INK2, lw=0.6, ls=(0, (2.5, 1.5)), label="prediction, mean over bins")
    ax.set_xlim(0, 511)
    ax.set_ylim(-4, 22)
    ax.set_xticks([0, 128, 256, 384, 511])
    labels(ax, "bin k", "|error|^2 (dB re 1 LSB^2)", fs=4.8)
    legend(ax, "upper center", fs=4.2)
    # ---- (d) error per pass -------------------------------------------------------------------------------------------------
    card(D, 2, 2, 48, ch, "(d)", "Where the error accrues", "mean-square error per word after each pass, six tones", plot=False, fs=7.4)
    for j, wl_ in enumerate((10, 16)):
        ax = D.inset(10.0 + 20.0 * j, 8.0, 15.5, 28.5, fs=4.4)
        re, im = fs.fft_inputs("six tones", wl_, 1)
        pe = fs.pass_error_power(re, im, wl_, rom)
        xs = np.arange(1, 4)
        ax.bar(xs - 0.2, np.maximum(pe["arith"][1:], 1e-4), 0.36, color=BLUE, label="rounding in the datapath")
        ax.bar(xs + 0.2, np.maximum(pe["rom"][1:], 1e-4), 0.36, color=ORANGE, label="ROM coefficient error")
        ax.set_yscale("log")
        ax.set_ylim(1e-3, 1e3)
        ax.set_xticks(xs)
        ax.set_xticklabels(["p0", "p1", "p2"])
        _title(ax, f"WL = {wl_}")
        if j == 0:
            labels(ax, None, "LSB^2 per word", fs=4.8)
            legend(ax, "upper left", fs=3.8)
        else:
            ax.set_yticklabels([])
    # ---- (e) SQNR vs level --------------------------------------------------------------------------------------------------
    card(D, 52, 2, 48, ch, "(e)", "SQNR against input level", "white noise, rounding; level = rms below full scale", plot=False, fs=7.4)
    ax = D.inset(60.0, 8.0, 36.5, 28.5, fs=4.6)
    for wl_, col in zip((10, 12, 14, 16), (SEQ_BLUE[2], SEQ_BLUE[3], SEQ_BLUE[4], SEQ_BLUE[6])):
        ax.plot(fs.LEVELS_DB, fs.sqnr_vs_level(wl_, rom=rom), "o-", color=col, lw=0.8, ms=2.2, label=f"WL = {wl_}")
    ax.axhline(np.mean(rom_db), color=ORANGE, lw=0.6, ls=(0, (2.5, 1.5)))
    ax.text(-65, np.mean(rom_db) + 1.0, "ROM ceiling", fontsize=4.4, color=ORANGE, va="bottom")
    ax.set_xlim(-68, -10)
    ax.set_ylim(-45, 55)
    ax.set_xticks(range(-66, -9, 12))
    labels(ax, "input level (dB re full scale)", "SQNR (dB)", fs=4.8)
    legend(ax, "lower right", fs=4.2)
    # ---- (f) noise model -----------------------------------------------------------------------------------------------------
    card(D, 102, 2, 46, ch, "(f)", "Noise model", "random part and rounding bias, per complex bin", plot=False, fs=7.4)
    var, bias, tot = fs.predicted_error_variance(), fs.predicted_bias_power(), fs.predicted_error_power()
    over, over_kind, over_wl = worst_overshoot(rows)
    lines_ = [
        ("every add/sub is followed by a rounded >> 1:", True),
        ("error 0 or +1/2 LSB (mean 1/4, variance 1/16);", False),
        ("noise injected after radix-2 stage s reaches", False),
        ("each output with weight 2^-(9 - s).", False),
        ("", False),
        (f"random part   2 x {var:.4f} = {2 * var:.3f} LSB^2", False),
        (f"bias part   {bias:.3f} LSB^2", False),
        (f"total   {tot:.3f} LSB^2   (rms {np.sqrt(tot):.2f} LSB)", True),
        ("", False),
        ("SQNR = 10 log10( sum |X|^2 / (512 x total) ),", True),
        ("combined with the ROM floor r by", False),
        ("1 / (10^(-q/10) + 10^(-r/10)).", False),
        ("", False),
        (f"Tone records measure up to {over:.1f} dB above the", False),
        (f"prediction ({over_kind}, WL = {over_wl}); noise-like", False),
        ("records agree to within 0.6 dB.", False),
    ]
    for j, (t, bold) in enumerate(lines_):
        D.text(104.5, 38.5 - 2.2 * j, t, fs=4.7, ha="left", color=INK if bold else INK2, bold=bold)
    D.save(FIGS / "fig28_fft_accuracy")
    with open(RESULTS / "model_fft_accuracy.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(fs.ACCURACY_HEADER)
        w.writerows(rows)


# ---------------------------------------------------------------------------------------------------------------
# fig26: the FFT core at signal level (drawn from the Verilog sources)
# ---------------------------------------------------------------------------------------------------------------
def _src(name: str) -> str:
    """Source of ``rtl/<name>.v`` with the comments removed."""
    t = (RTL_DIR / f"{name}.v").read_text()
    t = re.sub(r"/\*.*?\*/", "", t, flags=re.S)
    return re.sub(r"//[^\n]*", "", t)


def _expr(pattern: str) -> str:
    return " ".join(rtl_expr(pattern).split())


def core_facts():
    """The structure fig26 draws, read from ``rtl/fft512_mem.v``, ``r8_butterfly.v`` and ``btpsd_params.vh``; the figure asserts its drawing against it."""
    mods = parse_all()
    core = {i.name: i for i in mods["fft512_mem"].instances}
    bfly = {i.name: i for i in mods["r8_butterfly"].instances}
    consts = define_values(RTL_DIR / "btpsd_params.vh")
    core_src, bfly_src = _src("fft512_mem"), _src("r8_butterfly")
    q = {int(i): int(v) for i, v in re.findall(r"u_w(\d)\s*\([^;]*?\.q\(10'd(\d+)\)", bfly_src)}
    assert core["u_ram"].count == 8 and core["u_ram"].params == {"DW": "DW", "AW": "6"}
    assert core["u_meta1"].params == {"W": "9", "N": "1"} and core["u_meta2"].params == {"W": "9", "N": "LB"}
    assert core["u_bf"].module == "r8_butterfly" and core["u_bf"].count == 1
    assert bfly["u_tw"].count == 7 and all(bfly[f"u_w{i}"].module == "cru_rotate" for i in (1, 2, 3))
    assert bfly["u_dl_s"].params == {"W": "10*WL", "N": "LC"} and bfly["u_dl_y0"].params == {"W": "2*WL", "N": "LC"}
    assert bfly["u_dl_tq"].params == {"W": "70", "N": "LC + 2"}
    assert q == {1: 128, 2: 256, 3: 384}
    assert re.search(r"localparam\s+LB\s*=\s*`LAT_BFLY", core_src) and re.search(r"localparam\s+LAT\s*=\s*LB\s*\+\s*2", core_src)
    assert re.search(r"STALL\s*=\s*\(LAT\s*>\s*8\)\s*\?\s*\(LAT\s*-\s*8\)\s*:\s*0", core_src)
    assert re.search(r"localparam\s+LC\s*=\s*`LAT_CRU", bfly_src)
    lc, lb = consts["LAT_CRU"], consts["LAT_BFLY"]
    budget = fs.fft_cycle_budget()
    assert (budget["lat_cru"], budget["lat_bfly"], budget["lat"]) == (lc, lb, lb + 2) and budget["lat_bfly"] == 1 + lc + 1 + lc
    n_cru = sum(i.count for i in bfly.values() if i.module == "cru_rotate")
    n_dl = sum(i.count for i in list(core.values()) + list(bfly.values()) if i.module == "delay_line")
    return dict(core=core, bfly=bfly, q=q, lc=lc, lb=lb, lat=lb + 2, stall=budget["stall"], n_cru=n_cru, n_dl=n_dl, budget=budget)


def sim_done_cycle() -> int:
    """Cycle (0 = ``start`` high) in which ``done`` is first high in the RTL simulation stored in ``results/rtl_trace_fft512.csv`` (CI run)."""
    with open(RESULTS / "rtl_trace_fft512.csv") as f:
        rows = list(csv.DictReader(f))
    return next(int(r["cycle"]) for r in rows if r["done"] == "1")


def _tap(D, pts, label=None, **kw):
    """Control wire: thin dashed violet arrow."""
    D.arrow(pts, label=label, color=VIOLET, lw=0.85, ls=(0, (3.2, 2.0)), **kw)


def _wire(D, pts, label=None, bus=None, **kw):
    """Data bus: solid arrow."""
    D.arrow(pts, label=label, color=INK2, lw=1.35, bus=bus, **kw)


def _reg_line(D, x, ytop, ybot, label, fs=5.0):
    D.ax.plot([x, x], [ybot, ytop], color=INK2, lw=0.6, ls=(0, (3, 2)), zorder=2.5)
    D.text(x, ytop + 1.5, label, fs=fs, color=INK, bold=True)


def fig26_fft_core_blocks():
    F = core_facts()
    lc, lb, lat, stall = F["lc"], F["lb"], F["lat"], F["stall"]
    done = sim_done_cycle()
    assert done == F["budget"]["trace_done_cycle"]
    sig = fs.butterfly_signals(1, 19)
    W, H = 158, 187
    D = Diagram(W, H, "The FFT core at signal level: fft512_mem and its radix-8 butterfly",
                subtitle="Instance names, parameters and latencies are read from rtl/fft512_mem.v, r8_butterfly.v and cru_rotate.v (btpsd.rtlparse); clk and the asynchronous\n"
                         "active-low reset rst_n are omitted; bus widths are given in WL (word length), DW = 2 WL.")

    # ============================== (a) data loop and control pipeline ==============================================
    AY = 99.0
    card(D, 2, AY, 154, 75, "(a)", "Data loop (solid) and control pipeline (dashed)",
         "one butterfly per cycle: (1) issue in cycle t, (2) bank read and butterfly input in t + 1, (3) butterfly output and (4) write in t + 11; the written words can be read again from t + 12",
         plot=False, fs=7.6)
    y = lambda v: AY + v  # noqa: E731 - panel-relative height
    D.box(5, y(50), 36, 15.5, "issue controller\nst: ST_IDLE / ST_RUN / ST_DRAIN\niss_v, iss_p[1:0], iss_b[5:0]\nstall_cnt (STALL), drain_cnt (LAT - 2)", kind="control", fs=7.0, step=1)
    D.text(5.8, y(48.0), "start → busy … done", fs=5.0, color=INK2, ha="left", mono=True)
    m1, m2 = F["core"]["u_meta1"], F["core"]["u_meta2"]
    D.box(54, y(53), 26, 9, f"u_meta1\ndelay_line  W = {m1.params['W']}, N = {m1.params['N']}", kind="control", fs=6.8)
    D.box(94, y(53), 28, 9, f"u_meta2\ndelay_line  W = {m2.params['W']}, N = LB = {lb}", kind="control", fs=6.8)
    _tap(D, [(41, y(57.5)), (54, y(57.5))], "iss_v, iss_p, iss_b", loff=(0, 2.0), fs=5.4)
    _tap(D, [(80, y(57.5)), (94, y(57.5))], "v1, p1, b1", loff=(0, 2.0), fs=5.4)
    D.text(47.5, y(55.4), "cycle t", fs=4.8, color=MUTED)
    D.text(82.4, y(55.4), "t + 1", fs=4.8, color=MUTED, ha="left")
    _tap(D, [(122, y(57.5)), (153.5, y(57.5)), (153.5, y(9.5)), (130, y(9.5))])
    D.text(138.0, y(59.6), "vw, pw, bw  (t + 11)", fs=5.4, color=INK2)
    # control taps
    D.junction(46, y(57.5))
    _tap(D, [(46, y(57.5)), (46, y(37.5)), (52, y(37.5))])
    D.text(47.0, y(48.0), "iss_*", fs=4.8, color=VIOLET, ha="left", mono=True)
    D.text(46.9, y(39.0), "raddr", fs=4.6, color=VIOLET, ha="left", mono=True)
    D.junction(87, y(57.5))
    _tap(D, [(87, y(57.5)), (87, y(48.5)), (99, y(48.5)), (99, y(43.3))])
    D.junction(99, y(48.5))
    _tap(D, [(99, y(48.5)), (134, y(48.5)), (134, y(44.3))])
    D.text(100.4, y(46.0), "s1", fs=5.2, color=VIOLET, ha="left", mono=True, bold=True)
    D.text(116.5, y(50.2), "bf_tq = twq(p1, b1, m)", fs=5.0, color=VIOLET, mono=True)

    # host ports
    D.box(5, y(26), 36, 17.5, "load port\nld_en, ld_addr[8:0]\nld_re, ld_im  [WL]\nld_bank, ld_row from ld_addr", kind="io", fs=6.8)
    D.box(5, y(6), 36, 17.5, "unload port\nrd_k[8:0] → rd_re, rd_im  [WL]\nul_bank, ul_row from rd_k\nrd_w = rq_all[ul_bank_q]", kind="io", fs=6.8)
    # memory
    D.box(52, y(7), 28, 37, "g_bank[j].u_ram  × 8\n", kind="memory", fs=6.6, step=2)
    for j in range(8):
        yy = y(7.6) + 4.1 * (7 - j)
        D.ax.add_patch(Rectangle((52.6, yy), 26.8, 3.8, fc=tint(BANK_COL[j], 0.6), ec=BANK_COL[j], lw=0.5, zorder=3))
        D.text(54.0, yy + 1.9, f"bank {j}", fs=5.2, color=INK, bold=True, ha="left")
        D.text(78.6, yy + 1.9, "64 × DW", fs=4.8, color=INK2, ha="right")
    _wire(D, [(41, y(31.0)), (52, y(31.0))])
    D.text(46.5, y(32.7), "we, waddr, wdata", fs=4.4, color=INK2, mono=True)
    _wire(D, [(52, y(13.0)), (41, y(13.0))])
    D.text(46.5, y(14.7), "rq_all", fs=4.8, color=INK2, mono=True)
    # read rotate, butterfly, write rotate
    D.box(89, y(25), 20, 18, "read rotate\nslot m ← bank\n(m + s1) mod 8\ng_slot[m] × 8", kind="compute", fs=6.4)
    D.box(118, y(19), 32, 25, f"u_bf\nr8_butterfly\nWL, ROUND\nlatency LB = {lb} cycles\nxr, xi: 8 × WL\ntq: 7 × 10 bit", kind="compute", fs=7.0, step=3)
    D.box(89, y(8), 41, 8, "write rotate\nbank j ← slot (j - s_wr) mod 8", kind="compute", fs=6.4, step=4)
    _wire(D, [(80, y(34)), (89, y(34))], bus="8×2WL")
    _wire(D, [(109, y(34)), (118, y(34))], bus="8×WL")
    D.text(113.5, y(37.8), "bf_xr, bf_xi", fs=4.8, color=INK2, mono=True)
    _wire(D, [(140, y(19)), (140, y(13.5)), (130, y(13.5))])
    D.text(142.0, y(16.4), "bf_yr, bf_yi", fs=4.8, color=INK2, mono=True, ha="left")
    _wire(D, [(89, y(12.0)), (80, y(12.0))])
    D.text(84.5, y(13.7), "wdata", fs=4.6, color=INK2, mono=True)
    D.text(5.0, y(3.2), "run_mode = (st != ST_IDLE): the bank ports use pass addressing while a transform runs, the load / unload ports while idle", fs=4.9, color=INK2, ha="left", italic=True)
    lx, ly = 103.0, y(3.2)
    for i, (kind, lab) in enumerate((("memory", "memory"), ("compute", "arithmetic"), ("control", "control"), ("io", "host port"))):
        col = {"memory": PAL["blue"], "compute": PAL["aqua"], "control": PAL["violet"], "io": INK2}[kind]
        D.ax.add_patch(FancyBboxPatch((lx + 13.0 * i, ly - 0.8), 2.4, 1.6, boxstyle="round,pad=0,rounding_size=0.4",
                                      fc=tint(col, 0.78) if kind != "io" else "#e7e6e1", ec=col, lw=0.8))
        D.text(lx + 13.0 * i + 3.1, ly, lab, fs=4.8, ha="left")

    # ============================== (b) inside the butterfly ============================================================
    card(D, 2, 45, 154, 51, "(b)", "Inside u_bf (r8_butterfly): one butterfly per cycle through a 10-cycle pipeline",
         f"LAT_BFLY = 1 (stage 1 register) + LAT_CRU + 1 (stages 2 + 3 register) + LAT_CRU = 1 + {lc} + 1 + {lc} = {lb} cycles from x valid (t + 1) to y valid (t + 11)",
         plot=False, fs=7.6)
    bf, q = F["bfly"], F["q"]
    D.box(4, 60, 11, 26, "x\nxr, xi\n8 slots\n× WL", kind="io", fs=6.4)
    D.box(18, 60, 20, 26, "stage 1\nhadd / hsub, distance 4\ns_i = x_i + x_(i+4)\nd_i = x_i - x_(i+4)\neach result >> 1", kind="compute", fs=7.0)
    lane_y = [82.3, 74.3, 67.9, 61.5]
    D.box(43, 78.6, 28, 7.4, f"u_dl_s  carries s0..s3, d0\ndelay_line  W = 10·WL, N = LC = {lc}", kind="control", fs=6.2)
    for i, yy in zip((1, 2, 3), (71.4, 65.0, 58.6)):
        D.box(43, yy, 28, 5.8, f"u_w{i}  cru_rotate  q = {q[i]}  (d{i} · W8^{i})", kind="compute", fs=6.0)
    for yy in lane_y:
        D.arrow([(38, yy), (43, yy)], lw=1.0)
        D.arrow([(71, yy), (76, yy)], lw=1.0)
    D.text(40.5, 86.9, "s, d0", fs=4.4, color=INK2) if False else None
    D.box(76, 60, 26, 26, "stages 2 + 3\ngrp 0: s branch (u_dl_s)\ngrp 1: rotated d branch\nc0 .. c3: hadd / hsub\n(-j): swap and negate\ny3r, y3i: one register", kind="compute", fs=7.0)
    D.box(107, 78.6, 28, 7.4, f"u_dl_y0  slot 0 of y3\ndelay_line  W = 2·WL, N = LC = {lc}", kind="control", fs=6.2)
    D.box(107, 59.2, 28, 17.6, f"g_tw[k].u_tw  × {bf['u_tw'].count}\ncru_rotate, k = 1 .. 7\nq = tqd[(k-1)*10 +: 10]\nslot k of y3 → slot k of y", kind="compute", fs=7.0)
    D.arrow([(102, 82.3), (107, 82.3)], lw=1.0)
    D.arrow([(102, 68.0), (107, 68.0)], lw=1.0)
    D.arrow([(135, 82.3), (138, 82.3)], lw=1.0)
    D.arrow([(135, 68.0), (138, 68.0)], lw=1.0)
    D.box(138, 60, 15, 26, "y\nyr, yi\n8 slots\n× WL", kind="io", fs=6.4)
    # twiddle path
    tq_y = 56.2
    D.box(80, 54.0, 22, 4.4, f"u_dl_tq  N = LC + 2 = {lc + 2}", kind="control", fs=5.6)
    _tap(D, [(4, tq_y), (80, tq_y)])
    _tap(D, [(102, tq_y), (121, tq_y), (121, 59.2)])
    D.text(5.0, tq_y + 1.7, "tq[69:0]: 7 exponents, units of π/512", fs=4.8, color=VIOLET, ha="left")
    # time axis
    xs = [15.0, 40.0, 73.0, 104.0, 138.0]
    ticks = ["t + 1", "t + 2", "t + 6", "t + 7", "t + 11"]
    D.ax.plot([xs[0], xs[-1]], [50.2, 50.2], color=INK, lw=0.9, zorder=3)
    for xx, tk in zip(xs, ticks):
        D.ax.plot([xx, xx], [49.4, 51.0], color=INK, lw=0.9, zorder=3)
        D.text(xx, 47.4, tk, fs=5.6, color=INK, bold=True)
    for xx in xs[1:-1]:
        D.ax.plot([xx, xx], [51.4, 87.4], color=INK2, lw=0.6, ls=(0, (3, 2)), zorder=2.4)
    segs = [(0, "stage 1, R1: 1"), (1, f"u_w1..3, u_dl_s: LAT_CRU = {lc}"), (2, "stages 2 + 3: 1"), (3, f"g_tw, u_dl_y0: LAT_CRU = {lc}")]
    for k, lab in segs:
        D.text((xs[k] + xs[k + 1]) / 2, 52.4, lab, fs=5.0, color=INK2)
    D.text(4.0, 47.4, "cycle", fs=5.0, color=INK2, ha="left")

    # ============================== (c) worked example ==================================================================
    card(D, 2, 2, 80, 40, "(c)", "One butterfly through the signals",
         f"pass 1, b = {sig['b']} (hi = {sig['hi']}, lo = {sig['lo']}): s = (hi + lo) mod 8 = {sig['s']}, slot m reads bank (m + {sig['s']}) mod 8", plot=False, fs=7.6)
    cols = [("slot m", 8.0), ("bank", 19.0), ("row", 30.5), ("a2 a1 a0", 44.5), ("a", 57.5), ("tq", 67.5), ("angle", 77.0)]
    for lab, xx in cols:
        D.text(xx, 31.4, lab, fs=5.0, color=INK2, bold=True)
    D.ax.plot([4.0, 80.0], [29.9, 29.9], color=GRID, lw=0.7, zorder=3)
    for r, sl in enumerate(sig["slots"]):
        yy = 27.7 - 2.55 * r
        j = sl["bank"]
        D.ax.add_patch(Rectangle((16.9, yy - 1.05), 4.2, 2.1, fc=tint(BANK_COL[j], 0.55), ec=BANK_COL[j], lw=0.5, zorder=3))
        D.text(19.0, yy, str(j), fs=5.0, color=INK, bold=True)
        D.text(8.0, yy, str(sl["m"]), fs=5.4, color=INK)
        D.text(30.5, yy, f"{{{sig['hi']}, {sl['m']}}} = {sl['row']}", fs=5.2, color=INK)
        D.text(44.5, yy, f"{sl['a2']} {sl['a1']} {sl['a0']}", fs=5.4, color=INK)
        D.text(57.5, yy, str(sl["addr"]), fs=5.4, color=INK)
        D.text(67.5, yy, str(sl["tq"]), fs=5.4, color=INK)
        D.text(77.0, yy, f"{sl['angle_deg']:.1f}°", fs=5.2, color=INK)
    D.text(4.0, 7.0, f"tq = 16 lo m in pass 1; output k of the butterfly is rotated by exp(-j tq π/512), slot 0 is not rotated.", fs=4.6, color=INK2, ha="left")
    D.text(4.0, 4.4, f"The eight words differ in a1; they are written back to the same eight (bank, row) pairs {lat} cycles after the issue.", fs=4.6, color=INK2, ha="left")

    # ============================== (d) expressions and constants ======================================================
    card(D, 84, 2, 72, 40, "(d)", "Expressions and constants", "left: copied from rtl/fft512_mem.v by the figure script; right: values of this build", plot=False, fs=7.6)
    for r, pat in enumerate(CORE_EXPRESSIONS):
        D.text(86.0, 31.6 - 2.3 * r, _expr(pat), fs=4.4, color=INK, ha="left", mono=True)
    rows = [f"LB = LAT_BFLY = {lb}", f"LAT = LB + 2 = {lat}  (issue to readable)", f"STALL = LAT - 8 = {stall}  (idle after pass 0)", "DW = 2 WL",
            f"banks: 8 × 64 × DW = 1024 WL bit", f"{F['core']['u_ram'].count} ram_1r1w, 1 r8_butterfly,", f"{F['n_cru']} cru_rotate (3 + {bf['u_tw'].count}), {F['n_dl']} delay_line",
            f"first issue → done: 3·64 + STALL + LAT - 1 = {F['budget']['total_cycles']}", f"(CI RTL simulation: done sampled in trace cycle {done})"]
    for r, t in enumerate(rows):
        D.text(122.0, 31.6 - 2.9 * r, t, fs=5.0, color=INK, ha="left")
    D.save(FIGS / "fig26_fft_core_blocks")

FIGURES = [fig24_index_space, fig25_pass_snapshots, fig26_fft_core_blocks, fig27_twiddle_atlas, fig28_accuracy]


def main():
    for fn in FIGURES:
        fn()
        print("wrote", fn.__name__)


if __name__ == "__main__":
    main()
