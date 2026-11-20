"""Block diagrams 1-3: system overview, memory-based FFT architecture, radix-8 butterfly network.
Generated from code (matplotlib).  The paper's own figures are kept separately in docs/paper_figures."""
import matplotlib.pyplot as plt
import numpy as np

from _common import FIGS
from btpsd.schedule import bank_of, pass_addresses, row_of
from btpsd.signals import qpc
from btpsd.viz import INK, INK2, MUTED, PAL, SEQ_BLUE, Diagram, tint

LEGEND = [("io", "I/O stream"), ("memory", "memory"), ("compute", "arithmetic"), ("rom", "constant ROM"),
          ("control", "control / addressing"), ("accum", "accumulation")]


def _cmap():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("seqblue", ["#ffffff"] + SEQ_BLUE)


def _timeline(ax):
    """Overlapping chunks and the 1024-sample segments built from neighbouring chunk pairs."""
    ax.axis("off")
    ax.set_xlim(-0.1, 6.1)
    ax.set_ylim(-2.7, 2.9)
    for i in range(6):
        ax.add_patch(plt.Rectangle((i, 0), 1, 1, fc=tint(PAL["blue"], 0.84), ec=PAL["blue"], lw=0.9))
        ax.text(i + 0.5, 0.5, f"c{i}", ha="center", va="center", fontsize=6.6, color=INK, fontweight="bold")
    for k in range(5):
        up = k % 2 == 0
        y0, y1 = (1.15, 1.15 + 0.5 + 0.12 * (k // 2)) if up else (-0.15, -0.65 - 0.12 * (k // 2))
        ax.plot([k, k, k + 2, k + 2], [y0, y1, y1, y0], color=PAL["orange"] if up else PAL["violet"], lw=1.0)
        ax.text(k + 1, y1 + (0.32 if up else -0.34), f"S{k}", ha="center", va="center", fontsize=6.4, color=INK, fontweight="bold")
    ax.text(3.0, -2.35, "segment S_k = (c_k, c_k+1), hop 512, 50 % overlap", ha="center", fontsize=5.9, color=INK2)


def fig01_system_overview():
    from btpsd.bispectrum import BispecConfig, bispectrum_fixed, bispectrum_hw_to_float

    d = Diagram(118, 120, "Datapath overview: chunk-merge PSD and bispectrum estimator",
                subtitle="One memory-based radix-8 FFT serves both modes; WL, window, FD taps and ROM tolerance are parameters.")
    d.legend(3, 110.2, [e for e in LEGEND if e[0] != "rom"], dx=21)
    w, h = 17, 12
    x = [6 + i * (w + 5.2) for i in range(5)]
    y1, y2, y3 = 88, 68, 48
    b1 = d.box(x[0], y1, w, h, "Input stream\nx[n], WL bits\n512-sample chunks", "io", step=1)
    b2 = d.box(x[1], y1, w, h, "Chunk memory\nM1..M8\n8 banks x 64 words", "memory", step=2)
    b3 = d.box(x[2], y1, w, h, "Radix-8 memory FFT\n512 points, in place\n3 passes x 64 bfly", "compute", step=3)
    b4 = d.box(x[3], y1, w, h, "Spectrum read-out\ndigit-reversed address\nnew chunk A_s+1", "control", step=4)
    b5 = d.box(x[4], y1, w, h, "M9\nspectrum of the\nprevious chunk A_s", "memory", step=5)
    for a, b, bw in ((b1, b2, "WL"), (b2, b3, "2WL"), (b3, b4, "2WL"), (b4, b5, "2WL")):
        d.arrow([a.r, b.l], None, bus=bw)
    d.text((b4.r[0] + b5.l[0]) / 2, y1 + h / 2 + 3.7, "copy", fs=6.2)
    c6 = d.box(x[4], y2, w, h, "Segment merge\neven: (A+B)/2\nodd: FD6(A-B)/2", "compute", step=6)
    c7 = d.box(x[3], y2, w, h, "Window filter\n1 / 3 / 5 taps\non X(k)", "compute", step=7)
    c8 = d.box(x[2], y2, w, h, "Squarer\nsum of odd numbers\nre^2 + im^2", "compute", step=8)
    c9 = d.box(x[1], y2, w, h, "Save / Add\naccumulator\nover S segments", "accum", step=9)
    c10 = d.box(x[0], y2, w, h, "PSD P(k)\n>> avg shift\n1024 bins (symmetric)", "io", step=10)
    ym = (y1 + y2 + h) / 2
    d.arrow([(b5.b[0] + 2.5, y1), (c6.t[0] + 2.5, y2 + h)], None)
    d.text(b5.b[0] + 4.0, ym, "A_s(m)", fs=6.4, ha="left")
    d.arrow([(b4.b[0], y1), (b4.b[0], ym), (c6.t[0] - 2.5, ym), (c6.t[0] - 2.5, y2 + h)], None)
    d.text(b4.b[0] + 1.0, ym + 1.8, "A_s+1(m)", fs=6.4, ha="left")
    for a, b, lab, bw in ((c6, c7, "X(k)", "2WL"), (c7, c8, "Xw(k)", "2WL"), (c8, c9, "P", "2WL"), (c9, c10, "sum", "")):
        d.arrow([a.l, b.r], None, bus=bw or None)
        d.text((a.l[0] + b.r[0]) / 2, y2 + h / 2 - 2.2, lab, fs=6.4)
    e1 = d.box(x[0], y3, w, h, "Spectrum memory\nwindowed F(f)\n(banks M1..M8)", "memory", step=11)
    e2 = d.box(x[1], y3, w, h, "Triple product\nF(f1) F(f2) F*(f1+f2)\n2 complex mult", "compute", step=12)
    e3 = d.box(x[2], y3, w, h, "Accumulator\nover segments\n(wide register)", "accum", step=13)
    e4 = d.box(x[3], y3, w, h, "Bispectrum\nB(f1, f2)\nf1, f2 < N/4", "io", step=14)
    e5 = d.box(x[4], y3, w, h, "Bicoherence\n|B|^2/(P1 P2 P3)\n(optional)", "compute", step=15, dashed=True)
    for a, b in ((e1, e2), (e2, e3), (e3, e4)):
        d.arrow([a.r, b.l], None)
    d.arrow([e4.r, e5.l], None, ls="--")
    d.arrow([b2.t, (b2.t[0], 104.2), (1.8, 104.2), (1.8, e1.l[1]), e1.l], None)
    d.text(3.5, 106.0, "FFT result of the segment (bispectrum mode)", fs=6.6, ha="left")
    # ---- bottom strip: example results of the model -------------------------------------------------------
    d.text(59, 43.6, "Model examples on synthetic signals (not measured hardware)", fs=7.2, color=INK2, bold=True)
    d.card(2, 2.5, 38, 37)
    d.card(42, 2.5, 38, 37)
    d.card(82, 2.5, 34, 37)
    d.panel(4, 36.0, "(a)", "Chunks and segments")
    d.panel(44, 36.0, "(b)", "PSD, WL = 12")
    d.panel(84, 36.0, "(c)", "|B(f1, f2)| in dB, WL = 14")
    _timeline(d.inset(4.5, 5.0, 33, 26, grid=False))
    from _figdata import psd_example
    xs, ref, hw, snr = psd_example(wl=12)
    f = np.arange(513) * 256.0 / 1024
    ax = d.inset(50, 9.0, 27, 22)
    ax.plot(f, 10 * np.log10(ref[:513] + 1e-18), color=MUTED, lw=1.5, label="Welch, float")
    ax.plot(f, 10 * np.log10(hw[:513] + 1e-18), color=PAL["blue"], lw=0.9, label="model WL=12")
    ax.set_xlabel("frequency (Hz)", fontsize=6.0, labelpad=1)
    ax.set_ylabel("power (dB)", fontsize=6.0, labelpad=1)
    ax.legend(fontsize=5.4, frameon=False, loc="upper right", handlelength=1.4, borderpad=0.2)
    ax.text(0.04, 0.05, f"SQNR {snr:.1f} dB", transform=ax.transAxes, fontsize=5.8, color=INK)
    x, _ = qpc(512 * 17, coupled=True)
    re_, im_, S, sh = bispectrum_fixed(x, BispecConfig(wl=14, window="hann"))
    B = np.abs(bispectrum_hw_to_float(re_, im_, 14, S, sh))
    Bdb = 20 * np.log10(B / B.max() + 1e-12)
    ax = d.inset(90, 9.0, 22, 22, grid=False)
    ax.imshow(np.maximum(Bdb, -50).T, origin="lower", extent=[0, 128, 0, 128], cmap=_cmap(), aspect="equal",
              interpolation="nearest", vmin=-50, vmax=0)
    ax.set_xlabel("f1 (bin)", fontsize=6.0, labelpad=1)
    ax.set_ylabel("f2 (bin)", fontsize=6.0, labelpad=1)
    d.save(FIGS / "fig01_system_overview")


def fig02_memory_fft_architecture():
    d = Diagram(122, 93, "Memory-based radix-8 FFT: read, compute, write-back loop",
                subtitle="Eight single-port banks, one butterfly (8 words) per cycle; cyclic rotators keep every access conflict-free.")
    d.legend(3, 82.0, LEGEND, dx=19.5)
    d.box(26, 66, 54, 10, "Address generator (pass p, butterfly b)\nhi = b>>3, lo = b&7, s = (hi+lo) mod 8\nrow_j for j = 0..7, twiddle q_k for k = 1..7", "control", fs=6.8)
    d.box(88, 66, 32, 10, "Load (natural order)\none word per cycle\nbank, row from address a", "io", fs=6.6)
    wr = d.box(14, 51, 66, 8, "Write rotator\nbank j <- slot (j - s) mod 8   (s delayed by the pipeline latency)", "control", fs=6.8)
    bx0, bw, bg = 14.5, 7.5, 0.8
    by, bh = 26, 20
    d.rect(bx0 - 1.6, by - 1.6, 8 * bw + 7 * bg + 3.2, bh + 3.2)
    banks = [d.box(bx0 + j * (bw + bg), by, bw, bh, f"M{j + 1}\n ", "memory", fs=6.4) for j in range(8)]
    rd = d.box(14, 11, 66, 8, "Read rotator\nslot m <- bank (m + s) mod 8", "control", fs=6.8)
    hh = d._lh(6.4) + 2.0
    body = bh - hh
    for bk in banks:
        for r in range(0, 64, 4):
            yy = by + body - body * r / 64.0
            d.ax.plot([bk.x + 0.7, bk.x + bw - 0.7], [yy, yy], color=tint(PAL["blue"], 0.7), lw=0.35, zorder=2.3)
    b_ex = 19
    adr = pass_addresses(0)[b_ex]
    cells = sorted(zip(bank_of(adr).tolist(), row_of(adr).tolist()))
    for bnk, row in cells:
        yy = by + body - body * row / 64.0
        d.ax.add_patch(plt.Rectangle((banks[bnk].x + 0.7, yy - 0.45), bw - 1.4, 0.9, fc=PAL["orange"], ec="none", zorder=2.4))
    for j, bk in enumerate(banks):
        cx = bk.c[0]
        d.arrow([(cx, 51), (cx, by + bh + 1.7)], None, lw=0.8)
        d.arrow([(cx, by - 1.7), (cx, 19.0)], None, lw=0.8)
    d.text(83.0, 49.2, "write: 8 words", fs=6.2, ha="left")
    d.text(83.0, 21.6, "read: 8 words", fs=6.2, ha="left")
    bf = d.box(88, 9, 32, 12, "Radix-8 butterfly\n8-pt DFT, 3 radix-2 stages\n>>1 each; W8 via CRU", "compute", fs=6.6)
    d.box(88, 30, 32, 12, "7 x CRU twiddle rotators\nx W512^(j k) = e^(-j q pi/512)\nq_k from the address generator", "rom", fs=6.6)
    d.arrow([rd.r, bf.l], None, bus="8 x 2WL")
    d.arrow([(104, 21), (104, 30)], None)
    d.text(106.0, 25.5, "Y[k]", fs=6.4, ha="left")
    d.arrow([(104, 42), (104, 55), (80, 55)], None)
    d.text(92.0, 56.3, "rotated words", fs=6.2)
    d.arrow([(104, 66), (104, 62), (47, 62), (47, 59)], None)
    d.text(76.0, 63.3, "load path (write mux)", fs=6.2)
    d.box(0.8, 28, 8.4, 16, "Unload\nf -> (a2,a1,a0)\n= (f&7,\n(f>>3)&7,\nf>>6)", "io", fs=5.8)
    d.arrow([(bx0 - 1.6, 36), (9.2, 36)], None)
    d.text(5.0, 25.0, "to merge\nunit and M9", fs=5.8)
    bus = "#6250d6"
    d.arrow([(26, 71), (11.2, 71), (11.2, 15)], None, color=bus, ls="--", head=False, lw=0.9)
    for yy, lab in ((55, "s"), (42, "row_j"), (15, "s")):
        d.arrow([(11.2, yy), (14.0 if yy != 42 else bx0 - 1.6, yy)], None, color=bus, ls="--", lw=0.9)
        d.text(12.9, yy + 1.7, lab, fs=5.8, color=bus)
    cs = ", ".join(f"M{bnk + 1}:{row}" for bnk, row in cells)
    d.text(61, 6.2, "Address map:  a = 64 a2 + 8 a1 + a0,   bank(a) = (a2 + a1 + a0) mod 8,   row(a) = 8 a2 + a1   (a2, a1, a0 in 0..7)", fs=6.6, color=INK2)
    d.text(61, 2.6, f"Orange cells: the 8 words read by butterfly b = {b_ex} in pass 0, one per bank (bank:row  {cs}).  64 rows x 2 WL bits per bank.",
           fs=6.4, color=INK2)
    d.save(FIGS / "fig02_memory_fft_architecture")


def fig03_r8_butterfly_network():
    d = Diagram(118, 80, "Radix-8 butterfly: 8-point DFT as three radix-2 DIF stages")
    ys = [64 - 6.6 * i for i in range(8)]
    xs = [16, 42, 68]
    names_in = [f"x{i}" for i in range(8)]
    for i, yv in enumerate(ys):
        d.text(8.5, yv, names_in[i], fs=7.5, color=INK, ha="right", bold=True)
        d.arrow([(9.5, yv), (xs[0], yv)], None, head=False, lw=0.9)
    mult_lab = {0: {4: "1", 5: "W8^1", 6: "-j", 7: "W8^3"}, 1: {2: "1", 3: "-j", 6: "1", 7: "-j"}, 2: {}}
    dist = [4, 2, 1]
    row_state = list(range(8))  # row labels per stage for bookkeeping
    for s in range(3):
        xs_s = xs[s]
        xe = xs_s + 9
        pairs = []
        if s == 0:
            pairs = [(i, i + 4) for i in range(4)]
        elif s == 1:
            pairs = [(0, 2), (1, 3), (4, 6), (5, 7)]
        else:
            pairs = [(0, 1), (2, 3), (4, 5), (6, 7)]
        for (i, j) in pairs:
            yi, yj = ys[i], ys[j]
            for (a, b) in ((yi, yi), (yj, yj), (yi, yj), (yj, yi)):
                d.arrow([(xs_s, a), (xe, b)], None, head=False, lw=0.9, color=INK2)
            d.adder(xe, yi, "+", r=1.0)
            d.adder(xe, yj, "-", r=1.0)
            lab = mult_lab[s].get(j)
            # continuation wires
            nxt = xs_s + 26 if s < 2 else xs_s + 17
            d.arrow([(xe + 1.0, yi), (nxt, yi)], None, head=False, lw=0.9)
            if lab and lab != "1":
                d.box(xe + 4.2, yj - 1.5, 7.2, 3.0, lab.replace("W8^1", "W8¹").replace("W8^3", "W8³"), "rom", fs=6.2, round_=0.6)
                d.arrow([(xe + 1.0, yj), (xe + 4.2, yj)], None, head=False, lw=0.9)
                d.arrow([(xe + 11.4, yj), (nxt, yj)], None, head=False, lw=0.9)
            else:
                d.arrow([(xe + 1.0, yj), (nxt, yj)], None, head=False, lw=0.9)
    # outputs: stage-3 sums/diffs at xs[2]+9 ; wires continue to xs[2]+17
    out_order = [0, 4, 2, 6, 1, 5, 3, 7]
    xo = xs[2] + 17
    for r, yv in enumerate(ys):
        k = out_order[r]
        d.text(xo + 2.0, yv, f"Y[{k}]", fs=7.2, color=INK, ha="left", bold=True)
        if k == 0:
            d.arrow([(xo + 10.5, yv), (xo + 20.5, yv)], None, lw=0.9, ls="--")
            d.text(xo + 15.5, yv + 1.3, "x 1", fs=6.0)
        else:
            d.arrow([(xo + 10.5, yv), (xo + 13.5, yv)], None, head=False, lw=0.9)
            d.box(xo + 13.5, yv - 1.6, 8.2, 3.2, f"CRU q_{k}", "rom", fs=6.0, round_=0.6)
    d.text(xs[0] + 4.5, 71.2, "stage 1 (d = 4)", fs=7, color=INK)
    d.text(xs[1] + 4.5, 71.2, "stage 2 (d = 2)", fs=7, color=INK)
    d.text(xs[2] + 4.5, 71.2, "stage 3 (d = 1)", fs=7, color=INK)
    d.text(xo + 17.6, 71.2, "twiddle", fs=7, color=INK)
    d.text(59, 8.4, "Every add/subtract output is shifted right by one bit (rounded): a radix-8 butterfly scales by 1/8, the 512-point FFT by 1/512.\n"
           "W8^1 = e^(-j pi/4) and W8^3 are CRU rotations (ROM entry r = 128); -j is a free swap/negate.  Y[k] is the natural-order DFT output; the\n"
           "twiddle exponent of pass 0 is q = 2 b k, of pass 1 q = 16 lo k, pass 2 has none (units of pi/512).", fs=6.4, color=INK2)
    d.save(FIGS / "fig03_r8_butterfly_network")


if __name__ == "__main__":
    fig01_system_overview()
    fig02_memory_fft_architecture()
    fig03_r8_butterfly_network()
    print("diagrams 1-3 written")
