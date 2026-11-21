"""RTL-level figure sheets (15-18): module hierarchy, control state machines, pipeline and interface timing.

Everything here is derived from the Verilog sources (``btpsd.rtlparse``) and from the constants shared by model and RTL
(``btpsd.rtlmodel``); ``python/tests/test_rtl_structure.py`` checks the drawn hierarchy and the drawn transitions against
the sources.  Figures built from simulation output (cycle budget, issue trace) are in ``make_sim_figures.py``.
"""
from __future__ import annotations

import re

from matplotlib.patches import FancyBboxPatch, Rectangle

from _common import FIGS
from btpsd.rtlmodel import LAT, LAT_BFLY, LAT_CRU, stall_cycles, write_back_slack
from btpsd.rtlparse import RTL_DIR, fsm_edges, parse_all, reg_widths
from btpsd.viz import GRID, INK, INK2, KIND, MUTED, PAL, SURF, Diagram, tint

LEGEND = [("io", "top level / interface"), ("control", "control / delay"), ("compute", "arithmetic"),
          ("memory", "memory"), ("rom", "generated ROM"), ("accum", "accumulation")]


# ---------------------------------------------------------------------------------------------------------------
# drawing helpers
# ---------------------------------------------------------------------------------------------------------------
def _colors(kind):
    col = KIND.get(kind)
    edge = PAL[col] if col else INK2
    band = tint(PAL[col], 0.78) if col else "#e7e6e1"
    light = tint(PAL[col], 0.93) if col else "#f8f7f5"
    return edge, band, light


def card_box(d, x, y, w, h, kind, title, right=None, hh=4.6, lw=1.0, fs=7.8):
    """Rounded card with a coloured header band, the title at the left and an optional note at the right."""
    edge, band, light = _colors(kind)
    style = "round,pad=0,rounding_size=1.0"
    d.ax.add_patch(FancyBboxPatch((x + 0.3, y - 0.3), w, h, boxstyle=style, fc="#000000", ec="none", alpha=0.08, zorder=1))
    body = FancyBboxPatch((x, y), w, h, boxstyle=style, fc=light, ec="none", zorder=2)
    d.ax.add_patch(body)
    r = Rectangle((x, y + h - hh), w, hh, fc=band, ec="none", zorder=2.1)
    d.ax.add_patch(r)
    r.set_clip_path(body)
    d.ax.plot([x, x + w], [y + h - hh] * 2, color=edge, lw=0.6, alpha=0.55, zorder=2.15, solid_capstyle="butt")
    d.ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=style, fc="none", ec=edge, lw=lw, zorder=2.2))
    d.ax.text(x + 1.6, y + h - hh / 2, title, fontsize=fs, fontweight="bold", color=INK, ha="left", va="center",
              family="DejaVu Sans Mono", zorder=3)
    if right:
        d.ax.text(x + w - 1.4, y + h - hh / 2, right, fontsize=5.6, color=INK2, ha="right", va="center", zorder=3)


def chip(d, x, y, text, kind, fs=5.7, h=3.5):
    """Small rounded label (a leaf module or a state); returns its right edge."""
    w = 0.60 * len(text) * fs / 6.0 + 3.0
    edge, band, _ = _colors(kind)
    d.ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.7", fc=band, ec=edge, lw=0.8,
                                  zorder=3))
    d.ax.text(x + w / 2, y + h / 2 - 0.05, text, fontsize=fs, color=INK, ha="center", va="center", zorder=4)
    return x + w


def lines(d, x, y, rows, fs=6.2, color=INK2, dy=1.75, bold_first=False):
    for i, t in enumerate(rows):
        d.ax.text(x, y - i * dy, t, fontsize=fs, color=INK if (bold_first and i == 0) else color, ha="left", va="center",
                  zorder=3)


def _expr(s, env):
    return int(eval(s, {"__builtins__": {}}, dict(env)))  # noqa: S307 - expressions come from our own sources


def _kb(bits):
    return f"{bits / 1024:.1f} kbit" if bits >= 1024 else f"{bits} bit"


def _fmt(s):
    return s.replace("*", "·").replace(" ", "")


# ---------------------------------------------------------------------------------------------------------------
# fig15: module hierarchy
# ---------------------------------------------------------------------------------------------------------------
def storage_bits(mods, top, env):
    """Bits of ``ram_1r1w`` storage below ``top`` for the default parameters ``env`` (instances multiplied through)."""
    def walk(name, mult, e):
        m = mods[name]
        local = dict(e)
        if name == "fft512_mem":
            local["DW"] = 2 * local["WL"]
        total = 0
        for i in m.instances:
            if i.module == "ram_1r1w":
                words = 2 ** _expr(i.params["AW"], local)
                total += mult * i.count * words * _expr(i.params["DW"], local)
            else:
                total += walk(i.module, mult * i.count, local)
        return total
    return walk(top, 1, env)


def count_below(mods, top, leaf):
    def walk(name):
        n = 0
        for i in mods[name].instances:
            n += i.count * (1 if i.module == leaf else 0) + i.count * walk(i.module)
        return n
    return walk(top)


def fig15_hierarchy():
    mods = parse_all()
    n_lines = sum(m.lines for m in mods.values())
    D = Diagram(142, 114, "RTL module hierarchy and interfaces",
                subtitle=f"Read from the Verilog sources by btpsd.rtlparse: {len(mods)} modules, {n_lines} lines; counts are instances.")
    D.legend(3, 101.6, LEGEND, dx=23)

    def chips_for(mod, y, x0, kinds):
        x = x0
        for inst_mod, label, kind in kinds:
            n = sum(i.count for i in mods[mod].instances if i.module == inst_mod)
            x = chip(D, x, y, label.format(n=n), kind) + 1.4

    def top_text(x, y, h, rows):
        lines(D, x, y + h - 7.2, rows, bold_first=True)

    # ---- top levels ------------------------------------------------------------------------------------------
    t1 = (4, 78, 62, 19)
    t2 = (76, 78, 62, 19)
    card_box(D, *t1, "io", "bt_psd_top", f"rtl/bt_psd_top.v, {mods['bt_psd_top'].lines} lines")
    top_text(5.6, 78, 19, ["chunk-merge PSD estimator", "sample stream in, 1024-bin power spectrum out",
                           f"{len(mods['bt_psd_top'].states)}-state controller T_INIT .. T_SEG; two FFT cores form a ping-pong pair"])
    lines(D, 5.6, 84.0, ["in_valid, in_data[WL], in_ready, busy, avg_shift  |  res_addr[9:0] \u2192 res_data[ACCW]"], fs=5.4)
    card_box(D, *t2, "io", "bt_bispectrum_top", f"rtl/bt_bispectrum_top.v, {mods['bt_bispectrum_top'].lines} lines")
    top_text(77.6, 78, 19, ["direct bispectrum estimator", "sample stream in, B[f1][f2] with f1, f2 < 2^LOGB out",
                            f"{len(mods['bt_bispectrum_top'].states)}-state controller B_INIT .. B_ENG, hop-256 framing"])
    lines(D, 77.6, 84.0, ["in_valid, in_data[WL], in_ready, busy, avg_shift  |  res_idx \u2192 res_re, res_im"], fs=5.4)
    chips_for("bt_bispectrum_top", 79.6, 77.6, [("ram_1r1w", "ram_1r1w \u00d7{n}  half-frame buffer 256 \u00d7 WL", "memory")])

    # ---- engines and the shared FFT core -----------------------------------------------------------------------
    y0, h0 = 44, 27
    e1, c, e2 = (4, y0, 41, h0), (50, y0, 42, h0), (97, y0, 41, h0)
    card_box(D, *e1, "compute", "psd_engine", f"{mods['psd_engine'].lines} lines")
    lines(D, 5.6, y0 + h0 - 7.2, ["merge (6-tap FD), window,", "squarer, accumulation",
                                  f"{len(mods['psd_engine'].states)} states, one bin at a time"])
    chips_for("psd_engine", y0 + 11.6, 5.6, [("ram_1r1w", "ram_1r1w \u00d7{n}", "memory")])
    lines(D, 5.6, y0 + 9.6, ["segment 1024 \u00d7 2WL,", "accumulator 1024 \u00d7 ACCW"], fs=5.6, dy=1.45)
    chips_for("psd_engine", y0 + 1.2, 5.6, [("mag_sq_odd", "mag_sq_odd \u00d7{n}  sum-of-odds squarer", "compute")])

    card_box(D, *c, "compute", "fft512_mem", f"{mods['fft512_mem'].lines} lines")
    lines(D, 51.6, y0 + h0 - 7.2, ["512-point radix-8 FFT, in place,", "8 banks, one butterfly per cycle,",
                                   f"{len(mods['fft512_mem'].states)}-state issue controller"])
    chips_for("fft512_mem", y0 + 11.6, 51.6, [("ram_1r1w", "ram_1r1w \u00d7{n}", "memory"), ("delay_line", "delay_line \u00d7{n}", "control")])
    lines(D, 51.6, y0 + 9.6, ["banks 64 \u00d7 2WL each; control pipeline", "9 bit wide, 1 and LAT_BFLY stages"], fs=5.6, dy=1.45)

    card_box(D, *e2, "compute", "bispectrum_engine", f"{mods['bispectrum_engine'].lines} lines")
    lines(D, 98.6, y0 + h0 - 7.2, ["window, rounded product P12,", "P12 \u00b7 conj(F3), accumulation",
                                   f"{len(mods['bispectrum_engine'].states)} states, 4 cycles per (f1, f2)"])
    chips_for("bispectrum_engine", y0 + 11.6, 98.6, [("ram_1r1w", "ram_1r1w \u00d7{n}", "memory")])
    lines(D, 98.6, y0 + 9.6, ["2 \u00d7 windowed spectrum 512 \u00d7 2WL,", "2 \u00d7 accumulator 2^(2 LOGB) \u00d7 BACCW (re, im)"],
          fs=5.6, dy=1.45)

    # ---- butterfly and rotator ------------------------------------------------------------------------------------
    bf, cr = (50, 20, 42, 16), (50, 1.0, 42, 15)
    card_box(D, *bf, "compute", "r8_butterfly", f"{mods['r8_butterfly'].lines} lines")
    lines(D, 51.6, 20 + 16 - 7.2, ["scaled 8-point DFT (3 radix-2 stages,", f">>1 each) + output twiddles, latency {LAT_BFLY}"])
    chips_for("r8_butterfly", 21.6, 51.6, [("delay_line", "delay_line \u00d7{n}", "control")])
    card_box(D, *cr, "compute", "cru_rotate", f"{mods['cru_rotate'].lines} lines")
    lines(D, 51.6, 1.0 + 15 - 7.2, [f"y = x \u00b7 exp(-j q \u03c0/512), {LAT_CRU} register stages"])
    chips_for("cru_rotate", 3.0, 51.6, [("cru_rom", "cru_rom \u00d7{n}  generated, 129 entries", "rom")])

    # ---- edges with multiplicities -----------------------------------------------------------------------------------
    n_core = sum(i.count for i in mods["bt_psd_top"].instances if i.module == "fft512_mem")
    D.arrow([(24, 78), (24, 71)], None)
    D.text(25.2, 74.2, "\u00d71   go, clr \u2192 busy, done;  a/b read ports, result port", fs=5.4, ha="left")
    D.arrow([(56, 78), (56, 71)], None)
    D.text(54.6, 74.2, f"\u00d7{n_core}", fs=6.4, ha="right", color=INK, bold=True)
    D.arrow([(86, 78), (86, 71)], None)
    D.text(87.2, 74.2, "\u00d71", fs=6.4, ha="left", color=INK, bold=True)
    D.arrow([(118, 78), (118, 71)], None)
    D.text(119.2, 74.2, "\u00d71   go, clr \u2192 busy, done;  f read port", fs=5.4, ha="left")
    D.text(71.0, 76.4, "ld_en, ld_addr, ld_re/im  \u00b7  start \u2192 busy, done  \u00b7  rd_k \u2192 rd_re/im", fs=5.4, color=INK2)
    D.arrow([(71, 44), (71, 36)], None)
    D.text(72.2, 40.0, "\u00d71", fs=6.4, ha="left", color=INK, bold=True)
    D.arrow([(71, 20), (71, 16)], None)
    n_cru = count_below(mods, "r8_butterfly", "cru_rotate")
    D.text(69.6, 18.0, f"\u00d7{n_cru}  (3 \u00d7 W8, 7 \u00d7 twiddle)", fs=5.8, ha="right", color=INK, bold=True)

    # ---- storage table (computed from the parameters) -----------------------------------------------------------------
    env = {"WL": 12, "ACCW": 2 * 12 + 8, "BACCW": 2 * 12 + 10, "LOGB": 5}
    sp, sb = storage_bits(mods, "bt_psd_top", env), storage_bits(mods, "bt_bispectrum_top", env)
    D.card(4, 1.5, 41, 38)
    D.text(5.8, 36.0, "Instances and storage", fs=7.2, ha="left", color=INK, bold=True)
    D.text(5.8, 33.2, "WL = 12, ROUND = 1, LOGB = 5; RAM bits only", fs=5.2, ha="left")
    rows = [("", "bt_psd_top", "bt_bispectrum_top"),
            ("fft512_mem", str(count_below(mods, "bt_psd_top", "fft512_mem")), str(count_below(mods, "bt_bispectrum_top", "fft512_mem"))),
            ("cru_rotate", str(count_below(mods, "bt_psd_top", "cru_rotate")), str(count_below(mods, "bt_bispectrum_top", "cru_rotate"))),
            ("ram_1r1w", str(count_below(mods, "bt_psd_top", "ram_1r1w")), str(count_below(mods, "bt_bispectrum_top", "ram_1r1w"))),
            ("RAM bits", _kb(sp), _kb(sb))]
    D.ax.plot([5.8, 43], [30.4, 30.4], color=MUTED, lw=0.5, zorder=3)
    for i, (a, b, c2) in enumerate(rows):
        y = 28.2 - i * 3.6
        bold = i == 0
        D.text(5.8, y, a, fs=5.9, ha="left", color=INK, bold=bold)
        D.text(28.0, y, b, fs=5.9, ha="right", color=INK, bold=bold)
        D.text(43.0, y, c2, fs=5.9, ha="right", color=INK, bold=bold)
    D.text(5.8, 8.6, "Counts include every instance below the top.", fs=5.0, ha="left")
    D.text(5.8, 6.2, "RAM bits: sum over ram_1r1w of 2^AW \u00d7 DW.", fs=5.0, ha="left")

    D.card(97, 1.5, 41, 38)
    D.text(98.8, 36.0, "Notes", fs=7.2, ha="left", color=INK, bold=True)
    notes = ["Generated, do not edit: cru_rom.v,", "btpsd_consts.vh (python/scripts/gen_rtl.py).", "",
             "Not synthesised: no LUT, DSP, BRAM, timing", "or power figure exists for these sources.", "",
             "The tops are sequential references and differ", "from the paper's overlapped design",
             "(docs/ARCHITECTURE.md, section 6)."]
    lines(D, 98.8, 32.4, notes, fs=5.5, dy=2.35)
    D.save(FIGS / "fig15_rtl_hierarchy")


# ---------------------------------------------------------------------------------------------------------------
# finite-state-machine drawing (figs 16, 17)
# ---------------------------------------------------------------------------------------------------------------
def _poly_point(pts, f):
    """Point at fraction ``f`` of the length of the polyline ``pts``."""
    seg = [((x1 - x0) ** 2 + (y1 - y0) ** 2) ** 0.5 for (x0, y0), (x1, y1) in zip(pts, pts[1:])]
    target, acc = f * sum(seg), 0.0
    for (x0, y0), (x1, y1), L in zip(pts, pts[1:], seg):
        if acc + L >= target and L > 0:
            t = (target - acc) / L
            return x0 + t * (x1 - x0), y0 + t * (y1 - y0)
        acc += L
    return pts[-1]


DRAWN: dict[tuple[str, str], set] = {}  # (rtl file, state variable) -> drawn edges, filled by FsmDrawing.check


class FsmDrawing:
    """Nodes on a grid, edges as right-angle polylines; ``check`` compares the drawn edges with the Verilog."""

    def __init__(self, D, rtl_file):
        self.D, self.rtl_file, self.nodes, self.drawn = D, rtl_file, {}, {}

    def node(self, name, x, y, w, h, kind, text, fs=7.0):
        self.nodes[name] = self.D.box(x, y, w, h, f"{name}\n{text}", kind, fs=fs)
        return self.nodes[name]

    def edge(self, a, b, pts, label=None, lpos=0.5, loff=(0.0, 1.1), dashed=False, ha="center", fs=5.4):
        self.D.arrow(pts, None, ls="--" if dashed else "-")
        if label:
            x, y = _poly_point(pts, lpos)
            self.D.text(x + loff[0], y + loff[1], label, fs=fs, ha=ha, color=INK)
        self.drawn[(a, b)] = pts

    # straight connections between boxes -------------------------------------------------------------------------------
    def right(self, a, b, label=None, **kw):
        A, B = self.nodes[a], self.nodes[b]
        self.edge(a, b, [A.r, B.l], label, **kw)

    def left(self, a, b, label=None, **kw):
        A, B = self.nodes[a], self.nodes[b]
        self.edge(a, b, [A.l, B.r], label, **kw)

    def down(self, a, b, label=None, dx=0.0, **kw):
        A, B = self.nodes[a], self.nodes[b]
        self.edge(a, b, [(A.b[0] + dx, A.b[1]), (B.t[0] + dx, B.t[1])], label, **kw)

    def up(self, a, b, label=None, dx=0.0, **kw):
        A, B = self.nodes[a], self.nodes[b]
        self.edge(a, b, [(A.t[0] + dx, A.t[1]), (B.b[0] + dx, B.b[1])], label, **kw)

    def over(self, a, b, h, label=None, xa=0.0, xb=0.0, lpos=0.5, **kw):
        """From the top of ``a`` up by ``h``, sideways, down onto the top of ``b``."""
        A, B = self.nodes[a], self.nodes[b]
        ytop = max(A.y + A.h, B.y + B.h) + h
        pts = [(A.t[0] + xa, A.y + A.h), (A.t[0] + xa, ytop), (B.t[0] + xb, ytop), (B.t[0] + xb, B.y + B.h)]
        self.edge(a, b, pts, label, lpos=lpos, dashed=kw.pop("dashed", True), **kw)

    def under(self, a, b, h, label=None, xa=0.0, xb=0.0, lpos=0.5, y=None, **kw):
        """From the bottom of ``a`` down by ``h`` (or to the absolute ``y``), sideways, up onto the bottom of ``b``."""
        A, B = self.nodes[a], self.nodes[b]
        ybot = y if y is not None else min(A.y, B.y) - h
        pts = [(A.b[0] + xa, A.y), (A.b[0] + xa, ybot), (B.b[0] + xb, ybot), (B.b[0] + xb, B.y)]
        self.edge(a, b, pts, label, lpos=lpos, dashed=kw.pop("dashed", True), **kw)

    def check(self, var="st"):
        """The drawn edges must be exactly the transitions found in the Verilog."""
        rtl = set(fsm_edges(RTL_DIR / self.rtl_file, var))
        drawn = set(self.drawn)
        DRAWN[(self.rtl_file, var)] = drawn
        if rtl != drawn:
            raise AssertionError(f"{self.rtl_file}: missing in drawing {sorted(rtl - drawn)}, not in RTL {sorted(drawn - rtl)}")


FSM_LEGEND = [("io", "idle / end"), ("control", "clear accumulators"), ("memory", "load / copy samples"),
              ("compute", "FFT core running"), ("accum", "engine running")]


def _panel_card(D, x, y, w, h, tag, title, sub):
    D.card(x, y, w, h)
    D.panel(x + 1.8, y + h - 3.3, tag, title, fs=8.2)
    D.text(x + 1.8, y + h - 6.6, sub, fs=5.8, ha="left")


def fig16_top_fsms():
    mods = parse_all()
    D = Diagram(150, 140, "Control state machines of the top levels and the FFT core",
                subtitle="Every drawn transition is a transition of the Verilog (checked at generation time and by python/tests/test_rtl_structure.py).")
    D.legend(3, 128.6, FSM_LEGEND, dx=29)

    # ---- (a) bt_psd_top ---------------------------------------------------------------------------------------
    _panel_card(D, 2, 90, 146, 35, "(a)", "bt_psd_top: five states",
                "T_LOAD takes the 512 samples of a chunk into core[cur]; the segment (older core, newer core) goes to the engine from the second chunk on")
    F = FsmDrawing(D, "bt_psd_top.v")
    y, h, w, pitch, x0 = 99.0, 11.0, 22.0, 29.0, 8.0
    xs = [x0 + i * pitch for i in range(5)]
    F.node("T_INIT", xs[0], y, w, h, "io", "eng_clr <= 1")
    F.node("T_CLR", xs[1], y, w, h, "control", "engine zeroes 1024 bins")
    F.node("T_LOAD", xs[2], y, w, h, "memory", "in_ready = 1, ld_cnt++")
    F.node("T_FFT", xs[3], y, w, h, "compute", "core[cur] transforms")
    F.node("T_SEG", xs[4], y, w, h, "accum", "psd_engine, one segment")
    F.right("T_INIT", "T_CLR", "next cycle")
    F.right("T_CLR", "T_LOAD", "eng_done")
    F.right("T_LOAD", "T_FFT", "512th sample", loff=(0, 1.3))
    F.right("T_FFT", "T_SEG", "done &\nhave_prev", loff=(0, 2.0))
    F.over("T_FFT", "T_LOAD", 5.0, "done & !have_prev:  have_prev <= 1, cur <= ~cur", lpos=0.5, loff=(0, 1.1), xa=-3.0, xb=3.0)
    F.under("T_SEG", "T_LOAD", 5.0, "eng_done:  cur <= ~cur", lpos=0.5, loff=(0, 1.3), xa=-3.0, xb=3.0)
    F.check()
    D.text(8, 91.8, "in_ready = (st == T_LOAD),   busy = (st != T_LOAD);  start pulse to core[cur] with the 512th sample;  eng_go pulse on entering T_SEG",
           fs=5.6, ha="left")

    # ---- (b) bt_bispectrum_top ------------------------------------------------------------------------------------
    _panel_card(D, 2, 48, 146, 40, "(b)", "bt_bispectrum_top: six states",
                "frames of 512 samples, hop 256: the second half of a frame is kept in a 256-word RAM and copied to the start of the next frame (B_COPY)")
    G = FsmDrawing(D, "bt_bispectrum_top.v")
    y, h, w, pitch, x0 = 57.0, 11.0, 17.0, 24.5, 5.0
    xs = [x0 + i * pitch for i in range(6)]
    G.node("B_INIT", xs[0], y, w, h, "io", "eng_clr <= 1")
    G.node("B_CLR", xs[1], y, w, h, "control", "engine zeroes bins")
    G.node("B_COPY", xs[2], y, w, h, "memory", "256 words to core")
    G.node("B_EXT", xs[3], y, w, h, "memory", "take new samples")
    G.node("B_FFT", xs[4], y, w, h, "compute", "core transforms")
    G.node("B_ENG", xs[5], y, w, h, "accum", "bispectrum_engine")
    G.right("B_INIT", "B_CLR")
    G.right("B_CLR", "B_COPY", "!first", loff=(0, 1.3))
    G.right("B_COPY", "B_EXT", "256 copied", loff=(0, 1.3))
    G.right("B_EXT", "B_FFT", "512th\nsample", loff=(0, 2.0))
    G.right("B_FFT", "B_ENG", "core_done", loff=(0, 1.3))
    G.over("B_CLR", "B_EXT", 5.0, "eng_done, first frame", lpos=0.5, loff=(0, 1.1), xa=3.0, xb=-3.0)
    G.under("B_ENG", "B_COPY", 5.0, "eng_done:  first <= 0", lpos=0.5, loff=(0, 1.3), xa=-3.0, xb=3.0)
    G.check()
    D.text(5, 49.8, "in_ready = (st == B_EXT),   busy = (st != B_EXT);  core_start pulse with the 512th sample;  eng_go pulse on leaving B_FFT",
           fs=5.6, ha="left")

    # ---- (c) fft512_mem ---------------------------------------------------------------------------------------------
    _panel_card(D, 2, 4, 146, 42, "(c)", "fft512_mem: three states",
                "the issue controller sequences 3 passes of 64 butterflies; a stall after pass 0 covers the write-back latency")
    H = FsmDrawing(D, "fft512_mem.v")
    iss = 3 * 64 + stall_cycles(LAT)
    from btpsd.rtlmodel import Fft512Rtl
    sched = Fft512Rtl().issue_schedule(stall_cycles(LAT))
    last = sched[-1][0]
    yc = 27.5
    H.node("ST_IDLE", 6, yc - 5.5, 22, 11, "io", "load / unload ports")
    run = D.box(40, yc - 11.0, 74, 22, "", "compute")
    D.ax.add_patch(FancyBboxPatch((40, yc - 11.0), 74, 22, boxstyle="round,pad=0,rounding_size=1.0", fc="none", ec=PAL["aqua"], lw=1.0, zorder=2.2))
    D.text(42, yc + 8.0, "ST_RUN", fs=7.2, ha="left", color=INK, bold=True)
    D.text(112, yc + 8.0, f"issue controller, {iss} issue slots", fs=5.8, ha="right")
    first = {p: min(t for t, pp, _ in sched if pp == p) for p in range(3)}
    lastp = {p: max(t for t, pp, _ in sched if pp == p) for p in range(3)}
    segs = [("pass 0", first[0], lastp[0]), ("stall", lastp[0] + 1, first[1] - 1), ("pass 1", first[1], lastp[1]), ("pass 2", first[2], lastp[2])]
    xa, span, ybar = 42.0, 70.0, yc - 6.0
    for name, a, b in segs:
        n = b - a + 1
        wseg = span * n / iss
        hatch = name == "stall"
        D.ax.add_patch(Rectangle((xa, ybar), wseg, 10.0, fc=tint(PAL["aqua"], 0.55) if not hatch else "#ffffff",
                                 ec=PAL["aqua"] if not hatch else PAL["orange"], lw=0.8, hatch="////" if hatch else None, zorder=3))
        if not hatch:
            D.text(xa + wseg / 2, ybar + 6.2, name, fs=6.4, color=INK, bold=True)
            D.text(xa + wseg / 2, ybar + 2.8, f"b = 0..63,  {n} cycles", fs=5.6)
        else:
            D.text(xa + wseg / 2 + 1.5, ybar + 13.0, f"STALL = {n}", fs=5.6, color=INK)
            D.arrow([(xa + wseg / 2 + 1.5, ybar + 11.6), (xa + wseg / 2, ybar + 10.0)], None, lw=0.7)
        D.text(xa, ybar - 1.6, ("slot " if a == 0 else "") + str(a), fs=5.0, color=INK2, ha="left" if a == 0 else "center")
        xa += wseg
    D.text(xa, ybar - 1.6, str(last), fs=5.0, color=INK2)
    H.nodes["ST_RUN"] = run
    H.node("ST_DRAIN", 121, yc - 5.5, 22, 11, "compute", f"LAT - 2 = {LAT - 2} cycles")
    H.right("ST_IDLE", "ST_RUN", "start", loff=(0, 1.3))
    H.right("ST_RUN", "ST_DRAIN", "last issue", loff=(0, 1.3))
    H.under("ST_DRAIN", "ST_IDLE", 0, "drain_cnt = 0:  busy <= 0, done <= 1 (one-cycle pulse)", lpos=0.5, loff=(0, 1.3), y=yc - 14.5)
    H.check()
    D.text(6, 6.4, f"first issue to done: {last} + LAT = {last + LAT} cycles  (LAT_BFLY = {LAT_BFLY}, LAT = LAT_BFLY + 2 = {LAT};  "
           f"tb_fft512_mem checks this number in simulation)", fs=5.6, ha="left")
    D.save(FIGS / "fig16_fsm_top_levels")


ENGINE_LEGEND = [("io", "idle / end"), ("control", "clear / wait"), ("memory", "RAM read request"),
                 ("compute", "multiply-accumulate"), ("accum", "round, saturate, write"), ("rom", "squarer handshake")]

_GX = [6.0, 34.5, 63.0, 91.5, 120.0]   # five grid columns of the engine drawings
_GW = 23.0


def _engine_grid(D, F, y_top, h=10.5, pitch_y=19.5):
    """Row y positions (bottom edge) of the three grid rows."""
    return [y_top - i * pitch_y for i in range(3)]


def fig17_engine_fsms():
    mods = parse_all()
    D = Diagram(150, 173, "Control state machines of the two engines",
                subtitle="Same grid for both: idle / clear / end in the first column; solid arrows follow the main flow, dashed arrows are loops.")
    D.legend(3, 161.6, ENGINE_LEGEND, dx=24)
    h = 10.5

    # ---- (a) psd_engine ------------------------------------------------------------------------------------------
    _panel_card(D, 2, 84, 146, 74, "(a)", "psd_engine: 14 states, one segment per go pulse",
                "merge of the two chunk spectra (m = 0..511, six taps), then window, square and accumulate (k = 0..1023, five taps)")
    F = FsmDrawing(D, "psd_engine.v")
    r1, r2, r3 = 128.0, 108.5, 89.0
    X = _GX
    F.node("S_IDLE", X[0], r1, _GW, h, "io", "wait for clr / go")
    F.node("S_CLR", X[0], r2, _GW, h, "control", "acc[cnt] <= 0, 1024 cycles")
    F.node("S_END", X[0], r3, _GW, h, "io", "busy <= 0, done <= 1")
    F.node("M_A", X[1], r1, _GW, h, "memory", "a_addr = b_addr =\n(m + t - 2) mod 512")
    F.node("M_B", X[2], r1, _GW, h, "control", "RAM read latency")
    F.node("M_C", X[3], r1, _GW, h, "compute", "facc += h[t] (-1)^idx D\nt = 2: write X(2m)")
    F.node("M_F", X[4], r1, _GW, h, "accum", "write X(2m+1)\nfacc <= 0")
    F.node("W_A", X[1], r2, _GW, h, "memory", "x_raddr =\n(k + 2 - t) mod 1024")
    F.node("W_B", X[2], r2, _GW, h, "control", "RAM read latency")
    F.node("W_C", X[3], r2, _GW, h, "compute", "wacc += W[t] X(k - d)")
    F.node("W_F", X[4], r2, _GW, h, "accum", "sq_re, sq_im <=\nrsat(wacc)")
    F.node("W_G", X[4], r3, _GW, h, "rom", "sq_start")
    F.node("W_S", X[3], r3, _GW, h, "rom", "wait for sq_done")
    F.node("W_ACC", X[2], r3, _GW, h, "accum", "acc[k] += p")
    # frame states
    F.down("S_IDLE", "S_CLR", "clr", loff=(1.2, 0), ha="left")
    F.down("S_CLR", "S_END", "cnt == 1023", loff=(1.2, 0), ha="left")
    lane = X[0] - 3.0
    A, B = F.nodes["S_END"], F.nodes["S_IDLE"]
    F.edge("S_END", "S_IDLE", [A.l, (lane, A.l[1]), (lane, B.l[1]), B.l], None, dashed=True)
    D.text(lane + 1.5, (A.l[1] + B.l[1]) / 2, "done pulse", fs=5.2, rot=90, color=INK)
    F.right("S_IDLE", "M_A", "go", loff=(0, 1.1))
    # merge row
    F.right("M_A", "M_B")
    F.right("M_B", "M_C")
    F.right("M_C", "M_F", "t == 5", loff=(0, 1.1))
    F.over("M_C", "M_A", 4.5, "t != 5:  t++", xa=-3.0, xb=4.0, loff=(0, 1.1))
    F.over("M_F", "M_A", 9.0, "m != 511:  m++", xa=0.0, xb=-4.0, loff=(0, 1.1))
    # window row
    A, B = F.nodes["M_F"], F.nodes["W_A"]
    yg = r1 - 4.0
    F.edge("M_F", "W_A", [A.b, (A.b[0], yg), (B.t[0] - 4.0, yg), (B.t[0] - 4.0, B.y + B.h)], "m == 511", lpos=0.62, loff=(0, 1.1))
    F.right("W_A", "W_B")
    F.right("W_B", "W_C")
    F.right("W_C", "W_F", "t == 4", loff=(0, 1.1))
    A, B = F.nodes["W_C"], F.nodes["W_A"]
    yl = r2 + h + 2.2
    F.edge("W_C", "W_A", [(A.t[0], A.y + A.h), (A.t[0], yl), (B.t[0] + 4.0, yl), (B.t[0] + 4.0, B.y + B.h)], "t != 4:  t++", lpos=0.5,
           loff=(0, -1.3), dashed=True)
    # square and accumulate row (right to left)
    F.down("W_F", "W_G")
    F.left("W_G", "W_S")
    F.left("W_S", "W_ACC", "sq_done", loff=(0, 1.1))
    A, B = F.nodes["W_ACC"], F.nodes["S_END"]
    F.edge("W_ACC", "S_END", [A.l, B.r], "k == 1023", loff=(0, 1.1))
    A, B = F.nodes["W_ACC"], F.nodes["W_A"]
    yk = r3 + h + 3.8
    F.edge("W_ACC", "W_A", [(A.t[0], A.y + A.h), (A.t[0], yk), (B.t[0], yk), (B.t[0], B.y)], "k != 1023:  k++", lpos=0.55,
           loff=(0, 1.1), dashed=True)
    F.check()
    D.text(6, 85.8, "merge: 512 x (6 x 3 + 1) = 9 728 cycles;  window: 1 024 x (5 x 3 + W_F + W_G + W_S + W_ACC) = 1 024 x (18 + squarer wait) cycles;"
           "  the squarer wait is data dependent (max(|re|, |im|) + 2 cycles)", fs=5.4, ha="left")

    # ---- (b) bispectrum_engine ----------------------------------------------------------------------------------------
    _panel_card(D, 2, 3, 146, 78, "(b)", "bispectrum_engine: 14 states, one frame per go pulse",
                "window of the bins 0 .. 2 MAXB - 1 (k, five taps), then the triple loop over f1, f2 < MAXB = 2^LOGB with four cycles per pair")
    G = FsmDrawing(D, "bispectrum_engine.v")
    q1, q2, q3 = 49.0, 29.5, 10.0
    G.node("S_IDLE", X[0], q1, _GW, h, "io", "wait for clr / go")
    G.node("S_CLR", X[0], q2, _GW, h, "control", "acc[cnt] <= 0")
    G.node("S_END", X[0], q3, _GW, h, "io", "busy <= 0, done <= 1")
    G.node("W_A", X[1], q1, _GW, h, "memory", "f_addr =\n(k + 2 - t) mod 512")
    G.node("W_B", X[2], q1, _GW, h, "control", "RAM read latency")
    G.node("W_C", X[3], q1, _GW, h, "compute", "wacc += W[t] F(k - d)")
    G.node("W_F", X[4], q1, _GW, h, "accum", "write Fw(k) =\nrsat(wacc)")
    G.node("R_A", X[1], q2, _GW, h, "memory", "ra_addr = f1")
    G.node("R_B", X[2], q2, _GW, h, "control", "RAM read latency")
    G.node("R_C", X[3], q2, _GW, h, "accum", "latch Fw(f1), f2 <= 0")
    G.node("T_A", X[4], q2, _GW, h, "memory", "ra = f2, rb = f1 + f2")
    G.node("T_B", X[4], q3, _GW, h, "control", "RAM read latency")
    G.node("T_C", X[3], q3, _GW, h, "compute", "P12 = round(F1 F2 /\n2^(WL-1)), latch F3")
    G.node("T_D", X[2], q3, _GW, h, "accum", "acc[f1][f2] +=\nP12 conj(F3)")
    G.down("S_IDLE", "S_CLR", "clr", loff=(1.2, 0), ha="left")
    G.down("S_CLR", "S_END", "all bins", loff=(1.2, 0), ha="left")
    A, B = G.nodes["S_END"], G.nodes["S_IDLE"]
    G.edge("S_END", "S_IDLE", [A.l, (lane, A.l[1]), (lane, B.l[1]), B.l], None, dashed=True)
    D.text(lane + 1.5, (A.l[1] + B.l[1]) / 2, "done pulse", fs=5.2, rot=90, color=INK)
    G.right("S_IDLE", "W_A", "go", loff=(0, 1.1))
    G.right("W_A", "W_B")
    G.right("W_B", "W_C")
    G.right("W_C", "W_F", "t == 4", loff=(0, 1.1))
    G.over("W_C", "W_A", 4.5, "t != 4:  t++", xa=-3.0, xb=4.0, loff=(0, 1.1))
    G.over("W_F", "W_A", 9.0, "k != 2 MAXB - 1:  k++", xa=0.0, xb=-4.0, loff=(0, 1.1))
    A, B = G.nodes["W_F"], G.nodes["R_A"]
    yg = q1 - 4.0
    G.edge("W_F", "R_A", [A.b, (A.b[0], yg), (B.t[0], yg), (B.t[0], B.y + B.h)], "k == 2 MAXB - 1", lpos=0.62, loff=(0, 1.1))
    G.right("R_A", "R_B")
    G.right("R_B", "R_C")
    G.right("R_C", "T_A")
    G.down("T_A", "T_B")
    G.left("T_B", "T_C")
    G.left("T_C", "T_D")
    A, B = G.nodes["T_D"], G.nodes["S_END"]
    G.edge("T_D", "S_END", [A.l, B.r], "f1, f2 = MAXB - 1", loff=(0, 1.1))
    A, B = G.nodes["T_D"], G.nodes["R_A"]
    yk = q3 + h + 3.8
    G.edge("T_D", "R_A", [(A.t[0], A.y + A.h), (A.t[0], yk), (B.t[0], yk), (B.t[0], B.y)], "f2 == MAXB - 1, f1 != MAXB - 1:  f1++", lpos=0.55,
           loff=(0, 1.1), dashed=True)
    A, B = G.nodes["T_D"], G.nodes["T_A"]
    ylow, xr = q3 - 3.6, X[4] + _GW + 2.6
    G.edge("T_D", "T_A", [A.b, (A.b[0], ylow), (xr, ylow), (xr, B.y + B.h / 2), (B.x + B.w, B.y + B.h / 2)], "f2 != MAXB - 1:  f2++", lpos=0.38,
           loff=(0, 1.3), dashed=True)
    G.check()
    maxb = 32
    nb = 2 * maxb
    cyc = 16 * nb + 3 * maxb + 4 * maxb * maxb
    D.text(6, 4.8, f"MAXB = 32: window {nb} x (5 x 3 + 1) = {16 * nb}, f1 loop {maxb} x 3 = {3 * maxb}, pair loop {maxb * maxb} x 4 = {4 * maxb * maxb}"
           f"  ->  {cyc} cycles in the loops;  the testbench measures 5219 per frame in the engine phase (go and done handshakes included)", fs=5.4, ha="left")
    D.save(FIGS / "fig17_fsm_engines")


# ---------------------------------------------------------------------------------------------------------------
# fig18: pipeline timing
# ---------------------------------------------------------------------------------------------------------------
def gantt_bar(D, x, y, w, h, kind, rows, fs=4.6, hatch=False, lw=0.8):
    """One operation bar of a timing chart; ``rows`` are the text lines inside."""
    edge, band, _ = _colors(kind)
    D.ax.add_patch(Rectangle((x, y), w, h, fc="#ffffff" if hatch else band, ec=edge, lw=lw, hatch="////" if hatch else None,
                             zorder=3))
    n = len(rows)
    step = fs * 0.184
    for i, t in enumerate(rows):
        yy = y + h / 2 + ((n - 1) / 2 - i) * step
        D.ax.text(x + w / 2, yy, t, fontsize=fs, color=INK, ha="center", va="center", zorder=4,
                  bbox=dict(fc="#ffffff", ec="none", pad=0.6) if hatch else None)


def gantt_axis(D, x0, y, n, cw, first=0, fs=5.0, ytop=None, ybot=None):
    for c in range(n):
        D.text(x0 + (c + 0.5) * cw, y, str(first + c), fs=fs, color=INK)
    if ytop is not None:
        for c in range(n + 1):
            D.ax.plot([x0 + c * cw] * 2, [ybot, ytop], color=GRID, lw=0.45, zorder=1)


def reg_line(D, x, ytop, ybot, label, dy=1.8, fs=4.8):
    D.ax.plot([x, x], [ybot, ytop], color=INK2, lw=0.6, ls=(0, (3, 2)), zorder=2.5)
    if label:
        D.text(x, ytop + dy, label, fs=fs, color=INK)


def fig18_pipeline_timing():
    w = reg_widths(RTL_DIR / "cru_rotate.v", {"WL": 16})
    stage_regs = [("R1", "x, a, b, inv, quad, dir", ["x1r", "x1i", "a1", "b1", "inv1", "quad1", "ccw1"]),
                  ("R2", "y, inv, quad", ["y2r", "y2i", "inv2", "quad2"]),
                  ("R3", "z, quad", ["z3r", "z3i", "quad3"]),
                  ("R4", "yr, yi (output)", ["y4r", "y4i"])]
    assert len(stage_regs) == LAT_CRU
    LC, LB = LAT_CRU, LAT_BFLY
    D = Diagram(150, 189, "Pipeline timing: rotator, butterfly and write-back loop",
                subtitle="Cycle-level view read from the register stages of the Verilog; the slack of the write-back is computed from the issue schedule.")

    # ---- (a) cru_rotate ---------------------------------------------------------------------------------------------
    _panel_card(D, 2, 128, 146, 52, "(a)", f"cru_rotate: {LC} register stages, one rotation per cycle",
                "x and q are presented in cycle t, yr and yi are valid in cycle t + 4; the combinational octant reduction and ROM lookup belong to the first stage")
    x0, cw, n = 40.0, 9.0, 8
    stages = [("rom", ["octant,", "ROM"]), ("compute", ["rotate", "(a, b)"]), ("compute", ["gain x inv", "round, >> 14"]),
              ("compute", ["quadrant,", "saturate"])]
    ytop_a, ybot_a = 166.4, 140.6
    gantt_axis(D, x0, 167.6, n, cw, ytop=ytop_a, ybot=ybot_a)
    D.text(x0 - 2.0, 167.6, "cycle", fs=5.0, ha="right", color=INK2)
    for i in range(4):
        yb = 159.4 - 6.2 * i
        D.text(x0 - 2.0, yb + 2.6, "input n" + ("" if i == 0 else f" + {i}"), fs=5.4, ha="right", color=INK)
        for st, (kind, rows) in enumerate(stages):
            gantt_bar(D, x0 + (i + st) * cw + 0.3, yb, cw - 0.6, 5.2, kind, rows)
        xf = x0 + (i + LC) * cw + 0.3
        D.ax.add_patch(FancyBboxPatch((xf, yb), cw - 0.6, 5.2, boxstyle="round,pad=0,rounding_size=0.8", fc="none", ec=INK2,
                                      lw=0.7, ls=(0, (3, 2)), zorder=3))
        D.text(xf + (cw - 0.6) / 2, yb + 2.6, "y valid", fs=4.6, color=INK)
    for k, (name, _, _) in enumerate(stage_regs, start=1):
        reg_line(D, x0 + k * cw, ytop_a, ybot_a, name)
    yb = 136.6
    D.ax.plot([x0, x0 + LC * cw], [yb, yb], color=INK, lw=0.9)
    for xx in (x0, x0 + LC * cw):
        D.ax.plot([xx, xx], [yb - 0.8, yb + 0.8], color=INK, lw=0.9)
    D.text(x0 + LC * cw / 2, yb - 2.1, f"LAT_CRU = {LC}", fs=5.8, color=INK, bold=True)
    D.text(x0, 131.0, "inputs n .. n + 3 overlap: one result per cycle from cycle 4 on; the quadrant step swaps and negates, it needs no multiplier", fs=5.2, ha="left")
    # register table
    D.text(116.0, 166.4, "register stage", fs=5.2, ha="left", color=INK, bold=True)
    D.text(145.0, 166.4, "bits", fs=5.2, ha="right", color=INK, bold=True)
    for k, (name, what, regs) in enumerate(stage_regs):
        yy = 163.4 - 3.0 * k
        D.text(116.0, yy, f"{name}  {what}", fs=4.9, ha="left", color=INK2)
        D.text(145.0, yy, str(sum(w[r] for r in regs)), fs=5.2, ha="right", color=INK)
    D.text(116.0, 150.6, "WL = 16; widths read from", fs=4.7, ha="left")
    D.text(116.0, 148.4, "the reg declarations of cru_rotate.v", fs=4.7, ha="left")

    # ---- (b) r8_butterfly ----------------------------------------------------------------------------------------------
    _panel_card(D, 2, 62, 146, 64, "(b)", f"r8_butterfly: LAT_BFLY = {LB} cycles, one butterfly per cycle",
                f"stage 1, a bank of rotators, one register stage for the two remaining add/sub stages, a second bank of rotators (1 + {LC} + 1 + {LC} = {LB})")
    x0, cw, n = 52.0, 7.8, 12
    ytop_b, ybot_b = 113.0, 70.6
    gantt_axis(D, x0, 113.8, n, cw, ytop=ytop_b + 0.4, ybot=ybot_b)
    D.text(x0 - 2.0, 113.8, "cycle", fs=5.0, ha="right", color=INK2)
    lanes = ["stage 1: distance-4 half add/sub", "W8^i rotators, 3 x cru_rotate", f"delay line of s and d0, {LC} stages",
             "stages 2 + 3: distance-2 and -1", "twiddle rotators, 7 x cru_rotate", f"delay line of y0, {LC} stages",
             f"delay line of tq, {LC + 2} stages"]
    ys = [107.8 - 6.2 * j for j in range(len(lanes))]
    for lab, yb in zip(lanes, ys):
        D.text(5.0, yb + 2.4, lab, fs=5.3, ha="left", color=INK)
    hb = 4.8
    gantt_bar(D, x0 + 0.3, ys[0], cw - 0.6, hb, "compute", ["hadd,", "hsub"])
    for bank, c0 in ((1, 1), (4, 6)):
        for st, (kind, _) in enumerate(stages):
            gantt_bar(D, x0 + (c0 + st) * cw + 0.3, ys[bank], cw - 0.6, hb, kind, [("ROM", "rot", "gain", "quad")[st]])
    gantt_bar(D, x0 + 1 * cw + 0.3, ys[2], LC * cw - 0.6, hb, "control", [f"{LC} cycles"], hatch=True)
    gantt_bar(D, x0 + 5 * cw + 0.3, ys[3], cw - 0.6, hb, "compute", ["2 + 1", "stages"])
    gantt_bar(D, x0 + 6 * cw + 0.3, ys[5], LC * cw - 0.6, hb, "control", [f"{LC} cycles"], hatch=True)
    gantt_bar(D, x0 + 0 * cw + 0.3, ys[6], (LC + 2) * cw - 0.6, hb, "control", [f"{LC + 2} cycles, tq arrives with y3"], hatch=True)
    for k, lab in ((1, "s1, d1"), (6, "y3"), (LB, "y valid")):
        reg_line(D, x0 + k * cw, ytop_b, ybot_b - 0.4, lab, dy=2.6)
    D.ax.plot([x0 + LB * cw] * 2, [ybot_b - 0.4, ytop_b], color=INK, lw=1.0, zorder=2.6)
    D.text(x0, 66.0, f"y[0] leaves through the delay line, y[1] .. y[7] through the twiddle rotators with q = tq[k]; the words of the next "
           f"butterfly follow one cycle behind", fs=5.2, ha="left")

    # ---- (c) write-back hazard ----------------------------------------------------------------------------------------
    stall = stall_cycles(LAT)
    _panel_card(D, 2, 2, 146, 58, "(c)", f"fft512_mem: why the issue controller stalls {stall} cycles after pass 0",
                f"a butterfly issued in slot t reads the banks, is written in slot t + {LB + 1} and its rows can be read again in slot t + {LAT} = t + LAT")
    s0, ns, cw, x0 = 56, 21, 3.5, 22.0
    gantt_axis(D, x0, 49.4, ns, cw, first=s0, fs=4.0, ytop=47.2, ybot=20.6)
    D.text(x0 - 1.0, 49.4, "slot", fs=5.0, ha="right", color=INK2)
    ya, yb_, yc, yd, hh = 42.0, 36.0, 30.0, 21.4, 4.4
    for lab, yy in (("pass 0", ya), ("stall", yb_), ("pass 1", yc), ("pass 0, b = 56", yd)):
        D.text(x0 - 1.0, yy + hh / 2, lab, fs=5.2, ha="right", color=INK)
    for b in range(56, 64):
        gantt_bar(D, x0 + (b - s0) * cw + 0.2, ya, cw - 0.4, hh, "memory", [str(b)], fs=4.0)
    gantt_bar(D, x0 + (64 - s0) * cw + 0.2, yb_, stall * cw - 0.4, hh, "control", [f"stall {stall}"], hatch=True, fs=4.6)
    for b in range(0, 9):
        gantt_bar(D, x0 + (64 + stall + b - s0) * cw + 0.2, yc, cw - 0.4, hh, "memory", [str(b)], fs=4.0)
    t = 56 - s0
    gantt_bar(D, x0 + t * cw + 0.2, yd, cw - 0.4, hh, "memory", ["rd"], fs=3.8)
    gantt_bar(D, x0 + (t + 1) * cw + 0.2, yd, LB * cw - 0.4, hh, "compute", [f"r8_butterfly, {LB} cycles"], fs=4.8)
    gantt_bar(D, x0 + (t + 1 + LB) * cw + 0.2, yd, cw - 0.4, hh, "accum", ["wr"], fs=3.8)
    xr = x0 + (t + LAT) * cw
    D.arrow([(xr + 0.2, yd + hh), (xr + 0.2, yc)], None, ls="--")
    D.arrow([(x0 + t * cw + cw / 2, ya), (x0 + t * cw + cw / 2, yd + hh)], None, ls=":")
    D.text(xr + 2.0, yd - 2.2, f"readable from slot 56 + LAT = {56 + LAT}", fs=5.0, ha="left", color=INK)
    D.text(x0, 14.4, f"Pass-1 butterfly (hi, lo) needs the words written by pass-0 butterflies lo, lo + 8, ..., lo + 56. For butterflies 0 .. 7", fs=5.0, ha="left")
    D.text(x0, 12.0, f"the last producer is butterfly 56 + lo: issue slot {64 + stall} + lo is exactly the first slot in which the word can be read.", fs=5.0, ha="left")
    D.text(x0, 9.6, f"Without the stall the first butterfly of pass 1 would read its last row {stall} cycles before it is written; "
           f"no further stall is needed.", fs=5.0, ha="left")
    slack = write_back_slack()
    ax = D.inset(110.0, 9.0, 34.0, 32.0, title="slack of the write-back (cycles)", fs=4.6)
    b = range(64)
    ax.step(b, slack[1], where="mid", color=PAL["blue"], lw=0.9, label="pass 1")
    ax.step(b, slack[2], where="mid", color=PAL["aqua"], lw=0.9, label="pass 2")
    ax.set_xlim(0, 63)
    ax.set_ylim(-4, 62)
    ax.set_xticks([0, 16, 32, 48, 63])
    ax.set_yticks([0, 20, 40, 60])
    ax.set_xlabel("butterfly b", fontsize=4.8, labelpad=1.0)
    ax.axhline(0, color=PAL["red"], lw=0.6, ls="--")
    ax.text(62, 2.2, "hazard below", ha="right", va="bottom", fontsize=4.2, color=PAL["red"])
    ax.legend(fontsize=4.4, frameon=False, loc="center left", borderpad=0.2, handlelength=1.2)
    D.save(FIGS / "fig18_pipeline_timing")


FIGURES = [fig15_hierarchy, fig16_top_fsms, fig17_engine_fsms, fig18_pipeline_timing]


def main():
    for f in FIGURES:
        f()
        print("wrote", f.__name__)


if __name__ == "__main__":
    main()
