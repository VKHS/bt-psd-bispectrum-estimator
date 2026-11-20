"""Figure sheets 9-11: conflict-free address map, FFT schedule, fractional-delay kernel gallery.

Everything is computed from the Python model (btpsd.schedule, btpsd.fdfilter); values reported in the paper are
drawn as hollow black diamonds or marked "(paper)"; nothing here is measured hardware."""
import csv

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch, Rectangle

from _common import DATA, FIGS
from _figdata import read_csv
from _sheet import PAPER, card, cmap_blue, labels, legend
from btpsd.fdfilter import PAPER_FD6, fd_frequency_response, fd_group_delay, fd_kernel
from btpsd.schedule import bank_of, check_conflict_free, pass_addresses, row_of, schedule_cycles
from btpsd.viz import INK, INK2, MUTED, PAL, SERIES, Diagram, tint

BLUE, ORANGE, AQUA, VIOLET, YELLOW = PAL["blue"], PAL["orange"], PAL["aqua"], PAL["violet"], PAL["yellow"]
DASH = (0, (3.0, 1.6))
DIGIT = ("a2", "a1", "a0")
FIXED = ("a1, a0", "a2, a0", "a2, a1")
TWIDDLE = ("W512^(b k) with b = 8 hi + lo", "W512^(8 lo k)", "none (no twiddle)")


def fig09_address_map():
    assert check_conflict_free()
    W, H = 130, 137
    D = Diagram(W, H, "Conflict-free address map of the in-place radix-8 FFT",
                subtitle="a = 64 a2 + 8 a1 + a0,   bank = (a2 + a1 + a0) mod 8,   row = 8 a2 + a1.   "
                         "The 8 words of a butterfly differ in one digit, so they sit in 8 different banks.")
    b_ex = 19
    cm = cmap_blue()
    for p in range(3):
        adr = pass_addresses(p)
        slot = np.zeros((8, 64), dtype=int)
        who = np.zeros((8, 64), dtype=int)
        for b in range(64):
            for m in range(8):
                a = int(adr[b, m])
                slot[int(bank_of(a)), int(row_of(a))] = m
                who[int(bank_of(a)), int(row_of(a))] = b
        hh = (25.0, 25.0, 30.6)[p]
        y0 = (102.5, 76.1, 44.1)[p]
        ax = card(D, 2, y0, 126, hh, f"({'abc'[p]})", f"Pass {p}: a butterfly combines the 8 words that differ in {DIGIT[p]}",
                  f"slot m = {DIGIT[p]};  fixed digits ({FIXED[p]}) = (hi, lo) with b = 8 hi + lo;  twiddle after the butterfly: {TWIDDLE[p]}",
                  pad=(9.5, 7.2 if p == 2 else 1.8, 3.5, 8.4))
        ax.grid(False)
        ax.imshow(slot, cmap=cm, vmin=-1.5, vmax=7.5, aspect="auto", interpolation="nearest", extent=(-0.5, 63.5, 7.5, -0.5))
        for bk, rw in zip(*np.where(who == b_ex)):
            ax.add_patch(Rectangle((rw - 0.5, bk - 0.5), 1, 1, fc=ORANGE, ec="#ffffff", lw=0.5, zorder=3))
        ax.set_yticks(range(8))
        ax.set_yticklabels([f"M{i + 1}" for i in range(8)], fontsize=4.8)
        ax.set_xticks([0, 8, 16, 24, 32, 40, 48, 56, 63])
        if p < 2:
            ax.set_xticklabels([])
        else:
            labels(ax, "row in the bank (0 .. 63)")
        for sp in ("left", "bottom"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(length=0)
        D.text(124.5, y0 + hh - 3.4, f"orange: the 8 words of butterfly b = {b_ex}   |   shade: slot m (light 0 ... dark 7)", fs=5.4,
               ha="right", color=INK2)
    # ---- (d) rotation table ---------------------------------------------------------------------------------------------
    ax = card(D, 2, 2, 38, 40, "(d)", "Rotation s(b)", "slot m sits in bank (m + s) mod 8, in every pass")
    ax.grid(False)
    tab = (np.arange(8)[:, None] + np.arange(8)[None, :]) % 8
    ax.imshow(tab, cmap=cm, vmin=-2.5, vmax=7.5, aspect="auto", interpolation="nearest")
    for hi in range(8):
        for lo in range(8):
            ax.text(lo, hi, str(tab[hi, lo]), ha="center", va="center", fontsize=5.2, color="#ffffff" if tab[hi, lo] >= 5 else INK)
    ax.set_xticks(range(8))
    ax.set_yticks(range(8))
    ax.tick_params(length=0)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_visible(False)
    labels(ax, "lo = b & 7", "hi = b >> 3")
    # ---- (e) bank accessed by each slot -------------------------------------------------------------------------------------
    strip = np.full((26, 64), np.nan)
    for p in range(3):
        strip[9 * p: 9 * p + 8, :] = bank_of(pass_addresses(p)).T
    ax = card(D, 42, 2, 58, 40, "(e)", "Bank read by slot m of butterfly b", "rows: slot m = 0 .. 7;  colour: bank", pad=(8.8, 11.6, 3.0, 9.4))
    ax.grid(False)
    ax.imshow(np.ma.masked_invalid(strip), cmap=ListedColormap(SERIES[:8]), vmin=-0.5, vmax=7.5, aspect="auto",
              interpolation="nearest", extent=(-0.5, 63.5, 25.5, -0.5))
    ax.set_yticks([3.5, 12.5, 21.5])
    ax.set_yticklabels(["pass 0", "pass 1", "pass 2"], fontsize=5.0)
    ax.set_xticks([0, 8, 16, 24, 32, 40, 48, 56, 63])
    ax.tick_params(length=0)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_visible(False)
    labels(ax, "butterfly b (natural order)")
    ax.legend(handles=[Patch(fc=SERIES[i], label=f"M{i + 1}") for i in range(8)], ncol=8, fontsize=4.8, frameon=False,
              loc="upper center", bbox_to_anchor=(0.5, -0.26), handlelength=0.9, handleheight=0.9, columnspacing=0.7, handletextpad=0.35)
    # ---- (f) checks ----------------------------------------------------------------------------------------------------------------
    card(D, 102, 2, 26, 40, "(f)", "Checks", "btpsd.schedule", plot=False)
    lines = [("check_conflict_free() = True", True), ("3 passes x 64 butterflies", False), ("x 8 words each", False), ("", False),
             ("8 distinct banks per", False), ("butterfly; every word is", False), ("used once per pass", False), ("", False),
             ("in place: the 8 results", False), ("return to the 8 addresses", False), ("just read", False), ("", False),
             ("after pass 2, address", False), ("(a2, a1, a0) holds", False), ("X[a2 + 8 a1 + 64 a0]", True)]
    yy = 31.4
    for ln, bold in lines:
        if ln:
            D.text(103.8, yy, ln, fs=5.5, ha="left", color=INK if bold else INK2, bold=bold)
            yy -= 1.9
        else:
            yy -= 0.9
    D.save(FIGS / "fig09_bank_address_map")


def _gantt(ax, lat, ytop, hatch_col):
    sc = schedule_cycles(lat)
    st = sc["stalls_per_pass"]
    t = 0
    for p in range(3):
        if st[p]:
            ax.barh(ytop - p, st[p], left=t, height=0.7, fc="#ffffff", ec=hatch_col, hatch="////", lw=0.5, zorder=2)
            t += st[p]
        ax.barh(ytop - p, 64, left=t, height=0.7, color=AQUA, zorder=2)
        ax.text(t + 32, ytop - p, f"{t}..{t + 63}", ha="center", va="center", fontsize=4.9, color="#ffffff", zorder=3)
        t += 64
    ax.barh(ytop - 2, lat - 1, left=t, height=0.7, color=tint(BLUE, 0.55), zorder=2)
    return sc


def fig10_schedule():
    W, H = 130, 98
    D = Diagram(W, H, "Schedule of one 512-point FFT: 3 passes of 64 butterflies",
                subtitle="One butterfly issues per cycle and its 8 results are written back L cycles later; "
                         "a pass-1 butterfly cannot issue before its 8 input words are written.")
    ax = card(D, 2, 52, 126, 36.5, "(a)", "Issue timeline for two pipeline latencies L", None, pad=(10.0, 7.0, 3.0, 6.2))
    for lat, ytop in ((6, 8.0), (14, 3.5)):
        sc = _gantt(ax, lat, ytop, ORANGE)
        nst = sum(sc["stalls_per_pass"])
        ax.text(1.5, ytop + 0.78, f"L = {lat}:  {sc['total_cycles']} cycles" + (f",  {nst} stall cycles" if nst else ",  no stall"),
                ha="left", va="center", fontsize=5.8, color=INK, fontweight="bold")
    ax.set_yticks([8, 7, 6, 3.5, 2.5, 1.5])
    ax.set_yticklabels(["pass 0", "pass 1", "pass 2"] * 2, fontsize=5.0)
    ax.set_ylim(0.6, 9.0)
    ax.set_xlim(0, 262)
    ax.set_xticks(range(0, 221, 32))
    labels(ax, "clock cycle")
    ax.grid(axis="y", visible=False)
    ax.legend(handles=[Patch(fc=AQUA, label="issue (one butterfly per cycle)"), Patch(fc="#ffffff", ec=ORANGE, hatch="////", label="stall (data not yet written)"),
                       Patch(fc=tint(BLUE, 0.55), label="write-back drain of the last butterfly")], fontsize=5.2, frameon=False,
              loc="upper right", bbox_to_anchor=(1.0, 1.02), handlelength=1.3, handleheight=0.9, borderpad=0.1, labelspacing=0.25)
    # ---- (b) dependence matrix ------------------------------------------------------------------------------------------------------------
    dep = np.zeros((64, 64))
    a0 = pass_addresses(0)
    a1 = pass_addresses(1)
    producer = {int(a): b for b in range(64) for a in a0[b]}
    for bp in range(64):
        for a in a1[bp]:
            dep[bp, producer[int(a)]] = 1
    ax = card(D, 2, 2, 50, 47, "(b)", "Pass 1 reads all of pass 0", "dot: pass-0 butterfly b writes a word read by pass-1 butterfly b'",
              pad=(8.8, 7.4, 3.0, 12.0))
    D.text(3.8, 2 + 47 - 9.7, "orange: the 8 producers of b' = 0; the last (b = 56) writes at cycle 56 + L", fs=5.2, ha="left", color=INK)
    ax.grid(False)
    ax.imshow(dep, cmap=ListedColormap(["#ffffff", BLUE]), aspect="auto", interpolation="nearest", extent=(-0.5, 63.5, 63.5, -0.5))
    for b in np.where(dep[0] > 0)[0]:
        ax.add_patch(Rectangle((b - 0.5, -0.5), 1, 1, fc=ORANGE, ec=ORANGE, lw=0.8, zorder=3))
    ax.set_xticks([0, 16, 32, 48, 63])
    ax.set_yticks([0, 16, 32, 48, 63])
    labels(ax, "pass-0 butterfly b (producer)", "pass-1 butterfly b' (consumer)")
    # ---- (c) cycle count versus latency -------------------------------------------------------------------------------------------------------
    rows = read_csv("fft512_cycles_vs_latency.csv")
    lat = [int(r["pipeline_latency_cycles"]) for r in rows]
    tot = [int(r["total_cycles"]) for r in rows]
    iss = [int(r["issue_cycles"]) for r in rows]
    ax = card(D, 54, 2, 74, 47, "(c)", "Cycles per chunk versus pipeline latency", "model of btpsd.schedule.schedule_cycles")
    ax.plot(lat, tot, "-o", color=BLUE, lw=1.0, ms=3.0, label="total (first issue to last write)")
    ax.plot(lat, iss, "-s", color=AQUA, lw=0.9, ms=2.6, label="issue slots incl. stalls")
    ax.axvline(8, color=MUTED, lw=0.7, ls=DASH)
    ax.text(8.2, 217.5, "L = 8", fontsize=5.2, color=INK2, ha="left")
    ax.text(0.03, 0.60, "no stall for L <= 8:  total = 192 + L - 1\nL > 8:  stall = L - 8,  total = 192 + (L - 8) + L - 1",
            transform=ax.transAxes, fontsize=5.2, color=INK, va="center")
    ax.set_xlim(3.4, 16.6)
    ax.set_ylim(188, 220)
    ax.set_xticks(range(4, 17, 2))
    labels(ax, "pipeline latency L (cycles)", "cycles")
    legend(ax, "lower right")
    D.save(FIGS / "fig10_fft_schedule")


def fig11_fd_kernel_gallery():
    W, H = 130, 92
    D = Diagram(W, H, "Half-sample fractional-delay kernels with 3 to 8 taps",
                subtitle="Model of this repository (sinc design; 6 taps = Eq. 6 constants).  Even lengths are symmetric, hence linear phase; "
                         "counterpart of paper Figs. 6-11.")
    with open(DATA / "paper_table3_fd_taps_snr.csv") as f:
        pr = {int(r["fd_taps"]): float(r["snr_db"]) for r in csv.DictReader(f)}
    rows = read_csv("fd_merge_snr_vs_taps.csv")
    mod = {int(r["taps"]): float(r["snr_complex_spectrum_db"]) for r in rows if r["signal"] == "eeg_like" and r["kernel"] == "sinc"}
    x0, step, w = 2.0, 21.2, 20.4
    for c, L in enumerate(range(3, 9)):
        h = PAPER_FD6 if L == 6 else fd_kernel(L, "sinc")[0]
        col = BLUE if L % 2 == 0 else ORANGE
        x = x0 + c * step
        card(D, x, 2, w, 77, f"{L} taps", None, None, plot=False, fs=8.0)
        D.text(x + 1.8, 72.2, f"cf. paper Fig. {L + 3}", fs=5.4, ha="left", color=INK2)
        D.text(x + 1.8, 69.4, f"paper {pr[L]:.1f} dB (Table 3)", fs=5.2, ha="left", color=INK, italic=True)
        D.text(x + 1.8, 66.8, f"model {mod[L]:.1f} dB, odd bins", fs=5.2, ha="left", color=INK2)
        # impulse response
        ax = D.inset(x + 6.2, 47.5, w - 8.0, 14.2, title="impulse response", fs=5.0)
        ax.axhline(0, color=MUTED, lw=0.5)
        ax.vlines(np.arange(L), 0, h, color=col, lw=1.0)
        ax.plot(np.arange(L), h, "o", ms=2.4, color=col)
        ax.set_ylim(-0.4, 0.75)
        ax.set_xlim(-0.7, 7.7)
        ax.set_xticks(range(0, 8, 2))
        # group delay
        w_, gd = fd_group_delay(h)
        ax = D.inset(x + 6.2, 26.5, w - 8.0, 14.2, title="group delay (samples)", fs=5.0)
        ax.plot(w_ / np.pi, gd, color=col, lw=1.0)
        ax.set_ylim(0, 5)
        ax.set_xlim(0, 1)
        ax.set_xticks([0, 0.5, 1])
        ax.set_xticklabels(["0", ".5", "1"])
        # magnitude
        w_, Hr = fd_frequency_response(h)
        ax = D.inset(x + 6.2, 6.5, w - 8.0, 14.2, title="|H| (dB)", fs=5.0)
        ax.axhline(0, color=MUTED, lw=0.5, ls=DASH)
        ax.plot(w_ / np.pi, 20 * np.log10(np.abs(Hr) + 1e-9), color=col, lw=1.0)
        ax.set_ylim(-30, 5)
        ax.set_xlim(0, 1)
        ax.set_xticks([0, 0.5, 1])
        ax.set_xticklabels(["0", ".5", "1"])
        if c == 0:
            D.text(x + w / 2, 3.0, "frequency (x pi rad / sample)", fs=5.2, color=INK2)
    D.save(FIGS / "fig11_fd_kernel_gallery")


if __name__ == "__main__":
    fig09_address_map()
    fig10_schedule()
    fig11_fd_kernel_gallery()
    print("figures 9-11 written")
