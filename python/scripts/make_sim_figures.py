"""Figure sheets drawn from the RTL-simulation results in ``results/rtl_*.csv`` (19: interface timing, 20: cycle budget).

Every number on these sheets was printed by the testbenches during the Icarus Verilog run of the CI workflow
(``scripts/ingest_ci_traces.py`` stores it, ``results/rtl_sim_provenance.json`` says which run); the Python models enter
only as hollow markers that are compared with the simulation, never as data.
"""
from __future__ import annotations

import csv
import json
import math

import numpy as np
from matplotlib.patches import Polygon, Rectangle

from _common import FIGS, RESULTS
from btpsd import rtlmodel
from btpsd import vectors as V
from btpsd.rtlparse import RTL_DIR, define_values
from btpsd.viz import GRID, INK, INK2, MUTED, PAL, Diagram, tint
from make_rtl_figures import _colors, _panel_card

ROOT = RESULTS.parent
STATE_KIND = {"S_IDLE": "io", "S_END": "io", "S_CLR": "control", "M_A": "memory", "M_B": "control", "M_C": "compute",
              "M_F": "accum", "W_A": "memory", "W_B": "control", "W_C": "compute", "W_F": "accum", "W_G": "rom",
              "W_S": "rom", "W_ACC": "accum", "R_A": "memory", "R_B": "control", "R_C": "accum", "T_A": "memory",
              "T_B": "control", "T_C": "compute", "T_D": "accum"}


def rows(name):
    with (RESULTS / name).open() as f:
        return list(csv.DictReader(f))


def provenance():
    return json.loads((RESULTS / "rtl_sim_provenance.json").read_text())


# ---------------------------------------------------------------------------------------------------------------
# waveform helpers
# ---------------------------------------------------------------------------------------------------------------
def wave_bit(D, x0, y, h, cw, vals, color=INK):
    xs, ys = [], []
    for i, v in enumerate(vals):
        xs += [x0 + i * cw, x0 + (i + 1) * cw]
        ys += [y + (h if v else 0.0)] * 2
    D.ax.plot([x0, x0 + len(vals) * cw], [y, y], color=GRID, lw=0.5, zorder=2)
    D.ax.plot(xs, ys, color=color, lw=0.95, zorder=3, solid_joinstyle="miter")


def wave_bus(D, x0, y, h, cw, vals, fmt=str, fs=4.0, color=INK2, fill=None):
    """Bus waveform: one hexagon per run of equal values, the value in the middle."""
    runs, i = [], 0
    while i < len(vals):
        j = i
        while j + 1 < len(vals) and vals[j + 1] == vals[i]:
            j += 1
        runs.append((i, j, vals[i]))
        i = j + 1
    e = min(0.45, cw * 0.18)
    for a, b, v in runs:
        xa, xb = x0 + a * cw, x0 + (b + 1) * cw
        pts = [(xa + e, y), (xb - e, y), (xb, y + h / 2), (xb - e, y + h), (xa + e, y + h), (xa, y + h / 2)]
        D.ax.add_patch(Polygon(pts, closed=True, fc=fill or tint(PAL["blue"], 0.9), ec=color, lw=0.7, zorder=3))
        D.ax.text((xa + xb) / 2, y + h / 2 - 0.05, fmt(v), fontsize=fs, color=INK, ha="center", va="center", zorder=4)


def cycle_axis(D, x0, y, first, n, cw, ytop, ybot, fs=3.8, every=1):
    for c in range(n):
        xx = x0 + (c + 0.5) * cw
        if c % every == 0:
            D.text(xx, y, str(first + c), fs=fs, color=INK2)
    for c in range(n + 1):
        D.ax.plot([x0 + c * cw] * 2, [ybot, ytop], color=GRID, lw=0.4, zorder=1)


def window_break(D, x, y0, y1):
    """Axis break between two windows: a dotted guide and a short double slash at mid height."""
    D.ax.plot([x, x], [y0, y1], color=GRID, lw=0.7, ls=(0, (1, 2)), zorder=1)
    ym = (y0 + y1) / 2
    for dx in (-0.45, 0.45):
        D.ax.plot([x + dx - 0.45, x + dx + 0.45], [ym - 1.6, ym + 1.6], color=MUTED, lw=0.9, zorder=2)


# ---------------------------------------------------------------------------------------------------------------
# fig19: interface timing
# ---------------------------------------------------------------------------------------------------------------
def model_wait_histogram():
    """Histogram of the W_S residence (n + 3 cycles per bin, binned by floor(log2)) predicted by the model for the testbench data."""
    from btpsd.estimator import psd_fixed
    p = V.psd_vectors(12)
    res = psd_fixed(p["x"], p["cfg"], quantize=False)
    return rtlmodel.squarer_wait_histogram(res.X_re, res.X_im, p["cfg"].squarer_units)


def fig19_interface_timing():
    prov = provenance()
    fft = rows("rtl_trace_fft512.csv")
    sq = rows("rtl_trace_squarer.csv")
    col = lambda name, a, b: [int(r[name]) for r in fft if a <= int(r["cycle"]) <= b]  # noqa: E731
    D = Diagram(150, 168, "Handshakes of the FFT core and the squarer: waveforms from the RTL simulation",
                subtitle=f"Signals printed by the testbenches (Icarus Verilog, CI run {prov['run_id']}, commit {prov['commit'][:7]}); "
                         f"sampled at the rising edge, cycle 0 = start high. No Python model is involved.")

    # ---- (a) fft512_mem -----------------------------------------------------------------------------------------------
    top = 158.0
    _panel_card(D, 2, 88, 146, 70, "(a)", "fft512_mem, first transform of tb_fft512_mem (WL = 16)",
                "load, start, run, done: start is sampled in cycle 0, the first butterfly issues in cycle 1, done is a one-cycle pulse")
    names = ["ld_en", "start", "busy", "iss_v", "iss_p", "iss_b", "st", "done"]
    bits = {"ld_en", "start", "busy", "iss_v", "done"}
    ytop_rows = 139.5
    pitch, hh = 6.0, 3.6
    ys = [ytop_rows - pitch * i for i in range(len(names))]
    for nme, y in zip(names, ys):
        D.text(26.0, y + hh / 2, nme, fs=5.2, ha="right", color=INK, mono=True)
    cw = 2.9
    wins = [(-3, 4, 28.0, "start of the run"), (62, 71, 57.2, "pass 0 to pass 1"), (194, 211, 92.2, "last issue to done")]
    last = max(int(r["cycle"]) for r in fft if int(r["iss_v"]))
    for k, (a, b, x0, title) in enumerate(wins):
        n = b - a + 1
        cycle_axis(D, x0, 145.2, a, n, cw, ys[-1] - 1.0, ys[0] + hh + 0.4, every=1 if n <= 10 else 2)
        D.text(x0, 148.2, title, fs=5.6, ha="left", color=INK, bold=True)
        for nme, y in zip(names, ys):
            vals = col(nme if nme != "st" else "st", a, b)
            if nme in bits:
                wave_bit(D, x0, y, hh, cw, vals, color=PAL["blue"] if nme in ("start", "done") else INK)
            elif nme == "st":
                labels = {int(r["st"]): r["state"].replace("ST_", "") for r in fft}
                wave_bus(D, x0, y, hh, cw, vals, fmt=lambda v: labels[v], fs=3.6, fill=tint(PAL["violet"], 0.88))
            else:
                wave_bus(D, x0, y, hh, cw, vals, fs=3.6)
        if k:
            window_break(D, x0 - 2.6, ys[-1] - 1.0, ys[0] + hh + 0.4)
    # annotations
    ynote = 85.2
    t1 = [r for r in fft if int(r["busy"])][0]
    D.text(28.0, ynote + 8.0, f"start in cycle 0; busy and the first issue", fs=4.9, ha="left", color=INK)
    D.text(28.0, ynote + 5.8, f"(slot 0) follow in cycle {t1['cycle']}", fs=4.9, ha="left", color=INK)
    stall_cycles = [int(r["cycle"]) for r in fft if 60 <= int(r["cycle"]) <= 75 and not int(r["iss_v"]) and int(r["busy"])]
    D.text(57.2, ynote + 8.0, f"STALL = {len(stall_cycles)}: iss_v is low in cycles {stall_cycles[0]} .. {stall_cycles[-1]},", fs=4.9, ha="left", color=INK)
    D.text(57.2, ynote + 5.8, f"iss_p is already 1; pass 1 issues from cycle {stall_cycles[-1] + 1}", fs=4.9, ha="left", color=INK)
    done_c = [int(r["cycle"]) for r in fft if int(r["done"])][0]
    D.text(92.2, ynote + 8.0, f"last issue (slot {last - 1}) in cycle {last}; ST_DRAIN from cycle {last + 1}; done in cycle {done_c}", fs=4.9, ha="left", color=INK)
    D.text(92.2, ynote + 5.8, f"= 1 + {last - 1} + LAT ({rtlmodel.LAT}); first issue to done: {done_c - 1} cycles", fs=4.9, ha="left", color=INK)

    # ---- (b) mag_sq_odd ------------------------------------------------------------------------------------------------
    _panel_card(D, 2, 44, 146, 42, "(b)", "mag_sq_odd (serial sum-of-odds squarer), directed examples of tb_mag_sq_odd",
                "start with stable re and im; busy for n + 1 cycles, done n + 2 cycles after start with n = max(|re|, |im|); p keeps the previous result until done")
    exs = [(0, 3, -5), (1, 0, 0), (2, -4, 2)]
    names2 = ["start", "busy", "done", "p"]
    ys2 = [68.4 - 5.4 * i for i in range(4)]
    for nme, y in zip(names2, ys2):
        D.text(26.0, y + 1.6, nme, fs=5.2, ha="right", color=INK, mono=True)
    cw2 = 2.6
    for e, re, im in exs:
        x0 = 28.0 + e * 40.0
        data = [r for r in sq if int(r["example"]) == e]
        n = max(abs(re), abs(im))
        cycle_axis(D, x0, 73.0, 0, 14, cw2, ys2[-1] - 0.8, ys2[0] + 3.6, every=2)
        D.text(x0, 75.8, f"(re, im) = ({re}, {im}),  n = {n}", fs=5.2, ha="left", color=INK, bold=True)
        for nme, y in zip(names2, ys2):
            vals = [int(r[nme]) for r in data]
            if nme == "p":
                wave_bus(D, x0, y, 3.2, cw2, vals, fs=3.8)
            else:
                wave_bit(D, x0, y, 3.2, cw2, vals, color=PAL["blue"] if nme != "busy" else INK)
        dn = [int(r["cycle"]) for r in data if int(r["done"])][0]
        D.text(x0, 48.8 - 0.0, f"done in cycle {dn} = n + 2;  p = {int(data[dn]['p'])}", fs=4.9, ha="left", color=INK)

    # ---- (c) W_S residence histogram ---------------------------------------------------------------------------------------
    _panel_card(D, 2, 2, 146, 40, "(c)", "How long psd_engine waits in W_S: one visit per bin, n + 3 cycles each",
                "tb_psd_top, WL = 12, 3 segments x 1 024 bins, serial squarer; hollow markers: the model's prediction from the windowed spectra")
    hist = {int(r["log2_bin"]): int(r["bins"]) for r in rows("rtl_squarer_wait_hist.csv")}
    model = model_wait_histogram()
    ax = D.inset(14.0, 8.0, 62.0, 22.0, fs=4.6)  # the card header ends at y = 35
    ks = sorted(set(hist) | set(model))
    ax.bar(range(len(ks)), [hist.get(k, 0) for k in ks], color=PAL["orange"], width=0.62, label="simulation (RTL)")
    ax.plot(range(len(ks)), [model.get(k, 0) for k in ks], ls="none", marker="D", mfc="none", mec=INK, ms=4.0, mew=0.8, label="model")
    ax.set_xticks(range(len(ks)))
    ax.set_xticklabels([f"{2 ** k}-{2 ** (k + 1) - 1}" for k in ks], fontsize=4.0)
    ax.set_xlabel("cycles in W_S per bin", fontsize=4.8, labelpad=1.0)
    ax.set_ylabel("bins", fontsize=4.8, labelpad=1.0)
    ax.legend(fontsize=4.4, frameon=False, loc="upper right", borderpad=0.2)
    occ = {r["state"]: int(r["cycles"]) for r in rows("rtl_state_occupancy.csv") if r["testbench"] == "psd_engine"}
    d = define_values(ROOT / "sim" / "vectors" / "psd_wl12.vh")
    nbin = 1024 * d["PSD_S"]
    lines = [f"W_S in the whole run: {occ['W_S']:,} cycles over {nbin:,} visits".replace(",", " "),
             f"= {nbin:,} x 3 + {d['PSD_SQUARER_CYCLES']:,}  (the model's sum of n over all bins)".replace(",", " "),
             f"mean residence {occ['W_S'] / nbin:.2f} cycles per bin, longest visits in the {2 ** max(hist)} - {2 ** (max(hist) + 1) - 1} range",
             "the squarer is the only data-dependent part of the engine: with a multiplier W_S would last 2 cycles"]
    for i, t in enumerate(lines):
        D.text(84.0, 28.0 - 3.0 * i, t, fs=5.0, ha="left", color=INK if i < 2 else INK2)
    D.save(FIGS / "fig19_interface_timing")


# ---------------------------------------------------------------------------------------------------------------
# fig20: cycle budget
# ---------------------------------------------------------------------------------------------------------------
def _bars(ax, states, occ, want, per, fs=4.6):
    y = np.arange(len(states))[::-1]
    for yy, st in zip(y, states):
        edge, band, _ = _colors(STATE_KIND[st])
        ax.barh(yy, occ[st], height=0.62, color=band, ec=edge, lw=0.8, zorder=3)
        if st in want:
            ax.plot([want[st]], [yy], marker="D", ms=3.6, mfc="none", mec=INK, mew=0.8, zorder=5)
        ax.text(occ[st] + 0.012 * ax.get_xlim()[1], yy, f"{occ[st]:,}".replace(",", " ") + (f"  ({per[st]})" if st in per else ""), va="center",
                fontsize=fs, color=INK, zorder=6)
    ax.set_yticks(y)
    ax.set_yticklabels(states, fontsize=5.0, family="DejaVu Sans Mono")
    ax.tick_params(axis="y", length=0)


def fig20_cycle_budget():
    prov = provenance()
    occ_rows = rows("rtl_state_occupancy.csv")
    occ = lambda tb: {r["state"]: int(r["cycles"]) for r in occ_rows if r["testbench"] == tb}  # noqa: E731
    psd, pe = occ("psd_top"), occ("psd_engine")
    bsp, be = occ("bsp_top"), occ("bsp_engine")
    dp = define_values(ROOT / "sim" / "vectors" / "psd_wl12.vh")
    db = define_values(ROOT / "sim" / "vectors" / "bsp_wl12.vh")
    S, F = dp["PSD_S"], db["BSP_S"]
    want_p = rtlmodel.psd_engine_state_cycles(S, dp["PSD_SQUARER_CYCLES"])
    want_b = rtlmodel.bispectrum_engine_state_cycles(F, db["BSP_MAX_BIN"].bit_length() - 1)
    D = Diagram(150, 148, "Where the cycles go: state occupancy of the two estimators in the RTL simulation",
                subtitle=f"Cycles per state from reset release to the last result (CI run {prov['run_id']}, commit {prov['commit'][:7]}); "
                         f"hollow diamonds: the per-state formulas of btpsd.rtlmodel.")

    # ---- (a) psd_engine -------------------------------------------------------------------------------------------------
    _panel_card(D, 2, 70, 72, 70, "(a)", "psd_engine", f"tb_psd_top: one clear and {S} segments, WL = 12")
    states = [s for s in rtlmodel.PSD_STATES if s != "S_IDLE"]
    ax = D.inset(16.0, 75.0, 44.0, 53.0, fs=4.6)
    ax.set_xlim(0, 40000)
    per = {s: f"{pe[s] // S:,}".replace(",", " ") for s in states if s not in ("S_CLR", "S_END", "W_S")}
    per["W_S"] = f"{pe['W_S'] / S:,.0f}".replace(",", " ")
    _bars(ax, states, pe, want_p, per)
    ax.set_xticks([0, 10000, 20000, 30000])
    ax.set_xticklabels(["0", "10k", "20k", "30k"], fontsize=4.4)
    ax.set_xlabel("cycles (per segment in brackets)", fontsize=4.8, labelpad=1.0)
    ax.grid(axis="x", color=GRID, lw=0.5)

    # ---- (b) bispectrum_engine ---------------------------------------------------------------------------------------------
    _panel_card(D, 76, 70, 72, 70, "(b)", "bispectrum_engine", f"tb_bispectrum_top: one clear and {F} frames, WL = 12, 32 x 32 bins")
    states = [s for s in rtlmodel.BSP_STATES if s != "S_IDLE"]
    ax = D.inset(90.0, 75.0, 44.0, 53.0, fs=4.6)
    ax.set_xlim(0, 6400)
    per = {s: f"{be[s] // F:,}".replace(",", " ") for s in states if s not in ("S_CLR", "S_END")}
    _bars(ax, states, be, want_b, per)
    ax.set_xticks([0, 2000, 4000, 6000])
    ax.set_xticklabels(["0", "2k", "4k", "6k"], fontsize=4.4)
    ax.set_xlabel("cycles (per frame in brackets)", fontsize=4.8, labelpad=1.0)
    ax.grid(axis="x", color=GRID, lw=0.5)

    # ---- (c) top-level phases ------------------------------------------------------------------------------------------------
    _panel_card(D, 2, 2, 146, 66, "(c)", "Phases of the two top levels",
                "T_LOAD / B_EXT take the samples in, T_FFT / B_FFT run the core, T_SEG / B_ENG the engine; the engine dominates")
    D.legend(32.0, 57.6, [("control", "clear"), ("memory", "load / copy samples"), ("compute", "FFT core"), ("accum", "engine")], dx=28.0)
    ax = D.inset(30.0, 33.0, 108.0, 20.0, fs=4.6)
    ax.set_ylim(-0.75, 1.5)
    phases = [("bt_psd_top", psd, [("T_INIT", "control"), ("T_CLR", "control"), ("T_LOAD", "memory"), ("T_FFT", "compute"), ("T_SEG", "accum")]),
              ("bt_bispectrum_top", bsp, [("B_INIT", "control"), ("B_CLR", "control"), ("B_COPY", "memory"), ("B_EXT", "memory"), ("B_FFT", "compute"), ("B_ENG", "accum")])]
    for yi, (name, o, seq) in enumerate(phases[::-1]):
        left = 0
        for st, kind in seq:
            edge, band, _ = _colors(kind)
            ax.barh(yi, o[st], left=left, height=0.55, color=band, ec=edge, lw=0.8, zorder=3)
            if o[st] > 5000:
                ax.text(left + o[st] / 2, yi, f"{st}  {o[st]:,}".replace(",", " "), ha="center", va="center", fontsize=4.4, color=INK, zorder=4)
            left += o[st]
        ax.text(0, yi - 0.36, "   ".join(f"{st} {o[st]:,}".replace(",", " ") for st, _ in seq) + f"   (sum {left:,})".replace(",", " "),
                ha="left", va="top", fontsize=4.2, color=INK2, zorder=4)
    ax.set_yticks([0, 1])
    ax.set_yticklabels([phases[1][0], phases[0][0]], fontsize=4.8, family="DejaVu Sans Mono")
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, 125000)
    ax.set_xticks([0, 25000, 50000, 75000, 100000, 125000])
    ax.set_xticklabels(["0", "25k", "50k", "75k", "100k", "125k"], fontsize=4.4)
    ax.set_xlabel("cycles per phase, summed over the whole run (stacked by phase, not in time order)", fontsize=4.8, labelpad=1.0)
    ax.grid(axis="x", color=GRID, lw=0.5)
    seg_mean = psd["T_SEG"] / S
    fft_chunk = psd["T_FFT"] // dp["PSD_NCHUNKS"]
    load_chunk = psd["T_LOAD"] // dp["PSD_NCHUNKS"]
    steady = load_chunk + fft_chunk + seg_mean
    busy_b = sum(be[k] for k in rtlmodel.BSP_STATES if k not in ("S_IDLE", "S_CLR", "S_END")) // F
    sp = lambda v: f"{v:,.0f}".replace(",", " ")  # noqa: E731
    notes = [
        f"FFT phase: {fft_chunk} cycles per chunk in both top levels (first issue to done: 207 cycles, plus control).",
        f"PSD segment phase: {sp(seg_mean)} cycles per segment on this data, {sp(pe['W_S'] / S)} of them waiting for the squarer; one chunk costs about "
        f"{sp(steady)} cycles in steady state ({load_chunk} load + {fft_chunk} FFT + one segment).",
        f"Bispectrum engine phase: {sp(busy_b)} state cycles + 3 control cycles = {sp(bsp['B_ENG'] / F)} per frame, independent of the data.",
        "These are the sequential reference engines of this repository (one multiply-accumulate at a time), written to be checked against",
        "the bit-accurate model. Reference only, paper-reported (Table 9, proposed estimator, overlapped architecture): latency 1 252 cycles,",
        "205 Msps at 232 MHz, R = 0.88. No measurement of that design was made here and the two sets of numbers are not comparable.",
    ]
    for i, t in enumerate(notes):
        D.text(8.0, 24.0 - 3.3 * i, t, fs=5.0, ha="left", color=INK if i < 3 else INK2)
    D.save(FIGS / "fig20_cycle_budget")


FIGURES = [fig19_interface_timing, fig20_cycle_budget]


def main():
    for fn in FIGURES:
        fn()
        print("wrote", fn.__name__)


if __name__ == "__main__":
    main()
