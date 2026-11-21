"""Figure sheets 30 and 31: the Model Composer test bench and the record-and-replay wrapper with its test-signal generator.

Everything that names a block, a port, a type, a parameter or a state on the sheets is read from the files of ``modelcomposer/``
(``bt_psd_mc.v``, ``sig_gen.v``, ``bt_psd_selftest.v``, ``bt_psd_mc_config.m``, ``bt_mc_params.m``, ``bt_mc_build.m``) by
:func:`facts`, which also asserts that the pieces agree with each other.  The curves of sheet 31 come from the bit-accurate model
(``btpsd.siggen``, ``btpsd.mcmodel``, ``btpsd.estimator``) on the self-test record; the cycle counts on its timeline are the
model's PREDICTION (``btpsd.mcmodel.expected_cycles``), drawn hatched, and not a simulation result.  Nothing on either sheet is
taken from the paper, and nothing was executed in MATLAB, Simulink or Model Composer.
"""
from __future__ import annotations

import csv
import json
import re

import numpy as np
from matplotlib.patches import Rectangle

from _common import FIGS, RESULTS, ROOT
from _sheet import card, labels, legend
from btpsd import mcmodel as M
from btpsd import rtlparse, siggen
from btpsd.estimator import psd_exact_welch, psd_float_fd
from btpsd.rtlmodel import psd_engine_state_cycles
from btpsd.viz import INK, INK2, MUTED, PAL, Diagram, tint

MC = ROOT / "modelcomposer"
BLUE, ORANGE, AQUA, YELLOW, VIOLET = PAL["blue"], PAL["orange"], PAL["aqua"], PAL["yellow"], PAL["violet"]


# ---------------------------------------------------------------------------------------------------------------
# facts: what the sheets may say, read from the sources
# ---------------------------------------------------------------------------------------------------------------
def _strip_m(text: str) -> str:
    return "\n".join(re.sub(r"%.*", "", ln) for ln in text.splitlines())


def facts() -> dict:
    """Ports, types, parameters, blocks and states of the Model Composer pack, read from its files and cross-checked."""
    hdl = MC / "hdl"
    mc_ports = rtlparse.ports(hdl / "bt_psd_mc.v")
    gen_ports = rtlparse.ports(hdl / "sig_gen.v")
    st_ports = rtlparse.ports(hdl / "bt_psd_selftest.v")
    params_m = (MC / "matlab" / "bt_mc_params.m").read_text()
    p = {k: v.strip().strip("'") for k, v in re.findall(r"(?m)^p\.(\w+)\s*=\s*([^;]+);", params_m)}
    cfg = _strip_m((MC / "matlab" / "bt_psd_mc_config.m").read_text())
    in_order = re.findall(r"addSimulinkInport\('(\w+)'\)", cfg)
    out_order = re.findall(r"addSimulinkOutport\('(\w+)'\)", cfg)
    out_type = {n: t for n, t in re.findall(r"(\w+)_port\.setType\('(\w+)'\)", cfg)}
    in_width = {n: int(w) for n, w in re.findall(r"port\('(\w+)'\)\.width ~= (\d+)", cfg)}
    build = _strip_m((MC / "matlab" / "bt_mc_build.m").read_text())
    blocks = re.findall(r"add_blk\('([^']+)',\s*'(\w+)'", build)
    gates = re.findall(r"gw_in\('(\w+)',\s*rect\(\d, \d, \d+, \d+\),\s*'((?:[^']|'')+)',\s*([^,]+),\s*([^)]+)\)", build)
    wrapper = (hdl / "bt_psd_mc.v").read_text()
    edges = rtlparse.fsm_edges(hdl / "bt_psd_mc.v", "state")

    # consistency between the pieces ---------------------------------------------------------------------------------
    pm = {q.name: q for q in mc_ports}
    assert [q.name for q in mc_ports if q.direction == "input"][:2] == ["clk", "ce"]
    assert in_order == [q.name for q in mc_ports if q.direction == "input" and q.name not in ("clk", "ce")], in_order
    assert out_order == [q.name for q in mc_ports if q.direction == "output"], out_order
    assert in_width == {n: pm[n].width for n in in_width} and set(in_width) == set(in_order)
    assert pm["in_data"].width == int(p["wl"]) and pm["rd_data"].width == int(p["accw"]) and pm["rd_addr"].width == 10
    assert p["in_type"] == f"Fix_{p['wl']}_{int(p['wl']) - 1}" and p["out_type"] == f"UFix_{p['accw']}_0"
    assert out_type == {"busy": "UFix_1_0", "done": "UFix_1_0", "rd_data": p["out_type"]}, out_type
    assert int(p["nsamp"]) == 512 * int(p["nchunks"]) and int(p["segments"]) == int(p["nchunks"]) - 1
    assert {q.name: q.width for q in gen_ports if q.direction == "output"} == {"x": int(p["wl"]), "valid": 1}
    assert [q.name for q in st_ports if q.direction == "input"] == ["clk", "rst", "pinc", "noise_en", "noise_shift", "rd_addr"]
    states = rtlparse.state_values(rtlparse._strip_comments(wrapper))
    order = sorted(states, key=states.get)
    assert order == ["S_CAPTURE", "S_PRIME", "S_REPLAY", "S_DRAIN", "S_DONE"], order
    assert set(edges) == {("S_CAPTURE", "S_PRIME"), ("S_PRIME", "S_REPLAY"), ("S_REPLAY", "S_DRAIN"), ("S_DRAIN", "S_DONE")}, edges
    names = {n for _, n in blocks}
    assert {"sine", "noise", "sum", "siggen", "chirp", "stim", "one", "zero", "addr_src", "bt_psd_mc"} <= names, names
    gate = {g: (a, b.strip(), c.strip()) for g, a, b, c in gates}
    assert gate == {"gw_in_data": ("Signed  (2''s comp)", "p.wl", "p.wl - 1"), "gw_in_valid": ("Boolean", "1", "0"),
                    "gw_in_rst": ("Boolean", "1", "0"), "gw_in_addr": ("Unsigned", "10", "0")}, gate
    assert "xbsIndex_r4/Black Box" in {s for s, _ in blocks} and "xbsIndex_r4/Gateway Out" in build
    log_vars = re.findall(r"'(bt_\w+)'", " ".join(re.findall(r"taps = \{.*?\};", build, re.S)))
    assert log_vars == ["bt_busy", "bt_done", "bt_rd_data", "bt_in", "bt_valid", "bt_addr"], log_vars
    return dict(mc_ports=mc_ports, gen_ports=gen_ports, st_ports=st_ports, params=p, in_order=in_order, out_order=out_order,
                out_type=out_type, edges=edges, order=order, log_vars=log_vars, blocks=blocks)


# ---------------------------------------------------------------------------------------------------------------
# small drawing helpers
# ---------------------------------------------------------------------------------------------------------------
def _wire(D, pts, label=None, **kw):
    D.arrow(pts, label=label, color=INK2, lw=1.1, **kw)


def _tap(D, pts, label=None, **kw):
    D.arrow(pts, label=label, color=VIOLET, lw=0.85, ls=(0, (3.2, 2.0)), **kw)


def _port(D, x, y, name, left=True, fs=5.6):
    D.text(x + (1.2 if left else -1.2), y, name, fs=fs, color=INK, ha="left" if left else "right", mono=True, bold=True)


# ---------------------------------------------------------------------------------------------------------------
# fig30: the test bench
# ---------------------------------------------------------------------------------------------------------------
def fig30_test_bench():
    F = facts()
    p = F["params"]
    wl = int(p["wl"])
    W, H = 190, 158
    D = Diagram(W, H, "Model Composer test bench of the PSD estimator",
                subtitle="Block names, ports, types and parameters are read from modelcomposer/ (bt_mc_build.m, bt_psd_mc_config.m, bt_mc_params.m, bt_psd_mc.v) by python/scripts/make_mc_figures.py.\n"
                         "The Simulink side was not executed where this repository was written (no MATLAB, Simulink or Model Composer there); the HDL behind the Black Box was simulated in Icarus Verilog.")

    # ============================== (a) the bench =====================================================================
    AY, AH = 58.0, 84.0
    card(D, 2, AY, 186, AH, "(a)", "Test bench built by bt_mc_build.m: one sample period = one clock cycle",
         "the double-precision source is converted to fixed point by a Gateway In, recorded by the Black Box, and the result is read back through Gateway Out blocks into the MATLAB workspace",
         plot=False, fs=7.6)
    D.legend(7, AY + 3.4, [("io", "Simulink block, double precision"), ("rom", "Gateway: double ↔ fixed point"),
                           ("compute", "HDL in the Black Box"), ("memory", "workspace variable"), ("control", "clock source")], fs=5.4, dx=37)

    BH_ = 7.6  # height of the small blocks
    ys = dict(rst=128.2, in_valid=118.6, in_data=94.0, rd_addr=69.4)  # port heights, in the order of the Black Box ports
    BBX, BBW = 98.0, 40.0
    BBY0, BBY1 = 64.0, 131.5

    # ---- Black Box ---------------------------------------------------------------------------------------------------
    D.box(BBX, BBY0, BBW, BBY1 - BBY0, "", kind="compute")
    D.text(BBX + BBW / 2, 111.0, "Black Box", fs=7.6, color=INK, bold=True)
    D.text(BBX + BBW / 2, 107.6, "bt_psd_mc", fs=6.8, color=INK, mono=True, bold=True)
    for k, ln in enumerate(["records the first", f"{p['nsamp']} samples ({p['nchunks']} chunks),", "replays them into",
                            "bt_psd_top and raises done"]):
        D.text(BBX + BBW / 2, 88.5 - 2.3 * k, ln, fs=5.1, color=INK2)
    D.text(BBX + BBW / 2, 79.2, "Simulation mode:", fs=5.0, color=INK2)
    D.text(BBX + BBW / 2, 76.9, "Vivado Simulator", fs=5.0, color=INK, bold=True)
    D.text(BBX + BBW / 2, 73.9, f"WL = {p['wl']}  NCHUNKS = {p['nchunks']}", fs=5.0, color=INK2, mono=True)
    D.text(BBX + BBW / 2, 71.3, f"ROUND = {p['round']}  SQ_SERIAL = {p['sq_serial']}", fs=5.0, color=INK2, mono=True)
    D.text(BBX + BBW / 2, 66.6, "bt_psd_mc_config.m, bt_psd_mc_bundle.v", fs=4.4, color=MUTED)
    for name in F["in_order"]:
        _port(D, BBX, ys[name], name)
    out_y = dict(busy=120.0, done=109.0, rd_data=94.0)
    for name in F["out_order"]:
        _port(D, BBX + BBW, out_y[name], name, left=False)

    # ---- inputs: source -> Gateway In -> port --------------------------------------------------------------------------
    GX, GW = 52.0, 34.0

    def gateway(y, title, sub, fs=5.8):
        D.box(GX, y - BH_ / 2, GW, BH_, f"{title}\n{sub}", kind="rom", fs=fs)

    for name, const, val in (("rst", "Constant 0", "no reset: the block powers up armed"), ("in_valid", "Constant 1", "every sample is valid")):
        y = ys[name]
        D.box(5, y - BH_ / 2, 34, BH_, f"{const}\n{val}", kind="io", fs=5.6)
        gateway(y, "Gateway In", "Boolean")
        _wire(D, [(39, y), (GX, y)])
        _wire(D, [(GX + GW, y), (BBX, y)], label="Boolean", lpos=0.5, loff=(0, 1.5), fs=4.8)
    y = ys["rd_addr"]
    D.box(5, y - BH_ / 2, 34, BH_, "From Workspace  bt_addr_stim\naddress ramp 0 … 1023, repeating", kind="io", fs=5.4)
    gateway(y, "Gateway In", "UFix_10_0")
    _wire(D, [(39, y), (GX, y)])
    _wire(D, [(GX + GW, y), (BBX, y)], label="UFix_10_0", lpos=0.5, loff=(0, 1.5), fs=4.8)

    y = ys["in_data"]
    srcs = [("Sine Wave + Random Number", "Sum · bt_mc_build('sine')"), ("Signal Generator", "sine, square, sawtooth · 'signal'"),
            ("Chirp Signal", "0.01 → 0.45 cycles/sample · 'chirp'"), ("From Workspace  bt_stim", "file written by Python · 'file'")]
    pitch = 9.0
    top = y + 1.5 * pitch
    cy = []
    for k, (a, b) in enumerate(srcs):
        yy = top - pitch * k
        D.box(5, yy - BH_ / 2, 34, BH_, f"{a}\n{b}", kind="io", fs=5.4)
        cy.append(yy)
        _wire(D, [(39, yy), (43.0, yy)], head=False)
    D.ax.plot([43.0, 43.0], [cy[0], cy[-1]], color=INK2, lw=1.1, zorder=2.5)
    D.junction(43.0, y)
    D.text(43.0, cy[0] + 5.2, "choose one source", fs=4.8, color=INK2, italic=True)
    gateway(y, "Gateway In", f"{p['in_type']} · round half away from zero · saturate", fs=5.8)
    D.ax.texts[-1].set_fontsize(5.0)
    _wire(D, [(43.0, y), (GX, y)])
    _wire(D, [(GX + GW, y), (BBX, y)], label=p["in_type"], lpos=0.5, loff=(0, 1.5), fs=4.8)
    D.text(GX + GW / 2, y - 6.2, "−1 ≤ x < 1, code = x · 2^%d" % (wl - 1), fs=4.8, color=INK2)

    # ---- outputs: Gateway Out -> To Workspace --------------------------------------------------------------------------
    OX = 146.0
    for name, var in zip(F["out_order"], F["log_vars"][:3]):
        y = out_y[name]
        _wire(D, [(BBX + BBW, y), (OX, y)], label=F["out_type"][name], lpos=0.5, loff=(0, 1.5), fs=4.8)
        D.box(OX, y - 3.6, 15.0, 7.2, "Gateway\nOut", kind="rom", fs=5.4)
        _wire(D, [(OX + 15.0, y), (OX + 21.0, y)])
        D.box(OX + 21.0, y - 3.6, 18.0, 7.2, f"To Workspace\n{var}", kind="memory", fs=5.2)
    D.text(OX + 19.5, 86.4, "also logged the same way:", fs=4.8, color=INK2)
    D.text(OX + 19.5, 83.7, ", ".join(F["log_vars"][3:]), fs=5.0, color=INK, mono=True, bold=True)
    D.text(OX + 19.5, 81.1, "(outputs of the three Gateway In blocks)", fs=4.6, color=INK2)

    # ---- hub -----------------------------------------------------------------------------------------------------------
    D.box(OX, 66.0, 40.0, 9.0, "Model Composer Hub\nsystem period 1 · clk, ce implicit · single rate", kind="control", fs=5.4)
    _tap(D, [(OX, 70.5), (BBX + BBW, 70.5)], "clk, ce", lpos=0.5, loff=(0, 1.5), fs=4.8)

    # ============================== (b) check flow ====================================================================
    BY, BH = 2.0, 54.0
    card(D, 2, BY, 100, BH, "(b)", "Checking a run against the bit-exact model",
         "the expected spectrum is computed from the logged Gateway In codes: it holds for every source block and every gateway rounding",
         plot=False, fs=7.6)
    steps = [("export_mc_stimulus.py", "exact doubles\ncode / 2^%d" % (wl - 1), "bt_stim.mat\nbt_stim.csv", "io"),
             ("bt_mc_build.m", "builds bench (a)\nwith the chosen\nsource", "bt_psd_mc_demo\n.slx", "io"),
             ("bt_mc_run.m", "simulates until done\nhas been high for\n1100 cycles", "bt_log.mat\n(-v7)", "io"),
             ("compare_mc_result.py", "decodes rd_data,\nbit-exact check,\ncycle prediction", "PASS or FAIL\nexit status 0/1/2", "compute")]
    bw, gap, x0 = 20.0, 4.2, 5.0
    for k, (title, body, out, kind) in enumerate(steps):
        x = x0 + k * (bw + gap)
        D.box(x, 27.0, bw, 16.0, f"{title}\n{body}", kind=kind, fs=5.4, step=k + 1)
        D.box(x + 1.5, 15.5, bw - 3.0, 7.5, out, kind="memory", bold=False, fs=5.0)
        _wire(D, [(x + bw / 2, 27.0), (x + bw / 2, 23.0)], head=True)
        if k < 3:
            _wire(D, [(x + bw, 35.0), (x + bw + gap, 35.0)])
    D.text(x0, 11.0, "A run passes when all 1 024 words read through rd_addr / rd_data equal psd_fixed(logged codes), bit for bit.", fs=5.0, ha="left", color=INK)
    D.text(x0, 8.0, "The cycles from the first sample to done are compared with btpsd.mcmodel.expected_cycles and reported;", fs=5.0, ha="left", color=INK2)
    D.text(x0, 5.0, "a difference there is information, not a failure (the prediction is a model, see sheet 31 (b)).", fs=5.0, ha="left", color=INK2)

    # ============================== (c) interface table ===============================================================
    CX, CW = 104.0, 84.0
    card(D, CX, BY, CW, BH, "(c)", "Black Box interface", "port order = order of addSimulinkInport / addSimulinkOutport in bt_psd_mc_config.m",
         plot=False, fs=7.6)
    cols = [CX + 3.0, CX + 15.5, CX + 24.5, CX + 40.0]
    y0 = BY + BH - 11.5
    for c, t in zip(cols, ["port", "dir", "Simulink type", "meaning"]):
        D.text(c, y0, t, fs=5.4, ha="left", bold=True, color=INK)
    D.ax.plot([CX + 2.0, CX + CW - 2.0], [y0 - 1.6, y0 - 1.6], color=MUTED, lw=0.6)
    simty = {"rst": "Boolean", "in_valid": "Boolean", "in_data": p["in_type"], "rd_addr": "UFix_10_0", **F["out_type"]}
    mean = {"clk": "implicit pair with ce (Hub)", "ce": "qualifies in_valid only", "rst": "synchronous, active high",
            "in_valid": "1 for every sample", "in_data": f"{wl}-bit two's complement", "rd_addr": "bin 0 … 1023 (two-sided)",
            "busy": "from first sample to done", "done": "spectrum complete, stays", "rd_data": "word of bin rd_addr, 1 clock later"}
    rows = [("clk", "in"), ("ce", "in")] + [(n, "in") for n in F["in_order"]] + [(n, "out") for n in F["out_order"]]
    byname = {q.name: q for q in F["mc_ports"]}
    for r, (name, d) in enumerate(rows):
        yy = y0 - 4.3 - 3.7 * r
        assert byname[name].direction == ("input" if d == "in" else "output")
        implicit = name in ("clk", "ce")
        if r % 2 == 0:
            D.ax.add_patch(Rectangle((CX + 2.0, yy - 1.8), CW - 4.0, 3.6, fc="#f4f3ef", ec="none", zorder=0.7))
        D.text(cols[0], yy, name, fs=5.2, ha="left", mono=True, bold=True, color=INK)
        D.text(cols[1], yy, d, fs=5.2, ha="left", color=INK2)
        D.text(cols[2], yy, "(implicit)" if implicit else simty[name], fs=5.2, ha="left", mono=True, color=INK2)
        D.text(cols[3], yy, mean[name], fs=5.0, ha="left", color=INK2)
    D.save(FIGS / "fig30_mc_test_bench")
    return D


# ---------------------------------------------------------------------------------------------------------------
# fig31: wrapper, timeline, generator
# ---------------------------------------------------------------------------------------------------------------
def wrapper_facts() -> dict:
    """The conditions of the wrapper's state machine and the constants of the generator, read from the Verilog and cross-checked."""
    w = rtlparse._strip_comments((MC / "hdl" / "bt_psd_mc.v").read_text())
    g = rtlparse._strip_comments((MC / "hdl" / "sig_gen.v").read_text())
    sel = rtlparse._strip_comments((MC / "hdl" / "bt_psd_selftest.v").read_text())
    assert re.search(r"if \(wr_cnt == NSAMP - 1\) state <= S_PRIME", w)
    assert re.search(r"S_PRIME:\s+state <= S_REPLAY", w)
    assert re.search(r"S_REPLAY:\s+if \(accept && rp == NSAMP - 1\) state <= S_DRAIN", w)
    assert re.search(r"S_DRAIN:\s+if \(!core_busy\) state <= S_DONE", w)
    assert re.search(r"wire cap = \(state == S_CAPTURE\) && ce && in_valid", w)
    assert re.search(r"assign busy = \(state == S_CAPTURE\) \? \(wr_cnt != 0\) : \(state != S_DONE\)", w)
    assert re.search(r"assign done = \(state == S_DONE\)", w)
    assert re.search(r"reg \[3:0\] rst_sr = 4'b0000", w) and re.search(r"wire core_rst_n = &rst_sr", w)
    assert re.search(r"localparam \[4:0\]\s+AVG_SHIFT = clog2\(NCHUNKS - 1\)", w)
    assert re.search(r"reg \[WL-1:0\] ram \[0:NSAMP-1\]", w)
    taps = int(re.search(r"32'h([0-9A-Fa-f]{8}) : 32'h0", g).group(1), 16)
    seed = int(re.search(r"SEED = 32'h([0-9A-Fa-f]{8})", g).group(1), 16)
    assert (taps, seed) == (siggen.TAPS, siggen.SEED)
    assert re.search(r"phase <= phase \+ pinc", g) and re.search(r"nz_raw >>> noise_shift", g) and re.search(r"u_gen \(", sel)
    assert re.search(r"\.rst_n\(~rst\), \.en\(1'b1\)", sel) and re.search(r"\.in_valid\(valid\), \.in_data\(x\)", sel)
    rom = (MC / "hdl" / "sig_rom.v").read_text()
    assert len(re.findall(r"mem\[\d+\] = ", rom)) == 1 << siggen.ROM_BITS
    return dict(taps=taps, seed=seed)


def timeline(wl: int = 12, nchunks: int = 5, tone_bin: float = siggen.SELFTEST["bin_1024"], noise_shift: int = siggen.SELFTEST["noise_shift"]):
    """Predicted phases of one run of the Black Box on a generator record (cycles from the first sample to ``done``)."""
    cfg = M.McConfig(wl=wl, nchunks=nchunks)
    pinc = siggen.pinc_for(tone_bin)
    x = siggen.generate(cfg.nsamp, pinc, noise_shift, wl)
    res = M.expected(x, cfg)
    pred = M.expected_cycles(res, cfg)
    st = psd_engine_state_cycles(res.S, res.squarer_cycles, clears=0)
    merge = sum(st[k] for k in ("M_A", "M_B", "M_C", "M_F"))
    window = sum(st[k] for k in ("W_A", "W_B", "W_C", "W_F", "W_G", "W_ACC"))
    squarer = st["W_S"]
    other = st["S_END"] + res.S * M.T_SEG_EXTRA + 2
    phases = [("capture", cfg.nsamp), ("replay: load", nchunks * M.CHUNK), ("replay: FFT", nchunks * M.T_FFT_PER_CHUNK),
              ("segment merge", merge), ("window, accumulate", window), ("squarer", squarer), ("handshakes", other)]
    assert sum(v for _, v in phases) == pred, (sum(v for _, v in phases), pred)
    return dict(cfg=cfg, res=res, pred=pred, phases=phases, x=x, pinc=pinc)


def simulated_cycles():
    """The clock cycles that the CI simulation printed (``results/rtl_mc_cycles.csv``) and the id of the run that printed them."""
    with (RESULTS / "rtl_mc_cycles.csv").open() as f:
        rows = {(r["testbench"], int(r["run"])): int(r["cycles"]) for r in csv.DictReader(f)}
    run_id = json.loads((RESULTS / "rtl_sim_provenance.json").read_text())["run_id"]
    return rows, run_id


def _title(ax, s):
    ax.set_title(s, fontsize=5.2, color=INK2, loc="left", pad=2.4)


def _n(v: int) -> str:
    """Integer with a thin-space-like separator for thousands (a plain space, as in the text of the repository)."""
    return f"{v:,}".replace(",", " ")


def fig31_wrapper_generator():
    F = facts()
    WF = wrapper_facts()
    p = F["params"]
    wl = int(p["wl"])
    T = timeline(wl, int(p["nchunks"]))
    sim, run_id = simulated_cycles()
    prd = M.testbench_cycle_predictions()
    k0, k1, ks = ("tb_bt_psd_mc", 0), ("tb_bt_psd_mc", 1), ("tb_bt_psd_selftest", 0)
    W, H = 190, 210
    D = Diagram(W, H, "The record-and-replay wrapper, its run, and the test-signal generator",
                subtitle="State machine, conditions, constants and ports are read from modelcomposer/hdl/*.v; the curves are the bit-accurate model on synthetic records (no measurement).\n"
                         "The cycle counts in (b) are the model's PREDICTION (hatched), not a simulation result; tb/tb_bt_psd_mc.v prints the count that the RTL simulation takes.\n"
                         f"On the testbench records (Icarus Verilog, CI run {run_id}), cycles simulated / predicted: tb_bt_psd_mc {_n(sim[k0])} / {_n(prd[k0])} and "
                         f"{_n(sim[k1])} / {_n(prd[k1])}, tb_bt_psd_selftest {_n(sim[ks])} / {_n(prd[ks])}.")

    # ============================== (a) wrapper internals ===========================================================
    AY, AH = 134.0, 58.0
    card(D, 2, AY, 186, AH, "(a)", "bt_psd_mc: record the samples, replay them into the estimator core",
         "the Simulink source cannot be held off and the core stalls for its FFT and segment merge: the wrapper records first and feeds the core at the pace it accepts",
         plot=False, fs=7.6)
    # --- control row: the state machine ---------------------------------------------------------------------------
    fy = AY + 41.0
    bw, gap, x0 = 17.0, 10.5, 5.0
    cond = ["wr_cnt = NSAMP - 1", "always", "accept & rp = NSAMP - 1", "!core_busy"]
    sx = [x0 + k * (bw + gap) for k in range(len(F["order"]))]
    for k, n in enumerate(F["order"]):
        D.box(sx[k], fy, bw, 6.0, n, kind="control", fs=4.9)
    for k, c in enumerate(cond):
        xa, xb = sx[k] + bw, sx[k + 1]
        D.arrow([(xa, fy + 3.0), (xb, fy + 3.0)], color=VIOLET, lw=0.9)
        D.text((xa + xb) / 2, fy + 7.0, c, fs=4.0, color=INK2)
    D.text(138.0, AY + 45.6, "busy = CAPTURE ? (wr_cnt ≠ 0) : (state ≠ S_DONE)", fs=4.5, color=INK, ha="left", mono=True)
    D.text(138.0, AY + 43.2, "done = (state = S_DONE)", fs=4.5, color=INK, ha="left", mono=True)
    D.text(138.0, AY + 40.8, "rst (synchronous, active high): back to S_CAPTURE", fs=4.4, color=INK2, ha="left")
    # --- data row ------------------------------------------------------------------------------------------------------
    D.text(3.0, AY + 29.2, "in_valid, ce", fs=4.8, color=INK, ha="left", mono=True, bold=True)
    D.arrow([(3.0, AY + 26.5), (12.0, AY + 26.5)], color=INK2, lw=1.1)
    D.box(12.0, AY + 22.5, 23.0, 8.0, "capture gate\ncap = CAPTURE & ce & in_valid", kind="control", fs=4.8)
    _tap(D, [(sx[0] + bw / 2, fy), (sx[0] + bw / 2, AY + 30.5)], None)
    D.text(sx[0] + bw / 2 + 1.0, AY + 36.4, "state = S_CAPTURE", fs=4.2, color=VIOLET, ha="left")
    D.text(3.0, AY + 21.0, f"in_data [{wl}]", fs=4.8, color=INK, ha="left", mono=True, bold=True)
    D.arrow([(3.0, AY + 18.5), (44.0, AY + 18.5)], color=INK2, lw=1.1)
    D.arrow([(35.0, AY + 26.5), (44.0, AY + 26.5)], label="write enable", lpos=0.5, loff=(0, 1.5), fs=4.2, color=VIOLET, lw=0.85, ls=(0, (3.2, 2.0)))
    D.box(44.0, AY + 14.5, 32.0, 17.5, f"ram\n{int(p['nsamp'])} × {wl} bit\nwrite port: wr_cnt\nread port: ra, registered", kind="memory", fs=5.0)
    D.arrow([(76.0, AY + 23.0), (84.0, AY + 23.0)], color=INK2, lw=1.1)
    D.box(84.0, AY + 20.0, 14.0, 6.0, "rdata", kind="memory", fs=5.0, bold=False)
    D.arrow([(98.0, AY + 23.0), (112.0, AY + 23.0)], label=f"in_data [{wl}]", lpos=0.5, loff=(0, 1.5), fs=4.2, color=INK2, lw=1.1)
    D.text(60.0, AY + 11.6, "ra = (rp_next ≥ NSAMP) ? 0 : rp_next", fs=4.3, color=INK2)
    D.text(60.0, AY + 9.3, "rp_next = accept ? rp + 1 : rp", fs=4.3, color=INK2)
    # the core
    cx, cw = 112.0, 40.0
    D.box(cx, AY + 3.0, cw, 33.0, "", kind="compute")
    D.text(cx + cw / 2, AY + 33.2, "bt_psd_top  (u_core)", fs=6.6, color=INK, bold=True)
    D.text(cx + cw / 2, AY + 30.6, "chunk-merge Blackman-Tukey PSD", fs=4.6, color=INK2)
    D.text(cx + cw / 2, AY + 28.6, "radix-8 FFT-512, six-tap merge", fs=4.6, color=INK2)
    left = {"in_valid": 27.0, "in_data": 23.0, "rst_n": 17.0, "avg_shift": 10.5, "res_addr": 5.5}
    for lab, yy in left.items():
        _port(D, cx, AY + yy, lab, fs=4.8)
    for lab, yy in (("in_ready", 27.0), ("busy", 23.0), ("res_data", 5.5)):
        _port(D, cx + cw, AY + yy, lab, left=False, fs=4.8)
    # core_in_valid comes from the state machine
    xr = sx[2] + bw / 2
    _tap(D, [(xr, fy), (xr, AY + 37.6), (106.0, AY + 37.6), (106.0, AY + 27.0), (cx, AY + 27.0)])
    D.text(xr + 1.2, AY + 38.8, "core_in_valid = (state = S_REPLAY)", fs=4.2, color=VIOLET, ha="left")
    # reset generator and avg_shift
    D.box(80.0, AY + 13.0, 26.0, 6.5, "rst_sr[3:0] ← {rst_sr[2:0], ~rst}\ncore_rst_n = &rst_sr", kind="free", fs=4.4, bold=False)
    D.arrow([(106.0, AY + 17.0), (cx, AY + 17.0)], color=INK2, lw=1.0)
    D.box(84.0, AY + 8.4, 22.0, 4.0, f"avg_shift = clog2(NCHUNKS - 1) = {int(p['avg_shift'])}", kind="io", fs=4.0, bold=False)
    D.arrow([(106.0, AY + 10.5), (cx, AY + 10.5)], color=INK2, lw=1.0)
    # rd_addr straight through, outputs
    D.text(3.0, AY + 8.0, "rd_addr [10]", fs=4.8, color=INK, ha="left", mono=True, bold=True)
    D.arrow([(3.0, AY + 5.5), (cx, AY + 5.5)], color=INK2, lw=1.1)
    D.text(46.0, AY + 3.4, "straight through: rd_data follows rd_addr by one clock", fs=4.4, color=INK2)
    D.arrow([(cx + cw, AY + 5.5), (162.0, AY + 5.5)], color=INK2, lw=1.1)
    D.box(162.0, AY + 2.5, 22.0, 6.0, f"rd_data [{int(p['accw'])}]", kind="io", fs=5.0, bold=False)
    D.text(155.0, AY + 27.0, "accept = core_in_valid & in_ready", fs=4.3, color=VIOLET, ha="left", mono=True)
    D.text(155.0, AY + 23.0, "core_busy = busy (S_DRAIN waits for 0)", fs=4.3, color=VIOLET, ha="left", mono=True)
    D.text(155.0, AY + 17.0, "bundle: every module in one file,", fs=4.4, color=INK2, ha="left")
    D.text(155.0, AY + 15.0, "no `include (Black Box import)", fs=4.4, color=INK2, ha="left")

    # ============================== (b) timeline ====================================================================
    BY, BH = 88.0, 44.0
    pred = T["pred"]
    nsq = T["res"].squarer_cycles
    card(D, 2, BY, 186, BH, "(b)", f"One run with {int(p['nchunks'])} chunks: predicted clock cycles from the first sample to done",
         f"generator record (tone at bin {siggen.SELFTEST['bin_1024']:g}, noise_shift {siggen.SELFTEST['noise_shift']}): {_n(pred)} cycles; "
         f"the squarer needs 3 x 1 024 cycles per segment plus {_n(nsq)} data-dependent cycles for this record",
         plot=False, fs=7.6)
    ax = D.inset(16.0, BY + 9.0, 160.0, 24.5, fs=4.6, grid=False)
    tend = pred + 1500
    cols = [BLUE, AQUA, VIOLET, ORANGE, YELLOW, PAL["magenta"], MUTED]
    ybar, hbar = 6.3, 1.2
    t0 = 0
    for (name, n), col in zip(T["phases"], cols):
        ax.add_patch(Rectangle((t0, ybar), n, hbar, fc=tint(col, 0.55), ec=col, lw=0.7, hatch="////", zorder=3))
        t0 += n
    n_cap = T["cfg"].nsamp

    def trace(y, xs, ys_, label, col):
        ax.step(xs, [y + v * 0.7 for v in ys_], where="post", color=col, lw=0.9, zorder=3)
        ax.text(-0.012 * tend, y + 0.35, label, ha="right", va="center", fontsize=4.6, color=INK, family="DejaVu Sans Mono", fontweight="bold")

    trace(4.4, [0, n_cap, tend], [1, 0, 0], "in_valid", INK2)
    trace(3.2, [0, 1, pred, tend], [0, 1, 0, 0], "busy", BLUE)
    trace(2.0, [0, pred, tend], [0, 1, 1], "done", AQUA)
    trace(0.8, [0, pred, tend], [0, 0, 0], "rd_addr", VIOLET)
    ax.plot([pred, pred + 1024, pred + 1024, pred + 1100], [0.8, 1.5, 0.8, 0.87], color=VIOLET, lw=0.9, zorder=3)
    ax.set_xlim(0, tend)
    ax.set_ylim(0.3, 8.3)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ticks = [0, 20000, 40000, 60000, 80000, 100000, 120000, pred]
    ax.set_xticks(ticks)
    ax.set_xticklabels(["0"] + [f"{t // 1000} k" for t in ticks[1:-1]] + [_n(pred)], fontsize=4.4)
    ax.axvline(pred, color=AQUA, lw=0.6, ls=(0, (2, 1.5)), zorder=1)
    labels(ax, "clock cycles after the first sample (in_valid = 1 in every cycle)", None, fs=4.8)
    t0 = 0
    for (name, n), col in zip(T["phases"], cols):
        mid = t0 + n / 2
        t0 += n
        if n > 0.08 * tend:
            ax.text(mid, 7.85, f"{name}: {_n(n)}", ha="center", va="center", fontsize=4.6, color=INK)
    small = [(nm, n) for nm, n in T["phases"] if n <= 0.08 * tend]
    D.text(16.0, BY + 3.6, "short phases at the left: " + ", ".join(f"{nm} {_n(n)}" for nm, n in small), fs=4.6, color=INK, ha="left")
    D.text(176.0, BY + 3.6, "read-out: rd_addr counts 0 … 1 023 after done (at least 1 025 cycles)", fs=4.6, color=INK2, ha="right")
    D.text(176.0, BY + BH - 3.6, "traces idealised from the prediction", fs=4.4, color=MUTED, ha="right", italic=True)

    # ============================== (c) generator and self-test =====================================================
    CY, CH = 44.0, 42.0
    pinc = siggen.pinc_for(siggen.SELFTEST["bin_1024"])
    card(D, 2, CY, 186, CH, "(c)", "In-fabric test-signal generator and self-test (bt_psd_selftest = sig_gen + bt_psd_mc)",
         "plays the part of Simulink's Sine Wave + Random Number sources as plain Verilog; it is not the AMD DDS Compiler IP", plot=False, fs=7.6)
    D.ax.add_patch(Rectangle((4.0, CY + 2.0), 178.0, CH - 12.0, fc="none", ec=MUTED, lw=0.8, ls="--", zorder=0.9))
    D.text(6.0, CY + CH - 11.0, "bt_psd_selftest", fs=5.2, color=INK2, ha="left", mono=True, italic=True)
    gy = CY + 20.0
    D.box(8.0, gy, 26.0, 8.0, "phase accumulator\nphase += pinc\nPHW = 32", kind="control", fs=5.0)
    D.arrow([(34.0, gy + 4.0), (40.0, gy + 4.0)], label="phase[31:22]", lpos=0.5, loff=(0, 1.6), fs=4.2, color=INK2, lw=1.1)
    D.box(40.0, gy, 24.0, 8.0, f"sig_rom\n{1 << siggen.ROM_BITS} × {wl} bit\nsine, amplitude 0.5", kind="rom", fs=5.0)
    D.arrow([(64.0, gy + 4.0), (74.1, gy + 4.0)], label="tone", lpos=0.5, loff=(0, 1.6), fs=4.4, color=INK2, lw=1.1)
    D.adder(76.0, gy + 4.0, "+", r=1.9)
    D.box(8.0, gy - 10.5, 26.0, 8.0, f"LFSR  32 bit\nx^32+x^22+x^2+x+1\nseed {WF['seed']:08X}", kind="control", fs=5.0)
    D.arrow([(34.0, gy - 6.5), (40.0, gy - 6.5)], label=f"top {wl}", lpos=0.5, loff=(0, 1.6), fs=4.2, color=INK2, lw=1.1)
    D.box(40.0, gy - 10.5, 24.0, 8.0, "noise\n>>> noise_shift\nnoise_en gate", kind="compute", fs=5.0)
    D.arrow([(64.0, gy - 6.5), (76.0, gy - 6.5), (76.0, gy + 2.0)], label="noise", lpos=0.4, loff=(0, 1.6), fs=4.4, color=INK2, lw=1.1)
    D.box(82.0, gy, 22.0, 8.0, "saturating sum\nx [%d]\nvalid after 2 edges" % wl, kind="compute", fs=5.0)
    D.arrow([(77.9, gy + 4.0), (82.0, gy + 4.0)], color=INK2, lw=1.1)
    D.arrow([(104.0, gy + 4.0), (116.0, gy + 4.0)], label=f"x [{wl}], valid", lpos=0.5, loff=(0, 1.6), fs=4.4, color=INK2, lw=1.1)
    D.box(116.0, gy - 5.0, 30.0, 17.0, "bt_psd_mc\nrecord → replay → done\n%d chunks" % int(p["nchunks"]), kind="compute", fs=5.4)
    D.arrow([(146.0, gy + 4.0), (156.0, gy + 4.0)], label="busy, done", lpos=0.5, loff=(0, 1.6), fs=4.4, color=INK2, lw=1.1)
    D.text(169.0, gy + 4.7, "rd_addr → rd_data", fs=4.8, color=INK, ha="center", mono=True, bold=True)
    D.text(169.0, gy + 1.9, "(1 024 words, 1 clock later)", fs=4.4, color=INK2)
    D.text(8.0, CY + 6.6, f"tone = pinc / 2^32 cycles per sample, i.e. bin = pinc · 1 024 / 2^32; the self-test uses pinc = {_n(pinc)} (bin {siggen.tone_bin(pinc):.2f}), noise_shift = {siggen.SELFTEST['noise_shift']}",
           fs=4.7, color=INK, ha="left")
    D.text(8.0, CY + 4.0, "inputs: pinc [32], noise_en, noise_shift [4], rd_addr [10], rst, clk.  A one-cycle pulse on rst restarts generator and estimator; the estimator also starts by itself at power-up.",
           fs=4.7, color=INK2, ha="left")

    # ============================== (d), (e) the record and its spectrum ==============================================
    DY, DH = 2.0, 40.0
    rec = siggen.selftest_record()
    cfg4 = M.McConfig(wl=rec["wl"], nchunks=rec["nchunks"])
    x = rec["x"]
    tone_only = siggen.generate(len(x), rec["pinc"], None, rec["wl"])
    res = M.expected(x, cfg4)
    P = np.maximum(res.P_int.astype(float), 1.0)
    ref = psd_exact_welch(x.astype(float))
    fd = psd_float_fd(x.astype(float))
    f = siggen.tone_bin(rec["pinc"])
    card(D, 2, DY, 90, DH, "(d)", "The self-test record",
         f"first 160 samples, WL = {rec['wl']}: tone (amplitude 0.5 of full scale) plus LFSR noise shifted by {rec['noise_shift']}", plot=False, fs=7.6)
    a1 = D.inset(12.0, DY + 7.0, 76.0, 22.0, fs=4.6)
    n = np.arange(160)
    a1.plot(n, tone_only[:160], color=MUTED, lw=0.7, label="tone only", zorder=2)
    a1.plot(n, x[:160], color=BLUE, lw=0.8, label="tone + noise (the record)", zorder=3)
    a1.set_xlim(0, 159)
    a1.set_ylim(-1700, 1700)
    labels(a1, "sample n", "code", fs=4.8)
    legend(a1, "upper right", fs=4.4)
    D.card(94.0, DY, 94.0, DH)
    D.panel(95.8, DY + DH - 3.3, "(e)", "Its spectrum", fs=7.6)
    D.text(95.8, DY + DH - 6.6, f"{rec['nchunks']} chunks (3 segments), Hann window, dB relative to each curve's peak; one-sided half of the 1 024 bins", fs=5.6, ha="left", color=INK2)
    a2 = D.inset(106.0, DY + 7.0, 78.0, 22.0, fs=4.6)
    k = np.arange(512)
    db = lambda v: 10.0 * np.log10(np.maximum(v / np.max(v), 1e-12))  # noqa: E731
    a2.plot(k, db(ref[:512]), color=INK, lw=0.6, label="exact Welch average", zorder=2)
    a2.plot(k, db(fd[:512]), color=ORANGE, lw=0.9, label="float six-tap merge", zorder=3)
    a2.plot(k, db(P[:512]), color=BLUE, lw=0.7, label="model words (rd_data)", zorder=4)
    a2.axvline(f, color=MUTED, lw=0.6, ls=(0, (2.2, 1.6)), zorder=1)
    a2.text(f + 6, -2.0, f"tone, bin {f:.2f}", fontsize=4.3, color=INK2, va="top")
    a2.set_xlim(0, 511)
    a2.set_ylim(-70, 3)
    labels(a2, "bin of 1 024", "dB", fs=4.8)
    legend(a2, "upper right", fs=4.2)
    D.save(FIGS / "fig31_mc_wrapper_generator")
    return D


# ---------------------------------------------------------------------------------------------------------------
FIGURES = [fig30_test_bench, fig31_wrapper_generator]


def main():
    for fn in FIGURES:
        fn()
        print("wrote", fn.__name__)


if __name__ == "__main__":
    main()
