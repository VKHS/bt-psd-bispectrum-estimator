"""What the Model Composer sheets (fig30, fig31) draw and print is what the sources of ``modelcomposer/`` say.

``facts()`` and ``wrapper_facts()`` in ``make_mc_figures.py`` read the Black Box configuration, the build script, the parameters
and the Verilog and assert that they agree with each other; here that assertion is shown to bite when one of the pieces is
changed, and the numbers on the timeline are tied to ``btpsd.mcmodel``.
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python" / "scripts"))

import make_mc_figures as mf  # noqa: E402
from btpsd import mcmodel as M  # noqa: E402
from btpsd import siggen  # noqa: E402


@pytest.fixture()
def pack(tmp_path, monkeypatch):
    """A scratch copy of ``modelcomposer/`` that the figure code reads instead of the real one."""
    dst = tmp_path / "modelcomposer"
    shutil.copytree(ROOT / "modelcomposer", dst)
    monkeypatch.setattr(mf, "MC", dst)
    return dst


def _edit(path: Path, old: str, new: str):
    t = path.read_text()
    assert t.count(old) == 1, (path.name, old)
    path.write_text(t.replace(old, new))


def _texts(D):
    out = [t.get_text() for t in D.ax.texts]
    for ax in D.fig.axes:
        out += [t.get_text() for t in ax.texts] + [ax.get_xlabel(), ax.get_ylabel(), ax.get_title(loc="left")]
        if ax.get_legend() is not None:
            out += [t.get_text() for t in ax.get_legend().get_texts()]
    return "\n".join(out)


# ---------------------------------------------------------------------------------------------------------------
# facts
# ---------------------------------------------------------------------------------------------------------------
def test_facts_describe_the_committed_pack():
    F = mf.facts()
    assert F["in_order"] == ["rst", "in_valid", "in_data", "rd_addr"] and F["out_order"] == ["busy", "done", "rd_data"]
    assert F["out_type"] == {"busy": "UFix_1_0", "done": "UFix_1_0", "rd_data": "UFix_32_0"}
    assert F["order"] == ["S_CAPTURE", "S_PRIME", "S_REPLAY", "S_DRAIN", "S_DONE"]
    assert F["log_vars"] == ["bt_busy", "bt_done", "bt_rd_data", "bt_in", "bt_valid", "bt_addr"]
    p = F["params"]
    assert (p["wl"], p["nchunks"], p["nsamp"], p["segments"], p["avg_shift"], p["accw"]) == ("12", "5", "2560", "4", "2", "32")
    assert {q.name: q.width for q in F["mc_ports"]}["rd_data"] == 32
    assert [q.name for q in F["gen_ports"] if q.direction == "output"] == ["x", "valid"]
    wf = mf.wrapper_facts()
    assert (wf["taps"], wf["seed"]) == (siggen.TAPS, siggen.SEED)


def test_swapping_two_black_box_ports_in_the_configuration_is_detected(pack):
    cfg = pack / "matlab" / "bt_psd_mc_config.m"
    t = cfg.read_text()
    a, b = "  this_block.addSimulinkInport('in_valid');\n", "  this_block.addSimulinkInport('in_data');\n"
    assert a in t and b in t
    cfg.write_text(t.replace(a + b, b + a))
    with pytest.raises(AssertionError):
        mf.facts()


def test_a_changed_output_type_is_detected(pack):
    _edit(pack / "matlab" / "bt_psd_mc_config.m", "rd_data_port.setType('UFix_32_0')", "rd_data_port.setType('UFix_31_0')")
    with pytest.raises(AssertionError):
        mf.facts()


def test_a_parameter_file_that_disagrees_with_the_hdl_is_detected(pack):
    _edit(pack / "matlab" / "bt_mc_params.m", "p.nchunks = 5;", "p.nchunks = 6;")
    with pytest.raises(AssertionError):
        mf.facts()


def test_a_wider_data_port_than_the_gateway_type_is_detected(pack):
    _edit(pack / "hdl" / "bt_psd_mc.v", "input  signed [WL-1:0]     in_data,", "input  signed [WL:0]       in_data,")
    with pytest.raises(AssertionError):
        mf.facts()


def test_a_renamed_log_variable_is_detected(pack):
    _edit(pack / "matlab" / "bt_mc_build.m", "'bt_rd_data'", "'bt_rd_word'")
    with pytest.raises(AssertionError):
        mf.facts()


def test_a_changed_state_machine_condition_is_detected(pack):
    _edit(pack / "hdl" / "bt_psd_mc.v", "S_DRAIN:   if (!core_busy) state <= S_DONE;", "S_DRAIN:   if (core_busy) state <= S_DONE;")
    with pytest.raises(AssertionError):
        mf.wrapper_facts()


def test_a_changed_lfsr_polynomial_is_detected(pack):
    _edit(pack / "hdl" / "sig_gen.v", "32'h80200003 : 32'h0", "32'h80200005 : 32'h0")
    with pytest.raises(AssertionError):
        mf.wrapper_facts()


def test_a_truncated_sine_rom_is_detected(pack):
    rom = pack / "hdl" / "sig_rom.v"
    lines = rom.read_text().splitlines(keepends=True)
    keep = [ln for ln in lines if "mem[1023] =" not in ln]
    assert len(keep) == len(lines) - 1
    rom.write_text("".join(keep))
    with pytest.raises(AssertionError):
        mf.wrapper_facts()


# ---------------------------------------------------------------------------------------------------------------
# the timeline of sheet 31 (b)
# ---------------------------------------------------------------------------------------------------------------
def test_timeline_phases_add_up_to_the_prediction():
    T = mf.timeline()
    cfg, res = T["cfg"], T["res"]
    ph = dict(T["phases"])
    assert sum(ph.values()) == T["pred"] == M.expected_cycles(res, cfg)
    assert ph["capture"] == cfg.nsamp == 2560
    assert ph["replay: load"] == cfg.nchunks * 512 and ph["replay: FFT"] == cfg.nchunks * M.T_FFT_PER_CHUNK
    assert ph["segment merge"] == res.S * 512 * 19  # per segment 512 x (6 taps x 3 cycles + M_F)
    assert ph["window, accumulate"] == res.S * (5 * 3 * 1024 + 3 * 1024)  # 5 taps x 3 cycles, plus W_F, W_G, W_ACC
    assert ph["squarer"] == 3 * 1024 * res.S + res.squarer_cycles
    assert ph["handshakes"] == res.S + res.S * M.T_SEG_EXTRA + 2
    # the segment engine's data-independent part is the 31 233 cycles per segment of the CI measurement
    assert ph["segment merge"] + ph["window, accumulate"] + 3 * 1024 * res.S + res.S == 31233 * res.S
    assert (T["cfg"].nchunks, res.S) == (5, 4)


def test_timeline_record_is_the_generator_record():
    T = mf.timeline()
    want = siggen.generate(T["cfg"].nsamp, siggen.pinc_for(siggen.SELFTEST["bin_1024"]), siggen.SELFTEST["noise_shift"], 12)
    assert np.array_equal(T["x"], want)


def test_timeline_depends_on_the_data_only_through_the_squarer():
    quiet = mf.timeline(noise_shift=11)
    loud = mf.timeline(noise_shift=1)
    assert quiet["res"].squarer_cycles != loud["res"].squarer_cycles
    fixed = lambda T: T["pred"] - T["res"].squarer_cycles  # noqa: E731
    assert fixed(quiet) == fixed(loud)


# ---------------------------------------------------------------------------------------------------------------
# the sheets
# ---------------------------------------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def drawn(tmp_path_factory):
    out = tmp_path_factory.mktemp("figs")
    real = mf.FIGS
    mf.FIGS = out
    try:
        d30, d31 = mf.fig30_test_bench(), mf.fig31_wrapper_generator()
    finally:
        mf.FIGS = real
    return out, d30, d31


def test_sheets_are_written_as_png_and_svg(drawn):
    out, _, _ = drawn
    for stem in ("fig30_mc_test_bench", "fig31_mc_wrapper_generator"):
        assert (out / f"{stem}.png").stat().st_size > 100_000
        assert (out / f"{stem}.svg").stat().st_size > 50_000


def test_two_renderings_are_byte_identical(drawn, tmp_path):
    out, _, _ = drawn
    real = mf.FIGS
    mf.FIGS = tmp_path
    try:
        mf.fig30_test_bench()
        mf.fig31_wrapper_generator()
    finally:
        mf.FIGS = real
    for stem in ("fig30_mc_test_bench", "fig31_mc_wrapper_generator"):
        assert (tmp_path / f"{stem}.svg").read_bytes() == (out / f"{stem}.svg").read_bytes()


def test_fig30_names_the_blocks_ports_and_types_of_the_pack(drawn):
    _, d30, _ = drawn
    t = _texts(d30)
    for needle in ("Sine Wave + Random Number", "Signal Generator", "Chirp Signal", "From Workspace  bt_stim", "bt_addr_stim",
                   "Gateway In", "Gateway\nOut", "Fix_12_11", "UFix_10_0", "UFix_32_0", "UFix_1_0", "Boolean", "Vivado Simulator",
                   "Model Composer Hub", "bt_busy", "bt_done", "bt_rd_data", "bt_in, bt_valid, bt_addr", "export_mc_stimulus.py",
                   "bt_mc_build.m", "bt_mc_run.m", "compare_mc_result.py", "WL = 12  NCHUNKS = 5", "2560 samples (5 chunks)"):
        assert needle in t, needle
    for port in ("rst", "in_valid", "in_data", "rd_addr", "busy", "done", "rd_data", "clk", "ce"):
        assert port in t


def test_fig30_does_not_claim_a_simulink_execution(drawn):
    _, d30, _ = drawn
    t = _texts(d30)
    assert "was not executed where this repository was written" in t


def test_fig31_labels_the_timeline_as_a_prediction(drawn):
    _, _, d31 = drawn
    t = _texts(d31)
    T = mf.timeline()
    assert "PREDICTION" in t and "not a simulation result" in t and "predicted clock cycles" in t
    assert mf._n(T["pred"]) in t and "145 790" in t
    for name in ("segment merge", "window, accumulate", "squarer"):
        assert f"{name}: " in t
    assert "traces idealised from the prediction" in t
    assert mf._n(T["res"].squarer_cycles) in t


def test_fig31_names_the_states_conditions_and_generator_constants(drawn):
    _, _, d31 = drawn
    t = _texts(d31)
    for needle in ("S_CAPTURE", "S_PRIME", "S_REPLAY", "S_DRAIN", "S_DONE", "wr_cnt = NSAMP - 1", "accept & rp = NSAMP - 1", "!core_busy",
                   "core_rst_n = &rst_sr", "avg_shift = clog2(NCHUNKS - 1) = 2", "A5A51F3D", "x^32+x^22+x^2+x+1", "420 982 292",
                   "not the AMD DDS Compiler IP", "bin 100.37", "exact Welch average", "float six-tap merge", "model words (rd_data)"):
        assert needle in t, needle


def test_fig31_spectrum_panel_shows_the_tone_where_the_model_puts_it(drawn):
    _, _, d31 = drawn
    axes = [a for a in d31.fig.axes if a.get_xlabel() == "bin of 1 024"]
    assert len(axes) == 1
    a = axes[0]
    ys = {ln.get_label(): ln for ln in a.get_lines() if not ln.get_label().startswith("_")}
    assert set(ys) == {"exact Welch average", "float six-tap merge", "model words (rd_data)"}
    model = ys["model words (rd_data)"].get_ydata()
    assert len(model) == 512 and int(np.argmax(model)) == 100 and max(model) == pytest.approx(0.0)
