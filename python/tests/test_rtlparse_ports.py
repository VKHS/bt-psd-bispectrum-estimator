"""The port and parameter reader of ``btpsd.rtlparse`` against the headers of the Verilog sources.

The expected tables are written out by hand from the module headers; the reader is used to generate the Model Composer
Black Box configuration and to label the block diagrams, so a mistake here would propagate into both.
"""
from pathlib import Path

import pytest

from btpsd.rtlparse import RTL_DIR, Port, default_parameters, parse_all, ports


def _p(module, env=None):
    return {p.name: p for p in ports(RTL_DIR / f"{module}.v", env)}


def test_fft_core_ports():
    got = ports(RTL_DIR / "fft512_mem.v")
    assert [p.name for p in got] == ["clk", "rst_n", "ld_en", "ld_addr", "ld_re", "ld_im", "start", "busy", "done",
                                     "rd_k", "rd_re", "rd_im"]
    d = {p.name: p for p in got}
    assert (d["ld_addr"].direction, d["ld_addr"].width) == ("input", 9)
    assert (d["ld_re"].width, d["ld_re"].signed) == (16, True)
    assert (d["busy"].direction, d["busy"].is_reg, d["busy"].width) == ("output", True, 1)
    assert (d["rd_im"].direction, d["rd_im"].width, d["rd_im"].signed) == ("output", 16, True)
    assert not d["clk"].signed and d["clk"].width == 1


def test_parameter_override_changes_widths():
    assert _p("fft512_mem", {"WL": 12})["ld_re"].width == 12
    assert _p("r8_butterfly", {"WL": 12})["xr"].width == 8 * 12
    assert _p("r8_butterfly")["tq"].width == 70
    assert _p("cru_rotate", {"WL": 20})["yr"].width == 20


def test_dependent_parameters_follow_their_inputs():
    assert default_parameters(RTL_DIR / "bt_psd_top.v") == {"WL": 12, "ROUND": 1, "SQ_SERIAL": 1, "ACCW": 32}
    d = _p("bt_psd_top", {"WL": 10})
    assert d["in_data"].width == 10
    assert d["res_data"].width == 2 * 10 + 8  # ACCW = 2 WL + 8 follows WL
    assert _p("bt_psd_top", {"WL": 10, "ACCW": 40})["res_data"].width == 40


def test_macro_defaults_are_skipped_unless_supplied():
    assert "WM2" not in default_parameters(RTL_DIR / "bt_psd_top.v")  # default is a `define macro
    assert default_parameters(RTL_DIR / "bt_psd_top.v", {"WM2": 5})["WM2"] == 5


def test_top_level_ports():
    d = _p("bt_psd_top")
    assert list(d) == ["clk", "rst_n", "avg_shift", "in_valid", "in_data", "in_ready", "busy", "res_addr", "res_data"]
    assert (d["avg_shift"].width, d["res_addr"].width) == (5, 10)
    assert d["in_data"].signed and not d["res_data"].signed
    b = _p("bt_bispectrum_top")
    assert (b["res_re"].width, b["res_im"].width, b["res_re"].signed) == (34, 34, True)  # BACCW = 2 WL + 2 LOGB


def test_ram_and_delay_line_ports():
    r = _p("ram_1r1w")
    assert (r["wdata"].width, r["rdata"].width, r["waddr"].width, r["rdata"].is_reg) == (32, 32, 6, True)
    assert _p("ram_1r1w", {"DW": 40, "AW": 8})["raddr"].width == 8
    assert _p("delay_line", {"W": 70})["q"].width == 70


def test_every_module_header_is_readable_and_clocked():
    for name in parse_all():
        got = ports(RTL_DIR / f"{name}.v")
        assert got and all(isinstance(p, Port) and p.width >= 1 for p in got), name
        if name != "cru_rom":  # the generated ROM is combinational
            assert any(p.name == "clk" for p in got), name


def test_unparsable_header_is_reported(tmp_path):
    f = tmp_path / "bad.v"
    f.write_text("module bad (input a, output wire [3:0] b, c d);\nendmodule\n")
    with pytest.raises(ValueError):
        ports(f)


def test_shorthand_port_lists_inherit_the_declaration(tmp_path):
    f = tmp_path / "two.v"
    f.write_text("module two #(parameter N = 4) (input clk, input [N-1:0] a, b, output y);\nendmodule\n")
    got = {p.name: p for p in ports(f)}
    assert got["a"].width == got["b"].width == 4 and got["b"].direction == "input"
    assert got["y"].direction == "output"


def test_sized_literals_in_parameter_defaults_are_evaluated(tmp_path):
    from btpsd.rtlparse import _eval, default_parameters

    assert _eval("32'hA5A51F3D", {}) == 0xA5A51F3D
    assert _eval("4'b1010 + 1", {}) == 11 and _eval("12'o777", {}) == 511 and _eval("8'd200", {}) == 200
    assert _eval("8'sd200", {}) == -56 and _eval("16'sh7FFF", {}) == 32767 and _eval("16'shFFFF", {}) == -1
    assert _eval("8'hF_F", {}) == 255
    f = tmp_path / "lit.v"
    f.write_text("module lit #(parameter W = 12, parameter [31:0] SEED = 32'hA5A51F3D, parameter signed [7:0] K = -8'sd3) "
                 "(input clk, input [W-1:0] a);\nendmodule\n")
    v = default_parameters(f)
    assert v == {"W": 12, "SEED": 0xA5A51F3D, "K": -3}
    assert [p.width for p in ports(f)] == [1, 12]


def test_the_model_composer_modules_are_parsed():
    mc = Path(__file__).resolve().parents[2] / "modelcomposer" / "hdl"
    g = {p.name: p for p in ports(mc / "sig_gen.v")}
    assert (g["pinc"].width, g["noise_shift"].width, g["x"].width, g["valid"].width) == (32, 4, 12, 1) and g["x"].signed
    s = {p.name: p.width for p in ports(mc / "bt_psd_selftest.v")}
    assert s["rd_data"] == 32 and s["pinc"] == 32 and s["rd_addr"] == 10
    w = {p.name: p.width for p in ports(mc / "bt_psd_mc.v")}
    assert (w["in_data"], w["rd_addr"], w["rd_data"]) == (12, 10, 32)
