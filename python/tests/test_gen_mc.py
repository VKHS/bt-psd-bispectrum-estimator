"""Tests of the generator of the Model Composer files (python/scripts/gen_mc.py) and of the committed generated files."""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python" / "scripts"))
sys.path.insert(0, str(ROOT / "scripts"))

import gen_mc  # noqa: E402
import lint_rtl  # noqa: E402

from btpsd import mcmodel as M  # noqa: E402
from btpsd import rtlparse  # noqa: E402

MC = ROOT / "modelcomposer"
BUNDLE = (MC / "hdl" / "bt_psd_mc_bundle.v").read_text()


def _defined(text):
    return re.findall(r"(?m)^\s*module\s+(\w+)", text)


def test_committed_files_are_current(tmp_path):
    gen_mc.write_all(tmp_path, **{"wl": 12, "nchunks": 5, "rounding": 1, "sq_serial": 1})
    for rel in ("hdl/bt_psd_mc_bundle.v", "matlab/bt_psd_mc_config.m", "matlab/bt_mc_params.m", "hdl/sig_rom.v"):
        assert (tmp_path / rel).read_text() == (MC / rel).read_text(), f"{rel} is stale: run python/scripts/gen_mc.py"


def test_defaults_of_the_generator_are_the_defaults_of_the_files():
    d = gen_mc.DEFAULTS
    assert (d["wl"], d["nchunks"], d["rounding"], d["sq_serial"]) == (12, 5, 1, 1)
    params = rtlparse.default_parameters(MC / "hdl" / "bt_psd_mc.v")
    assert params["WL"] == 12 and params["NCHUNKS"] == 5 and params["ROUND"] == 1 and params["SQ_SERIAL"] == 1


def test_bundle_has_one_module_each_in_dependency_order():
    mods = _defined(BUNDLE)
    assert mods == ["delay_line", "ram_1r1w", "cru_rom", "cru_rotate", "r8_butterfly", "fft512_mem", "mag_sq_odd", "psd_engine",
                    "bt_psd_top", "bt_psd_mc"]
    # an instantiated module is always defined above its first use
    known = set(mods)
    seen = set()
    chunks = re.split(r"(?m)^\s*module\s+(\w+)", BUNDLE)
    for i in range(1, len(chunks), 2):
        name, body = chunks[i], chunks[i + 1]
        body = re.sub(r"//[^\n]*", "", body)
        for m in known - {name}:
            if re.search(rf"\b{m}\s*(#\s*\(|\w+\s*\()", body):
                assert m in seen, f"{name} instantiates {m} before it is defined"
        seen.add(name)


def test_bundle_has_no_directives_except_the_timescale_and_no_includes():
    uses = re.findall(r"`(\w+)", re.sub(r"//[^\n]*", "", BUNDLE))
    assert uses == ["timescale"], uses
    assert "include" not in re.sub(r"//[^\n]*", "", BUNDLE)


def test_bundle_is_lint_clean(capsys):
    assert lint_rtl.main([str(MC / "hdl" / "bt_psd_mc_bundle.v")]) == 0
    assert "0 problems" in capsys.readouterr().out


def test_sources_are_lint_clean_together(capsys):
    files = sorted((ROOT / "rtl").glob("*.v")) + sorted((MC / "hdl").glob("*.v")) + sorted((ROOT / "tb").glob("*.v"))
    assert lint_rtl.main([str(f) for f in files]) == 0
    assert "0 problems" in capsys.readouterr().out


def test_every_rtl_source_appears_in_the_bundle_with_macros_expanded():
    macros = gen_mc.header_macros()
    for f in gen_mc.BUNDLE_ORDER:
        flat = gen_mc.flatten((ROOT / "rtl" / f).read_text(), macros, f)
        flat = re.sub(r"(?m)^`timescale[^\n]*\n", "", flat)
        assert flat in BUNDLE, f"{f} is not contained verbatim in the bundle"
    assert "(-256)" in BUNDLE and "parameter signed [15:0] WM2 = 0," in BUNDLE


def test_flatten_rejects_unknown_includes_and_macros():
    macros = gen_mc.header_macros()
    with pytest.raises(ValueError, match="include"):
        gen_mc.flatten('`include "other.vh"\nmodule m; endmodule\n', macros, "x.v")
    with pytest.raises(ValueError, match="unknown macro"):
        gen_mc.flatten("module m; wire a = `NOSUCH; endmodule\n", macros, "x.v")


@pytest.mark.parametrize("wl,nch,rnd,sq", [(12, 5, 1, 1), (16, 9, 0, 0), (10, 3, 1, 0), (14, 17, 1, 1)])
def test_configurations_change_exactly_the_wrapper_defaults(wl, nch, rnd, sq):
    text = gen_mc.bundle_text(wl, nch, rnd, sq)
    wrapper = text[text.index("module bt_psd_mc"):]
    assert re.search(rf"parameter WL\s*=\s*{wl},", wrapper)
    assert re.search(rf"parameter ROUND\s*=\s*{rnd},", wrapper)
    assert re.search(rf"parameter SQ_SERIAL\s*=\s*{sq},", wrapper)
    assert re.search(rf"parameter NCHUNKS\s*=\s*{nch},", wrapper)
    assert text[: text.index("module bt_psd_mc")].count("module ") == BUNDLE[: BUNDLE.index("module bt_psd_mc")].count("module ")
    # the part before the wrapper is independent of the configuration (the core keeps its own defaults)
    assert text[: text.index("module bt_psd_mc")].split("// ===== modelcomposer")[0].split("\n", 3)[3] == \
        BUNDLE[: BUNDLE.index("module bt_psd_mc")].split("// ===== modelcomposer")[0].split("\n", 3)[3]


def test_set_defaults_requires_every_parameter():
    src = gen_mc.WRAPPER.read_text()
    with pytest.raises(ValueError, match="NCHUNKS"):
        gen_mc.set_defaults(src.replace("NCHUNKS   = 5", "NCHUNKX   = 5"), 12, 5, 1, 1)


@pytest.mark.parametrize("wl,nch", [(12, 5), (16, 9), (10, 3)])
def test_config_function_matches_the_ports_of_the_verilog(wl, nch):
    text = gen_mc.config_m(wl, nch, 1, 1)
    pl = rtlparse.ports(MC / "hdl" / "bt_psd_mc.v", dict(WL=wl, NCHUNKS=nch, ROUND=1, SQ_SERIAL=1))
    ins = [p for p in pl if p.direction == "input" and p.name not in ("clk", "ce")]
    outs = [p for p in pl if p.direction == "output"]
    assert re.findall(r"addSimulinkInport\('(\w+)'\)", text) == [p.name for p in ins] == ["rst", "in_valid", "in_data", "rd_addr"]
    assert re.findall(r"addSimulinkOutport\('(\w+)'\)", text) == [p.name for p in outs] == ["busy", "done", "rd_data"]
    for p in outs:
        assert f"{p.name}_port.setType('UFix_{p.width}_0');" in text
    for p in ins:
        assert f"this_block.port('{p.name}').width ~= {p.width}" in text
    assert [(p.name, p.width) for p in pl] == [("clk", 1), ("ce", 1), ("rst", 1), ("in_valid", 1), ("in_data", wl), ("rd_addr", 10),
                                              ("busy", 1), ("done", 1), ("rd_data", 2 * wl + 8)]
    assert "setEntityName('bt_psd_mc')" in text and "setTopLevelLanguage('Verilog')" in text
    assert "setup_as_single_rate(this_block,'clk','ce')" in text and "addClkCEPair(clkname,cename,theInputRate)" in text
    assert "bt_psd_mc_bundle.v" in text and "tagAsCombinational" not in text.replace("do not tag the block as combinational", "")
    assert text.count("function ") == 2 and text.count("\nreturn;") + text.count("return;") >= 2


def test_black_box_ports_are_all_single_rate_registered():
    """The wrapper must not have a combinational path from an input to an output (the config does not tag it as combinational)."""
    src = re.sub(r"//[^\n]*", "", (MC / "hdl" / "bt_psd_mc.v").read_text())
    assigns = re.findall(r"assign\s+(\w+)\s*=([^;]+);", src)
    outputs = {"busy", "done"}
    for name, rhs in assigns:
        if name in outputs:
            assert not re.search(r"\b(rst|in_valid|in_data|rd_addr|ce)\b", rhs), f"{name} depends combinationally on an input"


@pytest.mark.parametrize("wl,nch,rnd", [(12, 5, 1), (16, 9, 0), (10, 3, 1)])
def test_params_function_matches_the_python_model(wl, nch, rnd):
    text = gen_mc.params_m(gen_mc.build_params(wl, nch, rnd, 1))
    cfg = M.McConfig(wl=wl, nchunks=nch, round=rnd)
    vals = {k: float(v) for k, v in re.findall(r"(?m)^p\.(\w+) = ([0-9.e+-]+);", text)}
    assert vals["wl"] == wl and vals["nchunks"] == nch and vals["nsamp"] == cfg.nsamp and vals["segments"] == cfg.segments
    assert vals["avg_shift"] == cfg.avg_shift and vals["accw"] == cfg.accw and vals["bins"] == 1024
    assert vals["scale"] == cfg.scale
    assert f"p.in_type = 'Fix_{wl}_{wl - 1}';" in text and f"p.out_type = 'UFix_{cfg.accw}_0';" in text


def test_invalid_requests_are_rejected(tmp_path):
    with pytest.raises(ValueError):
        gen_mc.write_all(tmp_path, 12, 1, 1, 1)
    with pytest.raises(ValueError):
        gen_mc.write_all(tmp_path, 4, 5, 1, 1)


def test_cli_writes_into_the_given_directory(tmp_path, capsys, monkeypatch):
    monkeypatch.setattr(sys, "argv", ["gen_mc.py", "--wl", "14", "--nchunks", "9", "--out-dir", str(tmp_path)])
    gen_mc.main()
    out = capsys.readouterr().out
    assert out.count("wrote") == 4
    assert "WL = 14" in (tmp_path / "matlab" / "bt_psd_mc_config.m").read_text()
    assert "Fix_14_13" in (tmp_path / "matlab" / "bt_mc_params.m").read_text()
