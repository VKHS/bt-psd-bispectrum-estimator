"""Tests of scripts/lint_rtl.py: the default file set and the consistency check of modules defined in more than one file."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import lint_rtl  # noqa: E402

LEAF = "module leaf #(parameter W = 4) (input clk, input [W-1:0] d, output reg [W-1:0] q);\n always @(posedge clk) q <= d;\nendmodule\n"
TOP = "module top (input clk, input [3:0] a, output [3:0] y);\n leaf #(.W(4)) u (.clk(clk), .d(a), .q(y));\nendmodule\n"


def _write(tmp_path, **files):
    out = []
    for name, text in files.items():
        p = tmp_path / f"{name}.v"
        p.write_text(text)
        out.append(str(p))
    return out


def test_repository_sources_are_clean_and_the_default_set_includes_the_model_composer_hdl(capsys):
    assert lint_rtl.main([]) == 0
    out = capsys.readouterr().out
    n_files = len(list((ROOT / "rtl").glob("*.v")) + list((ROOT / "tb").glob("*.v")) + list((ROOT / "modelcomposer" / "hdl").glob("*.v")))
    assert f"lint_rtl: {n_files} files" in out and "0 problems" in out
    assert any((ROOT / "modelcomposer" / "hdl").glob("*.v")) and n_files >= 25


def test_a_module_repeated_with_the_same_interface_is_accepted(tmp_path, capsys):
    assert lint_rtl.main(_write(tmp_path, a=LEAF, b=LEAF + TOP)) == 0
    assert "0 problems" in capsys.readouterr().out


def test_a_repeated_module_with_another_port_list_is_reported(tmp_path, capsys):
    other = LEAF.replace("output reg [W-1:0] q", "output reg [W-1:0] q, output vld")
    assert lint_rtl.main(_write(tmp_path, a=LEAF, b=other + TOP)) == 1
    assert "module leaf is defined in a.v and in b.v with different parameters or ports" in capsys.readouterr().out


def test_a_repeated_module_with_another_parameter_list_is_reported(tmp_path, capsys):
    other = LEAF.replace("parameter W = 4", "parameter W = 4, parameter K = 1")
    assert lint_rtl.main(_write(tmp_path, a=LEAF, b=other + TOP)) == 1
    assert "different parameters or ports" in capsys.readouterr().out


def test_the_older_checks_still_work(tmp_path, capsys):
    bad = TOP.replace(".d(a)", ".dd(a)")
    assert lint_rtl.main(_write(tmp_path, a=LEAF + bad)) == 1
    assert "unknown port .dd" in capsys.readouterr().out
    unclosed = LEAF.replace("endmodule", "")
    assert lint_rtl.main(_write(tmp_path, a=unclosed)) == 1
