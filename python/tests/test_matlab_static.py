"""Static checks of the MATLAB files of the Model Composer pack (no MATLAB or Octave is available where this is developed).

The checker is deliberately simple: it strips comments and strings, then verifies that quotes close, brackets balance, every
block keyword (function, if, for, while, switch, try) has its ``end``, that no ``'`` is used as a transpose operator (the
files avoid transposes so that this check is exact), and that the names the scripts share with each other and with the Python
tools agree.  It finds typos; it does not prove that the scripts run.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
MATLAB = ROOT / "modelcomposer" / "matlab"
sys.path.insert(0, str(ROOT / "python" / "scripts"))

import compare_mc_result  # noqa: E402

OPENERS = {"function", "if", "for", "parfor", "while", "switch", "try"}


def scan(text: str, name: str = "file.m"):
    """Return (list of code lines without comments, strings emptied, plus a list of problems)."""
    problems, code = [], []
    stack = []
    cont = False
    for ln, raw in enumerate(text.splitlines(), 1):
        out, i, in_str = [], 0, False
        while i < len(raw):
            c = raw[i]
            if in_str:
                if c == "'":
                    if i + 1 < len(raw) and raw[i + 1] == "'":
                        i += 2
                        continue
                    in_str = False
                    out.append("'")
                i += 1
                continue
            if c == "%":
                break
            if raw.startswith("...", i):
                cont = True
                break
            if c == "'":
                prev = out[-1] if out else " "
                if prev.isalnum() or prev in ")]}_.'":
                    problems.append(f"{name}:{ln}: quote after '{prev}' (transpose?): not allowed in these files")
                in_str = True
                out.append("'")
                i += 1
                continue
            out.append(c)
            i += 1
        if in_str:
            problems.append(f"{name}:{ln}: unterminated string")
        line = "".join(out)
        for ch in line:
            if ch in "([{":
                stack.append((ch, ln))
            elif ch in ")]}":
                if not stack or "([{".index(stack[-1][0]) != ")]}".index(ch):
                    problems.append(f"{name}:{ln}: unbalanced {ch}")
                    stack = []
                else:
                    stack.pop()
        if stack and stack[-1][0] == "(" and not cont and line.strip():
            problems.append(f"{name}:{ln}: line ends inside ( without a continuation")
        cont = cont and raw.rstrip().endswith("...")
        code.append(line)
    if stack:
        problems.append(f"{name}: unclosed {stack[-1][0]} opened in line {stack[-1][1]}")
    return code, problems


def blocks(code, name="file.m", functions_end=True):
    """Check that every block keyword has its ``end``; ``functions_end=False`` for files whose functions end with ``return``."""
    problems, depth = [], 0
    for ln, line in enumerate(code, 1):
        # statements at bracket depth 0
        parts, cur, d = [], "", 0
        for ch in line:
            if ch in "([{":
                d += 1
            elif ch in ")]}":
                d -= 1
            if ch in ";," and d == 0:
                parts.append(cur)
                cur = ""
            else:
                cur += ch
        parts.append(cur)
        for st in parts:
            w = re.match(r"\s*([A-Za-z_]\w*)", st)
            if not w:
                continue
            word = w.group(1)
            if word in OPENERS and not re.match(r"\s*\w+\s*=", st) and (functions_end or word != "function"):
                depth += 1
            elif word == "end" and re.fullmatch(r"\s*end\s*", st):
                depth -= 1
                if depth < 0:
                    problems.append(f"{name}:{ln}: end without opener")
                    depth = 0
    if depth != 0:
        problems.append(f"{name}: {depth} block(s) not closed")
    return problems


FILES = sorted(MATLAB.glob("*.m"))


def test_the_expected_files_exist():
    assert [f.name for f in FILES] == ["bt_mc_build.m", "bt_mc_params.m", "bt_mc_run.m", "bt_psd_mc_config.m"]


@pytest.mark.parametrize("path", FILES, ids=lambda p: p.name)
def test_syntax_level_checks(path):
    code, problems = scan(path.read_text(), path.name)
    # the generated configuration is in the style of the Black Box wizard: functions end with return, not with end
    problems += blocks(code, path.name, functions_end=path.name != "bt_psd_mc_config.m")
    assert problems == []


def test_the_checker_itself_finds_typos():
    code, p = scan("x = 'abc;\n")
    assert any("unterminated" in q for q in p)
    code, p = scan("y = f(1, [2 3);\n")
    assert p
    code, p = scan("z = a';\n")
    assert any("transpose" in q for q in p)
    code, p = scan("if a\n  b = 1;\n")
    assert blocks(code)
    code, p = scan("function f(x)\n  y = x;\nreturn;\n")
    assert blocks(code) and blocks(code, functions_end=False) == []
    code, p = scan("if a, b = 1; end\nfor k = 1:3\n  s = s + k;\nend\n")
    assert p == [] and blocks(code) == []
    code, p = scan("x = 'it''s fine'; % comment with ' quote\nfprintf('%d\\n', x);\n")
    assert p == []
    code, p = scan("v = f(a, ...\n      b);\n")
    assert p == []


def test_script_functions_are_terminated_consistently():
    for name in ("bt_mc_build.m", "bt_mc_run.m"):
        text = (MATLAB / name).read_text()
        assert re.match(r"function ", text)
        assert text.rstrip().endswith("end")
    cfg = (MATLAB / "bt_psd_mc_config.m").read_text()
    assert cfg.count("\nreturn;") >= 2  # the generated configuration uses return, as the wizard output does


# ---------------------------------------------------------------------------------------------------------------------
# consistency between the scripts, the generated configuration and the Python checker
# ---------------------------------------------------------------------------------------------------------------------
def _config_ports():
    cfg = (MATLAB / "bt_psd_mc_config.m").read_text()
    return re.findall(r"addSimulinkInport\('(\w+)'\)", cfg), re.findall(r"addSimulinkOutport\('(\w+)'\)", cfg)


def test_black_box_port_numbers_used_in_the_build_script_follow_the_config():
    ins, outs = _config_ports()
    build = (MATLAB / "bt_mc_build.m").read_text()
    gate = {"rst": "gw_in_rst", "in_valid": "gw_in_valid", "in_data": "gw_in_data", "rd_addr": "gw_in_addr"}
    wired = dict((a, int(n)) for a, n in re.findall(r"wire_up\('(gw_in_\w+)/1', 'bt_psd_mc/(\d)'\)", build))
    assert sorted(wired) == sorted(gate.values())
    for port, g in gate.items():
        assert wired[g] == ins.index(port) + 1, f"{port}: Black Box input {ins.index(port) + 1}"
    taps = re.findall(r"'bt_psd_mc/(\d)', '(\w+)', '(\w+)'", build)
    assert [(int(n), name) for n, name, _ in taps] == [(outs.index(n) + 1, n) for n in outs]


def test_logged_variables_match_between_build_and_run_scripts_and_python():
    build = (MATLAB / "bt_mc_build.m").read_text()
    run = (MATLAB / "bt_mc_run.m").read_text()
    made = set(re.findall(r"'(bt_\w+)'\s*;?\s*(?:\n|\})", "\n".join(re.findall(r"taps = \{.*?\};", build, re.S))))
    made |= set(re.findall(r"'(bt_\w+)'", " ".join(re.findall(r"taps = \{.*?\};", build, re.S))))
    used = set(re.findall(r"from_base\('(bt_\w+)'\)", run)) | {"bt_done"}
    assert {"bt_busy", "bt_done", "bt_rd_data", "bt_in", "bt_valid", "bt_addr"} <= made
    assert used <= made
    # the fields written to bt_log are the signals that compare_mc_result reads
    fields = set(re.findall(r"res\.(\w+) =", run))
    assert set(compare_mc_result.SIGNALS) <= fields and {"wl", "nchunks"} <= fields
    # the address stimulus variable is the same in both scripts
    assert "'bt_addr_stim'" in build and "'bt_addr_stim'" in run


def test_parameters_used_by_the_scripts_exist_in_the_generated_params():
    params = (MATLAB / "bt_mc_params.m").read_text()
    have = set(re.findall(r"(?m)^p\.(\w+) =", params))
    for name in ("bt_mc_build.m", "bt_mc_run.m"):
        used = set(re.findall(r"\bp\.(\w+)", (MATLAB / name).read_text()))
        assert used <= have, f"{name} uses p.{sorted(used - have)}"


def test_the_stimulus_file_variables_are_the_ones_the_exporter_writes():
    exp = (ROOT / "python" / "scripts" / "export_mc_stimulus.py").read_text()
    build = (MATLAB / "bt_mc_build.m").read_text()
    run = (MATLAB / "bt_mc_run.m").read_text()
    assert '"bt_stim"' in exp and "'VariableName', 'bt_stim'" in build
    assert '"bt_stop_cycles"' in exp and "bt_stop_cycles" in run
    assert "modelcomposer/stimulus" in build and (ROOT / "python" / "scripts" / "export_mc_stimulus.py").exists()


def test_matlab_files_say_that_they_were_not_executed():
    for name in ("bt_mc_build.m", "bt_mc_run.m", "bt_psd_mc_config.m"):
        assert re.search(r"NOT EXECUTED|Not executed", (MATLAB / name).read_text()), name


def test_the_build_script_switches_the_black_box_to_the_vivado_simulator():
    build = (MATLAB / "bt_mc_build.m").read_text()
    # a Black Box starts in the mode 'Inactive' (it outputs zeros): the script must change it, or list it as a manual step
    assert re.search(r"set_sim_mode\('bt_psd_mc', 'Vivado Simulator'\)", build)
    body = build.split("function set_sim_mode")[1].split("function wire_up")[0]
    assert "'DialogParameters'" in body and "simulation\\s*mode" in body
    assert "notes{end + 1}" in body and "Inactive" in body
