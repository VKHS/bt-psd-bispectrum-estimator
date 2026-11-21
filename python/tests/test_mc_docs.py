"""``docs/MODEL_COMPOSER.md`` and ``modelcomposer/README.md`` against the sources, the models and the stored simulation results.

The guide quotes paths, command lines, ports, parameters and numbers.  Each of them is checked here against the file, the argparse
help, the Verilog header, the MATLAB source, the model or the CSV it comes from, so that the text cannot drift away from the code.
"""
from __future__ import annotations

import csv
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python" / "scripts"))

import export_mc_stimulus as ex  # noqa: E402
import make_mc_figures as mf  # noqa: E402
from btpsd import mcmodel as M  # noqa: E402
from btpsd import paperdata, rtlparse, signals, siggen  # noqa: E402
from btpsd.estimator import psd_exact_welch, psd_float_fd  # noqa: E402

GUIDE = ROOT / "docs" / "MODEL_COMPOSER.md"
PACK_README = ROOT / "modelcomposer" / "README.md"
MC_V = ROOT / "modelcomposer" / "hdl" / "bt_psd_mc.v"
SIG_GEN_V = ROOT / "modelcomposer" / "hdl" / "sig_gen.v"
MATLAB = ROOT / "modelcomposer" / "matlab"
GENERATED = ("modelcomposer/stimulus/", "modelcomposer/work/", "sim/work/")
SCRIPTS = ("gen_mc.py", "export_mc_stimulus.py", "compare_mc_result.py")


# ---------------------------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def text() -> str:
    return GUIDE.read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def flat(text) -> str:
    """The guide with every run of white space (line breaks of wrapped paragraphs included) as one space."""
    return re.sub(r"\s+", " ", text)


@pytest.fixture(scope="module")
def timelines() -> dict:
    """``mf.timeline`` (model prediction on the generator record) for the chunk counts the guide tabulates."""
    return {c: mf.timeline(12, c) for c in (2, 3, 4, 5, 9)}


def table(text: str, header: str) -> list[list[str]]:
    """Rows (backticks and bold marks removed) of the markdown table whose header line starts with ``header``."""
    lines = text.splitlines()
    start = next(i for i, ln in enumerate(lines) if ln.startswith(header))
    rows = []
    for ln in lines[start + 2:]:
        if not ln.startswith("|"):
            break
        rows.append([c.strip().replace("`", "").replace("**", "") for c in ln.strip().strip("|").split("|")])
    return rows


def num(cell: str) -> int:
    """An integer written with spaces between the thousands."""
    return int(cell.replace(" ", ""))


def code_spans(text: str) -> list[str]:
    """Command lines of the fenced bash blocks (continuations joined, comments dropped) and the single-line inline code spans."""
    out = []
    for block in re.findall(r"```bash\n(.*?)```", text, re.S):
        joined = re.sub(r"\\\n\s*", " ", block)
        for ln in joined.splitlines():
            ln = re.sub(r"\s+#.*$", "", ln.strip())
            if ln and not ln.startswith("#"):
                out.append(ln)
    out += re.findall(r"`([^`\n]+)`", re.sub(r"```.*?```", "", text, flags=re.S))
    return out


_help_cache: dict[str, str] = {}


def help_of(script: str) -> str:
    if script not in _help_cache:
        env = dict(os.environ, PYTHONPATH=str(ROOT / "python"))
        r = subprocess.run([sys.executable, str(ROOT / "python" / "scripts" / script), "--help"], capture_output=True, text=True,
                           env=env, cwd=ROOT, timeout=120)
        assert r.returncode == 0, (script, r.stderr)
        _help_cache[script] = r.stdout
    return _help_cache[script]


def tree_file_names() -> set[str]:
    names = set()
    for _, dirs, files in os.walk(ROOT):
        dirs[:] = [d for d in dirs if d not in {".git", "__pycache__", "work", ".pytest_cache"}]
        names.update(files)
    return names


def black_box_ports() -> tuple[list[str], list[str]]:
    cfg = (MATLAB / "bt_psd_mc_config.m").read_text()
    return re.findall(r"addSimulinkInport\('(\w+)'\)", cfg), re.findall(r"addSimulinkOutport\('(\w+)'\)", cfg)


# ---------------------------------------------------------------------------------------------------------------
# structure, paths, commands
# ---------------------------------------------------------------------------------------------------------------
def test_sections_are_numbered_in_order_and_the_references_to_them_resolve(text, flat):
    heads = dict(re.findall(r"(?m)^## (\d+)\. (.+)$", text))
    assert list(map(int, heads)) == list(range(1, 12))
    assert "block by block" in heads["4"] and "Reproduce" in heads["9"]
    assert re.findall(r"(?m)^### (6\.\d) ", text) == ["6.1", "6.2", "6.3"]
    for ref in re.findall(r"[Ss]ections? (\d+)", flat):
        assert 1 <= int(ref) <= 11, ref
    # the MATLAB build script points the reader to section 4 of the guide
    assert "docs/MODEL_COMPOSER.md section 4 gives the same bench block by block" in (MATLAB / "bt_mc_build.m").read_text()
    assert "block-by-block bench" not in text  # section 4 is the bench; nothing else claims to be


@pytest.mark.parametrize("doc", [GUIDE, PACK_README], ids=["guide", "pack-readme"])
def test_every_file_the_text_names_exists(doc):
    """Backticked paths exist; backticked bare file names exist somewhere in the tree (generated outputs are exempt)."""
    names = tree_file_names()
    pat = re.compile(r"[\w./<>*-]+\.(?:py|v|vh|m|sh|md|png|svg|csv|json|mem|yml)")
    missing = []
    for tok in re.findall(r"`([^`\n]+)`", re.sub(r"```.*?```", "", doc.read_text(encoding="utf-8"), flags=re.S)):
        if not pat.fullmatch(tok) or any(c in tok for c in "<>*") or tok.startswith(GENERATED):
            continue
        ok = ((ROOT / tok).exists() or (doc.parent / tok).exists()) if "/" in tok else tok in names
        if not ok:
            missing.append(tok)
    assert not missing, missing


def test_the_sheets_are_linked_and_exist(text):
    for fig in ("fig30_mc_test_bench", "fig31_mc_wrapper_generator"):
        assert f"](figures/{fig}.png)" in text
        assert (ROOT / "docs" / "figures" / f"{fig}.png").is_file() and (ROOT / "docs" / "figures" / f"{fig}.svg").is_file()
    for ref in re.findall(r"`(fig\d\d)`", text):
        assert list((ROOT / "docs" / "figures").glob(ref + "_*.png")), ref


def test_every_command_uses_options_that_its_script_accepts(text):
    flags_of: dict[str, set[str]] = {}
    for ln in code_spans(text):
        m = re.search(r"python3 (python/scripts/\w+\.py)(?: (.*))?$", ln)
        if not m:
            continue
        assert (ROOT / m.group(1)).is_file(), ln
        flags_of.setdefault(Path(m.group(1)).name, set()).update(re.findall(r"(?<![\w-])(--[\w-]+)", m.group(2) or ""))
    assert {"gen_mc.py", "export_mc_stimulus.py", "compare_mc_result.py", "smoke_test.py", "make_mc_figures.py"} <= set(flags_of)
    for script, flags in flags_of.items():
        if flags:
            assert script in SCRIPTS, f"{script} takes options in the guide but is not checked here"
            for f in flags:
                assert f in help_of(script), (script, f)


def test_option_names_quoted_in_the_prose_belong_to_one_of_the_three_tools(text):
    quoted = set(re.findall(r"`(--[\w-]+)`", re.sub(r"```.*?```", "", text, flags=re.S)))
    assert {"--wl", "--nchunks", "--round", "--sq-serial", "--lag"} <= quoted
    for f in quoted:
        assert any(f in help_of(s) for s in SCRIPTS), f


def test_make_and_xsim_targets_exist(text):
    phony = set(re.search(r"(?m)^\.PHONY:(.*)$", (ROOT / "sim" / "Makefile").read_text()).group(1).split())
    xsim = set(re.findall(r"\[(\w+)\]=tb_", (ROOT / "sim" / "xsim_smoke.sh").read_text()))
    saw_make = saw_xsim = 0
    for ln in code_spans(text):
        if ln.startswith("make -C sim "):
            saw_make += 1
            assert set(ln.split()[3:]) <= phony, ln
        m = re.match(r"bash sim/xsim_smoke\.sh ([\w ]+)$", ln)
        if m:
            saw_xsim += 1
            assert set(m.group(1).split()) <= xsim, ln
    assert saw_make >= 3 and saw_xsim >= 2
    assert {"mc", "mcb", "sg", "st"} <= phony and {"mc", "mcb", "sg", "st"} <= xsim


# ---------------------------------------------------------------------------------------------------------------
# the Black Box interface
# ---------------------------------------------------------------------------------------------------------------
def test_port_table_equals_the_verilog_header_and_the_configuration_function(text):
    rows = table(text, "| Port | Direction")
    ports = [p for p in rtlparse.ports(MC_V) if p.name not in ("clk", "ce")]
    assert [r[0] for r in rows] == [p.name for p in ports], "the table lists the ports in the order of the module header"
    cfg = (MATLAB / "bt_psd_mc_config.m").read_text()
    for (name, direction, typ, _), p in zip(rows, ports):
        assert direction == ("in" if p.direction == "input" else "out"), name
        if p.direction == "input":
            want = "Boolean" if p.width == 1 else (f"Fix_{p.width}_{p.width - 1}" if p.signed else f"UFix_{p.width}_0")
            assert re.search(rf"port\('{name}'\)\.width ~= {p.width}\)", cfg), name  # the configuration function checks the width
        else:
            want = f"UFix_{p.width}_0"
            assert f"{name}_port.setType('{want}')" in cfg, name
        assert typ == want, (name, typ, want)


def test_black_box_port_numbers_in_the_text_follow_the_configuration_order(text, flat):
    ins, outs = black_box_ports()
    assert ins == ["rst", "in_valid", "in_data", "rd_addr"] and outs == ["busy", "done", "rd_data"]
    assert "(`rst`, `in_valid`, `in_data`, `rd_addr`)" in flat and "(`busy`, `done`, `rd_data`)" in flat
    rows = table(text, "| Block | Parameters | Wired to |")
    build = (MATLAB / "bt_mc_build.m").read_text()
    wires = re.findall(r"wire_up\('(\w+)/1', 'bt_psd_mc/(\d)'\)", build)
    assert sorted(int(k) for _, k in wires) == [1, 2, 3, 4]
    for gateway, k in wires:
        row = next(r for r in rows if r[0].endswith(gateway))
        assert row[2] == f"Black Box input {k} ({ins[int(k) - 1]})", (gateway, row)
    bb = next(r for r in rows if r[0].startswith("Black Box"))
    assert bb[2] == "outputs " + ", ".join(f"{i} {n}" for i, n in enumerate(outs, 1))


def test_default_build_facts_equal_the_params_file_and_the_wrapper_parameters(text, flat):
    params = dict(re.findall(r"p\.(\w+) = ([^;]+);", (MATLAB / "bt_mc_params.m").read_text()))
    assert (params["wl"], params["nchunks"], params["round"], params["sq_serial"]) == ("12", "5", "1", "1")
    assert params["in_type"] == "'Fix_12_11'" and params["out_type"] == "'UFix_32_0'"
    assert float(params["scale"]) == 2.0**-12 and "2^-12" in text
    assert "word length 12, five chunks, rounding to nearest, serial squarer" in flat
    header = re.sub(r"//[^\n]*", "", MC_V.read_text())
    for name in ("WL", "ROUND", "SQ_SERIAL", "NCHUNKS", "ACCW", "WM2", "WM1", "W0", "WP1", "WP2"):
        assert re.search(rf"parameter[^,\n]*\b{name}\b", header), name
    assert "`WM2` to `WP2`" in flat and "ACCW      = 2 * WL + 8" in header and "`2 WL + 8`" in text
    cfg = M.McConfig(wl=12, nchunks=5)
    assert (cfg.nsamp, cfg.segments, cfg.avg_shift, cfg.accw, cfg.scale) == (2560, 4, 2, 32, 2.0**-12)


def test_the_state_diagram_and_the_reset_time_of_the_text_match_the_sources(text, flat):
    line = "rst --> CAPTURE --(NSAMP samples seen)--> PRIME --> REPLAY --(last sample accepted)--> DRAIN --(core idle)--> DONE"
    verilog = MC_V.read_text()
    assert line in text and line in verilog
    for state in ("S_CAPTURE", "S_PRIME", "S_REPLAY", "S_DRAIN", "S_DONE"):
        assert state in verilog
    with (ROOT / "results" / "rtl_state_occupancy.csv").open() as f:
        clear = [int(r["cycles"]) for r in csv.DictReader(f) if r["testbench"] == "psd_top" and r["state"] == "T_CLR"]
    assert clear == [1027] and "(1 027 cycles, hidden behind the capture)" in flat


def test_the_bench_table_equals_what_the_build_script_sets(text):
    build = (MATLAB / "bt_mc_build.m").read_text()
    rows = {r[0]: r for r in table(text, "| Block | Parameters | Wired to |")}
    ins, outs = black_box_ports()
    p = dict(re.findall(r"p\.(\w+) = ([^;]+);", (MATLAB / "bt_mc_params.m").read_text()))
    wl = int(p["wl"])

    sine = rows["Sine Wave (Simulink Sources)"][1]
    assert "'SineType', 'Time based'" in build and "sine type Time based" in sine
    assert "addParameter(ip, 'ToneAmp', 0.5)" in build and "'Amplitude', num2str(o.ToneAmp)" in build and "amplitude 0.5" in sine
    assert "addParameter(ip, 'ToneBin', 100.37)" in build and "sprintf('2*pi*%.10g/1024', o.ToneBin)" in build
    assert "frequency 2*pi*100.37/1024" in sine
    assert "'Phase', '0'" in build and "'SampleTime', '1'" in build and "phase 0, sample time 1" in sine
    noise = rows["Random Number (Simulink Sources)"][1]
    assert "addParameter(ip, 'NoiseRms', 0.05)" in build and "sprintf('%.10g', o.NoiseRms^2)" in build
    assert "'Mean', '0'" in build and "'Seed', '1'" in build and noise == "mean 0, variance 0.05^2, seed 1, sample time 1"
    assert "'Inputs', '++'" in build and rows["Sum"][1:] == ["++", "gw_in_data"]

    gates = {m[0]: (m[1], m[2].strip(), m[3].strip()) for m in
             re.findall(r"gw_in\('(\w+)', rect\([^)]*\), '((?:[^']|'')+)', ([^,]+), ([^)]+)\)", build)}
    assert set(gates) == {"gw_in_data", "gw_in_valid", "gw_in_rst", "gw_in_addr"}
    arith, nbits, binpt = gates["gw_in_data"]
    assert arith.startswith("Signed") and nbits == "p.wl" and binpt == "p.wl - 1"
    assert "set_blk(name, 'quantization', 'Round  (unbiased: +/- Inf)')" in build and "set_blk(name, 'overflow', 'Saturate')" in build
    assert "set_blk(name, 'period', '1')" in build
    assert rows["Gateway In gw_in_data"][1] == (f"signed (two's complement), {wl} bits, binary point {wl - 1}, quantization Round "
                                                "(unbiased: +/- Inf), overflow Saturate, sample period 1")
    assert gates["gw_in_valid"][0] == gates["gw_in_rst"][0] == "Boolean"
    assert rows["Constant 1, Gateway In gw_in_valid"][1] == rows["Constant 0, Gateway In gw_in_rst"][1] == "Boolean"
    assert gates["gw_in_addr"] == ("Unsigned", "10", "0")
    assert rows["From Workspace bt_addr_stim, Gateway In gw_in_addr"][1] == "unsigned, 10 bits, binary point 0"
    assert "'VariableName', 'bt_addr_stim'" in build and "'Value', '1'" in build and "'Value', '0'" in build

    bb = rows["Black Box bt_psd_mc"]
    assert "add_blk('xbsIndex_r4/Black Box', 'bt_psd_mc'" in build and "'config_file', fullfile(here, 'bt_psd_mc_config.m')" in build
    assert "set_sim_mode('bt_psd_mc', 'Vivado Simulator')" in build
    assert bb[1] == "configuration M-function modelcomposer/matlab/bt_psd_mc_config.m; Simulation mode: Vivado Simulator"

    taps = re.findall(r"'([\w/]+)', '(\w+)', '(bt_\w+)'", build)
    assert [t[1] for t in taps] == ["busy", "done", "rd_data", "in_data", "in_valid", "rd_addr"] and [t[0] for t in taps[:3]] == [
        f"bt_psd_mc/{i}" for i in (1, 2, 3)]
    assert [t[1] for t in taps[:3]] == outs
    log = rows["Gateway Out (six), To Workspace (six)"]
    assert log[1] == "one tap each on busy, done, rd_data, and on the gateway outputs of in_data, in_valid, rd_addr"
    assert log[2] == "variables " + ", ".join(t[2] for t in taps)
    assert ins == ["rst", "in_valid", "in_data", "rd_addr"]
    assert "'Model Composer Hub|System Generator'" in build and rows["Model Composer Hub"][1] == "sample period (Simulink system period) 1"


def test_the_solver_statement_equals_the_build_script(flat):
    build = (MATLAB / "bt_mc_build.m").read_text()
    cases = dict(re.findall(r"case '(\w+)'.*?solver = '(\w+)';", build, re.S))
    assert cases == {"sine": "FixedStepDiscrete", "signal": "ode1", "chirp": "ode1", "file": "FixedStepDiscrete"}
    assert "'SolverType', 'Fixed-step', 'Solver', solver, 'FixedStep', '1'" in build
    assert ("the solver is fixed-step with step size 1 (`FixedStepDiscrete` for the discrete sources, `ode1` for Signal Generator "
            "and Chirp Signal, which are continuous-time blocks)") in flat


def test_the_source_table_equals_the_build_script(text, flat):
    build = (MATLAB / "bt_mc_build.m").read_text()
    for opt in ("Name", "File", "ToneBin", "ToneAmp", "NoiseRms", "Waveform"):
        assert f"addParameter(ip, '{opt}'" in build and f"`'{opt}'`" in text, opt
    assert "'WaveForm', o.Waveform" in build and "sprintf('%.10g/1024', o.ToneBin), 'Units', 'Hertz'" in build
    assert "'f1', '0.01'" in build and "'f2', '0.45'" in build and "'T', num2str(p.nsamp)" in build
    assert "0.01 to 0.45 cycles per sample over the record (bins 10.24 to 460.8)" in flat
    assert (0.01 * 1024, 0.45 * 1024) == pytest.approx((10.24, 460.8))
    assert "'VariableName', 'bt_stim'" in build and "From Workspace" in text
    rows = table(text, "| `bt_mc_build(...)` |")
    assert [r[0] for r in rows] == ["'sine' (default)", "'signal'", "'chirp'", "'file'"]
    for case in ("sine", "signal", "chirp", "file"):
        assert f"case '{case}'" in build


def test_the_run_description_equals_the_run_script(text, flat):
    run = (MATLAB / "bt_mc_run.m").read_text()
    assert "[t mod(t, 1024)]" in run and "nnz(done ~= 0) >= 1100" in run and "for attempt = 1:6" in run and "stop = 2 * stop" in run
    assert "repeats every 1 024 clock cycles" in flat and "up to six attempts" in flat and "at least 1 100 cycles" in flat
    assert "bt_log.mat" in run and "modelcomposer/work/bt_log.mat" in text and "struct `bt_log`" in text
    logged = ("in_code", "in_valid", "rd_addr", "rd_data", "done", "busy", "wl", "nchunks")
    assert re.findall(r"res\.(\w+) =", run) == list(logged)
    after = flat[flat.index("struct `bt_log`"):]
    positions = [after.index(f"`{k}`") for k in logged]
    assert positions == sorted(positions), "the guide lists the fields of the log in the order the script writes them"


def test_checker_messages_in_the_text_exist_in_the_checker(text):
    sources = (ROOT / "python" / "btpsd" / "mcmodel.py").read_text() + (ROOT / "python" / "scripts" / "compare_mc_result.py").read_text()
    rows = table(text, "| Message | Meaning |")
    assert len(rows) == 5
    for stem in ('head = "PASS" if self.ok else "FAIL"', "of {self.n_bins} bins equal the bit-exact model",
                 "done after {self.cycles_to_done} cycles (predicted {self.predicted_cycles})",
                 "done was never high long enough to read a result", "bins never read", "two reads of bin"):
        assert stem in sources, stem
    assert rows[0][0].startswith("PASS: 1024 of 1024 bins equal the bit-exact model; done after <n> cycles (predicted <m>)")
    assert rows[1][0].startswith("FAIL: <k> of 1024 bins equal the bit-exact model")
    assert rows[2][0] == "done was never high long enough to read a result" and rows[3][0] == "<k> bins never read"
    assert rows[4][0] == "two reads of bin <b> differ"
    assert "keep the run going at least 1 026 cycles after done" in rows[3][1] and 1024 + 1 + 1 == 1026
    assert re.search(r"cycles after done\"|at least \{bins \+ lag \+ 1\} cycles after done", sources)
    # exit codes
    main = (ROOT / "python" / "scripts" / "compare_mc_result.py").read_text()
    assert main.count("return 2") == 1 and "return 0 if rc.ok else 1" in main and "return 1" in main
    assert "Exit status 0 means all words equal the model, 1 that they differ or the log is unusable, 2 that the file cannot be read." in re.sub(
        r"\s+", " ", text)


def test_stimulus_sources_in_the_text_are_those_of_the_exporter(text):
    generators = ("white", "pink", "ar2", "sines", "eeg_like", "eeg_coupled")
    assert set(generators) == set(signals.SIGNALS)
    assert set(ex.SOURCES) == set(generators) | {"tone", "siggen"}
    for s in ex.SOURCES:
        assert f"`{s}`" in text, s
    assert "expected_P_int" in (ROOT / "python" / "scripts" / "export_mc_stimulus.py").read_text() and "`expected_P_int`" in text


# ---------------------------------------------------------------------------------------------------------------
# generator numbers
# ---------------------------------------------------------------------------------------------------------------
def test_phase_increment_table_equals_the_model(text):
    rows = table(text, "| tone bin | `pinc` |")
    assert [r[0] for r in rows] == ["64", "100", "100.37 (the self-test)", "256"]
    for cell_bin, cell_pinc in rows:
        assert num(cell_pinc) == siggen.pinc_for(float(re.match(r"[\d.]+", cell_bin).group())), cell_bin
    assert siggen.pinc_for(siggen.SELFTEST["bin_1024"]) == 420982292
    assert siggen.tone_bin(siggen.pinc_for(100.37)) == pytest.approx(100.37, abs=1e-6)


def test_generator_constants_in_the_text_equal_the_sources(text, flat):
    verilog = SIG_GEN_V.read_text().upper()
    assert "80200003" in verilog and "A5A51F3D" in verilog
    assert siggen.TAPS == 0x80200003 and siggen.SEED == 0xA5A51F3D and "`0x80200003`" in text and "`0xA5A51F3D`" in text
    assert "x^32 + x^22 + x^2 + x + 1" in text and siggen.TONE_AMP == 0.5 and siggen.PHW == 32 and siggen.ROM_BITS == 10
    assert "phase[n] >> 22" in text and siggen.PHW - siggen.ROM_BITS == 22 and "the 10 top bits of the phase" in flat
    assert "round(0.5 * 2^(WL-1) * sin(2 pi k / 1024))" in text
    rom = siggen.rom_table(12)
    assert len(rom) == 1024 and rom[0] == 0 and rom.max() == 1024 and rom.min() == -1024
    peak = 2 ** (12 - 1 - siggen.SELFTEST["noise_shift"])
    assert peak == 256 and peak / 2 ** 11 == 0.125
    assert "peak `2^(WL-1-noise_shift)` codes: 256 codes, an eighth of full scale, for `WL` = 12 and `noise_shift` = 3" in flat
    assert "noise peak 0.125 of full scale" in flat
    assert "`valid` rises after two enabled clock edges" in flat


# ---------------------------------------------------------------------------------------------------------------
# spectrum and cycle numbers
# ---------------------------------------------------------------------------------------------------------------
def test_spectrum_floor_table_equals_the_model(text, flat):
    rec = siggen.selftest_record()
    cfg = M.McConfig(wl=rec["wl"], nchunks=rec["nchunks"])
    x = rec["x"]
    P = np.maximum(M.expected(x, cfg).P_int.astype(float), 1.0)
    curves = [psd_exact_welch(x.astype(float)), psd_float_fd(x.astype(float)), P]

    def db(v):
        return 10.0 * np.log10(np.maximum(v / np.max(v), 1e-12))

    rows = table(text, "| | exact Welch average |")
    assert [r[0] for r in rows] == ["bins 0 to 39", "bins 40 to 89"]
    for (lo, hi), row in zip(((0, 40), (40, 90)), rows):
        assert [float(c) for c in row[1:]] == [round(float(np.mean(db(c[:512])[lo:hi])), 1) for c in curves], row
    assert int(np.argmax(P[:512])) == 100 and 512 + int(np.argmax(P[512:])) == 924
    assert "the words peak at bin 100, and at its mirror bin 924 in the two-sided spectrum" in flat
    assert (rec["wl"], rec["nchunks"], rec["noise_shift"]) == (12, 4, 3) and "four chunks" in flat


def test_cycle_budget_table_equals_the_timeline(text, timelines):
    T = timelines[5]
    assert T["pred"] == 145790
    rows = table(text, "| phase | cycles |")
    assert [num(r[1]) for r in rows] == [v for _, v in T["phases"]] + [T["pred"]]
    assert sum(v for _, v in T["phases"]) == T["pred"]
    S, res = T["res"].S, T["res"]
    st = M.psd_engine_state_cycles(S, res.squarer_cycles, clears=0) if hasattr(M, "psd_engine_state_cycles") else None
    if st is None:
        from btpsd.rtlmodel import psd_engine_state_cycles

        st = psd_engine_state_cycles(S, res.squarer_cycles, clears=0)
    merge = sum(st[k] for k in ("M_A", "M_B", "M_C", "M_F"))
    window = sum(st[k] for k in ("W_A", "W_B", "W_C", "W_F", "W_G", "W_ACC"))
    assert merge == 9728 * S and window == 18432 * S and st["W_S"] == 3072 * S + res.squarer_cycles and res.squarer_cycles == 14683
    assert st["S_END"] + S * M.T_SEG_EXTRA + 2 == 14 and (st["S_END"] // S + M.T_SEG_EXTRA) == 3
    assert (M.CHUNK, M.T_FFT_PER_CHUNK) == (512, 209)
    labels = [r[0] for r in rows]
    assert labels == ["capture", "replay: load (512 per chunk)", "replay: FFT (209 per chunk)", "segment merge (9 728 per segment)",
                      "window and accumulate (18 432 per segment)",
                      "squarer (3 072 per segment plus 14 683 data-dependent cycles in all)",
                      "segment hand-offs (3 per segment) and 2", "to done"]


def test_chunk_table_equals_the_model_for_the_generator_record(text, flat, timelines):
    rows = table(text, "| chunks | segments |")
    assert [int(r[0]) for r in rows] == [2, 3, 4, 5, 9]
    for r in rows:
        T = timelines[int(r[0])]
        assert int(r[1]) == T["cfg"].segments
        assert num(r[2]) == T["pred"] == M.expected_cycles(T["res"], T["cfg"])
        assert num(r[3]) == M.recommended_stop(T["res"], T["cfg"])
    assert M.READOUT_MARGIN == 2048
    assert "the prediction plus 2 % plus a read-out window of 2 048 cycles" in flat


def test_gaussian_tone_record_prediction_in_the_text(flat):
    cfg = M.McConfig(wl=12, nchunks=5)
    pred = M.expected_cycles(M.expected(ex.make_record("tone", cfg), cfg), cfg)
    assert pred == 143274 and "gives 143 274 for the same build" in flat


def test_measured_cycle_table_equals_the_stored_results_and_the_predictions(text):
    rows = table(text, "| Testbench | Record | Cycles simulated | Cycles predicted |")
    with (ROOT / "results" / "rtl_mc_cycles.csv").open() as f:
        measured = {(r["testbench"], int(r["run"])): int(r["cycles"]) for r in csv.DictReader(f)}
    predicted = M.testbench_cycle_predictions()
    keys = [("tb_bt_psd_mc", 0), ("tb_bt_psd_mc", 1), ("tb_bt_psd_selftest", 0)]
    assert set(measured) == set(predicted) == set(keys) and len(rows) == 3
    for key, row in zip(keys, rows):
        assert row[0].startswith(key[0]) and (key[0] != "tb_bt_psd_mc" or row[0].endswith(f"run {key[1]}")), row
        assert num(row[2]) == measured[key] and num(row[3]) == predicted[key], row
    prov = json.loads((ROOT / "results" / "rtl_sim_provenance.json").read_text())
    flat = re.sub(r"\s+", " ", text).replace("`", "")
    assert f"CI run {prov['run_id']}" in flat and f"commit {prov['commit'][:7]}" in flat
    assert {"tb_bt_psd_mc", "tb_bt_psd_selftest", "tb_sig_gen"} <= set(prov["testbenches_passed"])


def test_the_guide_quotes_no_value_of_the_paper_s_comparison_table(text, flat):
    """The cycle counts are not compared with Table 9; the paper-reported figures stay in the files that label them."""
    assert text.count("Table 9") == 1 and "not comparable with Table 9 of the paper" in flat
    values = set()
    for r in paperdata.table9():
        for key in ("latency_cycles", "pdp", "lut", "ff"):
            if r.get(key) is not None and r[key] >= 1000:
                values.add(int(r[key]))
    spelled = {f"{v:,}".replace(",", " ") for v in values} | {str(v) for v in values}
    tokens = set(re.findall(r"(?<![\d.])\d{1,3}(?: \d{3})+(?![\d])|(?<![\d.])\d{4,}(?![\d])", text))
    assert not (tokens & spelled), tokens & spelled


# ---------------------------------------------------------------------------------------------------------------
# honesty of the status statements
# ---------------------------------------------------------------------------------------------------------------
def test_the_guide_states_what_was_not_executed(text, flat):
    assert "**not executed**" in text and "has **not been built or run here**" in flat
    for needle in ("MATLAB, Simulink, Model Composer and Vivado were not available where it was written",
                   "**The Simulink flow was not executed.**", "arith_type", "n_bits", "bin_pt", "config_file",
                   "the string of the simulation mode and the name of the Hub block", "simulated, not synthesised",
                   "(not executed here; `sim/xsim_smoke.sh`"):
        assert needle in flat, needle
    for path in ("sim/xsim_smoke.sh", "modelcomposer/matlab/bt_mc_build.m", "modelcomposer/matlab/bt_mc_run.m"):
        assert "NOT EXECUTED" in (ROOT / path).read_text().upper()


def test_the_unexecuted_parameter_names_are_those_the_build_script_sets(flat):
    build = (MATLAB / "bt_mc_build.m").read_text()
    for name in ("arith_type", "n_bits", "bin_pt", "quantization", "overflow", "period", "config_file"):
        assert f"'{name}'" in build and f"`{name}`" in flat or name == "config_file" and "config_file" in flat, name
    assert "(`arith_type`, `n_bits`, `bin_pt`, `quantization`, `overflow`, `period`)" in flat


def test_simulation_mode_statements_match_the_documentation_wording(flat):
    assert "In the mode *Inactive* a Black Box ignores its inputs and outputs zeros (UG958)" in flat
    assert "Select *Vivado Simulator* (UG1483, \"Configuring the HDL Simulator\")" in flat
    assert "The third option, *External co-simulator*, is for a Questa session." in flat
    assert "ug1483-model-composer-sys-gen-user-guide/Configuring-the-HDL-Simulator" in flat
    assert "ug958-vivado-sysgen-ref/Black-Box" in flat


def test_dds_recipe_is_marked_as_untried_and_claims_no_sfdr_figure(text, flat):
    section = text[text.index("### 6.3"):text.index("## 7.")]
    assert "**not been built or run here**" in re.sub(r"\s+", " ", section)
    assert "documentation only" in text[:text.index("## 1.")] and "no figure is claimed here" in re.sub(r"\s+", " ", section)
    assert not re.search(r"\d\s*dB", section), "no dB figure for the DDS Compiler is quoted"
    assert "f = f_clk * PINC / 2^32" in section and "`PINC = round(bin / 1024 * 2^32)`" in section
    assert "UG958 says to use it instead of the Sine Wave block when the design needs a phase generator, several channels or AXI4 ports" in flat


def test_requirements_statements(text, flat):
    assert "model_composer -configlib" in text and "UG1483 says it needs sudo access" in text
    assert "MATLAB and Simulink at a release that your Vitis version supports" in flat
    assert "licensed as an add-on to them" in flat
    assert "xbsIndex_r4" in text and "xbsIndex_r4" in (MATLAB / "bt_mc_build.m").read_text()
    assert "stops with a message if MATLAB was not started through Model Composer" in flat and "bt_mc_build:library" in (
        MATLAB / "bt_mc_build.m").read_text()


def test_pack_readme_lists_every_file_of_the_pack_and_points_to_the_guide():
    readme = PACK_README.read_text(encoding="utf-8")
    for p in sorted((ROOT / "modelcomposer").rglob("*")):
        rel = p.relative_to(ROOT / "modelcomposer").as_posix()
        if p.is_file() and rel != "README.md" and not rel.startswith(("stimulus/", "work/")):
            assert f"`{rel}`" in readme, rel
    assert "(../docs/MODEL_COMPOSER.md)" in readme and "**not executed**" in readme
    assert (ROOT / "docs" / "MODEL_COMPOSER.md").is_file()
