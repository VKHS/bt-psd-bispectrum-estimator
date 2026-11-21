"""Tests of the smoke test itself: it passes on a healthy checkout, and it fails when the things it guards go stale."""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from btpsd import smoke

ROOT = Path(__file__).resolve().parents[2]
FAST = ["environment", "layout", "generated_rtl", "schedule", "fft_model", "psd_model", "bispectrum_model", "model_composer_files",
        "model_composer_model", "verilog_lint", "figure_framework"]


def _by_name(results):
    return {r.name: r for r in results}


def test_fast_checks_pass_on_the_checkout():
    res = _by_name(smoke.run_checks(FAST))
    assert set(res) == set(FAST)
    for name, r in res.items():
        assert r.status in ("PASS", "SKIP"), (name, r.detail)
    assert res["schedule"].status == "PASS" and "207" in res["schedule"].detail
    assert res["fft_model"].status == "PASS" and "bit-exact" in res["fft_model"].detail
    assert res["model_composer_files"].status == "PASS" and "current" in res["model_composer_files"].detail
    assert res["model_composer_model"].status == "PASS" and "one wrong word is caught" in res["model_composer_model"].detail


def test_golden_vectors_check_passes():
    r = _by_name(smoke.run_checks(["golden_vectors"]))["golden_vectors"]
    assert r.status == "PASS", r.detail
    assert "identical" in r.detail


def test_every_check_is_registered_once_and_documented():
    names = [n for n, _, _ in smoke.CHECKS]
    assert len(names) == len(set(names))
    assert {"environment", "layout", "generated_rtl", "golden_vectors", "schedule", "fft_model", "psd_model", "bispectrum_model",
            "model_composer_files", "model_composer_model", "verilog_lint", "icarus_testbenches", "xsim_testbenches",
            "figure_framework", "figures_present"} <= set(names)
    for _, _, meta in smoke.CHECKS:
        assert meta["doc"], "a check without a docstring"


def test_quick_skips_slow_checks_and_xsim_is_opt_in():
    res = _by_name(smoke.run_checks(None, quick=True, xsim=False))
    assert "xsim_testbenches" not in res  # optional checks only run on request
    assert res["golden_vectors"].status == "SKIP" and res["icarus_testbenches"].status == "SKIP"
    assert "xsim_testbenches" in _by_name(smoke.run_checks(["xsim_testbenches"]))


def test_icarus_check_skips_without_the_tools_and_never_reports_a_pass_it_did_not_run():
    r = _by_name(smoke.run_checks(["icarus_testbenches"]))["icarus_testbenches"]
    if all(shutil.which(t) for t in ("iverilog", "vvp", "make")):
        assert r.status in ("PASS", "FAIL")
    else:
        assert r.status == "SKIP" and "not found" in r.detail


def test_xsim_check_fails_loudly_without_vivado():
    if shutil.which("xvlog"):
        pytest.skip("Vivado is installed here")
    r = _by_name(smoke.run_checks(["xsim_testbenches"]))["xsim_testbenches"]
    assert r.status == "FAIL" and "settings64.sh" in r.detail  # requested but impossible: a failure, not a skip


def test_a_stale_generated_file_is_detected(monkeypatch, tmp_path):
    g = smoke._script("gen_rtl")
    monkeypatch.setattr(g, "consts_vh", lambda: "// stale\n")
    r = _by_name(smoke.run_checks(["generated_rtl"]))["generated_rtl"]
    assert r.status == "FAIL" and "btpsd_consts.vh" in r.detail


def test_a_stale_vector_is_detected(monkeypatch):
    g = smoke._script("gen_vectors")
    real = g.write_all

    def tampered(out):
        s = real(out)
        p = Path(out) / "cru_wl16_out.mem"
        p.write_text(p.read_text().replace("\n", "\n0", 1))
        return s

    monkeypatch.setattr(g, "write_all", tampered)
    r = _by_name(smoke.run_checks(["golden_vectors"]))["golden_vectors"]
    assert r.status == "FAIL" and "cru_wl16_out.mem" in r.detail


def test_a_wrong_latency_constant_is_detected(monkeypatch):
    from btpsd import rtlmodel

    monkeypatch.setattr(rtlmodel, "LAT_BFLY", rtlmodel.LAT_BFLY + 1)
    r = _by_name(smoke.run_checks(["schedule"]))["schedule"]
    assert r.status == "FAIL" and "btpsd_params.vh" in r.detail


def test_exceptions_become_failures_not_crashes(monkeypatch):
    def boom():
        raise RuntimeError("boom")

    monkeypatch.setattr(smoke, "CHECKS", [("boom", boom, {"slow": False, "optional": False, "doc": "x"})])
    r = smoke.run_checks()[0]
    assert r.status == "FAIL" and "RuntimeError: boom" in r.detail


def test_command_line_exit_status_and_listing():
    p = subprocess.run([sys.executable, str(ROOT / "python" / "scripts" / "smoke_test.py"), "--list"], capture_output=True, text=True)
    assert p.returncode == 0 and "fft_model" in p.stdout and "[slow, optional]" in p.stdout
    p = subprocess.run([sys.executable, str(ROOT / "python" / "scripts" / "smoke_test.py"), "--only", "nonsense"], capture_output=True, text=True)
    assert p.returncode == 2 and "unknown check" in p.stdout
    p = subprocess.run([sys.executable, str(ROOT / "python" / "scripts" / "smoke_test.py"), "--only", "schedule,fft_model"],
                       capture_output=True, text=True)
    assert p.returncode == 0 and "2 passed, 0 failed" in p.stdout


def test_xsim_script_is_well_formed():
    s = (ROOT / "sim" / "xsim_smoke.sh").read_text()
    assert s.startswith("#!/usr/bin/env bash") and "NOT EXECUTED" in s
    for tb in ("tb_cru_rotate", "tb_r8_butterfly", "tb_mag_sq_odd", "tb_fft512_mem", "tb_psd_top", "tb_bispectrum_top",
               "tb_bt_psd_mc", "tb_sig_gen", "tb_bt_psd_selftest"):
        assert tb in s
    # the same source lists as sim/Makefile
    mk = (ROOT / "sim" / "Makefile").read_text()
    for token in ("rtl/ram_1r1w.v", "rtl/r8_butterfly.v", "rtl/fft512_mem.v", "rtl/mag_sq_odd.v", "rtl/psd_engine.v", "rtl/bt_psd_top.v",
                  "rtl/bispectrum_engine.v", "rtl/bt_bispectrum_top.v", "modelcomposer/hdl/bt_psd_mc.v",
                  "modelcomposer/hdl/bt_psd_mc_bundle.v", "modelcomposer/hdl/sig_gen.v", "modelcomposer/hdl/sig_rom.v",
                  "modelcomposer/hdl/bt_psd_selftest.v"):
        assert token in s and token in mk
    bash = shutil.which("bash")
    if bash:
        assert subprocess.run([bash, "-n", str(ROOT / "sim" / "xsim_smoke.sh")]).returncode == 0
        p = subprocess.run([bash, str(ROOT / "sim" / "xsim_smoke.sh")], capture_output=True, text=True, env={"PATH": "/nonexistent"})
        assert p.returncode == 2 and "settings64.sh" in p.stderr


def test_smoke_testbench_list_matches_the_makefile_and_the_xsim_script():
    import re

    mk = (ROOT / "sim" / "Makefile").read_text()
    all_line = re.search(r"(?m)^all:\s*(.+)$", mk).group(1).split()
    assert tuple(all_line) == smoke.TESTBENCHES
    for t in smoke.TESTBENCHES:
        assert re.search(rf"(?m)^{t}:\s", mk), f"no target {t} in sim/Makefile"
    xs = (ROOT / "sim" / "xsim_smoke.sh").read_text()
    default = re.search(r"\|\| TARGETS=\(([^)]*)\)", xs).group(1).split()
    assert tuple(default) == smoke.TESTBENCHES
    for t in smoke.TESTBENCHES:
        assert f"[{t}]=" in xs


def test_a_stale_model_composer_file_is_detected(monkeypatch):
    g = smoke._script("gen_mc")
    real = g.bundle_text
    monkeypatch.setattr(g, "bundle_text", lambda *a, **k: real(*a, **k) + "// stale\n")
    r = _by_name(smoke.run_checks(["model_composer_files"]))["model_composer_files"]
    assert r.status == "FAIL" and "bt_psd_mc_bundle.v" in r.detail and "gen_mc.py" in r.detail


def test_a_wrong_generator_constant_is_detected(monkeypatch):
    from btpsd import siggen

    monkeypatch.setattr(siggen, "TAPS", 0x80200002)  # not a maximal-length polynomial: the state sequence repeats early
    r = _by_name(smoke.run_checks(["model_composer_model"]))["model_composer_model"]
    assert r.status == "FAIL" and "LFSR repeats" in r.detail, r.detail
