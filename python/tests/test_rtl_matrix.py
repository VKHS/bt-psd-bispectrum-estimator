"""The configuration matrix of the RTL simulation: what is run, how it is generated, and what the CI run found."""
import csv
import json
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python" / "scripts"))

import run_rtl_matrix as M  # noqa: E402

RESULTS = ROOT / "results"


def _makefile_lists() -> dict[str, list[str]]:
    """Expand the file-list variables of sim/Makefile (``NAME := a $(OTHER) b``)."""
    raw = {m.group(1): m.group(2).split()
           for m in re.finditer(r"^(\w+)\s*:=\s*(.*)$", (ROOT / "sim" / "Makefile").read_text(), re.M)}

    def expand(tokens):
        out = []
        for t in tokens:
            ref = re.fullmatch(r"\$\((\w+)\)", t)
            out += expand(raw[ref.group(1)]) if ref else [t]
        return out

    return {k: expand(v) for k, v in raw.items() if k in ("CRU", "BFLY", "FFT", "PSD", "BSP")}


def test_driver_compiles_the_same_files_as_the_makefile():
    lists = _makefile_lists()
    assert M.FFT_FILES == lists["FFT"]
    assert M.PSD_FILES == ["tb/tb_psd_top.v", *lists["PSD"]]
    assert M.BSP_FILES == ["tb/tb_bispectrum_top.v", *lists["BSP"]]


def test_configurations_cover_the_rtl_parameters():
    cfgs = M.configurations()
    assert len({c.tag for c in cfgs}) == len(cfgs)
    for kind in ("psd", "bsp"):
        sel = [c for c in cfgs if c.kind == kind]
        assert {c.wl for c in sel} >= {10, 12, 14, 16}
        assert {c.window for c in sel} == {"rect", "hann", "hamming", "blackman"}
        assert {c.rnd for c in sel} == {0, 1}
    assert {c.serial for c in cfgs if c.kind == "psd"} == {0, 1}


@pytest.mark.parametrize("n", [1, 3, 4])
def test_shards_partition_the_matrix_and_balance_the_load(n):
    cfgs = M.configurations()
    parts = [M.shard(cfgs, i, n) for i in range(1, n + 1)]
    tags = [c.tag for part in parts for c in part]
    assert sorted(tags) == sorted(c.tag for c in cfgs) and len(set(tags)) == len(tags)
    assert all([c.tag for c in part] == [c.tag for c in cfgs if c in part] for part in parts)  # canonical order
    loads = [sum(M.cost(c) for c in part) for part in parts]
    assert max(loads) <= 1.1 * sum(loads) / n


@pytest.mark.parametrize("window", ["rect", "hann", "hamming", "blackman"])
def test_model_window_taps_equal_the_rtl_header(window):
    taps = M.window_taps(window, 10)  # asserts the equality with rtl/btpsd_consts.vh itself
    assert len(taps) == 5 and taps[2] > 0


def test_generated_include_and_vectors_of_one_configuration_per_testbench():
    for c, words in ((M.Config("psd", 10, "hamming", rnd=0, serial=0), 1024), (M.Config("bsp", 14, "blackman"), 32 * 32)):
        d, _ = M.make_vectors(c)
        inc = (M.WORK / f"{c.tag}.vh").read_text()
        macros = set(re.findall(r"`define (\w+)", inc))
        needed = {"PSD_S", "PSD_AVG_SHIFT", "PSD_NCHUNKS", "PSD_SQUARER_CYCLES", "TB_XFILE", "TB_PFILE"} if c.kind == "psd" \
            else {"BSP_S", "BSP_AVG_SHIFT", "BSP_MAX_BIN", "BSP_NSAMPLES", "TB_XFILE", "TB_REFILE", "TB_IMFILE"}
        assert needed <= macros and {"TB_WL", "TB_PARAMS"} <= macros
        assert f".WL({c.wl})" in inc and f".ROUND({c.rnd})" in inc
        assert (".SQ_SERIAL(0)" in inc) == (c.kind == "psd")
        lines = lambda name: (M.WORK / f"{c.tag}_{name}.mem").read_text().split()
        assert len(lines("x")) == len(d["x"])
        if c.kind == "psd":
            assert len(lines("P")) == words and all(len(v) == (2 * c.wl + 8 + 3) // 4 for v in lines("P"))
        else:
            assert len(lines("B_re")) == len(lines("B_im")) == words


def test_testbench_defaults_are_the_committed_configuration():
    for name, vh, pre in (("tb_psd_top.v", "psd_wl12.vh", "PSD"), ("tb_bispectrum_top.v", "bsp_wl12.vh", "BSP")):
        text = (ROOT / "tb" / name).read_text()
        assert f'`define TB_VH "{vh}"' in text and "`define TB_WL 12" in text
        assert (ROOT / "sim" / "vectors" / vh).exists()
        # nothing of the committed configuration is hard-wired any more: files and parameters come from macros
        assert "$readmemh(\"" not in text and "#(.WL(WL)" not in text


@pytest.mark.skipif(not (RESULTS / "rtl_config_matrix.csv").exists(), reason="no CI matrix result stored")
def test_stored_matrix_result_is_the_full_matrix_and_passes():
    rows = list(csv.DictReader((RESULTS / "rtl_config_matrix.csv").open()))
    assert [r["tag"] for r in rows] == [c.tag for c in M.configurations()]
    assert list(rows[0].keys()) == M.CSV_HEADER
    assert all(r["result"] == "PASS" for r in rows)
    for r, c in zip(rows, M.configurations()):
        assert int(r["words_checked"]) == 1024 and int(r["wl"]) == c.wl and r["window"] == c.window
        assert int(r["round"]) == c.rnd and r["squarer"] == c.squarer
        # the clamp of the half difference can only fire with rounding, and the vectors are made to reach it
        assert (int(r["model_clamped_diffs"]) > 0) == bool(c.rnd)
    prov = json.loads((RESULTS / "rtl_config_matrix_provenance.json").read_text())
    assert prov["configurations"] == prov["passed"] == len(rows)
