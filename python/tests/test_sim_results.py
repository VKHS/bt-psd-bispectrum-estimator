"""The committed RTL-simulation results (``results/rtl_*.csv``) against the model formulas.

The CSV files are what the Icarus Verilog run of the testbenches in CI printed (``scripts/ingest_ci_traces.py``).  The
tests check that the simulated machine does what the model of the RTL says: the issue schedule of the FFT core, the
handshake of the squarer, and the cycles per state of both engines.
"""
import csv
import json
from pathlib import Path

import pytest

from btpsd import rtlmodel
from btpsd.rtlparse import define_values

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"
VECTORS = ROOT / "sim" / "vectors"


def _rows(name):
    with (RESULTS / name).open() as f:
        return list(csv.DictReader(f))


def test_provenance():
    p = json.loads((RESULTS / "rtl_sim_provenance.json").read_text())
    assert p["testbenches_passed"] == ["tb_bispectrum_top", "tb_cru_rotate", "tb_fft512_mem", "tb_mag_sq_odd", "tb_psd_top",
                                       "tb_r8_butterfly"]
    assert len(p["commit"]) == 40 and p["run_url"].endswith(str(p["run_id"]))
    assert p["fft512_first_issue_to_done_cycles"] == [rtlmodel.LAT - 1 + 3 * 64 + rtlmodel.stall_cycles(rtlmodel.LAT)]


def test_fft_trace_follows_the_issue_schedule():
    rows = [{k: (int(v) if k != "state" else v) for k, v in r.items()} for r in _rows("rtl_trace_fft512.csv")]
    by_cycle = {r["cycle"]: r for r in rows}
    stall = rtlmodel.stall_cycles(rtlmodel.LAT)
    sched = rtlmodel.Fft512Rtl().issue_schedule(stall)
    assert by_cycle[0]["start"] == 1 and by_cycle[-1]["start"] == 0
    # the first butterfly issues in the cycle after the start cycle; busy follows the issue controller
    issued = {r["cycle"] - 1: (r["iss_p"], r["iss_b"]) for r in rows if r["iss_v"]}
    assert issued == {t: (p, b) for t, p, b in sched}
    last = sched[-1][0]
    assert [c for c, r in by_cycle.items() if r["busy"]] == list(range(1, 1 + last + rtlmodel.LAT))
    # done is a one-cycle pulse, LAT cycles after the last issue; the run takes 207 cycles from the first issue
    done = [c for c, r in by_cycle.items() if r["done"]]
    assert done == [1 + last + rtlmodel.LAT]
    assert done[0] - 1 == define_values(VECTORS / "fft_wl16.vh")["FFT_TOTAL_CYCLES"] == 207
    assert {r["state"] for r in rows if 1 <= r["cycle"] <= 207} == {"ST_RUN", "ST_DRAIN"}


@pytest.mark.parametrize("example, re, im", [(0, 3, -5), (1, 0, 0), (2, -4, 2)])
def test_squarer_handshake(example, re, im):
    rows = [{k: int(v) for k, v in r.items()} for r in _rows("rtl_trace_squarer.csv") if int(r["example"]) == example]
    n = max(abs(re), abs(im))
    assert [r["cycle"] for r in rows] == list(range(14))
    assert [r["cycle"] for r in rows if r["start"]] == [0]
    done = [r["cycle"] for r in rows if r["done"]]
    assert done == [n + 2]                                              # done is visible n + 2 cycles after start
    assert [r["cycle"] for r in rows if r["busy"]] == list(range(1, n + 2))
    assert rows[done[0]]["p"] == re * re + im * im


def test_psd_state_occupancy():
    d = define_values(VECTORS / "psd_wl12.vh")
    occ = {r["state"]: int(r["cycles"]) for r in _rows("rtl_state_occupancy.csv") if r["testbench"] == "psd_engine"}
    want = rtlmodel.psd_engine_state_cycles(d["PSD_S"], d["PSD_SQUARER_CYCLES"])
    assert {k: occ[k] for k in want} == want
    top = {r["state"]: int(r["cycles"]) for r in _rows("rtl_state_occupancy.csv") if r["testbench"] == "psd_top"}
    assert top["T_LOAD"] == d["PSD_NCHUNKS"] * 512
    fft_total = define_values(VECTORS / "fft_wl16.vh")["FFT_TOTAL_CYCLES"]
    assert top["T_FFT"] == d["PSD_NCHUNKS"] * (fft_total + 2)           # 209 cycles per chunk in the top-level phase
    # the whole segment phase = busy engine states + three control cycles per segment
    busy = sum(want[k] for k in want if k not in ("S_CLR", "S_END"))
    assert top["T_SEG"] == busy + 3 * d["PSD_S"]


def test_bispectrum_state_occupancy():
    d = define_values(VECTORS / "bsp_wl12.vh")
    occ = {r["state"]: int(r["cycles"]) for r in _rows("rtl_state_occupancy.csv") if r["testbench"] == "bsp_engine"}
    logb = d["BSP_MAX_BIN"].bit_length() - 1
    want = rtlmodel.bispectrum_engine_state_cycles(d["BSP_S"], logb)
    assert {k: occ[k] for k in want} == want
    top = {r["state"]: int(r["cycles"]) for r in _rows("rtl_state_occupancy.csv") if r["testbench"] == "bsp_top"}
    busy = sum(want[k] for k in want if k not in ("S_CLR", "S_END"))
    assert top["B_ENG"] == busy + 3 * d["BSP_S"] == 5219 * d["BSP_S"]
    assert top["B_FFT"] == (define_values(VECTORS / "fft_wl16.vh")["FFT_TOTAL_CYCLES"] + 2) * d["BSP_S"]


def test_squarer_wait_histogram_counts_every_bin():
    d = define_values(VECTORS / "psd_wl12.vh")
    assert sum(int(r["bins"]) for r in _rows("rtl_squarer_wait_hist.csv")) == 1024 * d["PSD_S"]


def test_squarer_wait_histogram_equals_model():
    """The residence histogram measured in the simulation is the one the bit-accurate model predicts, bin by bin."""
    from btpsd import vectors as V
    from btpsd.estimator import psd_fixed

    p = V.psd_vectors(12)
    res = psd_fixed(p["x"], p["cfg"], quantize=False)
    model = rtlmodel.squarer_wait_histogram(res.X_re, res.X_im, p["cfg"].squarer_units)
    simulated = {int(r["log2_bin"]): int(r["bins"]) for r in _rows("rtl_squarer_wait_hist.csv")}
    assert model == simulated
