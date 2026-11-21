"""The RTL structure that the figures and the cycle formulas rest on.

``btpsd.rtlparse`` reads the Verilog sources; these tests pin down what it must find there, check that the state-machine
drawings contain exactly the transitions of the sources, and tie the pipeline constants and the per-state cycle formulas
of the model to the RTL and to the numbers in ``sim/vectors``.
"""
import re
import sys
from pathlib import Path

import pytest

from btpsd import rtlmodel
from btpsd.rtlparse import RTL_DIR, define_values, fsm_edges, parse_all, reg_widths, state_values

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "python" / "scripts"


def _code(name: str) -> str:
    """Source of ``rtl/<name>`` without comments."""
    src = (RTL_DIR / name).read_text()
    src = re.sub(r"/\*.*?\*/", "", src, flags=re.S)
    return re.sub(r"//[^\n]*", "", src)


@pytest.fixture(scope="module")
def mods():
    return parse_all()


def test_module_inventory(mods):
    assert sorted(mods) == ["bispectrum_engine", "bt_bispectrum_top", "bt_psd_top", "cru_rom", "cru_rotate", "delay_line",
                            "fft512_mem", "mag_sq_odd", "psd_engine", "r8_butterfly", "ram_1r1w"]
    assert [m.name for m in mods.values() if m.generated] == ["cru_rom"]


def _count(mods, parent, child):
    return sum(i.count for i in mods[parent].instances if i.module == child)


def test_instance_counts(mods):
    assert _count(mods, "cru_rotate", "cru_rom") == 1
    # three rotators for W8^i on the difference branch, seven output twiddle rotators (generate loop k = 1..7)
    assert _count(mods, "r8_butterfly", "cru_rotate") == 3 + 7
    assert _count(mods, "r8_butterfly", "delay_line") == 3
    # eight banks (generate loop), one butterfly, two control delay lines
    assert _count(mods, "fft512_mem", "ram_1r1w") == 8
    assert _count(mods, "fft512_mem", "r8_butterfly") == 1
    assert _count(mods, "fft512_mem", "delay_line") == 2
    assert _count(mods, "psd_engine", "mag_sq_odd") == 1
    assert _count(mods, "bt_psd_top", "fft512_mem") == 2           # ping-pong pair
    assert _count(mods, "bt_psd_top", "psd_engine") == 1
    assert _count(mods, "bt_bispectrum_top", "fft512_mem") == 1
    assert _count(mods, "bt_bispectrum_top", "bispectrum_engine") == 1
    assert _count(mods, "bt_bispectrum_top", "ram_1r1w") == 1       # half-frame buffer


def test_state_machines(mods):
    names = {
        "bt_psd_top": ["T_INIT", "T_CLR", "T_LOAD", "T_FFT", "T_SEG"],
        "bt_bispectrum_top": ["B_INIT", "B_CLR", "B_COPY", "B_EXT", "B_FFT", "B_ENG"],
        "fft512_mem": ["ST_IDLE", "ST_RUN", "ST_DRAIN"],
        "psd_engine": list(rtlmodel.PSD_STATES),
        "bispectrum_engine": list(rtlmodel.BSP_STATES),
    }
    for mod, expect in names.items():
        states = mods[mod].states
        assert sorted(states, key=states.get) == expect, mod
        assert sorted(states.values()) == list(range(len(expect))), mod   # encodings 0 .. n-1
    for mod in ("mag_sq_odd", "ram_1r1w", "delay_line", "cru_rotate", "r8_butterfly", "cru_rom"):
        assert mods[mod].states == {}


def test_transition_counts():
    assert len(fsm_edges(RTL_DIR / "psd_engine.v")) == 19
    assert len(fsm_edges(RTL_DIR / "bispectrum_engine.v")) == 19
    assert len(fsm_edges(RTL_DIR / "bt_psd_top.v")) == 6
    assert len(fsm_edges(RTL_DIR / "bt_bispectrum_top.v")) == 7
    assert len(fsm_edges(RTL_DIR / "fft512_mem.v")) == 3


def test_drawn_state_machines_equal_the_sources(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(SCRIPTS))
    import make_rtl_figures as mrf

    monkeypatch.setattr(mrf, "FIGS", tmp_path)
    mrf.DRAWN.clear()
    for fn in mrf.FIGURES:
        fn()
    assert sorted(f.stem for f in tmp_path.glob("*.png")) == sorted(f.stem for f in tmp_path.glob("*.svg"))
    assert {k[0] for k in mrf.DRAWN} == {"bt_psd_top.v", "bt_bispectrum_top.v", "fft512_mem.v", "psd_engine.v",
                                         "bispectrum_engine.v"}
    for (rtl_file, var), drawn in mrf.DRAWN.items():
        assert drawn == set(fsm_edges(RTL_DIR / rtl_file, var)), rtl_file


def test_loop_bounds_in_the_engines():
    psd = _code("psd_engine.v")
    for needle in ("t == 3'd5",          # merge: six taps
                   "m == 9'd511",        # 512 merged pairs
                   "t == 3'd4",          # window: five taps
                   "k == 10'd1023",      # 1024 bins
                   "cnt == 10'd1023"):   # clear
        assert needle in psd, needle
    bsp = _code("bispectrum_engine.v")
    assert "t == 3'd4" in bsp and "k == NBINS - 1" in bsp
    assert "f2 == MAXB - 1" in bsp and "f1 == MAXB - 1" in bsp
    assert "NBINS = 2 * MAXB" in bsp and "MAXB  = 1 << LOGB" in bsp
    assert "cnt == {2*LOGB{1'b1}}" in bsp


def test_state_encodings_used_by_the_testbenches(mods):
    # the testbenches count cycles per state with literal encodings
    tb = (ROOT / "tb" / "tb_psd_top.v").read_text()
    assert f"dut.st == 3'd{mods['bt_psd_top'].states['T_FFT']}" in tb
    assert f"dut.st == 3'd{mods['bt_psd_top'].states['T_SEG']}" in tb
    assert f"dut.u_eng.st == 4'd{mods['psd_engine'].states['W_S']}" in tb
    tb = (ROOT / "tb" / "tb_bispectrum_top.v").read_text()
    assert f"dut.st == 3'd{mods['bt_bispectrum_top'].states['B_FFT']}" in tb
    assert f"dut.st == 3'd{mods['bt_bispectrum_top'].states['B_ENG']}" in tb


def test_pipeline_constants():
    d = define_values(RTL_DIR / "btpsd_params.vh")
    assert d["LAT_CRU"] == rtlmodel.LAT_CRU == 4
    assert d["LAT_BFLY"] == rtlmodel.LAT_BFLY == 10 == 1 + d["LAT_CRU"] + 1 + d["LAT_CRU"]
    assert rtlmodel.LAT == rtlmodel.LAT_BFLY + 2 == 12
    assert rtlmodel.stall_cycles(rtlmodel.LAT) == rtlmodel.LAT - 8 == 4
    bf = _code("r8_butterfly.v")
    assert "delay_line #(.W(10*WL), .N(LC))" in bf and "delay_line #(.W(70), .N(LC + 2))" in bf


def test_cru_register_stages():
    w = reg_widths(RTL_DIR / "cru_rotate.v", {"WL": 16})
    stages = {
        "R1": ["x1r", "x1i", "a1", "b1", "inv1", "quad1", "ccw1"],
        "R2": ["y2r", "y2i", "inv2", "quad2"],
        "R3": ["z3r", "z3i", "quad3"],
        "R4": ["y4r", "y4i"],
    }
    bits = {k: sum(w[n] for n in names) for k, names in stages.items()}
    assert bits == {"R1": 66, "R2": 69, "R3": 58, "R4": 32}
    assert len(stages) == rtlmodel.LAT_CRU


def test_fft_cycle_totals_match_the_vectors():
    d = define_values(ROOT / "sim" / "vectors" / "fft_wl16.vh")
    stall = rtlmodel.stall_cycles(rtlmodel.LAT)
    slots = 3 * 64 + stall
    assert d["FFT_STALL_CYCLES"] == stall
    assert d["FFT_TOTAL_CYCLES"] == slots - 1 + rtlmodel.LAT == 207


def test_psd_state_cycles_against_the_vectors():
    d = define_values(ROOT / "sim" / "vectors" / "psd_wl12.vh")
    c = rtlmodel.psd_engine_state_cycles(d["PSD_S"], d["PSD_SQUARER_CYCLES"])
    assert c["M_A"] + c["M_B"] + c["M_C"] + c["M_F"] == 3 * 9728
    # per bin: 5 taps x 3 + W_F + W_G + W_S (n + 3) + W_ACC = 21 + n
    assert c["W_A"] + c["W_B"] + c["W_C"] + c["W_F"] + c["W_G"] + c["W_S"] + c["W_ACC"] == \
        3 * 1024 * 21 + d["PSD_SQUARER_CYCLES"]
    multiplier = rtlmodel.psd_engine_state_cycles(d["PSD_S"], 0, serial=False)
    assert multiplier["W_S"] == 2 * 1024 * 3


def test_bispectrum_state_cycles():
    c = rtlmodel.bispectrum_engine_state_cycles(frames=1, logb=5, clears=0)
    assert sum(v for k, v in c.items() if k != "S_END") == 16 * 64 + 3 * 32 + 4 * 32 * 32 == 5216
    assert c["T_A"] == c["T_B"] == c["T_C"] == c["T_D"] == 32 * 32
