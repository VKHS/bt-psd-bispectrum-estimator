"""What the FFT sheets (fig24-fig28) draw and print is what the model, the schedule and the Verilog say.

The data behind the sheets is pinned in ``test_fft_study.py``; here the drawing helpers of ``make_fft_figures.py`` are checked:
the isometric index cube tiles its hexagon, the fibres and bank rotations equal the schedule, the lines marked on the
memory snapshots are where the spectral lines are, the Verilog expressions copied onto the sheets compute the maps of
``btpsd.schedule``, and the structure drawn in fig26 is the structure of the sources.
"""
import csv
import re
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python" / "scripts"))

import make_fft_figures as mf  # noqa: E402
from btpsd import fftstudy as fs  # noqa: E402
from btpsd.rtlmodel import Fft512Rtl  # noqa: E402
from btpsd.rtlparse import RTL_DIR  # noqa: E402
from btpsd.schedule import N, bank_of, pass_addresses, pass_twiddle_q, row_of, split  # noqa: E402


# ---------------------------------------------------------------------------------------------------------------
# fig24: index space
# ---------------------------------------------------------------------------------------------------------------
def test_cube_faces_tile_the_hexagon_exactly_once():
    faces = mf.cube_faces(1.0)
    assert len(faces) == 3 * 64
    corners = [mf.iso(0, 0, 8), mf.iso(8, 0, 8), mf.iso(8, 0, 0), mf.iso(8, 8, 0), mf.iso(0, 8, 0), mf.iso(0, 8, 8)]  # outline, one orientation
    rng = np.random.default_rng(3)
    pts = rng.uniform([-8.0, -9.0], [8.0, 9.0], size=(4000, 2))

    def inside(poly, p):
        """Strictly inside a convex polygon given in either orientation."""
        poly = np.asarray(poly)
        e = np.roll(poly, -1, axis=0) - poly
        d = p[:, None, :] - poly[None, :, :]
        cr = e[None, :, 0] * d[:, :, 1] - e[None, :, 1] * d[:, :, 0]
        return np.all(cr > 1e-9, axis=1) | np.all(cr < -1e-9, axis=1)

    in_hex = inside(corners, pts)
    count = np.zeros(len(pts), dtype=int)
    for poly, _ in faces:
        count += inside(poly, pts)
    assert in_hex.sum() > 1500  # the sample is not trivial
    assert np.all(count[in_hex] == 1)  # every interior point lies in exactly one visible cell
    assert np.all(count[~in_hex] == 0)


def test_cube_face_labels_sit_where_the_cell_is():
    """Each visible cell is labelled by its digits ``(a2, a1, a0)``; its centre is the projected centre of that face of the cube."""
    for poly, (a2, a1, a0) in mf.cube_faces(1.0):
        cx, cy = np.mean(poly, axis=0)
        candidates = []
        if a2 == 7:  # top face
            candidates.append(mf.iso(a0 + 0.5, a1 + 0.5, 8))
        if a0 == 7:  # right face
            candidates.append(mf.iso(8, a1 + 0.5, a2 + 0.5))
        if a1 == 7:  # left face
            candidates.append(mf.iso(a0 + 0.5, 8, a2 + 0.5))
        assert any(np.allclose((cx, cy), c) for c in candidates), (a2, a1, a0)


@pytest.mark.parametrize("p", [0, 1, 2])
def test_fibres_are_the_schedule_rotated_into_slot_order(p):
    for b in range(64):
        fib = mf.fibre(p, b)
        banks = [bk for _, bk in fib]
        hi, lo = b >> 3, b & 7
        s = (hi + lo) & 7
        assert banks == [(m + s) % 8 for m in range(8)]  # slot m reads bank (m + s) mod 8
        digits = np.array([d for d, _ in fib])
        varying = [i for i in range(3) if len(set(digits[:, i])) > 1]
        assert varying == [p]  # exactly the digit a2 / a1 / a0 of the pass changes within a butterfly
        assert list(digits[:, p]) == list(range(8))  # in slot order it counts 0 .. 7


def test_digit_flow_text_matches_the_pass_structure():
    """The sheet says: after pass p the digit of weight a2, a1, a0 (p = 0, 1, 2) has been replaced by k0, k1, k2."""
    table = fs.digit_table()
    assert [name for name, _ in table] == ["after load", "after pass 0", "after pass 1", "after pass 2"]
    assert table[0][1] == ("n2", "n1", "n0") and table[3][1] == ("k0", "k1", "k2")
    for p in range(3):
        before, after = table[p][1], table[p + 1][1]
        changed = [i for i in range(3) if before[i] != after[i]]
        assert changed == [p]  # pass p replaces exactly one digit, the one of weight a_(2-p)
        assert before[p].startswith("n") and after[p] == f"k{p}"


def test_mixed_radix_factorisation_identity():
    """W^(nk) = W8^(n2 k0) W^((8 n1 + n0) k0) W8^(n1 k1) W^(8 n0 k1) W8^(n0 k2) for all n and k (the equation of fig24 (f))."""
    w = np.exp(-2j * np.pi / N)
    w8 = w**64
    n = np.arange(N)
    n2, n1, n0 = n >> 6, (n >> 3) & 7, n & 7
    k = np.arange(N)
    k0, k1, k2 = k & 7, (k >> 3) & 7, k >> 6
    lhs = w ** np.outer(n, k)
    rhs = (w8 ** np.outer(n2, k0)) * (w ** np.outer(8 * n1 + n0, k0)) * (w8 ** np.outer(n1, k1)) * (w ** np.outer(8 * n0, k1)) * (w8 ** np.outer(n0, k2))
    assert np.allclose(lhs, rhs, atol=1e-9)


def _verilog_value(expr: str, env: dict[str, int], width: int) -> int:
    """Value of the slice / concatenation / sum expressions of the address maps (``{a[8:6], a[5:3]}``, ``a[8:6] + a[5:3] + a[2:0]``)."""
    rhs = expr.split("=", 1)[1].strip()
    if rhs.startswith("{"):
        value = 0
        for name, hi, lo in re.findall(r"(\w+)\[(\d+):(\d+)\]", rhs):
            w = int(hi) - int(lo) + 1
            value = (value << w) | ((env[name] >> int(lo)) & ((1 << w) - 1))
        return value
    py = re.sub(r"(\w+)\[(\d+):(\d+)\]", lambda m: f"(({env[m.group(1)]} >> {m.group(3)}) & {(1 << (int(m.group(2)) - int(m.group(3)) + 1)) - 1})", rhs)
    return eval(py, {"__builtins__": {}}) & ((1 << width) - 1)  # noqa: S307 - expressions come from our own sources


def test_copied_load_and_unload_maps_compute_the_schedule_maps():
    pats = {p.split()[-1].strip("\\]"): p for p in mf.INDEX_EXPRESSIONS[:4]}
    ld_bank, ld_row, ul_bank, ul_row = (mf.rtl_expr(p) for p in mf.INDEX_EXPRESSIONS[:4])
    assert ld_bank.startswith("ld_bank") and ld_row.startswith("ld_row") and ul_bank.startswith("ul_bank") and ul_row.startswith("ul_row")
    for a in range(N):  # load: natural order address a = (a2, a1, a0)
        assert _verilog_value(ld_bank, {"ld_addr": a}, 3) == bank_of(a)
        assert _verilog_value(ld_row, {"ld_addr": a}, 6) == row_of(a)
    for k in range(N):  # unload: frequency k, stored at a = 64 k0 + 8 k1 + k2
        a = 64 * (k & 7) + 8 * ((k >> 3) & 7) + (k >> 6)
        assert _verilog_value(ul_bank, {"rd_k": k}, 3) == bank_of(a)
        assert _verilog_value(ul_row, {"rd_k": k}, 6) == row_of(a)
    assert len(pats) == 4


def test_copied_rotation_expressions_compute_the_butterfly_rotation():
    ex = {name: mf._expr(p) for name, p in zip(("s_rd", "m_r", "raddr", "s1", "bk", "bf_tq", "s_wr", "m_w", "wrow", "we", "rd_w"), mf.CORE_EXPRESSIONS)}
    assert ex["s_rd"] == "s_rd = iss_b[5:3] + iss_b[2:0]" and ex["s1"] == "s1 = b1[5:3] + b1[2:0]" and ex["s_wr"] == "s_wr = bw[5:3] + bw[2:0]"
    assert ex["m_r"] == "m_r = JB - s_rd" and ex["m_w"] == "m_w = JB - s_wr" and ex["bk"] == "bk = MB + s1"
    for p in range(3):
        for b in range(64):
            sig = fs.butterfly_signals(p, b)
            for name, var in (("s_rd", "iss_b"), ("s1", "b1"), ("s_wr", "bw")):
                assert _verilog_value(ex[name], {var: b}, 3) == sig["s"]
            for j in range(8):  # bank j delivers (and takes back) slot m = j - s; the crossbar maps slot m to bank m + s
                m = (j - sig["s"]) & 7
                assert sig["slots"][m]["bank"] == j
    assert "run_mode ? row_fn(iss_p, m_r, iss_b[5:3], iss_b[2:0]) : ul_row" in ex["raddr"]
    assert "row_fn(pw, m_w, bw[5:3], bw[2:0])" in ex["wrow"]
    assert "(ld_en && (ld_bank == JB))" in ex["we"] and "rq_all[ul_bank_q*DW +: DW]" in ex["rd_w"]


def test_every_copied_pattern_finds_its_line_in_the_source():
    text = " ".join((RTL_DIR / "fft512_mem.v").read_text().split())
    for pat in mf.INDEX_EXPRESSIONS + mf.CORE_EXPRESSIONS:
        line = mf.rtl_expr(pat)
        assert line and " ".join(line.split()) in text, pat


# ---------------------------------------------------------------------------------------------------------------
# fig25: memory snapshots
# ---------------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("k", [37, 475, 150, 362, 0, 511])
def test_marked_cells_are_exactly_where_the_tone_sits(k):
    """The outlines of fig25 (``line_cells``) enclose the words that hold the spectral line after each pass, and only those."""
    x = np.exp(2j * np.pi * k * np.arange(N) / N)
    snaps = fs.pass_snapshots_float(x)
    for stage in (1, 2, 3):
        grid = np.abs(snaps[stage]).reshape(8, 64)  # row a2, column 8 a1 + a0
        marked = np.zeros((8, 64), dtype=bool)
        row, col, width = mf.line_cells(k, stage)
        marked[row, col:col + width] = True
        assert np.array_equal(grid > 1e-9, marked), (k, stage)
        assert width == {1: 64, 2: 8, 3: 1}[stage]
    row, col, _ = mf.line_cells(k, 3)
    assert 64 * row + col == fs.tone_addresses([k])[0]["addr"]


def test_two_tone_record_peaks_sit_in_the_marked_cells():
    wl = 12
    rom = mf.default_rom()
    re, im = fs.two_tone_record(wl, seed=1)
    snaps = fs.pass_snapshots_fixed(re, im, rom, wl)
    grid = np.abs(snaps[3]).reshape(8, 64)
    for k, _, _ in mf.LINE_BINS:
        row, col, _ = mf.line_cells(k, 3)
        assert grid[row, col] > 10 * np.median(grid)  # a spectral line stands far above the noise floor of the memory


# ---------------------------------------------------------------------------------------------------------------
# fig27: twiddles
# ---------------------------------------------------------------------------------------------------------------
def test_used_rom_entries_are_even_and_few():
    """The FFT only ever asks the ROM for even entries: every exponent is a multiple of two units of pi/512 (e = q / 2 is an integer)."""
    used = mf.used_rom_entries()
    assert set(used) == {"pass 0", "pass 1", "W8^i"}
    for name, r in used.items():
        assert np.all(r % 2 == 0), name
        assert np.all((r > 0) & (r <= 128)), name
    assert len(used["pass 0"]) == 64 and len(used["pass 1"]) == 8 and list(used["W8^i"]) == [128]
    assert len(set(used["pass 0"]) | set(used["pass 1"])) == 64


def test_pass_0_never_reaches_the_last_sector_of_the_circle():
    e = fs.twiddle_exponents_used()
    assert max(e["pass 0"]) == 63 * 7 == 441  # e = b k with b <= 63, k <= 7, so the exponents 442 .. 511 are never used
    assert len(e["pass 0"]) == 244 and len(e["pass 1"]) == 26 and len(e["W8"]) == 8
    q0 = pass_twiddle_q(0)
    assert set(int(v) // 2 for v in q0.ravel()) == set(e["pass 0"])  # exponents and rotator indices agree (q = 2 e)


# ---------------------------------------------------------------------------------------------------------------
# fig28: accuracy
# ---------------------------------------------------------------------------------------------------------------
def _csv_rows():
    with open(ROOT / "results" / "model_fft_accuracy.csv") as f:
        rd = csv.reader(f)
        header = next(rd)
        assert header == fs.ACCURACY_HEADER
        rows = []
        for r in rd:
            rows.append([r[0], int(r[1]), r[2], float(r[3]), float(r[4]) if r[4] else None, float(r[5]) if r[5] else None, float(r[6])])
        return rows


def test_worst_overshoot_of_the_noise_model_is_what_the_sheet_states():
    rows = [r for r in _csv_rows() if r[4] is not None]
    over, name, wl = mf.worst_overshoot(rows)
    assert 4.0 < over < 6.0  # the prediction is conservative: measured SQNR beats it by under 6 dB in every row
    assert name in fs.ACCURACY_INPUTS and wl in fs.ACCURACY_WLS
    assert all(r[3] - r[4] <= over + 1e-9 for r in rows if r[2] == "round")


def test_sqnr_rises_one_decibel_per_decibel_of_level():
    """White noise whose level rises by 6 dB gains 6 dB SQNR until the ROM coefficient error takes over."""
    levels = list(range(-66, -23, 6))
    sq = fs.sqnr_vs_level(12, levels)
    d = np.diff(sq) / np.diff(levels)
    assert np.all(d > 0.8) and np.all(d < 1.2), d
    full = fs.sqnr_vs_level(12, [-12])[0]
    assert full < 50  # the coefficient error caps the SQNR far below the 6 dB per bit line


# ---------------------------------------------------------------------------------------------------------------
# fig26: structure and numbers of the FFT core
# ---------------------------------------------------------------------------------------------------------------
def test_butterfly_signals_equal_the_schedule():
    for p in range(3):
        twq_table = pass_twiddle_q(p)
        for b in range(64):
            sig = fs.butterfly_signals(p, b)
            assert sig["hi"] == b >> 3 and sig["lo"] == b & 7 and sig["s"] == ((b >> 3) + (b & 7)) % 8
            addrs = {sl["addr"] for sl in sig["slots"]}
            assert addrs == {int(a) for a in pass_addresses(p)[b]}
            assert sorted(sl["bank"] for sl in sig["slots"]) == list(range(8))
            for sl in sig["slots"]:
                assert sl["bank"] == bank_of(sl["addr"]) and sl["row"] == row_of(sl["addr"])
                assert (sl["a2"], sl["a1"], sl["a0"]) == tuple(int(v) for v in split(sl["addr"]))
                assert sl["tq"] == int(twq_table[b, sl["m"]])
                assert sl["angle_deg"] == pytest.approx(sl["tq"] * 180 / 512)
            assert sig["slots"][0]["tq"] == 0  # slot 0 is never rotated


def test_the_worked_example_is_the_butterfly_of_fig24():
    sig = fs.butterfly_signals(1, 19)
    assert (sig["hi"], sig["lo"], sig["s"]) == (2, 3, 5)
    assert [sl["bank"] for sl in sig["slots"]] == [5, 6, 7, 0, 1, 2, 3, 4]
    assert [sl["tq"] for sl in sig["slots"]] == [48 * m for m in range(8)]  # 16 lo m with lo = 3
    assert [sl["row"] for sl in sig["slots"]] == [16 + m for m in range(8)]  # {hi, m}
    with pytest.raises(ValueError):
        fs.butterfly_signals(3, 0)
    with pytest.raises(ValueError):
        fs.butterfly_signals(0, 64)


def test_core_facts_come_from_the_sources():
    f = mf.core_facts()  # asserts every drawn instance, parameter, constant and coefficient against the Verilog
    assert (f["lc"], f["lb"], f["lat"], f["stall"]) == (4, 10, 12, 4)
    assert f["lb"] == 1 + f["lc"] + 1 + f["lc"]
    assert (f["n_cru"], f["n_dl"]) == (10, 5)
    assert f["core"]["u_ram"].count == 8 and f["bfly"]["u_tw"].count == 7
    assert f["q"] == {1: 128, 2: 256, 3: 384}  # W8^1, W8^2, W8^3 in units of pi/512


def test_cycle_budget_matches_simulation_and_model():
    b = fs.fft_cycle_budget()
    assert b["issue_slots"] == 3 * 64 + b["stall"] == 196
    assert b["total_cycles"] == b["issue_slots"] + b["lat"] - 1 == 207  # README / ARCHITECTURE / tb_fft512_mem
    assert b["trace_done_cycle"] == b["trace_first_issue"] + b["total_cycles"] == 208
    assert mf.sim_done_cycle() == b["trace_done_cycle"]  # RTL simulation (CI) stored in results/rtl_trace_fft512.csv
    assert Fft512Rtl(wl=12).run() == b["total_cycles"]  # the cycle model
    from btpsd.schedule import schedule_cycles

    assert schedule_cycles(b["lat"])["total_cycles"] == b["total_cycles"]
    with open(ROOT / "results" / "rtl_trace_fft512.csv") as f:
        rows = list(csv.DictReader(f))
    first = next(int(r["cycle"]) for r in rows if r["iss_v"] == "1")
    assert first == b["trace_first_issue"]  # the first issue is sampled one cycle after start


def test_bank_geometry_numbers_on_the_sheet():
    f = mf.core_facts()
    ram = f["core"]["u_ram"]
    assert ram.params["AW"] == "6" and 2 ** int(ram.params["AW"]) == 64  # 64 words per bank
    assert ram.count * 64 * 2 == 1024  # 8 banks x 64 words x 2 WL bit = 1024 WL bit of storage
    oc = fs.operation_counts()
    assert oc["complex_addsub"] == 24 * oc["butterflies"] == 4608  # 24 complex add/subtract per radix-8 butterfly, 192 butterflies
    assert oc["words_read"] == oc["words_written"] == 8 * oc["butterflies"]
