import pytest

from btpsd import paperdata as pd


def test_table3_literal():
    t = pd.table3()
    assert t[6] == 50.08 and t[3] == 8.14 and len(t) == 8
    assert all(t[n] > t[n - 1] for n in (4, 6, 8, 10))  # even taps beat the preceding odd count


def test_table9_internal_consistency():
    c = pd.table9_consistency()
    assert c == {"pdp_mismatch": [], "R_mismatch": []}
    rows = {r["method"]: r for r in pd.table9()}
    p = rows["Proposed Estimator"]
    assert (p["dsp"], p["lut"], p["ff"], p["power_mw"], p["latency_cycles"], p["pdp"]) == (0, 4834, 4420, 191, 1252, 239132)


def test_table8_and_10_shapes():
    t8 = pd.table8()
    assert len(t8) == 9 * 4 * 13
    assert all(r["sqnr_db"] is None for r in t8 if "Unrolled" in r["method"] or "Single" in r["method"])
    prop = [r for r in t8 if r["method"] == "Proposed hardware" and r["signal"] == "EEG-1"]
    assert [r["sqnr_db"] for r in prop] == [20.40, 20.48, 20.46, 20.42]
    t10 = pd.table10()
    assert len(t10) == 9 * 3 * 5
    assert [r["sqnr_db"] for r in t10 if r["signal"] == "EEG-1" and r["case"] == 1] == [11.14, 12.93, 15.08, 17.16, 19.40]


def test_table6_printed_is_verbatim_with_flags():
    t6 = pd.table6_printed()
    assert len(t6) == 128 and t6[0]["a"] == 127 and t6[0]["b"] == 1
    flagged = {r["row"] for r in t6 if r["flag"]}
    assert {8, 18} <= flagged and len(flagged) == 27


def test_table1_and_2_shape_and_pipeline_fill():
    t1, t2 = pd.table1(), pd.table2()
    assert len(t1) == 64 and len(t2) == 64 and "R0" in t1 and "R127" in t2
    assert all(len(v) == 16 for v in t1.values()) and sorted(next(iter(t1.values()))) == list(range(16))
    assert all(len(v) == 15 for v in t2.values()) and sorted(next(iter(t2.values()))) == list(range(8, 23))
    # Table 1: register R(8 g + i) is filled at clock i, so it starts with i dashes
    for g in range(8):
        for i in range(8):
            row = t1[f"R{8 * g + i}"]
            assert [c for c in range(16) if row[c] is None] == list(range(i))
    # Table 2: register R(64 + 8 g + i) starts with g dashes
    for g in range(8):
        for i in range(8):
            row = t2[f"R{64 + 8 * g + i}"]
            assert [c for c in range(8, 23) if row[c] is None] == list(range(8, 8 + g))


def test_table2_row_r113_as_printed_deviates_from_its_neighbours():
    t2 = pd.table2()
    first = lambda r: next(v for v in t2[r].values() if v is not None)
    last = lambda r: [v for v in t2[r].values() if v is not None][-1]
    assert (first("R112"), first("R113"), first("R114")) == (392, 340, 408)   # 8 i - 504 gives 392, 400, 408
    assert (last("R113"), last("R114"), last("R115")) == (1361, 1633, 1665)  # +32 steps give 1601


def test_table4_5_7_literal():
    assert pd.table4() == {"Rectangular": 1, "Triangular": 5, "Parzen": 3, "Hann": 3, "Hamming": 3, "Blackman": 5, "Kaiser": 5}
    t5 = pd.table5()
    assert [r["n"] for r in t5] == list(range(8))
    # the printed angle row matches atan(2**-n) for n = 0..4 and is off by 0.0016 to 0.0029 degrees from n = 5 on
    assert all(abs(r["printed_angle_deg"] - r["angle_deg"]) < 0.01 for r in t5[:5])
    assert all(0.001 < abs(r["printed_angle_deg"] - r["angle_deg"]) < 0.003 for r in t5[5:])
    # the printed tangent row is shifted for n = 1..4 and repeats 1/32
    assert [r["printed_tg"] for r in t5][:6] == [1, 0.25, 0.125, 0.0625, 0.03125, 0.03125]
    assert pd.table7()["EEG-3"] == "A10E" and len(pd.table7()) == 9
