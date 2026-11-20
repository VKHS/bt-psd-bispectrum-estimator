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
