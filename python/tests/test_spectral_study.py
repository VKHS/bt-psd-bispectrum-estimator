"""Tests of the spectral study behind fig29: the definitions, the statements printed on the sheet, and the table behind it."""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

from btpsd import spectralstudy as ss
from btpsd.estimator import HwConfig, psd_exact_welch, psd_fixed

ROOT = Path(__file__).resolve().parents[2]
CSV = ROOT / "results" / "model_spectral_study.csv"
FIG = ROOT / "docs" / "figures" / "fig29_spectral_study"


# ---------------------------------------------------------------------------------------------------------------
# records and definitions
# ---------------------------------------------------------------------------------------------------------------
def test_tone_is_a_cosine_of_the_requested_frequency():
    x = ss.tone(37.0)
    assert len(x) == ss.NCHUNKS * 512
    assert 0.89 < np.abs(x).max() <= ss.TONE_AMP + 1e-12
    assert int(np.argmax(np.abs(np.fft.rfft(x[:1024])))) == 37


def test_two_tones_amplitude_ratio():
    x = ss.two_tones(100.0, 150.0, ratio_db=-20.0)
    X = np.abs(np.fft.rfft(x[:1024]))
    assert X[150] / X[100] == pytest.approx(0.1, rel=1e-3)
    x0 = ss.two_tones(100.0, 150.0, ratio_db=0.0)
    X0 = np.abs(np.fft.rfft(x0[:1024]))
    assert X0[150] / X0[100] == pytest.approx(1.0, rel=1e-3)


def test_noise_record_is_seeded():
    a, b, c = ss.noise_record(5, 40), ss.noise_record(5, 40), ss.noise_record(6, 40)
    assert np.array_equal(a, b) and not np.array_equal(a, c)
    assert a.std() == pytest.approx(ss.AVG_SIGMA, rel=0.03)
    assert np.abs(a).max() < 1.0  # never clips at word length 16 (and at any word length)


def test_db_and_floor_definitions():
    assert ss.db_re(1.0, 1.0) == 0.0
    assert ss.db_re(0.01, 1.0) == pytest.approx(-20.0)
    assert ss.db_re(0.0, 1.0) == ss.FLOOR_DB
    p = np.full(512, 1e-9)
    p[100] = 1.0
    p[105] = 1e-1  # inside the excluded zone: ignored
    p[120] = 1e-3
    p[2] = 0.5  # below the analysis band: ignored
    p[510] = 0.5  # above the analysis band: ignored
    assert ss.leakage_floor_db(p, 100.0, 1.0) == pytest.approx(-30.0)
    assert ss.leakage_floor_db(p, 100.0, 1.0, exclude=3) == pytest.approx(-10.0)


def test_valley_depth_definition():
    p = np.full(512, 1e-6)
    p[100], p[106], p[103] = 1.0, 0.5, 1e-3
    p[101] = p[102] = p[104] = p[105] = 0.2
    assert ss.valley_depth_db(p, 100.0, 106.0) == pytest.approx(10 * np.log10(0.5 / 1e-3))
    assert np.isnan(ss.valley_depth_db(p, 100.0, 101.0))
    q = p.copy()
    q[103] = 0.0
    assert ss.valley_depth_db(q, 100.0, 106.0) == pytest.approx(-ss.FLOOR_DB)  # capped, not infinite
    # the peak value is the maximum over the tone's bin and its two neighbours
    r = np.full(512, 1e-6)
    r[101:107] = 0.5
    r[100], r[107], r[103] = 0.1, 1.0, 0.01
    r[99] = 1.0
    assert ss.valley_depth_db(r, 100.0, 107.0) == pytest.approx(10 * np.log10(1.0 / 0.01))


# ---------------------------------------------------------------------------------------------------------------
# what the sheet says about the merge
# ---------------------------------------------------------------------------------------------------------------
def test_three_estimates_agree_at_the_peak():
    P = ss.psds(ss.tone(ss.TONE_BIN), "hann")
    pk = {k: P[k][:512].max() for k in P}
    assert abs(10 * np.log10(pk["hw"] / pk["ref"])) < 0.6  # the merge removes a little of the peak's energy
    assert abs(10 * np.log10(pk["hw"] / pk["fd"])) < 0.2


def test_floor_between_bins_is_set_by_the_merge_not_by_the_window():
    i = ss.OFFSETS.index(0.5)
    fl = {w: ss.floor_vs_offset(w, [ss.OFFSETS[i]]) for w in ss.WINDOWS}
    hw = {w: float(fl[w]["hw"][0]) for w in ss.WINDOWS}
    ref = {w: float(fl[w]["ref"][0]) for w in ss.WINDOWS}
    fd = {w: float(fl[w]["fd"][0]) for w in ss.WINDOWS}
    for w in ss.WINDOWS:
        assert -32.0 < hw[w] < -24.0, (w, hw[w])  # about -28 dB for every window
        assert abs(hw[w] - fd[w]) < 0.3  # the float six-tap merge has the same floor: no quantisation involved
    assert max(hw.values()) - min(hw.values()) < 3.0
    # the exact average has window-dependent floors: sidelobe decay of the four windows
    assert ref["rect"] > -30.0 and ref["hamming"] < -44.0 and ref["hann"] < -70.0 and ref["blackman"] < -78.0
    assert hw["hann"] - ref["hann"] > 40.0 and hw["blackman"] - ref["blackman"] > 50.0


def test_tone_panel_numbers():
    t = ss.tone_panel()
    assert t["f"] == 300.5 and t["window"] == "hann"
    assert t["ref"][300] == pytest.approx(t["ref"][301], abs=1e-6)  # symmetric about the tone
    assert max(t["ref"][:513]) == pytest.approx(0.0, abs=1e-9)
    assert t["floor"]["hw"] == pytest.approx(t["floor"]["fd"], abs=0.2)
    assert t["floor"]["ref"] < -70.0
    assert -32.0 < t["floor"]["hw"] < -24.0
    # the whole one-sided band is returned
    assert len(t["bins"]) == 513 and t["bins"][0] == 0 and t["bins"][-1] == 512


def test_floor_on_bin_centres():
    off = [0.0, 1.0]
    f = ss.floor_vs_offset("hann", off)
    # even bin: the float merge is exact, so only the hardware (rotator ROM) has a floor
    assert f["fd"][0] <= ss.FLOOR_DB and f["ref"][0] <= ss.FLOOR_DB
    assert -65.0 < f["hw"][0] < -50.0
    # odd bin: the merge interpolates, so even the float model has a floor, and the hardware agrees with it
    assert -40.0 < f["fd"][1] < -30.0 and abs(f["hw"][1] - f["fd"][1]) < 0.3
    assert f["ref"][1] <= ss.FLOOR_DB


def test_even_bin_floor_follows_the_rom_and_the_between_bin_floor_does_not():
    def floors(f, tol):
        x = ss.tone(f)
        peak = psd_exact_welch(x, "hann")[:512].max()
        hw = psd_fixed(x, HwConfig(wl=ss.WL, window="hann", rom_tol=tol)).to_float()
        return ss.leakage_floor_db(hw, f, peak)

    coarse, fine = floors(300.0, 0.01), floors(300.0, 0.005)
    assert fine < coarse - 3.0  # a tighter rotator ROM lowers the floor of an on-bin tone
    assert abs(floors(300.5, 0.01) - floors(300.5, 0.005)) < 0.2  # but not the merge floor


def test_valley_is_limited_by_the_merge():
    sep = [8]
    v = {w: ss.valley_vs_separation(w, sep) for w in ss.WINDOWS}
    for w in ("hann", "hamming", "blackman"):
        assert v[w]["ref"][0] > 40.0  # exact: deep valleys for the tapered windows
        assert v[w]["hw"][0] < 36.0  # hardware: stops at the merge floor
        assert abs(v[w]["hw"][0] - v[w]["fd"][0]) < 2.5  # float six-tap merge sits close to the hardware model
    assert v["rect"]["ref"][0] < 30.0  # the rectangular window itself leaks


# ---------------------------------------------------------------------------------------------------------------
# averaging
# ---------------------------------------------------------------------------------------------------------------
def test_overlap_correlation_and_variance():
    assert ss.welch_overlap_rho("hann") == pytest.approx(1.0 / 6.0, abs=1e-12)
    assert ss.welch_overlap_rho("rect") == pytest.approx(0.5, abs=1e-12)
    assert ss.welch_overlap_variance("hann", 1) == pytest.approx(1.0)
    rho2 = (1.0 / 6.0) ** 2
    assert ss.welch_overlap_variance("hann", 8) == pytest.approx((1 + 2 * (1 - 1 / 8) * rho2) / 8)
    assert ss.welch_overlap_variance("rect", 64) > ss.welch_overlap_variance("hann", 64) > ss.welch_overlap_variance("blackman", 64)
    # at large S the overlap costs a factor 1 + 2 rho^2 in variance
    assert ss.welch_overlap_variance("hann", 1024) * 1024 == pytest.approx(1 + 2 * rho2, rel=1e-3)


def test_averaging_follows_the_analytic_curve():
    a = ss.averaging_study("hann", segments=(1, 4, 16), nreal=16)
    assert list(a["S"]) == [1, 4, 16]
    assert np.allclose(a["inv_sqrt"], [1.0, 0.5, 0.25])
    for key in ("hw", "ref"):
        assert np.allclose(a[key], a["analytic"], rtol=0.12), (key, a[key], a["analytic"])
    assert np.allclose(a["hw"], a["ref"], rtol=0.03)  # the merge error does not change the spread
    assert np.all(a["hw"][:-1] > a["hw"][1:])
    assert np.all(a["bias_db"] < 0.0) and np.all(a["bias_db"] > -0.4)


# ---------------------------------------------------------------------------------------------------------------
# the table behind the sheet
# ---------------------------------------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def table():
    with open(CSV, newline="") as f:
        rows = list(csv.reader(f))
    return rows[0], rows[1:]


def _find(rows, study, window, model, x):
    hit = [float(r[4]) for r in rows if r[0] == study and r[1] == window and r[2] == model and float(r[3]) == float(x)]
    assert len(hit) == 1, (study, window, model, x, hit)
    return hit[0]


def test_table_layout(table):
    header, rows = table
    assert header == ss.CSV_HEADER
    keys = [tuple(r[:4]) for r in rows]
    assert len(keys) == len(set(keys))
    n_tone = 3 * 513
    n_floor = len(ss.WINDOWS) * 3 * len(ss.OFFSETS)
    n_valley = len(ss.VALLEY_RATIOS_DB) * len(ss.WINDOWS) * 3 * len(ss.SEPARATIONS)
    n_avg = 5 * len(ss.AVG_SEGMENTS)
    assert len(rows) == n_tone + n_floor + n_valley + n_avg
    assert {r[0] for r in rows} == {"tone_psd_db_re_exact_peak", "floor_db_re_exact_peak_vs_offset", "valley_db_vs_separation_ratio0db",
                                    "valley_db_vs_separation_ratio-20db", "averaging_rel_std", "averaging_bias_db"}


def test_table_matches_the_model(table):
    _, rows = table
    t = ss.tone_panel()
    assert _find(rows, "tone_psd_db_re_exact_peak", "hann", ss.MODEL_LABELS["hw"], 300) == pytest.approx(float(t["hw"][300]), abs=2e-4)
    assert _find(rows, "tone_psd_db_re_exact_peak", "hann", "exact_welch", 450) == pytest.approx(float(t["ref"][450]), abs=2e-4)
    fl = ss.floor_vs_offset("hann", [0.5, 1.0])
    assert _find(rows, "floor_db_re_exact_peak_vs_offset", "hann", ss.MODEL_LABELS["hw"], 0.5) == pytest.approx(fl["hw"][0], abs=2e-4)
    assert _find(rows, "floor_db_re_exact_peak_vs_offset", "hann", "float_six_tap", 1.0) == pytest.approx(fl["fd"][1], abs=2e-4)
    vl = ss.valley_vs_separation("blackman", [8], ratio_db=-20.0)
    assert _find(rows, "valley_db_vs_separation_ratio-20db", "blackman", "exact_welch", 8) == pytest.approx(vl["ref"][0], abs=2e-4)
    assert _find(rows, "averaging_rel_std", "hann", "analytic_overlap", 8) == pytest.approx(np.sqrt(ss.welch_overlap_variance("hann", 8)), abs=2e-5)
    assert _find(rows, "averaging_rel_std", "hann", "inv_sqrt_S", 16) == pytest.approx(0.25, abs=1e-5)
    avg = ss.averaging_study("hann", segments=(2,), nreal=ss.AVG_NREAL)
    assert _find(rows, "averaging_rel_std", "hann", ss.MODEL_LABELS["hw"], 2) == pytest.approx(float(avg["hw"][0]), abs=2e-5)
    assert _find(rows, "averaging_bias_db", "hann", "hw_bias_db", 2) == pytest.approx(float(avg["bias_db"][0]), abs=2e-5)


def test_sheet_files_exist():
    for ext in (".svg", ".png"):
        assert FIG.with_suffix(ext).exists()
