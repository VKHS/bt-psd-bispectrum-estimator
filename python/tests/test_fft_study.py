"""The FFT study (btpsd.fftstudy) behind the sheets fig24-fig28: snapshots, counts, noise model and the committed accuracy table."""
import csv
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]

from btpsd import fftstudy as fs  # noqa: E402
from btpsd.cordic import PAPER_TOL, default_rom  # noqa: E402
from btpsd.r8fft import r8fft_fixed  # noqa: E402
from btpsd.schedule import N, bank_of, output_index, pass_addresses, pass_twiddle_q, row_of, split  # noqa: E402


@pytest.fixture(scope="module")
def rom():
    return default_rom()


# ---------------------------------------------------------------------------------------------------------------
# snapshots
# ---------------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("wl,rounding", [(10, "round"), (12, "round"), (12, "trunc"), (16, "round")])
def test_fixed_snapshots_end_in_the_model_fft(rom, wl, rounding):
    re, im = fs.fft_inputs("six tones", wl, seed=2)
    snaps = fs.pass_snapshots_fixed(re, im, rom, wl, rounding)
    assert len(snaps) == 4
    assert np.array_equal(snaps[0], re + 1j * im)  # after load: the input itself
    yr, yi = r8fft_fixed(re, im, rom, wl, rounding)
    assert np.array_equal(fs.natural_order(snaps[-1]), yr + 1j * yi)


def test_float_snapshots_end_in_the_scaled_dft():
    rng = np.random.default_rng(5)
    x = rng.standard_normal(N) + 1j * rng.standard_normal(N)
    snaps = fs.pass_snapshots_float(x)
    assert np.allclose(fs.natural_order(snaps[-1]), np.fft.fft(x) / N, atol=1e-12)
    # the ROM variant has no rounding: it differs from the exact one only by the coefficient error
    rom_snaps = fs.fft_float_rom(x)
    err = np.abs(fs.natural_order(rom_snaps[-1]) - np.fft.fft(x) / N)
    assert 0 < err.max() < 0.02 * np.abs(x).max()


@pytest.mark.parametrize("k", [0, 1, 37, 150, 255, 256, 400, 511])
def test_a_tone_is_localised_digit_by_digit(k):
    """After pass 0 a tone at bin k fills the 64 words with a2 = k0; after pass 1 the 8 words with (a2, a1) = (k0, k1); after pass 2 one word."""
    x = np.exp(2j * np.pi * k * np.arange(N) / N)
    snaps = fs.pass_snapshots_float(x)
    a = np.arange(N)
    a2, a1, a0 = split(a)
    k0, k1, k2 = k & 7, (k >> 3) & 7, (k >> 6) & 7
    for snap, mask in ((snaps[1], a2 == k0), (snaps[2], (a2 == k0) & (a1 == k1)), (snaps[3], (a2 == k0) & (a1 == k1) & (a0 == k2))):
        mag = np.abs(snap)
        assert np.all(mag[~mask] < 1e-12)
        assert np.all(mag[mask] > 1e-6)
    assert np.isclose(np.abs(snaps[3][64 * k0 + 8 * k1 + k2]), 1.0)  # a unit tone: X[k] / 512 = 1 sits in a single word


def test_tone_addresses_agree_with_the_schedule_maps():
    for row in fs.tone_addresses([0, 1, 37, 150, 511]):
        a = row["addr"]
        assert output_index(a) == row["k"]
        assert (row["bank"], row["row"]) == (int(bank_of(a)), int(row_of(a)))
        assert a == 64 * row["k0"] + 8 * row["k1"] + row["k2"]
    assert [r["addr"] for r in fs.tone_addresses([37, 150])] == [64 * 5 + 8 * 4 + 0, 64 * 6 + 8 * 2 + 2]


def test_digit_table_names_every_stage():
    names = [n for n, _ in fs.digit_table()]
    assert names == ["after load", "after pass 0", "after pass 1", "after pass 2"]
    assert fs.digit_table()[-1][1] == ("k0", "k1", "k2")


# ---------------------------------------------------------------------------------------------------------------
# operation counts
# ---------------------------------------------------------------------------------------------------------------
def _nontrivial(w):
    """True for a twiddle other than +-1, +-j (needs a real rotator)."""
    return min(abs(w.real), abs(w.imag)) > 1e-9


def test_rotation_counts_against_an_explicit_enumeration():
    rc = fs.rotation_counts()
    n = np.arange(N)
    # direct DFT: every product W^(nk)
    direct = sum(_nontrivial(w) for w in np.exp(-2j * np.pi * ((n[:, None] * n[None, :]) % N) / N).ravel())
    assert rc["direct"] == direct
    # in-place radix-2 DIF: stage s has 256 butterflies, butterfly j uses W_N^(2^(s-1) (j mod N/2^s))
    r2 = 0
    for s in range(1, 10):
        for j in range(N // 2):
            r2 += _nontrivial(np.exp(-2j * np.pi * ((j % (N >> s)) << (s - 1)) / N))
    assert rc["radix2"] == r2
    # this scheme: output twiddles of passes 0 and 1 and the two W8 rotations per butterfly
    tw = [sum(_nontrivial(np.exp(-1j * np.pi * int(q) / 512)) for q in pass_twiddle_q(p).ravel()) for p in range(3)]
    assert rc["twiddle_per_pass"] == tw
    assert tw[2] == 0
    assert rc["w8"] == 2 * 192
    assert rc["radix8"] == sum(tw) + 384


def test_operation_counts_are_consistent():
    oc = fs.operation_counts()
    assert oc["butterflies"] == 192 and oc["passes"] == 3
    assert oc["words_read"] == oc["words_written"] == 8 * 192
    assert oc["complex_addsub"] == oc["radix2_addsub"] == 4608  # same number of adders as the radix-2 FFT
    assert oc["radix2_words_read"] == 3 * oc["words_read"] and oc["radix2_passes"] == 3 * oc["passes"]


# ---------------------------------------------------------------------------------------------------------------
# noise model
# ---------------------------------------------------------------------------------------------------------------
def test_noise_model_matches_simulation_for_white_records(rom):
    """Mean error power per bin, split into the systematic part (mean over records) and the random part (variance), against the analytic prediction."""
    wl = 12
    errs = []
    for seed in range(120):
        re, im = fs.fft_inputs("white noise", wl, seed)
        got = fs.natural_order(fs.pass_snapshots_fixed(re, im, rom, wl)[-1])
        errs.append(got - np.fft.fft(re + 1j * im) / N)
    e = np.array(errs)
    random_part = float(e.var(axis=0).mean())  # complex variance per bin
    bias_part = float(np.mean(np.abs(e.mean(axis=0)) ** 2))
    assert random_part == pytest.approx(2 * fs.predicted_error_variance(), rel=0.08)
    assert bias_part == pytest.approx(fs.predicted_bias_power(), rel=0.2)  # 120 records: the mean itself still carries noise
    assert float(np.mean(np.abs(e) ** 2)) == pytest.approx(fs.predicted_error_power(), rel=0.08)


def test_predicted_sqnr_is_within_a_decibel_for_noise_like_records(rom):
    for kind in ("white noise", "EEG-like"):
        for wl in (10, 12, 14):
            re, im = fs.fft_inputs(kind, wl, seed=1)
            pred = fs.combined_sqnr(fs.predicted_sqnr(re, im), fs.rom_only_sqnr(re, im, rom))
            assert abs(fs.fft_sqnr(re, im, wl, rom) - pred) < 1.0, (kind, wl)


def test_rom_error_sets_the_ceiling_at_long_word_lengths(rom):
    """The coefficient error does not depend on the word length: at WL = 24 the fixed-point SQNR equals the ROM-only SQNR."""
    for kind in fs.ACCURACY_INPUTS:
        re, im = fs.fft_inputs(kind, 24, seed=1)
        assert fs.fft_sqnr(re, im, 24, rom) == pytest.approx(fs.rom_only_sqnr(re, im, rom), abs=0.05)
        assert 44.0 < fs.rom_only_sqnr(re, im, rom) < 50.0


def test_truncation_costs_decibels_at_short_word_lengths(rom):
    re, im = fs.fft_inputs("white noise", 12, seed=1)
    assert fs.fft_sqnr(re, im, 12, rom, "round") > fs.fft_sqnr(re, im, 12, rom, "trunc") + 1.0


def test_pass_error_power_splits_into_rom_and_arithmetic(rom):
    re, im = fs.fft_inputs("white noise", 14, seed=3)
    pe = fs.pass_error_power(re, im, 14, rom)
    assert pe["rom"][0] == pe["arith"][0] == pe["total"][0] == 0.0  # nothing happened yet after the load
    for t, r, a in zip(pe["total"][1:], pe["rom"][1:], pe["arith"][1:]):
        assert t == pytest.approx(r + a, rel=0.25)  # independent error sources add


def test_error_floor_per_bin_integrates_to_the_total(rom):
    re, im = fs.fft_inputs("six tones", 12, seed=1)
    floor = fs.error_floor_per_bin(re, im, 12, rom)
    assert floor.shape == (N,)
    ref = np.fft.fft(re + 1j * im) / N
    assert 10 * np.log10(np.sum(np.abs(ref) ** 2) / floor.sum()) == pytest.approx(fs.fft_sqnr(re, im, 12, rom), abs=1e-9)


# ---------------------------------------------------------------------------------------------------------------
# twiddle atlas
# ---------------------------------------------------------------------------------------------------------------
def test_rom_error_vectors(rom):
    v = fs.rom_error_vectors(rom)
    assert v["abs"].shape == (1024,)
    trivial = (v["q"] % 256) == 0
    assert np.all(v["abs"][trivial] < 1e-15)  # +-1, +-j are exact
    assert np.max(np.abs(v["angle"])) <= PAPER_TOL  # the rotator meets the tolerance that generated the ROM
    assert np.all(np.abs(v["gain"]) < 2.0**-9)
    assert np.isclose(v["abs"].max(), np.abs(v["err"]).max())


def test_twiddle_exponents_used():
    ex = fs.twiddle_exponents_used()
    assert ex["W8"] == [0, 64, 128, 192, 256, 320, 384, 448]
    assert set(ex["pass 1"]) == {(16 * lo * k % 1024) // 2 for lo in range(8) for k in range(8)}
    assert set(ex["pass 0"]) == {(2 * b * k % 1024) // 2 for b in range(64) for k in range(8)}


# ---------------------------------------------------------------------------------------------------------------
# committed table
# ---------------------------------------------------------------------------------------------------------------
def test_accuracy_csv_matches_the_model():
    with (ROOT / "results" / "model_fft_accuracy.csv").open() as f:
        rows = list(csv.reader(f))
    assert rows[0] == fs.ACCURACY_HEADER
    got = fs.accuracy_rows()
    assert len(rows) - 1 == len(got) == len(fs.ACCURACY_INPUTS) * len(fs.ACCURACY_WLS) * 2
    for line, ref in zip(rows[1:], got):
        assert line[:3] == [str(ref[0]), str(ref[1]), ref[2]]
        for a, b in zip(line[3:], ref[3:]):
            if b == "":
                assert a == ""
            else:
                assert float(a) == pytest.approx(float(b), abs=0.011)
