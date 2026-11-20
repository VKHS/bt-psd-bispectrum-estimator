import numpy as np
import pytest

from btpsd import schedule as S
from btpsd.cordic import default_rom
from btpsd.r8fft import (dft8_float, dft8_fixed, fft512_reference_scaled, r8fft_fixed, r8fft_float,
                         to_fixed)
from btpsd.metrics import snr_complex_db


def test_bank_mapping_is_conflict_free_and_bijective():
    assert S.check_conflict_free()
    a = np.arange(512)
    pairs = set(zip(S.bank_of(a).tolist(), S.row_of(a).tolist()))
    assert len(pairs) == 512
    np.testing.assert_array_equal(S.addr_of(S.bank_of(a), S.row_of(a)), a)


def test_pass_tables():
    assert S.pass_addresses(0)[1].tolist() == [1, 65, 129, 193, 257, 321, 385, 449]  # vary a2
    assert S.pass_addresses(2)[1].tolist() == list(range(8, 16))  # vary a0, (a2,a1) = (0,1)
    assert S.pass_twiddle_q(0)[1].tolist() == [0, 2, 4, 6, 8, 10, 12, 14]
    assert S.pass_twiddle_q(1)[9].tolist() == [0, 16, 32, 48, 64, 80, 96, 112]
    assert not S.pass_twiddle_q(2).any()
    assert S.pass_twiddle_q(0).max() == 2 * 63 * 7


def test_dft8_float_and_fixed():
    rng = np.random.default_rng(0)
    x = rng.standard_normal((5, 8)) + 1j * rng.standard_normal((5, 8))
    np.testing.assert_allclose(dft8_float(x), np.fft.fft(x, axis=-1), atol=1e-12)
    rom = default_rom()
    xr = rng.integers(-20000, 20000, (64, 8))
    xi = rng.integers(-20000, 20000, (64, 8))
    yr, yi = dft8_fixed(xr, xi, rom, 16)
    ref = np.fft.fft(xr + 1j * xi, axis=-1) / 8
    assert snr_complex_db(ref, yr + 1j * yi) > 45


def test_memory_fft_and_r8_float_equal_numpy():
    rng = np.random.default_rng(1)
    x = rng.standard_normal(512) + 1j * rng.standard_normal(512)
    np.testing.assert_allclose(r8fft_float(x), np.fft.fft(x), atol=1e-9)
    m = S.MemoryFFT()
    m.load(x)
    np.testing.assert_allclose(m.run(), np.fft.fft(x), atol=1e-9)
    assert len(m.trace) == 192
    for _, _, banks, _ in m.trace:
        assert sorted(banks) == list(range(8))


def test_fixed_fft_accuracy():
    rng = np.random.default_rng(2)
    wl = 16
    for tol, floor in ((0.01, 40.0), (0.0025, 48.0)):
        x = 0.5 * (rng.uniform(-1, 1, (4, 512)) + 1j * rng.uniform(-1, 1, (4, 512)))
        re, im = to_fixed(x, wl)
        fr, fi = r8fft_fixed(re, im, default_rom(tol=tol), wl)
        assert snr_complex_db(fft512_reference_scaled(re, im), fr + 1j * fi) > floor


def test_fixed_fft_no_overflow_for_full_scale_impulse_and_dc():
    rom = default_rom()
    wl = 16
    re = np.zeros(512, dtype=np.int64)
    im = np.zeros(512, dtype=np.int64)
    re[0] = 32767
    fr, fi = r8fft_fixed(re, im, rom, wl)
    assert np.max(np.abs(fr - 32767 / 512)) <= 3 and np.max(np.abs(fi)) <= 3
    fr, fi = r8fft_fixed(np.full(512, 32767), im, rom, wl)  # DC -> single bin (scale 1/512 * 512)
    assert abs(fr[0] - 32767) <= 2 and np.max(np.abs(fr[1:])) <= 3


def test_schedule_cycles_model():
    assert S.schedule_cycles(8)["stalls_per_pass"] == [0, 0, 0]
    r = S.schedule_cycles(12)
    assert r["stalls_per_pass"] == [0, 4, 0] and r["issue_cycles"] == 192 + 4
    assert r["total_cycles"] == r["issue_cycles"] + 12 - 1


def test_butterfly_clamps_the_extreme_difference():
    """(32767) - (-32768) halves to +32768 with round-half-up; the hardware clamps it to 32767, and so does the model."""
    from btpsd.cordic import default_rom
    from btpsd.r8fft import dft8_fixed

    wl = 16
    hi, lo = (1 << (wl - 1)) - 1, -(1 << (wl - 1))
    re = np.array([hi, lo, hi, lo, hi, lo, hi, lo], dtype=np.int64)
    im = np.array([lo, hi, lo, hi, lo, hi, lo, hi], dtype=np.int64)
    yr, yi = dft8_fixed(re, im, default_rom(), wl)
    assert yr.max() <= hi and yi.max() <= hi and yr.min() >= lo and yi.min() >= lo
    assert yr[4] == hi  # (A - B) / 2 would be +32768 without the clamp
