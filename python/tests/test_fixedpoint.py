import numpy as np
import pytest

from btpsd.fixedpoint import (QFormat, best_csd_int, csd_from_real, csd_value, naf, naf_weight,
                              naf_weight_array, rshift, saturate, wrap)


def test_qformat_ranges_and_roundtrip():
    q = QFormat(16, 15)
    assert (q.min_int, q.max_int) == (-32768, 32767)
    x = np.array([-1.0, -0.5, 0.0, 0.25, 0.999])
    c = q.quantize(x)
    assert c[0] == -32768 and c[2] == 0 and c[3] == 8192
    assert np.max(np.abs(q.to_float(c) - x)) <= 2.0**-16 + 1e-12
    assert q.quantize(2.0) == 32767 and q.quantize(-2.0) == -32768  # saturation
    assert q.quantize(1.0, sat=False) == -32768  # wrap-around


def test_rshift_modes():
    x = np.array([-7, -6, -5, -1, 0, 1, 5, 6, 7])
    np.testing.assert_array_equal(rshift(x, 1, "trunc"), x >> 1)
    np.testing.assert_array_equal(rshift(x, 1, "round"), np.floor(x / 2 + 0.5).astype(int))
    np.testing.assert_array_equal(rshift(x, 0, "round"), x)
    with pytest.raises(ValueError):
        rshift(x, -1)
    with pytest.raises(ValueError):
        rshift(x, 1, "banker")


def test_saturate_wrap():
    x = np.array([-40000, -32768, 0, 32767, 40000])
    np.testing.assert_array_equal(saturate(x, 16), [-32768, -32768, 0, 32767, 32767])
    np.testing.assert_array_equal(wrap(np.array([32768, 65535, -32769]), 16), [-32768, -1, 32767])


def test_naf_identity_nonadjacent_and_weight():
    for n in range(-2000, 2001):
        d = naf(n)
        assert csd_value(d) == n
        shifts = [e for _, e in d]
        assert all(b - a >= 2 for a, b in zip(shifts, shifts[1:]))  # non-adjacent
        assert naf_weight(n) == len(d)
    assert naf_weight(0) == 0 and naf_weight(7) == 2 and naf_weight(127) == 2
    np.testing.assert_array_equal(naf_weight_array([1, 2, 3, 7, 11]), [1, 1, 2, 2, 3])


def test_best_csd_int_is_optimal():
    rng = np.random.default_rng(1)
    pool = {k: [c for c in range(0, 6000) if naf_weight(c) <= k] for k in (1, 2, 3)}
    for _ in range(100):
        t = float(rng.uniform(0, 3000))
        k = int(rng.integers(1, 4))
        m = best_csd_int(t, k)
        ref = min(pool[k], key=lambda c: (abs(c - t), c))
        assert naf_weight(m) <= k and abs(m - t) == abs(ref - t)
    assert best_csd_int(-11585.2, 1) == -8192


def test_csd_from_real():
    m, digits = csd_from_real(0.6366, 10, 3)
    assert csd_value(digits) == m and len(digits) <= 3
    assert abs(m / 1024 - 0.6366) < 0.01
