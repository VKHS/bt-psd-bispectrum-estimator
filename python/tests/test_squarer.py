import numpy as np

from btpsd.squarer import mag_sq, spectrum_power_cycles, square_odds, square_serial


def test_sum_of_odds_equals_square():
    for n in range(-60, 61):
        assert square_odds(n) == n * n
    v, c = square_serial(np.array([-7, 0, 5, 32767]))
    assert v.tolist() == [49, 0, 25, 32767**2] and c.tolist() == [7, 0, 5, 32767]


def test_mag_sq_cycles():
    v, c = mag_sq(np.array([3, -4]), np.array([4, 3]), units=2)
    assert v.tolist() == [25, 25] and c.tolist() == [4, 4]
    _, c1 = mag_sq(np.array([3]), np.array([4]), units=1)
    assert c1.tolist() == [7]
    d = spectrum_power_cycles(np.array([1, 2, 3]), np.array([3, 2, 1]))
    assert d["total"] == 3 + 2 + 3 and d["max_per_bin"] == 3
