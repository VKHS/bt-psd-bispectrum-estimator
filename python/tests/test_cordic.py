import numpy as np
import pytest

from btpsd.cordic import (STEP, cordic_gain, cordic_rotate_float, cordic_rotate_int, default_rom,
                          generate_rom, octant_reduce, rotate_fixed, rotate_float, twiddle_q)


def test_octant_reduce_exhaustive():
    q = np.arange(1024)
    quad, r, ccw = octant_reduce(q)
    assert r.min() == 0 and r.max() == 128
    ideal = np.exp(-1j * q * STEP)
    res = np.where(ccw, np.exp(+1j * r * STEP), np.exp(-1j * r * STEP))
    np.testing.assert_allclose((-1j) ** quad * res, ideal, atol=1e-12)


def test_twiddle_q():
    assert list(twiddle_q([0, 1, 2, 255], 512)) == [0, 2, 4, 510]
    assert twiddle_q(3, 8) == 384  # W8^3 = exp(-j 3pi/4)


@pytest.mark.parametrize("tol", [0.01, 0.0025])
def test_rom_meets_tolerances(tol):
    rom = default_rom(tol=tol)
    assert rom.n == 129 and rom.a[0] == 1 and rom.b[0] == 0
    assert np.max(np.abs(rom.angle_error()[1:])) <= tol + 1e-12
    assert np.max(np.abs(rom.gain_error())) <= rom.gain_tol + 1e-12
    assert rom.a[128] == rom.b[128] == 1  # pi/4
    assert np.all(rom.b[1:] >= 1) and np.all(rom.a >= rom.b)


def test_finer_tolerance_costs_more():
    c_coarse = default_rom(tol=0.01).adders().mean()
    c_fine = default_rom(tol=0.0025).adders().mean()
    assert c_fine > c_coarse


def test_rotate_float_with_rom_is_close_to_ideal():
    rom = default_rom()
    x = np.array([1 + 0.5j])
    q = np.arange(1024)
    z = rotate_float(x, q, rom)
    ideal = rotate_float(x, q)
    assert np.max(np.abs(z - ideal)) < 0.012 + 0.002  # angle tol 0.01 rad + gain tol


def test_rotate_fixed_error_bound():
    rom = default_rom()
    rng = np.random.default_rng(0)
    wl = 16
    xr = rng.integers(-20000, 20000, 4096)
    xi = rng.integers(-20000, 20000, 4096)
    q = rng.integers(0, 1024, 4096)
    yr, yi = rotate_fixed(xr, xi, q, rom, wl)
    ref = (xr + 1j * xi) * np.exp(-1j * q * STEP)
    err = np.abs((yr + 1j * yi) - ref)
    mag = np.abs(xr + 1j * xi)
    assert np.all(err <= 0.0102 * mag + 0.0011 * mag + 2.0)  # tol + gain + rounding
    # identity and exact quadrants
    yr, yi = rotate_fixed(np.array([1234]), np.array([-567]), np.array([0]), rom, wl)
    assert (yr[0], yi[0]) == (1234, -567)
    yr, yi = rotate_fixed(np.array([1234]), np.array([-567]), np.array([256]), rom, wl)  # -j
    assert (yr[0], yi[0]) == (-567, -1234)


def test_classical_cordic():
    assert cordic_gain(16) == pytest.approx(1.6467602581, rel=1e-9)
    for ang in np.linspace(-3.1, 3.1, 13):
        z = cordic_rotate_float(1 + 0.5j, ang, 18)
        assert abs(z - (1 + 0.5j) * np.exp(1j * ang)) < 1e-4
    rng = np.random.default_rng(1)
    xr, xi = rng.integers(-20000, 20000, 512), rng.integers(-20000, 20000, 512)
    q = rng.integers(0, 1024, 512)
    yr, yi = cordic_rotate_int(xr, xi, q, 16)
    ref = (xr + 1j * xi) * np.exp(-1j * q * STEP)
    assert np.max(np.abs(yr + 1j * yi - ref)) < 12
