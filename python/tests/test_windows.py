import numpy as np
import pytest

from btpsd.windows import (TABLE4_SIGNIFICANT_TAPS, apply_window_freq, approx_spectral_kernel,
                           spectral_kernel, window)

L = 1024


@pytest.mark.parametrize("name", ["rect", "hann", "hamming", "blackman"])
def test_periodic_raised_cosine_is_exactly_sparse(name):
    w = window(name, L)
    W = np.fft.fft(w) / L
    ker = spectral_kernel(name)
    assert len(ker) == TABLE4_SIGNIFICANT_TAPS[name]
    mask = np.zeros(L, bool)
    for d, c in ker.items():
        mask[d % L] = True
        assert abs(W[d % L] - c) < 1e-14
    assert np.max(np.abs(W[~mask])) < 1e-14


def test_symmetric_window_is_not_exactly_sparse():
    W = np.fft.fft(np.hanning(L)) / L
    assert np.max(np.abs(W[2:-1])) > 1e-6  # leakage outside the 3 central bins


@pytest.mark.parametrize("name", ["rect", "hann", "hamming", "blackman"])
def test_frequency_domain_window_equals_time_domain(name):
    rng = np.random.default_rng(0)
    x = rng.standard_normal(L)
    lhs = apply_window_freq(np.fft.fft(x), spectral_kernel(name))
    rhs = np.fft.fft(x * window(name, L))
    np.testing.assert_allclose(lhs, rhs, atol=1e-10)


def test_known_kernels():
    assert spectral_kernel("hann") == {-1: -0.25, 0: 0.5, 1: -0.25}
    k = spectral_kernel("blackman")
    assert k[0] == pytest.approx(0.42) and k[1] == pytest.approx(-0.25) and k[2] == pytest.approx(0.04)
    assert spectral_kernel("hamming")[1] == pytest.approx(-0.23)


def test_approx_kernel_energy():
    ker, kept = approx_spectral_kernel("triangular", 5)
    assert set(ker) == {-2, -1, 0, 1, 2} and 0.9 < kept <= 1.0
    _, kept_h = approx_spectral_kernel("hann", 3)
    assert kept_h == pytest.approx(1.0)
