import numpy as np
import pytest

from btpsd.windows import (RAISED_COSINE, TABLE4_SIGNIFICANT_TAPS, apply_window_freq, approx_spectral_kernel,
                           kernel_exactness, spectral_kernel, window, window_metrics)

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


@pytest.mark.parametrize("name", ["rect", "hann", "hamming", "blackman"])
def test_window_metrics_closed_forms(name):
    """Coherent gain, noise bandwidth and main-lobe width follow from the three cosine coefficients."""
    a0, a1, a2 = RAISED_COSINE[name]
    m = window_metrics(name)
    assert m["cg"] == pytest.approx(a0, abs=1e-12)
    assert m["enbw"] == pytest.approx((a0**2 + (a1**2 + a2**2) / 2) / a0**2, abs=1e-9)
    assert m["null"] == pytest.approx({"rect": 1, "hann": 2, "hamming": 2, "blackman": 3}[name], abs=1e-3)


@pytest.mark.parametrize("name,psl", [("rect", -13.26), ("hann", -31.47), ("hamming", -42.67), ("blackman", -58.11)])
def test_peak_sidelobe_levels(name, psl):
    """Peak sidelobe levels of the four windows (the values of the standard window tables, to 0.05 dB)."""
    assert window_metrics(name)["psl_db"] == pytest.approx(psl, abs=0.05)


@pytest.mark.parametrize("name", ["hann", "hamming", "blackman"])
def test_kernel_exactness_periodic_against_symmetric(name):
    assert kernel_exactness(name, periodic=True) < 1e-12
    assert kernel_exactness(name, periodic=False) > 1e-4
