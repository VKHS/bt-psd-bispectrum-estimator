import numpy as np
import pytest
from scipy import signal

from btpsd.reference import (autocorrelation_biased, bt_correlogram_psd, bt_via_periodogram_smoothing,
                             frame_signal, welch_psd, welch_twosided)


def test_frame_signal():
    f = frame_signal(np.arange(10), 4, 2)
    assert f.shape == (4, 4) and list(f[1]) == [2, 3, 4, 5]
    assert frame_signal(np.arange(3), 4, 2).shape == (1, 4)  # zero-padded


def test_bt_equals_smoothed_periodogram():
    """Eq. (1): FT of the lag-windowed correlogram == periodogram circularly smoothed by W."""
    rng = np.random.default_rng(0)
    x = rng.standard_normal(300)
    for lw in ("bartlett", "hann", "hamming"):
        a = bt_correlogram_psd(x, M=40, lag_window=lw, nfft=1024)
        b = bt_via_periodogram_smoothing(x, M=40, lag_window=lw, nfft=1024)
        np.testing.assert_allclose(a, b, atol=1e-10)


def test_autocorrelation_matches_direct():
    rng = np.random.default_rng(1)
    x = rng.standard_normal(200)
    r = autocorrelation_biased(x, 5)
    for k in range(6):
        assert r[k] == pytest.approx(np.sum(x[: 200 - k] * x[k:]) / 200)


def test_welch_matches_scipy_density():
    rng = np.random.default_rng(2)
    x = rng.standard_normal(8192)
    P, K = welch_twosided(x, L=1024, hop=512, window="hann", scaling="density")
    _, ref = signal.welch(x, fs=1.0, window="hann", nperseg=1024, noverlap=512, nfft=1024,
                          detrend=False, return_onesided=False, scaling="density")
    np.testing.assert_allclose(P, ref, rtol=1e-10)
    assert K == 15


def test_welch_psd_one_sided_length():
    f, P = welch_psd(np.random.default_rng(3).standard_normal(4096), fs=256.0, N_fft=1024)
    assert len(f) == len(P) == 513 and f[-1] == pytest.approx(128.0, abs=0.3)
