import numpy as np
import pytest

from btpsd.fdfilter import (PAPER_FD6, PAPER_FD6_START, exact_half_bin, fd_group_delay, fd_kernel,
                            interp_half_bin, merge_chunk_spectra)


def test_exact_half_bin_matches_definition():
    rng = np.random.default_rng(0)
    N = 16
    a = rng.standard_normal(N) + 1j * rng.standard_normal(N)
    A = np.fft.fft(a)
    n = np.arange(N)
    direct = np.array([np.sum(a * np.exp(-2j * np.pi * (m + 0.5) * n / N)) for m in range(N)])
    np.testing.assert_allclose(exact_half_bin(A), direct, atol=1e-12)


def test_dirichlet_kernel_formula():
    rng = np.random.default_rng(1)
    N = 32
    A = np.fft.fft(rng.standard_normal(N) + 1j * rng.standard_normal(N))
    d = np.arange(N)
    K = 1j * np.exp(-1j * np.pi * (d - 0.5) / N) / (N * np.sin(np.pi * (d - 0.5) / N))
    out = np.array([np.sum(K * A[(m + d) % N]) for m in range(N)])
    np.testing.assert_allclose(out, exact_half_bin(A), atol=1e-12)


def test_merge_exact_equals_1024_point_dft():
    rng = np.random.default_rng(2)
    seg = rng.standard_normal(1024)
    A, B = np.fft.fft(seg[:512]), np.fft.fft(seg[512:])
    X = merge_chunk_spectra(A, B, mode="exact")
    np.testing.assert_allclose(X, np.fft.fft(seg), atol=1e-9)


def test_paper_taps_and_symmetric_kernel():
    h, st = fd_kernel(6, "paper6")
    assert st == PAPER_FD6_START == -2 and np.allclose(h, PAPER_FD6)
    hs, ss = fd_kernel(6, "sinc")
    assert ss == -2 and np.allclose(hs[1:5], PAPER_FD6[1:5], atol=1e-4)
    assert hs[0] == pytest.approx(2 / (5 * np.pi))  # the paper rounds this tap to 1/8
    assert np.allclose(hs, hs[::-1])
    assert abs(PAPER_FD6[0] - hs[0]) < 3e-3


def test_constant_group_delay_for_even_symmetric():
    for L in (4, 6, 8):
        h, _ = fd_kernel(L, "sinc")
        _, gd = fd_group_delay(h, 256)
        assert np.ptp(gd[2:-2]) < 1e-9 and gd[10] == pytest.approx((L - 1) / 2)
    h, _ = fd_kernel(5, "sinc")
    _, gd = fd_group_delay(h, 256)
    assert np.ptp(gd[2:-2]) > 1e-3  # odd length: frequency-dependent delay


def test_fd_error_decreases_with_taps_on_smooth_spectra():
    rng = np.random.default_rng(3)
    from scipy import signal
    x = signal.lfilter([1], [1, -1.8, 0.95], rng.standard_normal(1024 + 200))[200:]
    A, B = np.fft.fft(x[:512]), np.fft.fft(x[512:])
    ref = merge_chunk_spectra(A, B, mode="exact")[1::2]
    errs = []
    for L in (4, 6, 8, 12):
        h, st = fd_kernel(L, "sinc")
        out = interp_half_bin(A - B, h, st)
        errs.append(np.linalg.norm(out - ref) / np.linalg.norm(ref))
    assert errs == sorted(errs, reverse=True)
