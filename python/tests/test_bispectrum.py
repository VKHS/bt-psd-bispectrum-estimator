import numpy as np

from btpsd import bispectrum as B
from btpsd import signals as sg
from btpsd.metrics import snr_complex_db

NB = 512 * 100


def test_direct_bispectrum_symmetry():
    x = sg.pink(NB, 1)
    Bd = B.bispectrum_direct(x, L=512, hop=256, nfft=512, max_bin=64)
    np.testing.assert_allclose(Bd, Bd.T, atol=1e-12 * np.abs(Bd).max())


def test_bicoherence_detects_phase_coupling():
    k1, k2 = 32, 104
    xc, _ = sg.qpc(NB, coupled=True, seed=1)
    xu, _ = sg.qpc(NB, coupled=False, seed=1)
    bc = B.bicoherence(xc, L=512, hop=512, nfft=512, max_bin=128)
    bu = B.bicoherence(xu, L=512, hop=512, nfft=512, max_bin=128)
    assert bc[k1, k2] > 0.95 and bu[k1, k2] < 0.1
    assert 0.0 <= bc.min() and bc.max() <= 1.0 + 1e-9
    assert np.median(bc) < 0.05


def test_indirect_bispectrum_peak_location():
    x, _ = sg.qpc(512 * 40, coupled=True, seed=2)
    A = np.abs(B.bispectrum_indirect(x, maxlag=48, L=256, nfft=256, max_bin=64))
    A[:3, :] = 0
    A[:, :3] = 0
    i, j = np.unravel_index(np.argmax(A), A.shape)
    assert min(abs(i - 16) + abs(j - 52), abs(i - 52) + abs(j - 16)) <= 2


def test_fixed_bispectrum_vs_direct():
    x, _ = sg.qpc(512 * 40, coupled=True, seed=3)
    ref = B.bispectrum_direct(x, L=512, hop=256, nfft=512, max_bin=128)
    snr = {}
    for wl in (12, 16):
        re, im, S, sh = B.bispectrum_fixed(x, B.BispecConfig(wl=wl))
        snr[wl] = snr_complex_db(ref, B.bispectrum_hw_to_float(re, im, wl, S, sh))
    assert snr[16] > 40 and snr[16] > snr[12] + 10


def test_eeg_coupled_peak_at_alpha_harmonic():
    x = sg.eeg_coupled(512 * 80, seed=1)
    for win in ("rect", "hann"):
        Bd = np.abs(B.bispectrum_direct(x, L=512, hop=256, nfft=512, max_bin=64, window=win))
        Bd[:3, :] = 0
        Bd[:, :3] = 0
        i, j = np.unravel_index(np.argmax(Bd), Bd.shape)
        assert (i, j) == (20, 20), (win, i, j)


def test_third_order_cumulant_shape_and_symmetry():
    x, _ = sg.qpc(512 * 20, coupled=True, seed=4)
    lags, c3 = B.third_order_cumulant(x, maxlag=16, L=256)
    assert c3.shape == (len(lags), len(lags)) == (33, 33)
    assert lags[0] == -16 and lags[-1] == 16
    np.testing.assert_allclose(c3, c3.T, atol=1e-9 * np.abs(c3).max())
    xg = sg.white(512 * 20, seed=4)
    _, cg = B.third_order_cumulant(xg, maxlag=16, L=256)
    assert np.abs(cg).max() < 0.1 * np.abs(c3).max()
