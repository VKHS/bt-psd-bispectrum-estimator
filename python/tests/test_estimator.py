import numpy as np
import pytest

from btpsd import estimator as E
from btpsd import signals as sg
from btpsd.metrics import sqnr_db
from btpsd.fdfilter import fd_kernel

N = 512 * 32


def test_exact_merge_reproduces_welch():
    x = sg.pink(N, seed=1)
    for win in ("rect", "hann", "hamming", "blackman"):
        P = E.psd_float_fd(x, win, mode="exact")
        ref = E.psd_exact_welch(x, win)
        np.testing.assert_allclose(P, ref, rtol=1e-9, atol=1e-12 * ref.max())


def test_fd6_psd_level_accuracy_reflects_spectral_slope():
    ref_w, ref_b = E.psd_exact_welch(sg.white(N, 2), "hann"), E.psd_exact_welch(np.cumsum(sg.white(N, 2)), "hann")
    s_white = sqnr_db(ref_w, E.psd_float_fd(sg.white(N, 2), "hann"))
    s_brown = sqnr_db(ref_b, E.psd_float_fd(np.cumsum(sg.white(N, 2)), "hann"))
    assert 15 < s_white < 35  # honest limit of the 6-tap merge on noise-like data
    assert s_brown > s_white + 8


def test_more_taps_help_at_psd_level():
    x = sg.pink(N, seed=3)
    ref = E.psd_exact_welch(x, "hann")
    vals = []
    for L in (4, 6, 8, 12):
        h, st = fd_kernel(L, "sinc")
        vals.append(sqnr_db(ref, E.psd_float_fd(x, "hann", h=h, start=st)))
    assert vals == sorted(vals)


def test_fixed_matches_float_fd_at_16_bits():
    x = sg.pink(N, seed=4)
    ref_fd = E.psd_float_fd(x, "hann")
    res = E.psd_fixed(x, E.HwConfig(wl=16))
    assert res.S == 31 and res.avg_shift == 5
    s16 = sqnr_db(ref_fd, res.to_float())
    s20 = sqnr_db(ref_fd, E.psd_fixed(x, E.HwConfig(wl=20)).to_float())
    assert s16 > 45 and s20 > s16


def test_fixed_psd_scale_is_analytic_not_fitted():
    """A bin-centred full-scale-ish sine must land in the right bin with the right power (< 0.2 dB)."""
    n = np.arange(N)
    k0 = 100  # bin of the 1024-point DFT
    x = 0.5 * np.cos(2 * np.pi * k0 * n / 1024)
    res = E.psd_fixed(x, E.HwConfig(wl=20, window="rect"))
    P = res.to_float()
    ref = E.psd_exact_welch(x, "rect")
    assert np.argmax(P[:512]) == k0
    assert abs(10 * np.log10(P[k0] / ref[k0])) < 0.2


def test_config_integer_constants():
    cfg = E.HwConfig(fd_frac=12)
    np.testing.assert_array_equal(cfg.fd_int(), np.round(np.array(cfg.fd_h) * 4096).astype(int))
    assert E.HwConfig(window="hann").win_int() == {-1: -256, 0: 512, 1: -256}
    csd = E.HwConfig(fd_max_terms=2).fd_int()
    from btpsd.fixedpoint import naf_weight
    assert all(naf_weight(int(v)) <= 2 for v in csd)


def test_squarer_cycle_accounting():
    res = E.psd_fixed(sg.white(512 * 4, 5), E.HwConfig(wl=16))
    assert res.squarer_cycles > 0
    assert res.X_re.shape == (3, 1024)


def test_requires_two_chunks():
    with pytest.raises(ValueError):
        E.psd_fixed(np.zeros(512))
