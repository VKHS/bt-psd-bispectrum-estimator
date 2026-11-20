"""The full-scale segments of the top-level vectors must keep exercising the saturating paths of the FFT."""

import numpy as np
import pytest

import btpsd.r8fft as R
from btpsd import vectors as V


@pytest.fixture
def clamp_counter(monkeypatch):
    """Count the results of the half add/subtract that the clamp changes (only the rounded extreme difference)."""
    orig = R._addsub_half
    count = {"n": 0}

    def wrapped(ar, ai, br, bi, rounding, wl=None):
        free = orig(ar, ai, br, bi, rounding, None)
        out = orig(ar, ai, br, bi, rounding, wl)
        count["n"] += int(sum(np.sum(np.asarray(u) != np.asarray(v)) for u, v in zip(free, out)))
        return out

    monkeypatch.setattr(R, "_addsub_half", wrapped)
    return count


def test_psd_stress_segments_hit_the_butterfly_clamp(clamp_counter):
    V.psd_vectors(12, stress_from=99)
    assert clamp_counter["n"] == 0
    V.psd_vectors(12)
    assert clamp_counter["n"] > 0


def test_bispectrum_stress_frames_hit_the_butterfly_clamp(clamp_counter):
    V.bispec_vectors(12, stress_from=10**6)
    assert clamp_counter["n"] == 0
    V.bispec_vectors(12)
    assert clamp_counter["n"] > 0


def test_stress_input_reaches_both_extreme_codes():
    for x in (V.psd_vectors(12)["x"], V.bispec_vectors(12)["x"]):
        assert x.min() == -(1 << 11) and x.max() == (1 << 11) - 1
