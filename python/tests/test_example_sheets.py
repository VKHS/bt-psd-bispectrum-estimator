"""The numbers printed on the example sheets (fig21 windows, fig22 PSD examples) come from the model and match the committed CSV."""
import csv
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python" / "scripts"))

from _figdata import PSD_EXAMPLE_HEADER, psd_example_rows  # noqa: E402
from btpsd.estimator import HwConfig  # noqa: E402
from btpsd.windows import TABLE4_SIGNIFICANT_TAPS  # noqa: E402


def test_psd_example_numbers_match_the_committed_csv():
    with (ROOT / "results" / "model_psd_examples_wl10.csv").open() as f:
        rows = list(csv.reader(f))
    assert rows[0] == PSD_EXAMPLE_HEADER
    got = psd_example_rows()
    assert len(rows) - 1 == len(got) == 3
    for line, ref in zip(rows[1:], got):
        assert line[:5] == [str(v) for v in ref[:5]]
        assert float(line[5]) == pytest.approx(ref[5], abs=0.01)  # SQNR against the float Welch average, dB
        assert float(line[6]) == pytest.approx(ref[6], abs=1e-3)  # 1 LSB of the PSD word, dB (analytic)
        assert abs(int(line[7]) - ref[7]) <= 2  # bins whose integer result is 0


def test_psd_lsb_matches_the_analytic_scale():
    """One LSB of the accumulated power word is ``1024 * 2**avg_shift / (S * 2**(2 WL - 2))`` in the scale of the float reference."""
    import math

    wl, S, shift = 10, 9, 4
    expected = 10 * math.log10(1024 * 2**shift / (S * 2 ** (2 * wl - 2)))
    with (ROOT / "results" / "model_psd_examples_wl10.csv").open() as f:
        rows = list(csv.DictReader(f))
    assert all(float(r["psd_lsb_db"]) == pytest.approx(expected, abs=1e-3) for r in rows)


@pytest.mark.parametrize("name", ["rect", "hann", "hamming", "blackman"])
def test_integer_window_kernel_keeps_every_tap(name):
    """The 10-bit kernel constants of the window sheet: as many non-zero taps as the paper's Table 4 lists for the window."""
    ints = HwConfig(window=name).win_int()
    assert len(ints) == TABLE4_SIGNIFICANT_TAPS[name]
    assert all(v != 0 for v in ints.values())
