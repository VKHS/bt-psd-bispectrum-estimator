"""Figs. 13-15 of the paper show symmetric 128-point windows: a pixel measurement of the printed stems (data/paper_figs13_15_stems_read.csv)."""
import csv
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python" / "scripts"))

import measure_paper_window_stems as M  # noqa: E402


def _committed():
    with (ROOT / "data" / "paper_figs13_15_stems_read.csv").open() as f:
        return list(csv.DictReader(f))


def test_committed_csv_is_what_the_script_measures():
    rows = _committed()
    fresh = M.rows()
    assert [r["window"] for r in rows] == [r[0] for r in fresh]
    for r, f in zip(rows, fresh):
        assert float(r["read_centre_minus_neighbour"]) == pytest.approx(f[3], abs=0.01)
        assert r["source"].startswith("read from the printed panel")


@pytest.mark.parametrize("row", _committed(), ids=lambda r: r["window"])
def test_printed_stems_are_those_of_the_symmetric_window(row):
    read = float(row["read_centre_minus_neighbour"])
    assert abs(read - float(row["symmetric_window"])) < 0.15
    assert abs(read - float(row["periodic_window"])) > 0.3
