"""Which window form do the paper's Figs. 13-15 show?  A pixel measurement of the printed panels.

The left panel of each of Figs. 13-15 plots the DFT magnitude of a 128-point window as stems (hollow squares).  The
periodic window of this repository has the stems ``64, 32`` (Hann), ``69.12, 29.44`` (Hamming) and ``53.76, 32, 5.12``
(Blackman); the symmetric window (``numpy.hanning`` form, MATLAB default) has slightly different heights and a little
leakage outside the central bins.  The difference "central stem minus first neighbour" does not depend on the size of
the marker or on the zero of the axis, so it is read from the pixel rows of the stem tops and scaled with the
height of the axes frame (top tick label: 70, 70, 60).

Writes ``data/paper_figs13_15_stems_read.csv``.  The values are *read from the printed figure*, not printed in the paper.
"""
from __future__ import annotations

import csv

import numpy as np
from PIL import Image

from _common import DATA, ROOT
from btpsd.windows import window

PANELS = [("hann", 13, "paperfig13_window_hann.png", 70.0), ("hamming", 14, "paperfig14_window_hamming.png", 70.0),
          ("blackman", 15, "paperfig15_window_blackman.png", 60.0)]
HEADER = ["window", "paper_figure", "axis_top_value", "read_centre_minus_neighbour", "symmetric_window", "periodic_window", "source"]
SOURCE = "read from the printed panel by pixel measurement (python/scripts/measure_paper_window_stems.py); not a value printed in the paper"


def _longest_run(mask: np.ndarray) -> tuple[int, int]:
    best = cur = start = best_start = 0
    for i, v in enumerate(mask):
        if v:
            if cur == 0:
                start = i
            cur += 1
            if cur > best:
                best, best_start = cur, start
        else:
            cur = 0
    return best, best_start


def read_panel(path, top_value: float) -> float:
    """Centre stem minus first neighbour, in axis units, from the left panel of one figure."""
    dark = np.asarray(Image.open(path).convert("L")) < 160
    h, w = dark.shape
    runs = {x: _longest_run(dark[:, x]) for x in range(int(0.30 * w))}  # the left panel only
    tall = {x: r for x, r in runs.items() if r[0] > 150}
    top = min(s for _, s in tall.values())
    frame = sorted(x for x, (_, s) in tall.items() if abs(s - top) <= 2)
    left, right = frame[0], frame[-1]
    assert right - left > 300, "axes frame not found"
    bottom = tall[right][1] + tall[right][0] - 1
    px_per_unit = (bottom - top) / top_value
    inner = sorted(x for x in tall if left + 5 < x < right - 5)
    clusters, cur = [], [inner[0]]
    for x in inner[1:]:
        if x - cur[-1] <= 2:
            cur.append(x)
        else:
            clusters.append(cur)
            cur = [x]
    clusters.append(cur)
    assert len(clusters) == 3, f"expected three stems, found {len(clusters)}"
    start = [min(tall[x][1] for x in c) for c in clusters]
    assert abs(start[0] - start[2]) <= 2, "the two neighbours differ"
    centre = min(start)
    neighbour = (start[0] + start[2]) / 2.0
    return float((neighbour - centre) / px_per_unit)


def predictions(name: str, L: int = 128) -> tuple[float, float]:
    out = []
    for periodic in (False, True):
        F = np.abs(np.fft.fft(window(name, L, periodic=periodic)))
        out.append(float(F[0] - F[1]))
    return out[0], out[1]


def rows():
    res = []
    for name, fig, fn, top in PANELS:
        read = read_panel(ROOT / "docs" / "paper_figures" / fn, top)
        sym, per = predictions(name)
        res.append([name, fig, top, round(read, 2), round(sym, 2), round(per, 2), SOURCE])
    return res


def main():
    out = rows()
    with open(DATA / "paper_figs13_15_stems_read.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(HEADER)
        w.writerows(out)
    for r in out:
        print(f"Fig. {r[1]} {r[0]:9s} read {r[3]:6.2f}   symmetric {r[4]:6.2f}   periodic {r[5]:6.2f}")


if __name__ == "__main__":
    main()
