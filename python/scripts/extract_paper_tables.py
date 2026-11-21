"""Extract Tables 1, 2, 4, 5 and 7 of the paper from its PDF, exactly as printed.

    pip install pdfplumber
    PYTHONPATH=python python3 python/scripts/extract_paper_tables.py path/to/paper.pdf

The PDF is not part of the repository.  Cells are the printed tokens; ``-`` marks an empty cell (a dash in the paper).
Tables 3, 6, 8, 9 and 10 are in ``data/`` as well (see ``data/README.md``).

Table 1 / Table 2  register-file address sequences (64 registers x clock cycles), mapped to the clock columns by the
                   horizontal position of each token (the text layer drops the dashes)
Table 4            significant taps of the window frequency responses
Table 5            dyadic micro-rotations: the printed rows next to the values ``2**-n`` and ``atan(2**-n)``
Table 7            data codes of the test-bench records
"""

from __future__ import annotations

import collections
import csv
import math
import re
import sys
from pathlib import Path

from _common import ROOT

DATA = ROOT / "data"
DASH = "-"


def _words(page):
    return page.extract_words(x_tolerance=1.5, y_tolerance=2)


def _caption_page(pdf, n: int) -> int:
    """Index of the page with the caption ``Table n`` alone on a line at the left margin."""
    for i, page in enumerate(pdf.pages):
        rows = collections.defaultdict(list)
        for w in _words(page):
            if w["x0"] < 400:
                rows[round(w["top"])].append(w)
        for ws in rows.values():
            ws.sort(key=lambda w: w["x0"])
            text = "".join(w["text"] for w in ws)
            if text == f"Table{n}" and ws[0]["x0"] < 60 and len(ws) <= 2:
                return i
    raise LookupError(f"caption of Table {n} not found")


def _lines(words, top, bottom, xmax=10_000):
    rows = collections.defaultdict(list)
    for w in words:
        if top <= w["top"] <= bottom and w["x0"] < xmax:
            rows[round(w["top"])].append(w)
    return [sorted(rows[y], key=lambda w: w["x0"]) for y in sorted(rows)]


def register_table(pdf, n: int):
    """Table 1 or 2: ``(clock labels, [(register, [cell per clock])])``."""
    page = pdf.pages[_caption_page(pdf, n)]
    words = _words(page)
    header = next(ws for ws in _lines(words, 0, 10_000) if ws[0]["text"] == "Clock")[1:]
    centres = [(w["x0"] + w["x1"]) / 2 for w in header]
    clocks = [int(w["text"]) for w in header]
    rows = []
    for ws in _lines(words, 0, 10_000):
        if not re.fullmatch(r"R\d+", ws[0]["text"]):
            continue
        cells = [None] * len(clocks)
        for w in ws[1:]:
            xc = (w["x0"] + w["x1"]) / 2
            j = min(range(len(clocks)), key=lambda k: abs(centres[k] - xc))
            if cells[j] is not None:
                raise ValueError(f"two tokens in one cell: {ws[0]['text']} clock {clocks[j]}")
            cells[j] = w["text"]
        rows.append((ws[0]["text"], [DASH if c is None else c for c in cells]))
    return clocks, rows


def table4(pdf):
    page = pdf.pages[_caption_page(pdf, 4)]
    words = _words(page)
    names = next(ws for ws in _lines(words, 0, 10_000) if ws[0]["text"] == "Windowing")[1:]
    taps = next(ws for ws in _lines(words, 0, 10_000) if ws[0]["text"] == "F.")
    values = [w for w in taps if re.fullmatch(r"\d+", w["text"])]
    assert len(names) == len(values)
    return [(a["text"], b["text"]) for a, b in zip(names, values)]


def table5(pdf):
    """Printed rows of Table 5; the degree sign is a trailing ``0`` in the text layer."""
    page = pdf.pages[_caption_page(pdf, 5)]
    words = _words(page)
    top = next(w for w in words if w["text"] == "Table" and w["x0"] < 60)["top"]
    lines = _lines(words, top, top + 60, 600)
    n_line = next(ws for ws in lines if ws[0]["text"] == "n" and ws[1]["text"] == "0")
    tg_line = next(ws for ws in lines if ws[0]["text"].startswith("tg"))
    angle_line = next(ws for ws in lines if len(ws) == 8 and ws[0]["top"] > tg_line[0]["top"]
                      and all(re.fullmatch(r"\d+(\.\d+)?", w["text"]) for w in ws))
    n = [int(w["text"]) for w in n_line[1:]]
    tg = [w["text"] for w in tg_line[1:]]
    ang = [w["text"] for w in angle_line]
    assert len(n) == len(tg) == len(ang) == 8
    return n, tg, ang


def table7(pdf):
    page = pdf.pages[_caption_page(pdf, 7)]
    words = _words(page)
    top = next(w for w in words if w["text"] == "Table" and w["x0"] < 60 and 150 < w["top"] < 260)["top"]
    lines = _lines(words, top, top + 60)
    sig_line = next(ws for ws in lines if ws[0]["text"] == "Biomedical")
    sig = sig_line[2:]
    code = next(ws for ws in lines if ws[0]["text"] == "Data" and ws[0]["top"] > sig_line[0]["top"])
    code = [w for w in code if w["x0"] > sig[0]["x0"] - 5]
    assert len(sig) == len(code)
    return [(a["text"], b["text"]) for a, b in zip(sig, code)]


def write_csv(name: str, header, rows) -> None:
    DATA.mkdir(exist_ok=True)
    with open(DATA / name, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    print(f"wrote data/{name} ({len(rows)} rows)")


def source(n: int, extra: str = "") -> str:
    return f"paper-reported, Integration 99 (2024) 102241, Table {n}{extra}"


def main(pdf_path: str) -> None:
    import pdfplumber

    pdf = pdfplumber.open(pdf_path)
    for n, name in ((1, "paper_table1_register_sequence_1.csv"), (2, "paper_table2_register_sequence_2.csv")):
        clocks, rows = register_table(pdf, n)
        write_csv(name, ["register"] + [f"clock_{c}" for c in clocks] + ["source"],
                  [[r] + c + [source(n)] for r, c in rows])
    write_csv("paper_table4_window_significant_taps.csv", ["window", "significant_taps_printed", "source"],
              [[a, b, source(4)] for a, b in table4(pdf)])
    n, tg, ang = table5(pdf)
    rows = []
    for i, t, a in zip(n, tg, ang):
        deg = a[:-1] if a.endswith("0") else a  # the degree sign appears as a trailing 0
        rows.append([i, t, a, deg, f"{2.0 ** -i:.7g}", f"{math.degrees(math.atan(2.0 ** -i)):.4f}",
                     source(5, "; printed_* columns as printed, computed_* columns analytic, computed here")])
    write_csv("paper_table5_dyadic_micro_rotations.csv",
              ["n", "printed_tg_wi", "printed_angle_token", "printed_angle_deg", "computed_2^-n", "computed_atan_deg",
               "source"], rows)
    write_csv("paper_table7_testbench_data_codes.csv", ["signal", "data_code_printed", "source"],
              [[a, b, source(7)] for a, b in table7(pdf)])


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    main(sys.argv[1])
