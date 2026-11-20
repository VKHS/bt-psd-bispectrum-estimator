"""Loaders for the data tables stored in ``data/``.

Everything read here is labelled by its ``source`` column:

* ``paper-reported``  numbers printed in the paper (Integration 99 (2024) 102241),
  transcribed from the PDF text layer.  They are *reference data*, never produced
  by this repository's code.
* ``table6_printed.csv`` additionally carries consistency flags computed here.

Measured/modelled results of this repository live in ``results/`` and are
never mixed into these files.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

__all__ = [
    "DATA_DIR",
    "table3",
    "table6_printed",
    "table8",
    "table9",
    "table10",
    "table9_consistency",
]

DATA_DIR = Path(__file__).resolve().parents[2] / "data"


def _rows(name: str) -> list[dict]:
    with open(DATA_DIR / name, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def table3() -> dict[int, float]:
    """Paper Table 3: SNR (dB) versus number of FD-filter taps (3..10)."""
    return {int(r["fd_taps"]): float(r["snr_db"]) for r in _rows("paper_table3_fd_taps_snr.csv")}


def table6_printed() -> list[dict]:
    """Printed Table 6 (CRU coefficients) with numeric fields converted and consistency flags."""
    out = []
    for r in _rows("table6_printed.csv"):
        d = dict(r)
        for k in ("row", "a", "b", "n_mux", "n_add", "shift_size", "shift_num", "mem_buf",
                  "r_n_add", "r_shift_size", "r_shift_num", "r_radius"):
            d[k] = int(float(r[k]))
        for k in ("k_label", "k_coef", "err", "angle_steps_from_ab", "label_vs_ab_steps"):
            d[k] = float(r[k])
        d["k_nearest"] = int(float(r["k_nearest"]))
        out.append(d)
    return out


def table8() -> list[dict]:
    """Paper Table 8: SQNR (dB) at WL = 10 for 9 signals x 4 windows x methods (blank = not reported)."""
    out = []
    for r in _rows("paper_table8_sqnr_wl10.csv"):
        out.append({**r, "sqnr_db": float(r["sqnr_db"]) if r["sqnr_db"] != "" else None})
    return out


def table9() -> list[dict]:
    """Paper Table 9: FPGA resources / throughput / power / latency (Artix-7)."""
    out = []
    for r in _rows("paper_table9_fpga_comparison.csv"):
        out.append({
            "method": r["method"], "dsp": int(r["dsp"]), "lut": int(r["lut"]), "ff": int(r["ff"]),
            "throughput_msps": int(r["throughput_msps"]), "R": float(r["throughput_ratio_R"]),
            "freq_mhz": int(r["freq_mhz"]), "power_mw": int(r["power_mw"]),
            "latency_cycles": int(r["latency_cycles"]), "pdp": int(r["pdp"]),
        })
    return out


def table10() -> list[dict]:
    """Paper Table 10: bispectrum SQNR (dB) for WL = 10..18, cases 1-3."""
    return [{"signal": r["signal"], "case": int(r["case"]), "wl": int(r["wl"]),
             "sqnr_db": float(r["sqnr_db"])} for r in _rows("paper_table10_bispectrum_sqnr.csv")]


def table9_consistency() -> dict:
    """Internal consistency of Table 9: ``PDP = power x latency`` and ``R = throughput / frequency``."""
    bad_pdp, bad_r = [], []
    for r in table9():
        if r["pdp"] != r["power_mw"] * r["latency_cycles"]:
            bad_pdp.append(r["method"])
        if abs(r["R"] - r["throughput_msps"] / r["freq_mhz"]) > 0.011:
            bad_r.append(r["method"])
    return {"pdp_mismatch": bad_pdp, "R_mismatch": bad_r}
