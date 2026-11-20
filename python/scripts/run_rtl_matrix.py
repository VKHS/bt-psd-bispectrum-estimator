"""Check the RTL against the model for a matrix of configurations (Icarus Verilog).

The default testbenches check one configuration of each top level (WL = 12, Hann window, ROUND = 1, serial squarer).
This script generates the golden vectors of the model for further configurations, compiles the testbenches
``tb/tb_psd_top.v`` and ``tb/tb_bispectrum_top.v`` for each of them and compares the RTL with the model bit for bit:

    PYTHONPATH=python python3 python/scripts/run_rtl_matrix.py            # all configurations (needs iverilog)
    PYTHONPATH=python python3 python/scripts/run_rtl_matrix.py --list     # the configurations, no simulator needed
    PYTHONPATH=python python3 python/scripts/run_rtl_matrix.py --only psd_wl10_hann_r1_s1
    PYTHONPATH=python python3 python/scripts/run_rtl_matrix.py --shard 2/4    # one fourth of the work (CI runs 4 jobs)

What varies is what the RTL exposes as parameters: ``WL`` (word length), the window taps ``WM2 .. WP2``, ``ROUND``
(1: round half up, 0: truncate) and ``SQ_SERIAL`` (serial sum-of-odds squarer or multiplier).  The vectors are written to
``sim/work/matrix`` (not tracked).  Every finished configuration prints one ``MCSV`` line, so that the result can be read
back from a CI log (``ingest_ci_matrix.py``); the exit status is 0 only if all configurations pass.
"""
from __future__ import annotations

import argparse
import contextlib
import csv
import re
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from _common import ROOT

from btpsd import r8fft as R
from btpsd import vectors as V
from btpsd.windows import spectral_kernel

WORK = ROOT / "sim" / "work" / "matrix"

# Compile order as in sim/Makefile (tests/test_rtl_matrix.py compares the two).
FFT_FILES = ["rtl/ram_1r1w.v", "rtl/delay_line.v", "rtl/cru_rom.v", "rtl/cru_rotate.v", "rtl/r8_butterfly.v",
             "rtl/fft512_mem.v"]
PSD_FILES = ["tb/tb_psd_top.v", *FFT_FILES, "rtl/mag_sq_odd.v", "rtl/psd_engine.v", "rtl/bt_psd_top.v"]
BSP_FILES = ["tb/tb_bispectrum_top.v", *FFT_FILES, "rtl/bispectrum_engine.v", "rtl/bt_bispectrum_top.v"]

MAX_BIN = 32
CONSTS = ROOT / "rtl" / "btpsd_consts.vh"

CSV_HEADER = ["kind", "tag", "wl", "window", "round", "squarer", "segments", "words_checked", "model_clamped_diffs",
              "result"]


@dataclass(frozen=True)
class Config:
    kind: str  # 'psd' or 'bsp'
    wl: int
    window: str
    rnd: int = 1  # RTL parameter ROUND
    serial: int = 1  # RTL parameter SQ_SERIAL (psd only)

    @property
    def tag(self) -> str:
        base = f"{self.kind}_wl{self.wl}_{self.window}_r{self.rnd}"
        return base + (f"_s{self.serial}" if self.kind == "psd" else "")

    @property
    def rounding(self) -> str:
        return "round" if self.rnd else "trunc"

    @property
    def squarer(self) -> str:
        if self.kind != "psd":
            return "-"
        return "serial" if self.serial else "multiplier"


def configurations() -> list[Config]:
    psd = [Config("psd", wl, "hann") for wl in (10, 12, 14, 16)]
    psd += [Config("psd", wl, w) for wl in (12, 10) for w in ("rect", "hamming", "blackman")]
    psd += [Config("psd", 16, "blackman")]
    psd += [Config("psd", 12, "hann", rnd=0), Config("psd", 10, "blackman", rnd=0), Config("psd", 16, "hamming", rnd=0)]
    psd += [Config("psd", 12, "hann", serial=0), Config("psd", 16, "hann", serial=0)]
    bsp = [Config("bsp", wl, "hann") for wl in (10, 12, 14, 16)]
    bsp += [Config("bsp", 12, w) for w in ("rect", "hamming", "blackman")]
    bsp += [Config("bsp", 12, "hann", rnd=0), Config("bsp", 16, "blackman", rnd=0)]
    cfgs = psd + bsp
    assert len({c.tag for c in cfgs}) == len(cfgs)
    return cfgs


def cost(c: Config) -> float:
    """Rough simulation time in seconds on a CI runner: about 22 s per PSD run plus the cycles of the serial squarer
    (21 757 at WL = 12, four times as many per two more bits) at about 5 000 cycles/s; 4 s per bispectrum run."""
    if c.kind == "bsp":
        return 4.0
    return 22.0 + (21757 * 4 ** ((c.wl - 12) / 2) / 5000 if c.serial else 0.0)


def shard(cfgs: list[Config], i: int, n: int) -> list[Config]:
    """The configurations of shard ``i`` (1-based) of ``n``: longest-processing-time assignment by :func:`cost`."""
    assert 1 <= i <= n
    load, bins = [0.0] * n, [[] for _ in range(n)]
    for c in sorted(cfgs, key=lambda c: (-cost(c), c.tag)):
        k = min(range(n), key=lambda j: (load[j], j))
        bins[k].append(c.tag)
        load[k] += cost(c)
    keep = set(bins[i - 1])
    return [c for c in cfgs if c.tag in keep]  # canonical order


@contextlib.contextmanager
def clamp_counter():
    """Count the rounded half add/subtract results that the clamp changes (see tests/test_vectors_stress.py)."""
    orig = R._addsub_half
    count = {"n": 0}

    def wrapped(ar, ai, br, bi, rounding, wl=None):
        free = orig(ar, ai, br, bi, rounding, None)
        out = orig(ar, ai, br, bi, rounding, wl)
        count["n"] += int(sum(np.sum(np.asarray(u) != np.asarray(v)) for u, v in zip(free, out)))
        return out

    R._addsub_half = wrapped
    try:
        yield count
    finally:
        R._addsub_half = orig


def header_constants() -> dict[str, int]:
    """The ``define`` constants of rtl/btpsd_consts.vh (filter taps FD_H*, window taps WIN_*)."""
    return {m.group(1): int(m.group(2)) for m in re.finditer(r"`define (\w+) (-?\d+)\s*$", CONSTS.read_text(), re.M)}


def window_taps(window: str, win_frac: int) -> list[int]:
    """Integer kernel taps d = -2 .. 2 of the model, checked against the constants of the RTL header."""
    kernel = spectral_kernel(window)
    taps = [int(np.floor(kernel.get(d, 0.0) * (1 << win_frac) + 0.5)) for d in (-2, -1, 0, 1, 2)]
    hdr = header_constants()
    name = window.upper()
    in_header = [hdr[f"WIN_{name}_{k}"] for k in ("M2", "M1", "0", "P1", "P2")]
    assert taps == in_header, f"{window}: model taps {taps} differ from btpsd_consts.vh {in_header}"
    assert win_frac == hdr["WIN_FRAC"], "window fraction bits of model and RTL differ"
    return taps


def _params(c: Config, taps: list[int]) -> str:
    p = f".WL({c.wl}), .ROUND({c.rnd}), "
    p += f".SQ_SERIAL({c.serial}), " if c.kind == "psd" else f".LOGB(5), "
    return p + ".WM2({}), .WM1({}), .W0({}), .WP1({}), .WP2({})".format(*taps)


def make_vectors(c: Config) -> tuple[dict, int]:
    """Write the vectors and the include file of one configuration; return the model data and the clamp count."""
    WORK.mkdir(parents=True, exist_ok=True)
    rel = WORK.relative_to(ROOT).as_posix()
    with clamp_counter() as count:
        if c.kind == "psd":
            d = V.psd_vectors(c.wl, c.window, rounding=c.rounding)
        else:
            d = V.bispec_vectors(c.wl, c.window, max_bin=MAX_BIN, rounding=c.rounding)
    cfg = d["cfg"]
    taps = window_taps(c.window, cfg.win_frac)
    if c.kind == "psd":
        hdr = header_constants()
        assert [int(v) for v in cfg.fd_int()] == [hdr[f"FD_H{i}"] for i in range(6)], "filter taps of model and RTL differ"
    if c.kind == "psd":
        V.write_mem(WORK / f"{c.tag}_x.mem", d["x"], c.wl)
        V.write_mem(WORK / f"{c.tag}_P.mem", d["P"], 2 * c.wl + 8)
        lines = [f"`define PSD_S {d['S']}", f"`define PSD_AVG_SHIFT {d['avg_shift']}",
                 f"`define PSD_NCHUNKS {d['nchunks']}", f"`define PSD_SQUARER_CYCLES {d['squarer_cycles']}",
                 f"`define TB_XFILE \"{rel}/{c.tag}_x.mem\"", f"`define TB_PFILE \"{rel}/{c.tag}_P.mem\""]
    else:
        bw = 2 * c.wl + 10
        V.write_mem(WORK / f"{c.tag}_x.mem", d["x"], c.wl)
        V.write_mem(WORK / f"{c.tag}_B_re.mem", d["B_re"][:MAX_BIN, :MAX_BIN].ravel(), bw)
        V.write_mem(WORK / f"{c.tag}_B_im.mem", d["B_im"][:MAX_BIN, :MAX_BIN].ravel(), bw)
        lines = [f"`define BSP_S {d['S']}", f"`define BSP_AVG_SHIFT {d['avg_shift']}", f"`define BSP_MAX_BIN {MAX_BIN}",
                 f"`define BSP_NSAMPLES {len(d['x'])}",
                 f"`define TB_XFILE \"{rel}/{c.tag}_x.mem\"", f"`define TB_REFILE \"{rel}/{c.tag}_B_re.mem\"",
                 f"`define TB_IMFILE \"{rel}/{c.tag}_B_im.mem\""]
    lines += [f"`define TB_WL {c.wl}", f"`define TB_PARAMS {_params(c, taps)}"]
    (WORK / f"{c.tag}.vh").write_text(f"// generated by python/scripts/run_rtl_matrix.py for {c.tag}\n" + "\n".join(lines) + "\n")
    return d, count["n"]


def simulate(c: Config, iverilog: str, vvp: str) -> tuple[str, int]:
    """Compile and run the testbench of one configuration; return (result, words checked)."""
    files = PSD_FILES if c.kind == "psd" else BSP_FILES
    out = WORK / f"{c.tag}.vvp"
    cmd = [iverilog, "-g2005", "-Wall", "-I", "rtl", "-I", "sim/vectors", "-I", WORK.relative_to(ROOT).as_posix(),
           f'-DTB_VH="{c.tag}.vh"', "-o", out.relative_to(ROOT).as_posix(), *files]
    comp = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    if comp.returncode != 0:
        print(f"[{c.tag}] compile failed:\n{comp.stdout[-1500:]}{comp.stderr[-1500:]}", file=sys.stderr)
        return "ERROR", 0
    run = subprocess.run([vvp, out.relative_to(ROOT).as_posix()], cwd=ROOT, capture_output=True, text=True, timeout=900)
    (WORK / f"{c.tag}.log").write_text(run.stdout)
    m = re.search(r"TB_RESULT: (PASS|FAIL) \w+ \(([^)]*)\)", run.stdout)
    if not m:
        print(f"[{c.tag}] no TB_RESULT line:\n{run.stdout[-1500:]}{run.stderr[-800:]}", file=sys.stderr)
        return "ERROR", 0
    if m.group(1) == "FAIL":
        print(f"[{c.tag}] FAIL: " + " | ".join(ln for ln in run.stdout.splitlines() if "MISMATCH" in ln or "FAIL" in ln)[:1500],
              file=sys.stderr)
        return "FAIL", 0
    words = re.match(r"(\d+)(?: x (\d+))?", m.group(2))
    n = int(words.group(1)) * (int(words.group(2)) if words.group(2) else 1)
    return "PASS", n


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", action="store_true", help="print the configurations and exit")
    ap.add_argument("--only", nargs="*", default=None, help="run only these tags")
    ap.add_argument("--shard", default=None, metavar="I/N", help="run shard I of N (1-based), balanced by estimated time")
    ap.add_argument("--csv", type=Path, default=None, help="also write the result table to this CSV file")
    ap.add_argument("--iverilog", default="iverilog")
    ap.add_argument("--vvp", default="vvp")
    args = ap.parse_args()

    cfgs = configurations()
    if args.only:
        unknown = set(args.only) - {c.tag for c in cfgs}
        if unknown:
            print(f"unknown tags: {sorted(unknown)}", file=sys.stderr)
            return 2
        cfgs = [c for c in cfgs if c.tag in args.only]
    if args.shard:
        i, n = (int(v) for v in args.shard.split("/"))
        cfgs = shard(cfgs, i, n)
    if args.list:
        for c in cfgs:
            print(c.tag)
        return 0

    rows = []
    t0 = time.time()
    for c in cfgs:
        d, clamps = make_vectors(c)
        result, words = simulate(c, args.iverilog, args.vvp)
        row = [c.kind, c.tag, c.wl, c.window, c.rnd, c.squarer, d["S"], words, clamps, result]
        rows.append(row)
        print("MCSV," + ",".join(str(v) for v in row), flush=True)
    print("\n" + "".join(f"{h:>20}" if i != 1 else f"{h:<28}" for i, h in enumerate(CSV_HEADER)))
    for r in rows:
        print("".join(f"{str(v):>20}" if i != 1 else f"{str(v):<28}" for i, v in enumerate(r)))
    failed = [r[1] for r in rows if r[-1] != "PASS"]
    print(f"\n{len(rows) - len(failed)} of {len(rows)} configurations match the model bit for bit "
          f"({time.time() - t0:.0f} s)")
    if args.csv:
        with args.csv.open("w", newline="") as f:
            w = csv.writer(f, lineterminator="\n")
            w.writerow(CSV_HEADER)
            w.writerows(rows)
    if failed:
        print(f"not passing: {failed}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
