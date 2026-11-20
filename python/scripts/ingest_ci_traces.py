"""Turn the trace lines of a CI simulation log into the CSV files under ``results/``.

The self-checking testbenches in ``tb/`` print ``TRACE``, ``OCC`` and ``HIST`` lines next to their pass/fail output
(see the comments in the testbenches).  The simulator runs in GitHub Actions (job ``rtl-sim``); this script reads the
downloaded log, checks that every testbench passed and writes

    results/rtl_trace_fft512.csv         per-cycle signals of the first 512-point transform
    results/rtl_trace_squarer.csv        per-cycle signals of three squarer examples
    results/rtl_state_occupancy.csv      cycles per state of the top levels and engines (whole testbench run)
    results/rtl_squarer_wait_hist.csv    residence time of the PSD engine in W_S, binned by floor(log2(cycles))
    results/rtl_sim_provenance.json      where the numbers come from

Usage::

    gh run view RUN_ID --log --job JOB_ID > ci.log
    python scripts/ingest_ci_traces.py ci.log --run RUN_ID --commit SHA --created 2026-11-20T19:30:49Z

These files hold what the *RTL simulation* measured; nothing in them comes from the Python models.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

from _common import RESULTS
from btpsd.rtlparse import RTL_DIR, parse_all

LINE = re.compile(r"Z (TRACE|OCC|HIST),(.*?)\s*$")
EXPECTED_TESTBENCHES = ["tb_bispectrum_top", "tb_cru_rotate", "tb_fft512_mem", "tb_mag_sq_odd", "tb_psd_top",
                        "tb_r8_butterfly"]
MODULE_OF = {"psd_top": "bt_psd_top", "psd_engine": "psd_engine", "bsp_top": "bt_bispectrum_top",
             "bsp_engine": "bispectrum_engine"}


def _write(path: Path, header: list[str], rows: list[list]) -> None:
    with path.open("w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)
    print(f"wrote {path.relative_to(RESULTS.parent)} ({len(rows)} rows)")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("log", type=Path)
    ap.add_argument("--run", required=True, help="GitHub Actions run id")
    ap.add_argument("--commit", required=True, help="commit the run simulated")
    ap.add_argument("--created", required=True, help="time the run was created (UTC, ISO 8601)")
    ap.add_argument("--repo", default="VKHS/bt-psd-bispectrum-estimator")
    args = ap.parse_args()

    text = args.log.read_text(errors="replace")
    fft, sq, occ, hist = [], [], [], []
    for ln in text.splitlines():
        m = LINE.search(ln)
        if not m:
            continue
        kind, rest = m.group(1), m.group(2).split(",")
        if kind == "TRACE" and rest[0] == "fft":
            fft.append([int(v) for v in rest[1:]])
        elif kind == "TRACE" and rest[0] == "sq":
            sq.append([int(v) for v in rest[1:]])
        elif kind == "OCC":
            occ.append((rest[0], int(rest[1]), int(rest[2])))
        elif kind == "HIST":
            hist.append((int(rest[1]), int(rest[2])))

    results = sorted({m.group(1) for m in re.finditer(r"TB_RESULT: PASS (\w+)", text)})
    failed = re.findall(r"TB_RESULT: FAIL.*", text)
    if failed or results != EXPECTED_TESTBENCHES:
        print(f"testbench results: passed {results}, failed {failed}", file=sys.stderr)
        return 1
    if not (fft and sq and occ and hist):
        print("the log has no trace lines", file=sys.stderr)
        return 1

    mods = parse_all(RTL_DIR)
    fft_states = {v: k for k, v in mods["fft512_mem"].states.items()}
    _write(RESULTS / "rtl_trace_fft512.csv", ["cycle", "ld_en", "start", "busy", "done", "iss_v", "iss_p", "iss_b", "st", "state"],
           [r + [fft_states[r[-1]]] for r in fft])
    _write(RESULTS / "rtl_trace_squarer.csv", ["example", "cycle", "start", "busy", "done", "re", "im", "p"], sq)
    rows = []
    for tb, value, cycles in occ:
        names = {v: k for k, v in mods[MODULE_OF[tb]].states.items()}
        rows.append([tb, MODULE_OF[tb], value, names[value], cycles])
    _write(RESULTS / "rtl_state_occupancy.csv", ["testbench", "module", "state_value", "state", "cycles"], rows)
    _write(RESULTS / "rtl_squarer_wait_hist.csv", ["log2_bin", "cycles_from", "cycles_to", "bins"],
           [[k, 2 ** k, 2 ** (k + 1) - 1, n] for k, n in hist])

    took = sorted({int(v) for v in re.findall(r"transform took (\d+) cycles", text)})
    prov = {
        "what": "Measurements of the RTL simulation (Icarus Verilog) of the self-checking testbenches in tb/; "
                "nothing here comes from the Python models",
        "repository": args.repo,
        "workflow": "ci / job rtl-sim (make -C sim -k all)",
        "run_id": int(args.run),
        "run_url": f"https://github.com/{args.repo}/actions/runs/{args.run}",
        "commit": args.commit,
        "run_created_utc": args.created,
        "simulator": "Icarus Verilog, apt package of the ubuntu-latest runner (version not printed in the log)",
        "testbenches_passed": results,
        "fft512_first_issue_to_done_cycles": took,
        "files": {
            "rtl_trace_fft512.csv": "tb_fft512_mem, chunk 0: signals sampled at the rising edge, cycle 0 = start high",
            "rtl_trace_squarer.csv": "tb_mag_sq_odd, directed examples (3, -5), (0, 0), (-4, 2): cycle 0 = start high",
            "rtl_state_occupancy.csv": "tb_psd_top (WL = 12, 4 chunks, 3 segments, one clear) and tb_bispectrum_top "
                                       "(WL = 12, 5 frames, one clear), cycles per state from reset release to the last result",
            "rtl_squarer_wait_hist.csv": "tb_psd_top, bins by floor(log2(cycles in W_S)) per visit",
        },
    }
    out = RESULTS / "rtl_sim_provenance.json"
    out.write_text(json.dumps(prov, indent=2) + "\n")
    print(f"wrote {out.relative_to(RESULTS.parent)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
