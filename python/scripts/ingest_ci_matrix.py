"""Store the result table of the CI job ``rtl-matrix`` as ``results/rtl_config_matrix.csv`` with its provenance.

``run_rtl_matrix.py`` prints one ``MCSV`` line per configuration (the RTL against the model, bit for bit).  The simulator
runs in GitHub Actions, in several shard jobs; this script reads the downloaded logs of the shards, checks that the rows
are exactly the configurations of ``run_rtl_matrix.py`` (each once) and writes

    results/rtl_config_matrix.csv               one row per configuration
    results/rtl_config_matrix_provenance.json   run id, commit, simulator

Usage::

    for JOB_ID in ...; do gh run view RUN_ID --log --job $JOB_ID > matrix_$JOB_ID.log; done
    python scripts/ingest_ci_matrix.py matrix_*.log --run RUN_ID --commit SHA --created 2026-11-20T20:08:11Z

The table holds what the *RTL simulation* found; the model only supplies the expected values.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

from _common import RESULTS
from run_rtl_matrix import CSV_HEADER, configurations

LINE = re.compile(r"Z MCSV,(.*?)\s*$")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("logs", type=Path, nargs="+", help="log of each shard job")
    ap.add_argument("--run", required=True, help="GitHub Actions run id")
    ap.add_argument("--commit", required=True, help="commit the run simulated")
    ap.add_argument("--created", required=True, help="time the run was created (UTC, ISO 8601)")
    ap.add_argument("--repo", default="VKHS/bt-psd-bispectrum-estimator")
    args = ap.parse_args()

    found = [m.group(1).split(",") for log in args.logs for m in map(LINE.search, log.read_text(errors="replace").splitlines()) if m]
    expected = [c.tag for c in configurations()]
    tags = sorted(r[1] for r in found)
    if tags != sorted(expected):
        print(f"the logs have the configurations {tags}, expected each of {expected} once", file=sys.stderr)
        return 1
    rows = sorted(found, key=lambda r: expected.index(r[1]))
    if any(len(r) != len(CSV_HEADER) for r in rows):
        print("a row has the wrong number of columns", file=sys.stderr)
        return 1
    bad = [r[1] for r in rows if r[-1] != "PASS"]
    if bad:
        print(f"not passing: {bad}", file=sys.stderr)
        return 1

    out = RESULTS / "rtl_config_matrix.csv"
    with out.open("w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(CSV_HEADER)
        w.writerows(rows)
    prov = {
        "what": "RTL simulation (Icarus Verilog) of tb/tb_psd_top.v and tb/tb_bispectrum_top.v for the configurations of "
                "python/scripts/run_rtl_matrix.py; the Python model only supplies the expected values",
        "repository": args.repo,
        "workflow": "ci / job rtl-matrix",
        "run_id": int(args.run),
        "run_url": f"https://github.com/{args.repo}/actions/runs/{args.run}",
        "commit": args.commit,
        "run_created_utc": args.created,
        "simulator": "Icarus Verilog, apt package of the ubuntu-latest runner (version not printed in the log)",
        "configurations": len(rows),
        "passed": len(rows),
        "columns": {
            "round": "RTL parameter ROUND (1: round half up, 0: truncate)",
            "squarer": "serial sum-of-odds squarer or multiplier (RTL parameter SQ_SERIAL; psd only)",
            "segments": "PSD segments or bispectrum frames of the vectors",
            "words_checked": "result words compared with the model: 1024 PSD bins, or 32 x 32 bispectrum bins "
                             "(real and imaginary part compared together)",
            "model_clamped_diffs": "rounded butterfly differences that the model's clamp changed while building the "
                                   "vectors (reachable only with rounding)",
        },
    }
    (RESULTS / "rtl_config_matrix_provenance.json").write_text(json.dumps(prov, indent=2) + "\n")
    print(f"wrote {out.relative_to(RESULTS.parent)} ({len(rows)} rows) and its provenance")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
