"""Check the log of a Simulink run of the Model Composer design against the bit-exact model.

    PYTHONPATH=python python3 python/scripts/compare_mc_result.py modelcomposer/work/bt_log.mat
    PYTHONPATH=python python3 python/scripts/compare_mc_result.py --csv log.csv --wl 12 --nchunks 5

The ``.mat`` file written by ``modelcomposer/matlab/bt_mc_run.m`` holds the struct ``bt_log`` with one column vector per signal
(one entry per clock cycle): ``in_code`` (integer code at the output of the Gateway In), ``in_valid``, ``rd_addr``, ``rd_data``,
``done`` and the scalars ``wl``, ``nchunks``.  A CSV file needs the columns ``in_code,in_valid,rd_addr,rd_data,done`` (header
line).  The expected spectrum is computed from the logged input codes, so it does not depend on the source block or the rounding of
the gateway.  Exit status 0: all 1024 bins equal the model; 1: they do not (or the log is unusable); 2: the file cannot be read.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from btpsd import mcmodel as M

SIGNALS = ("in_code", "in_valid", "rd_addr", "rd_data", "done")


def _values(obj):
    """Numeric array out of a logged variable: a plain array, or a Simulink structure ``signals.values`` / ``values``."""
    if hasattr(obj, "signals"):
        sg = obj.signals
        obj = sg[0].values if isinstance(sg, np.ndarray) and sg.dtype == object else sg.values
    elif hasattr(obj, "values"):
        obj = obj.values
    return np.asarray(obj, dtype=float).ravel()


def load_mat(path: Path):
    from scipy.io import loadmat

    d = loadmat(str(path), squeeze_me=True, struct_as_record=False)
    log = d.get("bt_log", d)
    get = (lambda k: getattr(log, k)) if hasattr(log, "_fieldnames") else (lambda k: log[k])
    try:
        sig = {k: _values(get(k)) for k in SIGNALS}
    except (KeyError, AttributeError) as e:
        raise ValueError(f"{path}: missing signal {e}; expected the struct bt_log with {', '.join(SIGNALS)}") from e
    wl = nch = None
    for key, name in (("wl", "wl"), ("nchunks", "nchunks")):
        try:
            v = float(np.asarray(get(key)).ravel()[0])
        except (KeyError, AttributeError, IndexError):
            continue
        if name == "wl":
            wl = int(v)
        else:
            nch = int(v)
    return sig, wl, nch


def load_csv(path: Path):
    data = np.genfromtxt(path, delimiter=",", names=True)
    names = data.dtype.names or ()
    miss = [k for k in SIGNALS if k not in names]
    if miss:
        raise ValueError(f"{path}: missing column(s) {', '.join(miss)}; header must contain {', '.join(SIGNALS)}")
    return {k: np.asarray(data[k], dtype=float) for k in SIGNALS}, None, None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("log", nargs="?", type=Path, help="MAT file with the struct bt_log")
    ap.add_argument("--csv", type=Path, default=None, help="CSV log instead of a MAT file")
    ap.add_argument("--wl", type=int, default=None, help="word length (default: from the log, else 12)")
    ap.add_argument("--nchunks", type=int, default=None, help="chunks per run (default: from the log, else 5)")
    ap.add_argument("--round", type=int, default=1, choices=(0, 1), help="ROUND parameter of the block")
    ap.add_argument("--lag", type=int, default=1, help="clock cycles from rd_addr to rd_data in the log (1 for the block itself)")
    a = ap.parse_args(argv)
    if (a.log is None) == (a.csv is None):
        ap.error("give either a MAT file or --csv")
    try:
        sig, wl, nch = load_csv(a.csv) if a.csv else load_mat(a.log)
    except (OSError, ValueError) as e:
        print(f"cannot read the log: {e}", file=sys.stderr)
        return 2
    cfg = M.McConfig(wl=a.wl or wl or 12, nchunks=a.nchunks or nch or 5, round=a.round)
    try:
        rc = M.check_logs(sig["in_code"], sig["in_valid"], sig["rd_addr"], sig["rd_data"], sig["done"], cfg, lag=a.lag)
    except ValueError as e:
        print(f"FAIL: the log cannot be evaluated: {e}")
        return 1
    print(f"configuration: WL {cfg.wl}, {cfg.nchunks} chunks, {cfg.segments} segments, avg_shift {cfg.avg_shift}, Hann window, "
          f"rounding {'to nearest' if cfg.round else 'truncation'}")
    print(rc.summary())
    for b, got, want in rc.mismatches:
        print(f"  bin {b}: got {got}, expected {want}")
    if rc.ok and rc.peak_bin is not None:
        print(f"largest bin of the lower half: {rc.peak_bin} (P_welch = rd_data * {cfg.scale:.9g})")
    if not rc.ok and rc.lag == 1 and rc.message == "":
        # a hint for a wrong latency: does another lag reproduce the model?
        for lag in (0, 2, 3):
            try:
                alt = M.check_logs(sig["in_code"], sig["in_valid"], sig["rd_addr"], sig["rd_data"], sig["done"], cfg, lag=lag)
            except ValueError:
                continue
            if alt.ok:
                print(f"hint: the log matches the model with --lag {lag}; the address-to-data latency in the design differs from 1")
                break
    return 0 if rc.ok else 1


if __name__ == "__main__":
    sys.exit(main())
