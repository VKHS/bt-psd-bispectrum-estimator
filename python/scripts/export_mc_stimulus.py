"""Write a stimulus for the Model Composer / Simulink flow.

    PYTHONPATH=python python3 python/scripts/export_mc_stimulus.py --source pink --wl 12 --nchunks 5

Sources (all synthetic):

``white``, ``pink``, ``ar2``, ``sines``, ``eeg_like``, ``eeg_coupled``
    the generators of ``btpsd.signals`` (peak ``--peak`` of full scale);
``tone``
    ``--tone-amp * sin(2 pi f n + phase) + Gaussian noise`` with ``f = --tone-bin / 1024`` cycles per sample and standard
    deviation ``--noise-rms``: what the Simulink *Sine Wave* and *Random Number* blocks add up to;
``siggen``
    the first samples of the in-fabric generator ``sig_gen.v`` (``btpsd.siggen``) for ``--tone-bin`` and ``--noise-shift``
    (``--noise-shift -1`` for the tone alone), so that a Simulink run and ``bt_psd_selftest`` can be compared.

The samples are quantised to integer codes of ``Fix_<wl>_<wl-1>`` (round half up, saturating) and written as the doubles
``code / 2**(wl-1)``, which a Gateway In block of that type converts back to exactly the same code whatever its quantisation
mode.  Files written to ``--out-dir`` (default ``modelcomposer/stimulus``):

``<name>.mat``
    ``bt_stim`` (N x 2: time in clock cycles, value; the *From Workspace* data), ``bt_codes`` (N x 1 integer codes),
    ``bt_wl``, ``bt_nchunks``, ``bt_nsamp``, ``bt_stop_cycles`` (recommended simulation length in clock cycles);
``<name>.csv``
    ``index,code,value``;
``<name>.json``
    the configuration, the predicted cycle count and the expected spectrum words (``btpsd.mcmodel``).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.io import savemat

from _common import ROOT

from btpsd import mcmodel as M
from btpsd import signals, siggen
from btpsd.r8fft import to_fixed

DEFAULT_DIR = ROOT / "modelcomposer" / "stimulus"
SOURCES = sorted(signals.SIGNALS) + ["tone", "siggen"]


def make_record(source: str, cfg: M.McConfig, seed: int = 1, peak: float = 0.5, tone_bin: float = 100.37, tone_amp: float = 0.5,
                noise_rms: float = 0.05, noise_shift: int = 3) -> np.ndarray:
    """Integer codes of the first ``cfg.nsamp`` samples of ``source``."""
    n = cfg.nsamp
    if source in signals.SIGNALS:
        x = signals.SIGNALS[source](n, seed=seed, peak=peak)
        codes, _ = to_fixed(np.asarray(x, dtype=float), cfg.wl, "round")
        return codes
    if source == "tone":
        rng = np.random.default_rng(seed)
        k = np.arange(n)
        x = tone_amp * np.sin(2 * np.pi * (tone_bin / 1024.0) * k + rng.uniform(0, 2 * np.pi)) + noise_rms * rng.standard_normal(n)
        codes, _ = to_fixed(np.clip(x, -1.0, 1.0 - 2.0 ** -(cfg.wl - 1)), cfg.wl, "round")
        return codes
    if source == "siggen":
        return siggen.generate(n, siggen.pinc_for(tone_bin), None if noise_shift < 0 else noise_shift, cfg.wl)
    raise ValueError(f"unknown source {source!r}; choose one of {', '.join(SOURCES)}")


def write_stimulus(codes: np.ndarray, cfg: M.McConfig, out_dir: Path, name: str, extra: dict | None = None) -> list[Path]:
    """Write ``<name>.mat/.csv/.json`` for a record of integer codes; returns the paths."""
    out_dir.mkdir(parents=True, exist_ok=True)
    res = M.expected(codes, cfg)
    values = M.codes_to_gateway_values(codes, cfg.wl)
    n = len(codes)
    t = np.arange(n, dtype=float)
    stop = M.recommended_stop(res, cfg)
    mat = out_dir / f"{name}.mat"
    savemat(mat, {"bt_stim": np.column_stack([t, values]), "bt_codes": codes.astype(np.int32).reshape(-1, 1),
                  "bt_wl": float(cfg.wl), "bt_nchunks": float(cfg.nchunks), "bt_nsamp": float(n), "bt_stop_cycles": float(stop)})
    csv = out_dir / f"{name}.csv"
    with open(csv, "w") as f:
        f.write("index,code,value\n")
        for i, (c, v) in enumerate(zip(codes, values)):
            f.write(f"{i},{int(c)},{float(v)!r}\n")
    js = out_dir / f"{name}.json"
    meta = {
        "what": "stimulus for the Model Composer flow (synthetic); expected values are those of the bit-accurate Python model, "
                "not measurements of any hardware",
        "wl": cfg.wl, "nchunks": cfg.nchunks, "segments": cfg.segments, "avg_shift": cfg.avg_shift, "nsamp": n,
        "window": cfg.window, "round": cfg.round, "accw": cfg.accw, "power_scale": cfg.scale,
        "predicted_cycles_to_done": M.expected_cycles(res, cfg), "recommended_stop_cycles": stop,
        "squarer_cycles": res.squarer_cycles,
        "expected_P_int": [int(v) for v in res.P_int],
    }
    if extra:
        meta["source"] = extra
    js.write_text(json.dumps(meta, indent=1) + "\n")
    return [mat, csv, js]


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--source", default="pink", choices=SOURCES)
    ap.add_argument("--wl", type=int, default=12)
    ap.add_argument("--nchunks", type=int, default=5, help="512-sample chunks per run; segments = nchunks - 1 (must match the block)")
    ap.add_argument("--round", type=int, default=1, choices=(0, 1), help="ROUND parameter of the block (matters for the expected result)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--peak", type=float, default=0.5, help="peak of the btpsd.signals sources, fraction of full scale")
    ap.add_argument("--tone-bin", type=float, default=100.37, help="tone frequency in bins of a 1024-point spectrum (tone, siggen)")
    ap.add_argument("--tone-amp", type=float, default=0.5)
    ap.add_argument("--noise-rms", type=float, default=0.05, help="noise standard deviation, fraction of full scale (tone)")
    ap.add_argument("--noise-shift", type=int, default=3, help="noise amplitude 2^(wl-1-shift), -1: none (siggen)")
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_DIR)
    ap.add_argument("--name", default=None)
    a = ap.parse_args(argv)
    cfg = M.McConfig(wl=a.wl, nchunks=a.nchunks, round=a.round)
    codes = make_record(a.source, cfg, a.seed, a.peak, a.tone_bin, a.tone_amp, a.noise_rms, a.noise_shift)
    name = a.name or f"bt_stim_{a.source}_wl{a.wl}_c{a.nchunks}"
    extra = {"source": a.source, "seed": a.seed, "peak": a.peak, "tone_bin": a.tone_bin, "tone_amp": a.tone_amp,
             "noise_rms": a.noise_rms, "noise_shift": a.noise_shift}
    paths = write_stimulus(codes, cfg, a.out_dir, name, extra)
    res = M.expected(codes, cfg)
    for p in paths:
        print("wrote", p)
    print(f"{cfg.nsamp} samples, {cfg.segments} segments, avg_shift {cfg.avg_shift}; predicted cycles to done "
          f"{M.expected_cycles(res, cfg)}; recommended simulation length {M.recommended_stop(res, cfg)} clock cycles "
          f"(sample period 1 in the model)")


if __name__ == "__main__":
    main()
