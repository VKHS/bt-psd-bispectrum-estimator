"""Model Composer flow: expected results, cycle prediction and decoding of the logs of a Simulink run.

The Black Box ``bt_psd_mc`` (``modelcomposer/hdl/bt_psd_mc.v``) records the first ``nchunks * 512`` samples that arrive
with ``in_valid = 1``, replays them into the estimator and raises ``done``.  Everything a user of the Simulink flow needs
from Python is here:

* :class:`McConfig`          the parameters of the generated block (word length, chunks, window, rounding);
* :func:`gateway_codes`      what a Gateway In block does to a double-precision source (the three quantisation modes);
* :func:`expected`           the bit-exact spectrum the block must produce for a record of integer codes;
* :func:`expected_cycles`    a prediction of the clock cycles from the first sample to ``done`` (a *prediction*, see below);
* :func:`recommended_stop`   a Simulink stop time that leaves the read-out window after ``done``;
* :func:`decode_spectrum`    the 1024 result words out of the logged ``rd_addr`` / ``rd_data`` / ``done`` signals;
* :func:`check_logs`         the whole comparison of a logged run against the model.

Cycle prediction.  From the first captured sample to the first cycle with ``done = 1`` (continuous ``in_valid``)::

    N + C * (512 + T_FFT) + S * (T_SEG_FIXED + T_SEG_EXTRA) + n_sq + 2

with ``N = 512 C`` the record, ``C`` chunks, ``S = C - 1`` segments, ``T_FFT = 209`` cycles per chunk and
``T_SEG_FIXED = 31 233`` the data-independent cycles of the segment engine per segment
(:func:`btpsd.rtlmodel.psd_engine_state_cycles`), ``T_SEG_EXTRA = 2`` and ``n_sq`` the data-dependent cycles of the
sum-of-odds squarer (``PsdResult.squarer_cycles``).  ``T_FFT`` and ``T_SEG_EXTRA`` are read off the CI measurement of
``tb_psd_top`` (``results/rtl_state_occupancy.csv``: 836 cycles in ``T_FFT`` for 4 chunks, 115 462 cycles in ``T_SEG``
for 3 segments against 115 456 in the engine states); the handshake of the wrapper adds the final ``+ 2``.  The testbench
``tb/tb_bt_psd_mc.v`` prints the cycle count that the simulation really takes (``MC,run,0,cycles_from_start_to_done``).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .estimator import HwConfig, PsdResult, hw_power_scale, psd_exact_welch, psd_fixed
from .rtlmodel import psd_engine_state_cycles

__all__ = ["CHUNK", "BINS", "T_FFT_PER_CHUNK", "T_SEG_EXTRA", "McConfig", "gateway_codes", "codes_to_gateway_values", "expected",
           "expected_cycles", "recommended_stop", "decode_spectrum", "check_logs", "RunCheck", "captured_record"]

CHUNK = 512
BINS = 1024
T_FFT_PER_CHUNK = 209  # cycles of the top-level state T_FFT per chunk (836 for 4 chunks in the CI run of tb_psd_top)
T_SEG_EXTRA = 2  # cycles per segment in T_SEG beyond the engine's busy states (115 462 - 115 456 = 6 for 3 segments)
READOUT_MARGIN = 2048  # cycles kept after the predicted done: a free-running 10-bit address counter needs 1 025


@dataclass(frozen=True)
class McConfig:
    """Parameters of one generated Black Box (``python/scripts/gen_mc.py``)."""

    wl: int = 12
    nchunks: int = 5
    round: int = 1  # RTL parameter ROUND: 1 rounding to nearest (half up), 0 truncation
    window: str = "hann"  # the wrapper uses the window defaults of bt_psd_top (Hann)

    def __post_init__(self):
        if self.nchunks < 2:
            raise ValueError("nchunks must be at least 2")
        if self.wl < 8:
            raise ValueError("wl must be at least 8")
        if self.round not in (0, 1):
            raise ValueError("round must be 0 or 1")

    @property
    def nsamp(self) -> int:
        return self.nchunks * CHUNK

    @property
    def segments(self) -> int:
        return self.nchunks - 1

    @property
    def avg_shift(self) -> int:
        """``ceil(log2(segments))``, as ``bt_psd_mc`` computes it."""
        return (self.segments - 1).bit_length()

    @property
    def accw(self) -> int:
        return 2 * self.wl + 8

    @property
    def scale(self) -> float:
        """``P_welch = rd_data * scale`` (``btpsd.reference.welch_twosided(scaling='paper')`` scale)."""
        return hw_power_scale(self.wl, self.segments, self.avg_shift)

    def hw(self) -> HwConfig:
        return HwConfig(wl=self.wl, window=self.window, rounding="round" if self.round else "trunc")


# ---------------------------------------------------------------------------
# input side
# ---------------------------------------------------------------------------
def gateway_codes(x, wl: int = 12, mode: str = "round_inf") -> np.ndarray:
    """Integer codes a Gateway In block of type ``Fix_<wl>_<wl-1>`` produces from double-precision samples ``x``.

    ``mode``: ``"round_inf"`` Round (unbiased: +/- Inf), i.e. half away from zero; ``"round_even"`` Round (unbiased: Even
    Values); ``"trunc"`` Truncate (towards minus infinity).  Overflow saturates.  Samples that are exact multiples of
    ``2**-(wl-1)`` are unaffected by the mode, which is why the stimulus files of this repository hold exact codes.
    """
    v = np.asarray(x, dtype=float) * float(1 << (wl - 1))
    if mode == "round_inf":
        q = np.sign(v) * np.floor(np.abs(v) + 0.5)
    elif mode == "round_even":
        q = np.round(v)  # numpy rounds half to even
    elif mode == "trunc":
        q = np.floor(v)
    else:
        raise ValueError(f"unknown quantisation mode {mode!r}")
    return np.clip(q, -(1 << (wl - 1)), (1 << (wl - 1)) - 1).astype(np.int64)


def codes_to_gateway_values(codes, wl: int = 12) -> np.ndarray:
    """Doubles that a Gateway In quantises back to exactly ``codes`` (``code / 2**(wl-1)``, exact in binary floating point)."""
    return np.asarray(codes, dtype=np.int64).astype(float) / float(1 << (wl - 1))


def captured_record(codes, valid=None, nsamp: int | None = None) -> np.ndarray:
    """The samples ``bt_psd_mc`` records: the first ``nsamp`` codes whose ``in_valid`` is 1 (all of them when ``valid`` is None)."""
    c = np.asarray(codes).astype(np.int64).ravel()
    if valid is not None:
        v = np.asarray(valid).ravel() != 0
        if v.shape != c.shape:
            raise ValueError("codes and valid must have the same length")
        c = c[v]
    if nsamp is not None:
        if len(c) < nsamp:
            raise ValueError(f"only {len(c)} valid samples, the block records {nsamp}")
        c = c[:nsamp]
    return c


# ---------------------------------------------------------------------------
# expected result and timing
# ---------------------------------------------------------------------------
def expected(codes, cfg: McConfig | None = None) -> PsdResult:
    """Bit-exact result of the block for a record of integer codes (the first ``cfg.nsamp`` are used)."""
    cfg = cfg or McConfig()
    c = np.asarray(codes).astype(np.int64).ravel()
    if len(c) < cfg.nsamp:
        raise ValueError(f"need {cfg.nsamp} samples, got {len(c)}")
    lim = 1 << (cfg.wl - 1)
    if c.min() < -lim or c.max() > lim - 1:
        raise ValueError(f"codes outside the {cfg.wl}-bit range")
    return psd_fixed(c[: cfg.nsamp], cfg.hw(), quantize=False)


def expected_cycles(res: PsdResult, cfg: McConfig | None = None) -> int:
    """Predicted clock cycles from the first captured sample to the first cycle with ``done = 1`` (see the module docstring)."""
    cfg = cfg or McConfig()
    seg = sum(psd_engine_state_cycles(res.S, res.squarer_cycles, clears=0).values())  # busy states of the engine, no clear
    return cfg.nsamp + cfg.nchunks * (CHUNK + T_FFT_PER_CHUNK) + seg + res.S * T_SEG_EXTRA + 2


def recommended_stop(res: PsdResult, cfg: McConfig | None = None, margin: int = READOUT_MARGIN) -> int:
    """Simulation length in clock cycles that covers the predicted ``done`` plus 2 % and a read-out window of ``margin`` cycles."""
    p = expected_cycles(res, cfg)
    return int(np.ceil(p * 1.02)) + int(margin)


# ---------------------------------------------------------------------------
# output side
# ---------------------------------------------------------------------------
def decode_spectrum(addr, data, done, lag: int = 1, bins: int = BINS) -> np.ndarray:
    """The ``bins`` result words from per-clock logs of ``rd_addr``, ``rd_data`` and ``done``.

    ``rd_data`` follows ``rd_addr`` by ``lag`` clock cycles (1 for ``bt_psd_mc``).  A data sample counts when ``done`` is high in
    its cycle and in the cycle of its address.  Raises :class:`ValueError` when a bin is never read or two reads of one bin differ.
    """
    a = np.asarray(addr, dtype=float).ravel()
    d = np.asarray(data, dtype=float).ravel()
    n = np.asarray(done, dtype=float).ravel() != 0
    if not (len(a) == len(d) == len(n)):
        raise ValueError("rd_addr, rd_data and done must have the same length")
    if len(a) <= lag:
        raise ValueError("log too short")
    t = np.arange(lag, len(a))
    ok = n[t] & n[t - lag]
    if not ok.any():
        raise ValueError("done was never high long enough to read a result (run the simulation longer)")
    ai = np.rint(a[t[ok] - lag]).astype(np.int64)
    di = np.rint(d[t[ok]]).astype(np.int64)
    if ai.min() < 0 or ai.max() >= bins:
        raise ValueError("rd_addr outside 0 .. %d" % (bins - 1))
    out = np.full(bins, -1, dtype=np.int64)
    seen = np.zeros(bins, dtype=bool)
    for b, v in zip(ai, di):
        if seen[b] and out[b] != v:
            raise ValueError(f"two reads of bin {b} differ ({out[b]} and {v}): lag {lag} may be wrong or the result changed")
        out[b] = v
        seen[b] = True
    if not seen.all():
        miss = np.flatnonzero(~seen)
        raise ValueError(f"{len(miss)} bins never read (first: {miss[:8].tolist()}); keep the simulation running at least "
                         f"{bins + lag + 1} cycles after done")
    return out


@dataclass
class RunCheck:
    ok: bool
    n_bins: int
    n_mismatch: int
    mismatches: list = field(default_factory=list)  # (bin, got, expected), first few
    done_cycle: int | None = None  # index of the first cycle with done = 1 in the log
    first_valid_cycle: int | None = None
    cycles_to_done: int | None = None
    predicted_cycles: int | None = None
    lag: int = 1
    snr_db_vs_exact: float | None = None
    scale: float = 1.0
    peak_bin: int | None = None
    message: str = ""

    def summary(self) -> str:
        head = "PASS" if self.ok else "FAIL"
        t = f"{head}: {self.n_bins - self.n_mismatch} of {self.n_bins} bins equal the bit-exact model"
        if self.cycles_to_done is not None and self.predicted_cycles is not None:
            t += f"; done after {self.cycles_to_done} cycles (predicted {self.predicted_cycles})"
        if self.snr_db_vs_exact is not None:
            t += f"; signal-to-error ratio {self.snr_db_vs_exact:.1f} dB against the float Welch estimate"
        return t + (f" [{self.message}]" if self.message else "")


def _snr_db(est, ref) -> float:
    est = np.asarray(est, dtype=float)
    ref = np.asarray(ref, dtype=float)
    e = np.sum((est - ref) ** 2)
    return float("inf") if e == 0 else 10.0 * np.log10(np.sum(ref**2) / e)


def check_logs(in_code, in_valid, rd_addr, rd_data, done, cfg: McConfig | None = None, lag: int = 1) -> RunCheck:
    """Compare a logged Simulink run against the bit-exact model.

    The logs are per clock cycle (one entry per sample period of the design): ``in_code`` the integer code at the output of
    the Gateway In (``double value * 2**(wl-1)``), ``in_valid`` the valid input of the block, and ``rd_addr`` / ``rd_data`` /
    ``done`` its address input and outputs.  The expected spectrum is computed from the *logged input codes*, so it does not
    depend on which Simulink source produced them or how the gateway rounded them.
    """
    cfg = cfg or McConfig()
    codes = np.rint(np.asarray(in_code, dtype=float)).astype(np.int64).ravel()
    valid = np.asarray(in_valid).ravel()
    rec = captured_record(codes, valid, cfg.nsamp)
    res = expected(rec, cfg)
    pred = expected_cycles(res, cfg)
    dn = np.asarray(done).ravel() != 0
    vi = np.flatnonzero(np.asarray(valid) != 0)
    first_valid = int(vi[0]) if len(vi) else None
    done_cycle = int(np.flatnonzero(dn)[0]) if dn.any() else None
    rc = RunCheck(ok=False, n_bins=BINS, n_mismatch=BINS, done_cycle=done_cycle, first_valid_cycle=first_valid, predicted_cycles=pred,
                  lag=lag, scale=cfg.scale)
    if done_cycle is not None and first_valid is not None:
        rc.cycles_to_done = done_cycle - first_valid
    try:
        got = decode_spectrum(rd_addr, rd_data, done, lag=lag)
    except ValueError as e:
        rc.message = str(e)
        return rc
    want = res.P_int.astype(np.int64)
    bad = np.flatnonzero(got != want)
    rc.n_mismatch = int(len(bad))
    rc.mismatches = [(int(b), int(got[b]), int(want[b])) for b in bad[:10]]
    rc.ok = rc.n_mismatch == 0
    rc.peak_bin = int(np.argmax(got[: BINS // 2]))
    try:
        ref = psd_exact_welch(rec.astype(float) / float(1 << (cfg.wl - 1)), cfg.window)
        rc.snr_db_vs_exact = _snr_db(got.astype(float) * cfg.scale, ref)
    except Exception as e:  # information only
        rc.message = f"float comparison unavailable: {e}"
    return rc
