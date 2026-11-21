"""Executable desk-check of ``modelcomposer/hdl/bt_psd_mc.v``.

The wrapper's handshake with the estimator core (record, prime, replay, drain, done, re-arm) is emulated here statement by
statement from the Verilog, against a *protocol model* of ``rtl/bt_psd_top.v``: the model reproduces the top-level state machine
(T_INIT, T_CLR, T_LOAD, T_FFT, T_SEG), ``in_ready`` / ``busy`` and the durations measured in CI (``results/rtl_state_occupancy.csv``),
not the arithmetic.  What this proves: the pointer and prefetch logic hands every recorded sample to the core exactly once and in
order, ``done`` rises when the core is idle after the last segment, junk at the input is ignored, and the cycle count equals
``btpsd.mcmodel.expected_cycles``.  What it does not prove: that the Verilog text means what this emulation does; the
Icarus Verilog simulation of ``tb/tb_bt_psd_mc.v`` in CI decides that.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest

from btpsd import mcmodel as M
from btpsd import siggen
from btpsd.rtlmodel import psd_engine_state_cycles
from btpsd.squarer import mag_sq

ROOT = Path(__file__).resolve().parents[2]
WRAPPER = (ROOT / "modelcomposer" / "hdl" / "bt_psd_mc.v").read_text()
TOP = (ROOT / "rtl" / "bt_psd_top.v").read_text()

CAPTURE, PRIME, REPLAY, DRAIN, DONE = range(5)
INIT, CLR, LOAD, FFT, SEG = range(5)


class Core:
    """Protocol model of bt_psd_top: handshake and state durations only."""

    def __init__(self, nchunks, seg_cycles, fft=209, clr=1027):
        self.nchunks, self.seg_cycles, self.fft, self.clr = nchunks, list(seg_cycles), fft, clr
        self.accepted = []
        self.reset()

    def reset(self):
        self.st, self.cnt, self.ld_cnt, self.have_prev, self.nseg = INIT, 0, 0, False, 0

    @property
    def in_ready(self):
        return self.st == LOAD

    @property
    def busy(self):
        return self.st != LOAD

    def edge(self, in_valid, in_data):
        if self.st == INIT:
            self.st, self.cnt = CLR, 0
        elif self.st == CLR:
            self.cnt += 1
            if self.cnt == self.clr:
                self.st = LOAD
        elif self.st == LOAD:
            if in_valid:
                self.accepted.append(in_data)
                if self.ld_cnt == 511:
                    self.ld_cnt, self.st, self.cnt = 0, FFT, 0
                else:
                    self.ld_cnt += 1
        elif self.st == FFT:
            self.cnt += 1
            if self.cnt == self.fft:
                if self.have_prev:
                    self.st, self.cnt = SEG, 0
                else:
                    self.have_prev, self.st = True, LOAD
        else:  # SEG
            self.cnt += 1
            if self.cnt == self.seg_cycles[self.nseg]:
                self.nseg += 1
                self.st = LOAD


class Wrapper:
    """Statement-by-statement emulation of bt_psd_mc (parameter NSAMP = NCHUNKS * 512), power-up state of the Verilog."""

    def __init__(self, nchunks, seg_cycles, **core_kw):
        self.nsamp = nchunks * 512
        self.state, self.wr_cnt, self.rp = CAPTURE, 0, 0
        self.ram = [None] * self.nsamp
        self.rdata = None
        self.rst_sr = [0, 0, 0, 0]
        self.core = Core(nchunks, seg_cycles, **core_kw)

    def outputs(self):
        busy = (self.wr_cnt != 0) if self.state == CAPTURE else (self.state != DONE)
        return busy, self.state == DONE

    def step(self, rst, in_valid, in_data, ce=1):
        core = self.core
        core_in_valid = self.state == REPLAY
        accept = core_in_valid and core.in_ready
        rp_next = self.rp + 1 if accept else self.rp
        ra = 0 if rp_next >= self.nsamp else rp_next
        cap = self.state == CAPTURE and ce and in_valid
        core_busy = core.busy
        rdata_now = self.rdata
        rdata_next = self.ram[ra]  # reads the memory as it is before this edge
        # edge
        core_was_in_reset = not all(self.rst_sr)
        if cap:
            self.ram[self.wr_cnt] = in_data
        self.rdata = rdata_next
        self.rst_sr = self.rst_sr[1:] + [0 if rst else 1]
        if rst:
            self.state, self.wr_cnt, self.rp = CAPTURE, 0, 0
        else:
            self.rp = rp_next
            if self.state == CAPTURE:
                if cap:
                    self.wr_cnt += 1
                    if self.wr_cnt - 1 == self.nsamp - 1:
                        self.state = PRIME
            elif self.state == PRIME:
                self.state = REPLAY
            elif self.state == REPLAY:
                if accept and self.rp - (1 if accept else 0) == self.nsamp - 1:
                    self.state = DRAIN
            elif self.state == DRAIN:
                if not core_busy:
                    self.state = DONE
            else:
                self.state = DONE
        if core_was_in_reset:
            core.reset()  # asynchronous reset: the core stays in T_INIT while its reset is asserted
        else:
            core.edge(core_in_valid, rdata_now)
        if not all(self.rst_sr):
            core.reset()


def run(w, stream, tail=0, rst_at=()):
    """Drive ``stream`` = list of (valid, data) per cycle (then idle); returns (cycle of first done, per-cycle busy/done)."""
    done_cycle, hist = None, []
    n = len(stream) + tail
    for c in range(n):
        v, d = stream[c] if c < len(stream) else (0, 0)
        busy, done = w.outputs()
        hist.append((busy, done))
        if done and done_cycle is None:
            done_cycle = c
        w.step(1 if c in rst_at else 0, v, d)
    return done_cycle, hist


def seg_cycles_of(codes, cfg):
    res = M.expected(codes, cfg)
    fixed = sum(psd_engine_state_cycles(1, 0, clears=0).values())
    out = []
    for s in range(res.S):
        n = int(mag_sq(res.X_re[s], res.X_im[s], cfg.hw().squarer_units)[1].sum())
        out.append(fixed + n + M.T_SEG_EXTRA)
    return out


# ---------------------------------------------------------------------------------------------------------------------
# facts of the core that the emulation relies on, checked against the Verilog text
# ---------------------------------------------------------------------------------------------------------------------
def test_core_protocol_facts_hold_in_the_verilog():
    assert re.search(r"assign\s+in_ready\s*=\s*\(st\s*==\s*T_LOAD\)", TOP)
    assert re.search(r"assign\s+busy\s*=\s*\(st\s*!=\s*T_LOAD\)", TOP)
    assert "if (ld_cnt == 9'd511)" in TOP and "T_INIT: begin" in TOP
    assert re.search(r"T_CLR:\s*begin\s*if\s*\(eng_done\)\s*st\s*<=\s*T_LOAD;", TOP)
    # the samples are taken in T_LOAD whenever in_valid is high
    assert re.search(r"wire\s+ld_now\s*=\s*\(st\s*==\s*T_LOAD\)\s*&&\s*in_valid;", TOP)
    # the wrapper's own constants
    assert "S_CAPTURE = 3'd0, S_PRIME = 3'd1, S_REPLAY = 3'd2, S_DRAIN = 3'd3, S_DONE = 3'd4" in WRAPPER
    assert (CAPTURE, PRIME, REPLAY, DRAIN, DONE) == (0, 1, 2, 3, 4)
    assert "rd_addr" in WRAPPER and ".res_addr(rd_addr)" in WRAPPER and ".avg_shift(AVG_SHIFT)" in WRAPPER


# ---------------------------------------------------------------------------------------------------------------------
# behaviour (fast durations: the logic does not depend on them)
# ---------------------------------------------------------------------------------------------------------------------
FAST = dict(fft=21, clr=40)


def _record(n, seed=0):
    return [int(v) for v in np.random.default_rng(seed).integers(-2048, 2048, n)]


def test_every_sample_reaches_the_core_once_in_order_continuous_input():
    nch = 3
    rec = _record(nch * 512)
    w = Wrapper(nch, [57] * (nch - 1), **FAST)
    done_c, hist = run(w, [(1, v) for v in rec], tail=2500)
    assert w.core.accepted == rec
    assert done_c is not None and w.state == DONE
    assert hist[done_c][1] and all(h[1] for h in hist[done_c:])  # done stays high
    assert not any(h[1] for h in hist[:done_c])
    assert hist[0][0] is False and hist[1][0] is True  # busy from the first captured sample


def test_idle_cycles_and_junk_with_valid_low_are_not_recorded():
    nch = 3
    rec = _record(nch * 512, 1)
    stream = []
    for i, v in enumerate(rec):
        stream.append((0, -7 - i % 5))  # junk with valid low in front of every sample
        if i % 3 == 0:
            stream.append((0, 1234 % 2048))
        stream.append((1, v))
    stream += [(1, 99), (1, 98), (0, 5)] * 100  # samples offered after the record is complete: not captured
    w = Wrapper(nch, [57] * (nch - 1), **FAST)
    done_c, _ = run(w, stream, tail=2500)
    assert w.core.accepted == rec and done_c is not None


def test_ce_low_blocks_capture():
    nch = 2
    w = Wrapper(nch, [57], **FAST)
    for _ in range(50):
        w.step(0, 1, 5, ce=0)
    assert w.wr_cnt == 0 and w.state == CAPTURE
    rec = _record(nch * 512, 2)
    run(w, [(1, v) for v in rec], tail=1500)
    assert w.core.accepted == rec


def test_rearm_with_reset_pulse_during_replay_starts_a_clean_run():
    nch = 3
    first, second = _record(nch * 512, 3), _record(nch * 512, 4)
    w = Wrapper(nch, [57] * (nch - 1), **FAST)
    # let the first run get into the replay (more than half of the core's work), then reset
    done_c, _ = run(w, [(1, v) for v in first], tail=900)
    assert done_c is None and w.state == REPLAY  # not finished: the core is busy with the first segment
    w.step(1, 0, 0)  # one-cycle reset pulse
    assert w.state == CAPTURE and w.wr_cnt == 0 and w.rp == 0
    for _ in range(5):  # the core's reset is released three cycles after rst falls; idle meanwhile
        w.step(0, 0, 0)
    w.core.accepted.clear()
    done2, hist = run(w, [(1, v) for v in second], tail=2500)
    assert done2 is not None and w.core.accepted == second


def test_reset_after_done_rearms_and_second_result_comes_from_second_record():
    nch = 2
    first, second = _record(nch * 512, 5), _record(nch * 512, 6)
    w = Wrapper(nch, [57], **FAST)
    d1, _ = run(w, [(1, v) for v in first], tail=1500)
    assert d1 is not None and w.state == DONE and w.core.accepted == first
    w.step(1, 0, 0)
    for _ in range(5):
        w.step(0, 0, 0)
    w.core.accepted.clear()
    d2, _ = run(w, [(1, v) for v in second], tail=1500)
    assert d2 is not None and w.core.accepted == second


def test_the_core_may_be_slower_than_the_capture():
    """A long clear phase (core not ready when the replay starts) only delays the replay; nothing is lost."""
    nch = 2
    rec = _record(nch * 512, 7)
    w = Wrapper(nch, [57], fft=21, clr=3000)
    done_c, _ = run(w, [(1, v) for v in rec], tail=4000)
    assert w.core.accepted == rec and done_c is not None


def test_busy_is_high_between_the_first_sample_and_done_only():
    nch = 2
    rec = _record(nch * 512, 8)
    w = Wrapper(nch, [57], **FAST)
    done_c, hist = run(w, [(0, 0)] * 20 + [(1, v) for v in rec], tail=1500)
    assert not any(b for b, _ in hist[:21])  # idle before the first sample: busy low
    assert all(b for b, d in hist[21:done_c])
    assert not hist[done_c][0] and hist[done_c][1]


# ---------------------------------------------------------------------------------------------------------------------
# the cycle count with the measured durations equals the prediction of btpsd.mcmodel
# ---------------------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("gap_every", [0, 11])
def test_cycle_count_equals_the_model_prediction(gap_every):
    cfg = M.McConfig(wl=12, nchunks=3)
    codes = siggen.generate(cfg.nsamp, siggen.pinc_for(100.37), 3, 12)
    segs = seg_cycles_of(codes, cfg)
    w = Wrapper(cfg.nchunks, segs)  # fft = 209, clr = 1027
    stream, gaps = [], 0
    for i, v in enumerate(codes):
        if gap_every and i % gap_every == gap_every - 1:
            stream.append((0, 0))
            gaps += 1
        stream.append((1, int(v)))
    done_c, _ = run(w, stream, tail=sum(segs) + 4000)
    assert w.core.accepted == [int(v) for v in codes]
    assert done_c == M.expected_cycles(M.expected(codes, cfg), cfg) + gaps
