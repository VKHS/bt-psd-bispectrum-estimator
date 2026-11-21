"""Tests of the in-fabric test-signal generator: the model (btpsd.siggen), its Verilog (emulated cycle by cycle) and the generated table."""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python" / "scripts"))

from btpsd import siggen as G  # noqa: E402
from btpsd import vectors  # noqa: E402

import gen_mc  # noqa: E402

SRC = (ROOT / "modelcomposer" / "hdl" / "sig_gen.v").read_text()


# ---------------------------------------------------------------------------------------------------------------------
# table
# ---------------------------------------------------------------------------------------------------------------------
def test_rom_table_is_a_sine_of_the_documented_amplitude():
    for wl in (8, 12, 16):
        r = G.rom_table(wl)
        a = int(round(G.TONE_AMP * (1 << (wl - 1))))
        assert len(r) == 1024 and r[0] == 0 and r[512] == 0
        assert r[256] == a and r[768] == -a  # the extremes
        assert np.array_equal(r[1:512], -r[513:1024])  # odd symmetry (no ties for these amplitudes)
        assert np.array_equal(r[1:256], r[511:256:-1])  # mirror about the peak
        assert np.max(np.abs(r - G.TONE_AMP * (1 << (wl - 1)) * np.sin(2 * np.pi * np.arange(1024) / 1024))) <= 0.5 + 1e-9


def test_committed_rom_file_matches_the_table():
    text = (ROOT / "modelcomposer" / "hdl" / "sig_rom.v").read_text()
    assert text == gen_mc.sig_rom_text(12)
    ent = {int(k): int(v, 16) for k, v in re.findall(r"mem\[(\d+)\] = 12'h([0-9a-f]{3});", text)}
    assert sorted(ent) == list(range(1024))
    r = G.rom_table(12)
    assert all(((ent[k] - (int(r[k]) & 0xFFF)) == 0) for k in range(1024))
    assert "parameter WL = 12" in text and "always @(posedge clk) if (en) dout <= mem[addr];" in text


# ---------------------------------------------------------------------------------------------------------------------
# LFSR
# ---------------------------------------------------------------------------------------------------------------------
def _step_bitwise(s: int) -> int:
    """Galois step written differently from the model: shift out the LSB, then toggle the tap bits."""
    out = s & 1
    s >>= 1
    if out:
        for bit in (31, 21, 1, 0):  # taps 32, 22, 2, 1 of x^32 + x^22 + x^2 + x + 1
            s ^= 1 << bit
    return s


def test_lfsr_states_match_an_independent_bitwise_implementation():
    s = G.SEED
    got = G.lfsr_states(2000)
    for i in range(2000):
        assert int(got[i]) == s
        s = _step_bitwise(s)
    assert G.TAPS == (1 << 31) | (1 << 21) | 3


def test_lfsr_update_has_the_maximal_period_2_pow_32_minus_1():
    """The state map is a linear map over GF(2); its order is 2^32 - 1 iff the generator has the maximal period for every non-zero seed."""
    def apply(cols, v):
        out, i = 0, 0
        while v:
            if v & 1:
                out ^= cols[i]
            v >>= 1
            i += 1
        return out

    def mul(a, b):
        return [apply(a, c) for c in b]

    def power(m, e):
        r = [1 << i for i in range(32)]
        while e:
            if e & 1:
                r = mul(r, m)
            m = mul(m, m)
            e >>= 1
        return r

    t = [_step_bitwise(1 << i) for i in range(32)]
    ident = [1 << i for i in range(32)]
    n = (1 << 32) - 1
    assert power(t, n) == ident
    for q in (3, 5, 17, 257, 65537):  # prime factors of 2^32 - 1
        assert n % q == 0
        assert power(t, n // q) != ident


def test_lfsr_rejects_the_zero_seed_and_never_reaches_zero():
    with pytest.raises(ValueError):
        G.lfsr_states(10, 0)
    assert not np.any(G.lfsr_states(100000) == 0)
    assert len(set(G.lfsr_states(100000).tolist())) == 100000  # no repetition in 100 000 steps


# ---------------------------------------------------------------------------------------------------------------------
# the model against a literal emulation of sig_gen.v
# ---------------------------------------------------------------------------------------------------------------------
class SigGenRtl:
    """Cycle-by-cycle emulation of sig_gen.v (registers as in the Verilog, non-blocking updates)."""

    def __init__(self, wl=12, seed=G.SEED, amp=G.TONE_AMP):
        self.wl, self.seed = wl, seed
        self.rom = [int(v) & ((1 << wl) - 1) for v in G.rom_table(wl, amp)]
        self.phase, self.lfsr, self.v1, self.nz_q, self.x_r, self.v2 = 0, seed, 0, 0, 0, 0
        self.rom_q = 0

    @staticmethod
    def _signed(v, wl):
        return v - (1 << wl) if v >= (1 << (wl - 1)) else v

    def edge(self, rst_n, en, pinc, noise_en, noise_shift):
        wl = self.wl
        rom_next = self.rom[self.phase >> (G.PHW - 10)] if en else self.rom_q  # sig_rom: if (en) dout <= mem[addr]
        if not rst_n:
            self.phase, self.lfsr, self.v1, self.nz_q, self.x_r, self.v2 = 0, self.seed, 0, 0, 0, 0
        elif en:
            nz_raw = self._signed(self.lfsr >> (32 - wl), wl)
            nz_s = nz_raw >> noise_shift
            tone = self._signed(self.rom_q, wl)
            total = tone + self.nz_q
            x = max(-(1 << (wl - 1)), min((1 << (wl - 1)) - 1, total))
            self.nz_q = nz_s if noise_en else 0
            self.phase = (self.phase + pinc) & ((1 << G.PHW) - 1)
            self.lfsr = (self.lfsr >> 1) ^ (0x80200003 if self.lfsr & 1 else 0)
            self.x_r = x
            self.v2 = self.v1
            self.v1 = 1
        self.rom_q = rom_next

    @property
    def x(self):
        return self.x_r

    @property
    def valid(self):
        return self.v2


def _run_rtl(n, pinc, noise_shift, wl=12, gapped=False, noise_en=True):
    g = SigGenRtl(wl)
    g.edge(0, 0, pinc, noise_en, noise_shift or 0)
    g.edge(0, 0, pinc, noise_en, noise_shift or 0)
    out, ne, c = [], 0, 0
    while len(out) < n:
        en = (c % 3) != 1 if gapped else True
        g.edge(1, en, pinc, noise_en, noise_shift or 0)
        c += 1
        ne += 1 if en else 0
        assert g.valid == (1 if ne >= 2 else 0)
        if en and g.valid:
            out.append(g.x)
    return np.array(out)


@pytest.mark.parametrize("wl", [8, 12, 16])
@pytest.mark.parametrize("shift", [0, 2, 5])
def test_model_equals_the_rtl_emulation(wl, shift):
    if shift >= wl:
        pytest.skip("shift outside 0 .. wl-1")
    pinc = G.pinc_for(123.456)
    n = 600
    want = G.generate(n, pinc, shift, wl)
    assert np.array_equal(_run_rtl(n, pinc, shift, wl), want)
    assert np.array_equal(_run_rtl(n, pinc, shift, wl, gapped=True), want)  # holding while en = 0


def test_model_equals_the_rtl_emulation_for_random_increments_and_tone_only():
    rng = np.random.default_rng(5)
    for _ in range(6):
        pinc = int(rng.integers(0, 1 << 32))
        assert np.array_equal(_run_rtl(300, pinc, 3), G.generate(300, pinc, 3, 12))
        assert np.array_equal(_run_rtl(300, pinc, None, noise_en=False), G.generate(300, pinc, None, 12))


def test_saturation_is_exercised_and_exact():
    pinc = G.pinc_for(100.37)
    x = G.generate(2048, pinc, 0, 12)
    assert np.sum((x == 2047) | (x == -2048)) > 100
    assert np.array_equal(_run_rtl(2048, pinc, 0), x)


# ---------------------------------------------------------------------------------------------------------------------
# behaviour of the signal
# ---------------------------------------------------------------------------------------------------------------------
def test_tone_only_has_its_peak_at_the_requested_frequency():
    for b in (37.0, 100.37, 200.5, 311.25):
        x = G.generate(4096, G.pinc_for(b), None, 12).astype(float)
        spec = np.abs(np.fft.fft(x * np.hanning(4096)))[:2048]
        assert abs(np.argmax(spec) / 4.0 - b) <= 0.5  # 4096-point bins are 1/4 of a 1024-bin
        assert G.tone_bin(G.pinc_for(b)) == pytest.approx(b, abs=1e-6)
    assert G.pinc_for(0.0) == 0 and G.pinc_for(1024.0) == 0  # wraps


def test_noise_amplitude_follows_the_shift():
    base = G.generate(4096, 0, None, 12)  # pinc = 0: the table output stays at 0
    assert np.all(base == 0)
    for shift in (1, 3, 5):
        n = G.generate(4096, 0, shift, 12)
        assert abs(n).max() <= (1 << (11 - shift))
        assert abs(n).max() > (1 << (11 - shift)) // 2


def test_argument_checks():
    with pytest.raises(ValueError):
        G.generate(10, 1 << 32, None)
    with pytest.raises(ValueError):
        G.generate(10, 1, 12, 12)
    with pytest.raises(ValueError):
        G.generate(10, 1, -1, 12)


def test_verilog_constants_agree_with_the_model():
    assert "parameter [31:0] SEED = 32'hA5A51F3D" in SRC and G.SEED == 0xA5A51F3D
    assert "32'h80200003" in SRC and "parameter PHW  = 32" in SRC and G.PHW == 32
    assert "phase[PHW-1:PHW-10]" in SRC and G.ROM_BITS == 10
    assert "lfsr[31:32-WL]" in SRC and ">>> noise_shift" in SRC


# ---------------------------------------------------------------------------------------------------------------------
# vectors
# ---------------------------------------------------------------------------------------------------------------------
def test_vector_definitions():
    g = vectors.siggen_vectors(12)
    n = g["n"]
    a, c, d = g["x"][:n], g["x"][n:2 * n], g["x"][2 * n:]
    rec = G.selftest_record(12, 4)
    assert np.array_equal(a, rec["x"]) and g["pinc"] == rec["pinc"]
    assert g["sat_c"] > 100 and np.array_equal(d, G.generate(n, g["pinc"], None, 12))
    assert not np.array_equal(a, c) and not np.array_equal(a, d)
    t = vectors.selftest_vectors(12, 4)
    assert np.array_equal(t["x"], rec["x"]) and t["S"] == 3 and t["avg_shift"] == 2 and len(t["P"]) == 1024
    assert t["pinc"] == g["pinc"] and t["noise_shift"] == G.SELFTEST["noise_shift"]


def test_committed_vector_files_are_current(tmp_path):
    import gen_vectors

    gen_vectors.write_all(tmp_path)
    for name in ("siggen_wl12_x.mem", "siggen_wl12.vh", "selftest_wl12_P.mem", "selftest_wl12.vh"):
        assert (tmp_path / name).read_text() == (ROOT / "sim" / "vectors" / name).read_text(), name
