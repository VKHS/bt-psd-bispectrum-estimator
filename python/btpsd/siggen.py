"""Bit-exact model of ``modelcomposer/hdl/sig_gen.v``: a small in-fabric test-signal generator.

The generator is what the Simulink ``Sine Wave`` + ``Random Number`` sources and the AMD DDS Compiler provide in a Model
Composer design, written as a few lines of synthesisable Verilog so that the estimator can be exercised from inside the FPGA
(``bt_psd_selftest.v``).  It is not the AMD DDS Compiler IP and makes no claim about it.

Datapath, one sample per enabled clock cycle (sample ``n`` = 0, 1, ...)::

    phase_n  = (n * pinc) mod 2**PHW                       phase accumulator, value before the increment
    tone_n   = ROM[phase_n >> (PHW - 10)]                  1024-entry sine table, phase truncated, no interpolation
    state_n  = LFSR state after n updates from SEED        32-bit Galois LFSR, taps 0x80200003, shifted right
    noise_n  = (signed top WL bits of state_n) >> noise_shift     arithmetic shift, 0 when the noise is disabled
    x_n      = saturate(tone_n + noise_n)                  to WL bits, two's complement

``ROM[k] = clip(floor(TONE_AMP * 2**(WL-1) * sin(2 pi k / 1024) + 1/2))``.  All of it is integer arithmetic except the table,
whose entries are computed here and written into ``sig_rom.v`` by ``python/scripts/gen_mc.py`` (the test compares them).
"""
from __future__ import annotations

import numpy as np

__all__ = ["TONE_AMP", "ROM_BITS", "PHW", "SEED", "TAPS", "rom_table", "lfsr_states", "generate", "pinc_for", "tone_bin",
           "selftest_record"]

TONE_AMP = 0.5  # sine table amplitude, fraction of full scale
ROM_BITS = 10  # table address bits (1024 entries)
PHW = 32  # phase accumulator width
SEED = 0xA5A51F3D  # LFSR start state (non-zero)
TAPS = 0x80200003  # Galois feedback mask of x^32 + x^22 + x^2 + x + 1 (right-shifting form)


def rom_table(wl: int = 12, amp: float = TONE_AMP) -> np.ndarray:
    """The 1024-entry sine table as signed integers."""
    k = np.arange(1 << ROM_BITS)
    v = np.floor(amp * (1 << (wl - 1)) * np.sin(2.0 * np.pi * k / (1 << ROM_BITS)) + 0.5)
    return np.clip(v, -(1 << (wl - 1)), (1 << (wl - 1)) - 1).astype(np.int64)


def lfsr_states(n: int, seed: int = SEED) -> np.ndarray:
    """States ``state_0 = seed, state_1, ...`` (``n`` values, uint64 holding 32-bit words)."""
    out = np.empty(n, dtype=np.uint64)
    s = int(seed) & 0xFFFFFFFF
    if s == 0:
        raise ValueError("an LFSR seed of zero never leaves the all-zero state")
    for i in range(n):
        out[i] = s
        s = (s >> 1) ^ (TAPS if (s & 1) else 0)
    return out


def pinc_for(bin_1024: float, phw: int = PHW) -> int:
    """Phase increment of a tone at ``bin_1024`` (bins of a 1024-point spectrum, i.e. cycles per sample = bin / 1024)."""
    return int(round((bin_1024 / 1024.0) * (1 << phw))) % (1 << phw)


def tone_bin(pinc: int, phw: int = PHW) -> float:
    """Inverse of :func:`pinc_for`."""
    return (pinc / (1 << phw)) * 1024.0


def generate(n: int, pinc: int, noise_shift: int | None = None, wl: int = 12, seed: int = SEED, amp: float = TONE_AMP) -> np.ndarray:
    """First ``n`` samples of the generator as signed integer codes; ``noise_shift = None`` switches the noise off."""
    if not 0 <= pinc < (1 << PHW):
        raise ValueError("pinc must fit the phase accumulator")
    rom = rom_table(wl, amp)
    phase = (np.arange(n, dtype=np.uint64) * np.uint64(pinc)) & np.uint64((1 << PHW) - 1)
    tone = rom[(phase >> np.uint64(PHW - ROM_BITS)).astype(np.int64)]
    if noise_shift is None:
        nz = np.zeros(n, dtype=np.int64)
    else:
        if not 0 <= noise_shift < wl:
            raise ValueError("noise_shift must be in 0 .. wl-1")
        top = (lfsr_states(n, seed) >> np.uint64(32 - wl)).astype(np.int64)
        top = np.where(top >= (1 << (wl - 1)), top - (1 << wl), top)
        nz = top >> noise_shift  # arithmetic shift
    return np.clip(tone + nz, -(1 << (wl - 1)), (1 << (wl - 1)) - 1).astype(np.int64)


# the record used by the self-test testbench and by the golden vectors
SELFTEST = dict(wl=12, nchunks=4, bin_1024=100.37, noise_shift=3)


def selftest_record(wl: int = 12, nchunks: int = 4, bin_1024: float = SELFTEST["bin_1024"],
                    noise_shift: int | None = SELFTEST["noise_shift"]) -> dict:
    """The generator's first ``nchunks * 512`` samples and the parameters that produce them."""
    pinc = pinc_for(bin_1024)
    return dict(wl=wl, nchunks=nchunks, pinc=pinc, noise_shift=noise_shift, x=generate(nchunks * 512, pinc, noise_shift, wl))
