import numpy as np

from btpsd import rtlmodel as R
from btpsd import schedule as S
from btpsd import signals as sg
from btpsd.cordic import default_rom
from btpsd.estimator import HwConfig, psd_fixed
from btpsd.r8fft import r8fft_fixed, to_fixed

WL = 16


def _chunk(seed=0):
    rng = np.random.default_rng(seed)
    hi = (1 << (WL - 1)) - 1
    return rng.integers(-hi, hi, 512), rng.integers(-hi, hi, 512)


def test_row_bank_and_twiddle_functions_match_schedule():
    for p in range(3):
        adr = S.pass_addresses(p)
        tw = S.pass_twiddle_q(p)
        for b in range(64):
            hi, lo = b >> 3, b & 7
            s = (hi + lo) & 7
            for m in range(8):
                a = int(adr[b, m])
                assert R.row_fn(p, m, hi, lo) == int(S.row_of(a))
                assert (m + s) & 7 == int(S.bank_of(a))
            assert [R.twq(p, b, k) for k in range(8)] == tw[b].tolist()


def test_unload_address_map_is_the_inverse_of_output_index():
    for k in range(512):
        a2, a1, a0 = k & 7, (k >> 3) & 7, k >> 6
        a = 64 * a2 + 8 * a1 + a0
        assert int(S.output_index(a)) == k


def test_fft_rtl_matches_bit_accurate_model_and_schedule():
    re, im = _chunk(1)
    dut = R.Fft512Rtl(WL)
    dut.load(re, im)
    done = dut.run()
    yr, yi = dut.spectrum()
    er, ei = r8fft_fixed(re, im, default_rom(0.01), WL, "round")
    np.testing.assert_array_equal(yr, er)
    np.testing.assert_array_equal(yi, ei)
    assert done == S.schedule_cycles(R.LAT)["total_cycles"]
    assert R.stall_cycles(R.LAT) == S.schedule_cycles(R.LAT)["stalls_per_pass"][1]


def test_fft_rtl_needs_exactly_the_stall_it_has():
    re, im = _chunk(2)
    er, ei = r8fft_fixed(re, im, default_rom(0.01), WL, "round")
    st = R.stall_cycles(R.LAT)
    assert st > 0
    ok = R.Fft512Rtl(WL)
    ok.load(re, im)
    ok.run(stall=st)
    assert np.array_equal(ok.spectrum()[0], er)
    bad = R.Fft512Rtl(WL)
    bad.load(re, im)
    bad.run(stall=st - 1)
    assert not np.array_equal(bad.spectrum()[0], er)  # one cycle less reads a word before it is written


def test_psd_engine_loops_match_psd_fixed():
    cfg = HwConfig(wl=12, window="hann")
    x = sg.pink(4 * 512, 7, peak=0.5)
    codes, _ = to_fixed(x, cfg.wl, cfg.rounding)
    res = psd_fixed(codes, cfg, quantize=False)
    rom = default_rom(cfg.rom_tol)
    ch = codes.reshape(4, 512)
    Ar, Ai = r8fft_fixed(ch, np.zeros_like(ch), rom, cfg.wl, cfg.rounding)
    acc = np.zeros(1024, dtype=np.int64)
    for s in range(3):
        wr, wi = R.psd_segment_loops(Ar[s], Ai[s], Ar[s + 1], Ai[s + 1], cfg)
        np.testing.assert_array_equal(wr, res.X_re[s])
        np.testing.assert_array_equal(wi, res.X_im[s])
        acc += wr * wr + wi * wi
    np.testing.assert_array_equal(acc, res.acc)


def test_bispectrum_engine_loops_match_bispectrum_fixed():
    from btpsd import vectors as V
    from btpsd.bispectrum import BispecConfig
    from btpsd.fixedpoint import rshift
    from btpsd.reference import frame_signal

    v = V.bispec_vectors(wl=12, window="hann", max_bin=32)
    cfg = BispecConfig(wl=12, window="hann", max_bin=32)
    frames = frame_signal(v["x"].astype(np.int64), 512, cfg.hop)
    rom = default_rom(cfg.rom_tol)
    Fr, Fi = r8fft_fixed(frames, np.zeros_like(frames), rom, cfg.wl, cfg.rounding)
    accr = np.zeros((32, 32), dtype=np.int64)
    acci = np.zeros((32, 32), dtype=np.int64)
    for f in range(frames.shape[0]):
        tr, ti = R.bispectrum_frame_loops(Fr[f], Fi[f], cfg, 5)
        accr += tr
        acci += ti
    np.testing.assert_array_equal(rshift(accr, v["avg_shift"], cfg.rounding), v["B_re"])
    np.testing.assert_array_equal(rshift(acci, v["avg_shift"], cfg.rounding), v["B_im"])
