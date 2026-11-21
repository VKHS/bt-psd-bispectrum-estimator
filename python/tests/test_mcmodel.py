"""Tests of the Model Composer helpers: expected result, cycle prediction, log decoding and the two command-line tools."""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat, savemat

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python" / "scripts"))

from btpsd import mcmodel as M  # noqa: E402
from btpsd import siggen, vectors  # noqa: E402
from btpsd.estimator import HwConfig, hw_power_scale, psd_fixed  # noqa: E402

import compare_mc_result as CMP  # noqa: E402
import export_mc_stimulus as EXP  # noqa: E402


# ---------------------------------------------------------------------------------------------------------------------
# configuration and quantiser
# ---------------------------------------------------------------------------------------------------------------------
def test_config_derived_values_match_the_estimator():
    for nch in range(2, 12):
        cfg = M.McConfig(wl=12, nchunks=nch)
        x = np.random.default_rng(nch).integers(-100, 100, nch * 512)
        res = psd_fixed(x, cfg.hw(), quantize=False)
        assert cfg.segments == res.S
        assert cfg.avg_shift == res.avg_shift
        assert cfg.nsamp == nch * 512
        assert cfg.accw == 2 * cfg.wl + 8
        assert cfg.scale == pytest.approx(hw_power_scale(12, res.S, res.avg_shift))


def test_config_rejects_bad_parameters():
    for kw in (dict(nchunks=1), dict(wl=4), dict(round=2)):
        with pytest.raises(ValueError):
            M.McConfig(**kw)


def test_defaults_match_the_generator_defaults():
    import gen_mc

    d = M.McConfig()
    assert (d.wl, d.nchunks, d.round) == (gen_mc.DEFAULTS["wl"], gen_mc.DEFAULTS["nchunks"], gen_mc.DEFAULTS["rounding"])


def test_gateway_quantiser_modes_on_ties_and_saturation():
    lsb = 2.0 ** -11  # Fix_12_11
    ties = np.array([0.5, 1.5, 2.5, -0.5, -1.5, -2.5]) * lsb
    assert M.gateway_codes(ties, 12, "round_inf").tolist() == [1, 2, 3, -1, -2, -3]
    assert M.gateway_codes(ties, 12, "round_even").tolist() == [0, 2, 2, 0, -2, -2]
    assert M.gateway_codes(ties, 12, "trunc").tolist() == [0, 1, 2, -1, -2, -3]
    assert M.gateway_codes([1.0, -1.0, 5.0, -5.0], 12, "round_inf").tolist() == [2047, -2048, 2047, -2048]
    with pytest.raises(ValueError):
        M.gateway_codes([0.0], 12, "nearest")


def test_exact_codes_survive_every_gateway_mode():
    codes = np.arange(-2048, 2048)
    vals = M.codes_to_gateway_values(codes, 12)
    for mode in ("round_inf", "round_even", "trunc"):
        assert np.array_equal(M.gateway_codes(vals, 12, mode), codes)


def test_captured_record_follows_the_valid_flag():
    codes = np.arange(20)
    valid = np.array([0, 0, 1, 1, 0, 1] + [1] * 14)
    assert M.captured_record(codes, valid, 5).tolist() == [2, 3, 5, 6, 7]
    assert M.captured_record(codes).tolist() == codes.tolist()
    with pytest.raises(ValueError):
        M.captured_record(codes, valid, 40)
    with pytest.raises(ValueError):
        M.captured_record(codes, valid[:-1], 5)


# ---------------------------------------------------------------------------------------------------------------------
# expected result and cycles
# ---------------------------------------------------------------------------------------------------------------------
def test_expected_equals_psd_fixed_and_checks_its_input():
    p = vectors.psd_vectors(12)
    cfg = M.McConfig(wl=12, nchunks=4)
    res = M.expected(p["x"], cfg)
    assert np.array_equal(res.P_int, p["P"])  # the golden vectors of the PSD testbench
    assert M.expected(np.concatenate([p["x"], [1, 2, 3]]), cfg).P_int.tolist() == p["P"].tolist()  # extra samples are ignored
    with pytest.raises(ValueError):
        M.expected(p["x"][:-1], cfg)
    bad = p["x"].copy()
    bad[7] = 4096
    with pytest.raises(ValueError):
        M.expected(bad, cfg)


def test_truncating_block_uses_truncation_in_the_model():
    x = np.random.default_rng(3).integers(-1500, 1500, 4 * 512)
    r1 = M.expected(x, M.McConfig(nchunks=4, round=1)).P_int
    r0 = M.expected(x, M.McConfig(nchunks=4, round=0)).P_int
    assert not np.array_equal(r1, r0)
    assert np.array_equal(r0, psd_fixed(x, HwConfig(wl=12, rounding="trunc"), quantize=False).P_int)


def test_cycle_constants_agree_with_the_ci_measurement_of_the_psd_testbench():
    """results/rtl_state_occupancy.csv holds the measured occupancy of tb_psd_top (4 chunks, 3 segments, WL 12)."""
    occ = {}
    with open(ROOT / "results" / "rtl_state_occupancy.csv") as f:
        for r in csv.DictReader(f):
            if r["testbench"] == "psd_top":
                occ[r["state"]] = int(r["cycles"])
    p = vectors.psd_vectors(12)
    assert occ["T_FFT"] == 4 * M.T_FFT_PER_CHUNK
    assert occ["T_LOAD"] == 4 * M.CHUNK
    from btpsd.rtlmodel import psd_engine_state_cycles

    eng = sum(psd_engine_state_cycles(3, p["squarer_cycles"], clears=0).values())
    assert occ["T_SEG"] == eng + 3 * M.T_SEG_EXTRA


def test_cycle_prediction_is_composed_of_the_documented_terms():
    cfg = M.McConfig(wl=12, nchunks=4)
    rec = siggen.selftest_record(12, 4)
    res = M.expected(rec["x"], cfg)
    want = 2048 + 4 * (512 + 209) + 3 * 31233 + res.squarer_cycles + 3 * 2 + 2
    assert M.expected_cycles(res, cfg) == want
    assert M.recommended_stop(res, cfg) >= int(1.02 * want) + 2048
    assert M.recommended_stop(res, cfg) - want > 1025 + 1000  # a free-running address counter has room to cover all 1024 bins


def test_more_segments_take_proportionally_longer():
    c4, c8 = M.McConfig(nchunks=5), M.McConfig(nchunks=9)
    x = np.random.default_rng(1).integers(-300, 300, 9 * 512)
    r4, r8 = M.expected(x, c4), M.expected(x, c8)
    t4, t8 = M.expected_cycles(r4, c4), M.expected_cycles(r8, c8)
    assert t8 > t4 and (t8 - t4) > 4 * 31233


# ---------------------------------------------------------------------------------------------------------------------
# synthetic Simulink logs
# ---------------------------------------------------------------------------------------------------------------------
def synth_log(cfg, codes, idle_before=3, gap_every=0, lag=1, tail=2200, junk=True, readout="free"):
    """Logs in the format of bt_mc_run.m for a block that behaves as bt_psd_mc is meant to (timing from the model)."""
    rng = np.random.default_rng(0)
    in_code, in_valid = [], []
    for _ in range(idle_before):
        in_code.append(int(rng.integers(-5, 5)))
        in_valid.append(0)
    gaps = 0
    for i, c in enumerate(codes[: cfg.nsamp]):
        if gap_every and i % gap_every == gap_every - 1:
            in_code.append(int(rng.integers(-5, 5)))
            in_valid.append(0)
            gaps += 1
        in_code.append(int(c))
        in_valid.append(1)
    res = M.expected(codes, cfg)
    done_at = idle_before + M.expected_cycles(res, cfg) + gaps
    total = done_at + tail
    pad = total - len(in_code)
    in_code += [int(rng.integers(-5, 5)) for _ in range(pad)]
    in_valid += [1 if junk else 0] * pad  # valid samples after the record are not captured
    t = np.arange(total)
    addr = (t % 1024) if readout == "free" else np.zeros(total, dtype=int)
    done = (t >= done_at).astype(int)
    data = np.zeros(total, dtype=np.int64)
    for k in range(lag, total):
        data[k] = res.P_int[addr[k - lag]] if done[k] else int(rng.integers(0, 1000))
    return dict(in_code=np.array(in_code, dtype=float), in_valid=np.array(in_valid, dtype=float), rd_addr=addr.astype(float),
                rd_data=data.astype(float), done=done.astype(float)), res


@pytest.fixture(scope="module")
def cfg4():
    return M.McConfig(wl=12, nchunks=4)


@pytest.fixture(scope="module")
def codes4():
    return vectors.psd_vectors(12)["x"]


def test_decode_spectrum_roundtrip(cfg4, codes4):
    log, res = synth_log(cfg4, codes4)
    got = M.decode_spectrum(log["rd_addr"], log["rd_data"], log["done"])
    assert np.array_equal(got, res.P_int)


def test_decode_spectrum_failure_modes(cfg4, codes4):
    log, _ = synth_log(cfg4, codes4)
    never = np.zeros_like(log["done"])
    with pytest.raises(ValueError, match="never high"):
        M.decode_spectrum(log["rd_addr"], log["rd_data"], never)
    short = log["done"].copy()
    d0 = int(np.flatnonzero(short)[0])
    short[d0 + 300:] = 0  # only a part of the address range is read while done is high
    with pytest.raises(ValueError, match="never read"):
        M.decode_spectrum(log["rd_addr"], log["rd_data"], short)
    bad = log["rd_data"].copy()
    bad[d0 + 600] += 1
    bad[d0 + 600 + 1024] += 2  # the same bin twice with different words
    with pytest.raises(ValueError, match="differ"):
        M.decode_spectrum(log["rd_addr"], bad, log["done"])
    with pytest.raises(ValueError):
        M.decode_spectrum(log["rd_addr"][:-1], log["rd_data"], log["done"])
    with pytest.raises(ValueError):
        M.decode_spectrum([2000, 2000, 2000], [1, 1, 1], [1, 1, 1])


def test_check_logs_passes_a_correct_run_with_idle_cycles_and_gaps(cfg4, codes4):
    for idle, gap in ((0, 0), (3, 0), (5, 7)):
        log, _ = synth_log(cfg4, codes4, idle_before=idle, gap_every=gap)
        rc = M.check_logs(log["in_code"], log["in_valid"], log["rd_addr"], log["rd_data"], log["done"], cfg4)
        assert rc.ok, rc.summary()
        assert rc.n_mismatch == 0 and rc.mismatches == []
        assert rc.first_valid_cycle == idle
        assert rc.cycles_to_done == rc.predicted_cycles + (len(codes4) // gap if gap else 0)
        assert rc.snr_db_vs_exact is not None and rc.snr_db_vs_exact > 20
        assert "PASS" in rc.summary()


def test_check_logs_reports_a_wrong_bin_and_a_wrong_latency(cfg4, codes4):
    log, res = synth_log(cfg4, codes4)
    d0 = int(np.flatnonzero(log["done"])[0])
    data = log["rd_data"].copy()
    data[d0 + 100:] += (log["rd_addr"][d0 + 99:-1] == 17)[: len(data) - d0 - 100]  # bin 17 off by one in every read
    rc = M.check_logs(log["in_code"], log["in_valid"], log["rd_addr"], data, log["done"], cfg4)
    assert not rc.ok and rc.n_mismatch >= 1
    assert rc.mismatches[0][0] == 17 and rc.mismatches[0][1] == rc.mismatches[0][2] + 1
    assert "FAIL" in rc.summary()
    # a design with an extra register in the data path reads as a latency error, not as garbage
    log2, _ = synth_log(cfg4, codes4, lag=2)
    rc2 = M.check_logs(log2["in_code"], log2["in_valid"], log2["rd_addr"], log2["rd_data"], log2["done"], cfg4, lag=1)
    assert not rc2.ok
    rc3 = M.check_logs(log2["in_code"], log2["in_valid"], log2["rd_addr"], log2["rd_data"], log2["done"], cfg4, lag=2)
    assert rc3.ok


def test_check_logs_without_enough_valid_samples_is_an_error(cfg4, codes4):
    log, _ = synth_log(cfg4, codes4)
    with pytest.raises(ValueError, match="valid samples"):
        M.check_logs(log["in_code"][:1000], log["in_valid"][:1000], log["rd_addr"][:1000], log["rd_data"][:1000], log["done"][:1000], cfg4)


def test_check_logs_when_done_never_rises(cfg4, codes4):
    log, _ = synth_log(cfg4, codes4)
    rc = M.check_logs(log["in_code"], log["in_valid"], log["rd_addr"], log["rd_data"], np.zeros_like(log["done"]), cfg4)
    assert not rc.ok and "never high" in rc.message and rc.done_cycle is None


# ---------------------------------------------------------------------------------------------------------------------
# stimulus export
# ---------------------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("source", EXP.SOURCES)
def test_export_every_source(tmp_path, source):
    cfg = M.McConfig(wl=12, nchunks=3)
    codes = EXP.make_record(source, cfg, seed=2)
    assert codes.shape == (1536,) and codes.dtype == np.int64
    assert codes.min() >= -2048 and codes.max() <= 2047
    assert np.abs(codes).max() > 50  # not silent
    paths = EXP.write_stimulus(codes, cfg, tmp_path, "s", {"source": source})
    assert [p.name for p in paths] == ["s.mat", "s.csv", "s.json"]
    d = loadmat(paths[0])
    st = d["bt_stim"]
    assert st.shape == (1536, 2) and np.array_equal(st[:, 0], np.arange(1536))
    assert np.array_equal(np.rint(st[:, 1] * 2048).astype(int), codes)  # exact
    assert np.array_equal(st[:, 1] * 2048, np.rint(st[:, 1] * 2048))
    assert d["bt_codes"].ravel().tolist() == codes.tolist()
    assert d["bt_wl"].item() == 12 and d["bt_nchunks"].item() == 3 and d["bt_nsamp"].item() == 1536
    meta = json.loads(paths[2].read_text())
    res = M.expected(codes, cfg)
    assert meta["expected_P_int"] == [int(v) for v in res.P_int]
    assert meta["predicted_cycles_to_done"] == M.expected_cycles(res, cfg)
    assert d["bt_stop_cycles"].item() == meta["recommended_stop_cycles"] == M.recommended_stop(res, cfg)
    rows = list(csv.DictReader(open(paths[1])))
    assert len(rows) == 1536 and int(rows[5]["code"]) == int(codes[5]) and float(rows[5]["value"]) == st[5, 1]


def test_export_siggen_source_is_the_generator_record(tmp_path):
    cfg = M.McConfig(wl=12, nchunks=4)
    got = EXP.make_record("siggen", cfg, tone_bin=siggen.tone_bin(siggen.selftest_record()["pinc"]), noise_shift=3)
    assert np.array_equal(got, siggen.selftest_record(12, 4)["x"])
    tone_only = EXP.make_record("siggen", cfg, noise_shift=-1)
    assert np.array_equal(tone_only, siggen.generate(2048, siggen.pinc_for(100.37), None, 12))


def test_export_cli(tmp_path, capsys):
    EXP.main(["--source", "tone", "--nchunks", "3", "--out-dir", str(tmp_path), "--name", "t", "--tone-bin", "50.5", "--seed", "4"])
    out = capsys.readouterr().out
    assert "recommended simulation length" in out and (tmp_path / "t.mat").exists()
    with pytest.raises(SystemExit):
        EXP.main(["--source", "nosuch"])


# ---------------------------------------------------------------------------------------------------------------------
# result comparison tool
# ---------------------------------------------------------------------------------------------------------------------
def _write_mat(path, log, cfg, as_struct=True):
    d = dict(log, wl=float(cfg.wl), nchunks=float(cfg.nchunks))
    savemat(path, {"bt_log": d} if as_struct else d)


def test_compare_tool_passes_a_correct_log(tmp_path, cfg4, codes4, capsys):
    log, _ = synth_log(cfg4, codes4, idle_before=2)
    p = tmp_path / "bt_log.mat"
    _write_mat(p, log, cfg4)
    assert CMP.main([str(p)]) == 0
    out = capsys.readouterr().out
    assert "PASS: 1024 of 1024 bins" in out and "WL 12, 4 chunks, 3 segments" in out


def test_compare_tool_reads_separate_variables_and_simulink_structures(tmp_path, cfg4, codes4):
    log, _ = synth_log(cfg4, codes4)
    p = tmp_path / "flat.mat"
    _write_mat(p, log, cfg4, as_struct=False)
    assert CMP.main([str(p)]) == 0
    q = tmp_path / "swt.mat"  # 'Structure With Time' of To Workspace blocks
    d = {k: {"time": np.arange(len(v)), "signals": {"values": v.reshape(-1, 1), "dimensions": 1, "label": k}} for k, v in log.items()}
    savemat(q, {"bt_log": dict(d, wl=12.0, nchunks=4.0)})
    assert CMP.main([str(q)]) == 0


def test_compare_tool_csv_and_options(tmp_path, cfg4, codes4, capsys):
    log, _ = synth_log(cfg4, codes4, lag=2)
    p = tmp_path / "log.csv"
    with open(p, "w") as f:
        f.write("in_code,in_valid,rd_addr,rd_data,done\n")
        for row in zip(*(log[k] for k in CMP.SIGNALS)):
            f.write(",".join(str(int(v)) for v in row) + "\n")
    assert CMP.main(["--csv", str(p), "--wl", "12", "--nchunks", "4"]) == 1
    assert "hint: the log matches the model with --lag 2" in capsys.readouterr().out
    assert CMP.main(["--csv", str(p), "--wl", "12", "--nchunks", "4", "--lag", "2"]) == 0


def test_compare_tool_failure_and_error_exit_codes(tmp_path, cfg4, codes4, capsys):
    log, _ = synth_log(cfg4, codes4)
    d0 = int(np.flatnonzero(log["done"])[0])
    bad = dict(log)
    bad["rd_data"] = log["rd_data"].copy()
    bad["rd_data"][d0 + 5:] += 1  # every bin off by one
    p = tmp_path / "bad.mat"
    _write_mat(p, bad, cfg4)
    assert CMP.main([str(p)]) == 1
    assert "FAIL" in capsys.readouterr().out
    assert CMP.main([str(tmp_path / "missing.mat")]) == 2
    q = tmp_path / "partial.mat"
    savemat(q, {"bt_log": {"in_code": log["in_code"]}})
    assert CMP.main([str(q)]) == 2
    assert "missing signal" in capsys.readouterr().err
    short = {k: v[:500] for k, v in log.items()}
    r = tmp_path / "short.mat"
    _write_mat(r, short, cfg4)
    assert CMP.main([str(r)]) == 1
    assert "cannot be evaluated" in capsys.readouterr().out
    with pytest.raises(SystemExit):
        CMP.main([])
