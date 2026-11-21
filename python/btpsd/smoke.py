"""Smoke test of the repository: one command, a few seconds to half a minute, no simulator required.

    python3 python/scripts/smoke_test.py              # all checks; Icarus Verilog is used if it is installed
    python3 python/scripts/smoke_test.py --quick      # model and file checks only (no vector regeneration, no simulation)
    python3 python/scripts/smoke_test.py --xsim       # also run sim/xsim_smoke.sh with the Vivado simulator
    python3 python/scripts/smoke_test.py --list       # names of the checks
    python3 python/scripts/smoke_test.py --only fft_model,psd_model

(equivalently ``make smoke``).  The exit status is 0 if no check failed.  A check that cannot run because a tool is missing is
reported as SKIP, not as PASS.

What it establishes, in order:

* the Python environment and the files of the checkout are what the other tools expect;
* the generated RTL constants (``rtl/cru_rom.v``, ``rtl/btpsd_consts.vh``) and the golden vectors in ``sim/vectors`` are exactly
  what the models produce today (nothing stale);
* the FFT schedule is conflict-free, takes the cycle count the Verilog parameters imply, and the cycle-level model of the FFT core
  reproduces the bit-accurate FFT;
* the FFT, the PSD estimator and the bispectrum estimator agree with floating-point references to the accuracy documented in
  ``docs/``;
* the Model Composer pack (``modelcomposer/``): its generated files are current, the model of the test-signal generator and of
  the record-and-replay wrapper hang together, and (with a simulator) the wrapper, the bundle, the generator and the self-test
  pass their testbenches;
* the Verilog passes the structural lint, and (with a simulator) every self-checking testbench passes against the golden vectors;
* the figure framework renders deterministic files.

It is a smoke test: a pass says the repository is healthy, not that the Verilog is correct for every configuration (that is the job
of ``make -C sim matrix`` and the unit tests) and not that anything was synthesised.  The xsim step has not been executed in the
environment this repository was developed in (no Vivado there); it is a script for a machine that has one.
"""
from __future__ import annotations

import argparse
import importlib
import platform
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np

ROOT = Path(__file__).resolve().parents[2]

# thresholds: measured values (WL = 16, seeded records) are 48 dB, 66.8 dB, 45.4 dB and 48.4 dB; the margins catch breakage, not drift
FFT_SQNR_MIN_DB = 40.0
PSD_VS_FLOAT_MIN_DB = 50.0  # quantisation only (the floating-point six-tap merge is the reference)
PSD_VS_EXACT_MIN_DB = 35.0  # includes the error of the six-tap merge, see docs/PAPER_NOTES.md section 2.2
BSP_MIN_DB = 35.0
TESTBENCHES = ("cru", "bfly", "sq", "fft", "psd", "bsp", "mc", "mcb", "sg", "st")  # targets of sim/Makefile


class Skip(Exception):
    """Raised by a check that cannot run in this environment."""


@dataclass
class Result:
    name: str
    status: str  # PASS | FAIL | SKIP
    detail: str
    seconds: float


CHECKS: list[tuple[str, Callable[[], str], dict]] = []


def check(name: str, slow: bool = False, optional: bool = False):
    """Register a check.  ``slow`` checks are skipped by ``--quick``; ``optional`` ones run only on request (``--xsim``)."""
    def wrap(fn):
        CHECKS.append((name, fn, {"slow": slow, "optional": optional, "doc": (fn.__doc__ or "").strip().split("\n")[0]}))
        return fn
    return wrap


def _script(name: str):
    """Import a module of ``python/scripts`` (they are scripts, not part of the package)."""
    d = str(ROOT / "python" / "scripts")
    if d not in sys.path:
        sys.path.insert(0, d)
    return importlib.import_module(name)


def _run(cmd: list[str], cwd: Path = ROOT, timeout: int = 900) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout)


# ---------------------------------------------------------------------------------------------------------------
# environment and files
# ---------------------------------------------------------------------------------------------------------------
@check("environment")
def env_check() -> str:
    """Python, numpy and scipy are new enough; matplotlib is optional."""
    import scipy

    def ver(s: str) -> tuple[int, ...]:
        return tuple(int(p) for p in s.split(".")[:2])

    assert sys.version_info >= (3, 10), f"Python 3.10 or newer needed, found {platform.python_version()}"
    assert ver(np.__version__) >= (1, 24), f"numpy 1.24 or newer needed, found {np.__version__}"
    assert ver(scipy.__version__) >= (1, 10), f"scipy 1.10 or newer needed, found {scipy.__version__}"
    try:
        import matplotlib
        mpl = matplotlib.__version__
    except ImportError:
        mpl = "not installed (figure checks will be skipped)"
    return f"python {platform.python_version()}, numpy {np.__version__}, scipy {scipy.__version__}, matplotlib {mpl}"


@check("layout")
def layout_check() -> str:
    """The files of the checkout that the simulations and the figures read are present."""
    files = [ln.strip() for ln in (ROOT / "sim" / "filelist.f").read_text().splitlines() if ln.strip().endswith(".v")]
    need = files + ["sim/Makefile", "rtl/btpsd_params.vh", "rtl/btpsd_consts.vh", "README.md", "docs/ARCHITECTURE.md", "requirements.txt"]
    need += [f"tb/tb_{n}.v" for n in ("cru_rotate", "r8_butterfly", "mag_sq_odd", "fft512_mem", "psd_top", "bispectrum_top",
                                      "bt_psd_mc", "sig_gen", "bt_psd_selftest")]
    need += [f"modelcomposer/hdl/{n}" for n in ("bt_psd_mc.v", "bt_psd_mc_bundle.v", "sig_gen.v", "sig_rom.v", "bt_psd_selftest.v")]
    need += [f"modelcomposer/matlab/{n}" for n in ("bt_mc_build.m", "bt_mc_run.m", "bt_mc_params.m", "bt_psd_mc_config.m")]
    missing = [f for f in need if not (ROOT / f).is_file()]
    assert not missing, "missing: " + ", ".join(missing)
    mems = sorted((ROOT / "sim" / "vectors").glob("*.mem"))
    assert len(mems) >= 14, f"only {len(mems)} .mem files in sim/vectors"
    assert all(m.stat().st_size > 0 for m in mems), "empty .mem file"
    return f"{len(files)} RTL files, {len(mems)} vector files"


@check("generated_rtl")
def generated_rtl_check() -> str:
    """rtl/cru_rom.v and rtl/btpsd_consts.vh equal what python/scripts/gen_rtl.py writes today."""
    g = _script("gen_rtl")
    stale = []
    if (ROOT / "rtl" / "cru_rom.v").read_text() != g.cru_rom_v(0.01):
        stale.append("rtl/cru_rom.v")
    if (ROOT / "rtl" / "btpsd_consts.vh").read_text() != g.consts_vh():
        stale.append("rtl/btpsd_consts.vh")
    assert not stale, "stale (run `make -C sim vectors`): " + ", ".join(stale)
    return "cru_rom.v and btpsd_consts.vh are current"


@check("golden_vectors", slow=True)
def vectors_check() -> str:
    """Regenerating the golden vectors reproduces sim/vectors byte for byte."""
    g = _script("gen_vectors")
    with tempfile.TemporaryDirectory() as tmp:
        g.write_all(Path(tmp))
        new = {p.name: p.read_bytes() for p in Path(tmp).iterdir()}
    old = {p.name: p.read_bytes() for p in (ROOT / "sim" / "vectors").iterdir() if p.is_file()}
    assert set(new) == set(old), f"file sets differ: {sorted(set(new) ^ set(old))}"
    bad = sorted(n for n in new if new[n] != old[n])
    assert not bad, "stale (run `make -C sim vectors`): " + ", ".join(bad)
    return f"{len(new)} files identical"


# ---------------------------------------------------------------------------------------------------------------
# models
# ---------------------------------------------------------------------------------------------------------------
@check("schedule")
def schedule_check() -> str:
    """The FFT schedule is conflict-free and its cycle count follows from the Verilog latency parameters."""
    from . import schedule as S
    from .fftstudy import fft_cycle_budget
    from .rtlmodel import LAT, LAT_BFLY, LAT_CRU, stall_cycles
    from .rtlparse import define_values

    assert S.check_conflict_free(), "an FFT pass touches a bank twice"
    v = define_values(ROOT / "rtl" / "btpsd_params.vh")
    assert (v["LAT_CRU"], v["LAT_BFLY"]) == (LAT_CRU, LAT_BFLY), f"rtl/btpsd_params.vh {v} against the model ({LAT_CRU}, {LAT_BFLY})"
    cyc = S.schedule_cycles(LAT)
    assert cyc["stalls_per_pass"] == [0, stall_cycles(LAT), 0], cyc
    b = fft_cycle_budget()
    assert cyc["total_cycles"] == b["total_cycles"] == 207, (cyc, b)
    return f"conflict-free; LAT {LAT}, stall {stall_cycles(LAT)}, {cyc['total_cycles']} cycles first issue to done"


@check("fft_model")
def fft_check() -> str:
    """The cycle-level FFT core model is bit-exact with the vectorised FFT, which is accurate against the exact transform."""
    from .cordic import default_rom
    from .fftstudy import fft_inputs, fft_sqnr
    from .r8fft import r8fft_fixed
    from .rtlmodel import Fft512Rtl

    re, im = fft_inputs("EEG-like", 16, 7)
    dut = Fft512Rtl(16)
    dut.load(re, im)
    done = dut.run()
    yr, yi = dut.spectrum()
    er, ei = r8fft_fixed(re, im, default_rom(0.01), 16, "round")
    assert np.array_equal(yr, er) and np.array_equal(yi, ei), "Fft512Rtl differs from r8fft_fixed"
    assert done == 207, f"done after {done} cycles"
    q = fft_sqnr(re, im, 16)
    assert q >= FFT_SQNR_MIN_DB, f"FFT SQNR {q:.1f} dB < {FFT_SQNR_MIN_DB} dB"
    return f"bit-exact, done at cycle {done}, SQNR {q:.1f} dB (limit {FFT_SQNR_MIN_DB:.0f})"


@check("psd_model")
def psd_check() -> str:
    """The bit-accurate PSD estimator tracks the floating-point six-tap merge and the exact Welch average."""
    from . import signals
    from .estimator import HwConfig, psd_exact_welch, psd_fixed, psd_float_fd
    from .metrics import sqnr_db

    x = signals.eeg_like(512 * 9, fs=256.0, seed=11)
    hw = psd_fixed(x, HwConfig(wl=16)).to_float()
    q_fd = sqnr_db(psd_float_fd(x), hw)
    q_ex = sqnr_db(psd_exact_welch(x), hw)
    assert q_fd >= PSD_VS_FLOAT_MIN_DB, f"PSD against the float six-tap merge: {q_fd:.1f} dB < {PSD_VS_FLOAT_MIN_DB}"
    assert q_ex >= PSD_VS_EXACT_MIN_DB, f"PSD against the exact Welch average: {q_ex:.1f} dB < {PSD_VS_EXACT_MIN_DB}"
    assert int(np.argmax(hw[:512])) == int(np.argmax(psd_exact_welch(x)[:512])), "peak bin differs"
    return f"WL 16: {q_fd:.1f} dB against the float six-tap merge, {q_ex:.1f} dB against the exact average"


@check("bispectrum_model")
def bispectrum_check() -> str:
    """The bit-accurate bispectrum agrees with the direct floating-point estimator and finds the coupled pair."""
    from . import bispectrum as B
    from . import signals
    from .metrics import sqnr_db

    x = signals.qpc(512 * 12, coupled=True, seed=0)[0]
    ref = B.bispectrum_direct(x, L=512, hop=256, nfft=512, max_bin=128)
    re, im, S, sh = B.bispectrum_fixed(x, B.BispecConfig(wl=16))
    hw = B.bispectrum_hw_to_float(re, im, 16, S, sh)
    q = sqnr_db(np.abs(ref), np.abs(hw))
    pk_ref = np.unravel_index(int(np.argmax(np.abs(ref))), ref.shape)
    pk_hw = np.unravel_index(int(np.argmax(np.abs(hw))), hw.shape)
    assert q >= BSP_MIN_DB, f"bispectrum magnitude SQNR {q:.1f} dB < {BSP_MIN_DB}"
    assert pk_ref == pk_hw, f"peak at {pk_hw}, reference at {pk_ref}"
    return f"WL 16, {S} frames: magnitude SQNR {q:.1f} dB, peak at (f1, f2) = {tuple(int(v) for v in pk_hw)}"


# ---------------------------------------------------------------------------------------------------------------
# Model Composer pack
# ---------------------------------------------------------------------------------------------------------------
@check("model_composer_files")
def mc_files_check() -> str:
    """The include-free bundle, the Black Box configuration, the parameter file and the sine ROM equal what gen_mc.py writes today."""
    g = _script("gen_mc")
    names = ("hdl/bt_psd_mc_bundle.v", "matlab/bt_psd_mc_config.m", "matlab/bt_mc_params.m", "hdl/sig_rom.v")
    with tempfile.TemporaryDirectory() as tmp:
        g.write_all(Path(tmp), **g.DEFAULTS)
        new = {n: (Path(tmp) / n).read_text() for n in names}
    stale = [n for n in names if (ROOT / "modelcomposer" / n).read_text() != new[n]]
    assert not stale, "stale (run `python3 python/scripts/gen_mc.py`): " + ", ".join("modelcomposer/" + n for n in stale)
    return f"{len(names)} generated files are current ({g.DEFAULTS['wl']}-bit, {g.DEFAULTS['nchunks']} chunks)"


@check("model_composer_model")
def mc_model_check() -> str:
    """The test-signal generator, the wrapper model and the log decoder hang together on the self-test record."""
    from . import mcmodel as M
    from . import siggen

    rec = siggen.selftest_record()
    cfg = M.McConfig(wl=rec["wl"], nchunks=rec["nchunks"])
    x = rec["x"]
    assert len(x) == cfg.nsamp and np.abs(x).max() < (1 << (cfg.wl - 1)), "self-test record has the wrong length or range"
    assert len(set(int(v) for v in siggen.lfsr_states(4096))) == 4096, "LFSR repeats within 4096 states"
    f = siggen.tone_bin(rec["pinc"])
    assert abs(f - siggen.SELFTEST["bin_1024"]) < 1e-3, f"tone at bin {f}"
    res = M.expected(x, cfg)
    pk = int(np.argmax(res.P_int[:512]))
    assert pk == int(round(f)), f"spectral peak at bin {pk}, tone at {f:.2f}"
    # an ideal Black Box (the model's own result, read through a free-running address counter) passes the checker
    pred = M.expected_cycles(res, cfg)
    n = M.recommended_stop(res, cfg)
    code = np.zeros(n, dtype=np.int64)
    code[: cfg.nsamp] = x
    valid = np.zeros(n, dtype=np.int64)
    valid[: cfg.nsamp] = 1
    done = np.zeros(n, dtype=np.int64)
    done[pred:] = 1
    addr = np.arange(n) % M.BINS
    data = np.zeros(n, dtype=np.int64)
    data[1:] = res.P_int[addr[:-1]]
    rc = M.check_logs(code, valid, addr, data, done, cfg)
    assert rc.ok and rc.cycles_to_done == pred, rc.summary()
    wrong = res.P_int.copy()
    wrong[300] += 1  # one wrong result word must be found
    data[1:] = wrong[addr[:-1]]
    bad = M.check_logs(code, valid, addr, data, done, cfg)
    assert not bad.ok and bad.n_mismatch == 1, bad.summary()
    return f"tone at bin {f:.2f} found at bin {pk}; ideal log passes, one wrong word is caught; predicted {pred} cycles to done"


# ---------------------------------------------------------------------------------------------------------------
# Verilog
# ---------------------------------------------------------------------------------------------------------------
@check("verilog_lint")
def lint_check() -> str:
    """scripts/lint_rtl.py reports no problems (balanced constructs, defined modules, connected ports)."""
    r = _run([sys.executable, "scripts/lint_rtl.py"])
    last = (r.stdout.strip().splitlines() or ["(no output)"])[-1]
    assert r.returncode == 0, "\n".join(r.stdout.strip().splitlines()[-6:] + r.stderr.strip().splitlines()[-3:])
    return last


@check("icarus_testbenches", slow=True)
def icarus_check() -> str:
    """Every self-checking testbench passes against the golden vectors (Icarus Verilog, `make -C sim`)."""
    for tool in ("iverilog", "vvp", "make"):
        if shutil.which(tool) is None:
            raise Skip(f"{tool} not found")
    rows = []
    for t in TESTBENCHES:
        t0 = time.time()
        r = _run(["make", "-C", "sim", t], timeout=1200)
        log = ROOT / "sim" / "work" / f"{t}.log"
        text = log.read_text() if log.exists() else r.stdout
        line = next((ln for ln in text.splitlines() if "TB_RESULT" in ln), "no TB_RESULT line")
        assert r.returncode == 0 and "TB_RESULT: PASS" in line, f"{t}: {line}\n" + "\n".join(r.stderr.strip().splitlines()[-5:])
        rows.append(f"{t} {time.time() - t0:.1f}s")
    return "all pass (" + ", ".join(rows) + ")"


@check("xsim_testbenches", slow=True, optional=True)
def xsim_check() -> str:
    """Every self-checking testbench passes in the Vivado simulator (sim/xsim_smoke.sh).  Not executed in the development environment."""
    if shutil.which("xvlog") is None or shutil.which("xelab") is None or shutil.which("xsim") is None:
        raise AssertionError("xvlog/xelab/xsim not on PATH: source <Vivado install>/settings64.sh first")
    r = _run(["bash", "sim/xsim_smoke.sh"], timeout=3600)
    res = [ln for ln in r.stdout.splitlines() if "TB_RESULT" in ln]
    assert r.returncode == 0 and res and all("PASS" in ln for ln in res), "\n".join((res or r.stdout.splitlines())[-8:] + r.stderr.splitlines()[-4:])
    return f"{len(res)} testbenches pass in xsim"


# ---------------------------------------------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------------------------------------------
@check("figure_framework")
def figure_check() -> str:
    """The figure framework renders SVG and PNG, and two renderings of one diagram are byte-identical."""
    try:
        from .viz import Diagram
    except ImportError as e:  # matplotlib missing
        raise Skip(str(e)) from e

    def render(tmp: Path, tag: str) -> bytes:
        d = Diagram(40, 20, "smoke", subtitle="smoke test")
        d.box(4, 4, 14, 8, "compute\nblock", "compute")
        d.box(24, 4, 12, 8, "memory", "memory")
        d.arrow([(18, 8), (24, 8)], "x")
        d.save(tmp / tag)
        assert (tmp / f"{tag}.png").stat().st_size > 1000 and (tmp / f"{tag}.svg").stat().st_size > 1000
        return (tmp / f"{tag}.svg").read_bytes()

    with tempfile.TemporaryDirectory() as tmp:
        a, b = render(Path(tmp), "a"), render(Path(tmp), "b")
    assert a == b, "two renderings of the same diagram differ"
    return f"SVG and PNG written, {len(a)} bytes of SVG, deterministic"


@check("figures_present")
def figures_present_check() -> str:
    """Every PNG sheet in docs/figures has its SVG, and the figure catalogue lists it."""
    figs = ROOT / "docs" / "figures"
    pngs = sorted(p for p in figs.glob("fig*.png"))
    assert pngs, "no figure sheets found"
    missing = [p.name for p in pngs if not p.with_suffix(".svg").is_file() or p.stat().st_size < 5000]
    assert not missing, "missing SVG or empty PNG: " + ", ".join(missing)
    catalogue = (figs / "README.md").read_text()
    unlisted = [p.stem for p in pngs if p.stem.split("_")[0] not in catalogue]
    assert not unlisted, "not in docs/figures/README.md: " + ", ".join(unlisted)
    return f"{len(pngs)} sheets, SVG present, all catalogued"


# ---------------------------------------------------------------------------------------------------------------
# runner
# ---------------------------------------------------------------------------------------------------------------
def run_checks(names: list[str] | None = None, quick: bool = False, xsim: bool = False) -> list[Result]:
    results = []
    for name, fn, meta in CHECKS:
        if names is not None and name not in names:
            continue
        if meta["optional"] and not xsim and (names is None or name not in names):
            continue
        if quick and meta["slow"] and (names is None or name not in names):
            results.append(Result(name, "SKIP", "skipped by --quick", 0.0))
            continue
        t0 = time.time()
        try:
            detail, status = fn(), "PASS"
        except Skip as e:
            detail, status = str(e), "SKIP"
        except AssertionError as e:
            detail, status = str(e) or "assertion failed", "FAIL"
        except Exception as e:  # noqa: BLE001  a smoke test reports every failure, whatever its type
            detail, status = f"{type(e).__name__}: {e}", "FAIL"
        results.append(Result(name, status, detail, time.time() - t0))
    return results


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true", help="skip the slow checks (vector regeneration, simulation)")
    ap.add_argument("--xsim", action="store_true", help="also run the Vivado simulator testbenches (sim/xsim_smoke.sh)")
    ap.add_argument("--only", help="comma-separated check names")
    ap.add_argument("--list", action="store_true", help="list the checks and exit")
    a = ap.parse_args(argv)
    if a.list:
        for name, _, meta in CHECKS:
            tags = ", ".join(t for t in ("slow", "optional") if meta[t])
            print(f"{name:20s} {meta['doc']}" + (f"  [{tags}]" if tags else ""))
        return 0
    names = a.only.split(",") if a.only else None
    if names:
        unknown = [n for n in names if n not in {c[0] for c in CHECKS}]
        if unknown:
            print("unknown check(s): " + ", ".join(unknown) + "; see --list")
            return 2
    t0 = time.time()
    print(f"btpsd smoke test, repository {ROOT}")
    results = run_checks(names, quick=a.quick, xsim=a.xsim)
    for r in results:
        first, *rest = r.detail.split("\n")
        print(f"[{r.status}] {r.name:20s} {first}  ({r.seconds:.1f} s)")
        for ln in rest:
            print(f"       {ln}")
    n = {s: sum(r.status == s for r in results) for s in ("PASS", "FAIL", "SKIP")}
    print(f"{n['PASS']} passed, {n['FAIL']} failed, {n['SKIP']} skipped in {time.time() - t0:.1f} s")
    return 1 if n["FAIL"] else 0


if __name__ == "__main__":
    sys.exit(main())
