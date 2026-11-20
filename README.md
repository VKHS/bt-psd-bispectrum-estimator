# Memory-based radix-8 FFT power spectral and bispectral estimator

[![ci](https://github.com/VKHS/bt-psd-bispectrum-estimator/actions/workflows/ci.yml/badge.svg)](https://github.com/VKHS/bt-psd-bispectrum-estimator/actions/workflows/ci.yml)

Bit-accurate Python models, a reference RTL implementation, golden test vectors and figures for the estimator of

> A. Khalili Sadaghiani, B. Forouzandeh, *High-performance power spectral/bispectral estimator for biomedical signal
> processing applications using novel memory-based FFT processor*, Integration, the VLSI Journal 99 (2024) 102241.

The design computes a 1024-bin power spectrum from overlapping 512-sample chunks of EEG, EOG or EMG signals with one
memory-based radix-8 FFT, a 6-tap fractional-delay merge of the chunk spectra, a frequency-domain window and a
sum-of-odds squarer, and a direct bispectrum `B(f1, f2) = F(f1) F(f2) F*(f1 + f2)` from the same FFT core.

**Status.** The Python model is tested (119 tests). The RTL is simulated: `make -C sim all` compiles it with Icarus
Verilog (`iverilog -g2005 -Wall`) and runs six self-checking testbenches, which match the golden vectors of the model bit
for bit (rotator 9216 vectors, butterfly 2048, squarer 1024, FFT 4 chunks, PSD 1024 bins, bispectrum 32 x 32 bins).
The PSD and bispectrum inputs end in full-scale, hard-clipped segments, on which the butterfly clamp fires inside the
top-level runs. GitHub Actions repeats the simulation on every push. Twenty-two figure sheets are drawn from the models,
from the Verilog sources and from the simulation traces; the 32 figures of the paper are included with a manifest. Not
done: simulation with another simulator, synthesis, timing, power. The numbers of the paper are quoted as *paper*; everything measured here is *model* or *RTL simulation*
and uses synthetic signals.

| Part | Where | State |
|---|---|---|
| Floating-point references, fixed-point models, schedule model | `python/btpsd/` | tested |
| Reference RTL (rotator, butterfly, FFT core, PSD and bispectrum engines) | `rtl/` | simulated (Icarus Verilog), lint-clean; not synthesised |
| Testbenches and golden vectors | `tb/`, `sim/` | 6 of 6 pass (`make sim`, CI) |
| Experiments (model on synthetic signals; RTL simulation data with provenance) | `results/` (`results/README.md`) | `make experiments` |
| Figure sheets, 22 (architecture, RTL, simulation, windows, PSD examples, bispectra) | `docs/figures/` (`docs/figures/README.md`) | `make figures`; byte-identical on re-runs |
| The 32 figures of the paper, with the counterpart sheet of each | `docs/paper_figures/` | manifest and rights note in its `README.md` |
| Notes on the paper (findings, candidate errata) | `docs/PAPER_NOTES.md` | |
| Paper-to-code map: sections, equations (1)-(32), tables 1-10, figures 1-32 | `docs/PAPER_TO_CODE.md` | complete |
| The ten tables of the paper, as printed | `data/` | all ten present, loaders and tests |
| Architecture, arithmetic contract, differences from the paper | `docs/ARCHITECTURE.md` | |

## Quick start

```bash
pip install -r requirements.txt     # numpy, scipy, matplotlib, pillow, pytest (tested with numpy 2.5, scipy 1.18, Python 3.13)
make test          # unit tests of the model and of the RTL transliteration
make lint          # structural check of the Verilog
make sim           # RTL testbenches (Icarus Verilog); prints TB_RESULT: PASS/FAIL per testbench
make experiments   # results/*.csv from the models (about half a minute)
make figures       # docs/figures/* (22 sheets; fig19 and fig20 read results/rtl_*.csv)
```

```python
import sys; sys.path.insert(0, "python")
from btpsd import estimator as E, bispectrum as B, signals as sg
from btpsd.metrics import sqnr_db

x = sg.pink(20480, seed=0)                                    # synthetic test signal
res = E.psd_fixed(x, E.HwConfig(wl=12, window="hann"))        # bit-accurate hardware model, 12-bit words
print(sqnr_db(E.psd_exact_welch(x, window="hann"), res.to_float()))

x2, _ = sg.qpc(512 * 8, coupled=True, seed=0)                 # quadratic phase coupling at (f1, f2)
re, im, S, shift = B.bispectrum_fixed(x2, B.BispecConfig(wl=12, max_bin=32))
Bxy = B.bispectrum_hw_to_float(re, im, 12, S, shift)          # 32 x 32 bispectrum, analytic scale
```

The Vivado script `fpga/synth_ooc.tcl` (untested) produces out-of-context utilisation and timing reports.

## Results

The two blocks below are never compared digit by digit: the paper's records are not available here, and its SQNR
formula and test signals are not fully specified (`docs/PAPER_NOTES.md`, section 2.3).

**Paper-reported** (Table 9, Artix-7, WL = 10; reproduced here as `data/paper_table9_fpga_comparison.csv`)

| Design | DSP | LUT | FF | Throughput | Frequency | Power | Latency |
|---|---:|---:|---:|---:|---:|---:|---:|
| Proposed estimator | 0 | 4834 | 4420 | 205 Msps (R = 0.88) | 232 MHz | 191 mW | 1252 cycles |

**Measured in this repository** (model on synthetic signals, 5 seeds, `results/`; the RTL rows are cycle counts from the Icarus Verilog testbenches)

| Quantity | Value |
|---|---|
| 512-point FFT, RTL pipeline (rotator 4, butterfly 10, loop latency 12) | 196 issue cycles + 11 = 207 cycles, 4 stall cycles (schedule model); the testbench measures 207 cycles from the start edge to `done`, 209 per chunk in the FFT phase of the PSD and bispectrum top levels |
| Bispectrum engine, 32 x 32 bins (RTL simulation) | 5219 cycles per frame (FSM count 5216), plus 209 for the FFT |
| PSD segment engine (RTL simulation, serial squarer, vectors with a full-scale tail) | 115 462 cycles for 3 segments, 38.5 k per segment, of which 7.3 k are squarer cycles on average; a sequential reference engine, not comparable to the overlapped design of the paper |
| Constant-rotator ROM, tolerance 0.01 rad | 129 entries, maximum angle error 0.0097 rad, 7.6 adders per rotation on average, 512-point FFT SNR 46.6 dB at WL = 16 |
| Odd-bin accuracy of the 6-tap merge | 11.7 dB (analytic bound 11.74 dB); Table 3 of the paper reports 50.08 dB |
| PSD SQNR against the exact Welch average, WL = 10 / 12 / 16 | white 15.4 / 23.9 / 23.5 dB; pink 28.4 / 31.9 / 29.0 dB; EEG-like 28.7 / 40.0 / 44.2 dB |
| Bispectrum magnitude SQNR (coupled signal), WL = 10 / 12 / 14 / 16 / 18 | 22.7 / 34.6 / 43.1 / 52.2 / 59.6 dB |
| Serial sum-of-odds squarer, squarer pairs needed for R = 0.88 | 1.0 to 1.6 at WL = 10, 44 to 96 at WL = 16 |
| Periodic windows, equivalent noise bandwidth (bins) and peak sidelobe, rectangular / Hann / Hamming / Blackman (`fig21`) | 1.000 / 1.500 / 1.363 / 1.727 and -13.3 / -31.5 / -42.7 / -58.1 dB; the frequency-domain kernels of 1 / 3 / 3 / 5 bins are exact |
| PSD at WL = 10 on the 5120-sample records of `fig22`: SQNR against the float Welch average, EEG-like / AR(2) narrow-band / white | 27.8 / 18.1 / 12.9 dB; one LSB of the PSD word is -21.6 dB, so low-power bins quantise to zero (119 and 200 of 513 bins for the first two signals) |

Reading the table: from WL = 14 the PSD SQNR is set by the merge, not by the word length, because six taps keep
about 93 % of the interpolation kernel energy. `docs/PAPER_NOTES.md` explains this, lists what the paper leaves open,
and lists fifteen candidate errata with the printed and the expected form (one of them: the window panels of the paper's
Figs. 13-15 are drawn with symmetric windows, whose frequency-domain kernel is not the exact 3- or 5-bin kernel of Eq. 13).

## Figures

The 22 sheets in `docs/figures/` (PNG for reading, SVG for editing) are drawn by `python/scripts/make_*.py`; each one is
based on the model, on the Verilog sources or on the simulation traces, never on a mixture on one number
(`docs/figures/README.md` lists script and basis of every sheet). The paper's own figures are in `docs/paper_figures/`
(manifest and rights note in that folder's `README.md`).

**Architecture and datapath**

| | |
|---|---|
| ![System overview](docs/figures/fig01_system_overview.png) | ![Memory FFT architecture](docs/figures/fig02_memory_fft_architecture.png) |
| `fig01` system overview | `fig02` memory-based FFT architecture |
| ![Bank address map](docs/figures/fig09_bank_address_map.png) | ![FFT schedule](docs/figures/fig10_fft_schedule.png) |
| `fig09` conflict-free bank and row map | `fig10` schedule, latency and the pass 0 to 1 stall |

**Reference RTL** (structure read from the Verilog; timing and cycle counts from the Icarus Verilog simulation)

| | |
|---|---|
| ![RTL hierarchy](docs/figures/fig15_rtl_hierarchy.png) | ![Top-level state machines](docs/figures/fig16_fsm_top_levels.png) |
| `fig15` module hierarchy and interfaces | `fig16` state machines of the top levels and the FFT core |
| ![Engine state machines](docs/figures/fig17_fsm_engines.png) | ![Pipeline timing](docs/figures/fig18_pipeline_timing.png) |
| `fig17` state machines of the PSD and bispectrum engines | `fig18` rotator and butterfly pipelines, write-back hazard and slack |
| ![Interface timing](docs/figures/fig19_interface_timing.png) | ![Cycle budget](docs/figures/fig20_cycle_budget.png) |
| `fig19` simulated waveforms of the FFT core and the squarer | `fig20` simulated cycles per state and per phase |

**Model results** (synthetic signals)

| | |
|---|---|
| ![Word-length study](docs/figures/fig12_wordlength_study.png) | ![Bispectrum atlas](docs/figures/fig13_bispectrum_atlas.png) |
| `fig12` word-length and error budget | `fig13` 2-D bispectrum atlas, three windows |
| ![Window responses](docs/figures/fig21_window_responses.png) | ![PSD examples](docs/figures/fig22_psd_examples.png) |
| `fig21` the four windows: kernels, responses, figures of merit | `fig22` records and WL = 10 PSD against the float Welch average |

Also in `docs/figures/`: `fig03` radix-8 butterfly network, `fig04` segment merge, `fig05` six-tap filter, `fig06`
rotator, `fig07` window, squarer and accumulator, `fig08` bispectrum datapath, `fig11` fractional-delay kernels,
`fig14` bispectrum views.

**Bispectra from the paper** (Figs. 24-32, recorded EEG-1, EOG-1 and EMG-1; windows rectangular, Hann, Hamming; click for the full size)

| | Rectangular | Hann | Hamming |
|---|---|---|---|
| EEG-1 | [<img src="docs/paper_figures/paperfig24_bispectrum_eeg1_rectangular.jpg" width="300">](docs/paper_figures/paperfig24_bispectrum_eeg1_rectangular.jpg) | [<img src="docs/paper_figures/paperfig25_bispectrum_eeg1_hann.jpg" width="300">](docs/paper_figures/paperfig25_bispectrum_eeg1_hann.jpg) | [<img src="docs/paper_figures/paperfig26_bispectrum_eeg1_hamming.jpg" width="300">](docs/paper_figures/paperfig26_bispectrum_eeg1_hamming.jpg) |
| EOG-1 | [<img src="docs/paper_figures/paperfig27_bispectrum_eog1_rectangular.jpg" width="300">](docs/paper_figures/paperfig27_bispectrum_eog1_rectangular.jpg) | [<img src="docs/paper_figures/paperfig28_bispectrum_eog1_hann.jpg" width="300">](docs/paper_figures/paperfig28_bispectrum_eog1_hann.jpg) | [<img src="docs/paper_figures/paperfig29_bispectrum_eog1_hamming.jpg" width="300">](docs/paper_figures/paperfig29_bispectrum_eog1_hamming.jpg) |
| EMG-1 | [<img src="docs/paper_figures/paperfig30_bispectrum_emg1_rectangular.jpg" width="300">](docs/paper_figures/paperfig30_bispectrum_emg1_rectangular.jpg) | [<img src="docs/paper_figures/paperfig31_bispectrum_emg1_hann.jpg" width="300">](docs/paper_figures/paperfig31_bispectrum_emg1_hann.jpg) | [<img src="docs/paper_figures/paperfig32_bispectrum_emg1_hamming.jpg" width="300">](docs/paper_figures/paperfig32_bispectrum_emg1_hamming.jpg) |

## Repository layout

```text
python/btpsd/        models: fixedpoint, windows, fdfilter, reference (Welch), cordic (rotator ROM), r8fft + schedule,
                     squarer, estimator (PSD), bispectrum, rtlmodel (transliteration of the RTL), vectors, signals,
                     metrics, paperdata, rtlparse (reads the Verilog), viz
python/tests/        119 tests: model, RTL transliteration, structure of the Verilog, simulation results, figure data
python/scripts/      run_experiments; make_diagrams_a/b/c, make_plots, make_rtl_figures, make_sim_figures,
                     make_signal_figures (the figure sheets); ingest_ci_traces (simulation output -> results/);
                     measure_paper_window_stems; gen_cru_rom, gen_rtl, gen_vectors, extract_paper_tables
rtl/                 Verilog-2001 reference implementation (see rtl/README.md)
tb/                  six self-checking testbenches
sim/                 Makefile for Icarus Verilog, golden vectors (sim/vectors, see its README.md)
scripts/lint_rtl.py  structural check of the Verilog sources
fpga/                out-of-context Vivado synthesis script (untested)
data/                the ten tables of the paper as printed, values read from its figures, generated ROM tables
results/             model experiments (provenance.json) and RTL simulation data (rtl_*.csv, rtl_sim_provenance.json)
docs/                ARCHITECTURE.md, PAPER_NOTES.md, PAPER_TO_CODE.md, figures/ (22 sheets), paper_figures/ (32)
CITATION.cff         how to cite the article
legacy/              the original files of the repository, unchanged (legacy/README.md)
```

The original files of the repository are kept unchanged in `legacy/`; the maintained sources are the ones listed in
`rtl/README.md` and `python/btpsd/`.

## Further reading

1. `docs/ARCHITECTURE.md` - dataflow, arithmetic contract shared by model and RTL, FFT schedule, differences from the
   paper's overlapped design, verification map.
2. `docs/PAPER_NOTES.md` - findings on the paper's claims, SQNR definitions, squarer throughput, candidate errata.
3. `docs/PAPER_TO_CODE.md` - every section, equation, table and figure of the paper and where it lives here.
4. `docs/figures/README.md`, `docs/paper_figures/README.md` - the 22 sheets with their basis, the paper's 32 figures.
5. `results/README.md` - which file is model output, which is RTL simulation output, and how each is produced.
6. `rtl/README.md`, `sim/vectors/README.md`, `data/README.md` - RTL files, parameters, vector formats, tables.

## Citing

Please cite the article above; `CITATION.cff` holds the reference (GitHub shows it under "Cite this repository").
