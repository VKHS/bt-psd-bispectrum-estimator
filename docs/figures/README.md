# Figure sheets

Thirty-one sheets, each as PNG (for reading) and SVG (for editing). They are drawn by the scripts in `python/scripts/`
(`make figures` runs all of them) and share one style (`python/btpsd/viz.py`). The paper's own figures are in
`docs/paper_figures/`; the manifest there says which sheet is the counterpart of which paper figure.

**What each sheet is based on.** Three kinds of source are used and never mixed on one number:

| Basis | Meaning |
|---|---|
| *model* | computed by `python/btpsd` on synthetic signals (bit-accurate fixed-point models, floating-point references) |
| *Verilog* | read from `rtl/` by `btpsd.rtlparse` (hierarchy, states, transitions, register widths); the generator fails if a drawn transition is not in the source |
| *simulation* | cycle traces printed by the testbenches during the Icarus Verilog run of the CI workflow, stored in `results/rtl_*.csv` with the run id and commit in `results/rtl_sim_provenance.json` |
| *MATLAB source* | block names, ports, types and parameters read as text from `modelcomposer/matlab/*.m`; those scripts have not been run where this repository was written (no MATLAB there), so the sheet shows what they set, not what Simulink did |

Values reported in the paper appear only as hollow diamonds or in text marked *paper*.

## Architecture and datapath

| Sheet | Content | Script | Basis |
|---|---|---|---|
| `fig01_system_overview` | datapath of the PSD and bispectrum estimator, with a PSD example | `make_diagrams_a.py` | model |
| `fig02_memory_fft_architecture` | eight banks, rotators, butterfly, write-back loop | `make_diagrams_a.py` | model |
| `fig03_r8_butterfly_network` | radix-8 butterfly as three radix-2 DIF stages | `make_diagrams_a.py` | model |
| `fig04_segment_merge` | exact even bins, interpolated odd bins | `make_diagrams_b.py` | model |
| `fig05_fd6_filter` | folded six-tap fractional-delay filter, shift-add constants | `make_diagrams_b.py` | model |
| `fig06_cru_rotator` | constant rotator: octant reduction, ROM, rotation, gain, quadrant | `make_diagrams_b.py` | model |
| `fig07_window_squarer_accumulator` | frequency-domain window, sum-of-odds squarer, accumulator | `make_diagrams_b.py` | model |
| `fig08_bispectrum_datapath` | direct FFT triple product | `make_diagrams_b.py` | model |
| `fig09_bank_address_map` | conflict-free bank and row map of the in-place FFT | `make_diagrams_c.py` | model |
| `fig10_fft_schedule` | three passes of 64 butterflies, latency, the pass 0 to 1 stall | `make_diagrams_c.py` | model |
| `fig11_fd_kernel_gallery` | half-sample fractional-delay kernels with 3 to 8 taps | `make_diagrams_c.py` | model |

## RTL

| Sheet | Content | Script | Basis |
|---|---|---|---|
| `fig15_rtl_hierarchy` | module hierarchy and interfaces, instance counts | `make_rtl_figures.py` | Verilog |
| `fig16_fsm_top_levels` | state machines of `bt_psd_top`, `bt_bispectrum_top` and `fft512_mem` | `make_rtl_figures.py` | Verilog |
| `fig17_fsm_engines` | state machines of `psd_engine` and `bispectrum_engine` | `make_rtl_figures.py` | Verilog |
| `fig18_pipeline_timing` | rotator wavefront, butterfly lane diagram, write-back hazard and slack of all 64 butterflies | `make_rtl_figures.py` | Verilog, model |
| `fig19_interface_timing` | waveforms of `fft512_mem` and `mag_sq_odd`, residence time in `W_S` (model prediction as markers) | `make_sim_figures.py` | simulation |
| `fig20_cycle_budget` | cycles per state of both engines (per-state formulas as markers), phases of both top levels | `make_sim_figures.py` | simulation |
| `fig23_rtl_config_matrix` | which of the 25 configurations of the two top levels were simulated (word length, window, rounding, squarer), each bit-exact to the model, with the clamp count of the vectors | `make_sim_figures.py` | simulation |

## Results of the model

| Sheet | Content | Script | Basis |
|---|---|---|---|
| `fig12_wordlength_study` | SQNR against word length, error budget, window choice at WL = 10 | `make_plots.py` | model |
| `fig13_bispectrum_atlas` | 2-D bispectrum, signal class x window | `make_plots.py` | model |
| `fig14_bispectrum_views` | cumulant, magnitude and phase views | `make_plots.py` | model |
| `fig21_window_responses` | impulse response, DFT bins, response in dB and figures of merit of the four windows | `make_signal_figures.py` | model |
| `fig22_psd_examples` | records and PSD at WL = 10 against the float Welch average | `make_signal_figures.py` | model |

## The FFT in detail (`docs/FFT.md`)

| Sheet | Content | Script | Basis |
|---|---|---|---|
| `fig24_fft_index_space` | index cubes, address digits, bank map, factorisation of the 512-point transform | `make_fft_figures.py` | model, Verilog |
| `fig25_fft_pass_snapshots` | memory contents after the load and after each of the three passes | `make_fft_figures.py` | model |
| `fig26_fft_core_blocks` | signal-level block diagram of `fft512_mem`, the butterfly pipeline, the signals of one butterfly | `make_fft_figures.py` | Verilog, model, simulation |
| `fig27_twiddle_atlas` | the twiddle factors each pass uses, rotation counts, error of the rotator ROM entries | `make_fft_figures.py` | model |
| `fig28_fft_accuracy` | fixed-point FFT accuracy against word length, rounding, input level, per-pass error accrual | `make_fft_figures.py` | model |

## What the estimator resolves

| Sheet | Content | Script | Basis |
|---|---|---|---|
| `fig29_spectral_study` | a tone between two bins, leakage floor against tone position, two-tone resolution, spread of the averaged estimate; hardware model against the exact Welch average | `make_spectral_figures.py` | model |

## Model Composer test bench (`docs/MODEL_COMPOSER.md`)

| Sheet | Content | Script | Basis |
|---|---|---|---|
| `fig30_mc_test_bench` | the Simulink / Model Composer bench (sources, gateways, Black Box, logging), how a run is checked against the bit-exact model, the Black Box interface | `make_mc_figures.py` | MATLAB source, Verilog, model |
| `fig31_mc_wrapper_generator` | record-and-replay wrapper (state machine, conditions), predicted clock-cycle timeline of one run, the in-fabric test-signal generator and self-test, the self-test record and its spectrum | `make_mc_figures.py` | Verilog, model (cycle counts are a prediction) |

## Reading the RTL sheets

* The cycle counts on `fig19` and `fig20` belong to the **sequential reference engines** of this repository (one
  multiply-accumulate at a time) and to the testbench data of `sim/vectors`. They are not comparable with the figures of
  Table 9 of the paper, which describes an overlapped architecture and was not simulated here.
* The RTL is simulated, not synthesised: no sheet shows area, frequency or power.
* The clock-cycle timeline of `fig31` (b) is a **prediction** of the model, drawn hatched and labelled as such; the counts
  that the CI simulation prints are compared with it in `docs/MODEL_COMPOSER.md`. The Simulink side of `fig30` is read from
  the MATLAB scripts and has not been executed.
* `fig23` compares the RTL with the model: the expected values of every run come from the model, so it does not test the
  model against the paper.
* `python/tests/test_rtl_structure.py` and `python/tests/test_sim_results.py` tie the drawn structure and the printed
  numbers to the Verilog and to the model formulas.
