# Results

Everything in this folder was produced by this repository. Nothing is copied from the paper (paper-reported values are in
`data/paper_*.csv`) and nothing comes from recorded signals: the model files use synthetic signals
(`python/btpsd/signals.py`), the simulation files use the test vectors of `sim/vectors`.

## Model (synthetic signals)

Regenerate with `make experiments` (`python/scripts/run_experiments.py`, about half a minute); `provenance.json` records
the record length, the seeds and the library versions. SQNR = `10 log10(sum P_ref^2 / sum (P_ref - P_hw)^2)` in the linear
power domain with the analytic output scale (no fitting).

| File | Content | Produced by |
|---|---|---|
| `fd_merge_snr_vs_taps.csv` | odd-bin accuracy of the fractional-delay merge against the number of taps | `exp_table3_analog` |
| `sqnr_vs_wordlength.csv` | PSD SQNR against word length, ROM tolerance and signal | `exp_sqnr_vs_wl` |
| `cru_rom_tradeoff.csv` | rotator ROM size, angle error and adders against the angle tolerance | `exp_rom_tradeoff` |
| `bispectrum_sqnr_vs_wordlength.csv` | bispectrum SQNR against word length | `exp_bispectrum` |
| `fft512_cycles_vs_latency.csv` | cycles of one 512-point FFT against the pipeline latency (schedule model) | `exp_cycle_model` |
| `squarer_cycle_model.csv` | serial squarer cycles per segment and squarers needed for R = 0.88 | `exp_cycle_model` |
| `model_sqnr_wl10_synthetic.csv` | SQNR at WL = 10 for four signals and four windows | `exp_wl10_windows` |
| `model_psd_examples_wl10.csv` | the numbers printed on `docs/figures/fig22_psd_examples` (one record per signal class) | `python/scripts/make_signal_figures.py` |

## RTL simulation

Cycle traces and counts printed by the testbenches in `tb/` during the Icarus Verilog run of the CI workflow. They are
stored by `python/scripts/ingest_ci_traces.py`; `rtl_sim_provenance.json` names the workflow run and the commit. No Python
model is involved in producing them; `python/tests/test_sim_results.py` compares them with the model formulas.

| File | Content |
|---|---|
| `rtl_sim_provenance.json` | run id and URL, commit, testbenches that passed, description of each file |
| `rtl_trace_fft512.csv` | `fft512_mem`, first transform: `ld_en`, `start`, `busy`, `done`, issue valid, pass, butterfly, state, per cycle (cycle 0 = `start` high) |
| `rtl_trace_squarer.csv` | `mag_sq_odd`, three directed examples: `start`, `busy`, `done`, inputs, product, per cycle |
| `rtl_state_occupancy.csv` | cycles per FSM state of `bt_psd_top` / `psd_engine` and `bt_bispectrum_top` / `bispectrum_engine`, from reset release to the last result |
| `rtl_squarer_wait_hist.csv` | histogram of the time `psd_engine` waits in state `W_S` per bin |
| `rtl_config_matrix.csv` | one row per configuration of `python/scripts/run_rtl_matrix.py` (25 configurations (16 PSD, 9 bispectrum): word length 10, 12, 14 and 16, the four windows, `ROUND` 0 and 1, serial squarer and multiplier): result of the top-level testbench against the model, words compared, clamped differences in the vectors |
| `rtl_config_matrix_provenance.json` | the CI run and commit that produced `rtl_config_matrix.csv`, with the meaning of its columns |

The simulated engines are the sequential reference engines of this repository. Their cycle counts say nothing about the
overlapped architecture of the paper (Table 9) and are not to be compared with it; the RTL is not synthesised, so there is
no area, frequency or power result.
