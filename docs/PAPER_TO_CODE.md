# Paper to code

Where every section, equation, table and figure of

> A. Khalili Sadaghiani, B. Forouzandeh, *High-performance power spectral/bispectral estimator for biomedical signal
> processing applications using novel memory-based FFT processor*, Integration, the VLSI Journal 99 (2024) 102241

lives in this repository. Columns: **Model** = `python/btpsd`, **RTL** = `rtl/` (simulated), **Data** = `data/`
(paper-reported), **Figure** = `docs/figures/` (drawn here) or `docs/paper_figures/` (the paper's own images).

## Sections

| Section | Topic | Model | RTL | Figures and data |
|---|---|---|---|---|
| 1 | Introduction, taxonomy of spectral estimation | `reference.py` (references) | | `paperfig01` |
| 2 | Blackman-Tukey / Welch estimator, 512-point memory-based radix-8 FFT, R8MDC baseline | `r8fft.py`, `schedule.py`, `estimator.py` | `fft512_mem.v`, `r8_butterfly.v` | `fig02`, `fig03`, `fig09`, `fig10`, `paperfig02`, `paperfig03`, Tables 1, 2 |
| 2 | chunk merge with the fractional-delay filter | `fdfilter.py`, `estimator.py` | `psd_engine.v` | `fig04`, `fig05`, `fig11`, `paperfig04` to `paperfig11`, Table 3 |
| 2 | frequency-domain windowing | `windows.py`, `estimator.py` | `psd_engine.v` | `fig07`, `paperfig12` to `paperfig15`, Table 4 |
| 2 | central rotational unit (CRU) | `cordic.py`, `radix2_cordic_fft.py` | `cru_rotate.v`, `cru_rom.v` | `fig06`, `paperfig16`, `paperfig17`, Tables 5, 6 |
| 2 | squarer and averaging | `squarer.py`, `estimator.py` | `mag_sq_odd.v`, `psd_engine.v` | `fig07` |
| 2 | bispectrum | `bispectrum.py` | `bispectrum_engine.v`, `bt_bispectrum_top.v` | `fig08`, `fig13`, `fig14` |
| 3 | experiments: PSD, SQNR, resources, bispectrum | `python/scripts/run_experiments.py` (synthetic signals) | | `fig12`, `fig13`, `fig14`, `results/`, Tables 7 to 10, `paperfig18` to `paperfig32` |
| 4 | conclusion | | | |

## Equations

| Eq. | Content | Model | RTL | Remark |
|---:|---|---|---|---|
| 1 | Blackman-Tukey estimate, correlogram with a lag window = periodogram convolved with the window spectrum | `reference.bt_correlogram_psd`, `bt_via_periodogram_smoothing`, `bt_psd_reference` | | the implemented estimator is Welch's (`PAPER_NOTES.md`, item 1) |
| 2 | valid output from the chunk spectra with 512 final-stage butterflies (classic approach) | `fdfilter.exact_half_bin`, `merge_chunk_spectra` | | the exact reference the six-tap merge is compared with; erratum 1 (range of the outer index) |
| 3, 4 | even and odd bins of the 1024-bin segment spectrum from the two chunk spectra | `fdfilter.interp_half_bin`, `estimator._merge_fixed` | `psd_engine.v` (merge) | |
| 5 | error of the filter response `E = H - H_id` | `fdfilter.fd_frequency_response`, `fd_group_delay` | | analytic bound in `results/fd_merge_snr_vs_taps.csv` |
| 6 | six taps `{0.125, -0.2122, 0.6366, 0.6366, -0.2122, 0.125}` | `fdfilter.fd_kernel`, `HwConfig.fd_int` | `btpsd_consts.vh` (`FD_H0..5`) | integer taps with 12 fraction bits |
| 7 to 11 | periodogram, segmentation `x_k[n] = x[d + kL]`, average over segments | `reference.periodogram`, `welch_psd`, `frame_signal`; `estimator.psd_exact_welch`, `chunk_view` | `psd_engine.v` (accumulate) | |
| 12, 13 | windowing in the frequency domain, `F{x w} = X * W / 2 pi` | `windows.spectral_kernel`, `apply_window_freq`; `estimator._window_fixed` | `psd_engine.v` (window), `btpsd_consts.vh` (`WIN_*`) | |
| 14 to 16 | rotation as a product of `exp(j w_i)` | `cordic.rotate_float`, `cordic_rotate_float` | `cru_rotate.v` | |
| 17 to 21 | micro-rotation recursion with `tg(w)` | `cordic.cordic_rotate_int` | | errata 4 to 7 (text before Eq. 18, Eqs. 17, 19, 20, 21) |
| 22, 23 | dyadic steps `d_i 2^-i` | `cordic.cordic_rotate_int`, `radix2_cordic_fft.py` | | erratum 7 |
| 24, 25 | scale factor `k = prod cos(w_i)` | `cordic.cordic_gain` | `cru_rotate.v` (`inv` in the ROM) | |
| 26 | `sum (2 i - 1) = n^2`, squarer by adding odd numbers | `squarer.square_odds`, `square_serial`, `mag_sq` | `mag_sq_odd.v` | cycle count = `|code|`, `PAPER_NOTES.md` section 4 |
| 27 | averaging by a shift of `K - i` bits | `estimator.psd_fixed` (`avg_shift`) | `psd_engine.v`, `bt_psd_top.v` (`avg_shift`) | |
| 28 | `B(f1, f2) = F(f1) F(f2) F*(f1 + f2)` | `bispectrum.bispectrum_direct`, `bispectrum_fixed` | `bispectrum_engine.v` | |
| 29 | phase relation `psi = phi(f1) + phi(f2) - phi(f1 + f2)` | `bispectrum.bicoherence`, `signals.qpc` | | |
| 30, 31 | magnitude relation; `B` as the 2-D transform of the third-order cumulant | `bispectrum.bispectrum_indirect` | | erratum 8 for (30) |
| 32 | third-order cumulant | `bispectrum.third_order_cumulant` | | erratum 9 |

## Tables

| Table | Content | Where |
|---:|---|---|
| 1, 2 | register-file address sequences | `data/paper_table1_*`, `data/paper_table2_*`; erratum 11 (Table 2, R113), errata 13, 14 (text of Section 2) |
| 3 | SNR versus FD taps | `data/paper_table3_*`; compared with the model in `results/fd_merge_snr_vs_taps.csv` (not reproduced, `PAPER_NOTES.md` section 2.3) |
| 4 | window significant taps | `data/paper_table4_*`; `windows.py` (kernels with 1, 3 and 5 taps) |
| 5 | dyadic micro-rotations | `data/paper_table5_*`; errata 2, 3 |
| 6 | CRU coefficients | `data/table6_printed.csv`; erratum 10; the ROM of the repository is generated (`data/cru_rom_generated.csv`) |
| 7 | test-bench data codes | `data/paper_table7_*` |
| 8 | SQNR at WL = 10 | `data/paper_table8_*`; model counterpart `results/model_sqnr_wl10_synthetic.csv`; erratum 12 (empty columns) |
| 9 | FPGA comparison | `data/paper_table9_*`; no synthesis result of this repository exists |
| 10 | bispectrum SQNR | `data/paper_table10_*`; model counterpart `results/bispectrum_sqnr_vs_wordlength.csv` |

## Figures

The paper's 32 figures are in `docs/paper_figures/`; the manifest there names, for each figure, the page and the
counterpart drawn in this repository. The sheets of `docs/figures/` are drawn from the model and the RTL.
